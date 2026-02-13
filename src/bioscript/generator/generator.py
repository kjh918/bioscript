#!/usr/bin/env python3
import re
import yaml
import textwrap
import shlex
from pathlib import Path
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

# [패턴 통일] Generator와 Runner 모두 이 패턴을 사용
TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]+)\]")

@dataclass
class ToolMeta:
    name: str
    version: str
    profile: str

@dataclass
class ScriptConfig:
    tool: ToolMeta
    cmd_line: str
    inputs: Dict[str, Any]
    outputs: Dict[str, Any]
    params: Dict[str, Any]
    output_dir: str

# --- 내부 로직 함수들 ---

def _load_yaml(path: Path) -> Dict[str, Any]:
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f) or {}

def _parse_cfg(cfg: Dict[str, Any]) -> ScriptConfig:
    cmd_line = str(cfg.get("cmd_line", "")).strip()
    
    # [핵심 수정] YAML에 ${VAR}라고 적혀 있어도 [VAR]로 강제 변환하여 저장
    cmd_line = cmd_line.replace("${", "[").replace("}", "]")
    
    tb = cfg.get("tool", {})
    io = cfg.get("io", {})
    return ScriptConfig(
        tool=ToolMeta(
            name=str(tb.get("name", "tool")),
            version=str(tb.get("version", "unknown")),
            profile=str(tb.get("profile", "default"))
        ),
        cmd_line=cmd_line,
        inputs=io.get("inputs", {}),
        outputs=io.get("outputs", {}),
        params=cfg.get("params", {}),
        output_dir=str(cfg.get("output_dir", "./scripts"))
    )

def _get_all_keys(sc: ScriptConfig) -> List[str]:
    keys = set(sc.inputs.keys()) | set(sc.params.keys()) | set(sc.outputs.keys())
    for m in TOKEN_PAT.finditer(sc.cmd_line):
        keys.add(m.group(1))
    return sorted(list(keys))

# --- Runner 생성 함수 ---

def _python_runner_text(sc: ScriptConfig, outputs_meta: Dict, all_keys: List[str], required: Dict[str, bool], defaults: Dict[str, str]) -> str:
    # 생성된 .py 파일 안의 CMD_LINE이 [VAR] 형태가 되도록 보장됨
    return textwrap.dedent(f"""\
        #!/usr/bin/env python3
        import argparse, json, re, subprocess, os, shlex
        from pathlib import Path

        TOKEN_PAT = re.compile(r"\\[([A-Za-z0-9_]*)\\]")
        CMD_LINE = {sc.cmd_line!r}
        REQUIRED_KEYS = {[k for k, v in required.items() if v]!r}
        DEFAULTS = {defaults!r}
        OUTPUT_KEYS = {list(sc.outputs.keys())!r}

        def render(s, ctx):
            def repl(m):
                key = m.group(1)
                return str(ctx.get(key, m.group(0))).strip()
            res = s
            for _ in range(3):
                if "[" not in res: break
                res = TOKEN_PAT.sub(repl, res)
            return " ".join(res.split())

        def main():
            parser = argparse.ArgumentParser()
            for k in {all_keys!r}:
                parser.add_argument(f"--{{k}}", default=DEFAULTS.get(k, ""))
            parser.add_argument("--cwd", default=".")
            parser.add_argument("--emit-outputs", default="")
            
            args = parser.parse_args()
            ctx = vars(args)

            for rk in REQUIRED_KEYS:
                if not ctx.get(rk):
                    print(f"Error: Missing --{{rk}}")
                    exit(1)

            # Output 경로 렌더링
            for k in OUTPUT_KEYS:
                ctx[k] = render(ctx[k], ctx)

            final_cmd = render(CMD_LINE, ctx)
            print(f"\\n[RUNNING CMD]\\n{{final_cmd}}\\n")

            os.makedirs(args.cwd, exist_ok=True)
            proc = subprocess.run(final_cmd, shell=True, cwd=args.cwd)
            
            if args.emit_outputs:
                out_data = {{k: ctx[k] for k in OUTPUT_KEYS}}
                Path(args.emit_outputs).write_text(json.dumps(out_data, indent=2))
            
            exit(proc.returncode)

        if __name__ == "__main__":
            main()
    """)
def _bash_runner_text(sc: ScriptConfig, outputs_meta: Dict, all_keys: List[str], defaults: Dict[str, str]) -> str:
    # 1. 기본값 선언 (한 줄씩 할당)
    default_lines = []
    for k in all_keys:
        val = defaults.get(k, "")
        default_lines.append(f"{k}={shlex.quote(str(val))}")
    default_block = "\n".join(default_lines)

    # 2. CLI Case 문 (인자 파싱)
    case_lines = []
    for k in all_keys:
        case_lines.append(f'    --{k}) {k}="$2"; shift 2;;')
    case_block = "\n".join(case_lines)

    # 3. 핵심: [변수]를 실제 값으로 치환하는 Bash 로직
    # Bash의 ${s//pattern/string} 문법에서 [ ]는 특수문자이므로 이스케이프가 생명입니다.
    repl_lines = []
    for k in all_keys:
        # Bash에서 \[VAR\] 를 찾아서 $VAR(의 값)으로 바꿈
        repl_lines.append(f'    s="${{s//\\\[{k}\\\]/${{{k}}}}}"')
    repl_block = "\n".join(repl_lines)

    # 4. Output 변수들 렌더링 명령 (세미콜론으로 연결)
    out_render_lines = []
    for k in outputs_meta.keys():
        out_render_lines.append(f'{k}=$(render "${{{k}}}")')
    out_render_block = "\n".join(out_render_lines)

    # 5. 최종 템플릿 조립 (들여쓰기 꼬임 방지를 위해 왼쪽 정렬 유지)
    script_content = f"""#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
{default_block}
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
{case_block}
    --cwd) CWD="$2"; shift 2;;
    -h|--help) echo "Usage: $0 [options]"; exit 0;;
    *) echo "Unknown argument: $1" >&2; exit 1;;
  esac
done

# --- Template Engine ---
render() {{
  local s="$1"
  # 중첩된 변수 치환을 위해 3회 반복 (예: [A] -> [B] -> value)
  for i in {{1..3}}; do
{repl_block}
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}}

# --- Finalize Outputs & Command ---
{out_render_block}
CMD_LINE={shlex.quote(sc.cmd_line)}
CMD=$(render "$CMD_LINE")

echo -e "\\n[RUNNING CMD]\\n$CMD\\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
"""
    return script_content

# --- 메인 함수 ---

def make_runners(config_path: Path, outdir_override: Optional[Path] = None) -> None:
    """
    YAML 설정을 로드하여 run_{tool}_{profile}.py 및 .sh 파일을 생성합니다.
    """
    cfg_raw = _load_yaml(config_path)
    sc = _parse_cfg(cfg_raw)
    
    # [1] 출력 디렉토리 결정
    # CLI에서 넘어온 override가 있으면 우선 사용, 없으면 YAML 설정값 사용
    outdir = Path(outdir_override) if outdir_override else Path(sc.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # [2] 메타데이터 정리
    out_meta = {k: (v if isinstance(v, dict) else {"default": str(v)}) for k, v in sc.outputs.items()}
    keys = _get_all_keys(sc)
    req = {k: bool(v.get("required")) for k, v in sc.inputs.items() if isinstance(v, dict)}
    dfs = {k: str(v.get("default", "")) for k, v in {**sc.params, **out_meta}.items() if isinstance(v, dict)}
    
    # [3] 파일명 생성 (tool_name + profile)
    # 특수문자나 공백을 언더바(_)로 치환하여 안전한 파일명 생성
    safe_name = re.sub(r"[^A-Za-z0-9_]+", "_", sc.tool.name).strip("_")
    safe_profile = re.sub(r"[^A-Za-z0-9_]+", "_", sc.tool.profile).strip("_")
    
    base_filename = f"run_{safe_name}_{safe_profile}"
    py_path = outdir / f"{base_filename}.py"
    sh_path = outdir / f"{base_filename}.sh"
    
    # [4] 파일 쓰기
    py_path.write_text(_python_runner_text(sc, out_meta, keys, req, dfs), encoding="utf-8")
    sh_path.write_text(_bash_runner_text(sc, out_meta, keys, dfs), encoding="utf-8")
    sh_path.chmod(0o755)
    
    print(f"[*] Generated runners in: {outdir}")
    print(f"    - {py_path.name}")
    print(f"    - {sh_path.name}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        make_runners(Path(sys.argv[1]))