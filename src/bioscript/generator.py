from __future__ import annotations

import re
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

# [Var] only
TOKEN_PAT = re.compile(r"\[([A-Za-z_][A-Za-z0-9_]*)\]")


def _infer_option_value_template(cmd_line: str, key: str) -> Optional[str]:
    m = re.search(rf"(?:^|\s)--{re.escape(key)}\s+([^\s]+)", cmd_line)
    if not m:
        return None
    return m.group(1).strip()


@dataclass(frozen=True)
class ToolMeta:
    name: str
    version: str
    profile: str


@dataclass(frozen=True)
class IOSpec:
    inputs: Dict[str, Dict[str, Any]]
    outputs: Dict[str, Dict[str, Any]]


@dataclass(frozen=True)
class ScriptConfig:
    tool: ToolMeta
    cmd_line: str
    io: IOSpec
    params: Dict[str, Dict[str, Any]]
    output_dir: str


def _load_yaml(path: Path) -> Dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", (s or "").strip())


def _extract_tokens(s: str) -> List[str]:
    seen: List[str] = []
    for m in TOKEN_PAT.finditer(s or ""):
        key = m.group(1)
        if key and key not in seen:
            seen.append(key)
    return seen


def _parse_tool(cfg: Dict[str, Any]) -> ToolMeta:
    tb = cfg.get("tool") or {}
    name = (tb.get("name") or "").strip()
    if not name:
        raise ValueError("config.tool.name is required")
    version = str(tb.get("version") or "unknown").strip()
    profile = str(tb.get("profile") or "default").strip()
    return ToolMeta(name=name, version=version, profile=profile)


def _parse_cfg(cfg: Dict[str, Any]) -> ScriptConfig:
    cmd_line = str(cfg.get("cmd_line") or "").strip()
    if not cmd_line:
        raise ValueError("config.cmd_line is required")

    io = cfg.get("io") or {}
    inputs = io.get("inputs") or {}
    outputs = io.get("outputs") or {}
    if not isinstance(inputs, dict) or not isinstance(outputs, dict):
        raise ValueError("config.io.inputs and config.io.outputs must be dict")

    params = cfg.get("params") or {}
    if not isinstance(params, dict):
        raise ValueError("config.params must be dict")

    return ScriptConfig(
        tool=_parse_tool(cfg),
        cmd_line=cmd_line,
        io=IOSpec(inputs=inputs, outputs=outputs),
        params=params,
        output_dir=str(cfg.get("output_dir") or "./generated"),
    )


def _outdir(sc: ScriptConfig, outdir_override: Optional[Path]) -> Path:
    if outdir_override:
        outdir_override.mkdir(parents=True, exist_ok=True)
        return outdir_override
    d = Path(sc.output_dir)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _filenames(tool: ToolMeta) -> Tuple[str, str]:
    t = _safe_name(tool.name)
    v = _safe_name(tool.version)
    p = _safe_name(tool.profile)
    return f"run_{t}_{p}.py", f"run_{t}_{p}.sh"


def _normalize_outputs(sc: ScriptConfig) -> Dict[str, Dict[str, Any]]:
    out_meta: Dict[str, Dict[str, Any]] = {}
    for name, meta in sc.io.outputs.items():
        m = meta if isinstance(meta, dict) else {}
        norm = dict(m)
        if "default" not in norm or norm.get("default") in (None, ""):
            inferred = _infer_option_value_template(sc.cmd_line, name)
            if inferred:
                norm["default"] = inferred
        out_meta[name] = norm
    return out_meta


def _collect_defaults(sc: ScriptConfig, outputs_meta: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    defaults: Dict[str, str] = {}

    for k, meta in sc.params.items():
        if isinstance(meta, dict) and meta.get("default") is not None:
            defaults[k] = str(meta.get("default"))

    for k, meta in outputs_meta.items():
        if isinstance(meta, dict) and meta.get("default") is not None:
            defaults[k] = str(meta.get("default"))

    return defaults


def _collect_required(sc: ScriptConfig, outputs_meta: Dict[str, Dict[str, Any]]) -> Dict[str, bool]:
    required: Dict[str, bool] = {}
    for k, meta in sc.io.inputs.items():
        required[k] = bool(isinstance(meta, dict) and meta.get("required") is True)
    for k, meta in outputs_meta.items():
        required[k] = bool(isinstance(meta, dict) and meta.get("required") is True)
    return required


def _collect_all_arg_keys(sc: ScriptConfig, outputs_meta: Dict[str, Dict[str, Any]]) -> List[str]:
    keys: List[str] = []

    for d in (sc.io.inputs, sc.params, outputs_meta):
        for k in d.keys():
            if k not in keys:
                keys.append(k)

    for t in _extract_tokens(sc.cmd_line):
        if t not in keys:
            keys.append(t)

    for _, meta in outputs_meta.items():
        dv = meta.get("default") if isinstance(meta, dict) else None
        if isinstance(dv, str):
            for t in _extract_tokens(dv):
                if t not in keys:
                    keys.append(t)

    return keys


def _python_runner_text(
    sc: ScriptConfig,
    outputs_meta: Dict[str, Dict[str, Any]],
    all_keys: List[str],
    required: Dict[str, bool],
    defaults: Dict[str, str],
) -> str:
    return f"""#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
from pathlib import Path
from typing import Dict

TOKEN_PAT = re.compile(r"\\[([A-Za-z_][A-Za-z0-9_]*)\\]")

CMD_LINE = {sc.cmd_line!r}
ALL_KEYS = {all_keys!r}
REQUIRED = {required!r}
DEFAULTS = {defaults!r}
OUTPUTS_META = {outputs_meta!r}

def render(s: str, ctx: Dict[str, str]) -> str:
    def repl(m: re.Match) -> str:
        key = m.group(1)
        return (ctx.get(key) or "").strip()
    out = TOKEN_PAT.sub(repl, s)
    out = re.sub(r"\\s+", " ", out).strip()
    return out

def finalize_outputs(ctx: Dict[str, str]) -> Dict[str, str]:
    outputs: Dict[str, str] = {{}}
    for name, meta in OUTPUTS_META.items():
        meta = meta or {{}}
        cur = (ctx.get(name) or "").strip()
        if cur:
            outputs[name] = cur
            continue

        dv = meta.get("default")
        if isinstance(dv, str) and dv.strip():
            val = render(dv, ctx)
            outputs[name] = val
            ctx[name] = val
            continue

        if bool(meta.get("required") is True):
            raise SystemExit(f"Missing required output --{{name}} (no default/template inferred)")
    return outputs

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Runner: {sc.tool.name} {sc.tool.version} ({sc.tool.profile})",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    for k in ALL_KEYS:
        if REQUIRED.get(k, False) and k not in OUTPUTS_META:
            ap.add_argument(f"--{{k}}", required=True, default="")
        else:
            ap.add_argument(f"--{{k}}", required=False, default=DEFAULTS.get(k, ""))

    ap.add_argument("--cwd", default=".")
    ap.add_argument("--print-cmd", action="store_true")
    ap.add_argument("--print-outputs", action="store_true")
    ap.add_argument("--emit-outputs", default="", help="write outputs json to file path")

    args = ap.parse_args()
    ctx = {{k: ("" if v is None else str(v)) for k, v in vars(args).items()
            if k not in ("cwd","print_cmd","print_outputs","emit_outputs")}}

    for k, dv in DEFAULTS.items():
        if not ctx.get(k, "").strip():
            ctx[k] = str(dv)

    outputs = finalize_outputs(ctx)

    cmd = render(CMD_LINE, ctx)

    if args.print_cmd:
        print(cmd)
        return

    proc = subprocess.run(shlex.split(cmd), cwd=str(Path(args.cwd)))
    rc = proc.returncode

    if args.print_outputs:
        print(json.dumps(outputs, ensure_ascii=False))

    if args.emit_outputs:
        Path(args.emit_outputs).write_text(json.dumps(outputs, ensure_ascii=False, indent=2), encoding="utf-8")

    raise SystemExit(rc)

if __name__ == "__main__":
    main()
"""


def _bash_runner_text(
    sc: ScriptConfig,
    outputs_meta: Dict[str, Dict[str, Any]],
    all_keys: List[str],
    required: Dict[str, bool],
    defaults: Dict[str, str],
) -> str:
    lines_defaults: List[str] = []
    for k in all_keys:
        if k in defaults:
            lines_defaults.append(f"{k}={shlex.quote(str(defaults[k]))}")
        else:
            lines_defaults.append(f"{k}=''")

    lines_req_inputs: List[str] = []
    for k, req in required.items():
        if req and k not in outputs_meta:
            lines_req_inputs.append(f'[[ -n "${{{k}}}" ]] || {{ echo "Missing required --{k}" >&2; exit 2; }}')

    case_lines = [f"    --{k}) {k}=\"$2\"; shift 2;;" for k in all_keys]

    outputs_default_block: List[str] = []
    for out_name, meta in outputs_meta.items():
        dv = meta.get("default") if isinstance(meta, dict) else None
        req = bool(isinstance(meta, dict) and meta.get("required") is True)

        if isinstance(dv, str) and dv.strip():
            outputs_default_block.append(f"""
if [[ -z "${{{out_name}}}" ]]; then
  __tmp={shlex.quote(dv)}
  {out_name}="$(render "$__tmp")"
fi
""".strip())

        if req:
            outputs_default_block.append(f'[[ -n "${{{out_name}}}" ]] || {{ echo "Missing required output {out_name}" >&2; exit 2; }}')

    output_keys_arr = " ".join(outputs_meta.keys())

    replace_lines = []
    for t in all_keys:
        replace_lines.append(f's="${{s//[{t}]/${{{t}}}}}"')
    replace_block = "\n  ".join(replace_lines)

    emit_block = r"""
emit_outputs() {
  local out=""
  for k in "${OUTPUT_KEYS[@]}"; do
    out+="${k}\t${!k}\n"
  done
  if [[ "$PRINT_OUTPUTS" -eq 1 ]]; then
    printf "%b" "$out"
  fi
  if [[ -n "$EMIT_OUTPUTS" ]]; then
    printf "%b" "$out" > "$EMIT_OUTPUTS"
  fi
}
""".strip()

    return f"""#!/usr/bin/env bash
set -euo pipefail

# Runner: {sc.tool.name} {sc.tool.version} ({sc.tool.profile})

# defaults (override by CLI)
{chr(10).join(lines_defaults)}

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=({output_keys_arr})

usage() {{
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}}

while [[ $# -gt 0 ]]; do
  case "$1" in
{chr(10).join(case_lines)}
    --cwd) CWD="$2"; shift 2;;
    --print-cmd) PRINT_CMD=1; shift;;
    --print-outputs) PRINT_OUTPUTS=1; shift;;
    --emit-outputs) EMIT_OUTPUTS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

# required inputs
{chr(10).join(lines_req_inputs)}

render() {{
  local s="$1"
  {replace_block}
  s="$(echo "$s" | tr -s ' ' | sed 's/^ *//; s/ *$//')"
  echo "$s"
}}

{emit_block}

# finalize outputs
{chr(10).join(outputs_default_block)}

CMD_LINE={shlex.quote(sc.cmd_line)}
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \\"$CWD\\" && $CMD"

emit_outputs
"""


def make_runners(config_path: Path, outdir_override: Optional[Path] = None) -> None:
    cfg_raw = _load_yaml(config_path)
    sc = _parse_cfg(cfg_raw)

    outputs_meta = _normalize_outputs(sc)
    required = _collect_required(sc, outputs_meta)
    defaults = _collect_defaults(sc, outputs_meta)
    all_keys = _collect_all_arg_keys(sc, outputs_meta)

    outdir = _outdir(sc, outdir_override)
    py_name, sh_name = _filenames(sc.tool)

    py_path = outdir / py_name
    sh_path = outdir / sh_name

    py_path.write_text(_python_runner_text(sc, outputs_meta, all_keys, required, defaults), encoding="utf-8")
    sh_path.write_text(_bash_runner_text(sc, outputs_meta, all_keys, required, defaults), encoding="utf-8")
    sh_path.chmod(0o755)

    print("generated:", py_path, sh_path)
