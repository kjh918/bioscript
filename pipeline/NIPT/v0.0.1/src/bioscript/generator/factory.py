import yaml
import re
import os
from pathlib import Path
from types import SimpleNamespace
from .writers import PythonWriter, BashWriter, NextflowWriter, SnakemakeWriter

# --- [전역 설정 관리] ---
flags = SimpleNamespace(
    PYTHON_SHEBANG = "#!/usr/bin/env python3",
    BASH_SHEBANG   = "#!/bin/bash",
    META_TAG       = "[METADATA]"
)

class SafeDict(dict):
    """템플릿 치환 시 매핑되지 않은 키가 있어도 에러를 내지 않고 그대로 유지함"""
    def __missing__(self, key): return "{" + key + "}"

class BioPipelineFactory:
    def __init__(self, yaml_path):
        with open(yaml_path, 'r', encoding='utf-8') as f:
            self.data = yaml.safe_load(f)
        self.tmpl_dir = Path(__file__).parent / "templates"

    def _build_setup_sh(self, out_path, clean_tool, clean_profile, meta, setup_content):
        file_name = f"setup_{clean_tool}_{clean_profile}.sh"
        target_file = out_path / file_name
        
        script_lines = [
            "#!/bin/bash",
            "# ==========================================",
            f"# Tool Setup : {meta['name']} (v{meta['version']})",
            f"# Profile    : {meta['profile']}",
            f"# Writer     : {meta['writer_name']} ({meta['date']})",
            "# ==========================================",
            "",
            setup_content if setup_content else "# No specific setup required.",
            ""
        ]
        
        target_file.write_text("\n".join(script_lines), encoding='utf-8')
        target_file.chmod(0o755)
        print(f"[*] Created & Executable set: {target_file.name}")

    def render_smart(self, template_text, mapping, meta):
        """
        1단계: 메타데이터 치환 (shebang, threads 등)
        2단계: 스마트 인덴트 치환 (코드 블록 리스트)
        """
        # 1. 메타데이터 1차 치환
        content = template_text.format_map(SafeDict(**meta))

        # 2. 스마트 인덴트 치환
        final_lines = []
        for line in content.splitlines():
            # 라인 전체가 {BLOCK} 형태인 경우 들여쓰기 계산
            match = re.search(r"^(\s*)\{([A-Za-z0-9_]+)\}\s*$", line)
            if match:
                indent, key = match.group(1), match.group(2)
                if key in mapping and isinstance(mapping[key], list):
                    for l in mapping[key]:
                        final_lines.append(f"{indent}{l}")
                    continue
            
            # 인라인 치환 (cmd_line 등 단일 문자열 변수)
            temp_line = line
            for k, v in mapping.items():
                if isinstance(v, str) and f"{{{k}}}" in temp_line:
                    temp_line = temp_line.replace(f"{{{k}}}", v)
            final_lines.append(temp_line)
        
        return "\n".join(final_lines)

    def build_all(self, out_dir):
        """
        모든 언어 템플릿을 기반으로 실행 권한이 부여된 스크립트를 생성함
        파일명 형식: run_[tool]_[profile].[ext]
        """
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        
        # [MODIFIED] 변경 이유: info, tool 블록을 안전하게 가져오기 위해 변수 분리
        info_data = self.data.get('info', {})
        tool_data = self.data.get('tool', {})
        params_data = self.data.get('params', {})
        
        tool_raw = tool_data.get('name', 'Unknown_Tool')
        profile_raw = tool_data.get('profile', 'default')
        
        # [MODIFIED] 변경 이유: 대소문자 혼용(Threads vs threads) 방어 로직
        thread_node = params_data.get('Threads') or params_data.get('threads') or {}
        threads_val = thread_node.get('default', 1) if isinstance(thread_node, dict) else (thread_node or 1)

        # [MODIFIED] 변경 이유: 템플릿의 주석 블록 치환을 위해 info(user, date)와 description 추가
        common_meta = {
            "meta_tag": flags.META_TAG,
            "name": tool_raw,
            "version": tool_data.get('version', '0.1.0'),
            "profile": profile_raw,
            "description": tool_data.get('description', '').strip(),
            "writer_name": info_data.get('user', 'Unknown'),
            "date": info_data.get('date', 'Unknown'),
            "threads": threads_val
        }

        # 파일명 정규화 (소문자화 및 공백 제거)
        clean_tool = tool_raw.lower().replace(" ", "_")
        clean_profile = profile_raw.lower().replace(" ", "_")

        formats = {
            "python": (PythonWriter, "py", flags.PYTHON_SHEBANG),
            "bash": (BashWriter, "sh", flags.BASH_SHEBANG),
            "nextflow": (NextflowWriter, "nf", ""),
            "snakemake": (SnakemakeWriter, "smk", "") # 추가된 부분
        }

        for lang, (WriterClass, ext, shebang) in formats.items():
            tmpl_file = self.tmpl_dir / f"{lang}.template.{ext}"
            if not tmpl_file.exists():
                continue

            meta = {**common_meta, "shebang": shebang}
            writer = WriterClass(self.data)
            blocks = writer.get_blocks()
            
            # 코드 생성
            final_code = self.render_smart(tmpl_file.read_text(), blocks, meta)
            
            # 파일명 결정: run_[tool]_[profile].[ext]
            file_name = f"run_{clean_tool}_{clean_profile}.{ext}"
            target_file = out_path / file_name
            
            # 파일 저장
            target_file.write_text(final_code, encoding='utf-8')
            
            # --- [실행 권한 추가: 755 (rwxr-xr-x)] ---
            target_file.chmod(0o755)
            
            print(f"[*] Created & Executable set: {target_file.name}")
            
        setup_content = tool_data.get('setup', '').strip()
        self._build_setup_sh(out_path, clean_tool, clean_profile, common_meta, setup_content)