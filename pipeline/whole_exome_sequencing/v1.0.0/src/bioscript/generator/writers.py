import re
import os

class BaseWriter:
    def __init__(self, data):
        self.data = data
        self.io = data.get('io', {})
        self.params = data.get('params', {})
        self.cmd = data.get('cmd_line', "").strip().replace('"',"'").replace('\t', '\\t')

    def _get_merged_vars(self, include_outputs=True):
        """치환을 위한 전체 변수 맵 (inputs + params + optional outputs)"""
        v = {}
        sources = [self.io.get('inputs', {}), self.params]
        if include_outputs:
            sources.append(self.io.get('outputs', {}))
            
        for d in sources:
            for k, val in d.items():
                v[k] = val.get('default', '') if isinstance(val, dict) else val
        return v

    def _get_full_info(self):
        """help, required 정보를 포함한 맵 (argparse/help/case문 생성용)"""
        info = {}
        # 도움말과 인자 선언에는 outputs를 제외한 inputs와 params만 사용
        for d in [self.io.get('inputs', {}), self.io.get('outputs', {}), self.params]:
            for k, val in d.items():
                info[k] = val if isinstance(val, dict) else {"default": val}
        return info

    def safe_replace(self, cmd, vars_map, prefix="", suffix=""):
        """정의된 변수만 [Key] -> {Key} 등으로 치환하여 samtools [NM] 등 보호"""
        def sub_func(m):
            var_name = m.group(1)
            if var_name in vars_map:
                return f"{prefix}{var_name}{suffix}"
            return f"[{var_name}]"
        
        pattern = r"\[([A-Za-z0-9_]+)\]"
        return re.sub(pattern, sub_func, cmd)

# ---------------------------------------------------------

class BashWriter(BaseWriter):
    def get_blocks(self):
        vars_map = self._get_merged_vars(include_outputs=True)
        full_info = self._get_full_info()
        
        h_req, h_opt = [], []
        decl_lines, case_lines, validation_lines = [], [], []
        
        # 1. Inputs & Params 처리
        for k, info in full_info.items():
            h = info.get('help', 'No description')
            is_req = info.get('required', False)
            df = info.get('default', '')
            
            # Help Block
            if is_req:
                h_req.append(f'echo "  --{k:15} {h}"')
                validation_lines.append(f'if [[ -z "${{{k}:-}}" ]]; then echo "Error: --{k} is required"; usage; fi')
            else:
                h_opt.append(f'echo "  --{k:15} {h} (Default: {df})"')
            
            # Decls & Case
            decl_lines.append(f'{k}="{df}"')
            case_lines.append(f'--{k}) {k}="$2"; shift 2 ;;')

        # 2. Output Decls (내부 계산 경로)
        out_decls = []
        for k, v in self.io.get('outputs', {}).items():
            val = v.get('default', '') if isinstance(v, dict) else v
            val_replaced = self.safe_replace(str(val), vars_map, "${", "}")
            out_decls.append(f'{k}="{val_replaced}"')

        mkdir_targets = list(set([k for k in vars_map.keys() if "Dir" in k or "Path" in k.lower()] + list(self.io.get('outputs', {}).keys())))
        
        mkdir_lines = []
        for k in mkdir_targets:
            # bash 문법: 변수가 존재할 때만 실행, 확장자가 있으면 dirname 추출
            mkdir_lines.append(f'if [[ -n "${{{k}:-}}" ]]; then')
            mkdir_lines.append(f'  if [[ "${{{k}}}" == *.* ]]; then mkdir -p "$(dirname "${{{k}}}")"; else mkdir -p "${{{k}}}"; fi')
            mkdir_lines.append(f'fi')

        return {
            "help_block_req": h_req or ['echo "  (None)"'],
            "help_block_opt": h_opt or ['echo "  (None)"'],
            "default_decl_block": decl_lines,
            "case_block": case_lines,
            "validation_block": validation_lines,
            "output_decl_block": out_decls,
            "mkdir_block": mkdir_lines,
            "cmd_line": self.safe_replace(self.cmd, vars_map, "${", "}")
        }

# ---------------------------------------------------------
class PythonWriter(BaseWriter):
    def get_blocks(self):
        vars_map = self._get_merged_vars(include_outputs=True)
        full_info = self._get_full_info()
        output_keys = list(self.io.get('outputs', {}).keys())
        
        # 1. Argparse 블록
        arg_lines = []
        for k, info in full_info.items():
            h = str(info.get('help', 'No description')).replace("'", "\\'")
            req = info.get('required', False)
            df = info.get('default', '')
            arg_lines.append(f"parser.add_argument('--{k}', required={req}, default='{df}', help='{h} (Default: {df})')")

        # 2. Output Fallback 블록 (CLI로 안 들어왔을 경우 템플릿 치환)
        out_lines = []
        for k, v in self.io.get('outputs', {}).items():
            val = v.get('default', '') if isinstance(v, dict) else v
            if val:
                val_replaced = self.safe_replace(str(val), vars_map, "{", "}")
                # args.k가 비어있으면 템플릿으로 계산된 경로 할당
                out_lines.append(f'if not {k}:\n    {k} = f"{val_replaced}"')

        dir_keys = [k for k in vars_map.keys() if "Dir" in k or "Path" in k.lower()]
        mkdir_lines = []
        for k in set(dir_keys + output_keys):
            # os.path.splitext를 통해 파일인지 폴더인지 구분하여 부모 디렉토리 확보
            mkdir_lines.append(f"if {k}:")
            mkdir_lines.append(f"    _tgt = os.path.dirname({k}) if os.path.splitext({k})[1] else {k}")
            mkdir_lines.append(f"    if _tgt: os.makedirs(_tgt, exist_ok=True)")

        # 4. cmd_line 치환 로직
        replaced_cmd = self.safe_replace(self.cmd, vars_map, "{", "}")

        return {
            "argparse_block": arg_lines,
            "var_decl_block": [f"{k} = args.{k}" for k in full_info.keys()],
            "output_decl_block": out_lines,
            "mkdir_block": mkdir_lines,  # 스마트 들여쓰기를 위해 리스트 반환
            "cmd_line": replaced_cmd
        }
# ---------------------------------------------------------

class NextflowWriter(BaseWriter):
    def get_blocks(self):
        vars_map = self._get_merged_vars(include_outputs=True)
        output_keys = self.io.get('outputs', {}).keys()
        
        in_elements = ["val(sample_id)"]
        for k in self.io.get('inputs', {}).keys():
            in_elements.append(f"path({k})" if any(x in k.lower() for x in ["fastq", "bam", "vcf"]) else f"val({k})")
        
        return {
            "process_name": str(self.data['tool']['name']).upper(),
            "publish_path": "${params.outdir}/${sample_id}",
            "input_block": [f"tuple {', '.join(in_elements)}"],
            "output_block": [f'path "*", emit: {k}' for k in output_keys],
            "script_vars": [f'def {k} = params.{k} ?: "{v}"' for k, v in vars_map.items() if k not in output_keys],
            "cmd_line": [self.safe_replace(self.cmd, vars_map, "${", "}")]
        }

# ---------------------------------------------------------

class SnakemakeWriter(BaseWriter):
    def get_blocks(self):
        input_keys = self.io.get('inputs', {}).keys()
        output_keys = self.io.get('outputs', {}).keys()
        
        def snakemake_sub(m):
            var_name = m.group(1)
            if var_name.lower() == "threads": return "{threads}"
            if var_name in input_keys: return f"{{input.{var_name}}}"
            if var_name in output_keys: return f"{{output.{var_name}}}"
            return f"{{params.{var_name}}}"

        final_cmd = re.sub(r"\[([A-Za-z0-9_]+)\]", snakemake_sub, self.cmd)
        
        return {
            "rule_name": str(self.data['tool']['name']).lower().replace(" ", "_"),
            "input_block": [f'{k} = "{v.get("default","") if isinstance(v,dict) else v}"' for k, v in self.io.get('inputs', {}).items()],
            "output_block": [f'{k} = "{v.get("default","") if isinstance(v,dict) else v}"' for k, v in self.io.get('outputs', {}).items()],
            "param_block": [f'{k} = "{v.get("default","") if isinstance(v,dict) else v}"' for k, v in self.params.items()],
            "cmd_line": [line.strip() for line in final_cmd.splitlines() if line.strip()]
        }