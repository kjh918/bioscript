{shebang}
# {meta_tag}
# TOOL_NAME = {name}
# VERSION = {version}
# THREADS = {threads}

# Tool Info: {name} ({version})
# Profile: {profile}

usage() {{
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    {help_block_req}
    echo ""
    echo "Optional Parameters:"
    {help_block_opt}
    echo "  -h, --help       Show this help message"
    exit 1
}}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
{default_decl_block}

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        {case_block}
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
{validation_block}

# --- [Output Paths] ---
{output_decl_block}

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="{cmd_line}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
{mkdir_block}

eval "$cmd"