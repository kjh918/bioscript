// [METADATA]
// TOOL_NAME = {name}
// THREADS = {threads}

process {process_name} {{
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "{publish_path}", mode: 'copy'

    input:
    {input_block}

    output:
    {output_block}

    script:
    // 로컬 변수 정의 (YAML params 기반)
    {script_vars}
    """
    {cmd_line}
    """
}}