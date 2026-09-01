# {meta_tag}
# TOOL_NAME = {name}
# VERSION = {version}
# THREADS = {threads}

rule {rule_name}:
    input:
        {input_block}
    output:
        {output_block}
    params:
        {param_block}
    threads: {threads}
    shell:
        """
        {cmd_line}
        """