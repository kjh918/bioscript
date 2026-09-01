{shebang}
# {meta_tag}
# TOOL_NAME = {name}
# VERSION = {version}
# THREADS = {threads}
# PROFILE = {profile}

"""
Tool: {name} ({version})
Profile: {profile}
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="{name} Analysis Runner")
    
    # --- [Argument Parsing] ---
    {argparse_block}
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    {var_decl_block}

    # --- [Output Paths] ---
    {output_decl_block}

    # --- [Command Execution] ---
    cmd = f"{cmd_line}"
    
    print(f"\\n[RUNNING]\\n{{cmd}}\\n")
    
    {mkdir_block}
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()