import os

def log(msg): 
    print(f"[INFO] {msg}", flush=True)

def ensure_dir(path): 
    os.makedirs(path, exist_ok=True)