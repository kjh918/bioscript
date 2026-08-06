import os
import pysam
import pandas as pd
import numpy as np

def log(msg): 
    print(f"[INFO] {msg}", flush=True)

def ensure_dir(path): 
    os.makedirs(path, exist_ok=True)

def chrom_key(c):
    s = str(c).replace("chr", "")
    if s == "X": return 23
    if s == "Y": return 24
    try: return int(s)
    except: return 99

def sort_chroms(chroms):
    return sorted(chroms, key=chrom_key)


def safe_log2fc(val, ref):
    return np.log2((val + 1e-9) / (ref + 1e-9))