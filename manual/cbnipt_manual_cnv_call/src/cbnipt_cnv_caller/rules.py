import yaml
from pathlib import Path

# rules.yaml 경로 동적 탐색 (main.py가 있는 폴더 기준)
yaml_path = Path(__file__).parent / "config.yaml"

with open(yaml_path, "r", encoding="utf-8") as f:
    CFG = yaml.safe_load(f)['CFG']

CALL_COLORS = {"NORMAL": "#4CAF50", "SUSPICIOUS": "#FF9800", "ABNORMAL": "#F44336"}