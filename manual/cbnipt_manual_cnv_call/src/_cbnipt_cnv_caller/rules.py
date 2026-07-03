import yaml
from pathlib import Path

# config.yaml 경로 동적 탐색 (이 파일이 있는 폴더 기준)
yaml_path = Path(__file__).parent / "config.yaml"

with open(yaml_path, "r", encoding="utf-8") as f:
    _RAW = yaml.safe_load(f)

CFG = _RAW["CFG"]
CALL_COLORS = _RAW["CALL_COLORS"]