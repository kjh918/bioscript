"""
GRCh38 chromosome reference sizes and display lists.
No third-party imports – pure data.
"""

CHROM_SIZES: dict[str, int] = {
    "1":  248_956_422, "2":  242_193_529, "3":  198_295_559,
    "4":  190_214_555, "5":  181_538_259, "6":  170_805_979,
    "7":  159_345_973, "8":  145_138_636, "9":  138_394_717,
    "10": 133_797_422, "11": 135_086_622, "12": 133_275_309,
    "13": 114_364_328, "14": 107_043_718, "15": 101_991_189,
    "16":  90_338_345, "17":  83_257_441, "18":  80_373_285,
    "19":  58_617_616, "20":  64_444_167, "21":  46_709_983,
    "22":  50_818_468, "X":  156_040_895, "Y":   57_227_415,
}

ALL_CHROMS: list[str] = [str(i) for i in range(1, 23)] + ["X", "Y"]
FEMALE_CHROMS: list[str] = [str(i) for i in range(1, 23)] + ["X"]
MALE_CHROMS: list[str] = [str(i) for i in range(1, 23)] + ["X", "Y"]

# Dash-Bio Ideogram annotation track definitions (stable, shared by both apps)
ANNOTATION_TRACKS: list[dict] = [
    {"id": "cnv",  "displayName": "CNV",  "color": "#E53E3E", "shape": "rectangle"},
    {"id": "gene", "displayName": "Gene", "color": "#6B46C1", "shape": "circle"},
]

LEGEND: list[dict] = [{"name": "Legend", "rows": [
    {"name": "Trisomy / Gain",  "color": "#FC8181", "shape": "rectangle"},
    {"name": "Monosomy / Loss", "color": "#90CDF4", "shape": "rectangle"},
    {"name": "Partial Gain",    "color": "#FBD38D", "shape": "rectangle"},
    {"name": "Gene",            "color": "#B794F4", "shape": "circle"},
]}]
