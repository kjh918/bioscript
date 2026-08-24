"""
Built-in demo SampleData used when no --sample TSV is supplied.
Mimics a female NIPT-positive case with multiple CNV findings.
"""

from .core.models import SampleData, CnvEvent, GeneAnnotation
from .core.reference import CHROM_SIZES


def make_demo_sample() -> SampleData:
    return SampleData(
        id="DEMO_001",
        sex="female",
        events=[
            CnvEvent(chr="21", type="trisomy",      cn=3, iscn="+21",
                     start=1,          stop=CHROM_SIZES["21"],  color="#FC8181"),
            CnvEvent(chr="X",  type="monosomy",     cn=1, iscn="-X",
                     start=1,          stop=CHROM_SIZES["X"],   color="#90CDF4"),
            CnvEvent(chr="5",  type="partial_loss", cn=1, iscn="del(5p)",
                     start=1,          stop=45_368_512,          color="#90CDF4"),
            CnvEvent(chr="17", type="partial_gain", cn=3, iscn="dup(17q)",
                     start=43_044_000, stop=CHROM_SIZES["17"],  color="#FBD38D"),
        ],
        genes=[
            GeneAnnotation(chr="21", name="DYRK1A", start=37_700_000, stop=37_865_335, color="#B794F4"),
            GeneAnnotation(chr="21", name="RUNX1",  start=34_787_801, stop=36_004_954, color="#B794F4"),
            GeneAnnotation(chr="17", name="TP53",   start=7_661_779,  stop=7_687_538,  color="#68D391"),
            GeneAnnotation(chr="17", name="BRCA1",  start=43_044_295, stop=43_125_483, color="#68D391"),
            GeneAnnotation(chr="17", name="ERBB2",  start=39_687_914, stop=39_730_426, color="#68D391"),
            GeneAnnotation(chr="5",  name="CTNND2", start=11_080_000, stop=11_690_000, color="#68D391"),
            GeneAnnotation(chr="X",  name="MECP2",  start=154_021_573,stop=154_137_103,color="#68D391"),
            GeneAnnotation(chr="X",  name="DMD",    start=31_119_221, stop=33_339_609, color="#68D391"),
            GeneAnnotation(chr="13", name="BRCA2",  start=32_315_086, stop=32_400_266, color="#68D391"),
            GeneAnnotation(chr="7",  name="CFTR",   start=117_480_025,stop=117_668_665,color="#68D391"),
        ],
    )
