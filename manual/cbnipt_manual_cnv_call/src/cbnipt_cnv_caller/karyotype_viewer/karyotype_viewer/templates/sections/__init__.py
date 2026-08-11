from .s01_header      import SECTION_HEADER
from .s02_banner      import SECTION_BANNER
from .s03_patient     import SECTION_PATIENT
from .s04_sample      import SECTION_SAMPLE
from .s05_events      import SECTION_EVENTS
from .s06_karyotype   import SECTION_KARYOTYPE
from .s07_chrom_detail import SECTION_CHROM_DETAIL
from .s08_findings    import SECTION_FINDINGS
from .s09_footer      import SECTION_FOOTER

ALL_SECTIONS = [
    SECTION_HEADER,
    SECTION_PATIENT,
    SECTION_SAMPLE,
    SECTION_EVENTS,
    SECTION_KARYOTYPE,
    SECTION_FINDINGS,
    SECTION_FOOTER,
]

__all__ = [
    "SECTION_HEADER", "SECTION_BANNER", "SECTION_PATIENT",
    "SECTION_SAMPLE", "SECTION_EVENTS", "SECTION_KARYOTYPE",
    "SECTION_CHROM_DETAIL", "SECTION_FINDINGS", "SECTION_FOOTER",
    "ALL_SECTIONS",
]