
import pysam


class Fragment:
    """
    2개의 Paired-end read를 묶어 관리하는 Domain Object.
    """
    def __init__(self, read_a: pysam.AlignedSegment, read_b: pysam.AlignedSegment):
        # [MODIFIED] 객체 생성 시 무조건 2개의 리드가 입력되어야 하며, 유효하지 않으면 Exception 발생
        if read_a is None or read_b is None:
            raise ValueError("Fragment 객체는 반드시 2개의 read(AlignedSegment)를 가져야 합니다.")
        if read_a.query_name != read_b.query_name:
            raise ValueError(f"Read ID가 일치하지 않습니다: {read_a.query_name} vs {read_b.query_name}")

        self.read_id = read_a.query_name
        
        # Read1, Read2 명확히 분리하여 할당
        if read_a.is_read1:
            self.read1, self.read2 = read_a, read_b
        else:
            self.read1, self.read2 = read_b, read_a

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def midpoint(self) -> int:
        # [MODIFIED] Fragment의 중심점을 반환 (Coverage 기준점을 Read start가 아닌 Midpoint로 전환하기 위함)
        return (self.start + self.end) // 2

    def get_gc_content(self) -> float:
        """Fragment 전체 구간의 GC 비율 계산"""
        seq = self.fasta.fetch(self.chrom, self.start, self.end).upper()
        if not seq:
            return 0.0
        gc_count = seq.count('G') + seq.count('C') + seq.count('S')
        return gc_count / len(seq)

    def get_end_motifs(self, motif_len: int) -> tuple:
        """5' 및 3' End Motif 추출 (예: 4-mer)"""
        # [MODIFIED] Reference 서열을 기반으로 Motif 추출. (실제 데이터에 따라 Read 자체 서열을 쓸 수도 있음)
        motif_5 = self.fasta.fetch(self.chrom, self.start, self.start + motif_len).upper()
        motif_3 = self.fasta.fetch(self.chrom, self.end - motif_len, self.end).upper()
        return motif_5, motif_3

    def is_valid_midpoint(self, bed_start: int, bed_end: int) -> bool:
        """Fragment의 중심점이 지정된 BED 구간 내에 존재하는지 확인"""
        return bed_start <= self.midpoint < bed_end