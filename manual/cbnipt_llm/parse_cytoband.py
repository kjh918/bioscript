import json
from pathlib import Path


def add_cytobands_to_json(
    chromosome_json_path: str,
    cytoband_txt_path: str,
    output_json_path: str,
):
    """
    기존 chromosome JSON에 UCSC cytoband 정보를 추가한다.

    입력 cytoband TXT 형식:
        chr1    0          2300000    p36.33    gneg
        chr1    2300000    5300000    p36.32    gpos25
        ...

    결과:
        chromosomes
          └─ chr1
              ├─ p_arm
              ├─ centromere
              ├─ q_arm
              └─ cytobands
                   ├─ p36.33
                   ├─ p36.32
                   └─ ...
    """

    # ---------------------------------------------------------
    # 1. 기존 chromosome JSON 로드
    # ---------------------------------------------------------
    json_path = Path(chromosome_json_path)

    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    chromosomes = data.setdefault("chromosomes", {})

    # ---------------------------------------------------------
    # 2. cytoband TXT 파싱
    # ---------------------------------------------------------
    txt_path = Path(cytoband_txt_path)

    band_count = 0

    with txt_path.open("r", encoding="utf-8") as f:

        for line_no, line in enumerate(f, start=1):

            line = line.strip()

            # 빈 줄 / comment 무시
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 5:
                print(
                    f"[warning] line {line_no}: "
                    f"필드 부족 → {line}"
                )
                continue

            chrom, start, end, band, stain = parts[:5]

            # 우리가 관리하는 chromosome만 사용
            if chrom not in chromosomes:
                print(
                    f"[skip] {chrom}: "
                    f"chromosome JSON에 없음"
                )
                continue

            try:
                start = int(start)
                end = int(end)

            except ValueError:
                print(
                    f"[warning] line {line_no}: "
                    f"좌표 변환 실패 → {line}"
                )
                continue

            # -------------------------------------------------
            # UCSC BED-style:
            #   start = 0-based
            #   end   = exclusive
            #
            # JSON:
            #   1-based inclusive
            #
            # ex)
            #   UCSC: 0 ~ 2300000
            #   JSON: 1 ~ 2300000
            # -------------------------------------------------

            start_1based = start + 1
            end_1based = end

            cytobands = chromosomes[chrom].setdefault(
                "cytobands",
                {}
            )

            cytobands[band] = {
                "start": start_1based,
                "end": end_1based,
                "stain": stain,
            }

            band_count += 1

    # ---------------------------------------------------------
    # 3. 저장
    # ---------------------------------------------------------
    output_path = Path(output_json_path)

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(
            data,
            f,
            ensure_ascii=False,
            indent=2,
        )

    print(f"완료: {output_path}")
    print(f"cytoband {band_count}개 추가")


if __name__ == "__main__":

    add_cytobands_to_json(
        chromosome_json_path="/storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/nipt_pipeline/references/hg38.json",
        cytoband_txt_path="/storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Resources/reference/a.txt",
        output_json_path="hg38_with_cytobands.json",
    )