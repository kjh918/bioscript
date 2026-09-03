# karyotype_viewer

Interactive clinical karyotype report viewer built with **FastAPI + Dash-Bio**.

---

## Package layout

```
karyotype_viewer/
├── core/               # UI-독립 데이터 레이어
│   ├── reference.py    # GRCh38 염색체 크기 + Ideogram 설정 상수
│   ├── models.py       # SampleData / CnvEvent / GeneAnnotation 데이터클래스
│   ├── parser.py       # TSV/CSV 파일 로더 & 컬럼 alias 정규화
│   ├── annot.py        # Ideogram rawAnnots / ISCN / Dropdown options 빌더
│   └── data.py         # 데모 CN/BAF DataFrame 생성기
│
├── ui/                 # Plotly figure + CSS-in-Python 공유 컴포넌트
│   ├── styles.py       # PAGE / CARD / dashbar() / card() 헬퍼
│   └── figures.py      # region_fig() – CN/BAF scatter 서브플롯
│
├── apps/               # Dash 앱 팩토리 (각 앱 = 독립 WSGI 인스턴스)
│   ├── overview.py     # create_overview_app() – 전체 핵형 Ideogram
│   └── detail.py       # create_detail_app()   – 단일 염색체 + brush + CN/BAF
│
├── server.py           # create_app() – FastAPI 호스트 + 두 Dash 앱 마운트
├── demo_sample.py      # make_demo_sample() – 내장 예시 샘플 (DEMO_001)
├── main.py             # CLI 진입점 (argparse + uvicorn)
└── __main__.py         # python -m karyotype_viewer 지원
```

---

## 설치

```bash
pip install -e .
# 또는
pip install fastapi uvicorn dash dash-bio plotly pandas numpy
```

---

## 실행

### 1. 내장 데모 데이터로 미리보기

```bash
python -m karyotype_viewer
# → http://localhost:8050/
```

### 2. 샘플 TSV 파일 지정

```bash
python -m karyotype_viewer --sample data/example_sample.tsv
```

### 3. 초기 염색체 + 포트 지정

```bash
python -m karyotype_viewer --sample data/example_sample.tsv --chrom 17 --port 8080
```

### 4. `create_app()` 직접 호출 (스크립트)

```python
import uvicorn
from karyotype_viewer.demo_sample import make_demo_sample
from karyotype_viewer.server import create_app

sample = make_demo_sample()
app = create_app(sample, initial_chrom="21")

uvicorn.run(app, host="0.0.0.0", port=8050)
```

---

## 샘플 TSV 포맷

`data/example_sample.tsv` 참조. 탭 구분, 헤더 필수.

| 컬럼 | 설명 |
|------|------|
| sample_id | 샘플 ID |
| sex | female / male |
| chr | 염색체 번호 (1-22, X, Y) |
| type | trisomy / monosomy / partial_gain / partial_loss / **gene** |
| cn | Copy number (정수) |
| iscn | ISCN 표기 (CNV 행) |
| start | 시작 position (bp) |
| stop | 종료 position (bp) |
| color | 어노테이션 색상 (hex) |
| gene | 유전자명 (type=gene 행에서 사용) |

---

## CN/BAF 데이터 TSV 포맷 (UI 업로드용)

| 컬럼 | alias | 설명 |
|------|-------|------|
| chrom | chromosome, chr | 염색체 (선택) |
| pos | position, start | genomic position (bp) **필수** |
| cn | copy_number, ratio, log2 | copy number **필수** |
| baf | vaf, b_allele_freq | B-allele frequency (선택) |

```bash
# 예시 데이터 생성
python data/generate_example_cnv.py
```

---

## 아키텍처 요약

```
http://localhost:8050/          ← FastAPI (host shell HTML)
  └── /overview/ [iframe]       ← Dash App 1 (WSGIMiddleware)
        Ideogram (vertical, 전체 핵형)
        click → postMessage → host → detail iframe

  └── /detail/   [iframe]       ← Dash App 2 (WSGIMiddleware)
        Ideogram (horizontal, 단일 염색체, brush)
        dcc.Dropdown (annotation 선택)
        dcc.Graph (CN / BAF scatter)
        dcc.Upload (TSV/CSV 업로드)
```

두 Dash 앱은 **독립 WSGI 인스턴스**로 분리되어 Ideogram.js DOM을 공유하지 않습니다.
염색체 선택은 `postMessage` → `dash_clientside.set_props` 로 서버 라운드트립 없이 전달됩니다.
