import subprocess
from pathlib import Path

# 검색할 루트 디렉터리
ROOT_DIR = "/storage/home/jhkim/scripts/bioscript/manual/cbnipt_manual_cnv_call/src/cbnipt_cnv_caller/karyotype_viewer"

# 찾고 싶은 문자열
KEYWORDS = [
    "overall"]


def find_files(root_dir):
    """
    find 명령으로 검색 대상 파일 목록 수집
    """
    result = subprocess.run(
        [
            "find",
            root_dir,
            "-type", "f",
            "(",
            "-name", "*.py",
            "-o", "-name", "*.js",
            "-o", "-name", "*.html",
            "-o", "-name", "*.json",
            ")",
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    return [Path(p) for p in result.stdout.splitlines() if p.strip()]


def search_file(file_path, keywords):
    """
    파일 안에서 keyword가 포함된 줄을 찾음.
    """
    matches = []

    try:
        with file_path.open("r", encoding="utf-8", errors="replace") as f:
            for line_no, line in enumerate(f, start=1):
                for keyword in keywords:
                    if keyword in line:
                        matches.append(
                            {
                                "keyword": keyword,
                                "line_no": line_no,
                                "line": line.rstrip(),
                            }
                        )
    except Exception as e:
        print(f"[ERROR] {file_path}: {e}")

    return matches


def main():
    files = find_files(ROOT_DIR)

    print(f"검색 파일 수: {len(files)}")
    print("=" * 100)

    found_count = 0

    for file_path in files:
        matches = search_file(file_path, KEYWORDS)

        if not matches:
            continue

        found_count += 1

        print(f"\nFILE: {file_path}")

        for m in matches:
            print(
                f"  L{m['line_no']:>5} "
                f"[{m['keyword']}] "
                f"{m['line']}"
            )

    print("\n" + "=" * 100)
    print(f"keyword 발견 파일 수: {found_count}")


if __name__ == "__main__":
    main()