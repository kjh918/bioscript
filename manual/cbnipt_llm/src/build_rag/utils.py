import os
import sys
from pathlib import Path
from typing import TextIO, Optional


class Tee:
    """출력을 여러 stream에 동시에 기록합니다."""

    def __init__(self, *streams: TextIO):
        self.streams = streams

    def write(self, message: str) -> int:
        for stream in self.streams:
            stream.write(message)
            stream.flush()
        return len(message)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()

    def isatty(self) -> bool:
        """일부 라이브러리가 stdout의 터미널 여부를 확인할 때 사용합니다."""
        return any(
            hasattr(stream, "isatty") and stream.isatty()
            for stream in self.streams
        )


class TeeLogger:
    """
    stdout과 stderr를 터미널 및 로그 파일에 동시에 연결합니다.

    close() 호출 시 stdout/stderr를 원래 상태로 복원합니다.
    """

    def __init__(
        self,
        log_path: str,
        mode: str = "a",
        encoding: str = "utf-8",
    ):
        self.log_path = Path(log_path)
        self.mode = mode
        self.encoding = encoding

        self.original_stdout: TextIO = sys.stdout
        self.original_stderr: TextIO = sys.stderr
        self.log_file: Optional[TextIO] = None

    def start(self) -> "TeeLogger":
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

        self.log_file = open(
            self.log_path,
            self.mode,
            encoding=self.encoding,
            buffering=1,
        )

        sys.stdout = Tee(self.original_stdout, self.log_file)
        sys.stderr = Tee(self.original_stderr, self.log_file)

        return self

    def close(self) -> None:
        sys.stdout = self.original_stdout
        sys.stderr = self.original_stderr

        if self.log_file is not None and not self.log_file.closed:
            self.log_file.flush()
            self.log_file.close()

    def __enter__(self) -> "TeeLogger":
        return self.start()

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


def get_default_log_path(output_path: str) -> str:
    """
    출력 파일 경로를 기준으로 로그 파일 경로를 생성합니다.

    예:
        medgen_results.json -> medgen_results.log
        result/out.json     -> result/out.log
    """
    output = Path(output_path)
    return str(output.with_suffix(".log"))


def setup_tee_logging(
    log_path: str,
    mode: str = "a",
) -> TeeLogger:
    """TeeLogger를 생성하고 즉시 시작합니다."""
    return TeeLogger(
        log_path=log_path,
        mode=mode,
    ).start()