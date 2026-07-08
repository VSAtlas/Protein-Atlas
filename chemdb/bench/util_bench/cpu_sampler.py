from __future__ import annotations

import csv
import os
import shutil
import subprocess
import threading
import time
from pathlib import Path
from typing import Any


class CPUSampler:
    backend: str

    def start(self) -> None:
        raise NotImplementedError

    def stop(self) -> None:
        raise NotImplementedError

    def join(self, timeout: float | None = None) -> None:
        raise NotImplementedError


class ProcStatSampler(CPUSampler):
    def __init__(self, out_csv: Path, interval_sec: float, alloc_cpus: int) -> None:
        self.backend = "procstat"
        self._out_csv = out_csv
        self._interval_sec = max(0.05, float(interval_sec))
        self._alloc_cpus = max(1, int(alloc_cpus))
        self._stop_evt = threading.Event()
        self._thread = threading.Thread(target=self._run, name="procstat-sampler", daemon=True)

    def start(self) -> None:
        self._thread.start()

    def stop(self) -> None:
        self._stop_evt.set()

    def join(self, timeout: float | None = None) -> None:
        self._thread.join(timeout)

    @staticmethod
    def _read_proc_stat() -> dict[int, list[int]]:
        values: dict[int, list[int]] = {}
        with Path("/proc/stat").open("r", encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith("cpu"):
                    continue
                parts = line.split()
                name = parts[0]
                if name == "cpu":
                    continue
                if not name[3:].isdigit():
                    continue
                cpu_idx = int(name[3:])
                nums = [int(x) for x in parts[1:11]]
                values[cpu_idx] = nums
        return values

    def _run(self) -> None:
        self._out_csv.parent.mkdir(parents=True, exist_ok=True)
        with self._out_csv.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["t_wall_ns", "cpu_idx", "user", "system", "idle", "iowait", "utilization"])

            prev = self._read_proc_stat()
            while not self._stop_evt.wait(self._interval_sec):
                now = self._read_proc_stat()
                now_ts = time.time_ns()

                cpus = sorted(set(prev.keys()) & set(now.keys()))
                if not cpus:
                    prev = now
                    continue

                max_cpus = min(self._alloc_cpus, len(cpus))
                for cpu_idx in cpus[:max_cpus]:
                    a = prev[cpu_idx]
                    b = now[cpu_idx]
                    d = [max(0, b_i - a_i) for a_i, b_i in zip(a, b)]

                    user = float(d[0] + d[1])
                    system = float(d[2] + d[5] + d[6])
                    idle = float(d[3])
                    iowait = float(d[4])
                    total = user + system + idle + iowait + float(d[7])
                    if total <= 0:
                        continue

                    user_f = user / total
                    system_f = system / total
                    idle_f = idle / total
                    iowait_f = iowait / total
                    util_f = max(0.0, min(1.0, 1.0 - idle_f))

                    writer.writerow([
                        now_ts,
                        cpu_idx,
                        f"{user_f:.8f}",
                        f"{system_f:.8f}",
                        f"{idle_f:.8f}",
                        f"{iowait_f:.8f}",
                        f"{util_f:.8f}",
                    ])

                fh.flush()
                prev = now


class MpstatSampler(CPUSampler):
    def __init__(self, raw_path: Path, interval_sec: float) -> None:
        self.backend = "mpstat"
        self._raw_path = raw_path
        self._interval_sec = max(1, int(round(float(interval_sec))))
        self._proc: subprocess.Popen[str] | None = None
        self._fh: Any = None

    def start(self) -> None:
        self._raw_path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = self._raw_path.open("w", encoding="utf-8", buffering=1)
        cmd = ["mpstat", "-P", "ALL", str(self._interval_sec)]
        env = dict(os.environ)
        env.setdefault("LC_ALL", "C")
        self._proc = subprocess.Popen(
            cmd,
            stdout=self._fh,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
        )

    def stop(self) -> None:
        if self._proc is None:
            return
        if self._proc.poll() is None:
            self._proc.terminate()
            try:
                self._proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self._proc.kill()
                self._proc.wait(timeout=5)

    def join(self, timeout: float | None = None) -> None:
        if self._proc is not None and self._proc.poll() is None:
            self._proc.wait(timeout=timeout)
        if self._fh is not None:
            self._fh.flush()
            self._fh.close()
            self._fh = None


def start_sampler(
    *,
    run_dir: str | Path,
    interval_sec: float,
    alloc_cpus: int,
    force_procstat: bool = False,
) -> CPUSampler:
    run_path = Path(run_dir)
    if (not force_procstat) and shutil.which("mpstat"):
        sampler: CPUSampler = MpstatSampler(
            run_path / "cpu_mpstat.log", interval_sec=interval_sec
        )
        sampler.start()
        return sampler

    sampler = ProcStatSampler(
        run_path / "cpu.log", interval_sec=interval_sec, alloc_cpus=alloc_cpus
    )
    sampler.start()
    return sampler
