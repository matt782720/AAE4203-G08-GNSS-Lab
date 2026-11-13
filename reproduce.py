"""Regenerate and validate both G08 submission CSVs from the lab inputs."""
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run(*arguments):
    subprocess.run([str(arg) for arg in arguments], cwd=ROOT, check=True)


def main():
    for name in ("data", "figures", "benchmark"):
        (ROOT / name).mkdir(exist_ok=True)

    rover_obs = ROOT / "data/G08.obs"
    rover_nav = ROOT / "data/G08.nav"
    run(ROOT / "tools/convbin.exe", "-r", "ubx", "-v", "3.04", "-od", "-os",
        "-oi", "-ot", "-f", "2", "-o", rover_obs, "-n", rover_nav,
        ROOT / "raw/G08_owndata_20261030.ubx")
    run(ROOT / "tools/rnx2rtkp.exe", "-k", ROOT / "code/spp.conf", "-p", "0",
        "-e", "-o", ROOT / "data/RTKLIB_GPS.pos", rover_obs, rover_nav)
    run(sys.executable, "code/process_g08.py")
    run(sys.executable, "course/WLSE_spp.py", "--mode", "gps", "--input-dir",
        "course/benchmark/gps_csv", "--output-dir", "benchmark",
        "--submission", "G08_benchmark.csv")
    run(sys.executable, "course/tools/check_submission.py", "G08_owndata.csv",
        "--week", "2442", "--epoch-file", "data/epochs.csv")
    run(sys.executable, "course/tools/check_submission.py", "benchmark/G08_benchmark.csv",
        "--week", "2158", "--epoch-file", "course/benchmark/gps_csv/epochs.csv")

    for line in (ROOT / "submission.sha256").read_text(encoding="ascii").splitlines():
        expected, filename = line.split(maxsplit=1)
        actual = hashlib.sha256((ROOT / filename).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError(f"Regenerated CSV differs from the submission copy: {filename}")
        print(f"Byte-for-byte match: {filename}")


if __name__ == "__main__":
    main()
