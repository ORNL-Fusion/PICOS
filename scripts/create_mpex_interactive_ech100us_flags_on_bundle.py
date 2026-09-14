#!/usr/bin/env python3
"""Create an interactive-only NERSC bundle with restart/ECH flags visible."""

from __future__ import annotations

import argparse
import re
import shutil
import stat
import tarfile
import textwrap
from pathlib import Path


RUN_ROOT_NAME = "PICOS_NERSC_MPEX_interactive_1ms_then_wellmin_ech100us_flags_on"
STEADY_TAG = "mpex_scenario14_ex8_steady_interactive_1ms_p32768_coll_mpexprof_nersc_nonrel"
ECH_SOURCE_TAG = "mpex_scenario14_ex8_well_min_65p588ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel"
ECH_TAG = "mpex_scenario14_ex8_well_min_65p588ghz_target8_p32768_100us_coll_mpexprof_nersc_nonrel"

DEFAULT_STEADY_HDF5 = (
    f"/pscratch/sd/a/atul19/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}/"
    "steady_interactive_1ms/run/picosFILES/outputFiles/HDF5"
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def write_executable(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip())
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def replace_key(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^({re.escape(key)}\s+).*$", re.MULTILINE)
    updated, count = pattern.subn(rf"\g<1>{value}", text)
    if count == 0:
        raise RuntimeError(f"Could not replace key {key}")
    return updated


def referenced_input_files(*texts: str) -> set[str]:
    refs: set[str] = set()
    pattern = re.compile(r"^\s*[A-Za-z0-9_]*fileName[A-Za-z0-9_]*\s+(\S+)\s*$")
    for text in texts:
        for line in text.splitlines():
            match = pattern.match(line)
            if match:
                name = match.group(1)
                if name.lower() != "none":
                    refs.add(name)
    return refs


def read_input(input_dir: Path, tag: str) -> str:
    return (input_dir / f"input_file_{tag}.input").read_text()


def read_ions(input_dir: Path, tag: str) -> str:
    return (input_dir / f"ions_properties_{tag}.ion").read_text()


def steady_input(input_dir: Path) -> str:
    text = read_input(input_dir, STEADY_TAG)
    return replace_key(text, "Poisson_BCModel", "2")


def steady_ions(input_dir: Path) -> str:
    return read_ions(input_dir, STEADY_TAG)


def ech_input(input_dir: Path) -> str:
    text = read_input(input_dir, ECH_SOURCE_TAG)
    text = replace_key(text, "Poisson_BCModel", "2")
    text = replace_key(text, "restart_enabled", "1")
    text = replace_key(text, "restart_path", DEFAULT_STEADY_HDF5)
    text = replace_key(text, "restart_snapshot", "-1")
    text = replace_key(text, "restart_continueTime", "1")
    return text


def ech_ions(input_dir: Path) -> str:
    text = read_ions(input_dir, ECH_SOURCE_TAG)
    text = replace_key(text, "NPC1", "41")
    text = replace_key(text, "NPC2", "165")
    return text


def write_case(case_dir: Path, input_text: str, ion_text: str, source_input_dir: Path) -> None:
    input_dir = case_dir / "inputFiles"
    input_dir.mkdir(parents=True, exist_ok=True)
    (input_dir / "input_file.input").write_text(input_text)
    (input_dir / "ions_properties.ion").write_text(ion_text)
    for name in sorted(referenced_input_files(input_text, ion_text)):
        src = source_input_dir / name
        if not src.is_file():
            raise FileNotFoundError(src)
        shutil.copy2(src, input_dir / name)
    write_executable(case_dir / "run_case_common_nersc.sh", common_script())
    write_executable(case_dir / "run_case_nersc_interactive.sh", interactive_case_script(case_dir.name))


def common_script() -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    run_picos_case() {{
      : "${{CASE_TAG:?CASE_TAG is required}}"
      : "${{CASE_SUBDIR:?CASE_SUBDIR is required}}"
      : "${{RESTART_REQUIRED:?RESTART_REQUIRED is required}}"

      PICOS_LOAD_MODULES=${{PICOS_LOAD_MODULES:-1}}
      PICOS_GCC_MODULE=${{PICOS_GCC_MODULE:-gcc-native/14}}
      if [ "${{PICOS_LOAD_MODULES}}" = "1" ] && command -v module >/dev/null 2>&1; then
        module load cpu
        module load PrgEnv-gnu
        module load "${{PICOS_GCC_MODULE}}"
        module load cray-hdf5
      fi

      if command -v g++ >/dev/null 2>&1; then
        LIBSTDCXX_DIR=$(dirname "$(g++ -print-file-name=libstdc++.so.6)")
        if [ -f "${{LIBSTDCXX_DIR}}/libstdc++.so.6" ]; then
          export LD_LIBRARY_PATH="${{LIBSTDCXX_DIR}}:${{LD_LIBRARY_PATH:-}}"
        fi
      fi

      PICOS_ROOT=${{PICOS_ROOT:-${{HOME}}/myRepos/PICOS_ECH}}
      PICOS_BUILD_DIR=${{PICOS_BUILD_DIR:-${{PICOS_ROOT}}/build}}
      PICOS_BIN=${{PICOS_BIN:-${{PICOS_BUILD_DIR}}/picosFILES/src/xpicos}}
      RUN_ROOT=${{RUN_ROOT:-${{SCRATCH}}/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}}}
      MPI_RANKS=${{MPI_RANKS:-128}}
      CPUS_PER_TASK=${{CPUS_PER_TASK:-1}}

      export OMP_NUM_THREADS=${{OMP_NUM_THREADS:-1}}
      export OMP_PLACES=${{OMP_PLACES:-threads}}
      export OMP_PROC_BIND=${{OMP_PROC_BIND:-spread}}
      export HDF5_USE_FILE_LOCKING=FALSE

      if [ ! -x "${{PICOS_BIN}}" ]; then
        echo "Missing PICOS++ executable: ${{PICOS_BIN}}" >&2
        exit 2
      fi
      if [ $((MPI_RANKS % 2)) -ne 0 ]; then
        echo "MPI_RANKS must be even for PICOS++." >&2
        exit 2
      fi

      CASE_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
      STAGE_ROOT="${{RUN_ROOT}}/${{CASE_SUBDIR}}/run"
      RUN_PICOS_FILES="${{STAGE_ROOT}}/picosFILES"
      mkdir -p "${{RUN_PICOS_FILES}}/inputFiles" "${{STAGE_ROOT}}/logs"
      rsync -a "${{CASE_DIR}}/inputFiles/" "${{RUN_PICOS_FILES}}/inputFiles/"
      git -C "${{PICOS_ROOT}}" log --oneline -1 > "${{STAGE_ROOT}}/commitHash.txt" || true

      input_file="${{RUN_PICOS_FILES}}/inputFiles/input_file.input"
      if [ "${{RESTART_REQUIRED}}" = "1" ]; then
        STEADY_HDF5=${{STEADY_HDF5:-${{RUN_ROOT}}/steady_interactive_1ms/run/picosFILES/outputFiles/HDF5}}
        if [ ! -f "${{STEADY_HDF5}}/PARTICLES_FILE_0.h5" ]; then
          echo "Missing steady-state restart files under: ${{STEADY_HDF5}}" >&2
          exit 2
        fi
        tmp_file="${{input_file}}.tmp"
        awk -v restart_path="${{STEADY_HDF5}}" '
          /^restart_enabled[[:space:]]/ {{ print "restart_enabled             1"; next }}
          /^restart_path[[:space:]]/ {{ print "restart_path                " restart_path; next }}
          /^restart_snapshot[[:space:]]/ {{ print "restart_snapshot            -1"; next }}
          /^restart_continueTime[[:space:]]/ {{ print "restart_continueTime        1"; next }}
          {{ print }}
        ' "${{input_file}}" > "${{tmp_file}}"
        mv "${{tmp_file}}" "${{input_file}}"
      fi

      if [ -d "${{RUN_PICOS_FILES}}/outputFiles/HDF5" ]; then
        if [ "${{OVERWRITE_OUTPUTS:-0}}" = "1" ]; then
          rm -rf "${{RUN_PICOS_FILES}}/outputFiles"
        else
          echo "Output exists: ${{RUN_PICOS_FILES}}/outputFiles/HDF5" >&2
          echo "Set OVERWRITE_OUTPUTS=1 to replace it." >&2
          exit 2
        fi
      fi
      ABS_OUTPUT_DIR="${{RUN_PICOS_FILES}}/outputFiles"

      echo "[$(date)] Starting ${{CASE_TAG}}"
      echo "RUN_ROOT=${{RUN_ROOT}}"
      echo "PICOS_BIN=${{PICOS_BIN}}"
      echo "MPI_RANKS=${{MPI_RANKS}}"
      echo "ABS_OUTPUT_DIR=${{ABS_OUTPUT_DIR}}"
      if [ "${{RESTART_REQUIRED}}" = "1" ]; then
        echo "STEADY_HDF5=${{STEADY_HDF5}}"
      fi

      (
        cd "${{RUN_PICOS_FILES}}"
        srun -n "${{MPI_RANKS}}" -c "${{CPUS_PER_TASK}}" --cpu-bind=cores \\
          "${{PICOS_BIN}}" 1-D "${{ABS_OUTPUT_DIR}}"
      ) > "${{STAGE_ROOT}}/logs/${{CASE_TAG}}.log" 2>&1

      echo "[$(date)] Finished ${{CASE_TAG}}"
      echo "Output: ${{RUN_PICOS_FILES}}/outputFiles"
    }}
    """


def interactive_case_script(subdir: str) -> str:
    if subdir == "steady_interactive_1ms":
        case_tag = STEADY_TAG
        restart = "0"
    else:
        case_tag = ECH_TAG
        restart = "1"
    return f"""
    #!/bin/bash
    set -euo pipefail

    if [ -z "${{SLURM_JOB_ID:-}}" ]; then
      echo "Start an interactive allocation first:" >&2
      echo "  salloc -A m77 -C cpu -q interactive -N 1 -t 04:00:00 --ntasks-per-node=128" >&2
      exit 2
    fi

    CASE_TAG={case_tag}
    CASE_SUBDIR={subdir}
    RESTART_REQUIRED={restart}
    source "$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)/run_case_common_nersc.sh"
    run_picos_case
    """


def top_level_run_script() -> str:
    return """
    #!/bin/bash
    set -euo pipefail

    if [ -z "${SLURM_JOB_ID:-}" ]; then
      echo "Start an interactive allocation first:" >&2
      echo "  salloc -A m77 -C cpu -q interactive -N 1 -t 04:00:00 --ntasks-per-node=128" >&2
      exit 2
    fi

    ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    cd "${ROOT_DIR}/steady_interactive_1ms"
    ./run_case_nersc_interactive.sh
    cd "${ROOT_DIR}/well_min_interactive_100us"
    ./run_case_nersc_interactive.sh
    """


def readme_text() -> str:
    return f"""
    # Interactive MPEX 1 ms Steady Then 100 us ECH Restart

    This bundle is interactive-only and contains only:

    - `steady_interactive_1ms`: kinetic D+ plus kinetic electrons, collisions
      and pair source on, ECH off, 1 ms.
    - `well_min_interactive_100us`: restart from the steady case, electron ECH
      on, ion RF off, 100 us.

    Unlike the larger production bundle, the ECH input deck already has the
    visible restart flags on:

    - `restart_enabled = 1`
    - `restart_snapshot = -1`
    - `restart_continueTime = 1`
    - `SW_RFheating = 1`
    - `SW_RFheatingElectrons = 1`
    - `SW_RFheatingIons = 0`

    The run wrapper still patches `restart_path` at runtime so it matches the
    actual `RUN_ROOT`.

    ## Run on NERSC

    ```bash
    cd /pscratch/sd/a/atul19/picosRuns/MPEX_ECH_runs
    tar -xzf {RUN_ROOT_NAME}.tar.gz
    cd {RUN_ROOT_NAME}

    salloc -A m77 -C cpu -q interactive -N 1 -t 04:00:00 --ntasks-per-node=128

    export PICOS_ROOT=$HOME/myRepos/PICOS_ECH
    export PICOS_BUILD_DIR=$PICOS_ROOT/build
    export RUN_ROOT=$SCRATCH/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}
    export OVERWRITE_OUTPUTS=1
    export MPI_RANKS=128
    export OMP_NUM_THREADS=1

    ./run_interactive_1ms_then_wellmin_ech100us.sh
    ```

    Check:

    ```bash
    find $RUN_ROOT -maxdepth 8 -type f -name '*.h5' -print | head
    tail -80 $RUN_ROOT/well_min_interactive_100us/run/logs/*.log
    ```
    """


def validate(root: Path) -> None:
    steady = (root / "steady_interactive_1ms/inputFiles/input_file.input").read_text()
    ech = (root / "well_min_interactive_100us/inputFiles/input_file.input").read_text()
    ions = (root / "well_min_interactive_100us/inputFiles/ions_properties.ion").read_text()
    required = {
        "steady": ["Poisson_BCModel             2", "SW_RFheating                0", "restart_enabled             0"],
        "ech": [
            "Poisson_BCModel             2",
            "SW_RFheating                1",
            "SW_RFheatingElectrons       1",
            "SW_RFheatingIons            0",
            "restart_enabled             1",
            "restart_continueTime        1",
            DEFAULT_STEADY_HDF5,
        ],
        "ions": ["NPC1                          41", "NPC2                          165", "Z2                            -1"],
    }
    for label, text in (("steady", steady), ("ech", ech), ("ions", ions)):
        missing = [item for item in required[label] if item not in text]
        if missing:
            raise RuntimeError(f"{label} missing {missing}")


def create_archive(root: Path, archive: Path) -> None:
    if archive.exists():
        archive.unlink()
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(root, arcname=root.name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, default=repo_root() / "validation" / RUN_ROOT_NAME)
    parser.add_argument("--archive", type=Path, default=None)
    args = parser.parse_args()

    root = args.out_root.resolve()
    archive = (args.archive if args.archive else root.with_suffix(".tar.gz")).resolve()
    input_dir = repo_root() / "picosFILES" / "inputFiles"

    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)

    write_case(root / "steady_interactive_1ms", steady_input(input_dir), steady_ions(input_dir), input_dir)
    write_case(root / "well_min_interactive_100us", ech_input(input_dir), ech_ions(input_dir), input_dir)
    write_executable(root / "run_interactive_1ms_then_wellmin_ech100us.sh", top_level_run_script())
    (root / "README.md").write_text(textwrap.dedent(readme_text()).lstrip())

    validate(root)
    create_archive(root, archive)
    print(root)
    print(archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
