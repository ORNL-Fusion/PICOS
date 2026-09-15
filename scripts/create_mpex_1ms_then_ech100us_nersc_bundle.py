#!/usr/bin/env python3
"""Create a NERSC bundle for 1 ms kinetic steady state then 100 us ECH.

The steady stage runs PICOS_ECH with kinetic D+ ions and kinetic electrons,
collisions on, z=0 pair source on, and ECH off.  The ECH stage restarts from
the final steady HDF5 snapshot and applies electron ECH for another 100 us.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import stat
import tarfile
import textwrap
from dataclasses import dataclass
from pathlib import Path


RUN_ROOT_NAME = "PICOS_NERSC_MPEX_kinetic_1ms_then_ech100us"

STEADY_2MS_TAG = "mpex_scenario14_ex8_steady_2ms_p262144_coll_mpexprof_nersc_nonrel"
STEADY_1MS_TAG = "mpex_scenario14_ex8_steady_1ms_p262144_coll_mpexprof_nersc_nonrel"
SMOKE_STEADY_TAG = "mpex_scenario14_ex8_steady_interactive_1ms_p32768_coll_mpexprof_nersc_nonrel"

GYROPERIODS_1MS = "1.0430367059868449e+04"
GYROPERIODS_100US = "1.0430367059868452e+03"
OUTPUT_CADENCE_1MS_20 = "5.2151835299342247e+02"
MPEX_NE = "5.0000000000000000e+19"
MPEX_TEMP = "1.5000000000000000e+01"
MPEX_SOURCE_RATE = "6.0000000000000000e+22"
MPEX_SOURCE_CENTER = "0.0000000000000000e+00"
MPEX_SOURCE_SIGMA = "2.9999999999999999e-01"
MPEX_PROFILE_N = 200
MPEX_ZMIN = -2.0
MPEX_ZMAX = 8.0


@dataclass(frozen=True)
class Case:
    subdir: str
    tag: str
    job_name: str
    walltime: str
    restart_required: bool
    steady_subdir: str
    source_input_tag: str
    source_ion_tag: str


CASES = (
    Case(
        "steady_1ms",
        STEADY_1MS_TAG,
        "mpex_s1ms",
        "18:00:00",
        False,
        "",
        STEADY_2MS_TAG,
        STEADY_2MS_TAG,
    ),
    Case(
        "left_resonance_70ghz",
        "mpex_scenario14_ex8_left_resonance_70ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "ech_l70_1m",
        "11:30:00",
        True,
        "steady_1ms",
        "mpex_scenario14_ex8_left_resonance_70ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "mpex_scenario14_ex8_left_resonance_70ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
    ),
    Case(
        "right_resonance_70ghz",
        "mpex_scenario14_ex8_right_resonance_70ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "ech_r70_1m",
        "11:30:00",
        True,
        "steady_1ms",
        "mpex_scenario14_ex8_right_resonance_70ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "mpex_scenario14_ex8_right_resonance_70ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
    ),
    Case(
        "well_min_65p588ghz",
        "mpex_scenario14_ex8_well_min_65p588ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "ech_w65_1m",
        "11:30:00",
        True,
        "steady_1ms",
        "mpex_scenario14_ex8_well_min_65p588ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "mpex_scenario14_ex8_well_min_65p588ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
    ),
    Case(
        "steady_interactive_1ms",
        SMOKE_STEADY_TAG,
        "mpex_i1ms",
        "04:00:00",
        False,
        "",
        SMOKE_STEADY_TAG,
        SMOKE_STEADY_TAG,
    ),
    Case(
        "well_min_interactive_100us",
        "mpex_scenario14_ex8_well_min_65p588ghz_target8_p32768_100us_coll_mpexprof_nersc_nonrel",
        "ech_iw65",
        "04:00:00",
        True,
        "steady_interactive_1ms",
        "mpex_scenario14_ex8_well_min_65p588ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "mpex_scenario14_ex8_well_min_65p588ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
    ),
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def write_executable(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip())
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def replace_key(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^({re.escape(key)}\s+).*$", re.MULTILINE)
    replacement = rf"\g<1>{value}"
    updated, count = pattern.subn(replacement, text)
    if count == 0:
        raise RuntimeError(f"Could not replace key {key}")
    return updated


def append_or_replace_key(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^({re.escape(key)}\s+).*$", re.MULTILINE)
    lines = []
    count = 0
    for line in text.splitlines():
        match = pattern.match(line)
        if match:
            count += 1
            if count == 1:
                lines.append(f"{match.group(1)}{value}")
        else:
            lines.append(line)
    updated = "\n".join(lines) + ("\n" if text.endswith("\n") else "")
    if count == 0:
        updated = text.rstrip() + f"\n{key:<30} {value}\n"
    return updated


def normalize_mpex_input(text: str) -> str:
    replacements = {
        "quietStart": "1",
        "CV_ne": MPEX_NE,
        "CV_Te": MPEX_TEMP,
        "CV_Tpar": MPEX_TEMP,
        "CV_Tper": MPEX_TEMP,
        "IC_ne": MPEX_NE,
        "IC_Te": MPEX_TEMP,
        "IC_Tpar": MPEX_TEMP,
        "IC_Tper": MPEX_TEMP,
        "pairSource_rate": MPEX_SOURCE_RATE,
        "pairSource_mean_x": MPEX_SOURCE_CENTER,
        "pairSource_sigma_x": MPEX_SOURCE_SIGMA,
        "pairSource_positionMode": "0",
        "pairSource_weightMode": "1",
        "pairSource_Ti_birth": MPEX_TEMP,
        "pairSource_Te_birth": MPEX_TEMP,
        "filtersPerIterationFields": "3",
        "filtersPerIterationIons": "3",
    }
    for key, value in replacements.items():
        text = append_or_replace_key(text, key, value)
    return text


def normalize_mpex_ions(text: str) -> str:
    replacements = {
        "BC_type_1": "1",
        "BC_type_2": "1",
        "BC_T_1": MPEX_TEMP,
        "BC_T_2": MPEX_TEMP,
        "BC_mean_x_1": MPEX_SOURCE_CENTER,
        "BC_mean_x_2": MPEX_SOURCE_CENTER,
        "BC_sigma_x_1": MPEX_SOURCE_SIGMA,
        "BC_sigma_x_2": MPEX_SOURCE_SIGMA,
        "BC_G_1": MPEX_SOURCE_RATE,
        "BC_G_2": MPEX_SOURCE_RATE,
        "IC_weightScale_1": "1.0000000000000000e+00",
        "IC_weightScale_2": "1.0000000000000000e+00",
    }
    for key, value in replacements.items():
        text = append_or_replace_key(text, key, value)
    return text


def paper_profile_text(name: str) -> str | None:
    if name.endswith(("_ne_norm.txt", "_te_norm.txt", "_ti_norm.txt")):
        return "\n".join(["1.0000000000000000e+00"] * MPEX_PROFILE_N) + "\n"
    if name.endswith("_pair_source_norm.txt"):
        values = []
        for ii in range(MPEX_PROFILE_N):
            x = MPEX_ZMIN + (MPEX_ZMAX - MPEX_ZMIN) * ii / max(MPEX_PROFILE_N - 1, 1)
            values.append(math.exp(-0.5 * ((x - 0.0) / 0.3) ** 2))
        vmax = max(values) if values else 1.0
        return "\n".join(f"{value / vmax:.16e}" for value in values) + "\n"
    return None


def source_input_text(input_dir: Path, case: Case) -> str:
    path = input_dir / f"input_file_{case.source_input_tag}.input"
    text = path.read_text()
    text = replace_key(text, "Poisson_BCModel", "2")
    text = normalize_mpex_input(text)
    if case.tag == STEADY_1MS_TAG:
        text = text.replace("2 ms", "1 ms")
        text = replace_key(text, "simulationTime", GYROPERIODS_1MS)
        text = replace_key(text, "RF_ion_t_OFF", GYROPERIODS_1MS)
        text = replace_key(text, "RF_electron_t_OFF", GYROPERIODS_1MS)
        text = replace_key(text, "outputCadence", OUTPUT_CADENCE_1MS_20)
        text = text.replace("40 output intervals", "20 output intervals")
    if case.subdir == "well_min_interactive_100us":
        text = replace_key(text, "simulationTime", GYROPERIODS_100US)
        text = replace_key(text, "outputCadence", "2.0860734119736904e+01")
    return text


def source_ion_text(input_dir: Path, case: Case) -> str:
    path = input_dir / f"ions_properties_{case.source_ion_tag}.ion"
    text = normalize_mpex_ions(path.read_text())
    if case.subdir == "well_min_interactive_100us":
        text = replace_key(text, "NPC1", "41")
        text = replace_key(text, "NPC2", "165")
    return text


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


def populate_case(root: Path, input_dir: Path, case: Case) -> None:
    case_dir = root / case.subdir
    dst = case_dir / "inputFiles"
    dst.mkdir(parents=True, exist_ok=True)

    input_text = source_input_text(input_dir, case)
    ion_text = source_ion_text(input_dir, case)
    (dst / "input_file.input").write_text(input_text)
    (dst / "ions_properties.ion").write_text(ion_text)

    for name in sorted(referenced_input_files(input_text, ion_text)):
        src = input_dir / name
        if not src.is_file():
            raise FileNotFoundError(f"Referenced input file is missing: {src}")
        profile_text = paper_profile_text(name)
        if profile_text is None:
            shutil.copyfile(src, dst / name)
        else:
            (dst / name).write_text(profile_text)

    write_executable(case_dir / "run_case_common_nersc.sh", common_script())
    write_executable(case_dir / "run_case_nersc.sh", batch_script(case))
    write_executable(case_dir / "run_case_nersc_interactive.sh", interactive_script(case))


def common_script() -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    run_picos_case() {{
      : "${{CASE_TAG:?CASE_TAG is required}}"
      : "${{CASE_SUBDIR:?CASE_SUBDIR is required}}"
      : "${{CASE_DIR:?CASE_DIR is required}}"
      : "${{RESTART_REQUIRED:?RESTART_REQUIRED is required}}"
      STEADY_SUBDIR=${{STEADY_SUBDIR:-}}

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
        echo "Set PICOS_ROOT, PICOS_BUILD_DIR, or PICOS_BIN." >&2
        exit 2
      fi
      if [ $((MPI_RANKS % 2)) -ne 0 ]; then
        echo "MPI_RANKS must be even for PICOS++." >&2
        exit 2
      fi

      STAGE_ROOT="${{RUN_ROOT}}/${{CASE_SUBDIR}}/run"
      RUN_PICOS_FILES="${{STAGE_ROOT}}/picosFILES"
      mkdir -p "${{RUN_PICOS_FILES}}/inputFiles" "${{STAGE_ROOT}}/logs"
      rsync -a "${{CASE_DIR}}/inputFiles/" "${{RUN_PICOS_FILES}}/inputFiles/"
      git -C "${{PICOS_ROOT}}" log --oneline -1 > "${{STAGE_ROOT}}/commitHash.txt" || true

      input_file="${{RUN_PICOS_FILES}}/inputFiles/input_file.input"
      if [ ! -f "${{input_file}}" ] || [ ! -f "${{RUN_PICOS_FILES}}/inputFiles/ions_properties.ion" ]; then
        echo "Missing untagged input_file.input or ions_properties.ion under ${{RUN_PICOS_FILES}}/inputFiles." >&2
        exit 2
      fi

      if [ "${{RESTART_REQUIRED}}" = "1" ]; then
        : "${{STEADY_SUBDIR:?STEADY_SUBDIR is required for restart cases}}"
        STEADY_HDF5=${{STEADY_HDF5:-${{RUN_ROOT}}/${{STEADY_SUBDIR}}/run/picosFILES/outputFiles/HDF5}}
        if [ ! -f "${{STEADY_HDF5}}/PARTICLES_FILE_0.h5" ]; then
          echo "Missing steady-state restart files under: ${{STEADY_HDF5}}" >&2
          echo "Run the matching 1 ms steady case first or set STEADY_HDF5 explicitly." >&2
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
      echo "SLURM_JOB_ID=${{SLURM_JOB_ID:-interactive}}"
      echo "PICOS_ROOT=${{PICOS_ROOT}}"
      echo "PICOS_BIN=${{PICOS_BIN}}"
      echo "RUN_ROOT=${{RUN_ROOT}}"
      echo "STAGE_ROOT=${{STAGE_ROOT}}"
      echo "MPI_RANKS=${{MPI_RANKS}}"
      echo "OMP_NUM_THREADS=${{OMP_NUM_THREADS}}"
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


def batch_script(case: Case) -> str:
    restart = "1" if case.restart_required else "0"
    return f"""
    #!/bin/bash
    #SBATCH --account=m77
    #SBATCH -N 1
    #SBATCH -C cpu
    #SBATCH -q regular
    #SBATCH -J {case.job_name}
    #SBATCH -t {case.walltime}
    #SBATCH --ntasks-per-node=128
    #SBATCH --cpus-per-task=1
    #SBATCH -o slurm-%x-%j.out
    #SBATCH -e slurm-%x-%j.err

    set -euo pipefail

    CASE_TAG={case.tag}
    CASE_SUBDIR={case.subdir}
    STEADY_SUBDIR={case.steady_subdir}
    RESTART_REQUIRED={restart}
    CASE_DIR=${{SLURM_SUBMIT_DIR:-$(pwd)}}
    if [ ! -f "${{CASE_DIR}}/inputFiles/input_file.input" ]; then
      echo "Submit from this case directory:" >&2
      echo "  cd <bundle>/{case.subdir} && sbatch ./run_case_nersc.sh" >&2
      exit 2
    fi

    source "${{CASE_DIR}}/run_case_common_nersc.sh"
    run_picos_case
    """


def interactive_script(case: Case) -> str:
    restart = "1" if case.restart_required else "0"
    return f"""
    #!/bin/bash
    set -euo pipefail

    if [ -z "${{SLURM_JOB_ID:-}}" ]; then
      echo "Start an interactive allocation first, for example:" >&2
      echo "  salloc -A m77 -C cpu -q interactive -N 1 -t 04:00:00 --ntasks-per-node=128" >&2
      exit 2
    fi

    CASE_TAG={case.tag}
    CASE_SUBDIR={case.subdir}
    STEADY_SUBDIR={case.steady_subdir}
    RESTART_REQUIRED={restart}
    CASE_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
    source "${{CASE_DIR}}/run_case_common_nersc.sh"
    run_picos_case
    """


def submit_all_script() -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    ROOT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
    export RUN_ROOT=${{RUN_ROOT:-${{SCRATCH}}/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}}}
    echo "RUN_ROOT=${{RUN_ROOT}}"

    steady_job=$(cd "${{ROOT_DIR}}/steady_1ms" && sbatch --parsable ./run_case_nersc.sh)
    echo "Submitted 1 ms kinetic steady-state job: ${{steady_job}}"
    echo "Submitting 100 us ECH restart jobs with dependency afterok:${{steady_job}}"

    for case_dir in left_resonance_70ghz right_resonance_70ghz well_min_65p588ghz; do
      job_id=$(cd "${{ROOT_DIR}}/${{case_dir}}" && sbatch --parsable --dependency=afterok:${{steady_job}} ./run_case_nersc.sh)
      echo "Submitted ${{case_dir}}: ${{job_id}}"
    done

    echo "Restart HDF5 path:"
    echo "  ${{RUN_ROOT}}/steady_1ms/run/picosFILES/outputFiles/HDF5"
    echo "Monitor with:"
    echo "  squeue -u $USER"
    echo "  sacct -j ${{steady_job}} --format=JobID,JobName,State,Elapsed,ExitCode"
    """


def submit_ech_only_script() -> str:
    return """
    #!/bin/bash
    set -euo pipefail

    : "${STEADY_HDF5:?Set STEADY_HDF5=/path/to/1ms/steady/outputFiles/HDF5 first.}"
    ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    export RUN_ROOT=${RUN_ROOT:-${SCRATCH}/picosRuns/MPEX_ECH_runs/PICOS_NERSC_MPEX_kinetic_1ms_then_ech100us}
    export STEADY_HDF5

    for case_dir in left_resonance_70ghz right_resonance_70ghz well_min_65p588ghz; do
      job_id=$(cd "${ROOT_DIR}/${case_dir}" && sbatch --parsable ./run_case_nersc.sh)
      echo "Submitted ${case_dir}: ${job_id}"
    done
    """


def run_interactive_smoke_script() -> str:
    return """
    #!/bin/bash
    set -euo pipefail

    if [ -z "${SLURM_JOB_ID:-}" ]; then
      echo "Start an interactive allocation first, for example:" >&2
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
    # MPEX Scenario 14: 1 ms Kinetic Steady State Then 100 us ECH

    This bundle uses only the PICOS_ECH branch.  The restart source must be a
    kinetic run with both `spp_1` D+ ions and `spp_2` electrons.  Do not restart
    these ECH decks from the main-vs-ECH hybrid compatibility archive; that
    archive contains only `spp_1`.

    Production chain:

    - `steady_1ms`: kinetic D+ plus kinetic electrons, collisions on, z=0 pair
      source on, flat `n/T` initial profiles, `n_e=n_i=5e19 m^-3`,
      `T_e=T_i=15 eV`, `G=6e22 s^-1`, reformulated Poisson field solve,
      ECH off, 1 ms.
    - `left_resonance_70ghz`: restart from `steady_1ms`, electron ECH on,
      ion RF off, 100 us.
    - `right_resonance_70ghz`: same, right-side resonance, 100 us.
    - `well_min_65p588ghz`: same, well-minimum case, 100 us.

    Interactive smoke chain:

    - `steady_interactive_1ms`: lower-particle 1 ms kinetic steady state.
    - `well_min_interactive_100us`: 100 us ECH restart from that smoke steady
      state.

    The restart scripts set:

    - `restart_enabled = 1`
    - `restart_snapshot = -1`
    - `restart_continueTime = 1`
    - `restart_path = $STEADY_HDF5` or the matching local steady HDF5 directory

    ## Production on NERSC

    ```bash
    cd /pscratch/sd/a/atul19/picosRuns/MPEX_ECH_runs
    tar -xzf {RUN_ROOT_NAME}.tar.gz
    cd {RUN_ROOT_NAME}

    export PICOS_ROOT=$HOME/myRepos/PICOS_ECH
    export PICOS_BUILD_DIR=$PICOS_ROOT/build
    export RUN_ROOT=$SCRATCH/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}
    ./submit_1ms_then_ech100us_nersc.sh
    ```

    ## Interactive smoke test

    ```bash
    cd /pscratch/sd/a/atul19/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}
    salloc -A m77 -C cpu -q interactive -N 1 -t 04:00:00 --ntasks-per-node=128

    export PICOS_ROOT=$HOME/myRepos/PICOS_ECH
    export PICOS_BUILD_DIR=$PICOS_ROOT/build
    export RUN_ROOT=$SCRATCH/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}_interactive
    export OVERWRITE_OUTPUTS=1
    ./run_interactive_1ms_then_wellmin_ech100us.sh
    exit
    ```

    ## ECH only after an existing 1 ms kinetic steady state

    ```bash
    export PICOS_ROOT=$HOME/myRepos/PICOS_ECH
    export PICOS_BUILD_DIR=$PICOS_ROOT/build
    export RUN_ROOT=$SCRATCH/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}_ech_only
    export STEADY_HDF5=/path/to/steady_1ms/run/picosFILES/outputFiles/HDF5
    ./submit_ech100us_after_existing_1ms_nersc.sh
    ```

    Monitor:

    ```bash
    squeue -u $USER
    sacct -j <jobid> --format=JobID,JobName,State,Elapsed,ExitCode
    find $RUN_ROOT -maxdepth 8 -type f -name '*.h5' -print | head
    tail -f $RUN_ROOT/steady_1ms/run/logs/{STEADY_1MS_TAG}.log
    ```

    Keep the same `MPI_RANKS` for the steady and ECH restart stages, because
    each particle MPI rank loads its corresponding `PARTICLES_FILE_<rank>.h5`.
    """


def validate_bundle(root: Path) -> None:
    checks = {
        root / "steady_1ms/inputFiles/input_file.input": [
            "numberOfParticleSpecies     2",
            f"simulationTime              {GYROPERIODS_1MS}",
            "Poisson_BCModel             2",
            "quietStart                  1",
            f"CV_ne                       {MPEX_NE}",
            f"IC_ne                       {MPEX_NE}",
            f"pairSource_rate             {MPEX_SOURCE_RATE}",
            f"pairSource_sigma_x          {MPEX_SOURCE_SIGMA}",
            "SW_RFheating                0",
            "SW_pairSource               1",
            "restart_enabled             0",
        ],
        root / "well_min_65p588ghz/inputFiles/input_file.input": [
            "numberOfParticleSpecies     2",
            f"simulationTime              {GYROPERIODS_100US}",
            "Poisson_BCModel             2",
            "SW_RFheating                1",
            "SW_RFheatingIons            0",
            "SW_RFheatingElectrons       1",
            "SW_Collisions               1",
            f"CV_ne                       {MPEX_NE}",
            f"IC_ne                       {MPEX_NE}",
            f"pairSource_rate             {MPEX_SOURCE_RATE}",
        ],
        root / "steady_interactive_1ms/inputFiles/ions_properties.ion": [
            "Z2                            -1",
            "BC_mean_x_2                   0.0000000000000000e+00",
            f"BC_sigma_x_2                  {MPEX_SOURCE_SIGMA}",
            f"BC_G_1                        {MPEX_SOURCE_RATE}",
            f"BC_G_2                        {MPEX_SOURCE_RATE}",
        ],
    }
    for path, required in checks.items():
        text = path.read_text()
        missing = [item for item in required if item not in text]
        if missing:
            raise RuntimeError(f"{path} is missing required settings: {missing}")


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

    for case in CASES:
        populate_case(root, input_dir, case)

    write_executable(root / "submit_1ms_then_ech100us_nersc.sh", submit_all_script())
    write_executable(root / "submit_ech100us_after_existing_1ms_nersc.sh", submit_ech_only_script())
    write_executable(root / "run_interactive_1ms_then_wellmin_ech100us.sh", run_interactive_smoke_script())
    (root / "README.md").write_text(textwrap.dedent(readme_text()).lstrip())

    validate_bundle(root)
    create_archive(root, archive)

    print(root)
    print(archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
