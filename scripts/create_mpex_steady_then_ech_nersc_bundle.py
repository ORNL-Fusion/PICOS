#!/usr/bin/env python3
"""Create the NERSC MPEX steady-state then ECH restart bundle.

The bundle contains one 2 ms RF-off steady-state case and three dependent
100 us ECH restart cases: left 70 GHz resonance, right 70 GHz resonance, and
the magnetic-well-minimum resonance.  The generated Slurm wrappers stage only
the input files required for each case.
"""

from __future__ import annotations

import argparse
import os
import shutil
import stat
import tarfile
import textwrap
from dataclasses import dataclass
from pathlib import Path


STEADY_TAG = "mpex_scenario14_ex8_steady_2ms_p262144_coll_mpexprof_nersc_nonrel"
SMOKE_TAG = "mpex_scenario14_ex8_steady_interactive_1ms_p32768_coll_mpexprof_nersc_nonrel"
RUN_ROOT_NAME = "PICOS_NERSC_MPEX_steady_then_ech_restart2ms_seeded_z0source"


@dataclass(frozen=True)
class BundleCase:
    subdir: str
    tag: str
    job_name: str
    walltime: str
    restart_required: bool


SMOKE_CASE = BundleCase("steady_interactive_1ms", SMOKE_TAG, "mpex_i1ms", "04:00:00", False)


PRODUCTION_CASES = (
    BundleCase("steady_2ms", STEADY_TAG, "mpex_steady2", "24:00:00", False),
    BundleCase(
        "left_resonance_70ghz",
        "mpex_scenario14_ex8_left_resonance_70ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "ech_l70_r2",
        "11:30:00",
        True,
    ),
    BundleCase(
        "right_resonance_70ghz",
        "mpex_scenario14_ex8_right_resonance_70ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "ech_r70_r2",
        "11:30:00",
        True,
    ),
    BundleCase(
        "well_min_65p588ghz",
        "mpex_scenario14_ex8_well_min_65p588ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel",
        "ech_w65_r2",
        "11:30:00",
        True,
    ),
)


ALL_CASES = (SMOKE_CASE, *PRODUCTION_CASES)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def case_input_files(tag: str) -> list[str]:
    return [
        f"input_file_{tag}.input",
        f"ions_properties_{tag}.ion",
        f"{tag}_B_norm.txt",
        f"{tag}_ne_norm.txt",
        f"{tag}_te_norm.txt",
        f"{tag}_ti_norm.txt",
        f"{tag}_pair_source_norm.txt",
        f"{tag}_rf_power_norm.txt",
    ]


def write_executable(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip())
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def copy_case_inputs(root: Path, case: BundleCase, source_dir: Path) -> None:
    case_dir = root / case.subdir
    input_dir = case_dir / "inputFiles"
    input_dir.mkdir(parents=True, exist_ok=True)
    for name in case_input_files(case.tag):
        src = source_dir / name
        if not src.is_file():
            raise FileNotFoundError(f"Missing generated input file: {src}")
        shutil.copy2(src, input_dir / name)


def assert_contains(path: Path, required: list[str]) -> None:
    text = path.read_text()
    missing = [item for item in required if item not in text]
    if missing:
        raise RuntimeError(f"{path} is missing required settings: {missing}")


def validate_case_inputs(root: Path) -> None:
    steady_input = root / "steady_2ms" / "inputFiles" / f"input_file_{STEADY_TAG}.input"
    assert_contains(
        steady_input,
        [
            "SW_EfieldSolve              1",
            "SW_fieldSolveModel          2",
            "SW_Collisions               1",
            "SW_RFheating                0",
            "SW_pairSource               1",
            "IC_weightScale              1.0000000000000000e+00",
            "CV_ne                       1.0000000000000000e+19",
            "CV_Te                       1.5000000000000000e+01",
            "IC_ne                       1.0000000000000000e+19",
            "IC_Te                       1.5000000000000000e+01",
            "restart_enabled             0",
            "pairSource_mean_x           0.0000000000000000e+00",
            "pairSource_positionMode     0",
            "pairSource_weightMode       1",
            "pairSource_Ti_birth         1.5000000000000000e+01",
            "pairSource_Te_birth         1.5000000000000000e+01",
        ],
    )

    smoke_input = root / "steady_interactive_1ms" / "inputFiles" / f"input_file_{SMOKE_TAG}.input"
    smoke_ion = root / "steady_interactive_1ms" / "inputFiles" / f"ions_properties_{SMOKE_TAG}.ion"
    assert_contains(
        smoke_input,
        [
            "SW_EfieldSolve              1",
            "SW_fieldSolveModel          2",
            "SW_Collisions               1",
            "SW_RFheating                0",
            "SW_pairSource               1",
            "IC_weightScale              1.0000000000000000e+00",
            "CV_ne                       1.0000000000000000e+19",
            "CV_Te                       1.5000000000000000e+01",
            "IC_ne                       1.0000000000000000e+19",
            "IC_Te                       1.5000000000000000e+01",
            "restart_enabled             0",
            "pairSource_mean_x           0.0000000000000000e+00",
            "pairSource_positionMode     0",
            "pairSource_weightMode       1",
            "pairSource_Ti_birth         1.5000000000000000e+01",
            "pairSource_Te_birth         1.5000000000000000e+01",
        ],
    )
    assert_contains(
        smoke_ion,
        [
            "BC_mean_x_1                   0.0000000000000000e+00",
            "BC_mean_x_2                   0.0000000000000000e+00",
            "IC_weightScale_1              1.0000000000000000e+00",
            "IC_weightScale_2              1.0000000000000000e+00",
        ],
    )

    for case in PRODUCTION_CASES[1:]:
        input_path = root / case.subdir / "inputFiles" / f"input_file_{case.tag}.input"
        ion_path = root / case.subdir / "inputFiles" / f"ions_properties_{case.tag}.ion"
        assert_contains(
            input_path,
            [
                "SW_EfieldSolve              1",
                "SW_fieldSolveModel          2",
                "SW_Collisions               1",
                "SW_RFheating                1",
                "SW_RFheatingIons            0",
                "SW_RFheatingElectrons       1",
                "SW_pairSource               1",
                "IC_weightScale              1.0000000000000000e+00",
                "CV_ne                       1.0000000000000000e+19",
                "CV_Te                       1.5000000000000000e+01",
                "IC_ne                       1.0000000000000000e+19",
                "IC_Te                       1.5000000000000000e+01",
                "restart_enabled             0",
                "pairSource_mean_x           0.0000000000000000e+00",
                "pairSource_positionMode     0",
                "pairSource_weightMode       1",
                "pairSource_Ti_birth         1.5000000000000000e+01",
                "pairSource_Te_birth         1.5000000000000000e+01",
                "RF_electron_Prf                      3.0000000000000000e+05",
            ],
        )
        assert_contains(
            ion_path,
            [
                "BC_mean_x_1                   0.0000000000000000e+00",
                "BC_mean_x_2                   0.0000000000000000e+00",
                "IC_weightScale_1              1.0000000000000000e+00",
                "IC_weightScale_2              1.0000000000000000e+00",
            ],
        )


def common_script() -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    run_picos_case() {{
      : "${{CASE_TAG:?CASE_TAG is required}}"
      : "${{CASE_SUBDIR:?CASE_SUBDIR is required}}"
      : "${{RESTART_REQUIRED:?RESTART_REQUIRED is required}}"
      : "${{STEADY_TAG:?STEADY_TAG is required}}"
      : "${{CASE_DIR:?CASE_DIR is required}}"

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

      HDF5_INSTALL=${{HDF5_INSTALL:-${{PICOS_ROOT}}/HDF5/lib}}
      ARMADILLO_INSTALL=${{ARMADILLO_INSTALL:-${{PICOS_ROOT}}/arma_libs/lib64}}
      if [ -d "${{HDF5_INSTALL}}" ]; then
        export LD_LIBRARY_PATH="${{LD_LIBRARY_PATH:-}}:${{HDF5_INSTALL}}"
      fi
      if [ -d "${{ARMADILLO_INSTALL}}" ]; then
        export LD_LIBRARY_PATH="${{LD_LIBRARY_PATH:-}}:${{ARMADILLO_INSTALL}}"
      fi

      if [ ! -x "${{PICOS_BIN}}" ]; then
        echo "Missing PICOS++ executable: ${{PICOS_BIN}}" >&2
        echo "Set PICOS_ROOT, PICOS_BUILD_DIR, or PICOS_BIN to the NERSC build." >&2
        exit 2
      fi
      if [ $((MPI_RANKS % 2)) -ne 0 ]; then
        echo "MPI_RANKS must be even for PICOS++." >&2
        exit 2
      fi

      STAGE_ROOT="${{RUN_ROOT}}/${{CASE_SUBDIR}}/run"
      RUN_PICOS_FILES="${{STAGE_ROOT}}/picosFILES"
      mkdir -p "${{RUN_PICOS_FILES}}/inputFiles" "${{RUN_PICOS_FILES}}/outputFiles" "${{STAGE_ROOT}}/logs"
      rsync -a "${{CASE_DIR}}/inputFiles/" "${{RUN_PICOS_FILES}}/inputFiles/"
      git -C "${{PICOS_ROOT}}" log --oneline -1 > "${{STAGE_ROOT}}/commitHash.txt" || true

      input_file="${{RUN_PICOS_FILES}}/inputFiles/input_file_${{CASE_TAG}}.input"
      ion_file="${{RUN_PICOS_FILES}}/inputFiles/ions_properties_${{CASE_TAG}}.ion"
      if [ ! -f "${{input_file}}" ] || [ ! -f "${{ion_file}}" ]; then
        echo "Missing input deck for ${{CASE_TAG}} under ${{RUN_PICOS_FILES}}/inputFiles." >&2
        exit 2
      fi

      if [ "${{RESTART_REQUIRED}}" = "1" ]; then
        STEADY_HDF5=${{STEADY_HDF5:-${{RUN_ROOT}}/steady_2ms/run/picosFILES/outputFiles/${{STEADY_TAG}}/HDF5}}
        if [ ! -d "${{STEADY_HDF5}}" ]; then
          echo "Missing steady-state restart directory: ${{STEADY_HDF5}}" >&2
          echo "Run the 2 ms steady-state job first or set STEADY_HDF5 explicitly." >&2
          exit 2
        fi
        tmp_file="${{input_file}}.tmp"
        awk -v restart_path="${{STEADY_HDF5}}" '
          /^restart_enabled[[:space:]]/ {{ print "restart_enabled             1"; next }}
          /^restart_path[[:space:]]/ {{ print "restart_path                " restart_path; next }}
          /^restart_snapshot[[:space:]]/ {{ print "restart_snapshot            -1"; next }}
          /^restart_continueTime[[:space:]]/ {{ print "restart_continueTime        0"; next }}
          {{ print }}
        ' "${{input_file}}" > "${{tmp_file}}"
        mv "${{tmp_file}}" "${{input_file}}"
      fi

      if [ -d "${{RUN_PICOS_FILES}}/outputFiles/${{CASE_TAG}}" ]; then
        if [ "${{OVERWRITE_OUTPUTS:-0}}" = "1" ]; then
          rm -rf "${{RUN_PICOS_FILES}}/outputFiles/${{CASE_TAG}}"
        else
          echo "Output exists: ${{RUN_PICOS_FILES}}/outputFiles/${{CASE_TAG}}" >&2
          echo "Set OVERWRITE_OUTPUTS=1 to replace it." >&2
          exit 2
        fi
      fi

      echo "[$(date)] Starting ${{CASE_TAG}}"
      echo "SLURM_JOB_ID=${{SLURM_JOB_ID:-interactive}}"
      echo "PICOS_ROOT=${{PICOS_ROOT}}"
      echo "PICOS_BIN=${{PICOS_BIN}}"
      echo "RUN_ROOT=${{RUN_ROOT}}"
      echo "STAGE_ROOT=${{STAGE_ROOT}}"
      echo "MPI_RANKS=${{MPI_RANKS}}"
      echo "CPUS_PER_TASK=${{CPUS_PER_TASK}}"
      echo "OMP_NUM_THREADS=${{OMP_NUM_THREADS}}"
      echo "LIBSTDCXX_DIR=${{LIBSTDCXX_DIR:-unset}}"
      if [ "${{RESTART_REQUIRED}}" = "1" ]; then
        echo "STEADY_HDF5=${{STEADY_HDF5}}"
      fi

      (
        cd "${{RUN_PICOS_FILES}}"
        srun -n "${{MPI_RANKS}}" -c "${{CPUS_PER_TASK}}" --cpu-bind=cores \\
          "${{PICOS_BIN}}" 1-D outputFiles "${{CASE_TAG}}"
      ) > "${{STAGE_ROOT}}/logs/${{CASE_TAG}}.log" 2>&1

      echo "[$(date)] Finished ${{CASE_TAG}}"
      echo "Output: ${{RUN_PICOS_FILES}}/outputFiles/${{CASE_TAG}}"
    }}
    """


def batch_script(case: BundleCase) -> str:
    restart = "1" if case.restart_required else "0"
    return f"""
    #!/bin/bash
    #SBATCH --account=m77
    #SBATCH -N 1
    #SBATCH -C cpu
    #SBATCH -q regular
    #SBATCH -J {case.job_name}
    #SBATCH --mail-user=kumara@ornl.gov
    #SBATCH --mail-type=END,FAIL
    #SBATCH -t {case.walltime}
    #SBATCH --ntasks-per-node=128
    #SBATCH --cpus-per-task=1
    #SBATCH -o slurm-%x-%j.out
    #SBATCH -e slurm-%x-%j.err

    set -euo pipefail

    CASE_TAG={case.tag}
    CASE_SUBDIR={case.subdir}
    STEADY_TAG={STEADY_TAG}
    RESTART_REQUIRED={restart}

    CASE_DIR=${{SLURM_SUBMIT_DIR:-$(pwd)}}
    if [ ! -f "${{CASE_DIR}}/inputFiles/input_file_${{CASE_TAG}}.input" ]; then
      echo "Submit this script from its case directory:" >&2
      echo "  cd <bundle>/{case.subdir} && sbatch ./run_case_nersc.sh" >&2
      exit 2
    fi

    source "${{CASE_DIR}}/run_case_common_nersc.sh"
    run_picos_case
    """


def interactive_script(case: BundleCase) -> str:
    restart = "1" if case.restart_required else "0"
    return f"""
    #!/bin/bash
    set -euo pipefail

    CASE_TAG={case.tag}
    CASE_SUBDIR={case.subdir}
    STEADY_TAG={STEADY_TAG}
    RESTART_REQUIRED={restart}
    CASE_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"

    source "${{CASE_DIR}}/run_case_common_nersc.sh"
    run_picos_case
    """


def submit_script() -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    ROOT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
    export RUN_ROOT=${{RUN_ROOT:-${{SCRATCH}}/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}}}

    echo "RUN_ROOT=${{RUN_ROOT}}"
    steady_job=$(cd "${{ROOT_DIR}}/steady_2ms" && sbatch --parsable ./run_case_nersc.sh)
    echo "Submitted 2 ms steady-state job: ${{steady_job}}"
    echo "Submitting ECH restart jobs with dependency afterok:${{steady_job}}"
    echo "Restart path for ECH jobs:"
    echo "  ${{RUN_ROOT}}/steady_2ms/run/picosFILES/outputFiles/{STEADY_TAG}/HDF5"

    for case_dir in left_resonance_70ghz right_resonance_70ghz well_min_65p588ghz; do
      job_id=$(cd "${{ROOT_DIR}}/${{case_dir}}" && sbatch --parsable --dependency=afterok:${{steady_job}} ./run_case_nersc.sh)
      echo "Submitted ${{case_dir}}: ${{job_id}}"
    done

    echo "Monitor with:"
    echo "  squeue -u $USER"
    echo "  sacct -j ${{steady_job}} --format=JobID,State,Elapsed,ExitCode"
    """


def run_interactive_1ms_script() -> str:
    return """
    #!/bin/bash
    set -euo pipefail

    ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    cd "${ROOT_DIR}/steady_interactive_1ms"
    ./run_case_nersc_interactive.sh
    """


def submit_ech_only_script() -> str:
    return """
    #!/bin/bash
    set -euo pipefail

    : "${STEADY_HDF5:?Set STEADY_HDF5=/path/to/steady/output/HDF5 before running this script.}"
    ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    export STEADY_HDF5
    export RUN_ROOT=${RUN_ROOT:-${SCRATCH}/picosRuns/MPEX_ECH_runs/PICOS_NERSC_MPEX_steady_then_ech_restart2ms_seeded_z0source}

    for case_dir in left_resonance_70ghz right_resonance_70ghz well_min_65p588ghz; do
      job_id=$(cd "${ROOT_DIR}/${case_dir}" && sbatch --parsable ./run_case_nersc.sh)
      echo "Submitted ${case_dir}: ${job_id}"
    done
    """


def readme_text() -> str:
    return f"""
    # MPEX Scenario-14 NERSC Steady-Then-ECH Restart Bundle

    This bundle runs a two-stage PICOS++ workflow on Perlmutter CPU nodes.

    Stage 1 is a 2 ms plasma steady-state run initialized with finite
    `n_e = n_i = 1e19 m^-3` and `T_e = T_i = 15 eV`, with collisions,
    coupled ion-electron pair source, and the reformulated Poisson field solve
    enabled, but ECH off.
    Stage 2 restarts from the last steady-state HDF5 snapshot and runs three
    100 us ECH cases:

    - `left_resonance_70ghz`
    - `right_resonance_70ghz`
    - `well_min_65p588ghz`

    The bundle also includes `steady_interactive_1ms`, an RF-off steady-state
    interactive test with the same physics switches, 32768 particles, and 1 ms
    physical time. Use it first inside an interactive allocation to validate
    modules, MPI launch, HDF5 output, source behavior, and early profile
    evolution before launching the full 2 ms batch stage.

    Important physics switches in these decks:

    - `SW_EfieldSolve = 1`
    - `SW_fieldSolveModel = 2`
    - `ReformulatedPoisson_quasiNeutral = 1`
    - `SW_Collisions = 1`
    - `SW_pairSource = 1`
    - `IC_weightScale = 1`
    - `pairSource_mean_x = 0`
    - `pairSource_positionMode = 0`
    - `pairSource_weightMode = 1`
    - `BC_mean_x_1 = BC_mean_x_2 = 0`
    - Stage 1: `SW_RFheating = 0`
    - Stage 2: `SW_RFheating = 1`, `SW_RFheatingIons = 0`, `SW_RFheatingElectrons = 1`

    The top-level submit script submits the steady-state job first and submits
    the three ECH jobs with `afterok` dependency. The ECH run scripts patch
    `restart_enabled`, `restart_path`, `restart_snapshot`, and
    `restart_continueTime` in their staged input decks at runtime.

    ## NERSC usage

    First smoke-test the steady-state path:

    ```bash
    cd /pscratch/sd/a/atul19/picosRuns/MPEX_ECH_runs
    tar -xzf {RUN_ROOT_NAME}.tar.gz
    cd {RUN_ROOT_NAME}

    salloc -A m77 -C cpu -q interactive -N 1 -t 04:00:00 --ntasks-per-node=128
    export PICOS_ROOT=/global/homes/a/atul19/myRepos/PICOS_ECH
    export PICOS_BUILD_DIR=${{PICOS_ROOT}}/build
    export RUN_ROOT=${{SCRATCH}}/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}_interactive_1ms
    export OVERWRITE_OUTPUTS=1
    ./run_steady_interactive_1ms_first.sh
    exit
    ```

    Then submit the production dependency chain:

    ```bash
    cd {RUN_ROOT_NAME}

    export PICOS_ROOT=/global/homes/a/atul19/myRepos/PICOS_ECH
    export PICOS_BUILD_DIR=${{PICOS_ROOT}}/build
    ./submit_steady_then_ech_nersc.sh
    ```

    Monitor:

    ```bash
    squeue -u $USER
    sacct -j <jobid> --format=JobID,JobName,State,Elapsed,ExitCode
    tail -f $SCRATCH/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}/steady_2ms/run/logs/{STEADY_TAG}.log
    ```

    Check the 1 ms interactive steady profiles after the test run:

    ```bash
    cd $PICOS_ROOT
    python3 scripts/check_mpex_steady_state_profiles.py \\
      $SCRATCH/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}_interactive_1ms \\
      --out-dir $SCRATCH/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}_interactive_1ms/steady_profile_check
    ```

    Check the full 2 ms steady profiles after the batch run:

    ```bash
    cd $PICOS_ROOT
    python3 scripts/check_mpex_steady_state_profiles.py \\
      $SCRATCH/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}/steady_2ms/run \\
      --out-dir $SCRATCH/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}/steady_2ms/steady_profile_check
    ```

    To rerun a case, set `OVERWRITE_OUTPUTS=1` before submission.

    If the steady state already exists, submit only ECH:

    ```bash
    export STEADY_HDF5=/pscratch/sd/a/atul19/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}/steady_2ms/run/picosFILES/outputFiles/{STEADY_TAG}/HDF5
    ./submit_ech_after_existing_steady_nersc.sh
    ```
    """


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
    source_dir = repo_root() / "picosFILES" / "inputFiles"

    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)

    for case in ALL_CASES:
        copy_case_inputs(root, case, source_dir)
        case_dir = root / case.subdir
        write_executable(case_dir / "run_case_common_nersc.sh", common_script())
        write_executable(case_dir / "run_case_nersc.sh", batch_script(case))
        write_executable(case_dir / "run_case_nersc_interactive.sh", interactive_script(case))

    write_executable(root / "submit_steady_then_ech_nersc.sh", submit_script())
    write_executable(root / "run_steady_interactive_1ms_first.sh", run_interactive_1ms_script())
    write_executable(root / "submit_ech_after_existing_steady_nersc.sh", submit_ech_only_script())
    (root / "README.md").write_text(textwrap.dedent(readme_text()).lstrip())

    validate_case_inputs(root)
    create_archive(root, archive)

    print(root)
    print(archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
