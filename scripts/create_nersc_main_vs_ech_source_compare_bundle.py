#!/usr/bin/env python3
"""Create a NERSC bundle comparing PICOS main and PICOS_ECH MPEX hybrid runs."""

from __future__ import annotations

import argparse
import shutil
import stat
import tarfile
import textwrap
from dataclasses import dataclass
from pathlib import Path


RUN_ROOT_NAME = "PICOS_NERSC_main_vs_ech_MPEX_scenario14_mainformat_1ms"
PROFILE_FILES = (
    "MPEX_B_norm_PICOS_scenario_14.txt",
    "MPEX_Tpar_norm_scenario_14.txt",
    "MPEX_Tper_norm_scenario_14.txt",
    "MPEX_n_norm_scenario_14.txt",
)


@dataclass(frozen=True)
class CompareCase:
    subdir: str
    tag: str
    job_name: str
    root_env: str
    build_env: str
    bin_env: str
    default_root: str


CASES = (
    CompareCase(
        "main_mpex_scenario14_1ms",
        "picos_main_mpex_scenario14_mainformat_1ms",
        "picos_m14",
        "PICOS_MAIN_ROOT",
        "PICOS_MAIN_BUILD_DIR",
        "PICOS_MAIN_BIN",
        "${HOME}/myRepos/PICOS_main",
    ),
    CompareCase(
        "ech_branch_mpex_scenario14_1ms",
        "picos_ech_mpex_scenario14_mainformat_1ms",
        "picos_e14",
        "PICOS_ECH_ROOT",
        "PICOS_ECH_BUILD_DIR",
        "PICOS_ECH_BIN",
        "${HOME}/myRepos/PICOS_ECH",
    ),
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def write_executable(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip())
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_vector(path: Path, values: list[float]) -> None:
    path.write_text("".join(f"{value:.16e}\n" for value in values))


def input_deck(tag: str) -> str:
    return textwrap.dedent(
        f"""
        // PICOS main-vs-PICOS_ECH MPEX Scenario 14 source comparison.
        // Common physics: one D+ guiding-center kinetic species, fluid
        // electrons, Ohm-law electric field, collisions on, RF off, and the
        // legacy warm source centered at z=0.
        // This deck intentionally uses only main-branch input keywords, with
        // SW_Bohm and Bohm_* included, so the main and PICOS_ECH branches run
        // the same hybrid model.
        // =============================================================================
        mpisForFields               2
        quietStart                  1
        numberOfParticleSpecies     1
        numberOfTracerSpecies       0
        advanceParticleMethod       1

        // Characteristic values:
        // =============================================================================
        CV_ne                       1.0000000000000000e+19
        CV_Te                       1.5000000000000000e+01
        CV_B                        1.3079603575279184e+00
        CV_Tpar                     1.5000000000000000e+01
        CV_Tper                     1.5000000000000000e+01

        // Simulation time is in background-ion gyroperiod units. For the MPEX
        // Scenario 14 field this is about 1 ms.
        // =============================================================================
        DTc                         5.0000000000000003e-01
        simulationTime              1.2000000000000000e+04

        // Switches:
        // =============================================================================
        SW_EfieldSolve               1
        SW_BfieldSolve               0
        SW_Collisions                1
        SW_RFheating                 0
        SW_advancePos                1
        SW_linearSolve               0
        SW_Bohm                      0

        // Magnetic field initial conditions:
        // =============================================================================
        IC_uniformBfield            0
        IC_BX                       1.3079603575279184e+00
        IC_BY                       0.0
        IC_BZ                       0.0
        IC_BX_NX                    200
        IC_BX_fileName              MPEX_B_norm_PICOS_scenario_14.txt

        // Geometry:
        // =============================================================================
        dp                          4.9774922103303126e-01
        r1                          0.0
        r2                          5.0000000000000003e-02
        LX_min                      -2.0000000000000000e+00
        LX_max                      8.0000000000000000e+00

        // Electron fluid/profile initial conditions:
        // =============================================================================
        IC_ne                       1.0000000000000000e+19
        IC_Te                       1.5000000000000000e+01
        IC_Te_NX                    200
        IC_Te_fileName              MPEX_Tper_norm_scenario_14.txt

        // RF operator is disabled, but both branches still parse these keys.
        // =============================================================================
        RF_Prf                      0.0
        RF_n_harmonic               1
        RF_freq                     8.3850000000000000e+06
        RF_x1                       4.0000000000000000e+00
        RF_x2                       6.5000000000000000e+00
        RF_t_ON                     1.2000000000000000e+04
        RF_t_OFF                    1.2000000000000000e+04
        RF_kpar                     2.0000000000000000e+01
        RF_kper                     1.0000000000000000e+02
        RF_handedness               -1
        RF_Prf_fileName             Prf_profile.txt
        RF_Prf_NS                   200

        // Bohm boundary condition block. SW_Bohm=0 here, but PICOS_main expects
        // these legacy keys to exist in the input file.
        // =============================================================================
        Bohm_type                    2
        Bohm_edgeCells               4
        Bohm_t_ON                    50
        Bohm_gamma_i                 3

        // Output variables:
        // =============================================================================
        outputCadence               5.0000000000000000e+02
        outputs_variables           {{X_p,V_p,a_p,BX_p,BX_m,n_m,Tpar_m,Tper_m,Te_m,u_m,EX_m}}

        // Data smoothing:
        // =============================================================================
        smoothingParameter          1.0000000000000000e-04
        filtersPerIterationFields   2
        filtersPerIterationIons     2
        """
    ).strip() + "\n"


def ion_deck() -> str:
    return textwrap.dedent(
        """
        // PICOS main-vs-PICOS_ECH MPEX Scenario 14 species deck.
        // Species 1 is D+. The source is the legacy warm plasma source
        // centered at z=0 for this hybrid main-branch compatibility test.
        // =============================================================================
        SPECIES1                      1
        NPC1                          2500
        pctSupPartOutput1             100
        Z1                            1
        M1                            2.0000000000000000e+00

        IC_type_1                     1
        IC_Tper_1                     1.5000000000000000e+01
        IC_Tper_fileName_1            MPEX_Tper_norm_scenario_14.txt
        IC_Tper_NX_1                  200
        IC_Tpar_1                     1.5000000000000000e+01
        IC_Tpar_fileName_1            MPEX_Tpar_norm_scenario_14.txt
        IC_Tpar_NX_1                  200
        IC_densityFraction_1          1.0
        IC_densityFraction_fileName_1 MPEX_n_norm_scenario_14.txt
        IC_densityFraction_NX_1       200

        // Boundary conditions:
        // 1 = warm plasma source, 2 = NBI, 3 = periodic, 4 = simple reinjection.
        // =============================================================================
        BC_type_1                     1
        BC_T_1                        1.5000000000000000e+01
        BC_E_1                        0.0
        BC_eta_1                      45
        BC_mean_x_1                   0.0000000000000000e+00
        BC_sigma_x_1                  4.0000000000000002e-01
        BC_G_1                        1.0000000000000000e+22
        BC_G_fileName_1               G_profile.txt
        BC_G_NS_1                     200
        """
    ).strip() + "\n"


def common_script(case: CompareCase) -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    run_compare_case() {{
      : "${{CASE_TAG:?CASE_TAG is required}}"
      : "${{CASE_SUBDIR:?CASE_SUBDIR is required}}"
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

      PICOS_ROOT=${{{case.root_env}:-{case.default_root}}}
      PICOS_BUILD_DIR=${{{case.build_env}:-${{PICOS_ROOT}}/build}}
      PICOS_BIN=${{{case.bin_env}:-${{PICOS_BUILD_DIR}}/picosFILES/src/xpicos}}
      RUN_ROOT=${{RUN_ROOT:-${{SCRATCH}}/picosRuns/{RUN_ROOT_NAME}}}
      MPI_RANKS=${{MPI_RANKS:-128}}
      CPUS_PER_TASK=${{CPUS_PER_TASK:-1}}

      export OMP_NUM_THREADS=${{OMP_NUM_THREADS:-1}}
      export OMP_PLACES=${{OMP_PLACES:-threads}}
      export OMP_PROC_BIND=${{OMP_PROC_BIND:-spread}}
      export HDF5_USE_FILE_LOCKING=FALSE

      if [ ! -x "${{PICOS_BIN}}" ]; then
        echo "Missing PICOS++ executable: ${{PICOS_BIN}}" >&2
        echo "Set {case.root_env}, {case.build_env}, or {case.bin_env}." >&2
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

      if [ ! -f "${{RUN_PICOS_FILES}}/inputFiles/input_file_${{CASE_TAG}}.input" ]; then
        echo "Missing input_file_${{CASE_TAG}}.input under ${{RUN_PICOS_FILES}}/inputFiles." >&2
        exit 2
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
      echo "PICOS_ROOT=${{PICOS_ROOT}}"
      echo "PICOS_BIN=${{PICOS_BIN}}"
      echo "RUN_ROOT=${{RUN_ROOT}}"
      echo "STAGE_ROOT=${{STAGE_ROOT}}"
      echo "MPI_RANKS=${{MPI_RANKS}}"
      echo "OMP_NUM_THREADS=${{OMP_NUM_THREADS}}"

      (
        cd "${{RUN_PICOS_FILES}}"
        srun -n "${{MPI_RANKS}}" -c "${{CPUS_PER_TASK}}" --cpu-bind=cores \\
          "${{PICOS_BIN}}" 1-D outputFiles "${{CASE_TAG}}"
      ) > "${{STAGE_ROOT}}/logs/${{CASE_TAG}}.log" 2>&1

      echo "[$(date)] Finished ${{CASE_TAG}}"
      echo "Output: ${{RUN_PICOS_FILES}}/outputFiles/${{CASE_TAG}}"
    }}
    """


def batch_script(case: CompareCase) -> str:
    return f"""
    #!/bin/bash
    #SBATCH --account=m77
    #SBATCH -N 1
    #SBATCH -C cpu
    #SBATCH -q regular
    #SBATCH -J {case.job_name}
    #SBATCH -t 04:00:00
    #SBATCH --ntasks-per-node=128
    #SBATCH --cpus-per-task=1
    #SBATCH -o slurm-%x-%j.out
    #SBATCH -e slurm-%x-%j.err

    set -euo pipefail

    CASE_TAG={case.tag}
    CASE_SUBDIR={case.subdir}
    CASE_DIR=${{SLURM_SUBMIT_DIR:-$(pwd)}}
    if [ ! -f "${{CASE_DIR}}/inputFiles/input_file_${{CASE_TAG}}.input" ]; then
      echo "Submit from this case directory:" >&2
      echo "  cd <bundle>/{case.subdir} && sbatch ./run_case_nersc.sh" >&2
      exit 2
    fi

    source "${{CASE_DIR}}/run_case_common_nersc.sh"
    run_compare_case
    """


def interactive_script(case: CompareCase) -> str:
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
    CASE_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
    source "${{CASE_DIR}}/run_case_common_nersc.sh"
    run_compare_case
    """


def submit_script() -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    ROOT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
    export RUN_ROOT=${{RUN_ROOT:-${{SCRATCH}}/picosRuns/{RUN_ROOT_NAME}}}
    echo "RUN_ROOT=${{RUN_ROOT}}"

    main_job=$(cd "${{ROOT_DIR}}/{CASES[0].subdir}" && sbatch --parsable ./run_case_nersc.sh)
    ech_job=$(cd "${{ROOT_DIR}}/{CASES[1].subdir}" && sbatch --parsable ./run_case_nersc.sh)

    echo "Submitted main branch job: ${{main_job}}"
    echo "Submitted PICOS_ECH branch job: ${{ech_job}}"
    echo "Monitor with:"
    echo "  squeue -u $USER"
    echo "  sacct -j ${{main_job}},${{ech_job}} --format=JobID,JobName,State,Elapsed,ExitCode"
    """


def run_both_interactive_script() -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    if [ -z "${{SLURM_JOB_ID:-}}" ]; then
      echo "Start an interactive allocation first, for example:" >&2
      echo "  salloc -A m77 -C cpu -q interactive -N 1 -t 04:00:00 --ntasks-per-node=128" >&2
      exit 2
    fi

    ROOT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
    cd "${{ROOT_DIR}}/{CASES[0].subdir}"
    ./run_case_nersc_interactive.sh
    cd "${{ROOT_DIR}}/{CASES[1].subdir}"
    ./run_case_nersc_interactive.sh
    """


def readme_text() -> str:
    return f"""
    # PICOS main vs PICOS_ECH hybrid source comparison

    This bundle runs the same master-compatible MPEX scenario-14 hybrid source
    deck with two executables:

    - `{CASES[0].subdir}`: `~/myRepos/PICOS_main/build/picosFILES/src/xpicos`
    - `{CASES[1].subdir}`: `~/myRepos/PICOS_ECH/build/picosFILES/src/xpicos`

    The physics is intentionally restricted to features common to both
    branches: one D+ guiding-center kinetic species, fluid electrons, Ohm-law
    electric field solve, collisions on, RF off, and the legacy z=0 warm source.
    This is the right comparison to check that the PICOS_ECH branch still
    reproduces the main-branch steady-source behavior before enabling kinetic
    electrons, explicit pair source, Poisson solve, or ECH.

    ## Batch run

    ```bash
    cd $SCRATCH/picosRuns
    tar -xzf {RUN_ROOT_NAME}.tar.gz
    cd {RUN_ROOT_NAME}

    export PICOS_MAIN_ROOT=$HOME/myRepos/PICOS_main
    export PICOS_ECH_ROOT=$HOME/myRepos/PICOS_ECH
    ./submit_both_nersc.sh
    ```

    ## Interactive run

    ```bash
    cd $SCRATCH/picosRuns/{RUN_ROOT_NAME}
    salloc -A m77 -C cpu -q interactive -N 1 -t 04:00:00 --ntasks-per-node=128
    export PICOS_MAIN_ROOT=$HOME/myRepos/PICOS_main
    export PICOS_ECH_ROOT=$HOME/myRepos/PICOS_ECH
    export OVERWRITE_OUTPUTS=1
    ./run_both_interactive.sh
    exit
    ```

    ## Compare outputs

    After both jobs finish:

    ```bash
    python analyze_nersc_main_vs_ech_source_compare.py \\
      --run-root $SCRATCH/picosRuns/{RUN_ROOT_NAME}
    ```

    Expected outputs:

    - `analysis/main_vs_ech_hybrid_summary.csv`
    - `analysis/main_vs_ech_hybrid_final_profiles.png`
    - `analysis/main_vs_ech_hybrid_histories.png`
    """


def populate_case(root: Path, case: CompareCase, profile_dir: Path) -> None:
    case_dir = root / case.subdir
    input_dir = case_dir / "inputFiles"
    input_dir.mkdir(parents=True, exist_ok=True)
    (input_dir / f"input_file_{case.tag}.input").write_text(input_deck(case.tag))
    (input_dir / f"ions_properties_{case.tag}.ion").write_text(ion_deck())
    for name in PROFILE_FILES:
        src = profile_dir / name
        if not src.is_file():
            raise FileNotFoundError(src)
        shutil.copy2(src, input_dir / name)
    write_vector(input_dir / "G_profile.txt", [1.0] * 200)
    write_vector(input_dir / "Prf_profile.txt", [0.0] * 200)
    write_executable(case_dir / "run_case_common_nersc.sh", common_script(case))
    write_executable(case_dir / "run_case_nersc.sh", batch_script(case))
    write_executable(case_dir / "run_case_nersc_interactive.sh", interactive_script(case))


def create_archive(root: Path, archive: Path) -> None:
    if archive.exists():
        archive.unlink()
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(root, arcname=root.name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, default=repo_root() / "validation" / RUN_ROOT_NAME)
    parser.add_argument("--archive", type=Path, default=None)
    parser.add_argument("--profile-dir", type=Path, default=repo_root() / "templateFILES")
    args = parser.parse_args()

    root = args.out_root.resolve()
    archive = (args.archive if args.archive else root.with_suffix(".tar.gz")).resolve()
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)

    for case in CASES:
        populate_case(root, case, args.profile_dir.resolve())

    write_executable(root / "submit_both_nersc.sh", submit_script())
    write_executable(root / "run_both_interactive.sh", run_both_interactive_script())
    shutil.copy2(repo_root() / "scripts" / "analyze_nersc_main_vs_ech_source_compare.py", root)
    (root / "README.md").write_text(textwrap.dedent(readme_text()).lstrip())
    create_archive(root, archive)

    print(root)
    print(archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
