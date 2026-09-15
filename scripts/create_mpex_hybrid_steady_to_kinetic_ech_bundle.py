#!/usr/bin/env python3
"""Create a NERSC bundle for hybrid steady state followed by kinetic ECH.

The fully kinetic RF-off steady stage is not a good MPEX preconditioning step
on the present coarse 1D electrostatic grid: kinetic electrons drain before a
device-length column is established.  This bundle therefore runs the RF-off
steady stage in the original PICOS++ hybrid mode, converts the final hybrid
mesh profiles to kinetic-ion/electron IC profiles, and then runs a kinetic
electron ECH stage from those profiles.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import stat
import tarfile
import textwrap
from pathlib import Path


RUN_ROOT_NAME = "PICOS_NERSC_MPEX_hybrid_steady_to_kinetic_ech100us"
STEADY_SOURCE_TAG = "mpex_scenario14_ex8_steady_interactive_1ms_p32768_coll_mpexprof_nersc_nonrel"
ECH_SOURCE_TAG = "mpex_scenario14_ex8_well_min_65p588ghz_target8_p262144_100us_coll_mpexprof_nersc_nonrel"
STEADY_TAG = "mpex_scenario14_ex8_hybrid_steady_1ms_p32768_coll_mpexprof_nersc_nonrel"
ECH_TAG = "mpex_scenario14_ex8_kinetic_from_hybrid_well_min_65p588ghz_100us_coll_mpexprof_nersc_nonrel"
MPEX_NE = "5.0000000000000000e+19"
MPEX_TEMP = "1.5000000000000000e+01"
MPEX_SOURCE_RATE = "6.0000000000000000e+22"
MPEX_SOURCE_CENTER = "0.0000000000000000e+00"
MPEX_SOURCE_SIGMA = "2.9999999999999999e-01"
MPEX_PROFILE_N = 200
MPEX_ZMIN = -2.0
MPEX_ZMAX = 8.0


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def write_executable(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip())
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def replace_key(text: str, key: str, value: str, *, required: bool = True) -> str:
    pattern = re.compile(rf"^({re.escape(key)}\s+).*$", re.MULTILINE)
    updated, count = pattern.subn(rf"\g<1>{value}", text)
    if count == 0 and required:
        raise RuntimeError(f"Could not replace key {key}")
    return updated


def append_or_replace(text: str, key: str, value: str) -> str:
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
        text = append_or_replace(text, key, value)
    return text


def normalize_mpex_ions(text: str) -> str:
    replacements = {
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
        text = append_or_replace(text, key, value)
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


def read_input(input_dir: Path, tag: str) -> str:
    return (input_dir / f"input_file_{tag}.input").read_text()


def read_ions(input_dir: Path, tag: str) -> str:
    return (input_dir / f"ions_properties_{tag}.ion").read_text()


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


def truncate_after_species_one(ion_text: str) -> str:
    marker = "// Species 2:"
    idx = ion_text.find(marker)
    if idx >= 0:
        return ion_text[:idx].rstrip() + "\n"
    marker = "SPECIES2"
    idx = ion_text.find(marker)
    if idx >= 0:
        return ion_text[:idx].rstrip() + "\n"
    return ion_text


def steady_input(input_dir: Path) -> str:
    text = normalize_mpex_input(read_input(input_dir, STEADY_SOURCE_TAG))
    text = replace_key(text, "numberOfParticleSpecies", "1")
    text = replace_key(text, "SW_fieldSolveModel", "0")
    text = replace_key(text, "Poisson_BCModel", "2")
    text = replace_key(text, "SW_pairSource", "0")
    text = replace_key(text, "SW_RFheating", "0")
    text = replace_key(text, "SW_RFheatingIons", "0")
    text = replace_key(text, "SW_RFheatingElectrons", "0")
    text = replace_key(text, "restart_enabled", "0")
    text = replace_key(text, "restart_path", "none")
    text = replace_key(
        text,
        "outputs_variables",
        "{X_p,V_p,a_p,BX_p,EX_p,BX_m,dBX_m,ddBX_m,n_m,Tpar_m,Tper_m,Te_m,u_m,EX_m}",
    )
    return text.replace("kinetic D+ plus kinetic electrons", "hybrid D+ plus fluid electrons")


def steady_ions(input_dir: Path) -> str:
    text = truncate_after_species_one(read_ions(input_dir, STEADY_SOURCE_TAG))
    # Use main-branch optimized source rate scale as the hybrid preconditioner.
    text = replace_key(text, "BC_G_1", MPEX_SOURCE_RATE)
    text = replace_key(text, "BC_T_1", MPEX_TEMP)
    text = replace_key(text, "BC_type_1", "1")
    text = replace_key(text, "BC_mean_x_1", MPEX_SOURCE_CENTER)
    text = replace_key(text, "BC_sigma_x_1", MPEX_SOURCE_SIGMA)
    text = append_or_replace(text, "IC_weightScale_1", "1.0000000000000000e+00")
    return text


def ech_input(input_dir: Path) -> str:
    text = normalize_mpex_input(read_input(input_dir, ECH_SOURCE_TAG))
    text = replace_key(text, "Poisson_BCModel", "2")
    text = replace_key(text, "SW_RFheating", "1")
    text = replace_key(text, "SW_RFheatingIons", "0")
    text = replace_key(text, "SW_RFheatingElectrons", "1")
    text = replace_key(text, "SW_pairSource", "0")
    text = replace_key(text, "restart_enabled", "0")
    text = replace_key(text, "restart_path", "none")
    text = replace_key(text, "restart_continueTime", "0")
    text = replace_key(
        text,
        "IC_BX_fileName",
        f"{ECH_TAG}_B_norm_from_hybrid.txt",
    )
    text = replace_key(
        text,
        "IC_Te_fileName",
        f"{ECH_TAG}_te_norm_from_hybrid.txt",
    )
    text = replace_key(
        text,
        "outputs_variables",
        "{X_p,V_p,a_p,mu_p,BX_p,EX_p,BX_m,dBX_m,ddBX_m,n_m,Tpar_m,Tper_m,u_m,EX_m,Phi_m}",
    )
    return text


def ech_ions(input_dir: Path) -> str:
    text = normalize_mpex_ions(read_ions(input_dir, ECH_SOURCE_TAG))
    replacements = {
        "IC_Tper_fileName_1": f"{ECH_TAG}_ti_norm_from_hybrid.txt",
        "IC_Tpar_fileName_1": f"{ECH_TAG}_ti_norm_from_hybrid.txt",
        "IC_densityFraction_fileName_1": f"{ECH_TAG}_ne_norm_from_hybrid.txt",
        "IC_Tper_fileName_2": f"{ECH_TAG}_te_norm_from_hybrid.txt",
        "IC_Tpar_fileName_2": f"{ECH_TAG}_te_norm_from_hybrid.txt",
        "IC_densityFraction_fileName_2": f"{ECH_TAG}_ne_norm_from_hybrid.txt",
        # Preserve computational weights during the short kinetic ECH stage.
        "BC_type_1": "4",
        "BC_type_2": "4",
    }
    for key, value in replacements.items():
        text = replace_key(text, key, value)
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
        if src.is_file():
            profile_text = paper_profile_text(name)
            if profile_text is None:
                shutil.copyfile(src, input_dir / name)
            else:
                (input_dir / name).write_text(profile_text)
    write_executable(case_dir / "run_case_common_nersc.sh", common_script())
    write_executable(case_dir / "run_case_nersc_interactive.sh", interactive_case_script(case_dir.name))


def converter_script() -> str:
    return f"""
    #!/usr/bin/env python3
    from __future__ import annotations

    import argparse
    from pathlib import Path

    import h5py
    import numpy as np

    TAG = "{ECH_TAG}"

    def numeric_steps(path: Path) -> list[str]:
        with h5py.File(path, "r") as handle:
            return sorted([k for k in handle.keys() if k.isdigit()], key=int)

    def read(path: Path, dataset: str) -> np.ndarray:
        with h5py.File(path, "r") as handle:
            return np.asarray(handle[dataset], dtype=float).reshape(-1)

    def interp_to_profile(z: np.ndarray, values: np.ndarray, n: int, zmin: float, zmax: float) -> np.ndarray:
        xp = np.linspace(zmin, zmax, n)
        order = np.argsort(z)
        out = np.interp(xp, z[order], values[order], left=values[order][0], right=values[order][-1])
        return out

    def normalized(values: np.ndarray, base: float, floor: float = 1.0e-8) -> np.ndarray:
        out = np.asarray(values, dtype=float) / max(base, 1.0e-300)
        out[~np.isfinite(out)] = floor
        return np.maximum(out, floor)

    def main() -> int:
        parser = argparse.ArgumentParser()
        parser.add_argument("--steady-hdf5", required=True, type=Path)
        parser.add_argument("--ech-input-dir", required=True, type=Path)
        parser.add_argument("--n", type=int, default=200)
        parser.add_argument("--cv-ne", type=float, default=5.0e19)
        parser.add_argument("--ti-base", type=float, default=15.0)
        parser.add_argument("--te-base", type=float, default=15.0)
        parser.add_argument("--b-base", type=float, default=1.3584640814139803)
        args = parser.parse_args()

        main_h5 = args.steady_hdf5 / "main.h5"
        particles = args.steady_hdf5 / "PARTICLES_FILE_0.h5"
        fields = args.steady_hdf5 / "FIELDS_FILE_0.h5"
        if not main_h5.is_file() or not particles.is_file() or not fields.is_file():
            raise FileNotFoundError(f"Missing HDF5 restart output under {{args.steady_hdf5}}")

        with h5py.File(main_h5, "r") as handle:
            z = np.asarray(handle["geometry/x_m"], dtype=float).reshape(-1)
        step = numeric_steps(particles)[-1]
        fstep = numeric_steps(fields)[-1]

        n = read(particles, f"/{{step}}/ions/spp_1/n_m")
        ti = read(particles, f"/{{step}}/ions/spp_1/Tpar_m")
        b = read(fields, f"/{{fstep}}/fields/BX_m/x")
        try:
            te = read(particles, f"/{{step}}/ions/spp_1/Te_m")
        except Exception:
            te = np.full_like(n, args.te_base)

        zmin = float(np.nanmin(z))
        zmax = float(np.nanmax(z))
        # HDF5 x_m is cell-centered; input profiles span LX_min..LX_max.
        lx_min = -2.0
        lx_max = 8.0
        n_prof = interp_to_profile(z, n, args.n, lx_min, lx_max)
        ti_prof = interp_to_profile(z, ti, args.n, lx_min, lx_max)
        te_prof = interp_to_profile(z, te, args.n, lx_min, lx_max)
        b_prof = interp_to_profile(z, b, args.n, lx_min, lx_max)

        args.ech_input_dir.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.ech_input_dir / f"{{TAG}}_ne_norm_from_hybrid.txt", normalized(n_prof, args.cv_ne), fmt="%.16e")
        np.savetxt(args.ech_input_dir / f"{{TAG}}_ti_norm_from_hybrid.txt", normalized(ti_prof, args.ti_base), fmt="%.16e")
        np.savetxt(args.ech_input_dir / f"{{TAG}}_te_norm_from_hybrid.txt", normalized(te_prof, args.te_base), fmt="%.16e")
        np.savetxt(args.ech_input_dir / f"{{TAG}}_B_norm_from_hybrid.txt", normalized(b_prof, args.b_base), fmt="%.16e")
        print(f"Wrote kinetic ECH IC profiles from hybrid steady snapshot {{step}} to {{args.ech_input_dir}}")
        return 0

    if __name__ == "__main__":
        raise SystemExit(main())
    """


def common_script() -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    run_picos_case() {{
      : "${{CASE_TAG:?CASE_TAG is required}}"
      : "${{CASE_SUBDIR:?CASE_SUBDIR is required}}"

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

      CASE_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
      STAGE_ROOT="${{RUN_ROOT}}/${{CASE_SUBDIR}}/run"
      RUN_PICOS_FILES="${{STAGE_ROOT}}/picosFILES"
      mkdir -p "${{RUN_PICOS_FILES}}/inputFiles" "${{STAGE_ROOT}}/logs"
      rsync -a "${{CASE_DIR}}/inputFiles/" "${{RUN_PICOS_FILES}}/inputFiles/"
      git -C "${{PICOS_ROOT}}" log --oneline -1 > "${{STAGE_ROOT}}/commitHash.txt" || true

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
    case_tag = STEADY_TAG if subdir == "hybrid_steady_1ms" else ECH_TAG
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
    source "$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)/run_case_common_nersc.sh"
    run_picos_case
    """


def top_level_run_script() -> str:
    return f"""
    #!/bin/bash
    set -euo pipefail

    if [ -z "${{SLURM_JOB_ID:-}}" ]; then
      echo "Start an interactive allocation first:" >&2
      echo "  salloc -A m77 -C cpu -q interactive -N 1 -t 04:00:00 --ntasks-per-node=128" >&2
      exit 2
    fi

    ROOT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
    RUN_ROOT=${{RUN_ROOT:-${{SCRATCH}}/picosRuns/MPEX_ECH_runs/{RUN_ROOT_NAME}}}
    export RUN_ROOT

    cd "${{ROOT_DIR}}/hybrid_steady_1ms"
    ./run_case_nersc_interactive.sh

    STEADY_HDF5="${{RUN_ROOT}}/hybrid_steady_1ms/run/picosFILES/outputFiles/HDF5"
    ECH_STAGE="${{RUN_ROOT}}/kinetic_well_min_ech100us/run/picosFILES"
    mkdir -p "${{ECH_STAGE}}/inputFiles"
    rsync -a "${{ROOT_DIR}}/kinetic_well_min_ech100us/inputFiles/" "${{ECH_STAGE}}/inputFiles/"

    python "${{ROOT_DIR}}/make_kinetic_ic_from_hybrid_steady.py" \\
      --steady-hdf5 "${{STEADY_HDF5}}" \\
      --ech-input-dir "${{ECH_STAGE}}/inputFiles"

    cd "${{ROOT_DIR}}/kinetic_well_min_ech100us"
    ./run_case_nersc_interactive.sh
    """


def readme_text() -> str:
    return f"""
    # MPEX Hybrid Steady To Kinetic ECH 100 us

    This bundle is meant to replace the fully kinetic RF-off pre-run that
    collapsed to the z=0 source region.

    Workflow:

    1. `hybrid_steady_1ms`: one kinetic D+ species with fluid electrons,
       flat `n/T` initial profiles, `n_e=n_i=5e19 m^-3`,
       `T_e=T_i=15 eV`, `G=6e22 s^-1` at `z=0`, generalized Ohm
       field solve, collisions on, RF off.
    2. `make_kinetic_ic_from_hybrid_steady.py`: converts the final hybrid mesh
       density, Ti, Te, and B profiles into normalized IC files.
    3. `kinetic_well_min_ech100us`: kinetic D+ plus kinetic electrons, ECH on
       for electrons only. Boundary type is set to simple same-weight reinjection
       for this short kinetic ECH stage so the initial column is not destroyed by
       explicit-rate electron recycling.

    ## Run on NERSC interactive node

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

    ./run_hybrid_steady_then_kinetic_ech100us.sh
    ```

    Check:

    ```bash
    find $RUN_ROOT -maxdepth 8 -type f -name '*.h5' -print | head
    tail -80 $RUN_ROOT/hybrid_steady_1ms/run/logs/*.log
    tail -80 $RUN_ROOT/kinetic_well_min_ech100us/run/logs/*.log
    ```
    """


def validate(root: Path) -> None:
    steady = (root / "hybrid_steady_1ms/inputFiles/input_file.input").read_text()
    steady_ion = (root / "hybrid_steady_1ms/inputFiles/ions_properties.ion").read_text()
    ech = (root / "kinetic_well_min_ech100us/inputFiles/input_file.input").read_text()
    ech_ion = (root / "kinetic_well_min_ech100us/inputFiles/ions_properties.ion").read_text()
    checks = {
        "steady input": [
            "numberOfParticleSpecies     1",
            "quietStart                  1",
            f"CV_ne                       {MPEX_NE}",
            f"IC_ne                       {MPEX_NE}",
            "SW_fieldSolveModel          0",
            "SW_pairSource               0",
            "SW_RFheating                0",
            "Te_m",
        ],
        "steady ion": [
            f"BC_G_1                        {MPEX_SOURCE_RATE}",
            f"BC_sigma_x_1                  {MPEX_SOURCE_SIGMA}",
            "BC_mean_x_1                   0.0000000000000000e+00",
        ],
        "ech input": [
            "quietStart                  1",
            f"CV_ne                       {MPEX_NE}",
            f"IC_ne                       {MPEX_NE}",
            "SW_RFheating                1",
            "SW_RFheatingElectrons       1",
            "SW_RFheatingIons            0",
            "SW_pairSource               0",
            f"{ECH_TAG}_te_norm_from_hybrid.txt",
        ],
        "ech ion": [
            "BC_type_1                     4",
            "BC_type_2                     4",
            f"{ECH_TAG}_ne_norm_from_hybrid.txt",
        ],
    }
    texts = {
        "steady input": steady,
        "steady ion": steady_ion,
        "ech input": ech,
        "ech ion": ech_ion,
    }
    for label, needles in checks.items():
        missing = [needle for needle in needles if needle not in texts[label]]
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

    write_case(root / "hybrid_steady_1ms", steady_input(input_dir), steady_ions(input_dir), input_dir)
    write_case(root / "kinetic_well_min_ech100us", ech_input(input_dir), ech_ions(input_dir), input_dir)
    write_executable(root / "make_kinetic_ic_from_hybrid_steady.py", converter_script())
    write_executable(root / "run_hybrid_steady_then_kinetic_ech100us.sh", top_level_run_script())
    (root / "README.md").write_text(textwrap.dedent(readme_text()).lstrip())

    validate(root)
    create_archive(root, archive)
    print(root)
    print(archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
