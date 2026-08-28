#!/usr/bin/env python3
"""Cross-platform driver for this project's two WarpX builds.

Standard library only, and targets Python 3.8+: this runs *before* .venv exists,
so it cannot import anything from it.

    warpx info                 show what was detected, and build nothing
    warpx root                 print the project root
    warpx sync [uv args]       build pywarpx into .venv
    warpx rebuild [uv args]    force a pywarpx rebuild
    warpx build [cmake args]   build the standalone solver
    warpx run [warpx args]     run the standalone solver

Configuration resolves in one direction, most explicit first:

    command-line flag  >  environment variable  >  autodetection

so `--compute cuda`, then `WARPX_COMPUTE=CUDA`, then whatever the machine looks
like. `warpx info` prints the result of that resolution without building.
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

DIMS = "3"          # what WarpX_DIMS / $WARPX_DIMS take
DIM_SUFFIX = "3d"   # what WarpX appends to the binary name for those dims
COMPUTE_CHOICES = ("auto", "NOACC", "OMP", "CUDA", "SYCL", "HIP")

WINDOWS = platform.system() == "Windows"
MACOS = platform.system() == "Darwin"


# --------------------------------------------------------------------- errors


class WarpXError(RuntimeError):
    """Something the user has to fix; reported without a traceback."""


# ----------------------------------------------------------------------- root


def find_root(explicit: str | None = None) -> Path:
    """Locate the project root: $WARPX_PROJECT, then upwards from cwd, then
    upwards from this file (so the driver works from outside the project)."""
    if explicit:
        return Path(explicit).resolve()
    env = os.environ.get("WARPX_PROJECT")
    if env:
        return Path(env).resolve()

    for start in (Path.cwd(), Path(__file__).resolve().parent):
        for d in (start.resolve(), *start.resolve().parents):
            if (d / "vendor" / "warpx").is_dir() and (d / "pyproject.toml").is_file():
                return d

    raise WarpXError(
        "not inside a WarpX project; set $WARPX_PROJECT to the project root"
    )


def require_submodule(root: Path) -> None:
    if not (root / "vendor" / "warpx" / "CMakeLists.txt").is_file():
        raise WarpXError(
            "vendor/warpx is empty — run: git submodule update --init --recursive"
        )


# ------------------------------------------------------------------ toolchain


def have(prog: str) -> bool:
    return shutil.which(prog) is not None


def cuda_available() -> bool:
    """True when there is a CUDA *compiler*, not merely a CUDA driver.

    nvidia-smi is deliberately not enough: plenty of machines can run CUDA
    binaries but cannot compile them, and guessing CUDA there swaps a working
    build for a confusing failure.
    """
    if have("nvcc") or os.environ.get("CUDACXX"):
        return True
    for var in ("CUDA_HOME", "CUDA_PATH"):
        home = os.environ.get(var)
        if home:
            exe = "nvcc.exe" if WINDOWS else "nvcc"
            if (Path(home) / "bin" / exe).is_file():
                return True
    return False


def is_omp_prefix(d: Path) -> bool:
    """True if d looks like a prefix holding an OpenMP runtime."""
    if not (d / "include" / "omp.h").is_file():
        return False
    libs = ("libomp.dylib", "libomp.so", "libomp.a", "libgomp.dylib", "libgomp.so")
    return any((d / "lib" / lib).is_file() for lib in libs)


def omp_candidates() -> list[Path]:
    """Prefixes that might hold libomp, most specific first. No package manager
    is required — Homebrew is one entry here, not the mechanism."""
    out: list[str] = []
    conda = os.environ.get("CONDA_PREFIX")
    if conda:
        out.append(conda)
    brew = os.environ.get("HOMEBREW_PREFIX")
    if brew:
        out += [os.path.join(brew, "opt", "libomp"), brew]
    if have("brew"):
        try:
            out.append(
                subprocess.run(
                    ["brew", "--prefix", "libomp"],
                    capture_output=True, text=True, check=True,
                ).stdout.strip()
            )
        except (subprocess.CalledProcessError, OSError):
            pass
    out += [
        "/opt/homebrew/opt/libomp",
        "/usr/local/opt/libomp",
        "/home/linuxbrew/.linuxbrew/opt/libomp",
        "/opt/local",
        "/usr/local",
        "/usr",
    ]
    return [Path(p) for p in out if p]


def find_omp_root() -> Path | None:
    """Only macOS actually needs this. GCC ships libgomp, and MSVC has /openmp
    built in, so on Linux and Windows the compiler already knows."""
    env = os.environ.get("OpenMP_ROOT")
    if env:
        return Path(env)
    if not MACOS:
        return None
    for d in omp_candidates():
        if is_omp_prefix(d):
            return d
    return None


def prefix_candidates() -> list[Path]:
    """Prefixes where CMake might find FFTW, HDF5, ADIOS and friends."""
    out: list[str] = []
    conda = os.environ.get("CONDA_PREFIX")
    if conda:
        # conda-forge puts Windows headers and libs under Library/.
        out.append(os.path.join(conda, "Library") if WINDOWS else conda)
    if WINDOWS:
        vcpkg = os.environ.get("VCPKG_ROOT")
        if vcpkg:
            out.append(os.path.join(vcpkg, "installed", "x64-windows"))
    else:
        brew = os.environ.get("HOMEBREW_PREFIX")
        if brew:
            out.append(brew)
        if have("brew"):
            try:
                out.append(
                    subprocess.run(
                        ["brew", "--prefix"],
                        capture_output=True, text=True, check=True,
                    ).stdout.strip()
                )
            except (subprocess.CalledProcessError, OSError):
                pass
        out += [
            "/opt/homebrew",
            "/home/linuxbrew/.linuxbrew",
            "/opt/local",
            "/usr/local",
        ]
    return [Path(p) for p in out if p]


def find_prefix_path() -> Path | None:
    for d in prefix_candidates():
        if (d / "lib").is_dir() and (d / "include").is_dir():
            return d
    return None


def resolve_compute(flag: str | None) -> str:
    """Pick a compute backend: flag, then $WARPX_COMPUTE, then the machine."""
    if flag and flag != "auto":
        return flag.upper()
    env = os.environ.get("WARPX_COMPUTE")
    if env:
        return env.upper()

    if cuda_available():
        return "CUDA"
    if WINDOWS:
        # WarpX's own Windows CI builds NOACC, and MSVC only implements
        # OpenMP 2.0. Serial is slower but it actually works; opt in to OMP
        # explicitly if your Windows toolchain handles it.
        return "NOACC"
    if MACOS:
        # AppleClang ships no OpenMP runtime, so OMP is only real if libomp is.
        return "OMP" if find_omp_root() else "NOACC"
    return "OMP"


def resolve_fft(flag: str | None, compute: str) -> str:
    """FFT needs FFTW on CPU; CUDA brings its own cuFFT."""
    if flag:
        return "ON" if flag.upper() in ("1", "ON", "TRUE", "YES") else "OFF"
    env = os.environ.get("WARPX_FFT")
    if env:
        return "ON" if env.upper() in ("1", "ON", "TRUE", "YES") else "OFF"
    if compute == "CUDA":
        return "ON"
    if WINDOWS:
        # FFTW is rarely present on a stock Windows toolchain; a missing
        # dependency should not be the first thing a new checkout hits.
        return "OFF"
    return "ON"


def resolve_mpi(flag: str | None) -> str:
    if flag:
        return "ON" if flag.upper() in ("1", "ON", "TRUE", "YES") else "OFF"
    env = os.environ.get("WARPX_MPI")
    if env:
        return "ON" if env.upper() in ("1", "ON", "TRUE", "YES") else "OFF"
    return "OFF"


def build_env(compute: str, mpi: str) -> dict:
    """Environment for pywarpx's setup.py, which is configured by env var."""
    env = os.environ.copy()
    env["WARPX_COMPUTE"] = compute
    env["WARPX_MPI"] = mpi
    env["WARPX_DIMS"] = DIMS

    omp = find_omp_root()
    if omp:
        env["OpenMP_ROOT"] = str(omp)

    prefix = find_prefix_path()
    if prefix:
        existing = env.get("CMAKE_PREFIX_PATH", "")
        parts = existing.split(os.pathsep) if existing else []
        if str(prefix) not in parts:
            env["CMAKE_PREFIX_PATH"] = os.pathsep.join([str(prefix)] + parts)
    return env


def ncpu() -> int:
    return os.cpu_count() or 4


# -------------------------------------------------------------------- running


def run(cmd: list[str], env: dict | None = None, cwd: Path | None = None) -> None:
    print("+ " + " ".join(str(c) for c in cmd), file=sys.stderr)
    try:
        proc = subprocess.run([str(c) for c in cmd], env=env, cwd=cwd)
    except FileNotFoundError:
        raise WarpXError(f"{cmd[0]}: not found on PATH")
    if proc.returncode != 0:
        raise SystemExit(proc.returncode)


def solver_path(root: Path) -> Path:
    """Find the built solver.

    CMake symlinks the fully-qualified binary to a short `warpx.3d`, but that
    step needs privileges Windows does not grant by default, so fall back to
    globbing and taking the newest match.
    """
    bindir = root / "vendor" / "warpx" / "build" / "bin"
    for name in (f"warpx.{DIM_SUFFIX}.exe", f"warpx.{DIM_SUFFIX}"):
        p = bindir / name
        if p.is_file():
            return p
    matches = [
        p for p in bindir.glob(f"warpx.{DIM_SUFFIX}.*")
        if p.is_file() and os.access(p, os.X_OK)
    ]
    if matches:
        return max(matches, key=lambda p: p.stat().st_mtime)
    raise WarpXError("solver not built — run: warpx build")


# ------------------------------------------------------------------- commands


def cmd_root(args) -> None:
    print(find_root(args.project))


def cmd_info(args) -> None:
    root = find_root(args.project)
    compute = resolve_compute(args.compute)
    mpi = resolve_mpi(args.mpi)
    fft = resolve_fft(args.fft, compute)
    omp = find_omp_root()
    prefix = find_prefix_path()

    try:
        solver: object = solver_path(root)
    except WarpXError as exc:
        solver = f"({exc})"

    rows = [
        ("project root", root),
        ("platform", f"{platform.system()} {platform.machine()}"),
        ("python", f"{sys.version.split()[0]} ({sys.executable})"),
        ("submodule", "present" if (root / "vendor/warpx/CMakeLists.txt").is_file()
                      else "MISSING — git submodule update --init --recursive"),
        ("compute", compute + ("" if args.compute or os.environ.get("WARPX_COMPUTE")
                               else "  (autodetected)")),
        ("cuda compiler", "yes" if cuda_available() else "no"),
        ("dims", DIMS),
        ("mpi", mpi),
        ("fft", fft),
        ("OpenMP_ROOT", omp or "(not needed / not found)"),
        ("CMAKE_PREFIX_PATH", prefix or "(none prepended)"),
        ("generator", "Ninja" if have("ninja") else "(CMake default)"),
        ("parallelism", ncpu()),
        ("solver", solver),
    ]
    width = max(len(k) for k, _ in rows)
    for key, value in rows:
        print(f"{key.rjust(width)}  {value}")


def cmd_sync(args) -> None:
    root = find_root(args.project)
    require_submodule(root)
    compute = resolve_compute(args.compute)
    mpi = resolve_mpi(args.mpi)
    if not have("uv"):
        raise WarpXError("uv not found on PATH — see https://docs.astral.sh/uv/")
    print(f"pywarpx: compute={compute} mpi={mpi} dims={DIMS}", file=sys.stderr)
    run(["uv", "sync", "--project", str(root)] + args.rest,
        env=build_env(compute, mpi))


def cmd_rebuild(args) -> None:
    args.rest = ["--reinstall-package", "pywarpx"] + args.rest
    cmd_sync(args)


def cmd_build(args) -> None:
    root = find_root(args.project)
    require_submodule(root)
    compute = resolve_compute(args.compute)
    mpi = resolve_mpi(args.mpi)
    fft = resolve_fft(args.fft, compute)

    if not have("cmake"):
        raise WarpXError("cmake not found on PATH")

    src = root / "vendor" / "warpx"
    build = src / "build"
    print(f"warpx.{DIM_SUFFIX}: compute={compute} mpi={mpi} fft={fft}", file=sys.stderr)

    configure = ["cmake", "-S", str(src), "-B", str(build)]
    if have("ninja"):
        configure += ["-G", "Ninja"]
    configure += [
        f"-DWarpX_MPI={mpi}",
        f"-DWarpX_COMPUTE={compute}",
        f"-DWarpX_DIMS={DIMS}",
        f"-DWarpX_FFT={fft}",
        "-DCMAKE_BUILD_TYPE=Release",
    ] + args.rest
    run(configure)

    # --config is required by multi-config generators (Visual Studio, Xcode)
    # and ignored by single-config ones, so it is safe to always pass.
    run(["cmake", "--build", str(build),
         "--parallel", str(ncpu()), "--config", "Release"])


def cmd_run(args) -> None:
    root = find_root(args.project)
    run([str(solver_path(root))] + args.rest)


# ----------------------------------------------------------------------- main


COMMANDS = {
    "info": (cmd_info, "show detected configuration"),
    "root": (cmd_root, "print the project root"),
    "sync": (cmd_sync, "build pywarpx into .venv"),
    "rebuild": (cmd_rebuild, "force a pywarpx rebuild"),
    "build": (cmd_build, "build the standalone solver"),
    "run": (cmd_run, "run the standalone solver"),
}

GLOBAL_FLAGS = ("--project", "--compute", "--mpi", "--fft")


def split_argv(argv: list[str]) -> tuple[list[str], list[str]]:
    """Split at the subcommand: warpx's own flags before it, pass-through after.

    argparse.REMAINDER cannot do this — it refuses to start collecting on an
    option-like token, so `warpx sync --reinstall` would be rejected instead of
    forwarded. Splitting by hand also keeps the order of the forwarded
    arguments exactly as typed.
    """
    for i, tok in enumerate(argv):
        if tok in COMMANDS:
            head, tail = argv[: i + 1], argv[i + 1 :]
            # `warpx build --compute cuda` would otherwise hand --compute to
            # cmake, which quietly ignores it and builds the wrong backend.
            stray = [t for t in tail if t.split("=")[0] in GLOBAL_FLAGS]
            if stray:
                raise WarpXError(
                    f"{stray[0]} is a warpx option, not a pass-through argument — "
                    f"put it before '{tok}': warpx {stray[0]} ... {tok} ..."
                )
            return head, tail
    return argv, []


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    parser = argparse.ArgumentParser(
        prog="warpx",
        description="Build and run WarpX for this project.",
        epilog="Anything after the subcommand is passed straight through to uv, "
               "cmake or the solver, so warpx's own options go before it: "
               "`warpx --compute cuda build -DAMReX_CUDA_ARCH=8.6`.",
    )
    parser.add_argument("--project", metavar="DIR",
                        help="project root (default: autodetect)")
    parser.add_argument("--compute", choices=COMPUTE_CHOICES,
                        help="compute backend (default: autodetect)")
    parser.add_argument("--mpi", metavar="ON|OFF", help="build with MPI (default: OFF)")
    parser.add_argument("--fft", metavar="ON|OFF",
                        help="build with FFT support (default: autodetect)")

    subs = parser.add_subparsers(dest="command", required=True)
    for name, (fn, helptext) in COMMANDS.items():
        subs.add_parser(name, help=helptext, add_help=False).set_defaults(func=fn)

    try:
        head, tail = split_argv(argv)
        args = parser.parse_args(head)
        args.rest = tail
        args.func(args)
    except WarpXError as exc:
        print(f"warpx: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
