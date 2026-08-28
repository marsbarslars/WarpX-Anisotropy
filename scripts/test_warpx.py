#!/usr/bin/env python3
"""Tests for the parts of warpx.py that cannot be exercised on one machine.

Backend selection and solver lookup differ per platform, so they are tested by
faking the platform rather than by having three machines. Standard library only,
same as the driver:

    python3 scripts/test_warpx.py
"""

import importlib.util
import os
import pathlib
import sys
import tempfile
import time
import unittest

_spec = importlib.util.spec_from_file_location(
    "warpx_driver", pathlib.Path(__file__).resolve().parent / "warpx.py"
)
warpx = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(warpx)


class PlatformCase(unittest.TestCase):
    """Fakes OS, CUDA and libomp, and clears the env vars that would override."""

    def setUp(self):
        self._saved = {k: os.environ.pop(k, None) for k in
                       ("WARPX_COMPUTE", "WARPX_FFT", "WARPX_MPI", "OpenMP_ROOT")}
        self._orig = (warpx.WINDOWS, warpx.MACOS,
                      warpx.cuda_available, warpx.find_omp_root)

    def tearDown(self):
        for k, v in self._saved.items():
            if v is not None:
                os.environ[k] = v
        (warpx.WINDOWS, warpx.MACOS,
         warpx.cuda_available, warpx.find_omp_root) = self._orig

    def fake(self, osname, cuda=False, libomp=False):
        warpx.WINDOWS = osname == "Windows"
        warpx.MACOS = osname == "Darwin"
        warpx.cuda_available = lambda: cuda
        warpx.find_omp_root = lambda: "/fake/libomp" if libomp else None


class TestComputeSelection(PlatformCase):
    def test_cuda_wins_on_every_os(self):
        for osname in ("Darwin", "Linux", "Windows"):
            self.fake(osname, cuda=True)
            self.assertEqual(warpx.resolve_compute(None), "CUDA", osname)

    def test_linux_defaults_to_omp(self):
        self.fake("Linux")
        self.assertEqual(warpx.resolve_compute(None), "OMP")

    def test_macos_needs_libomp_for_omp(self):
        self.fake("Darwin", libomp=True)
        self.assertEqual(warpx.resolve_compute(None), "OMP")
        self.fake("Darwin", libomp=False)
        self.assertEqual(warpx.resolve_compute(None), "NOACC")

    def test_windows_defaults_to_serial(self):
        # MSVC implements only OpenMP 2.0; WarpX's own Windows CI builds NOACC.
        self.fake("Windows", libomp=True)
        self.assertEqual(warpx.resolve_compute(None), "NOACC")

    def test_every_combination_is_buildable(self):
        valid = set(warpx.COMPUTE_CHOICES) - {"auto"}
        for osname in ("Darwin", "Linux", "Windows"):
            for cuda in (False, True):
                for libomp in (False, True):
                    self.fake(osname, cuda, libomp)
                    self.assertIn(warpx.resolve_compute(None), valid)


class TestPrecedence(PlatformCase):
    def test_flag_beats_env_beats_detection(self):
        self.fake("Linux")                       # would detect OMP
        os.environ["WARPX_COMPUTE"] = "NOACC"
        self.assertEqual(warpx.resolve_compute(None), "NOACC")
        self.assertEqual(warpx.resolve_compute("CUDA"), "CUDA")

    def test_auto_is_not_treated_as_a_backend(self):
        self.fake("Linux")
        self.assertEqual(warpx.resolve_compute("auto"), "OMP")

    def test_flag_is_normalised(self):
        self.fake("Linux")
        self.assertEqual(warpx.resolve_compute("cuda"), "CUDA")


class TestFft(PlatformCase):
    def test_cuda_gets_fft_even_on_windows(self):
        self.fake("Windows", cuda=True)
        self.assertEqual(warpx.resolve_fft(None, "CUDA"), "ON")

    def test_windows_cpu_skips_fft(self):
        self.fake("Windows")
        self.assertEqual(warpx.resolve_fft(None, "NOACC"), "OFF")

    def test_posix_gets_fft(self):
        self.fake("Linux")
        self.assertEqual(warpx.resolve_fft(None, "OMP"), "ON")

    def test_explicit_off(self):
        self.fake("Linux")
        self.assertEqual(warpx.resolve_fft("OFF", "CUDA"), "OFF")


class TestSolverLookup(unittest.TestCase):
    @staticmethod
    def _tree(names):
        root = pathlib.Path(tempfile.mkdtemp())
        bindir = root / "vendor" / "warpx" / "build" / "bin"
        bindir.mkdir(parents=True)
        for i, name in enumerate(names):
            p = bindir / name
            p.write_text("")
            p.chmod(0o755)
            os.utime(p, (time.time() + i, time.time() + i))
        return root

    def test_prefers_the_short_alias(self):
        root = self._tree(["warpx.3d.NOMPI.OMP.DP", "warpx.3d"])
        self.assertEqual(warpx.solver_path(root).name, "warpx.3d")

    def test_windows_exe_alias(self):
        root = self._tree(["warpx.3d.exe"])
        self.assertEqual(warpx.solver_path(root).name, "warpx.3d.exe")

    def test_globs_when_the_symlink_is_missing(self):
        # Windows withholds the privileges CMake's create_symlink needs.
        root = self._tree(["warpx.3d.NOMPI.NOACC.DP.exe"])
        self.assertEqual(warpx.solver_path(root).name, "warpx.3d.NOMPI.NOACC.DP.exe")

    def test_newest_build_wins(self):
        root = self._tree(["warpx.3d.NOMPI.OMP.DP", "warpx.3d.NOMPI.CUDA.DP"])
        self.assertEqual(warpx.solver_path(root).name, "warpx.3d.NOMPI.CUDA.DP")

    def test_reports_when_nothing_is_built(self):
        with self.assertRaises(warpx.WarpXError):
            warpx.solver_path(self._tree([]))


class TestArgvSplit(unittest.TestCase):
    def test_forwards_option_like_arguments(self):
        head, tail = warpx.split_argv(["sync", "--reinstall", "-v"])
        self.assertEqual(head, ["sync"])
        self.assertEqual(tail, ["--reinstall", "-v"])

    def test_keeps_global_flags_in_the_head(self):
        head, tail = warpx.split_argv(["--compute", "CUDA", "build", "-DX=1"])
        self.assertEqual(head, ["--compute", "CUDA", "build"])
        self.assertEqual(tail, ["-DX=1"])

    def test_preserves_argument_order(self):
        _, tail = warpx.split_argv(["run", "inputs.txt", "max_step=5", "-a", "b"])
        self.assertEqual(tail, ["inputs.txt", "max_step=5", "-a", "b"])

    def test_rejects_a_global_flag_after_the_subcommand(self):
        # Otherwise cmake silently ignores it and builds the wrong backend.
        for argv in (["build", "--compute", "CUDA"], ["build", "--compute=CUDA"]):
            with self.assertRaises(warpx.WarpXError):
                warpx.split_argv(argv)

    def test_no_subcommand(self):
        self.assertEqual(warpx.split_argv(["--help"]), (["--help"], []))


class TestOmpDetection(unittest.TestCase):
    def test_requires_header_and_library(self):
        d = pathlib.Path(tempfile.mkdtemp())
        (d / "include").mkdir()
        (d / "lib").mkdir()
        self.assertFalse(warpx.is_omp_prefix(d))     # nothing yet
        (d / "include" / "omp.h").write_text("")
        self.assertFalse(warpx.is_omp_prefix(d))     # header alone is not enough
        (d / "lib" / "libomp.dylib").write_text("")
        self.assertTrue(warpx.is_omp_prefix(d))


if __name__ == "__main__":
    unittest.main(verbosity=2)
