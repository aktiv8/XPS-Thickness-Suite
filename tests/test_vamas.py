"""
VAMAS reader/writer tests. The synthetic round trip always runs; the test against the real K-Alpha+ tilt series runs
only when that folder exists on this machine (those files are not in the repository).

    python -m unittest discover -s tests -v
"""
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import synthetic as S  # noqa: E402
import vamas as VM  # noqa: E402

REAL = r"C:\Users\david\Documents\Working Directory\20250626\Dave Morgan\X-Ray014 400um - FG  ON"
REAL_SI = os.path.join(REAL, "Stage Tilt ARXPS Si Wafer", "i tilt.vms")
REAL_GAAS = os.path.join(REAL, "Stage Tilt ARXPS GaAs Wafer", "gaas titl.vms")


def _text(path):
    with open(path, encoding="latin-1", newline="") as fh:     # newline="" keeps the CRLF line ends the writer used
        return fh.read()


def _roundtrip(blocks, **kw):
    fd, path = tempfile.mkstemp(suffix=".vms")
    os.close(fd)
    try:
        VM.write(path, blocks, **kw)
        return VM.read(path)
    finally:
        os.remove(path)


class RoundTrip(unittest.TestCase):
    def setUp(self):
        self.blocks = S.two_state_series(azimuths=(0.0, 90.0), noise=True)
        self.file = _roundtrip(self.blocks)

    def test_blocks_energy_axis_and_counts_survive(self):
        self.assertEqual(len(self.file.blocks), len(self.blocks))
        for a, b in zip(self.blocks, self.file.blocks):
            self.assertEqual(a.block_id, b.block_id)
            self.assertTrue(np.allclose(a.be_ev, b.be_ev, atol=1e-6))
            self.assertTrue(np.allclose(a.counts, b.counts))
            self.assertTrue(np.allclose(a.transmission, b.transmission))
            self.assertEqual(b.scans, 4)
            self.assertAlmostEqual(b.dwell_s, 0.05)
            self.assertAlmostEqual(b.hv_ev, 1486.68)

    def test_binding_energy_is_ascending_and_is_hv_minus_kinetic(self):
        be = self.file.blocks[0].be_ev
        self.assertTrue(np.all(np.diff(be) > 0))
        self.assertAlmostEqual(be[0], 90.0, places=6)
        self.assertAlmostEqual(be[-1], 115.0, places=6)

    def test_regions_and_series_group_by_azimuth_and_sort_by_angle(self):
        self.assertEqual(self.file.regions(), ["Si2p Scan"])
        series = self.file.series("Si2p Scan", "exp1", "exp2")
        self.assertEqual(sorted(series), [0.0, 90.0])
        for az, rows in series.items():
            self.assertEqual([a for a, _ in rows], list(S.DEFAULT_TILTS))

    def test_default_sources_pick_tilt_then_azimuth(self):
        self.assertEqual(self.file.default_sources("Si2p Scan"), ("exp1", "exp2"))

    def test_unspecified_values_are_none(self):
        f = _roundtrip(S.two_state_series())                     # azimuth None -> 1e37 in the file
        self.assertIsNone(f.blocks[0].exp_values[1])
        self.assertEqual(f.sources_present("Si2p Scan"), ["exp1"])
        self.assertEqual(f.default_sources("Si2p Scan"), ("exp1", None))
        self.assertEqual(list(f.series("Si2p Scan", "exp1", None)), [None])

    def test_value_count_is_the_total_across_interleaved_columns(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "a.vms")
            VM.write(path, self.blocks[:1])
            text = _text(path).splitlines()
            n = len(self.blocks[0].be_ev)
            self.assertIn(str(2 * n), text)                      # Counts + Transmission


class Errors(unittest.TestCase):
    def _parse(self, text):
        return VM.parse(text)

    def test_not_a_vamas_file(self):
        with self.assertRaises(VM.VamasError):
            self._parse("hello\nworld\n")

    def test_truncated_file(self):
        fd, path = tempfile.mkstemp(suffix=".vms")
        os.close(fd)
        try:
            VM.write(path, S.two_state_series()[:2])
            lines = _text(path).splitlines()
            with self.assertRaises(VM.VamasError):
                self._parse("\n".join(lines[:len(lines) // 2]))
        finally:
            os.remove(path)

    def test_map_experiment_mode_is_refused(self):
        fd, path = tempfile.mkstemp(suffix=".vms")
        os.close(fd)
        try:
            VM.write(path, S.two_state_series()[:1])
            text = _text(path).replace("\r\nNORM\r\n", "\r\nMAP\r\n", 1)
            with self.assertRaisesRegex(VM.VamasError, "MAP"):
                self._parse(text)
        finally:
            os.remove(path)

    def test_missing_file(self):
        with self.assertRaises(VM.VamasError):
            VM.read(os.path.join(tempfile.gettempdir(), "no-such-file.vms"))


@unittest.skipUnless(os.path.isfile(REAL_SI), "the real K-Alpha+ tilt series is not on this machine")
class RealKAlphaFiles(unittest.TestCase):
    def test_si_wafer_series(self):
        f = VM.read(REAL_SI)
        self.assertEqual(len(f.blocks), 30)
        self.assertEqual(f.regions(), ["C1s Scan", "O1s Scan", "Si2p Scan", "Si2s Scan", "XPS Survey"])
        self.assertEqual(f.default_sources("Si2p Scan"), ("exp1", None))
        self.assertEqual([a for a, _ in f.series("Si2p Scan", "exp1")[None]], [0, 15, 30, 45, 60, 75])
        b = f.blocks_of("Si2p Scan")[0]
        self.assertEqual(len(b.be_ev), 251)                      # 502 values, two interleaved columns
        self.assertAlmostEqual(b.be_ev[0], 90.08, places=2)
        self.assertAlmostEqual(b.be_ev[-1], 115.08, places=2)
        self.assertIsNotNone(b.transmission)

    @unittest.skipUnless(os.path.isfile(REAL_GAAS), "GaAs file not present")
    def test_gaas_series(self):
        f = VM.read(REAL_GAAS)
        self.assertEqual(len(f.blocks), 36)
        self.assertEqual(f.regions(), ["As2p3 Scan", "C1s Scan", "Ga2p Scan", "Ga3d As3d Scan", "O1s Scan", "XPS Survey"])
        self.assertEqual(f.blocks_of("Ga3d As3d Scan")[0].scans, 6)


if __name__ == "__main__":
    unittest.main()
