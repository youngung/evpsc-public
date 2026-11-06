import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GEN_DIM = ROOT / "pyscripts" / "gen_dim.py"


class GenerateDimensionFileTests(unittest.TestCase):
    def run_gen_dim(self, *args):
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "parameters.dim"
            subprocess.run(
                [sys.executable, str(GEN_DIM), "-o", str(output), *args],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            return output.read_text()

    def assert_parameter(self, contents, name, value):
        match = re.search(
            rf"^\s*parameter\({name}=(\d+)\)\s*$",
            contents,
            flags=re.IGNORECASE | re.MULTILINE,
        )
        self.assertIsNotNone(match, f"missing parameter {name}")
        self.assertEqual(int(match.group(1)), value)

    def test_default_dimensions(self):
        contents = self.run_gen_dim()

        for name, value in {
            "nphmx": 1,
            "ngrmx": 10000,
            "nmodmx": 6,
            "ntwmmx": 2,
            "nsysmx": 48,
            "ntwsmx": 24,
            "nneimx": 1,
            "ndifmx": 12,
            "iregime": 2,
        }.items():
            with self.subTest(parameter=name):
                self.assert_parameter(contents, name, value)

    def test_custom_dimensions_and_regime(self):
        contents = self.run_gen_dim(
            "--nphmx", "3",
            "--ngrmx", "2500",
            "--nmodmx", "4",
            "--ntwmmx", "1",
            "--nsysmx", "18",
            "--ntwsmx", "6",
            "--nneimx", "2",
            "--ndifmx", "8",
            "--iregime", "1",
        )

        for name, value in {
            "nphmx": 3,
            "ngrmx": 2500,
            "nmodmx": 4,
            "ntwmmx": 1,
            "nsysmx": 18,
            "ntwsmx": 6,
            "nneimx": 2,
            "ndifmx": 8,
            "iregime": 1,
        }.items():
            with self.subTest(parameter=name):
                self.assert_parameter(contents, name, value)


if __name__ == "__main__":
    unittest.main()
