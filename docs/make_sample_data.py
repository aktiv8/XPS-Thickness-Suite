"""
Regenerate docs/sample/si-native-oxide-arxps.vms, the synthetic Si 2p tilt series used by the ARXPS worked example.

    python docs/make_sample_data.py

The file is made by synthetic.py (a 1.5 nm SiO2 layer on Si, charging drift, Poisson noise, seeded), so it is
reproducible and contains no measurement. tests/test_docs.py checks that the committed file is what this script writes.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import synthetic  # noqa: E402

SAMPLE = os.path.join(HERE, "sample", "si-native-oxide-arxps.vms")
TRUE_THICKNESS_NM = 1.5


def main():
    os.makedirs(os.path.dirname(SAMPLE), exist_ok=True)
    synthetic.write_sample(SAMPLE, d_nm=TRUE_THICKNESS_NM)
    print("wrote", SAMPLE)


if __name__ == "__main__":
    main()
