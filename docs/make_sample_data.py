"""
Regenerate the synthetic Si 2p tilt series used by the ARXPS worked examples:
docs/sample/si-native-oxide-arxps.vms (spectra only) and si-native-oxide-arxps-fitted.vms (the same spectra with a
CasaXPS-style peak fit stored in each block).

    python docs/make_sample_data.py

The files are made by synthetic.py (a 1.5 nm SiO2 layer on Si, charging drift, Poisson noise, seeded), so they are
reproducible and contain no measurement. tests/test_docs.py checks that the committed files are what this script writes.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import synthetic  # noqa: E402

SAMPLE = os.path.join(HERE, "sample", "si-native-oxide-arxps.vms")
SAMPLE_FITTED = os.path.join(HERE, "sample", "si-native-oxide-arxps-fitted.vms")
TRUE_THICKNESS_NM = 1.5


def main():
    os.makedirs(os.path.dirname(SAMPLE), exist_ok=True)
    synthetic.write_sample(SAMPLE, d_nm=TRUE_THICKNESS_NM)
    synthetic.write_sample(SAMPLE_FITTED, fitted=True, d_nm=TRUE_THICKNESS_NM)
    print("wrote", SAMPLE, "and", SAMPLE_FITTED)


if __name__ == "__main__":
    main()
