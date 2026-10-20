"""
Regenerate the synthetic Si 2p tilt series used by the ARXPS worked examples:
docs/sample/si-native-oxide-arxps.vms (spectra only) and si-native-oxide-arxps-fitted.vms (the same spectra with a
CasaXPS-style peak fit stored in each block), plus the two dual-anode files daxpes-tio2-al-ka.vms and
daxpes-tio2-ag-la.vms (a 3 nm TiO2-like film covering 60 % of the area, with a thick reference, for the DAXPES tab).

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
DAXPES_SOFT = os.path.join(HERE, "sample", "daxpes-tio2-al-ka.vms")
DAXPES_HARD = os.path.join(HERE, "sample", "daxpes-tio2-ag-la.vms")
DAXPES_THICKNESS_NM, DAXPES_COVERAGE = 3.0, 0.6


def main():
    os.makedirs(os.path.dirname(SAMPLE), exist_ok=True)
    synthetic.write_sample(SAMPLE, d_nm=TRUE_THICKNESS_NM)
    synthetic.write_sample(SAMPLE_FITTED, fitted=True, d_nm=TRUE_THICKNESS_NM)
    synthetic.write_daxpes_samples(DAXPES_SOFT, DAXPES_HARD, d_nm=DAXPES_THICKNESS_NM, cov=DAXPES_COVERAGE)
    print("wrote", SAMPLE, "and", SAMPLE_FITTED, "and", DAXPES_SOFT, "and", DAXPES_HARD)


if __name__ == "__main__":
    main()
