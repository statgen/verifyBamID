# Overlapping-mate regression

After building verifyBamID, run:

```sh
python3 tests/test_overlap_mates.py
```

The test requires Python 3 and `samtools` on `PATH`. Use `--verifybamid` or
`--samtools` to select another executable. It creates temporary SAM/BAM and VCF
fixtures, runs the actual program, and checks the retained read counts. Fixtures
cover low-quality, N and deleted bases in either mate order, valid overlapping
mates, distinct read names, and `--ignoreOverlapPair`.
