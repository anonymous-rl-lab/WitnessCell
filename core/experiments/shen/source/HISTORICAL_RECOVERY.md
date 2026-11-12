# Historical source recovery

Official repository: https://github.com/ucsd-ccbb/mali-dual-crispr-pipeline

Data commit: `7c05a5a551dc5cd3a3c841d2f22215a450b353a3` (2017-03-06).
Historical count paths were removed in the package reorganization commit `c646ecb7899e106bdbcff46a2f005a79b796d31c` (2017-03-31).

Read the four exact download URLs, SHA256 and Git blob identifiers in `../measurement/SOURCE_PROVENANCE.json`. These are original bytes recovered from public Git history, not synthetic examples or scores converted back into counts.

The discovery inventory `historical_counts_receipts.json` lists ten source candidates. Only the four CV4 files in SOURCE_PROVENANCE.json were used and included here. Nonmatching experimental variants were not pooled into the assay.

The published repository also contains an A549 four-timepoint raw count fixture. Its eight selected sample columns agree exactly for all 24,309 common constructs. That fixture omits dual-NTC controls, which is why the historical full count files were needed.

The author demonstration notebook uses days 3,14,20,28 despite d21 count headers. WitnessCell follows the paper, headers and pre-existing protocol: 3,14,21,28. This is documented before model fitting in SOURCE_AND_FIT_LOCK.json.

ORIGINAL_PIPELINE_LICENSE.txt is copied verbatim from the author-provided ZIP. It is retained for that source software and is not a new licence for all data or WitnessCell modules.
