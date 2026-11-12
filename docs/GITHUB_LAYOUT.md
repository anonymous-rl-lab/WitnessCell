# GitHub layout and publication

The archive expands to one repository root, `witnesscell/`. It can be committed directly after rights-holder review; no Git history, credentials or external repository URL is included.

- `README.md`, `README_CN.md`, pinned requirements, `.python-version`, `CITATION.cff` and licence-scope documents are at the root.
- `.github/workflows/reproduce.yml` provides smoke checks on pushes/PRs and a manual full replay. Hosted execution has not been performed by this packaging task.
- `tests/` tests output-path safety, task declarations and controller algebra/leakage boundaries through existing interface checks.
- Scientific inputs remain local and lossless. Packed payloads and an alias map avoid committing multiple byte-identical copies. The preparation command reconstructs the logical layout and generated paths are ignored.
- All physical files are below GitHub's 100 MiB hard limit; the release audit also checks the 50 MiB warning threshold. Use Git via CLI for initial publication, not uploading the archive as one repository source file.

Reference documentation consulted for packaging: https://docs.github.com/en/actions/tutorials/build-and-test-code/python and https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github .

A GitHub-ready layout is separate from permission to publish all scientific modules. LICENSE.md records the remaining author decision without inventing a licence.
