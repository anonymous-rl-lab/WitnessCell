# Anonymous GitHub upload checklist

## Before upload

Run:

```bash
bash scripts/run_smoke.sh
python scripts/audit_repository.py
```

Both commands must pass. Do not upload the original 430 MiB evidence archive;
it contains redundant assets and machine-specific historical metadata.

## Create the repository

Use a neutral account and repository name during double-blind review. Extract
the delivered ZIP, enter its root, then run:

```bash
git init
git add .
git commit -m "Anonymous WitnessCell reproducibility release v18"
git branch -M main
git remote add origin <anonymous-repository-url>
git push -u origin main
```

No Git LFS configuration is required for this compact release; the repository
audit rejects any file above 50 MiB.

## Repository settings

- keep Actions enabled so both CI workflows run;
- enable branch protection after the first successful run;
- do not add author names, affiliations, personal account links, ORCID, email,
  project DOI, funding text or preregistration-account metadata;
- do not rewrite frozen `NOT_ADJUDICATED` decisions;
- do not commit downloaded h5ad files or generated prediction directories.

After unblinding, author metadata and the canonical repository/DOI may be added
in a separate release commit.
