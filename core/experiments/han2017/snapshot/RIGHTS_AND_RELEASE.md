# Code access and release status

The complete original scientific archive is supplied inside paper package v34 as `reproduction/grn_repo_v27.zip`. It is byte-identical to the original scientific-code release v27; its size and SHA256 are recorded in `reproduction/ARCHIVE_PROVENANCE.json`. Paper and original-code version numbers refer to different objects. The runnable transcriptomic and functional adapters are in `software/`.

The MIT file in `software/LICENSE-ADAPTER.txt` applies only to the explicitly named new adapters and tests. The original scientific modules, archived upstream material and dependencies retain their existing rights. No unified licence file was present in the supplied original archive. Including source code for review does not establish a new open-source licence for that archive.

Before external submission/release, the authors should:

1. Approve a licence covering the scientific modules they own, identify third-party exceptions and preserve their notices.
2. Supply the complete reviewer archive or an accessible reviewer repository, and test access from outside the authors' account.
3. Establish the public persistent release and record the actual repository URL, release tag and archive DOI when available.
4. Ensure external raw-data and large-weight access instructions remain valid.

No public URL, DOI, anonymous external access test or author-approved upstream licence is asserted by this package. These release actions were not performed during manuscript preparation. No scientific results depend on a future release identifier.

## Component-level scope

| Component | Current status | Required handling |
| --- | --- | --- |
| Named adapters, example generator and interface tests | MIT scope explicitly listed in `software/LICENSE-ADAPTER.txt` | Preserve that notice and its stated file scope |
| Seven RNA and three functional native modules | Unchanged originals, hashes retained | Authors/rightsholders determine the applicable use licence; adapter MIT does not extend to them |
| `reproduction/grn_repo_v27.zip` | Complete original archive with historical records and third-party material | Preserve original notices and review ownership before any unified release licence |
| Revision analysis and figure/document scripts | Included for review and reproducibility | Include in the author-owned licence decision where ownership permits |
| Public data, annotations, dependencies and external weights | Original accessions, provenance and asset indexes retained | Follow source-specific terms and do not relicense third-party assets |

Submission-stage code can be provided as a ZIP or accessible link. A public DOI is not asserted or fabricated to fill an unavailable field. Publication-stage persistent access and an actual version identifier remain author release responsibilities. `AUTHOR_ACTIONS.md` lists the outstanding author declarations separately from completed computational checks.

The new `partner_control/` scientific code, interface and analyses are included for review and reproduction. Their public licence is part of the author-owned release decision above; the older adapter-only MIT file does not automatically extend to them. The Shen archive preserves its source notices.
