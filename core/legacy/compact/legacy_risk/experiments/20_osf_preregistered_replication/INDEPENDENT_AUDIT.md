# Independent audit

The returned AutoDL archive was independently inspected after upload.

## Integrity checks

- The gzip/tar archive opened without error.
- The result archive SHA-256 was `511860f645158a9af536a3944e03d4cc440674e23b55d826bf615e02d96bea0c`.
- The recorded source hash was the frozen `955d446a0cfe65a6f4f191d74c7a7d16a6083488b9beb40a5e155ab94d528bd7`.
- The result copy of `FROZEN_PROTOCOL.json` matched the registered-package copy with SHA-256 `fa54095c35fa369b4c0519439df5dabb4d722b6eab7624aaa18349e48693fb21`.
- The run started at `2026-08-10T10:56:15Z` and completed at `2026-08-10T11:13:03Z`.
- No traceback, exception, numerical failure, non-finite output or repair attempt was found in the returned logs.

## Independent numerical comparison

The 12 regenerated `deploy_predictions.npz` archives were compared directly against the frozen archives in the source package with `allow_pickle=False`. For every array, the key, shape, dtype and value were checked independently.

- Archives: 12
- Arrays per archive: 13
- Total comparisons: 156
- Exact comparisons: 156
- Global maximum absolute difference: `0.0`
- Paired condition–seed rows: 654

The independently recomputed paired effects were identical to the registered values:

- all-gene MSE gain: `0.001079217249737252`;
- top-100 MSE gain: `0.0030453806218846413`;
- top-100 Pearson delta: `0.0013006984799760163`.

The audit therefore confirms the packaged status `PASS_PREREGISTERED_COMPUTATIONAL_REPLICATION` within the registered computational scope.

