#!/usr/bin/env bash

set -euo pipefail

test -n "${FRESH_PACKET_ROOT:-}"
test -n "${FRESH_TASK_MANIFEST:-}"
packet_root=$(readlink -f -- "${FRESH_PACKET_ROOT}")
test -d "${packet_root}"
(
  cd -- "${packet_root}"
  sha256sum --check --strict packet.sha256
)
manifest_path=$(readlink -f -- "${FRESH_TASK_MANIFEST}")
case "${manifest_path}" in
  "${packet_root}/tasks.small.tsv" | "${packet_root}/tasks.large.tsv") ;;
  *)
    echo "Unfrozen task manifest: ${manifest_path}" >&2
    exit 64
    ;;
esac
line_number=$((SLURM_ARRAY_TASK_ID + 2))
IFS=$'\t' read -r \
  task_id \
  sample_id \
  stack \
  tier \
  attempt_id \
  viralscan_cache_path \
  viralscan_cache_manifest_sha256 \
  read1_path \
  read2_path \
  read1_storage_bytes \
  read2_storage_bytes \
  read1_storage_sha256 \
  read2_storage_sha256 \
  viralscan_path \
  output_path \
  index_path \
  t2g_path \
  whitelist_path \
  cores \
  status_path \
  stdout_path \
  stderr_path \
  < <(sed -n "${line_number}p" "${manifest_path}")
test -n "${task_id}"
test -n "${attempt_id}"
test -n "${viralscan_cache_path}"

cache_manifest="${viralscan_cache_path}/data/manifest.json"
test -f "${cache_manifest}"
observed_cache_sha256=$(sha256sum -- "${cache_manifest}" | cut -d' ' -f1)
if [ "${observed_cache_sha256}" != "${viralscan_cache_manifest_sha256}" ]; then
  echo "Viral-data cache manifest drifted: ${cache_manifest}" >&2
  exit 65
fi
export VIRALSCAN_CACHE="${viralscan_cache_path}"

"${packet_root}/../env_full/bin/python" \
  "${packet_root}/source/run_fresh_control.py" \
  --sample-id "${sample_id}" \
  --stack "${stack}" \
  --viralscan "${viralscan_path}" \
  --output "${output_path}" \
  --read1 "${read1_path}" \
  --read2 "${read2_path}" \
  --read1-storage-bytes "${read1_storage_bytes}" \
  --read2-storage-bytes "${read2_storage_bytes}" \
  --read1-storage-sha256 "${read1_storage_sha256}" \
  --read2-storage-sha256 "${read2_storage_sha256}" \
  --index "${index_path}" \
  --t2g "${t2g_path}" \
  --whitelist "${whitelist_path}" \
  --cores "${cores}" \
  --status "${status_path}" \
  --stdout "${stdout_path}" \
  --stderr "${stderr_path}" \
  --attempt-id "${attempt_id}" \
  --viralscan-cache "${viralscan_cache_path}"
