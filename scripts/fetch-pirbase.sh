#!/usr/bin/env bash
# Mirror piRBase v3.0 into a directory the pipeline can read instead of the
# site. Their server drops any plain GET of a file above ~100MB but still
# answers bounded range requests, so every file is fetched in 50MB pieces.
# Rerunnable: finished files are skipped, partial ones resume.
set -euo pipefail

base='http://bigdata.ibp.ac.cn/piRBase/download/v3.0'
out=${1:?usage: fetch-pirbase.sh <output dir>}
chunk=$((50 * 1024 * 1024))
mkdir -p "$out/fasta"

"$(dirname "$0")/../bin/rnac" pirbase urls "$base" | cut -d, -f3 | while read -r url; do
  file="$out/fasta/$(basename "$url")"
  total=$(curl -s -m 60 -r 0-0 -D - -o /dev/null "$url" | tr -d '\r' | sed -n 's/^Content-Range: bytes 0-0\///p')
  if [ -z "$total" ]; then
    echo "SKIP $(basename "$url"): no response from server" >&2
    continue
  fi
  have=$(stat -c %s "$file" 2>/dev/null || stat -f %z "$file" 2>/dev/null || echo 0)
  while [ "$have" -lt "$total" ]; do
    end=$((have + chunk - 1)); [ "$end" -ge "$total" ] && end=$((total - 1))
    curl -sf -m 1800 -r "$have-$end" "$url" >> "$file" || { echo "retry $(basename "$url") at $have" >&2; sleep 30; }
    have=$(stat -c %s "$file" 2>/dev/null || stat -f %z "$file")
  done
  gzip -t "$file" && echo "OK   $(basename "$url") $total bytes"
done
