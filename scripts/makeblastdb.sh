#!/usr/bin/env bash
# Build a nucleotide BLAST database in a temporary directory, then move it to
# OUT_PREFIX.
#
# Usage:
#   bash makeblastdb.sh IN.fna OUT_PREFIX
#
# makeblastdb re-opens the new database by absolute path when it finishes and
# exits with an error if that path contains spaces (e.g. a Dropbox folder),
# even though the files are fine. $TMPDIR has no spaces, and blastn itself
# handles relative -db paths with spaces in the working directory.

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: bash makeblastdb.sh IN.fna OUT_PREFIX" >&2
    exit 1
fi

IN=$1
OUT=$2
NAME=$(basename "$OUT")

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

makeblastdb -in "$IN" -dbtype nucl -out "$tmp/$NAME"
mkdir -p "$(dirname "$OUT")"
mv "$tmp/$NAME".* "$(dirname "$OUT")/"
