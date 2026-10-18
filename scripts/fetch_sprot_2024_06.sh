#!/usr/bin/env bash
# UniProtKB/Swiss-Prot release 2024_06 (27 Nov 2024): the last release before 2025-01-01.
set -euo pipefail
cd "$(dirname "$0")/.." && mkdir -p data/raw/sprot_2024_06 && cd data/raw/sprot_2024_06
curl -sSL --retry 5 -C - -o sprot.tar.gz \
  "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2024_06/knowledgebase/uniprot_sprot-only2024_06.tar.gz"
tar -xzf sprot.tar.gz uniprot_sprot.dat.gz && rm sprot.tar.gz
