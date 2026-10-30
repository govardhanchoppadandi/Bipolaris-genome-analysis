#!/usr/bin/env bash
set -Eeuo pipefail
LINEAGE="pleosporales_odb12"
DOWNLOAD_PATH="${HOME}/busco_downloads"

command -v busco >/dev/null || {
    echo "BUSCO not found. Activate busco_latest first."
    exit 1
}
busco --version
mkdir -p "$DOWNLOAD_PATH"

if [[ ! -d "$DOWNLOAD_PATH/lineages/$LINEAGE" ]]; then
    busco --download "$LINEAGE" --download_path "$DOWNLOAD_PATH"
else
    echo "Already downloaded: $DOWNLOAD_PATH/lineages/$LINEAGE"
fi

test -d "$DOWNLOAD_PATH/lineages/$LINEAGE"
echo "Lineage ready: $LINEAGE"
