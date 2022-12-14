#!/usr/bin/env bash

## Author: Francesco Tabaro francesco.tabaro@embl.it

usage() { echo "Usage: $0 [-n]" 1>&2; exit 1; }

DRY_RUN=0
VERBOSE=0
while getopts ":nhv" o; do
    case "${o}" in
        n)
            DRY_RUN=1
            ;;
        v)
            VERBOSE=1
            ;;
        h)
            usage
            ;;
        *)
    esac
done
shift $((OPTIND-1))

SOURCE_HTML=$(find ../Rmd/ -name '*.html')
TARGET_FOLDER="_posts"
PAGES=$(find $TARGET_FOLDER -type f -printf "%f\n" | sed -r 's/[0-9]{4}-[0-9]{2}-[0-9]{2}-(.*)/\1/')
TODAY=$(date -I)

for F in $SOURCE_HTML; do
  BN=$(basename $F)

  # Try to extract date from HTML - Rmarkdown
  matched_date=$(grep '<meta name="date" content=".*" />' "$F")
  if [ -n "$matched_date" ]; then
      DATE=$(echo "$matched_date" | sed -r 's,<meta name="date" content="(.*)" />,\1,' | sed -r 's/,//' | xargs -I{} date -d "{}" -I)
    else
      DATE=$TODAY
    fi

  # Compute output path
  BN="${BN//_/-}"
  TARGET_FILE="$DATE-$BN"
  TARGET_PATH=$TARGET_FOLDER/$TARGET_FILE

  # Copy HTML to _pages only if output file does not exist yet
  if ! echo "$PAGES" | grep -q "^$BN"; then
      if [ $DRY_RUN -eq 1 ]; then
          echo "[DRY_RUN] $F -> $TARGET_PATH"
        else
          cp -v "$F" "$TARGET_PATH"
        fi
    elif [ $VERBOSE -eq 1 ]; then
        echo "$TARGET_PATH exists. Skipping" >&2
    fi
  done