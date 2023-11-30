#! /bin/bash

# Halt on error.
set -euo pipefail

/usr/local/bin/R --no-restore --file=scripts/main.R > scripts/main.log