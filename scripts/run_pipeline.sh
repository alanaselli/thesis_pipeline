#! /bin/bash

# Halt on error.
set -euo pipefail

/usr/local/bin/R --no-restore --file=scripts/main.R > scripts/main.log \
&& /usr/local/bin/R --no-restore --file=scripts/convert_SNP_pos.R > scripts/SNP_pos.log \
&& /usr/local/bin/R --no-restore --file=scripts/plot_LD_decay.R > scripts/plot_LD_decay.log \
&& /usr/local/bin/R --no-restore --file=scripts/ROH.R > scripts/ROH.log \
&& /usr/local/bin/R --no-restore --file=scripts/PCA.R > scripts/PCA.log