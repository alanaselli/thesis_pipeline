#! /bin/bash

/usr/local/bin/R --no-restore --file=main.R > main.log \
&& /usr/local/bin/R --no-restore --file=convert_SNP_pos.R > SNP_pos.log \
&& /usr/local/bin/R --no-restore --file=plot_LD_decay.R > plot_LD_decay.log \
&& /usr/local/bin/R --no-restore --file=ROH.R > ROH.log \
&& /usr/local/bin/R --no-restore --file=PCA.R > PCA.log