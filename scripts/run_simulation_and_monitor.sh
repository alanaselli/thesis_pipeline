#!/bin/bash

# Halt on error.
set -euo pipefail

# Run monitor_memory in the background
scripts/monitor_memory.sh  &

# Capture the background process ID (PID)
script1_pid=$!

# Configure for BLUPF90
ulimit -s unlimited
export OMP_STACKSIZE=64M

# Run simulation
/usr/local/bin/R --no-restore --file=scripts/main.R --args -g '10' -m '3 2' -f '20 15 10 5' -c '50' -n '10_50' > scripts/10_50.log
/usr/local/bin/R --no-restore --file=scripts/main.R --args -g '50' -m '3 2' -f '20 15 10 5' -c '50' -n '50_50' > scripts/50_50.log
/usr/local/bin/R --no-restore --file=scripts/main.R --args -g '10' -m '6 4' -f '40 30 20 10' -c '100' -n '10_100' > scripts/10_100.log
/usr/local/bin/R --no-restore --file=scripts/main.R --args -g '50' -m '6 4' -f '40 30 20 10' -c '100' -n '50_100' > scripts/50_100.log
/usr/local/bin/R --no-restore --file=scripts/main.R --args -g '10' -m '6 4' -f '80 60 40 20' -c '200' -n '10_200' > scripts/10_200.log
/usr/local/bin/R --no-restore --file=scripts/main.R --args -g '50' -m '6 4' -f '80 60 40 20' -c '200' -n '50_200' > scripts/50_200.log

# Kill the background process of script1
kill $script1_pid
