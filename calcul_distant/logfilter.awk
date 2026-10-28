# Log filter for porousGasificationFoamMR runs (2026-09-25).
# Usage: mpirun ... porousGasificationFoamMR -parallel 2>&1 | awk -v N=20 -f logfilter.awk > log
#
# The raw solver output is ~4.3 kB per time step (~100 GB for 400 simulated
# minutes). This keeps:
#   - everything before the first "Time =" (case setup, model selection);
#   - the full block of every N-th time step (default N=20);
#   - on EVERY step, any line signalling a problem (FOAM Warning/FATAL,
#     janafThermo, SIGFPE/Floating point, nan, error), so no alert is lost.
# The monitoring greps ("^Time =", "^deltaT", "T gas min/max",
# "pyrolysis min/max(T)", "laminarEDM:", "janafThermo") keep working.
BEGIN { if (N == "") N = 20; step = -1 }
/^Time = / { step++ }
{
    if (step < 0 || step % N == 0) {
        if ($0 ~ /^Time = / && nwarn > 20000) print "logfilter: warning lines so far = " nwarn
        if ($0 ~ /attempt to use|From function|FOAM Warning|in file/ && nwarn > 20000) { nwarn++; next }
        print; fflush(); next
    }
    if ($0 ~ /FOAM Warning|FOAM FATAL|From function|janafThermo|Floating point|sigFpe|[Nn]a[Nn]|rror|DBG/) {
        # cap (2026-09-25): a diverging run printed >20 GB of janafThermo
        # warnings; keep the first 20000 warning lines, then only a counter
        # on each printed step. FATAL/Floating point are always printed.
        if ($0 ~ /FOAM FATAL|Floating point|sigFpe|DBG/) { print; fflush(); next }
        nwarn++
        if (nwarn <= 20000) { print; fflush() }
    }
}
