#!/bin/bash
# HANDOFF steps a and b of PLAN.md. Restartable (finished trials are skipped). Commits every 5 min and per stage.
#   nohup bash tro/results/continuous_baselines/run_handoff.sh > /tmp/handoff.log 2>&1 &
cd "$(dirname "$0")/../.."            # -> tro/
commit() {
  (cd .. && git add tro/results && git commit -q -m "$1

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TByxUTGsR8uPvFhDBetfgi" || true
   for d in 2 4 8 16; do git push -q origin claude/tro-2027-paper-review-wg3qzv && break; sleep $d; done)
}
( while sleep 300; do commit "Handoff runs: partial results (in progress)"; done ) &
LOOP=$!
# a: isolated timing of ours (nothing else running)
python prototype_3d/patch_size_timing.py --fields 12mm 6mm && commit "Isolated timing of our filter, 12 mm and 6 mm (handoff a)"
# b: ours 6 mm without sampled data (rows at the samples, like the baselines), outcomes
python prototype_3d/certificate_upgrades.py --variant free --field 6mm --qp daqp --tag fixed --trials 0 15 &
python prototype_3d/certificate_upgrades.py --variant free --field 6mm --qp daqp --tag fixed --trials 15 30 &
wait %2 %3 2>/dev/null; wait $(jobs -p | grep -v "^$LOOP$") 2>/dev/null
kill $LOOP
python prototype_3d/continuous_baselines_report.py > /dev/null
commit "Ours 6 mm without sampled data (handoff b) and updated comparison table"
echo HANDOFF_DONE
