#!/bin/bash
# Box = patch with one parameter (the SDF cell size). 8 mm is the predicted finest size within the 20-Hz budget:
# fitting p95 = a * h^-b to the isolated 12-mm (23.1 ms) and 6-mm (77.8 ms) runs gives p95(8 mm) ~ 47 ms.
#   bash tro/results/continuous_baselines/run_8mm.sh
cd "$(dirname "$0")/../.."            # -> tro/
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
set -o pipefail
python prototype_3d/refit_native_fields.py 8 || { echo "REFIT_FAILED"; exit 1; }
for r in "0 5" "5 10" "10 15" "15 20" "20 25" "25 30"; do
  SUMMED_REFINE_DEPTH=0 python prototype_3d/certificate_upgrades.py --variant free --field 8mm --qp daqp --tag fixed --trials $r &
done
wait
n=$(ls results/certificate_upgrades/free_8mm_fixed/franka_*.json 2>/dev/null | wc -l)
[ "$n" = "30" ] || { echo "OUTCOMES_INCOMPLETE $n"; exit 1; }
(cd prototype_3d && python continuous_baselines_starts.py) || echo "STARTS_FAILED"
# isolated timing: nothing else may run
SUMMED_REFINE_DEPTH=0 python prototype_3d/patch_size_timing.py --variant free --fields 8mm || { echo "TIMING_FAILED"; exit 1; }
python prototype_3d/continuous_baselines_report.py > /dev/null && echo RUN_8MM_DONE
