#!/bin/bash
# Box = patch (6-mm SDF cells) with the vertex certificate and the unit lift (variant unit; 8 rows per box, the
# paper's Proposition coef). Same conditions as every other row of the fair table.
#   bash tro/results/continuous_baselines/run_6mm_vertex.sh
cd "$(dirname "$0")/../.."            # -> tro/
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for r in "0 5" "5 10" "10 15" "15 20" "20 25" "25 30"; do
  SUMMED_REFINE_DEPTH=0 python prototype_3d/certificate_upgrades.py --variant unit --field 6mm --qp daqp --tag fixed --trials $r &
done
wait
n=$(ls results/certificate_upgrades/unit_6mm_fixed/franka_*.json 2>/dev/null | wc -l)
[ "$n" = "30" ] || { echo "OUTCOMES_INCOMPLETE $n"; exit 1; }
n=$(powershell -NoProfile -Command "@(Get-CimInstance Win32_Process | Where-Object { \$_.Name -like 'python*' -and \$_.CommandLine -match 'certificate_upgrades|continuous_baselines|patch_size_timing|refit_native' }).Count")
[ "${n//[^0-9]/}" = "0" ] || { echo "another experiment is running; not timing"; exit 1; }
SUMMED_REFINE_DEPTH=0 python prototype_3d/patch_size_timing.py --variant unit --fields 6mm || { echo "TIMING_FAILED"; exit 1; }
python prototype_3d/continuous_baselines_report.py > /dev/null && echo RUN_6MM_VERTEX_DONE
