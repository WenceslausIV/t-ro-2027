#!/bin/bash
# Box = patch with the vertex certificate (variant unit) at a refit SDF cell size: refit, outcomes, start validity,
# isolated timing, report. Same conditions as every other row of the fair table.
#   bash tro/results/continuous_baselines/run_vertex_size.sh 5
mm=$1
cd "$(dirname "$0")/../.."            # -> tro/
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python prototype_3d/refit_native_fields.py $mm || { echo "REFIT_FAILED"; exit 1; }
for r in "0 5" "5 10" "10 15" "15 20" "20 25" "25 30"; do
  SUMMED_REFINE_DEPTH=0 python prototype_3d/certificate_upgrades.py --variant unit --field ${mm}mm --qp daqp --tag fixed --trials $r &
done
wait
n=$(ls results/certificate_upgrades/unit_${mm}mm_fixed/franka_*.json 2>/dev/null | wc -l)
[ "$n" = "30" ] || { echo "OUTCOMES_INCOMPLETE $n"; exit 1; }
(cd prototype_3d && python continuous_baselines_starts.py) || echo "STARTS_FAILED"
n=$(powershell -NoProfile -Command "@(Get-CimInstance Win32_Process | Where-Object { \$_.Name -like 'python*' -and \$_.CommandLine -match 'certificate_upgrades|continuous_baselines|patch_size_timing|refit_native' }).Count")
[ "${n//[^0-9]/}" = "0" ] || { echo "another experiment is running; not timing"; exit 1; }
SUMMED_REFINE_DEPTH=0 python prototype_3d/patch_size_timing.py --variant unit --fields ${mm}mm || { echo "TIMING_FAILED"; exit 1; }
python prototype_3d/continuous_baselines_report.py > /dev/null && echo RUN_VERTEX_${mm}MM_DONE
