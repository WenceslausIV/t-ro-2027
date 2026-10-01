#!/bin/bash
# Fair comparison (PLAN.md, "LOCAL CONTINUATION"), on the user's machine (AMD Ryzen 9 3900X, Windows).
# Restartable: every runner skips trials whose json already exists.
#   bash tro/results/continuous_baselines/run_fair.sh outcomes   # phase A, parallel queue (max 10 jobs)
#   bash tro/results/continuous_baselines/run_fair.sh timing     # phase B, one process at a time, nothing else running
cd "$(dirname "$0")/../.."            # -> tro/
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
CU=prototype_3d/certificate_upgrades.py
CB=prototype_3d/continuous_baselines.py
LOG=results/continuous_baselines/logs
mkdir -p $LOG

outcomes() {
  # Every method: our fitted obstacle fields, rows at the 10-ms samples, DAQP with exact check, slack fallback,
  # mesh audit at EVERY state (+ final).
  {
    for r in "0 15" "15 30"; do
      echo "SUMMED_REFINE_DEPTH=0 python $CU --variant free --qp daqp --tag fixed --trials $r"
      echo "SUMMED_REFINE_DEPTH=0 python $CU --variant free --field 6mm --qp daqp --tag fixed --trials $r"
      echo "SUMMED_REFINE_DEPTH=1 python $CU --variant free --qp daqp --tag fixed_refine1 --trials $r"
      echo "SUMMED_REFINE_DEPTH=1 python $CU --variant free --field 6mm --qp daqp --tag fixed_refine1 --trials $r"
      echo "python $CB --methods points_delta --obstacle fitted --audit-every 1 --tag _audit1 --trials $r"
      echo "python $CB --methods points_delta --obstacle fitted --edge-mm 5 --audit-every 1 --tag _e5 --trials $r"
      echo "python $CB --methods capsule --obstacle fitted --audit-every 1 --tag _audit1 --trials $r"
      echo "python $CB --methods spheres_enclosing --obstacle fitted --audit-every 1 --tag _audit1 --trials $r"
      echo "python $CB --methods spheres_kmeans --spheres 16 --obstacle fitted --audit-every 1 --tag _audit1 --trials $r"
      echo "python $CB --methods spheres_kmeans --spheres 64 --obstacle fitted --audit-every 1 --tag _audit1 --trials $r"
      echo "python $CB --methods spheres_kmeans --spheres 128 --obstacle fitted --audit-every 1 --tag _audit1 --trials $r"
    done
  } | xargs -P 10 -I{} bash -c '{} >> '"$LOG"'/outcomes_$(echo "{}" | md5sum | cut -c1-8).log 2>&1'
  echo OUTCOMES_DONE
}

timing() {
  # Isolated: one process at a time; check that no other experiment runs before starting.
  n=$(powershell -NoProfile -Command "@(Get-CimInstance Win32_Process | Where-Object { \$_.Name -like 'python*' -and \$_.CommandLine -match 'certificate_upgrades|continuous_baselines|patch_size_timing' }).Count")
  if [ "${n//[^0-9]/}" != "0" ]; then echo "another experiment is running ($n); not timing"; exit 1; fi
  SUMMED_REFINE_DEPTH=0 python prototype_3d/patch_size_timing.py --variant free --fields 12mm 6mm
  SUMMED_REFINE_DEPTH=1 python prototype_3d/patch_size_timing.py --variant free --fields 12mm 6mm
  SUMMED_REFINE_DEPTH=0 python prototype_3d/patch_size_timing.py --variant sampled --fields 12mm 6mm
  for m in "points_delta" "points_delta --edge-mm 5" "capsule" "spheres_enclosing" \
           "spheres_kmeans --spheres 16" "spheres_kmeans --spheres 64" "spheres_kmeans --spheres 128"; do
    t=_timing_local; case "$m" in *edge-mm*) t=_e5_timing_local;; esac
    python $CB --methods $m --obstacle fitted --audit-every 0 --tag $t
  done
  # exact-SDF reference rows, re-timed on this machine
  for m in "points_delta" "capsule" "spheres_enclosing" "spheres_kmeans --spheres 16" "spheres_kmeans --spheres 64"; do
    python $CB --methods $m --obstacle exact --audit-every 0 --tag _timing_local
  done
  echo TIMING_DONE
}

"$@"
