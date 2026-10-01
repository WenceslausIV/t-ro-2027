#!/bin/bash
# Star tube, 30 random goals: ours (vertex, 6-mm box = patch) vs full sphere decomposition (links and tube),
# three sphere settings (the best within the 20-Hz budget is reported). Restartable (finished trials are skipped).
#   bash tro/results/star_random_goals/run_star.sh outcomes   # every state audited, queue of max 10 jobs
#   bash tro/results/star_random_goals/run_star.sh timing     # isolated: one process, no audit, nothing else running
cd "$(dirname "$0")/../.."            # -> tro/
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
P=prototype_3d/star_random_goals.py
LOG=results/star_random_goals/logs
mkdir -p $LOG
SPH=("--link-spheres 64 --tube-spheres 256 --voxel-mm 10 --tag _L64_T256"
     "--link-spheres 64 --tube-spheres 1024 --voxel-mm 5 --tag _L64_T1024"
     "--link-spheres 128 --tube-spheres 1024 --voxel-mm 5 --tag _L128_T1024")

outcomes() {
  {
    for r in "0 5" "5 10" "10 15" "15 20" "20 25" "25 30"; do echo "python -u $P run ours --trials $r"; done
    for s in "${SPH[@]}"; do for r in "0 10" "10 20" "20 30"; do echo "python -u $P run spheres $s --trials $r"; done; done
  } | xargs -P 10 -I{} bash -c '{} >> '"$LOG"'/out_$(echo "{}" | md5sum | cut -c1-8).log 2>&1'
  echo STAR_OUTCOMES_DONE
}

timing() {
  n=$(powershell -NoProfile -Command "@(Get-CimInstance Win32_Process | Where-Object { \$_.Name -like 'python*' -and \$_.CommandLine -match 'star_random_goals|certificate_upgrades|continuous_baselines|patch_size_timing' }).Count")
  [ "${n//[^0-9]/}" = "0" ] || { echo "another experiment is running; not timing"; exit 1; }
  python -u $P run ours --no-audit --tag _timing
  for s in "${SPH[@]}"; do python -u $P run spheres ${s/--tag _/--tag _timing_} --no-audit; done
  echo STAR_TIMING_DONE
}

"$@"
