#!/bin/bash
# Star tube, 30 random goals: ours (vertex, 6-mm box = patch) vs full sphere decomposition (links and tube),
# three sphere settings. Restartable (finished trials and audits are skipped).
#   bash tro/results/star_random_goals/run_star.sh sims     # isolated, one process at a time, no audit -> timing
#   bash tro/results/star_random_goals/run_star.sh audits   # certified mesh-distance bounds of every saved state
#   bash tro/results/star_random_goals/run_star.sh report
cd "$(dirname "$0")/../.."            # -> tro/
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
P=prototype_3d/star_random_goals.py
SPH=("--link-spheres 64 --tube-spheres 256 --voxel-mm 10 --tag _L64_T256"
     "--link-spheres 64 --tube-spheres 1024 --voxel-mm 5 --tag _L64_T1024"
     "--link-spheres 128 --tube-spheres 1024 --voxel-mm 5 --tag _L128_T1024")

sims() {
  n=$(powershell -NoProfile -Command "@(Get-CimInstance Win32_Process | Where-Object { \$_.Name -like 'python*' -and \$_.CommandLine -match 'star_random_goals|certificate_upgrades|continuous_baselines|patch_size_timing' }).Count")
  [ "${n//[^0-9]/}" = "0" ] || { echo "another experiment is running; not timing"; exit 1; }
  python -u $P run ours --no-audit
  for s in "${SPH[@]}"; do python -u $P run spheres $s --no-audit; done
  echo STAR_SIMS_DONE
}

audits() {
  for f in ours spheres_L64_T256 spheres_L64_T1024 spheres_L128_T1024; do python -u $P audit --folder $f --workers 12; done
  echo STAR_AUDITS_DONE
}

report() { python prototype_3d/star_random_goals_report.py && echo STAR_REPORT_DONE; }

joint() {
  # joint Bernstein certificate (w = 1), box = patch at 6 and 8 mm: isolated sims, then audits and report
  n=$(powershell -NoProfile -Command "@(Get-CimInstance Win32_Process | Where-Object { \$_.Name -like 'python*' -and \$_.CommandLine -match 'star_random_goals|certificate_upgrades|continuous_baselines|patch_size_timing' }).Count")
  [ "${n//[^0-9]/}" = "0" ] || { echo "another experiment is running; not timing"; exit 1; }
  for f in 6mm 8mm; do python -u $P run ours --variant joint --field $f --no-audit --tag _joint_$f; done
  for f in 6mm 8mm; do python -u $P audit --folder ours_joint_$f --workers 12; done
  report
  echo STAR_JOINT_DONE
}

"$@"
