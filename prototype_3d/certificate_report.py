"""Summarize completed pilot runs without changing their stored measurements."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    records = [(p, json.loads(p.read_text(encoding='utf-8'))) for p in sorted(args.directory.glob('*.json'))]
    records = [(p, d) for p, d in records if 'metadata' in d and 'metrics' in d]
    lines = ['# Certificate pilot results', '',
             'Local single-thread measurements; simulated time and CPU time are distinct.',
             'The original paper result files are not replaced by this comparison.', '',
             '| Run | Reached [s] | Distance bound [mm] | Slack steps | Input modification integral | Median/p95 [ms] | Max rows |',
             '|---|---:|---:|---:|---:|---:|---:|']
    hashes = {}
    for path, d in records:
        m, meta = d['metrics'], d['metadata']
        hashes.setdefault(meta['scene'], set()).add(meta['geometry_sha256'])
        reached = '--' if m['reached_s'] is None else f"{m['reached_s']:.2f}"
        lines.append(f"| {path.stem} | {reached} | {m['min_gap_bound_mm']:.4f} | {m['slack_steps']} | "
                     f"{m['control_change_integral']:.5f} | {m['t_median_ms']:.2f}/{m['t_p95_ms']:.2f} | {m['rows_max']} |")
    lines += ['', '## Shared geometry check', '']
    for scene, values in hashes.items():
        if len(values) != 1:
            raise ValueError(f'Geometry changed within scene {scene}: {values}')
        lines.append(f'- {scene}: one common SHA-256 across its completed variants (`{next(iter(values))}`).')
    lines += ['', '## Interpretation', '',
              '- Distance bounds below zero are inconclusive, not confirmed collisions.',
              '- The 3D pilot uses the existing statewise distance routines, not a continuous-time audit.',
              '- All comparisons use the same unit-lift screening at a given pose.',
              '- Fixed multipliers change the certificate only. They do not remove the boxes.',
              '- Read README.md for the meaning of the screening h bound and goal flags.', '']
    target = args.directory.parent / 'REPORT.md'
    target.write_text('\n'.join(lines), encoding='utf-8')
    print(target)


if __name__ == '__main__':
    main()
