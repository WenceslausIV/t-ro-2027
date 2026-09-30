# Measured results (2026-09-26)

All five configurations (three method families) completed six fixed-horizon simulations. No collision or unresolved clearance interval was found. No deadlock or goal-success ranking was computed.

| Method | Matched-pose median [ms] | Matched-pose p95 [ms] | Minimum sampled physical gap [mm] | Mean input TV [m/s] |
|---|---:|---:|---:|---:|
| Ours | 0.106 | 0.161 | 29.42 | 0.89 |
| Sampling-64 | 0.159 | 0.247 | 118.75 | 1.07 |
| Sampling-256 | 0.164 | 0.259 | 39.98 | 0.85 |
| Sampling-1024 | 0.207 | 0.308 | 29.34 | 0.78 |
| Ball-world | 0.217 | 0.336 | 22.06 | 8.57 |

The field and Sampling-1024 have similar minimum observed physical gaps. The field has lower measured query cost in this implementation; Sampling-1024 has lower input TV. Ball-world runs attain a smaller gap while retaining an interval-wide lower bound above 20 mm. These observations are not a general performance ordering.

Preparation: common geometry 0.0174 s; field fitting plus contact/band/regularity checks 4.3859 s. Field runtime polynomial arrays: 460800 bytes. Sampling coordinate storage and setup timings are in comparison.json; tree allocation is excluded.

For baseline adaptations, limitations, verification details, and reproduction commands, see [README.md](README.md). Closed-loop timings are separately available in [table.md](table.md). The paper table uses matched-pose timings.
