# Certificate pilot results

Local single-thread measurements; simulated time and CPU time are distinct.
The original paper result files are not replaced by this comparison.

| Run | Reached [s] | Distance bound [mm] | Slack steps | Input modification integral | Median/p95 [ms] | Max rows |
|---|---:|---:|---:|---:|---:|---:|
| dock_bernstein-norm | -- | 1.6498 | 0 | 0.65255 | 10.97/13.44 | 3753 |
| dock_bernstein-one | -- | 1.6788 | 0 | 0.65277 | 10.62/12.80 | 3753 |
| dock_bernstein-proj | -- | 1.6497 | 0 | 0.65289 | 11.02/13.54 | 3753 |
| dock_vertex-norm | -- | 1.6569 | 0 | 0.65542 | 8.60/11.20 | 1672 |
| dock_vertex-one | -- | 1.6810 | 0 | 0.65562 | 8.33/10.16 | 1672 |
| dock_vertex-proj | -- | 1.6569 | 0 | 0.65530 | 8.59/10.89 | 1672 |
| franka_00_bernstein-one | -- | 2.2763 | 0 | 5.47085 | 82.84/161.47 | 95418 |
| franka_00_vertex-one | -- | 2.2741 | 0 | 5.53833 | 45.89/66.17 | 25096 |
| franka_01_bernstein-one | 6.63 | 4.3670 | 0 | 3.53343 | 72.62/249.63 | 139239 |
| franka_01_vertex-one | 7.38 | 4.5798 | 0 | 4.24282 | 46.00/115.53 | 42408 |
| franka_02_bernstein-one | 4.79 | 16.4684 | 0 | 0.97425 | 6.26/39.33 | 20520 |
| franka_02_vertex-one | 4.83 | 16.5548 | 0 | 1.03954 | 6.16/27.07 | 5432 |
| franka_03_bernstein-one | -- | 1.5663 | 0 | 9.23301 | 253.11/282.71 | 146691 |
| franka_03_vertex-one | -- | 2.0434 | 0 | 9.28317 | 114.57/133.42 | 41776 |
| franka_04_bernstein-one | 4.77 | 30.4680 | 0 | 0.07062 | 2.86/12.25 | 2700 |
| franka_04_vertex-one | 4.77 | 30.4680 | 0 | 0.07074 | 2.66/12.24 | 800 |
| franka_05_bernstein-one | 4.01 | 17.6877 | 0 | 0.98453 | 8.06/20.43 | 11988 |
| franka_05_vertex-one | 4.04 | 18.9586 | 0 | 1.03415 | 7.81/16.74 | 3152 |
| star_bernstein-norm | 3.92 | 3.4095 | 0 | 1.85718 | 52.19/165.61 | 89424 |
| star_bernstein-one | 3.93 | 3.4232 | 0 | 1.87049 | 52.17/165.72 | 89424 |
| star_bernstein-proj | 3.95 | 3.4477 | 0 | 1.89570 | 50.99/166.47 | 89424 |
| star_vertex-norm | 3.99 | 3.5102 | 0 | 1.95149 | 33.90/85.29 | 26432 |
| star_vertex-one | 3.99 | 3.5046 | 0 | 1.94767 | 34.65/83.73 | 26440 |
| star_vertex-proj | 4.00 | 3.5213 | 0 | 1.96451 | 34.73/83.25 | 26432 |

## Shared geometry check

- dock: one common SHA-256 across its completed variants (`1bf159f020226c65a7c0c0394a6d4bf3293709d833592de8f84df9ccc53f8c75`).
- franka: one common SHA-256 across its completed variants (`a676a509db7e6b45cbe145150d2beda6bafb055c64638093ed21ccd0333404e0`).
- rounded star tube: one common SHA-256 across its completed variants (`2c3b6b0976223209eaa5107b2483b72660c616b550fb7e6066f247044bcc252f`).

## Interpretation

- Distance bounds below zero are inconclusive, not confirmed collisions.
- The 3D pilot uses the existing statewise distance routines, not a continuous-time audit.
- All comparisons use the same unit-lift screening at a given pose.
- Fixed multipliers change the certificate only. They do not remove the boxes.
- Read README.md for the meaning of the screening h bound and goal flags.
