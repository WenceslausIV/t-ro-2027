# Fig. 1 body with 100 SDF patches

Reproduce with `python prototype_3d/patch_resolution_check.py` from the project root.

The shape is GT[3] from five_robot_setup(seed=3, scale=1.0). Both fits use the same physical domain, 5-mm training grid, approximate polygon-distance targets, cubic degree, and ridge parameter 1e-4. The existing 119x90 grid is compared with exactly 10x10 rectangular patches. The thin z representation is retained. Normalized-coordinate evaluation of the fine fit matched the production body_field implementation exactly on 1,000 independent probes.

| Quantity | 10,710 patches | 100 patches |
| --- | ---: | ---: |
| Patch size [mm] | 10 x 10 | 119 x 90 |
| Effective planar controls | 11,346 | 169 |
| Coefficient solve median [ms] | 6.53 | 1.50 |
| Data generation + fitting median [s] | 1.05 | 1.16 |
| Certified boundary level [mm] | 0.271 | 16.828 |
| SDF RMSE within 30 mm of boundary [mm] | 0.0467 | 4.801 |
| Sampled certified-contour gap median / max [mm] | 0.270 / 0.429 | 17.380 / 36.865 |

The 100-patch fit retains the overall shape, but its errors around concave features require a larger enclosing level. The level is certified on continuous polygon edges by branch-and-bound, using a global physical planar Hessian bound from Bernstein coefficients and a 0.1-mm termination tolerance. As in the production implementation, this uses double precision without outward rounding.

Error metrics use exact signed point-to-segment distances on a separate 481x361 diagnostic grid. Contour gaps use vertices extracted by contour interpolation. These two diagnostics are sampled measurements, not continuous bounds. No control simulation was performed and no claim about goal success follows from this experiment.

Timings exclude imports, use one BLAS thread and three repetitions; background activity was not excluded. The fixed distance-training data dominate total time, so the total-time difference does not establish a slowdown caused by fewer patches. Reducing the training-data count is a separate experiment.

See comparison.png (or comparison.pdf) for the true polygon, zero contour, certified level contour and knot grid. Full measurements are in results.json. Production paper figures and simulation results are unchanged.
