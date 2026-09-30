# Fig. 1 with approximately 1-cm geometric tolerance

The selection target is the maximum sampled distance from the certified level contour to the true polygon, not just SDF RMSE or the numerical value of the certified level. These distances do not constitute a continuous 1-cm gap certificate.

Tested square cell sizes, coarsest first: 8 cm, then 4 cm. Both shapes met the target at 4 cm. This is a choice among tested resolutions, not a minimum-patch solution.

| Body | Patches | Certified level [mm] | Sampled maximum contour gap [mm] |
| --- | ---: | ---: | ---: |
| B_i | 30 x 23 = 690 | 2.834 | 7.879 |
| B_j | 25 x 26 = 650 | 1.854 | 3.043 |

The same seed-3 shapes, 5-mm fitting sample spacing, cubic basis and ridge parameter are retained. The domain uses the original 10-cm padding rounded up to whole cells. Boundary levels use the production continuous edge certificate. Contours are extracted on a 1-mm grid and compared to exact polygon segments. The 5-mm Bernstein surface cover of B_i has 718 boxes.

In the updated visualization, dashed blue squares explicitly show the 4-cm parent SDF patches near the surface. Small solid/faint boxes show the separate 5-mm surface cover. A legend distinguishes these two sizes; changing the SDF resolution does not automatically change the surface-cover box size.

Reproduce selection: `python prototype_3d/overview_resolution.py`.
Reproduce the paper figure: `python planar/make_paper_figs.py overview`.

`overview_before.png` preserves the previous paper figure; `overview.png` is the replacement, also saved to `paper/figs/overview.png`. The caption in `paper/main.tex` records the illustration's new resolution. Other experiments retain their existing resolutions. The manuscript compiled successfully to 11 pages.
