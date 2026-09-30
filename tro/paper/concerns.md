# Open concerns

## 1. Non-smooth control input from box activation (2026-09-30)

**Concern.** A cover box enters the QP only when its certified lower bound of the lifted barrier drops below the
activation threshold η (5 cm planar, 3 cm spatial). At entry, its constraint roughly requires an approach speed of at
most γη (≈ 0.25 m/s planar). If the body is approaching faster, the new row is violated by the current input, and the
QP solution u* jumps. Online refinement, where a box is replaced by its children, switches rows in the same way. The
barrier h and the safe set stay smooth and continuous; the non-smoothness is in the feedback u*(x).

**Why it matters.** The standard CBF theorem (Ames et al.) assumes a locally Lipschitz feedback. That assumption gives
existence and uniqueness of closed-loop solutions, and then invariance follows. A discontinuous feedback is not
forbidden, but it falls outside that theorem, and the existence of solutions is not guaranteed. Nonsmooth-barrier work
(Glotfelter, Cortés, Egerstedt 2017, already cited) treats such feedbacks through Filippov solutions.

**Current manuscript.** Theorem `thm:invariance` is stated conditionally. It says that *if* an absolutely continuous
trajectory exists that satisfies the certificates almost everywhere, it stays safe. Limitations says that existence and
uniqueness are not established for the hard-activation QP. This is correct but weaker than the standard CBF guarantee,
and a reviewer may ask whether such a solution exists.

**Options.**
1. **Filippov strengthening (recommended; theory only, no code or experiment change).**
   - Existence: u* takes values in the compact set U, so the closed-loop right-hand side is locally bounded, and
     Filippov solutions exist (Filippov 1988), provided u* is measurable in the state.
   - Safety of every Filippov solution: for a fixed surface point x with h(x) < η₀, the true condition
     F(x, pose, u) = ḣ + γh ≥ 0 is continuous in the pose and affine in u. Every admissible input near the pose
     satisfies it through the certificates (Prop. `prop:coef` or the joint certificate). This holds for parent rows and
     child rows alike. The Filippov set is built from limits and convex combinations of such inputs, so every element
     satisfies F ≥ 0. The per-point comparison argument of `thm:invariance` then applies.
   - Remaining checks: measurability of u*(pose), since the QP data are piecewise continuous on finitely many regions;
     QP feasibility at all states, because slack voids the guarantee; the epigraph form A u + T|E u| + C ≥ 0 with
     T ≤ 0 has a convex feasible set in u.
   - Result: existence plus safety of all solutions without a Lipschitz assumption, the same strength as the standard
     theorem.
2. **Lipschitz feedback by redundancy-based activation (implementation change).** Activate a box only once its row
   holds for every u ∈ U, i.e. γ·(lower bound) ≥ max over U of the velocity terms. Pruning then never changes the QP
   solution, and with a fixed cover (native default, no online refinement) the CBF-QP is locally Lipschitz under the
   usual strict-feasibility and constraint-qualification conditions. Cost: η grows to decimeters, and the number of
   active rows and the filter time grow accordingly.

**Status.** Not yet addressed in the manuscript. The review record is `review_claude_codex.md` (C-06 and C-14), where
the a.e. / absolutely continuous reformulation was adopted but the Filippov existence argument was not.
