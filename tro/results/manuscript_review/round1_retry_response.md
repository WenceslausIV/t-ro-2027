# Advisory review (manuscript e20f1c7d…f326)

This is based on reading the supplied text only; nothing was run. Arithmetic cross-checks were done by hand.

## 1. Definite errors

**D1. Only the boundary is certified, but solid enclosure is claimed (Sec. IV-B, Thm. 1).**
Prop. 1 proves $\partial\mathcal G_i\subset\{\phi_i<l_i\}$. The text then asserts $\mathcal B_i\supseteq\mathcal G_i$ and "the fit affects only the size of $l_i$, never its validity." A fitted spline has no maximum principle, so this does not follow.

- *Counterexample:* take $\mathcal G$ the closed unit disc and $\phi(\mathbf x)=\norm{\mathbf x}-1+2e^{-\norm{\mathbf x}^2/0.01}$. Then $\max_{\partial\mathcal G}\phi\approx0$, so $l\approx\epsilon$, yet $\phi(\mathbf 0)\approx1>l$ and $\mathbf 0\in\mathcal G\setminus\mathcal B$. Interior sign errors of this kind are routine for "approximate signed distances" of meshes.
- *Consequence:* Thm. 1's initial hypothesis $\{\phi_A\le l_A\}\cap\{\phi_B<l_B\}=\emptyset$ does not imply $\mathcal G_A\cap\mathcal G_B=\emptyset$ at $t=0$. A small $\mathcal G_A$ placed in the spurious interior hole of $\mathcal B_B$ satisfies every hypothesis and $h\ge0$ on $\mathcal S_A$, while intersecting $\mathcal G_B$ for all time. The proof's last step applies Lemma 1 to $\mathcal G_A,\mathcal G_B$ without the disjoint start that lemma needs.
- *Minimal repair:* add "$\mathcal G_A\cap\mathcal G_B=\emptyset$ at $t=0$" to Thm. 1 and Cor. 1. The existing proof is then correct with boundary-only certification: a first ground-truth contact lies on both boundaries, hence in $\{\phi_A\le l_A\}\cap\{\phi_B<l_B\}$.
- *Wording:* replace "$\mathcal B_i\supseteq\mathcal G_i$" and "certifies enclosure of the true geometry" by boundary enclosure in Sec. IV-B, Contribution 1, and the Conclusion.
- *Alternative:* if full-solid enclosure is wanted, run the branch-and-bound over a cover of the solid (octree of $\{d\le0\}$), certifying $\max_{\mathcal G}\phi<l$.

**D2. Prop. 1 termination is not proved for octree pieces.**
The proof shows $\mathrm{UB}\le\bar l_i+2G\rho+\tfrac12\kappa\rho^2$ for boxes, but discarding needs $\mathrm{UB}<\ell+\epsilon$, and nothing forces $\ell\to\bar l_i$.

- *Counterexample:* if $\ell$ comes from a fixed sample set with $\bar l_i-\ell\ge\epsilon$, every box containing the maximizer $\mathbf p^\star$ has $\mathrm{UB}\ge\phi(\mathbf p^\star)\ge\ell+\epsilon$ at all depths. The procedure never terminates.
- *Repair:* for each retained box, update $\ell$ with $\phi$ at an exact boundary point within $c\rho$ of its center. Then $\mathrm{UB}\le\ell+(1+c)G\rho+\tfrac12\kappa\rho^2$ and termination follows.
- *Same gap, milder, for triangles and edges:* the text says "vertices, centers, *or* projected points." The proof needs the center value always included in $\ell$.
- *Validity versus tightness:* enclosure only needs $\mathrm{UB}<l$ on every retained piece. The tightness claim $l_i\le\bar l_i+\epsilon$ additionally needs the points defining $\ell$ to lie exactly on $\partial\mathcal G_i$. An approximate projection breaks tightness, not enclosure; say so.
- *Box exclusion:* the "may contain boundary points" test must be stated and conservative (for an exact SDF, keep a box unless $|d(\mathbf c)|>\rho$).

**D3. Abstract sentence 1 and L1: convexity and closest-point derivatives.**
"Need convex level sets" and "holds only for convex level sets" are wrong in both directions.

- *Not sufficient:* a disc robot (convex level sets) beside an L-shaped obstacle has a closest point that jumps when its center crosses the obstacle's medial axis. Even with both bodies convex, two squares with parallel faces have a segment of minimizers, and the selected point jumps between corners as the relative angle passes through zero. The relevant global condition is *strict* convexity of the pair.
- *Not necessary:* a locally unique, nondegenerate minimizer (LICQ plus second-order sufficiency) gives a smooth closest point by the implicit function theorem on the KKT system. Non-convex pairs satisfy this at generic poses.
- *The derivative does not require differentiating the closest point at all.* For $H(\text{pose})=\min_{\mathbf x\in\mathcal S_A}\phi_B(\mathbf y(\mathbf x))-l_B$ with a unique minimizer, the envelope theorem gives $\dot H=\nabla\phi_B(\mathbf y^\star)^\top\dot{\mathbf y}(\mathbf x^\star)$ with $\mathbf x^\star$ held fixed. With several minimizers, $H$ is locally Lipschitz and its directional derivative is the minimum over minimizers (Danskin). A valid nonsmooth CBF constrains all (ε-)active minimizers.
- *What actually fails* is the single-closest-pair implementation, exactly as Sec. VII-D already concedes. The abstract and L1 should say the same.
- *Novelty framing:* $\mathcal H_{AB}=\{H\ge0\}$ is this same min-function's safe set; the Franka KKT baseline uses exactly it. The contribution is semi-infinite enforcement over all surface points, not a new safe set.
- *Overgeneralization:* bounding circles and the heuristic margin are properties of [jo2026geometry], not of the class of geometry-aware CBFs.
- *Suggested wording:* "The closest-point CBF of [prior] differentiates the barrier at a single minimizer, valid only where that minimizer is unique and nondegenerate (globally guaranteed for strictly convex pairs); it handles non-convex shapes through bounding circles and relies on a heuristic margin."
- *"We show that collision avoidance reduces to…"* Lemma 1 is elementary and gives a *sufficient* condition, conditional on a collision-free start. Tone it down.

**D4. The locally Lipschitz hypothesis is incompatible with the implemented filter (Thm. 2, Cor. 1).**
Hard activation at $\eta$ inserts a row that the current input generally violates. The row requires velocity terms $\ge-\gamma\eta=-0.25$ m/s (planar), while approach speeds reach 1 m/s or more, so $\mathbf u^\ast$ jumps. Online refinement at $\tfrac12M_Br^2>\vartheta$ and the slack fallback also switch the constraint set. The care taken to make $M_B$ Lipschitz is undone by these switches, so Cor. 1's hypothesis is false in general, not merely unverified.

- *Repair A (keeps the results):* the proof of Thm. 2 is pointwise in $\mathbf x$ and needs only an absolutely continuous closed-loop solution with the constraints holding for almost every $t$. Then $h_{\mathbf x}(t)$ is absolutely continuous and $\dot h_{\mathbf x}\ge-\gamma h_{\mathbf x}$ a.e. whenever $h_{\mathbf x}<\eta$, which suffices. State the theorem for Carathéodory/Filippov solutions. Filippov convexification is harmless here: each row is affine in $\mathbf u$ with continuous data, activation is an open condition, and both parent and child rows imply $\Psi\ge0$ on the surface.
- *Repair B (regularity for free):* choose $\eta\ge\gamma^{-1}\sup_{\mathbf u\in\mathcal U}|\text{velocity terms}|$. Pruned rows are then redundant and pruning does not change the minimizer. This needs $\eta$ of order decimetres with current settings.

## 2. Missing hypotheses

- **Radius restriction in Prop. 1.** $\kappa_i(\mathbf c)$ bounds the Hessian only for $\rho\le h$. The UB is therefore invalid for coarse initial triangles, long edges, or root octree boxes unless pieces are pre-split to $\rho\le h$. Also require the $\rho$-neighbourhood of $\partial\mathcal G_i$ to lie in $\Omega_i$. The $C^2$ property of $\kappa_i$ is never used there.
- **Domains.** Thm. 1 needs $\{\phi\le l\}$ compact in the *interior* of $\Omega$; otherwise its boundary includes pieces of $\partial\Omega$ where $\phi<l$. It is assumed, not certified, though a Bernstein sign check on the domain's boundary cells would do it. $\phi_B(\mathbf y(\mathbf x))$ is undefined when $\mathcal S_A$ maps outside $\Omega_B$; state the extension and show that it preserves the lower bounds used in pruning.
- **Thm. 2 and Cor. 1 are conditional on QP feasibility for all time.** Prop. 5 gives feasibility only where constant terms are nonnegative, which can fail inside $\mathcal H_{AB}$. Cor. 1's phrase "boxes enter the QP only with positive lower bounds" is ambiguous: if it means $\eta>0$, say that; if it means nonnegative constant terms, it is a strong state restriction.
- **"Falls on the inputs only."** The formal safe set is box-independent, but the set of states from which the filter is feasible depends on the boxes. "The boxes only restrict the inputs" understates this; add the qualifier "where the QP is feasible."
- **Prop. 6.** The argument is sound by compactness (retained boxes converge to $\mathcal S_A$ without needing $\nabla\phi_A\ne0$). However, it is per pose and per input, with no uniform depth and no statement about QP feasibility.
- **Frames.** State once that $(\mathbf V,\boldsymbol\Omega)$ and $\mathbf d$ are expressed in $B$'s frame, and that the Frobenius bound is used as an operator-norm bound on the third-derivative tensor.

I found no error in Lemma 1, the window argument for $\kappa_i$ and $M_B$, eq. (psi-lb), Prop. 3, the joint certificate, or the Bernstein ranges in Prop. 4.

## 3. Native-patch interpretation

- **The Taylor model of $\phi_A$ is unnecessary on a native patch.** $\phi_A-l_A$ is exactly one tensor-cubic there, with known Bernstein coefficients. Using them gives $M_A=0$ at the price of degree-elevating $q_B\circ\mathbf y$ to $4^n$ coefficients, which makes "native" mean something mathematically. At minimum, $M_A$ should be taken over the single cell: the present window of half-width $r=\sqrt n\,s/2$ pulls in neighbouring cells the box never touches.
- **"By an amount that shrinks with the box size" is vacuous for the declared default.** There $s=h_A$ is fixed. Shrinkage is proved only for the constant term at rest (Prop. 4) and for a fixed strictly feasible input (Prop. 6); there is no bound for the input-dependent terms.
- **Prop. 4's first-order loss explains the native results.** The term $\tfrac s2\norm{\mathbf b}_1$ is 5 to 20 mm for 1 to 4 cm cells, matching the 9.1 mm (dock) and 21.4 mm (crossing) gaps against 1.7 and 4.5 mm.
- **The "default" is not the evaluated method.** Natively, the Franka arm reaches the goal in 6 of 30 trials, and 2 trials use slack over 1497 steps, so the hard constraints and hence the guarantee are lost there. All headline numbers come from a non-default variant with offline halving (6 mm from 12 mm cells) plus online refinement. That variant also uses heterogeneous rules: a 1-cm restriction in the "original" runs and a distance-independent rule for the tube and Sec. VII-H. Either make subdivision the primary method or lead with the native numbers.
- **Retention test.** "All one strict sign" in floating point can drop a cell whose coefficient is within rounding of zero. Use a tolerance: retain if $\min\le\tau$ and $\max\ge-\tau$.
- **Counts to reconcile (unverified).** Fig. 1 reports 86 patches at $l=2.0$ mm; Sec. VII-B reports 87 native patches for the Fig. 2 fields, whose levels are 2.8 and 1.9 mm. Confirm these are different fields, or reconcile.

## 4. Experiment–claim consistency

- **The abstract's baseline comparison is selective.** With the 1-cm margin used in the cited work, spheres and samples have 0 collisions and minimum distances of 0.9 and 1.4 mm against 1.1 mm for the proposed method. They also reach more goals (20 and 23 against 17) and are 4 to 7 times faster. The abstract reports only the margin-free variants (17 and 19 collisions). Report both.
- **Success rates are omitted.** The subdivided method completes 12 of 20 planar instances (circles: 20 of 20) and 17 of 30 Franka trials. The abstract attributes stalling only to native patches.
- **"Certified" is overloaded.** Table III's "certified guarantee: yes" and "certified min. distance" describe a discrete-time, double-precision run whose guarantee Sec. VIII itself disclaims. The audits are spatially exhaustive but temporally sampled. With 1 rad/s joints, 10-ms steps and 1.1-mm clearance, between-step contact is not excluded by the audit. Rename the row to something like "continuous-time certificate (exact arithmetic)" and the distances to "audited lower bound at saved states."
- **The dual-arm evidence is weak.** The lower bound is evaluated at only the 20 steps per trial with the smallest *sampled* distance. "Two Pandas kept their meshes apart" needs either an all-step audit or this qualifier.
- **Scope of "with optional subdivision" in the abstract.** The Franka sentence is grammatically detached from it; attach the qualifier explicitly.
- **"Min. $h$" in Tables I and II compares different quantities.** Ours is a lower bound of $\tilde h$ on active boxes; the baseline's is its own barrier value. The boldface implies a like-for-like comparison.
- **Timing.** 13.3/46.5 ms (Franka) and 77/103 ms (tube) exceed the 10-ms step. The tube paragraph says so; the Franka paragraph and abstract should too.
- **Sec. VII-H.** One scene plus six hand-selected trials, acknowledged as exploratory; keep the conclusions at that strength.

The row and box counts I checked are internally consistent: 401×4, 255×4, 126×4, 192×4, 4374×8, 6730×8, 3961×8, 258×8, and the 8-versus-27 rows per box in VII-H. The docking depths match the gaps, and 1.7 mm equals 1.2 + 0.5 mm.

## 5. Numerical versus mathematical guarantees

- **Floating point.** The limitation is stated, but it touches three certificate steps differently.
  - Cover retention: sign tests near zero (see Sec. 3).
  - Level certificate: $\epsilon=0.1$ mm dwarfs rounding error, so a one-line error bound would make it rigorous.
  - Subdivision: $S=M_B^{-1}TM_B$ goes through the power basis. De Casteljau uses only convex combinations, is stable, and is trivially outward-roundable.
- **QP.** "Solved exactly" by NNLS with constraint generation still means rows hold to solver tolerance. Report the maximum row violation, and confirm that every active row is checked at termination. "Min. $h$ = 0.0" needs a stated tolerance.
- **Sampled data.** Between steps $h_{\mathbf x}$ may become negative. The constant terms then go negative, and zero input is no longer certified. A tightening such as $\Psi\ge L\,\Delta t\,\norm{\mathbf u}$, with $L$ built from the existing $M_B$ and $\norm{\mathbf H}$ bounds, would close the gap using quantities already computed.

## 6. Priority edits for Codex

1. Thm. 1 and Cor. 1: add initial ground-truth disjointness and replace solid-enclosure language (D1).
2. Prop. 1: require $\rho\le h$ and fix the octree termination argument and the $\ell$ update rule (D2).
3. Rewrite abstract sentence 1 and L1 (D3). Add margin-inclusive baseline results and goal-reaching rates to the abstract.
4. Replace the Lipschitz hypothesis with a solution-concept statement (D4).
5. Relabel "certified" in Table III and the dual-arm text. Resolve the default-versus-evaluated mismatch.
6. Optional strengthening: exact Bernstein form of $\phi_A$ on native patches.