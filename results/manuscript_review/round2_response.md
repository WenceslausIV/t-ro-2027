Nothing below was run; these are proofs and counterexamples from the supplied text only.

## 1. Revised first-contact theorem

**Verdict: valid after one added hypothesis and one convention; as written there is a counterexample.**

**Counterexample (domain face).** "Compact subset of its domain" lets $\{\phi_A\le l_A\}$ touch $\partial\Omega_A$, and then $\partial\{\phi_A\le l_A\}\not\subset\mathcal S_A$.
- Take $\Omega_A=[0,1]^2$, $\phi_A(\mathbf x)=x_1-\tfrac12$ (an exact tensor cubic spline, by linear precision), $l_A=0$, and $\mathcal G_A=[0.1,0.4]\times[0.1,0.9]$. Then $\partial\mathcal G_A\subset\{\phi_A<0\}$ and $\mathcal S_A=\{x_1=\tfrac12\}$.
- Let $\mathcal G_B$ be a disc of radius 0.05 with $\phi_B=\norm{\mathbf y}-0.06$, $l_B=0$, starting at $(-1,0.5)$ and translating in $+x_1$ to $(0.25,0.5)$.
- All hypotheses hold and $h\ge0$ on $\mathcal S_A$ throughout, since the disc never reaches $x_1=\tfrac12$. Yet $\mathcal G_B$ enters $\mathcal G_A$ through the face $x_1=0$.

**Repair.**
- Require $\{\phi_A\le l_A\}\subset\operatorname{int}\Omega_A$, or equivalently $\phi_A>l_A$ on $\partial\Omega_A$. This is certifiable from the Bernstein coefficients of the boundary cells.
- State that $h(\mathbf x)\ge0$ is imposed only where $\mathbf y(\mathbf x)\in\Omega_B$. Currently $h$ is undefined elsewhere.

**Proof check with these additions.**
1. If $\phi_A(\mathbf z)<l_A$ at $\mathbf z\in\operatorname{int}\Omega_A$, continuity makes $\mathbf z$ interior to $\{\phi_A\le l_A\}$. Hence $\partial\{\phi_A\le l_A\}\subset\mathcal S_A$.
2. Lemma 1 applied to $\{\phi_A\le l_A\}$ and $\{\phi_B\le l_B-\epsilon\}$ gives a first-contact point $\mathbf z\in\mathcal S_A$ with $\phi_B(\mathbf z)\le l_B-\epsilon$. It lies in $\Omega_B$, so $h(\mathbf z)<0$ is a genuine violation.
3. Only $\mathbf z\in\partial K_A$ and $\mathbf z\in K_B^\epsilon$ are used, so no interior condition on $B$ is needed for this theorem. The proof's "$\phi_B=l_B-\epsilon$" should read "$\le$".
4. Taking the union over $\epsilon$ gives $\{\phi_A\le l_A\}\cap\{\phi_B<l_B\}=\emptyset$ for all time.
5. The ground-truth step is correct with the new initial disjointness: a first contact lies on $\partial\mathcal G_A\cap\partial\mathcal G_B\subset\{\phi_A<l_A\}\cap\{\phi_B<l_B\}$.

"Certified regions may touch" in the corollary is consistent with this.

## 2. Measurable-input replacement for the Lipschitz theorem

**Verdict: true, with the proof below. Null sets depending on $\mathbf x$ are harmless.**

**Standing assumptions.**
- $\gamma>0$; $\phi_B\in C^1$ (the spline is $C^2$).
- $\mathbf y(\mathbf x,\cdot)\in\Omega_B$ on the interval considered (item 3 removes this).
- The certificate at time $t$ uses the same $\mathbf u(t)$ that appears in the dynamics at $t$, with $\mathbf a(t)\ge|(\mathbf V,\boldsymbol\Omega)(t)|$.

**Proof.** Fix $\mathbf x\in\mathcal S_A$.
1. $h_{\mathbf x}(t)=\phi_B(\mathbf y(\mathbf x,\text{pose}(t)))-l_B$ is a $C^1$ function composed with an absolutely continuous curve, hence absolutely continuous.
2. Wherever the pose is differentiable and the dynamics hold, $\dot h_{\mathbf x}=\nabla\phi_B(\mathbf y)^\top(\mathbf J\mathbf u+\mathbf W(\mathbf u)\mathbf d)$.
3. By hypothesis and Prop. 3 (or the joint certificate with $g_A(\mathbf x)=0$), $\dot h_{\mathbf x}\ge-\gamma h_{\mathbf x}$ for $t\in\{h_{\mathbf x}<\eta\}\setminus N_{\mathbf x}$, with $N_{\mathbf x}$ null.
4. Suppose $h_{\mathbf x}(t_1)<0$. Let $t_0=\sup\{t\le t_1: h_{\mathbf x}(t)\ge0\}$, which exists because $h_{\mathbf x}(0)\ge0$. Continuity gives $h_{\mathbf x}(t_0)=0$ and $h_{\mathbf x}<0<\eta$ on $(t_0,t_1]$.
5. $g=e^{\gamma t}h_{\mathbf x}$ is absolutely continuous with $g'\ge0$ a.e. on $(t_0,t_1)$. So $g(t_1)-g(t_0)=\int g'\ge0$, giving $h_{\mathbf x}(t_1)\ge0$, a contradiction. ∎

The argument is run separately for each $\mathbf x$, so the uncountable union $\bigcup_{\mathbf x}N_{\mathbf x}$ is never formed.

**Remarks for the text.**
- **Absolute continuity cannot be weakened.** "Continuous and a.e. differentiable" fails: a negated Cantor function has $\dot h=0\ge-\gamma h$ a.e. wherever $h<0$, yet decreases from 0.
- **Only $\{h_{\mathbf x}<0\}$ is used.** So $\eta>0$ is a robustness margin, not a logical requirement.
- **The theorem is conditional.** It asserts nothing about existence and is vacuous if no such trajectory exists.
- **It does not cover the implementation.** A zero-order hold evaluates the certificate at the sampled pose and applies the input over the whole step, so the hypothesis holds only at the sample instants, a null set. Steps that use slack violate the hypothesis outright.
- **A simpler sufficient hypothesis.** "For a.e. $t$, every active box satisfies its rows" uses one null set and implies the per-point version.

## 3. Boxes dropped near $\partial\Omega_B$

**Verdict: you are right that compactness is insufficient.** It gives a collar of unknown positive width, while a dropped box reaches a width set by $r$.

**Geometry of the drop.** Let $\rho_{\rm drop}\ge r$ be the margin the code actually tests. Audit this: it may be $r$, or $r$ plus the majorant window, with a Euclidean or cube test. A box is dropped when $\operatorname{dist}(\mathbf y_{\mathbf c},\R^n\setminus\Omega_B)<\rho_{\rm drop}$. Every point of a dropped box then lies outside $\Omega_B$ or within $w=\rho_{\rm drop}+r_{\max}$ of $\partial\Omega_B$. Here $r_{\max}$ is the largest half-diagonal among all covers paired with $B$. If clusters are dropped by the same test, use their radius $R$.

**Minimal collar condition (C).**
$$\phi_B(\mathbf y)-l_B\ge0\quad\text{for all }\mathbf y\in\Omega_B\text{ with }\operatorname{dist}(\mathbf y,\partial\Omega_B)\le w .$$
Use $\ge\eta$ if the text keeps "pruned boxes have $h\ge\eta$."

**Why (C) suffices.** Extend item 2 by treating $h_{\mathbf x}$ as unconstrained while $\mathbf y\notin\Omega_B$.
1. If $h_{\mathbf x}(t_1)<0$, then $\mathbf y(t_1)$ is in $\Omega_B$ and outside the collar.
2. Let $t_0=\sup\{t\le t_1:\mathbf y\notin\Omega_B\text{ or }h_{\mathbf x}\ge0\}$. On $(t_0,t_1]$ the point stays outside the collar, so its box is not dropped and the certificate applies.
3. At $t_0$, either $h_{\mathbf x}(t_0)\ge0$ directly, or $\mathbf y(t_0)\in\partial\Omega_B$, which lies in the collar. Either way $h_{\mathbf x}(t_0)=0$, and the comparison argument concludes.

**Why (C) is minimal.** If $\phi_B<l_B$ somewhere within $w$ of $\partial\Omega_B$, a pose exists in which a dropped box carries a surface point there.

**Offline certificate.** Every Bernstein coefficient of $\phi_B-l_B$ (or $-\eta$) is nonnegative on every cell meeting the collar, i.e. the outer $\lceil w/h_B\rceil$ layers. This is a property of the fitted coefficients, not of the fit target. The 15-cm padding makes it plausible, but edge control points are weakly constrained by the fit.

**Explicit skip certificate (alternative).**
- Precompute per cell $m_{\rm cell}=\min$ Bernstein coefficient of $\phi_B-l_B$.
- Skip a domain-violating box only if $m_{\rm cell}\ge\eta$ on every cell of $\Omega_B$ meeting the cube of half-width $r$ about $\mathbf y_{\mathbf c}$.
- If the test fails, the step is uncertified and must be reported as such, not silently dropped.
- Under (C) this test never fails; without (C) it can.

**What must be audited before claiming implementation safety.**
1. The actual drop test and its margin, including the norm used.
2. $r_{\max}$ over all partners. Native patches are the worst case: $r=\sqrt n\,h_A/2$, about 2.8 cm for 4-cm planar cells.
3. The coefficient certificate for each $B$ at that $w$, with a floating-point margin.
4. That out-of-domain evaluation cannot return clamped or extrapolated values that pass screening.
5. A runtime count of dropped boxes and of failed skip tests.

Item 3 is the theorem's hypothesis. Items 1, 2 and 4 establish that the code matches it.

## 4. Min-of-rounded-box SDF in the level certificate

**Properties of $d=\min_i d_i$.**
- $d$ is 1-Lipschitz and vanishes on $\partial\mathcal G$.
- Outside the union, $d$ is the exact distance.
- Inside, $|d|$ is only a lower bound on the distance to $\partial\mathcal G$. With crossing slabs or tangential overlaps it can be far smaller than the true distance, with no uniform constant.
- $\{d=0\}\supseteq\partial\mathcal G$. The inclusion is strict if two primitives abut exactly, because an internal seam has $d=0$.

**Validity of a completed enclosure: holds, independently of $\ell$.**
- If a box contains $\mathbf p\in\partial\mathcal G$, then $|d(\mathbf c)|\le\norm{\mathbf c-\mathbf p}\le r$, so it is retained. Interior inexactness is irrelevant.
- Every retained box is discarded only when $\mathrm{UB}<\ell+\epsilon\le l$. So $\partial\mathcal G\subset\{\phi<l\}$ for any value of $\ell$, attained or not.
- Conditions to audit:
  - $r$ is the half-diagonal, not the half-side.
  - The same $r$ enters the UB.
  - $\rho\le h$, or a radius-dependent Hessian window is used.
  - The root boxes cover $\partial\mathcal G$.
  - The queue is truly empty; no piece is dropped at a depth cap.
  - Each primitive SDF is exactly 1-Lipschitz.

**Tightness ($l\le\bar l+\epsilon$): not established by the current code.**
- It needs $\ell\le\bar l$, so every value entering $\ell$ must be a rigorous lower bound of $\phi$ at a true boundary point.
- A projected point with residual below $10^{-9}$ is not on $\partial\mathcal G$.
- If it landed slightly inside, the residual does not even bound its distance to $\partial\mathcal G$, by the interior property above.

**Termination: not guaranteed.**
- The revised text's "$C=1$ when $|\mathrm{SDF}(\mathbf c)|\le\rho$" is true for outside centers only.
- For inside centers the oracle gives no boundary point within $C\rho$. Near an internal seam or a thin exterior cusp no outside center may exist at a given level, so $\ell$ is never raised there.
- If $\phi$ on an internal seam exceeds $\ell+\epsilon$, refinement never stops.

**Minimal sound repair.**
1. **Deflated lower bound.** Accept a projected point $\mathbf p$ only if $0\le d(\mathbf p)\le\tau$, on the exact side. A true boundary point $\mathbf b$ then exists with $\norm{\mathbf b-\mathbf p}=d(\mathbf p)$. Update $\ell\leftarrow\max(\ell,\ \phi(\mathbf p)-G_{\rm loc}\,d(\mathbf p)-\delta_{\rm fp})$, with $G_{\rm loc}$ a certified gradient bound from the first-derivative Bernstein coefficients of the cells meeting the ball. This gives $\ell\le\bar l$ rigorously, without treating the residual as zero.
2. **Depth cap with absorption.** Stop at a fixed depth and return $l=\max(\ell+\epsilon,\ \max_{\text{leftover}}\mathrm{UB})$ plus a rounding margin. The procedure terminates by construction and validity is unconditional.
3. **Restate Prop. 1 for implicit bodies.** Claim "returns a valid $l$ with the a posteriori gap $l-\ell\ge l-\bar l$ reported" instead of "$l\le\bar l+\epsilon$ and terminates for every $\epsilon$." If the theoretical termination statement is kept, its hypothesis should be the following:
   - every retained box supplies a point of $\{d=0\}$ within $\omega(\rho)$ of its center, with $\omega\to0$ (a modulus suffices; linear $C\rho$ is not needed);
   - the target is $\max_{\{d=0\}}\phi$, which equals $\bar l$ only when no primitives abut exactly.