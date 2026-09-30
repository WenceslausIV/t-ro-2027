Codex response to your review: I accepted the boundary/solid issue and independently proved it with an exact C2 tensor cubic spline counterexample. I added initial ground-truth disjointness, removed the full-solid claim, fixed the closed-loop corollary (certified boundaries may touch), and added an octree boundary-oracle condition. I reject strict convexity alone as a blanket sufficient regularity statement, and will not assert Filippov existence/admissibility without proof. Native patches remain the user-selected design; the true 6-mm field-refit pilot was added (one trial, no general performance claim). You are advisory only. Please answer four narrow questions with at most 1800 words:
1. Validate or find a counterexample to the revised first-contact theorem below.
2. We can replace the Lipschitz-feedback theorem by: given an absolutely continuous pose trajectory and measurable locally integrable inputs satisfying the dynamics a.e., the active-box certificate holds for each fixed surface point for a.e. time when h_x<eta; then all h_x>=0 if initially nonnegative. Null sets may depend on x. Prove or refute this without claiming existence of a closed-loop solution or arbitrary Filippov selection.
3. Finite-domain code drops boxes whose balls are not contained in B's spline domain. Compactness of {phi_B<=l_B} alone does NOT certify this pruning for large boxes. Formulate the minimal domain-collar condition or an explicit conservative skip certificate, and explain which must be numerically audited before claiming implementation safety.
4. The code uses min of rounded-box SDFs (exact distance outside the union, generally not exact signed distance inside). Existing implicit-level code retains |sdf(center)|<=r and updates LB from outside-center finite-difference projections with residual <1e-9. Distinguish validity of completed upper-bound enclosure from tightness and guaranteed termination. Suggest a minimal sound repair without trusting a nearly-zero residual as an exact attained value.
Do not suggest a generic sampled-data tightening without a proved uniform bound. Do not claim any tests were run.

CURRENT REVISED THEORY:
\section{Certified Regions and First Contact}\label{sec:regions}

\subsection{First Contact}
\begin{lemma}[First contact]\label{lem:contact}
Let $\mathcal K_A(s),\mathcal K_B(s)\subset\R^n$, $s\in[0,T]$, be images of compact sets under rigid
motions that are continuous in $s$, disjoint at $s=0$ and intersecting at $s=T$. Then at
$s^\star=\inf\{s\mid\mathcal K_A(s)\cap\mathcal K_B(s)\neq\emptyset\}>0$, some point lies on both
boundaries, $\partial\mathcal K_A(s^\star)\cap\partial\mathcal K_B(s^\star)\neq\emptyset$.
\end{lemma}
\begin{proof}
The set of intersecting parameters is closed, so the bodies intersect at $s^\star>0$; pick $\mathbf z$
in both. If $\mathbf z$ were interior to $\mathcal K_B(s^\star)$, the point of $\mathcal K_A$ that
is at $\mathbf z$ at $s^\star$ would, by continuity, lie in $\mathcal K_B(s)$ for $s<s^\star$ close
to $s^\star$, a contradiction. Symmetrically, $\mathbf z$ is not interior to $\mathcal K_A(s^\star)$,
since the point of $\mathcal K_B$ at $\mathbf z$ would lie in $\mathcal K_A$ before $s^\star$.
\end{proof}
The lemma requires neither convexity nor smoothness. It reduces collision avoidance to excluding
coincident boundary points, provided that the trajectory starts collision-free.

\subsection{Certified Regions}\label{sec:cert-body}
For each body we construct a sublevel region $\mathcal B_i$ of a fitted SDF and certify
$\partial\mathcal G_i\subset\{\phi_i<l_i\}$. This boundary certificate suffices for the
first-contact argument below. It does not by itself establish $\mathcal G_i\subseteq\mathcal B_i$:
an unconstrained fitted field can have an interior positive bump. We do not assume full-solid
enclosure without an additional interior certificate.

\emph{Level sets of SDFs.} Let $\phi_i$ be a tensor-product cubic B-spline fitted to (approximate)
signed distances of $\mathcal G_i$ on a box $\Omega_i$, and let $\kappa_i(\mathbf c)$ be a certified
bound on $\norm{\nabla^2\phi_i}$ over the ball of radius $\rho\le h$ (the cell size) around $\mathbf c$: a
cubic B-spline whose control value $m$ is the maximum of the cellwise bounds from the Bernstein
coefficients of the second derivatives over cells $m-4,\dots,m+1$ per axis, which is $C^2$ and, its
weights being a partition of unity, bounds every cell adjacent to that of $\mathbf c$. The smallest level that encloses the boundary is
$\bar l_i=\max_{\partial\mathcal G_i}\phi_i$, and we certify it by branch-and-bound over pieces of the
boundary: triangles for a mesh, edges for a polygon, and, for a body given by an exact SDF, boxes of an
octree that may contain boundary points. For a piece $\mathcal T$ with center $\mathbf c$ and radius
$\rho$ (every point within $\rho$ of $\mathbf c$), Taylor's theorem gives
\begin{equation}\label{eq:level-ub}
\max_{\mathcal T}\phi_i\le\mathrm{UB}(\mathcal T)=\phi_i(\mathbf c)+\norm{\nabla\phi_i(\mathbf c)}\rho
+\tfrac12\kappa_i(\mathbf c)\rho^2 ,
\end{equation}
and $\phi_i$ at vertices, centers, or certified projected points on $\partial\mathcal G_i$ gives
attained values whose maximum $\ell$ lower-bounds $\bar l_i$. For each retained octree box we
require a boundary point within $C\rho$ of its center, with a constant $C$ independent of
refinement; for an exact SDF, a nearest boundary point gives $C=1$ when
$|\mathrm{SDF}(\mathbf c)|\le\rho$. This is a requirement on the boundary oracle, not a
consequence of evaluating a fitted field. Pieces with $\mathrm{UB}<\ell+\epsilon$ are discarded;
the others are split (a triangle into four, an edge or a box into halves or octants) and re-bounded.
\begin{proposition}\label{prop:level}
When no piece is left, $l_i=\ell+\epsilon$ satisfies $\partial\mathcal G_i\subset\{\phi_i<l_i\}$ and
$l_i\le\bar l_i+\epsilon$. With uniformly bounded derivative majorants and the boundary-point
condition above, the procedure terminates for every $\epsilon>0$.
\end{proposition}
\begin{proof}
Every boundary point lies in a discarded piece, where $\phi_i\le\mathrm{UB}<\ell+\epsilon$ with $\ell$
nondecreasing, and $\ell\le\bar l_i$ because it is attained. For mesh triangles and polygon
edges, the center is a boundary evaluation, so $\phi_i(\mathbf c)\le\ell$.
For an octree box, its boundary point $\mathbf b$ gives
$\phi_i(\mathbf c)\le\phi_i(\mathbf b)+GC\rho\le\ell+GC\rho$, where $G$ is a uniform gradient bound
on the relevant neighborhood. Hence in both cases
$\mathrm{UB}-\ell\le G(1+C)\rho+\tfrac12\kappa\rho^2$, taking $C=0$ for on-boundary centers.
Uniform subdivision makes this bound smaller than $\epsilon$ after finitely many levels.
Each level has finitely many pieces, which proves termination.
\end{proof}
The bounds certify the whole continuous boundary, with tightness up to $\epsilon$ under the
conditions above. Boundary evaluations supply lower bounds; Taylor bounds control the points
between them. The fit affects the resulting level, while the enclosure inequality follows from
the certificate, so the training targets may be approximate. We set $\mathcal B_i=\{\phi_i\le l_i\}$ and write
$\mathcal S_i=\{\phi_i=l_i\}$ for the \emph{certified surface}.

% -----------------------------------------------------------------------------
\section{Surface-Cover Barriers}\label{sec:cover}
We keep one SDF per body and enforce the first-contact condition on the whole certified surface; the
construction is the same for planar and spatial bodies, rigid or articulated. Let body $A$ (e.g., a
robot or a link) have the certified surface $\mathcal S_A=\{\phi_A=l_A\}$ in its frame, and let body $B$
(e.g., another robot or an obstacle) have the field $\phi_B$ and level $l_B$ in its frame, both from
Proposition~\ref{prop:level}. We write $\mathbf y(\mathbf x)$ for a point $\mathbf x$ of $A$ mapped into
$B$'s frame by the current relative pose. For planar bodies, the fields are bivariate and the boxes
below are squares.

\subsection{Safe Set on the Continuous Surface}
The safe set is the first-contact condition on every point of the certified surface,
\begin{equation}\label{eq:safe}
\mathcal H_{AB}=\big\{\text{poses}\ \big|\ h(\mathbf x)=\phi_B(\mathbf y(\mathbf x))-l_B\ge0\ \ \forall\mathbf x\in\mathcal S_A\big\};
\end{equation}
no box, sample, or bound enters its definition. We preserve this continuous spline-defined set
and construct finite Bernstein certificates for its admissible controls. Surface separation alone
does not exclude pre-existing containment of one body inside the other; collision avoidance uses
the initialization in the following theorem.
\begin{theorem}\label{thm:cover}
Suppose $\{\phi_A\le l_A\}$ and $\{\phi_B\le l_B\}$ are compact subsets of their domains and
the strict boundary certificates of Proposition~\ref{prop:level} hold. Suppose initially
$\mathcal G_A\cap\mathcal G_B=\emptyset$ and
$\{\phi_A\le l_A\}\cap\{\phi_B<l_B\}=\emptyset$.
If a continuous trajectory stays in $\mathcal H_{AB}$, then $\mathcal G_A$ and $\mathcal G_B$
never intersect. All set intersections are understood in a common spatial frame.
\end{theorem}
\begin{proof}
For $\epsilon>0$, if $\{\phi_A\le l_A\}$ met $\{\phi_B\le l_B-\epsilon\}$, Lemma~\ref{lem:contact} would
give a first point on both boundaries, where $\phi_A=l_A$ and $\phi_B=l_B-\epsilon$, contradicting
$h\ge0$ on $\mathcal S_A$. Hence $\{\phi_A\le l_A\}\cap\{\phi_B<l_B\}=\emptyset$ throughout. A first
contact of $\mathcal G_A$ and $\mathcal G_B$ would lie on $\partial\mathcal G_A\cap\partial\mathcal G_B$
(Lemma~\ref{lem:contact}), where $\phi_A<l_A$ and $\phi_B<l_B$ by Proposition~\ref{prop:level}.
\end{proof}

\subsection{Certified Cover and Lifted Barrier}
On each cell of $\phi_A$, the Bernstein coefficients of $\phi_A-l_A$ follow from~\eqref{eq:b2b}. If
they all have one strict sign on a box, the box contains no point of $\mathcal S_A$ by the
convex-hull property. The converse need not hold: a retained box may contain no surface point.
The default retains each native SDF cell whose coefficient interval contains zero, without
subdividing it. These patches are the boxes $\{\mathcal Q_k\}$, with side $s=h_A$ fixed by the SDF
grid, centers $\mathbf c_k$, and half-diagonal $r=\sqrt n\,s/2$; their union contains all of
$\mathcal S_A$. There is no independent cover-size parameter. Optionally, retained patches can be
halved using~\eqref{eq:subdiv}, discarding children whose coefficient intervals exclude zero.
The cover is computed once per body; it only partitions
the surface for the enforcement below and does not enter~\eqref{eq:safe}. Subdivision restricts the
same polynomial to smaller domains; it neither samples nor refits the surface.

By the comparison lemma, $\mathcal H_{AB}$ stays invariant if $\dot h(\mathbf x)\ge-\gamma h(\mathbf x)$
at every surface point. The surface is implicit, so we impose this condition on whole boxes, after
lifting the barrier to reduce the tightening caused by points off the surface: for $\mathbf x$ in a
box,
\begin{equation}\label{eq:lift}
\tilde h(\mathbf x)=\phi_A(\mathbf x)-l_A+\phi_B(\mathbf y(\mathbf x))-l_B ,
\end{equation}
which equals $h$ on $\mathcal S_A$. Since $\mathbf x$ is fixed in $A$'s frame, $\phi_A(\mathbf x)$ does
not change, and $\dot{\tilde h}(\mathbf x)=\nabla\phi_B(\mathbf y)^\top\dot{\mathbf y}(\mathbf x)$. A
box point displaced from the surface by $\delta$ toward $B$ has $\phi_A-l_A\approx\delta$ and
$\phi_B-l_B\approx h-\delta$ for SDF-like fields, so $\tilde h\approx h$: the lifting cancels the
first-order effect only to the extent that the opposing gradients cancel. The box points move with
$\dot{\mathbf y}(\mathbf x)=\mathbf J\mathbf u+\mathbf W(\mathbf u)(\mathbf y-\mathbf y_{\mathbf c})$,
affine in the inputs $\mathbf u$, where $\mathbf J\mathbf u$ is the velocity of the box center in
$B$'s frame and $\mathbf W(\mathbf u)$ the relative angular velocity as a skew matrix
(Sec.~\ref{sec:cover-der}).

\subsection{Bernstein and Vertex Constraints}
A box in $A$'s frame can cross several knot cells of $B$ after a rigid transformation.
Certified Taylor models avoid constructing their intersections online. These models are used
only to bound the control inequality; the represented surfaces and certified levels remain
those of the original splines.
Cellwise bounds on the Frobenius norm of the third-derivative tensor follow from the Bernstein
coefficients of the third partial derivatives. A tensor-product cubic B-spline is $C^2$ with bounded,
piecewise-polynomial third
derivatives, so its Hessian is Lipschitz with these bounds on any region they cover, which is what the
third-order Taylor remainders below require. On each (fixed) box of $A$, $M_A$ is their maximum over
all cells meeting the bounding cube of half-width $r$ around its center, including for a full
native patch. For $\phi_B$, the bound must vary continuously with the pose:
$M_B$ is piecewise multilinear on a grid of spacing $h/2$, with node values the maxima over the cells meeting
the cube of half-width $h/2+r$ around the node; interpolation is a convex combination of nodes whose
cubes all contain the ball of radius $r$ around $\mathbf y_{\mathbf c}$, so $M_B$ bounds the tensor on
that ball and is Lipschitz. Let $q_A$ and $q_B$ be the second-order Taylor polynomials of
$\phi_A$ at $\mathbf c$ and of $\phi_B$ at $\mathbf y_{\mathbf c}$. Taylor's theorem gives, on a box,
$\phi_A\ge q_A-\tfrac16M_Ar^3$, $\phi_B\ge q_B-\tfrac16M_Br^3$, and
$\norm{\nabla\phi_B-\nabla q_B}\le\tfrac12M_Br^2$. The box points move relative to $B$ with the twist
$(\mathbf V,\boldsymbol\Omega)$ of $A$ relative to $B$ at $A$'s origin $\mathbf p_A$, linear in $\mathbf u$,
so $\norm{\dot{\mathbf y}(\mathbf x)}\le\norm{\mathbf V}_1+\rho\norm{\boldsymbol\Omega}_1$ with
$\rho=\norm{\mathbf c}+r$, since $\mathbf c$ is expressed in $A$'s body frame.
With $\mathbf d=\mathbf y-\mathbf y_{\mathbf c}$, $\norm{\mathbf d}\le r$,
$\mathbf g=\nabla\phi_B(\mathbf y_{\mathbf c})$, $\mathbf H=\nabla^2\phi_B(\mathbf y_{\mathbf c})$, and
$\dot{\mathbf y}_{\mathbf c}=\mathbf J\mathbf u$, the models give on the box
\begin{equation}\label{eq:psi-lb}
\begin{aligned}
\Psi\ge{}&\mathbf g^\top\dot{\mathbf y}_{\mathbf c}+\mathbf d^\top\mathbf H\dot{\mathbf y}_{\mathbf c}
+\mathbf g^\top\mathbf W\mathbf d-r^2\norm{\mathbf H}\norm{\boldsymbol\Omega}\\
&+\gamma\big(q_A-l_A+q_B-l_B-\tfrac16(M_A+M_B)r^3\big)\\
&-\tfrac12M_Br^2\big(\norm{\mathbf V}_1+\rho\norm{\boldsymbol\Omega}_1\big),
\end{aligned}
\end{equation}
where $\Psi(\mathbf x,\mathbf u)=\nabla\phi_B(\mathbf y)^\top\dot{\mathbf y}(\mathbf x)+\gamma\tilde h(\mathbf x)$;
the first three terms are affine in $\mathbf d$ and in $\mathbf u$.

\emph{Joint coefficient certificate.} Before bounding the terms separately as
in~\eqref{eq:psi-lb}, one can retain the full quadratic model. More generally, write
$g_A=\phi_A-l_A$ and $F=\dot h+\gamma h$. For any scalar $w$ held constant over a box,
$F+\gamma w g_A=F$ on $\mathcal S_A$. Define
\begin{equation}\label{eq:joint-model}
\begin{aligned}
p_w={}&(\mathbf g+\mathbf H\mathbf d)^\top
              (\mathbf J\mathbf u+\mathbf W(\mathbf u)\mathbf d)\\
 &+\gamma\big(w(q_A-l_A)+q_B-l_B\big),\\
E_w={}&\tfrac\gamma6(|w|M_A+M_B)r^3\\
 &+\tfrac12M_Br^2(\mathbf1^\top\mathbf a_V+
                       \rho\mathbf1^\top\mathbf a_\Omega).
\end{aligned}
\end{equation}
The two-sided Taylor bounds give $F+\gamma w g_A\ge p_w-E_w$.
Consequently, if $b_k(p_w)$ are its $3^n$ tensor Bernstein coefficients,
\begin{equation}\label{eq:joint-coef}
b_k(p_w)\ge E_w\quad\text{for every }k,\qquad
\mathbf a\ge|(\mathbf V,\boldsymbol\Omega)|,
\end{equation}
then the CBF condition holds on $\mathcal S_A\cap\mathcal Q$. This certificate keeps
the value and velocity terms in the same polynomial, preserving their cancellation.
For fixed $w$, its constraints are affine in $(\mathbf u,\mathbf a)$.
The unit choice $w=1$ recovers the summed-field residual. Other choices are equality
multipliers for the surface condition, not new barriers to be differentiated;
they change neither the surface nor~\eqref{eq:safe}.

\emph{Vertex reduction.} The implementation used in the main experiments reduces
the number of rows by bounding the value and velocity terms separately:
\begin{proposition}\label{prop:coef}
Let $\beta_k$, $k=1,\dots,3^n$, be the Bernstein coefficients of degree two per coordinate of
$q_A-l_A+q_B(\mathbf y(\cdot))-l_B$ on a box $\mathcal Q$, let $\mathbf d_j$, $j=1,\dots,2^n$, be the
vertices of $\mathcal Q$ relative to its center in $B$'s frame, and let
$\sigma_V=\frac12M_Br^2$ and $\sigma_\Omega=\frac12\rho M_Br^2+r^2\norm{\mathbf H}$. If, for all $j$,
\begin{equation}\label{eq:coef}
\begin{aligned}
&(\mathbf g+\mathbf H\mathbf d_j)^\top\mathbf J\mathbf u+\mathbf g^\top\mathbf W(\mathbf u)\mathbf d_j
-\sigma_V\mathbf 1^{\!\top}\mathbf a_V-\sigma_\Omega\mathbf 1^{\!\top}\mathbf a_\Omega\\
&\quad+\gamma\big(\min_k\beta_k-\tfrac16(M_A+M_B)r^3\big)\ge0,\\
&\mathbf a=(\mathbf a_V,\mathbf a_\Omega)\ge|(\mathbf V,\boldsymbol\Omega)|,
\end{aligned}
\end{equation}
then $\Psi(\mathbf x,\mathbf u)\ge0$ for every $\mathbf x\in\mathcal Q$.
\end{proposition}
\begin{proof}
In~\eqref{eq:psi-lb}, the terms affine in $\mathbf d$ attain their minimum over the (rotated) box at a
vertex, and the value model is at least $\min_k\beta_k$ by the convex-hull property; the twist terms
are bounded by $\mathbf a$.
\end{proof}
The constraints, $2^n$ per box, are affine in $(\mathbf u,\mathbf a)$. The epigraph variables $\mathbf a$ (six per pair of
bodies in space, three in the plane) are shared by all boxes of the pair, and since the left side
of~\eqref{eq:coef} decreases in $\mathbf a$, the feasible set in $\mathbf u$ is exactly that of the same
constraints with $|(\mathbf V,\boldsymbol\Omega)|$; the filter remains a QP. Joint-by-joint speed bounds
would be looser, since joint velocities partly cancel in the twist.
The remainder terms tighten only the admissible inputs, and the one
proportional to the twist vanishes at rest. Online, $\phi_B$, its first two derivatives, and $M_B$ are
evaluated at box centers; the Taylor data of $\phi_A$ are computed offline. At rest, the loss of the
constraints is bounded as follows:
\begin{proposition}\label{prop:tight}
With $\mathbf b=\nabla\phi_A(\mathbf c)+\Rot^\top\nabla\phi_B(\mathbf y_{\mathbf c})$ and
$\mathbf Q=\nabla^2\phi_A(\mathbf c)+\Rot^\top\nabla^2\phi_B(\mathbf y_{\mathbf c})\Rot$, where
$\Rot$ maps $A$'s frame to $B$'s, the constant terms of~\eqref{eq:coef} are at least
$\gamma\big(\tilde h(\mathbf c)-\tfrac s2\norm{\mathbf b}_1-\tfrac{s^2}8\sum_{ij}|Q_{ij}|
-\tfrac16(M_A+M_B)r^3\big)$.
\end{proposition}
\begin{proof}
With $\mathbf x=\mathbf c+s\boldsymbol\tau$, $\boldsymbol\tau\in[-\frac12,\frac12]^n$, the degree-two
coefficients of $\tau_i$, $\tau_i\tau_j$, and $\tau_i^2$ lie in $[-\frac12,\frac12]$,
$[-\frac14,\frac14]$, and $[-\frac14,\frac14]$.
\end{proof}
The first-order term vanishes only where the two gradients cancel, $\nabla\phi_A=-\Rot^\top\nabla\phi_B$:
exact SDFs satisfy this at a smooth external tangency where both gradients exist; fitted fields
need not satisfy it, since
their gradients may differ in norm and in direction. A bound on $\phi_B$ alone over a box loses
$\norm{\nabla\phi_B}r$ even there. The proposition concerns the constant terms at rest; the velocity
terms of~\eqref{eq:coef} add their own, input-dependent loss.

\subsection{Derivatives and Pruning}\label{sec:cover-der}
For planar $SE(2)$ single integrators, a point of $A$ at world position $\mathbf w$ moves in $B$'s frame
with $\Rot(\theta_B)^{\!\top}\big(\mathbf v_A+\omega_AJ(\mathbf w-\mathbf p_A)-\mathbf v_B-\omega_BJ(\mathbf w-\mathbf p_B)\big)$,
$J=\left[\begin{smallmatrix}0&-1\\1&0\end{smallmatrix}\right]$, and
$\mathbf W=(\omega_A-\omega_B)J$. For a link $A$ of an articulated robot, a revolute joint $j$ before
the link contributes the column $\mathbf z_j\times(\mathbf w-\mathbf o_j)$ and the angular velocity
$\mathbf z_j\dot q_j$; for two moving arms, the joints of $B$'s arm enter with the opposite sign. In
every case, $\mathbf J$ and $\mathbf W$ are linear in the inputs and are obtained directly from
the kinematics, without a closest-point optimization.

\emph{Pruning.} A box needs no constraint if the constant terms of~\eqref{eq:coef}, divided by
$\gamma$, are at least an activation threshold $\eta>0$: then $h\ge\eta$ on its part of the surface.
We group the boxes of each body into clusters with center $\bar{\mathbf c}$ and radius $R$ (covering
the boxes' balls) and skip a cluster if $\phi_B(\bar{\mathbf c})-GR-l_B+\min(\phi_A-l_A)\ge\eta$,
where $G$ bounds $\norm{\nabla\phi_B}$ on the cells the cluster can reach and the minimum of
$\phi_A-l_A$ over the cluster's boxes is certified offline; remaining boxes are first screened with
$\phi_B(\mathbf y_{\mathbf c})-\norm{\nabla\phi_B(\mathbf y_{\mathbf c})}r-\frac12\kappa_B(\mathbf y_{\mathbf c})r^2-l_B+\min_{\mathcal Q}(\phi_A-l_A)\ge\eta$.
If $r>h_B$, the Hessian majorant uses the same radius-dependent window construction as $M_B$.

\emph{Optional refinement.} The default uses native patches and performs no online refinement.
Where the field varies rapidly, $\frac12M_Br^2$ can exceed
$\norm{\nabla\phi_B}$: the bound then fails to resolve the gradient direction and can reject useful
motion, causing a stall. An active box
with $\frac12M_Br^2>\vartheta$ is therefore replaced online by those of its $2^n$ half-size boxes whose
bounds of $\phi_A-l_A$ bracket zero; they cover the same part of $\mathcal S_A$.
For unchanged derivative bounds, halving $r$ quarters the gradient remainder and reduces the
value remainder by a factor of eight. The subdivided implementation uses $\vartheta=0.1$ and at most two halvings.
The original experiments additionally limit refinement to boxes whose lower bound of $\tilde h$
is below 1~cm, except for the tube scene; the common distance-independent rule is evaluated
in the certificate comparison below.
\begin{theorem}\label{thm:invariance}
Let a locally Lipschitz feedback satisfy~\eqref{eq:coef} or~\eqref{eq:joint-coef} on every box, or on the kept half-size boxes of
a refined box, whose certified lower bound of $\tilde h$ is below $\eta$. Then $\mathcal H_{AB}$ is
forward invariant.
\end{theorem}
\begin{proof}
Fix $\mathbf x\in\mathcal S_A$ with $h(\mathbf x)\ge0$ initially. Whenever $h(\mathbf x)<\eta$, every
box containing $\mathbf x$ has a lower bound of $\tilde h$ below $\tilde h(\mathbf x)=h(\mathbf x)<\eta$,
so it is constrained. Proposition~\ref{prop:coef}, or the joint certificate with $g_A(\mathbf x)=0$, gives
$\dot h(\mathbf x)\ge-\gamma h(\mathbf x)$. Hence
$h(\mathbf x)$ cannot cross zero.
\end{proof}

\begin{proposition}[Certification of strictly feasible inputs]\label{prop:complete}
Fix a pose and input $\mathbf u$. Suppose $\mathcal S_A$ is compact, its mapped neighborhood lies
inside the field domains, and $F(\mathbf x,\mathbf u)\ge\epsilon>0$ on $\mathcal S_A$.
Assume uniformly bounded derivative majorants on that neighborhood. If subdivision is allowed
to arbitrary depth, with consistent exclusion bounds for $g_A$, then a finite subdivision
certifies this input using either~\eqref{eq:coef} or~\eqref{eq:joint-coef} with $w=1$.
\end{proposition}
\begin{proof}
The unit lifted residual equals $F$ on the surface. Continuity and compactness give a neighborhood
where it is at least $\epsilon/2$. As the maximum box diameter tends to zero, consistent exclusion
bounds discard boxes bounded away from the surface. The remaining boxes eventually lie in that
neighborhood. Bernstein enclosure errors, the losses from separating the terms, and the Taylor
remainders converge uniformly to zero for the fixed input. Eventually their sum is smaller than
$\epsilon/2$, so every retained box passes its certificate.
\end{proof}
This result concerns optional subdivision and strictly feasible inputs, not completeness of the
fixed native-patch certificate or finite certification of inputs with $F=0$.
The two-level subdivided implementation has a finite computational budget. Re-centering Taylor models can
also prevent successive finite certificates from being nested.

% -----------------------------------------------------------------------------
\section{Safety Filter}\label{sec:filter}
Collecting the constraints~\eqref{eq:coef} of the active boxes of all pairs, the filter is
\begin{equation}\label{eq:qp}
\begin{aligned}
(\mathbf u^\ast,\mathbf a^\ast)=\argmin_{\mathbf u\in\mathcal U,\,\mathbf a}\ &\norm{\mathbf u-\mathbf u^{\rm nom}}^2+\varepsilon_a\norm{\mathbf a}^2\\
\text{s.t. }&\eqref{eq:coef}\ \text{for all active boxes},
\end{aligned}
\end{equation}
with $\mathbf u$ the stacked planar velocities or the joint velocities, $\mathbf a$ the epigraph
variables of all active pairs, and $\varepsilon_a=10^{-3}$. With
polytopic input bounds, \eqref{eq:qp} is a QP, which we reformulate as a least-distance program and solve numerically by
non-negative least squares~\cite{lawson1995solving}, warm-started by constraint generation.
\begin{proposition}[Feasibility]\label{prop:feas}
For driftless dynamics---planar single integrators, unicycles, or velocity-controlled joints---and
any $\mathcal U\ni\mathbf 0$, $(\mathbf u,\mathbf a)=\mathbf 0$ is feasible for~\eqref{eq:qp} whenever
the constant terms of~\eqref{eq:coef} are nonnegative on all active boxes; by
Proposition~\ref{prop:tight}, this holds wherever $\tilde h$ at the box centers exceeds the bound given
there, which shrinks with the box size.
\end{proposition}
\begin{proof}
All velocity terms of~\eqref{eq:coef} are linear in $\mathbf u$ without drift, so
$(\mathbf u,\mathbf a)=\mathbf 0$ leaves the constant terms.
\end{proof}
\begin{corollary}\label{cor:closed}
Suppose each pair satisfies the initial conditions of Theorem~\ref{thm:cover},
\eqref{eq:qp} stays feasible and enforces the active-box conditions of Theorem~\ref{thm:invariance},
and its minimizer is locally Lipschitz. Then, along the resulting continuous-time trajectory,
no pair of ground-truth shapes intersects. The certified sublevel regions may touch at their
level boundaries.
\end{corollary}
\begin{proof}
Theorem~\ref{thm:invariance} keeps each pose in $\mathcal H_{AB}$, and Theorem~\ref{thm:cover}
concludes.
\end{proof}
Stopping is certified when the constant terms are nonnegative. A finite certificate may fail to
certify stopping near the boundary even though zero input preserves the exact safe set for the
driftless model. Proposition~\ref{prop:tight} bounds the constant-term loss; it does not give a
general physical width of this region without additional assumptions on the fields.
No step needed slack in the reported subdivided experiments, and the safe set itself is unaffected.
Feasibility does not imply progress: exact non-convex geometry can hook and block
bodies~\cite{reis2020undesirable}, and we do not resolve such deadlocks. In the experiments, a slack
variable is added only if~\eqref{eq:qp} is infeasible.

% -----------------------------------------------------------------------------
