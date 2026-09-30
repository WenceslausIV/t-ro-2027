# 이론 검토: 건전성·견고성·우아한 업그레이드 (2026-09-30, Claude)

범위: `main.tex`의 이론 부분(II–VI절, 증명 전부)만 본다. 수치·코드·결과 대조는 이 문서의 범위가 아니다.
결론 먼저: **치명적 오류는 없다.** 모든 보조정리·명제·정리의 증명을 한 줄씩 다시 유도했고, 성립한다.
다만 (1) 필요 이상으로 강한 가정, (2) 부정확한 서술, (3) 기호 충돌이 있었고, (4) 이론을 한 단계 끌어올릴
수 있는 지점이 여러 곳 있다. (1)–(3)은 `main.tex`에 반영했고(§3), (4) 중 방법·실험을 바꾸는 것은
결정이 필요하므로 LaTeX 초안과 함께 제안한다(§4).

---

## 1. 결과별 건전성 판정

| 결과 | 판정 | 확인한 핵심 |
|---|---|---|
| Lemma 1 (first contact) | ✅ | 교차 시각 집합은 닫힘(거리함수 연속). 내부점 z의 물체 좌표 역상이 열린 공 안에 머무름 → 강체운동 연속성으로 모순. 볼록성·매끄러움 불필요. |
| Prop. 1 (certified level) | ✅ | κ 스플라인: 점 c가 span k에 있을 때 활성 제어값 m=k..k+3 각각의 창 [m−4,m+1]이 모두 [k−1,k+1]을 포함 → 반경 ρ≤Δ 공을 덮음. UB = Taylor 2차 + κ. 종료 증명의 UB−ℓ ≤ G(1+C)ρ + κρ²/2 정확. 작업 한도 시 강한 포함 유지. |
| Thm 1 (cover ⇒ 비충돌) | ✅ | ∂{φ≤l} ⊂ {φ=l} (연속성). ε-축소 후 Lemma 1 → 모순. 두 번째 단계 정확. **한 방향(S_A vs φ_B)만으로 충분**함도 정확. |
| (psi-lb) | ✅ | ∇φ_B(y)=g+Hd+e, ‖e‖≤½M_B r², dᵀHWd ≥ −r²‖H‖‖Ω‖ (‖[Ω]×‖₂=‖Ω‖), ‖ẏ‖≤‖V‖₁+r̄‖Ω‖₁. 3차 Taylor 잔차(Hessian Lipschitz, Frobenius ≥ 연산자 노름). |
| M_B 다중선형 보간 | ✅ | 노드와 y_c의 ∞-거리 ≤ Δ/2 → 노드 창(반폭 Δ/2+r)이 y_c의 r-정육면체 포함 → 볼록결합도 상한. |
| Joint certificate | ✅ | F+γwg_A ≥ p_w − E_w; p_w는 x에 대해 전차수 2 → 3ⁿ Bernstein 계수가 정확; S_A 위에서 g_A=0. |
| Prop. 2 (vertex) | ✅ | 속도항은 d에 대해 affine → 회전된 상자의 꼭짓점에서 최소. 값/속도 분리는 보수적(최소의 합 ≤ 합의 최소). |
| Prop. 3 (tight) | ✅ | τ, τᵢτⱼ, τᵢ²의 2차 Bernstein 계수가 각각 [−½,½], [−¼,¼], [−¼,¼] → 식 정확. |
| Pruning (cluster/box/collar) | ✅ | min(φ_A−l_A)를 더하는 것은 보수적(표면을 포함하는 상자에서 ≤0). collar 논증(공이 Ω 밖 → 모든 점이 2r collar 또는 Ω 밖) 정확. |
| Thm 2 (invariance) | ✅ | 점별 AC 논증. **가정이 필요 이상으로 강했음 → §3.1** |
| Prop. 4 (complete) | ✅ | Bernstein 범위 → 실제 범위 수렴, 재유지 상자 → S_A 근방, 잔차 → 0. |
| Prop. 5 / Cor. 1 | ✅ | driftless이면 (u,a)=0에서 속도·잔차항이 모두 0. |

**견고성 측면의 장점(논문에 더 강조할 가치가 있음):** 적합(fitting) 품질은 보수성에만 영향을 주고
건전성에는 영향을 주지 않는다. 내부 bump(G ⊄ B), 먼 곳의 가짜 sublevel 섬, 경사 노름 ≠ 1 등은 모두
보수성(혹은 가짜 장애물)으로만 나타난다. 경계 포함 인증(Prop. 1)과 첫 접촉 논증 덕분이다. → 본문 Thm 1
뒤에 한 문장으로 반영(§3.2).

---

## 2. 발견한 약점 (중요도 순)

1. **[중] Thm 2의 활성 임계값 η₀>0 가정은 불필요.** 증명은 h_x<0인 구간 (t₀,t₁]에서만 부등식을 쓴다.
   그 구간에서는 h≥0 인증을 가진 상자가 x를 포함할 수 없으므로 η₀=0(“h≥0 인증”)이면 충분하다.
   따라서 `η₀=min(η,δ_B)`라는 부수 장치가 통째로 사라진다. η의 진짜 역할은 연속시간 증명이 아니라
   샘플 사이 운동의 버퍼다(§4.2에서 정량화).
2. **[중] Corollary 1이 QP 실현가능성을 가정 — 사실 불필요.** driftless 시스템에서는 u=0이면 자세가
   고정되어 모든 h_x가 상수다. 연속시간 증명은 “h_x<0일 때 ḣ_x≥0”만 쓰므로, QP가 불가능할 때 slack 대신
   **u=0으로 fallback하면 실현가능성 가정 없이 안전이 보장**된다. 현재 실험은 slack을 쓰므로 slack이 쓰인
   구간(native 12-mm 2회, KKT 7회)에서만 보장이 끊겼던 것이다.
3. **[중] “fitted field의 경사는 노름과 방향이 모두 다를 수 있다”(Prop. 3 뒤)는 접점에서는 틀림.** S_A 위
   h의 최소점(수준집합의 외접점)에서 라그랑주 조건 Rᵀ∇φ_B = −λ∇φ_A (λ=‖∇φ_B‖/‖∇φ_A‖)가 **적합된 장에서도**
   성립한다. 즉 방향은 항상 반평행이고, 노름 불일치 (1−λ)만 1차항에 남는다. 그리고 이 λ가 바로 논문의
   w_N(및 접점에서 w_P)이다 → 고정 multiplier 실험의 “개선 없음”이 이론적으로 설명된다(적합 장의 λ≈1).
4. **[소] M_B “must vary continuously”는 불필요.** 불변성 증명은 시각별 부등식만 쓴다. 연속성은
   폐루프 정칙성(해의 존재)용일 뿐 안전 증명에는 쓰이지 않는다. κ의 C² 성질도 마찬가지.
5. **[소] 사용되지 않는 Lemma (lem:cbf)** 가 “locally Lipschitz feedback”을 가정 → 구현 필터가 만족하지
   못하는 가정을 서론에 굳이 제시하는 셈. Thm 2 증명은 Prop. 2를 인용하지 않았음.
6. **[소] 기호 충돌:** M_B(Bernstein→power 행렬 vs 3차 도함수 상한), J(속도 Jacobian 𝐉 vs 평면 회전생성자 J
   vs 입력 변형 적분 J), h(barrier vs 셀 크기 h_A), ρ(Prop. 1 반경 vs ‖c‖+r), Q(B-spline 행렬 vs 상자 𝒬 vs
   Hessian 합 𝐐), a/b(epigraph 𝐚 vs 비교절의 경사 𝐚,𝐛), s(Lemma 1 매개변수 vs 상자 변).
7. **[구조적 한계] 보장이 연속시간 한정 + 폐루프 해 존재 미증명.** 구현은 10-ms ZOH. 두 한계는 §4.2의
   sampled-data 정리 하나로 동시에 해소된다(ZOH이면 폐루프가 차분방정식이라 존재 문제가 사라짐).
8. **[구조적 한계] 보수성이 1차(O(s))로만 정량화됨(Prop. 3).** multiplier를 결정변수로 두면 2차 O(s²)가
   된다(§4.1). native patch(s=Δ_A)의 손실이 장 해상도의 제곱으로 줄어든다는 뜻이라 “native가 기본값”이라는
   설계를 이론적으로 뒷받침한다.

---

## 3. `main.tex`에 반영한 수정 (이론 텍스트만, 방법·실험 불변)

1. **Preliminaries:** 미사용 Lemma(lem:cbf) 삭제 → 한 문장(“고전 결과는 Lipschitz feedback 가정;
   Thm 2가 여기서 쓰는 형태”). Bernstein→power 행렬 M_B → Λ, 부분구간 [ξ₀, ξ₁]로 표기.
2. **Lemma 1:** 매개변수 s → t.
3. **Certified regions:** 셀 크기 h → Δ_i, κ 설명 단순화(“인접 셀 상한; 구현은 C² 스플라인으로 평활”).
   내부 bump 서술 2문장 → 1문장. primitive 합집합 단락 압축.
4. **Safe set (IV-A):** Danskin 정리로 재해석 추가 — H_AB = {min_{S_A} h ≥ 0}; closest-point barrier는 한
   최소점만 제약하므로 최소점이 유일할 때만 충분(L1, L2의 엄밀한 원인); 모든 x에 부과하면 선택 없이 모든
   최소점을 덮는다. `danskin1966theory` 인용 추가(SIAM J. Appl. Math. 14(4):641–664, 1966).
5. **Thm 1:** 결론에 {φ_A≤l_A}∩{φ_B<l_B}=∅ 유지를 명시, 그리고 역방향 한 문장 → “연속 운동에서 H_AB는
   정확히 인증된 비겹침 조건이며, 내부 적합 오차는 건전성에 영향 없음”.
6. **Joint certificate:** w는 상자마다·시간마다 달라도 됨; F+γwg_A = 표면 조건의 라그랑지안, w = 등식 multiplier.
   “summed-field residual”(코드 용어) → “lifted barrier (lift)”.
7. **Prop. 3 뒤 단락:** 접점에서의 라그랑주 조건으로 교체(방향은 항상 반평행, 노름 불일치만 남고 w=λ로 제거).
8. **Pruning / domain collar:** “Thm 2는 생략 상자에 h≥0만 요구; η는 샘플 간 버퍼”; η₀ 문장 삭제.
9. **Thm 2:** 가정 최소화(h≥0 인증, driftless & u=0 허용), 증명에서 Prop. 2 인용, h_x<0이면 y가 φ_B가 C²인
   영역에 있음을 명시. 증명 뒤 한 문장: 연속시간에서는 ḣ≥0 (h<0일 때)만 필요, γ는 샘플링에서 의미.
10. **Safety filter:** QP 제약에 (joint-coef) 포함. Cor. 1을 “feasible QP 해 또는 u=0”으로 일반화.
    u=0 fallback 문단 추가(실험은 slack을 썼고 모든 slack step을 보고한다고 명시). Limitations 반영.
11. **기호:** 평면 J → 𝐄=R(π/2); ρ(IV절) → r̄=‖c‖+r; Prop. 3 행렬 𝐐 → 𝐍; 비교절 𝐚,𝐛 → 𝛎_A,𝛎_B; 입력 변형 J → D.
12. **Certificate comparison:** w_N, w_P가 외접점에서 모두 λ와 같음 → “개선 없음”은 λ≈1로 예상된 결과라고 명시.

빌드: pdfLaTeX(Times), 경고·미정의 참조 없음, **13쪽 유지**(수정 전과 동일).

---

## 4. 우아한 업그레이드 제안 (결정 필요)

### 4.1 [가장 추천] 자유 multiplier = 라그랑지안 인증서, 그리고 2차 정확도 정리

**관찰.** (joint-coef)의 Bernstein 계수는 w에 대해서도 affine이다(g_A는 입력과 무관하므로 w·g_A 항이
bilinear가 되지 않는다). 따라서 상자마다 w_Q를 **QP 결정변수**로 두어도 필터는 QP로 남는다
(|w|는 epigraph 하나 추가, 또는 외접 상황에 맞게 w_Q≥0로 제한).

- 실현가능 입력집합 = ⋃_w (고정 w 집합)의 사영 → 볼록, 그리고 **모든 고정 w(특히 w=1)를 지배**.
- Thm 2는 그대로 성립(w는 시간에 따라 바뀌어도 됨 — 이번에 본문에 명시해 둠).
- 고정 w_N, w_P가 못 하는 일 — 입력에 의존하는 CBF 조건의 **법선 방향 미분 상쇄** — 를 QP가 자동으로 한다.
  실험에서 관찰된 “접촉 근처에서 미끄러지는 링크가 느려짐”의 1차 원인(곡률×속도 항)을 제거한다.
- 비용: 활성 상자당 변수 1개. 대안: 클러스터당 공유 w, 또는 직전 step의 최적 w 사용(변수 추가 없음, 여전히 건전).

**새 명제(2차 정확도).** 초안:

```latex
\begin{proposition}[Second-order accuracy]\label{prop:second}
Suppose $\norm{\nabla\phi_A}\ge\mu>0$ on the retained boxes and the derivative majorants of
$\phi_A$ and $\phi_B$ up to third order are bounded on the relevant neighborhoods. For a pose and
an input $\mathbf u$, let $F^\ast(\mathbf u)=\min_{\mathbf x\in\mathcal S_A}F(\mathbf x,\mathbf u)$.
There are constants $s_0,C_0,C_1$, independent of the box size, such that for boxes of side
$s\le s_0$ and $F^\ast(\mathbf u)\ge(C_0+C_1\norm{\mathbf u})s^2$, every retained box satisfies
\eqref{eq:joint-coef} for some $w$ with $|w|\le C_2(1+\norm{\mathbf u})$.
\end{proposition}
\begin{proof}[Proof sketch]
At the box center, $\nabla_{\mathbf x}p_w=\nabla_{\mathbf x}F+\gamma w\nabla g_A$. The choice
$w=-\hat{\mathbf n}^\top\nabla_{\mathbf x}F/(\gamma\norm{\nabla g_A})$, $\hat{\mathbf n}=\nabla g_A/\norm{\nabla g_A}$,
leaves the tangential part $\mathbf t$, so the degree-two Bernstein coefficients exceed
$p_w(\mathbf c)-\frac s2\norm{\mathbf t}_1-O(s^2)$, and $E_w=O((1+\norm{\mathbf u})s^2)$.
A retained box lies within $O(s)$ of some $\bar{\mathbf x}\in\mathcal S_A$, where
$p_w(\mathbf c)=F(\bar{\mathbf x})-\mathbf t^\top(\bar{\mathbf x}-\mathbf c)+O(s^2)$ and
$\norm{\mathbf t}\le\norm{\nabla_{\mathcal S}F(\bar{\mathbf x})}+O(s)$. Following $\mathcal S_A$ from
$\bar{\mathbf x}$ for a length $cs$ against $\nabla_{\mathcal S}F$ reaches a point where
$F\le F(\bar{\mathbf x})-cs\norm{\nabla_{\mathcal S}F(\bar{\mathbf x})}+O(s^2)$, and $F\ge F^\ast$ there.
Hence every coefficient exceeds $F^\ast-O((1+\norm{\mathbf u})s^2)$.
\end{proof}
```

검증 메모: (i) 재유지 상자는 Bernstein 계수–Greville 점 오차 O(s²) 때문에 S_A에서 O(s²)/μ 이내; (ii) F는
∇φ_B를 포함하므로 C^{1,1}(Hessian Lipschitz)이면 충분, 곡선 논법은 국소적(음함수 정리)이라 전역 측지선 불필요;
(iii) 모든 상수는 ‖u‖에 affine. w=1이면 같은 논법에서 법선 성분 n̂ᵀ(∇F+γ∇g_A)가 남아
**1차 손실 C s(|1−λ| + ‖u‖)** 가 된다 — 현재 구현과 제안의 차이를 정확히 보여 준다.

**해석(논문 서사에 유리):**
- 정지(u=0): min_{S_A} h ≥ C₀s²/γ 이면 정지가 인증된다. “인증 불가 층”의 두께가 O(s²).
- native patch(s=Δ_A)이면 유한 인증서의 손실이 **장 해상도의 제곱**에 비례 → 별도 cover 해상도가
  없는 native 기본값의 이론적 근거. (현재 가장 큰 리뷰 리스크인 “기본값이 약하다”에 대한 이론적 반론.)
- IV절 전체를 “semi-infinite CBF의 다면체 내부 근사” U_s ⊆ U_ex = {u : F(x,u)≥0 ∀x∈S_A}로 재구성하면
  Prop. 4(수렴) + 새 명제(속도 O(s²))가 한 문장으로 정리된다. 사용자의 제약(“보수성은 입력에만”)을 그대로
  수학적 구조로 표현하는 틀이다.

필요 작업: 코드(상자당 w 변수 또는 직전 w 재사용), 최소한 native Franka 30회 재실험. 반영 여부 결정 필요.

### 4.2 [강력 추천] Sampled-data 정리: 이산 구현 자체를 보장

ZOH 입력 u_k를 [t_k, t_k+Δt)에 유지. 고정된 표면점 x에 대해 e(τ)=h_x(t_k+τ)는 구간 내 C²이고
|ë| ≤ (H c₁² + G c₂)‖u_k‖² =: κ‖u_k‖² (H, G: 한 스텝 sweep 영역의 ‖∇²φ_B‖, ‖∇φ_B‖ 상한; c₁, c₂: 고정입력에서
상대 점속도·가속도 상수 — 평면 SE(2)와 관절속도 모두 가속도가 속도의 제곱에 비례).
ė(0) ≥ −γe(0) + ½κΔt‖u_k‖² 이고 γΔt≤1 이면
e(τ) ≥ e(0)(1−γτ) + ½κ‖u_k‖²τ(Δt−τ) ≥ 0 (τ∈[0,Δt]).

```latex
\begin{theorem}[Sampled-data enforcement]\label{thm:sampled}
Let the inputs be held on $[t_k,t_{k+1})$, $t_{k+1}-t_k\le\Delta t$, $\gamma\Delta t\le1$, with
driftless kinematics. Let $\norm{\dot{\mathbf y}}\le c_1\norm{\mathbf u}$ and
$\norm{\ddot{\mathbf y}}\le c_2\norm{\mathbf u}^2$ under a held input, let $G$ and $H$ bound
$\norm{\nabla\phi_B}$ and $\norm{\nabla^2\phi_B}$ on the region swept in one step, and set
$\kappa=Hc_1^2+Gc_2$. Suppose that at every $t_k$, either $\mathbf u_k=\mathbf 0$, or each cover box
has a certificate $h\ge Gc_1\bar u\Delta t$ on its surface portion ($\bar u\ge\norm{\mathbf u}$ on
$\mathcal U$), or satisfies~\eqref{eq:coef} or~\eqref{eq:joint-coef} with right side
$\frac12\kappa\Delta t\norm{\mathbf u_k}^2$ instead of zero. If $h\ge0$ on $\mathcal S_A$ at $t_0$,
then $h\ge0$ on $\mathcal S_A$ for all $t\ge t_0$, between samples included.
\end{theorem}
\begin{proof}
Fix $\mathbf x\in\mathcal S_A$ and $k$, and let $e(\tau)=h_{\mathbf x}(t_k+\tau)$, so
$|\ddot e|\le\kappa\norm{\mathbf u_k}^2$. For an active box,
$e(\tau)\ge e(0)+\dot e(0)\tau-\frac12\kappa\norm{\mathbf u_k}^2\tau^2
\ge e(0)(1-\gamma\tau)+\frac12\kappa\norm{\mathbf u_k}^2\tau(\Delta t-\tau)\ge0$ on $[0,\Delta t]$
whenever $e(0)\ge0$. For a skipped box, $e(\tau)\ge e(0)-Gc_1\bar u\tau\ge0$. For
$\mathbf u_k=\mathbf 0$, $e$ is constant. Induction over $k$ concludes.
\end{proof}
```

왜 우아한가:
- 추가 항이 **‖u‖²에 비례 → 정지 상태에서 0**. Prop. 5(정지 실현가능성)와 u=0 fallback이 그대로 살아 있다.
  → Thm 1 + 이 정리 + fallback = **구현된 이산 제어기에 대해 무조건적(실현가능성·해 존재 가정 없는) 비충돌 보장**.
  현재 Limitations의 “Closed-loop solutions”와 “Discrete time” 두 문단이 통째로 해소된다.
- η의 역할이 정량화된다: η ≥ G c₁ ū Δt. (지금은 “버퍼”라는 정성적 설명뿐.)
- QP 유지: ‖u‖² ≤ ū‖u‖₁ 로 상계하면 epigraph 𝐬≥|𝐮| 하나로 affine. (정확한 볼록 2차 제약으로 두면 SOCP.)
- 주의: 곡률이 큰 곳(H ~ 1/필렛반경)에서 affine 상계는 접근 속도를 상당히 제한할 수 있다. 상자별 국소 상수
  (H ≈ ‖H(y_c)‖ + M_B·(r+sweep))를 쓰면 평평한 곳에서는 거의 0이다. Δt에 선형으로 줄어든다.
- 계산 지연 τ_c도 같은 틀에서 [t_k, t_k+τ_c+Δt] 구간으로 확장 가능.

필요 작업: 행에 항 하나 추가 + 재실험(적어도 대표 장면). 사용자의 제약과 충돌 없음(안전집합은 그대로 연속,
보수성은 입력에만).

### 4.3 견고성 통합 remark (텍스트만, 짧게 가능)

모든 교란이 **같은 행의 affine tightening**으로 들어간다는 점을 한 단락으로 정리할 수 있다:
- 정적 불확실성(자세 추정 오차 δ_p, δ_R; 메시–실물 오차) → 레벨 상향 l_B + G(δ_p + r̄δ_R), 또는 ground truth 팽창.
- 속도 추종 오차 ‖d‖≤d̄ → 상수항 −G d̄ (정지 시 h ≥ G d̄/γ 필요 — ISSf 형태).
- 샘플링 → §4.2의 ‖u‖² 항. 부동소수점 → 최종 부등식의 outward-rounded 재평가(오프라인 Prop. 1은 ε=0.1 mm
  여유가 반올림 오차보다 수십 자릿수 커서 사후 구간 재평가만으로 엄밀화 가능).
리뷰어가 흔히 묻는 “perception/모델 오차에 대한 견고성”에 대한 답이 된다.

### 4.4 범위 확장 remark (선택)

가속도 입력(2차 동역학)으로의 확장: ψ_x = ḣ_x + αh_x에 같은 표면-족 논법 적용. ḧ_x = ẏᵀ∇²φ_B ẏ + ∇φ_Bᵀÿ에서
ÿ가 가속도 입력에 affine이고, 3차 스플라인의 Hessian이 Lipschitz이므로 Taylor 모델이 그대로 선다.
실험 없이 “구조가 유지된다”는 remark 수준으로만 권장.

---

## 5. 반영하지 않은 것과 이유

- native patch를 기본값으로 두는 설계, 실험 구성: 사용자 결정 사항(C-01) — 손대지 않음. 대신 §4.1이 이 설계를
  이론적으로 뒷받침하는 가장 좋은 경로다.
- §4.1, §4.2의 새 명제/정리: 방법·구현이 바뀌고 재실험이 필요하므로 초안만 제시.
- 그림, 안전집합 정의, SDF 사용: CLAUDE.md의 제약대로 불변.

---

## 6. 후속: 제안 §4.1, §4.2 구현 (2026-09-30, 사용자 요청 “다 고쳐”)

- 원고: Prop. `prop:second`(2차 정확도), 자유 multiplier 문단, Thm `thm:sampled`, Cor. `cor:sampled` 추가.
  Limitations의 두 문단을 “Closed loop and sampling” 하나로 통합. 인증서 비교 절은 고정 multiplier 4행을
  한 문장(λ≈1 해석)으로 압축.
- 코드: `summed.py`(mult_mode `free`, `sampled_tightening`, `solve(..., Wc)`), `summed3d.py`(`pair_rows(..., sd)`,
  `stack`), 새 파일 `sampled_data.py`, `certificate_upgrades.py`, `certificate_upgrades_report.py`.
  기본 경로(`vertex`/`one`)는 바뀌지 않았고, 기존 native 결과가 비트 단위로 재현됨을 확인.
- 수치 검증(증명 아님): 행의 실현가능 경계에서 7.2e5개 상자 점 모두 F+γw g_A > 0; sampled 곡률 상한은
  관측 |e''|의 8.5배 이상, 이동 상한은 관측 이동의 2.5배 이상.
- 첫 관찰: 역사적 native trial 18(1000 스텝 전부 slack)은 단위 lift의 상자 하한이 시작부터 음수(−1.66 mm)였기
  때문이며, 자유 multiplier로는 slack 0 — 2차 정확도 명제가 예측한 효과.
- 30회 결과: `results/certificate_upgrades/summary.md`.
- 최종 결과(`results/certificate_upgrades/README.md`): 12-mm unit 6/30(slack 2회) → free 8/30(slack 0) →
  sampled 7/30(u=0 대체 0). 6-mm sampled 13/28(대체 0, 최소 4.3 mm). 6-mm의 trial 16, 19는 계산 비용
  (dense NNLS QP, 활성 상자 1260–1490개) 때문에 사용자 결정으로 중단했고, 논문 표 IV 주석에 명시.
  두 trial은 어떤 방법도 도착하지 못한 trial이라 같은 28회 기준 비교(subdivided 17)는 불변.
