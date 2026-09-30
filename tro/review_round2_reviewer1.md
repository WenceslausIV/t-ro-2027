# Reviewer 1 — 기하·증명·CBF 안전성 재검토

검토일: 2026-09-26. 현재 `tro/main.tex` 전체와 접촉 인증, regularity, activation, collision 검사 구현을 읽었다. 원고·시뮬레이션 코드는 수정하지 않았다. 아래 줄 번호는 이번 검토 시점 기준이다.

## 총평

**수학적 핵심은 이전 버전보다 상당히 명료해졌고, 현재 first-contact lemma 및 그에 기반한 조건부 collision-avoidance theorem에서 반례를 찾지 못했다. 다만 실제 switching QP가 corollary의 feedback 가정을 만족하는지는 아직 입증되지 않았다. 이 부분을 보강하고 앞부분의 보장 범위를 정리하는 major revision을 권한다.** 이는 핵심 접촉 인증 알고리즘이 틀렸다는 판정이 아니다.

좋아진 부분은 다음과 같다.

- `main.tex:470`의 first-contact 증명은 compactness와 연속 경로만으로 양쪽 접점이 경계에 있음을 올바르게 보인다. 비볼록성이나 경계 미분 가능성을 요구하지 않는다.
- `main.tex:372`의 상대 좌표 시간미분은 회전 프레임 항과 두 로봇의 angular coefficient 부호가 맞다.
- `main.tex:586`의 strict pruning 및 nondecreasing lower bound를 통한 인증 논리는 타당하다. `main.tex:598`은 이전의 잘못된 “image width가 정확히 반감” 주장을 고쳤다.
- `main.tex:458`의 connected ground truth와 winding argument, `main.tex:418`의 중앙집중 driftless feasibility, `main.tex:626`의 level-set regularity 검사는 명시된 전제에서 타당하다.
- contact exclusion과 obstacle 전체 sublevel enclosure, 연속시간과 sampled controller, exact arithmetic과 double precision을 구분한 것은 적절하다.

## R1-1. 활성화 경계에서 QP feedback의 locally Lipschitz 가정이 자동으로 성립하지 않는다

**중요도: P2 / 이론과 구현 사이의 조건부 공백. 정리 자체의 오류는 아님.**

근거: `tro/main.tex:411`, `tro/main.tex:434`, `tro/main.tex:528`, `cspace_sdf_cbf_compare.py:662`의 활성 pair 처리.

Corollary는 QP minimizer가 locally Lipschitz라고 명시적으로 가정한다. 따라서 조건부 정리는 이 점을 숨기지 않는다. 그러나 구현은 `|t|_inf < D-delta`에서 constraint를 갑자기 넣고 뺀다. activation band에서 `h>0`인 사실만으로 해당 constraint가 redundant라는 결론은 나오지 않는다. 따라서 band 및 regularity certificate를 모두 통과해도 실제 feedback은 활성화 면에서 불연속일 수 있다.

간단한 허용 가능한 예:

- 정지한 원판과 이동하는 원판의 반지름 합을 `R=2.7`로 둔다.
- `D=R+0.3=3`, `delta=0.125`, 활성화 위치 `a=2.875`.
- periodic theta-independent quadratic field `phi(t)=||t||^2/(2R)`, certified level `l*=R/2+0.002=1.352`. 접촉 집합에서 `max phi=R/2=1.35`이므로 strict contact certificate이며, 2 mm tolerance와도 부합한다. 이 quadratic은 cubic spline 공간에 표현 가능하다.
- 활성화 band에서 `h>0`이고 level set에서 gradient도 nonzero다.
- `gamma=5`, 입력 상한 `1`, `t=(a,0)` 근처에서 nominal radial velocity가 `-1`이라고 하자.

활성화 면의 바깥에서는 constraint가 없으므로 `u_x=-1`. 안쪽에서는 `2a/(2R) u_x >= -gamma h(a)`에 의해 `u_x >= -0.839` 정도로 제한된다. 따라서 면을 건너면서 feedback이 점프한다. 모든 certificate를 통과했다는 사실만으로 corollary의 globally stated locally-Lipschitz hypothesis를 결론낼 수 없다.

**수정 제안:** (a) switching collar 전체에서 해당 constraint가 모든 허용 입력에 대해 redundant임을 추가로 인증하거나, (b) 적절한 연속적인 활성화/전역 barrier extension을 구성하거나, (c) absolutely continuous trajectory에서 active interval별 differential inequality와 band separation을 사용하는 안전성 정리를 제시한다. (c)를 택하면 solution existence와 사용할 solution notion을 명시해야 한다. 실제 실험에서 이 불연속이 발생했다는 주장은 하지 않는다. 현재 인증 항목만으로 발생하지 않음을 보이지 못한다는 지적이다.

## R1-2. “any collision-free start”는 실제 보장보다 강하다

**중요도: P2 / 명확한 보장 범위 불일치.**

근거: `tro/main.tex:163` 대 `tro/main.tex:516`, `tro/main.tex:529`; 관련 요약 `tro/main.tex:504`, `tro/main.tex:1129`.

기여 항목은 “from any collision-free start”라고 쓰지만 corollary는 초기 fitted-body separation에 더해, field domain 안에서는 `h(0)>=0`을 요구한다. strict contact exclusion과 field continuity 때문에 실제 contact 부근의 collision-free 상태 중에도 `h<0`인 상태가 일반적으로 존재한다. outward boundary enclosure까지 고려하면 ground truth끼리 떨어져 있어도 fitted bodies가 겹칠 수 있다.

따라서 physical collision-free, fitted-body collision-free, certified superlevel membership은 서로 다른 초기 조건이다. 현재 정리는 이들 중 필요한 조건을 명확히 요구하므로 정리를 바꿀 필요는 없고, 기여·요약 문장을 “from initially separated fitted bodies with nonnegative active barriers” 정도로 맞추면 된다. 일반 collision-free 초기 상태를 모두 수용하려면 별도의 initialization/recovery 결과가 필요하다.

## R1-3. certified level은 boundary overestimate의 정확한 최대값이 아니다

**중요도: P2 / 수학적 해석 과장.**

근거: `tro/main.tex:763`–`767`, `tro/main.tex:781`–`783`, `tro/main.tex:744`–`748`.

`l*`가 “the largest amount by which the field overestimates the distance at the obstacle boundary”라고 단정하는 것은 정확하지 않다. 증명된 관계는

`max_boundary(C) phi <= max_M phi < l* <= max_M phi + epsilon`

이다. 첫 번째 부등식이 equality일 이유가 없고, 마지막에는 numerical branch-and-bound tolerance가 들어간다. 원고가 올바르게 설명했듯 `M`에는 obstacle 내부에서 양쪽 경계가 교차하는 configuration도 포함된다. arbitrary fitted smooth field는 그런 내부 configuration에서 최대를 가질 수 있다.

예를 들어 같은 크기의 원판 쌍에서는 한 orientation의 `M`이 translation collision disk 전체를 채운다. obstacle boundary에서는 정확하지만 내부에서 양의 smooth bump를 갖는 field를 생각하면 `max_M phi`는 boundary error가 아니라 내부 bump가 결정한다. 따라서 “fit accuracy at the obstacle boundary가 certificate를 결정한다”는 일반 정리는 나오지 않는다.

**수정 제안:** “a certified upper bound on boundary overestimation, obtained over the larger contact set and including at most epsilon optimization slack”로 바꾼다. 실험 field의 최대점이 실제로 외곽 collision boundary 근처인지 별도로 확인했다면 그 결과에 한정해서 empirical interpretation을 유지할 수 있다. “a smaller level does not imply larger safe set”이라는 현재 caveat는 좋은 수정이며 그대로 유지하면 된다.

## R1-4. strict convexity만으로 closest-point trajectory의 smoothness를 보장한다는 문장은 과하다

**중요도: P3 / 기존 방법의 한계 설명 정교화.**

근거: `tro/main.tex:126`–`129`, `tro/main.tex:322`–`324`.

“smoothness … guaranteed only for convex level sets”는 convexity가 필요조건인 것처럼 읽히고, “This holds for strictly convex level sets”는 strict convexity만으로 충분하다는 인상을 준다. 비볼록 body라도 특정 configuration의 unique nondegenerate closest pair는 국소적으로 smooth하게 움직일 수 있다. 반대로 strict convexity만으로 curvature nondegeneracy가 보장되지는 않는다.

구체적으로 smooth strictly convex body `x^4+y^4<=1` 두 개의 중심을 `(0,0)`, `(3,0)`에 두고 두 번째 body를 작은 theta만큼 회전시키자. 접촉 후보를 `a=(f(y_a),y_a)`, `b=(3,0)+R(theta)(-f(y_b),y_b)`, `f(y)=(1-y^4)^(1/4)`로 표현하면 최소거리의 stationarity는 `(b-a) dot a'=0`, `(b-a) dot b'=0`이다. theta=0 부근에서 이 식은 `y_a^3+y_b^3~theta`, `y_b-y_a=O(theta)`를 주므로 `y_a~(theta/2)^(1/3)`. closest-point 위치는 theta에 대해 미분 가능하지 않다. 두 곡선은 해당 지점에서 곡률이 0이지만 여전히 smooth하고 strictly convex다.

bounded numerical stationarity check도 이 scaling을 확인했다: theta `1e-3, 1e-5, 1e-7`에서 `y_a/theta`가 각각 약 `78.6, 1709, 36840`으로 증가하고 residual은 `1e-17` 이하였다. 이 확인은 논문의 핵심 알고리즘에 대한 수치 실험이 아니라 설명 문장의 반례 점검이다.

**수정 제안:** global uniqueness/regularity가 확보된 경우와 nonconvex feature switching을 구분하고, closest-point smoothness에는 적절한 nondegeneracy 조건을 명시한다. 이 지적은 제안한 C-space field approach의 장점을 부정하지 않는다.

## 다른 reviewer와의 분담 및 남는 한계

`step_gap`가 residual-motion tolerance나 max-depth에 도달한 unresolved interval을 별도로 표시하지 않는 문제는 Reviewer 2가 실제 반례와 함께 다룬다. 여기서는 중복 항목으로 세지 않았다. finite-precision certificate와 sampled-data guarantee 부재는 이미 원고가 인정한 한계이며, 이번 검토에서 새 오류인 것처럼 다시 주장하지 않았다.

이론 쪽의 가장 중요한 다음 작업은 **수학적으로 올바른 contact certificate와, 실제로 쓰는 switching sampled QP가 정확히 어떤 조건에서 연결되는지**를 한 문장 이상의 결과로 완결하는 것이다. 현재 수준에서 “핵심 theorem이 틀렸다”는 판단은 부당하지만, “certificate가 통과했으므로 구현된 controller의 조건도 모두 충족된다”는 판단도 아직 이르다.
