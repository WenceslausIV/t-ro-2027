# T-RO 수정 원고: 독립 리뷰어 3명 검토 종합

검토일: 2026-09-26. 이는 실제 저널 심사가 아니라 이론, 실험, 기여도에 대한 독립 모의 리뷰다. 원고·실험 코드는 수정하지 않았다.

검토 원고: `tro/main.tex`, 1189줄, SHA-256 `4bbdd6f3c8f02a82ad8fbe4a30422a15d22a195d0b9ec7dad199b20a4b166e27`.
검토 중 별도 프로세스가 실험 코드와 로그를 갱신하고 있었다. 아래 수치와 줄 번호는 확인한 파일 상태 기준이며, 이후 수정본까지 잘못되었다는 뜻은 아니다.

## 판단

**이전 원고보다 이론의 가정과 보장 범위가 훨씬 정확해졌다. T-RO급 연구로 발전시킬 중심 기여는 남아 있지만, 현재 원고 그대로의 투고는 권하지 않는다.**

핵심 최초 접촉 논리와 상대 자세 미분식에서 새로 확인한 치명적 반례는 없다. 가장 확실한 투고 차단 요인은 도킹 실험의 코드·결과·본문·그림이 서로 다른 버전을 설명한다는 점이다. 그다음은 구현 제어기가 이론의 정칙성 가정을 만족하는지, 수치 검증이 어떤 수준까지 보장하는지, 강한 비교법에 대해서도 이점이 남는지다.

각 리뷰의 상세 근거:

- [Reviewer 1 — 이론·증명](review_round2_reviewer1.md)
- [Reviewer 2 — 실험·구현·재현성](review_round2_reviewer2.md)
- [Reviewer 3 — 기여도·관련 연구·저널 적합성](review_round2_reviewer3.md)

| 관점 | 모의 판정 | 주된 이유 |
|---|---|---|
| Reviewer 1: 이론 | Major revision | 핵심 조건부 증명은 타당하나 실제 switching feedback에 대한 적용 조건 보강 필요 |
| Reviewer 2: 실험 | Major revision | 도킹 버전 불일치, intersample 검증의 미결정 처리, 측정 범위 문제 |
| Reviewer 3: 기여도 | Weak reject / 큰 수정 | 실질적인 확장은 인정하지만 강한 비교법과 신규성 설명 부족 |

## 1. [높음, 확인된 불일치] 도킹 대표 실험을 하나의 버전으로 통일해야 한다

근거: `tro/main.tex:905–973`, `tro/main.tex:1139–1142`, `dock_shapes.py:14–18`, `results/dock2.json`, `tro/figs/dock_paper.png`.

| 항목 | 원고·그림 | 현재 코드·저장 결과 |
|---|---:|---:|
| full seating 깊이 | 0.32 m | 0.58 m |
| 설계 clearance 파라미터 | 0.03 m | 0.025 m |
| B-spline cells | 75 × 75 × 120 | 78 × 78 × 144 |
| B-spline 계수 수 | 약 730,000 | 944,784 |
| Bézier 계수 수 | B-spline과 같다고 설명 | 745,290 |
| B-spline seating shortfall | 5.0 cm | 7.073 cm |
| B-spline final ground-truth gap | 42.0 mm | 37.218 mm |
| Bézier final ground-truth gap | +8.8 mm | −2.150 mm |
| Bézier minimum ground-truth gap | 본문은 접촉하지 않았다고 설명 | −2.631 mm |
| Bézier minimum h | −30.4 mm | −36.050 mm |

본문 969줄은 두 비교법의 실제 형상이 닿지 않았다고 하는데, 한계 절 1142줄은 Bézier가 충돌했다고 한다. 현재 JSON은 후자와 일치한다. 그림도 0.32 m seating 및 약 −30 mm Bézier h를 표시하므로 단순 표 오타가 아니다.

**수정:** 대표 실험의 geometry/config/result를 고정한 다음 그 결과에서 표·본문·그림을 함께 생성한다. 현재 설정으로 비교하면 B-spline은 Bézier보다 약 26.8% 많은 계수를 쓰므로 equal-size 비교라는 문구도 바꿔야 한다. 이전 설정을 유지하려면 그 설정의 결과와 생성 근거를 별도로 보존해야 한다.

## 2. [중간, 구현 적용 조건] hard activation과 locally Lipschitz feedback 사이가 연결되지 않았다

근거: `tro/main.tex:411–442`, `tro/main.tex:528–537`.

Corollary는 QP 최적해가 locally Lipschitz라고 명시적으로 가정한다. 그러나 `||t||∞ < D−δ` 경계에서 제약을 갑자기 추가하는 구현에서는, band에서 h>0라는 사실만으로 제어입력의 연속성이 따라오지 않는다. h>0여도 접근 속도가 크면 CBF 제약이 제어입력을 바꿀 수 있다.

예를 들어 φ(t)=||t||², l*=r²이고 activation 경계가 t=(a,0), a>r이면, 경계 안쪽에서는 inward velocity가 `v_x >= −γ(a²−r²)/(2a)`로 제한되지만 바깥쪽에서는 nominal input이 그대로 허용된다. 두 값이 다르면 jump가 생긴다. 접촉 반경을 r보다 작게 잡으면 contact exclusion 및 level regularity와도 양립한다.

이는 **명시된 가정 아래의 정리가 틀렸다는 뜻이 아니다.** 구현이 그 가정을 만족한다는 설명이 부족하다는 뜻이다.

**수정 선택지:** activation collar에서 제약이 모든 허용 입력에 대해 redundant임을 추가 검증하거나, 실제 switched filter에 대해 absolutely continuous trajectory와 a.e. differential inequality를 사용한 안전성 논리를 제시한다. 단순히 h>0만 확인하는 것과 구분해야 한다.

## 3. [중간, 재현된 검증 한계] intersample 검사에서 '미결정'과 '안전'을 구분해야 한다

근거: `cspace_sdf_cbf_compare.py:574–610`의 `step_gap`, 특히 604줄; `tro/main.tex:700–708`.

이 함수는 endpoint distance가 이동 상계보다 큰 구간을 제외하는 올바른 논리를 포함한다. 그러나 이동량이 0.1 mm보다 작거나 최대 깊이에 도달하면, 아직 안전을 증명하지 못한 구간도 종료하고 양의 sampled distance를 반환한다.

리뷰어 2의 작은 반례: 한 변 20 μm인 두 정사각형 중 하나를 중심 +40 μm에서 −40 μm로 이동시키면 endpoint gap은 양수지만 중간에는 겹친다. 기본 tolerance에서는 양의 gap을 반환할 수 있다. 이는 검사 함수의 반례이며 **원고의 실제 B-spline 실험에서 충돌이 있었다는 증거는 아니다.**

또한 이 함수의 반환값은 연속 시간의 정확한 minimum distance가 아니다. 멀리 떨어진 pair에서는 bounding-radius 하계를, 근접 구간에서는 평가한 시점의 최소치를 사용한다.

**수정:** `collision / verified_free / unresolved`와 미결정 구간 수를 별도로 기록한다. tolerance 종료를 안전 판정으로 합치지 않는다. 거리 표는 sampled minimum인지, certified lower bound인지 명시한다. 유니사이클은 선형 보간 경로와 zero-order-hold 실제 원호 궤적도 구분한다.

긍정적 확인: 저장된 `five_random_trajectory.npz`의 1500개 구간은 실제 polygon 반경과 pose increment를 이용한 별도 이동 상계 검사에서 약 **14.0 mm의 양의 여유**가 남았다. 따라서 해당 선형 보간 궤적의 안전성을 의심할 근거는 발견하지 못했다. 이 값은 원고의 sampled minimum 18.8 mm와 다른 종류의 지표다.

## 4. [중간, 이론보다 강한 표현] 초기 조건과 certified level의 의미

근거: `tro/main.tex:163–164`, `tro/main.tex:529–532`, `tro/main.tex:752–783`.

- “from any collision-free start”는 너무 넓다. 안전하게 분리되어 있어도 fitted barrier가 h<0일 수 있다. Corollary처럼 fitted bodies의 초기 비충돌 및 domain 안에서는 h>=0 조건을 함께 명시해야 한다.
- l*는 `max_M φ`의 ε-정확 상계다. M에는 C-obstacle 내부의 boundary-crossing configuration도 들어간다. 따라서 765줄의 “the largest amount … at the obstacle boundary”는 등식처럼 읽히며 정확하지 않다. **장애물 경계에서의 overestimate에 대한 상계**라고 해야 한다.
- 서로 다른 field의 l* 크기만으로 보수성을 순위 매길 수 없다. φ와 l*에 같은 상수를 더하면 h는 그대로인데 l*만 달라진다. 767줄에서는 이미 인정하지만 755–759줄의 “less conservative / no certificate … less conservative”에는 같은 주의가 적용되지 않았다. 실제 safe-set geometry, clearance, attainable insertion depth로 비교하는 편이 정확하다.
- 622줄의 “rounding errors … small compared with 2 mm”는 rounding bound가 없다는 뒤 문장과 맞지 않는다. bound나 검증 실험이 없다면 크기가 확인된 사실처럼 쓰지 않는다.
- 126–129줄 및 322–324줄의 closest-point smoothness 설명도 정교화할 수 있다. Strict convexity만으로는 곡률 비퇴화가 보장되지 않고, 반대로 비볼록 body라도 unique nondegenerate closest pair가 있는 국소 구간에서는 smoothness가 가능하다. Reviewer 1은 smooth strictly convex superellipse를 이용한 반례를 제시했다. 이는 제안한 핵심 알고리즘의 오류가 아니라 기존 방법의 가정 설명 문제다.

## 5. [중간, 관련 연구·기여 표현] '제약 하나'와 'certified C-space' 자체를 독점적 차별점으로 삼지 않는다

근거: `tro/main.tex:195–241`, Reviewer 3의 원문 대조.

기존 Boolean CBF composition도 여러 primitive 조건을 하나의 smooth scalar barrier로 합칠 수 있다. 따라서 convex decomposition이 반드시 여러 QP row를 요구하는 것처럼 일반화하면 반박받을 수 있다. 제약 개수와 온라인 평가에 필요한 형상 계산량을 구분해야 한다.

- [Molnar & Ames, Composing Control Barrier Functions for Complex Safety Specifications](https://arxiv.org/abs/2309.06647): smooth composition과의 차이를 설명한다.
- [Dai et al., Certified polyhedral decompositions of collision-free configuration space](https://doi.org/10.1177/02783649231201437): C-IRIS 계열의 기존 C-space certificate와의 관계를 논의한다. 직접 baseline을 반드시 실행해야 한다는 뜻은 아니다.
- [Wu et al., Optimization-free Smooth Control Barrier Function for Polygonal Collision Avoidance](https://arxiv.org/abs/2502.16293): 기본 polygon 가정의 convexity와 composite/non-convex 확장을 구분해서 서술한다.

방어 가능한 중심 기여는 **회전하는 비볼록 planar body pair의 연속 접촉 집합에 대한 검증으로, 부정확한 학습 target과 독립적으로 smooth field의 안전 수준을 정하고, 온라인에서는 local-support field를 평가하는 구조**다.

ICRA 2027 원고와 공유하는 winding/enclosure/B-spline 변환을 별도 기여로 주장하지 않는 현재 설명은 적절하다. 공통 도구를 쓴다는 사실만으로 중복 투고라고 판단할 근거는 없다. 관계 설명과 관련 원고 제공은 유지한다.

## 6. [높은 심사 위험, 수학 오류 아님] 강한 baseline 및 정량적인 task success가 필요하다

근거: `tro/main.tex:329–338`, `tro/main.tex:962–973`, `tro/main.tex:1032–1053`, `tro/main.tex:1163–1167`.

현재는 단일 closest feature와 bounding circle만 비교한다. 한계를 솔직하게 적은 것은 개선이지만, reviewer가 요구할 비교 실험 자체를 대신하지는 못한다. multi-feature/nonsmooth 또는 convex-decomposition 기반의 타당한 비교법 하나가 가장 설득력 있는 다음 실험이다.

Random swaps에서는 안전한 완료 수가 ours 13/20, circle 20/20이다. 논문도 이를 인정한다. 그러므로 무조건적인 전체 성능 우위보다 **좁은 형상 상호작용에서의 유용성**을 입증해야 한다. 대표 docking도 full seating 전에 정지하므로, success tolerance를 사전에 정하고 성공률·잔여 오차를 함께 평가해야 한다. 목표가 의도적으로 infeasible인 lobe insertion은 safety/contact-stop 실험으로 해석한다.

검토 중 `results/run_dock2dt.log`에는 10 ms와 5 ms 실험이 진행되고 있었다. timestep sweep 코드가 없다는 지적은 하지 않는다. 완료 후 결과와 원고에 반영하고, 이것만으로 feature-switching 문제를 분리했다고 단정하지 않는다.

Hardware는 T-RO의 무조건적 필수 요건이라고 판단하지 않았다. 현재 가장 먼저 할 일은 이미 있는 증거의 일관성 확보와 적절한 비교법이다.

## 7. [중간, 측정 정의] certification timing의 포함 범위

근거: `tro/main.tex:916–917`, `cspace_experiments.py`의 `dock2_setup`, `results/dock2.json`.

원고는 band 및 regularity를 포함한 인증 시간이라고 설명한다. 저장되는 `t_cert`는 `cert.certify()`가 반환하는 시간이며 band/regularity 호출은 그 이후다. 정확한 숫자의 크기와 별개로 포함 범위가 다르다.

**수정:** contact / band / regularity / total 시간을 분리하거나, 전체를 둘러싼 wall-clock timer를 사용한다. 형상 fitting 및 target 생성 시간의 포함 여부도 같은 방식으로 명시한다. 런타임 수치는 재실행 환경에 따라 바뀔 수 있으므로 configuration/run 식별자와 묶는다.

## 수정된 것으로 인정하는 항목

- Ground-truth connectedness 및 복수 boundary/component 설명.
- First-contact exclusion과 full C-obstacle enclosure의 구분.
- Driftless safe-set에서 u=0에 의한 공동 QP feasibility.
- 잘못된 'subdivision마다 image width 절반' 대신 개선된 종료 논리.
- Exact-arithmetic theorem과 검증되지 않은 floating-point 구현의 구분.
- C0 cellwise gradient 검사와 C1 regularity의 구분; global polynomial periodic seam 한계.
- Random trials min h −2.4 mm, goal reaching과 collision-free completion 구분.
- CPU/GPU baseline timing 경로 공개.
- Closest-point 비교법의 제한된 범위 및 deadlock 한계 공개.
- 익명 모드, bibliography placeholder 정리, 새 관련 연구 추가.

## 정적 검사와 검토 한계

- Abstract는 공백 기준 190단어.
- 미정의 citation/reference, 중복 label, 누락된 includegraphics 파일은 발견하지 못했다.
- 도킹 및 5-robot 그림을 직접 확인했다. 도킹 그림은 구버전 수치다.
- `pdflatex`, `latexmk`, `tectonic`이 PATH에 없어 PDF 컴파일, 최종 페이지 수, float 배치 검사는 하지 못했다.
- 전체 실험을 다시 실행하지 않았다. 저장 JSON/log 대조, 코드 분석, 작은 수치 반례, 저장 5-robot 궤적의 이동 상계 검사를 사용했다.
- 본문과 관련 코드가 동시에 수정될 수 있으므로, 투고 전에는 하나의 고정된 원고·설정·결과 버전으로 최종 일치 검사를 해야 한다.

## 권장 수정 순서

1. 도킹 geometry/config/result/표/그림의 버전을 통일한다.
2. 초기 조건, level 의미, activation feedback 가정을 정리한다.
3. collision audit의 unresolved 상태와 정확한 거리·시간 측정 정의를 저장한다.
4. 적절한 강한 baseline 및 완료 기준을 갖춘 대표 실험을 보강한다.
5. 기존 certified C-space 및 Boolean composition과의 차이를 분명히 쓴다.
6. PDF를 컴파일해 투고 형식과 그림 가독성을 확인한다.
