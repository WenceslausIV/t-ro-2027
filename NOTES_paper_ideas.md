# 논문 개선 아이디어 메모 (tro/main.tex)

논문을 읽으면서 떠오른 "더 좋은 방법"을 모아 둔다. 본문에 바로 반영한 것은 [반영]으로 표시한다. 나머지는 사용자 판단이 필요하거나 실험이 더 필요한 것이다.

## A. 방법 자체

1. **cover의 최악 step 시간 (p95 22 ms, 5대 swap)**
   - 원인 후보:
     - (i) 활성 모서리 수천 개 → 행 수가 큰 NNLS
     - (ii) Python에서 쌍마다 field 평가
   - 개선 방향:
     - 상자 모서리 대신 **상자 단위 barrier 하나** 쓰기: min over corners는 C¹이 아니다. 대신 Prop. taylor의 lower bound φ(c) − |∇φ(c)| r − ½κr² − l을 쓴다. 이것은 매끄럽고 corner 4/8개를 1개로 줄이지만, 더 보수적이다(+|∇φ|r).
       - 행 수가 1/2ⁿ로 줄고, 보수성은 r 정도(5 mm 상자에서 약 3.5 mm) 늘어난다.
       - 실험 한 번으로 표에 넣을 가치가 있다.
     - QP를 활성 집합 warm start로 풀기(이미 constraint generation 있음).
     - C 구현으로 옮기기.
2. **적응형 상자 크기**
   - 곡률이 큰 곳만 잘게, 평평한 곳은 크게 한다.
   - Prop. tight의 보수성 항 |∇φ|r + ½κr²를 목표 허용치 δ 이하로 맞추는 방식이다.
   - 상자 수가 줄고, 보수성은 균일하게 δ로 제어된다.
3. **공간적으로 변하는 level l(x)** (compiled field)
   - 본문에 이미 future work로 적혀 있다.
   - cover 쪽에는 이득이 작다(level이 이미 0.3–1.4 mm).
4. **이산시간 보장**
   - compiled field는 −2.4 mm까지 내려간다. cover는 ACT = 5 cm로 미리 켜져서 0 이상을 유지했다.
   - 5 cm 선활성화가 sampled-data 보장과 어떻게 연결되는지 한 문단으로 정리할 수 있다: 한 step당 h 감소량 ≤ Lipschitz × 이동량.
5. **데드락**
   - cover와 compiled 모두 같은 instance에서 데드락이 난다.
   - 비볼록 hook이 원인이라 필터만으로는 안 된다.
   - 별 튜브 역방향 실험도 같은 이유로 실패했다.
   - future work로 "planner reference + CBF"를 명시하면 좋다.

## B. 실험·표현

6. **Table III 시간 열**
   - 지금은 "instance 중앙값의 중앙값 / 최대"다.
   - Franka 표는 "전체 step의 median / p95"다.
   - 통일하면 좋다. 값은 `five_cover_stats_5.json`의 t_med_pooled, t_p95_pooled에 이미 있다(cover 3.6/7.0, compiled 0.63/1.18).
   - 단, closest와 circle은 pooled 값이 없어서 재실행이 필요하다.
7. **Docking 5 mm 대 2.5 mm 상자** 결과를 Prop. tight 예측선과 함께 작은 그래프로 보여주면 이론과 실험 연결이 선명해진다(페이지가 허용하면).
8. **Fig 1(b)** 확대 창의 위치를 (a)에 작은 사각형으로 표시하면 연결이 더 분명하다.
9. **두 팔 실험의 거리 평가**
   - 지금은 각 trial에서 "샘플 거리가 가장 작은 20개 step"에서만 정확한 점-삼각형 거리로 하한을 구한다.
   - 샘플로 step을 고르는 셈이라 "샘플 기반 평가 금지" 원칙에 약하게 걸린다.
   - 개선: 모든 step에서 하한을 계산하거나, 평면 실험처럼 step 사이 구간 감사(rigid-motion bound + bisection)를 적용한다. 계산량은 크지만 오프라인이라 가능하다.
10. **Docking 문단의 "Both C⁰ fields" (15.8, 9.6 mm)**
    - 두 번째 C⁰ field(30³ cells로 추정)는 표에 열이 없어서 독자가 헷갈린다.
    - 표 각주에 한 줄로 설명하거나, 문장을 "the $C^0$ field (and a coarser one)"처럼 풀어 쓰면 좋다.

## C. 이번에 다듬으며 반영한 것 [반영]

- Abstract를 문장 단위로 나눴다(offline cover / online corner barrier / compiled 변형 / 결과).
- Intro "this minimum" → "the minimum of that SDF over that boundary"로 대상을 명확히 했다.
- Certified Regions: 주 방법이 쓰는 SDF level 문단을 앞으로, compiled field만 쓰는 spline 경계 문단을 뒤로 옮겼다.
  - B_i의 두 의미(cover: {φ ≤ l}, compiled: spline 내부)도 명시했다.
- cover 상자 중심 기호 c̄_k → c_k. cluster 중심 c̄와 겹치지 않게 했다.
- Planar Setup: 한 문장이던 설정을 제어 / cover / compiled / 계산 환경으로 나눴다.
- 5대 swap: cover 시간(median 1.8 ms, p95 22 ms)을 숨기지 않고 본문에 적었다.
- 미로: 로켓의 64 제어점·4 mm offset은 compiled 전용이라 compiled 문장으로 옮겼다.
- "Our field"를 "Our compiled field"로 바꿨다(sampling/ball-world 비교).
