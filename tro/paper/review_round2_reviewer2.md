# Reviewer 2 — 실험·구현·재현성 재검토

검토일: 2026-09-26. 현재 원고와 코드·JSON·로그를 읽고, 작은 충돌 검사 반례와 저장된 five-robot trajectory의 motion bound를 실행했다. 전체 시뮬레이션은 재실행하지 않았으며 원고·코드·결과는 수정하지 않았다. 검토 중 다른 프로세스가 `cspace_sdf_cbf_compare.py`와 timestep sweep 로그를 갱신했으므로 아래 스냅샷 기준이다.

**판정: 실험 부분은 major revision. 가장 급한 문제는 docking 실험의 새 코드/결과와 이전 버전 표·설명이 섞인 것이다.** 연구 아이디어를 부정하는 오류는 아니지만, 현재 PDF에 상충하는 충돌 판정이 들어가 있어 제출 전 반드시 정리해야 한다.

## 1. [P1, 확인된 불일치] Docking의 형상·모델·성과가 서로 다른 버전이다

근거: `tro/main.tex:906`, `tro/main.tex:916`, `tro/main.tex:924`, `tro/main.tex:945`, `tro/main.tex:968`, `tro/main.tex:1141`; `dock_shapes.py:14`; `results/dock2.json:7`.

| 항목 | 원고 | 현재 코드/저장 결과 |
|---|---:|---:|
| full seating 깊이 | 0.32 m | 0.58 m (`0.88-0.30`) |
| 설계 좌우 clearance | 3 cm | 2.5 cm |
| B-spline cells | 75×75×120 | 78×78×144 |
| B-spline coefficients | 약 730,000 | 944,784 |
| Bézier coefficients | B-spline과 동일하다고 설명 | 745,290 |
| certified level, B-spline / Bézier | 12.3 / 28.4 mm | 12.491 / 27.081 mm |
| B-spline seating shortfall | 5.0 cm | 7.073 cm |
| B-spline final gap | 42.0 mm | 37.218 mm |
| Bézier final gap | +8.8 mm | −2.150 mm |
| Bézier minimum gap | 본문: true shapes do not touch | −2.631 mm |
| Bézier minimum h | −30.4 mm | −36.050 mm |
| closest minimum h | −13.6 mm | −16.385 mm |

`main.tex:969`는 두 대안의 실제 형상이 닿지 않았다고 하지만 `main.tex:1142`는 Bézier의 실제 충돌을 인정한다. 현재 `dock2.json:55`와 `:56`은 후자와 일치한다. 이는 표현 취향이 아니라 핵심 안전성 결과의 직접적인 모순이다. 최종 gap과 trajectory 전체 minimum gap도 분리해서 제시해야 한다.

같은 coefficient budget이라는 주장도 현재 docking에는 성립하지 않는다. B-spline이 Bézier보다 **26.77% 많은 계수**를 사용한다. 기본 C-shape representation 실험의 110,592 대 115,248 비교는 여전히 별개로 유효하다.

수정: docking 형상과 결과 버전을 하나로 고정하고 본문·표·그림·caption·limitations를 일괄 재생성한다. 같은 budget 비교를 유지하려면 실제 모델 크기를 맞추거나, 서로 다른 budget이라고 명확히 쓰고 별도 matched-budget ablation을 제시한다. angular-resolution sweep의 72/96/120-cell 수치도 바뀐 형상으로 얻은 결과인지 다시 확인해야 한다.

## 2. [P2, 재현 확인] 새 intersample 검사는 개선됐지만, unresolved interval을 안전 판정과 구별하지 않는다

근거: `cspace_sdf_cbf_compare.py:574`의 `step_gap`, 특히 `:604`; `tro/main.tex:700`.

`max(endpoint distances)>motion bound`일 때 interval이 안전하다는 판정은 타당하다. 그러나 `w<tol` 또는 `max_depth`에 도달했을 때도 동일하게 종료하여, 아직 충돌을 배제하지 못한 interval을 양수 값으로 반환할 수 있다. 본문은 0.1 mm tolerance를 밝히지만, 그 잔여 구간을 미확정 상태로 분류하지 않는 점까지는 설명하지 않는다.

현재 함수로 실행한 반례:

```python
P = np.array([[-1e-5,-1e-5], [1e-5,-1e-5],
              [1e-5, 1e-5], [-1e-5,1e-5]])
x0 = np.array([[0.,0.,0.], [4e-5,0.,0.]])
x1 = np.array([[0.,0.,0.], [-4e-5,0.,0.]])
rho = np.linalg.norm(P,axis=1).max() * np.ones(2)
C.step_gap([P,P],rho,x0,x1)           # +0.000020 m
C.true_distance(P,P)                # 0.0 at the midpoint: overlap/contact
C.step_gap([P,P],rho,x0,x1,tol=1e-6)  # 0.0
```

이는 작은 형상의 단위 반례이며 **논문의 저장된 trajectory에서 누락된 충돌을 발견했다는 뜻은 아니다.** tolerance 아래의 grazing contact에도 같은 논리적 한계가 있다. 함수 이름/docstring의 “smallest distance over one integration step”도 정확한 연속시간 최소값이라는 뜻으로는 성립하지 않는다. 안전하다고 판정한 구간은 중간 거리의 최솟값을 더 찾지 않으며, 먼 pair는 bounding-circle lower bound를 반환한다.

수정: `safe / collision / unresolved`를 구분하고, 모든 실험의 unresolved 수를 저장한다. 최소 거리 측정값과 interval별 검증된 하한도 구별한다. 무충돌을 주장하려면 unresolved=0을 보이거나 남은 구간에 보수적인 처리를 적용한다.

## 3. [P2, 확인된 계측 범위 불일치] Certification time이 band/regularity를 포함한다는 설명은 코드와 다르다

근거: `tro/main.tex:916`; `cspace_experiments.py:887`–`:891`. `t_cert`는 `cert.certify()`에서 받은 시간이며, `certify_band()`와 `certify_regular()`는 그 뒤에 실행한다. 그런데 본문은 “certification including band and regularity in 14.6 s”라고 쓴다. 현재 JSON의 14.767 s도 같은 필드이다.

`repr`의 `cspace_experiments.py:182`와 static의 `:631`도 같은 구조다. 해당 표의 certification time은 contact-level B&B 시간이라고 구체화하거나, 별도 wrapper timer로 세 검증 전체를 재측정해야 한다. Field fitting과 raster generation 포함 여부도 표마다 구별해야 한다.

## 4. [P2, 평가 범위] Timestep과 비매끄러운 baseline의 차이를 아직 원고 결과로 충분히 분리하지 않았다

현재 원고는 단일 closest feature의 한계와 10 ms discretization을 명시하고 있어 이전보다 공정하다 (`tro/main.tex:1033`). 그러나 C²라는 성질 자체가 충돌/채터링 개선의 원인이라는 강한 결론에는 approximation accuracy, grid resolution, certificate level, discrete integration이 함께 바뀌는 문제가 남는다. 특히 새 docking의 unequal coefficient budgets는 이 해석을 더 어렵게 한다.

**dt sweep은 코드에 이미 구현되어 있고 실행 중이므로 ‘없는 실험’이라고 판정하지 않는다.** 검토 시 `results/run_dock2dt.log`에는 10 ms와 5 ms 결과가 있으며, Bézier min gap이 각각 −2.631 / −2.466 mm, min h가 −36.050 / −35.614 mm였다. 5 ms에서도 문제가 남는다는 유용한 증거다. 최종 2/1 ms 결과까지 확정한 뒤 이를 원고에 포함하면 시간 이산화와 representation의 영향을 더 잘 구분할 수 있다.

남아 있는 실험적 약점은 저자도 인정한 all-active-feature nonsmooth 또는 decomposition-based baseline의 부재이다 (`tro/main.tex:1164`). 단순히 limitation 문장을 추가했다고 강한 journal 비교가 완성되는 것은 아니다. 최소 하나의 강한 geometric baseline을 같은 input bounds/step size에서 비교하고, activation 기준·CPU/GPU 비용·안전 여유·성공 조건을 명시하는 것이 좋다. 이는 발견된 수학적 오류가 아니라 T-RO 수준 설득력의 문제다.

## 5. [P3, 해석 제한] Static obstacle 수 증가 실험에서 활성 장애물 수는 최대 0/0/1/2이다

근거: `results/static.json`의 각 method `n_constraints`, `cspace_experiments.py:643`, `tro/main.tex:1113`.

1/3 obstacles에서는 closest baseline constraint가 전혀 활성화되지 않고, 5에서는 최대 1개, 7에서는 최대 2개다. 따라서 0.2→2.8 ms 차이는 총 장애물 수 scaling뿐 아니라 no-interaction→interaction의 변화도 포함한다. 하나의 합집합 field lookup이 장애물 수에 직접 의존하지 않는 구조적 설명은 타당하지만, 이 표만으로 dense-many-obstacle 상황의 실증적 scaling까지 보여준 것은 아니다.

수정: 각 경우 활성 constraint 평균/최댓값을 표기하고, 동일한 local interaction 조건에서 obstacle 수만 변화시키는 실험을 추가하거나 결론을 현재 네 배치로 제한한다.

## 실제로 고쳐진 부분과 독립 확인

- 이전 edge-crossing 반례인 수직으로 겹친 긴 직사각형에서 현재 `true_distance`와 `closest_pair` 모두 0.0을 반환했다. 이전의 +1.8 오류는 해결됐다.
- Random stats 표는 현재 JSON과 일치한다. B-spline 최저 h=−2.386566 mm, 충돌 0/20, 모든 goal 도달 13/20; closest 충돌 5/20, goal 도달 12/20, 충돌 없이 완료 11/20; circle 완료 20/20이다. 공통 무충돌 완료 11개에서 6.193 / 6.185 / 6.252 s도 맞는다.
- five-robot GIF 경로의 기본 audit는 recorded poses만 검사한다 (`polygon_validation.py:24`, `cspace_cbf_5robots.py:332`). 그러나 저장된 전체 경로에 별도의 충분조건을 적용하니 **1500개 interval 모두 안전하다는 양의 하한**을 얻었다. 각 step의 저장된 모든 pair 최소거리에서 두 로봇의 가장 큰 motion bound 합을 뺀 최솟값은 **+14.0371448 mm**, 미확정 interval은 0이었다. 따라서 이 경로의 선형 pose 보간에서 충돌이 있을 것이라고 의심할 근거는 없다. 18.7934 mm는 기록된 pose 최소값이며, 여기서 얻은 전체 보간 경로의 보수적인 하한은 14.0371 mm다.
- Ground-truth polygon의 음수 지표가 translational penetration depth가 아니라 deepest vertex depth라는 설명, CPU/GPU baseline 계측 범위, circle barrier의 squared units, continuous-time 보장과 sampled control의 차이가 이제 명시되어 있다. 모두 유의미한 수정이다.

위 five-robot 확인은 다음 계산으로 재현했다:

```python
z = np.load('results/five_random_trajectory.npz')
x = np.concatenate([z['x'], z['final'][None]])
rho = np.array([np.linalg.norm(z[f'GT_{i}'],axis=1).max() for i in range(5)])
dx = np.diff(x,axis=0)
m = np.linalg.norm(dx[:,:,:2],axis=2) + abs(dx[:,:,2])*rho
lower = z['d'] - np.sort(m,axis=1)[:,-2:].sum(axis=1)
print(lower.min(), sum(lower<=0))  # 0.01403714481148401, 0
```

## 제출 전 우선순위

1. Docking 결과/형상 버전을 고정하고 모든 수치와 충돌 판정을 통일한다.
2. Intersample audit에 미확정 분류와 검증된 거리 하한을 추가하고 실행 결과에 남긴다.
3. Certification timer 범위를 바로잡고, 진행 중 dt sweep을 최종 보고한다.
4. 실험 확장이 가능하다면 matched-budget comparison과 강한 geometric baseline을 우선한다. 다양한 shape seed와 hardware는 그다음 확장이다.

스냅샷 SHA-256 (검토 종료 시; 이후 병행 수정은 별도):

```text
tro/main.tex                    4BBDD6F3C8F02A82AD8FBE4A30422A15D22A195D0B9EC7DAD199B20A4B166E27
cspace_sdf_cbf_compare.py       A9B327AD1FD051E11A05BFABCEDEE7050064C7FF3697B757AD6EBE6373B6D1F5
cspace_experiments.py           BCFC1EADA0F674860CF94F855CA0C3A1F956A437FD39ED4AC5E48E436E6762B8
dock_shapes.py                  A87BF088021CEBE3E4B073AECBF0F0A56E5D2A8BAFF8BBEFC619818CD4DA7156
results/dock2.json              EE1B59AA59088237700BFB93910CC040DD317092AD37F347A6B16B123CB01C44
results/stats.json              7E13D67CFB9B50CAABF207FA78FF4FC13B171B8B18DDA57C038252B1B6554E54
```
