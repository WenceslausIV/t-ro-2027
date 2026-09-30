검토일: 2026-09-26. Reviewer 3 — 신규성, 관련 연구, 선행/동시 투고 원고와의 관계, T-RO 적합성.

**모의 판정: 현재는 weak reject / 큰 수정 후 재심사 쪽이다. 다만 ICRA 논문의 분량만 늘린 원고라는 이유로 거절할 논문은 아니다.** 핵심 연구 기여는 남아 있다. 이번 수정으로 related work, exact-arithmetic 한계, 안전성과 목표 도달의 구분, ICRA 2027 원고의 기여 귀속이 상당히 좋아졌다. 현재의 가장 큰 문제는 그 기여가 최신 대안 대비 얼마나 유리한지를 보여 주는 증거와, 본문 전체의 주장·실험 동기화다. 아래의 신규성 평가는 reviewer 판단이며, 수학적 오류 판정과 구분한다.

검토 범위는 현재 `tro/main.tex` 전체와 `references.bib`, 로컬 ICRA 2026 PDF의 방법 부분, ICRA 2027 PDF의 representation/enclosure/CBF 증명, 아래에 연결한 일차 문헌과 T-RO 공식 지침이다. 실험 전체를 재실행하거나 PDF를 컴파일하지 않았다. 줄 번호는 검토 시점 파일 기준이다.

**1. [중요·신규성 논증] ‘한 개 scalar constraint’와 ‘형상 복잡도에 독립적인 온라인 계산’을 구분해야 한다.**

위치: `main.tex:196–207`, `main.tex:312–315`, `main.tex:444–450`, `main.tex:1113–1117`.

Convex decomposition을 이용하는 통상적인 직접 구현에서는 piece pair마다 제약을 만들 수 있다. 그러나 이것이 반드시 많은 QP 제약을 필요로 한다는 일반론은 아니다. [Molnar and Ames, Composing Control Barrier Functions for Complex Safety Specifications](https://arxiv.org/abs/2309.06647)는 Boolean 조합을 하나의 미분 가능한 CBF로 표현한다. 이미 인용한 [Wu et al.](https://arxiv.org/html/2502.16293v1)의 Section V도 single CBF를 사용한다. 따라서 ‘상대 논문은 많은 제약, 우리는 하나’만으로 차별성을 방어하면 반박을 받기 쉽다.

차별성은 더 구체적으로 적는 편이 강하다. 원고의 cubic field는 pose마다 로컬 계수만 평가하므로, 고정된 field를 준비한 이후 온라인 평가량이 원래 boundary span 수와 직접 비례하지 않는다. 반면 smooth composition은 QP row를 하나로 줄여도 constituent geometry 값을 계산해야 한다. **제약 수, geometry evaluation 비용, 오프라인 시간/메모리를 별개의 축으로 비교하라.** 이는 제안 방법과 composition이 동등하다는 뜻이 아니며, 모든 composition이 bounded-input 환경에서 자동으로 feasible하다는 뜻도 아니다.

현재 `main.tex:204–207`의 Wu 설명은 대체로 맞지만 ‘polygons’를 ‘convex polygons’로 한정하면 더 정확하다. 해당 원문의 Assumption 1은 convexity를 명시한다.

**2. [중요·관련 연구] C-space의 연속 영역을 수학적으로 검증하는 기존 연구와의 관계가 빠져 있다.**

위치: `main.tex:215–241`; `references.bib` 전체.

신경망 C-space barrier와 고전 C-space construction을 추가한 것은 적절하다. 그러나 이 논문의 핵심을 ‘continuous-domain certificate’로 옮겼으므로, 그에 가까운 certification 문헌도 논의할 필요가 있다. [Dai et al., Certified polyhedral decompositions of collision-free configuration space, IJRR 2024](https://journals.sagepub.com/doi/10.1177/02783649231201437)는 C-IRIS를 통해 다항식/SOS 기반으로 collision-free configuration regions를 검증한다. 원고의 citation list에 없다.

그 방법이 이 논문을 선취했다는 뜻은 아니다. C-IRIS는 인증된 convex region들의 구성이고, 여기서는 fitted smooth field의 contact-set exclusion을 인증하여 reactive CBF를 만든다. 이 차이를 짧게 설명하면 왜 ‘기존의 certified C-space + controller’가 곧바로 대체하지 않는지 명확해진다. C-IRIS 전체를 실험 baseline으로 구현하는 것을 필수 요구하지는 않는다.

**3. [중요·reviewer 판단] 남아 있는 실험은 최신 형상 기반 안전 필터보다 낫다는 근거가 부족하다.**

위치: `main.tex:329–338`, `main.tex:693–698`, `main.tex:997–1053`, `main.tex:1163–1167`.

이제 단일 closest-pair baseline의 한계와 구현별 timing 차이를 솔직하게 명시했다. 이것은 해석 오류를 고쳤지만 비교의 공백 자체를 채우지는 않는다. 저자 스스로 nonconvex switch에서 유효한 nonsmooth constraint들을 처리하지 않는다고 인정한 baseline의 실패는, 제안 field를 사용할 필요성을 입증하는 강한 비교가 되기 어렵다. 원래 ICRA 방법의 알려진 한계를 재현한 결과로는 유용하다.

T-RO 기여를 설득하려면 적어도 하나의 제대로 된 형상 대응 대안을 같은 모델·정확도·control budget 아래 비교하는 것이 가장 효과적이다. 예를 들어 convex decomposition + 적절한 nonsmooth/feature-aware constraints, 또는 convex pieces의 smooth composition을 사용할 수 있다. [Chen et al.의 최신 exact signed-distance/nonsmooth CBF 원고](https://arxiv.org/abs/2608.02886)는 이미 관련 연구에 포함되어 있으므로, 그 계열을 실제 비교하지 않은 이유와 예상 tradeoff를 최소한 설명해야 한다. 제약 수만이 아니라 certified/observed clearance, goal completion, preprocessing time, online time을 함께 평가하라.

이전 검토에서도 baseline 부족은 지적되었다. 이번 재검토의 판단은 **limitation에 적은 것만으로 이 출판 수준의 증거 부족까지 해소되지는 않았다**는 것이다. 반드시 특정 논문 하나를 구현해야 한다는 규정은 아니다.

**4. [구체적 확인] ICRA 2026에 대한 실질적 확장은 있고, ICRA 2027과도 기여를 구분할 수 있다.**

위치: `main.tex:179–190`, `main.tex:235–241`, `main.tex:457–467`.

로컬 ICRA 2026 PDF pp. 4–5는 workspace closest-point BP-SDF barrier와 bounding-circle inflation에 의한 translation-only C-space를 다룬다. 이에 비해 현재 원고는 rotating shape pair의 SE(2) field, contact-set branch-and-bound, activation-band/regularity verification으로 바뀐다. 이것은 단순 추가 실험이나 증명 보충 이상의 변화다.

로컬 ICRA 2027 PDF pp. 3–6과의 공통 내용은 periodic cubic boundary fit, B-spline-to-Bernstein 변환, hull/winding enclosure 논리다. ICRA 2027의 중심은 disk robot의 boundary-dependent semantic margin과 온라인 degree-six coefficient constraints다. 현재 T-RO의 중심은 두 비볼록 강체에 대한 오프라인 3-parameter contact certification과 한 개 smooth pair field다. 현재 Relation 문단은 공통 기반을 숨기지 않고 출처를 밝히므로 개선되었다.

다만 ‘first-contact lemma’, 상대 좌표 derivative, driftless zero-input feasibility를 각각 대형 독립 신규성처럼 나열하면 contribution이 잘게 부풀려졌다는 인상을 줄 수 있다. 이들은 유용한 구성 요소이며, **핵심 기여는 이들을 결합한 인증 알고리즘과 그 실용적 비용/보수성**이라고 정리하는 편이 강하다. lemma가 새로 번호 붙었다는 사실과 연구적으로 새로운 원리가 생겼다는 사실은 다르다.

동시 심사 중인 ICRA 2027 원고와 도구를 공유한다는 사실만으로 중복 투고를 판정할 수 없다. Cover letter에는 두 원고가 공유하는 결과와 별개 결과를 구체적으로 나누고, 편집자·심사자가 관련 원고를 확인할 수 있게 준비하는 것이 좋다. 단순 인용으로 실제 중복 결과가 허용되는 것은 아니다.

**5. [주장 정밀성·재수정 필요] 안전성/보수성에 관한 몇 문장이 여전히 본문의 정확한 조건보다 강하다.**

- `main.tex:163–166`: ‘any collision-free start’는 Corollary의 `main.tex:529–532`보다 강하다. 초기 fitted bodies separation과 active barrier의 `h>=0`가 필요하다. 물리적으로 collision-free인 모든 pose가 fitted safe set에 포함된다는 증명은 없다. contribution 문장에 초기 safe-set 조건을 넣어야 한다.
- `main.tex:755`는 서로 다른 field의 낮은 level을 바로 ‘less conservative’라고 부른다. 바로 뒤 `main.tex:763–767`에서는 작은 level만으로 큰 safe set을 결론 낼 수 없다고 정확하게 썼다. 수준 차이를 이야기할 때는 ‘smaller certified offset’로 쓰고, 실제 보수성은 safe-volume/clearance/task outcome으로 별도 제시하라. `main.tex:765`의 ‘is ... the largest amount’도 maximum 자체가 아니라 그것의 upper bound라는 표현이 정확하다.
- `main.tex:1139–1142`는 C0 docking에서 실제 충돌과 약 36 mm 위반을 보고하지만, `main.tex:948–951`, `main.tex:969`는 다른 값과 비충돌을 말한다. 어느 결과가 최신인지 raw results와 그림을 동기화해야 한다. 이 모순은 독립적으로 확인했고 실험 담당 reviewer에게 전달했다.
- `main.tex:619–624`는 rounding error가 ‘small compared with tolerance’라고 단정하면서 동시에 unquantified/unbounded라고 한다. 일반적인 기대와 검증된 수치 경계를 구분하라. 정확산술 조건 명시는 잘 수정되었으며, 여기서는 남은 단정만 문제다.

**형식과 정책 확인.** [현재 T-RO 공식 안내](https://www.ieee-ras.org/publications/t-ro/t-ro-information-for-authors/)는 2025년 1월부터 evolved-paper 별도 분류를 없애고 substantive new results를 요구한다. 최초 투고 18쪽 이하, regular-paper abstract 200단어 이하, double-anonymous 규칙이 적용된다. 현재 abstract는 공백 기준 190단어이며 author anonymization은 활성화되어 있다. 컴파일하지 못했으므로 18쪽 충족과 최종 float 가독성은 미확인이다. T-RO가 하드웨어 실험을 보편적으로 의무화한다고 해석할 근거는 찾지 못했다. 시뮬레이션만으로 게재가 불가능하다는 주장은 하지 않는다.

`references.bib:62–66`의 Sailing Through Point Clouds 서지는 권/호/쪽이 빠졌다. [저자 공식 구현 저장소의 권장 citation](https://github.com/BolunDai0216/SailingThroughPointClouds)에 따르면 RA-L 9(9), 7731–7738 (2024)다. 새로 추가된 Long/Chi/Wu/Chen 문헌의 존재와 현재 제목은 확인했으며, Chen 2026의 v1과 v2 제목이 다르므로 예전 v1 HTML 제목만 보고 현 bibliography를 오류로 판단하면 안 된다.

**강점과 게재 전망.** 학습 target의 오차와 안전성 인증을 분리한 설계, nonconvex pair rotation을 직접 처리하는 표현, 같은 Bernstein machinery로 band/regularity까지 검사하는 일관성이 강점이다. 현 원고는 T-RO를 목표로 할 만한 연구 주제를 갖추었다. 그러나 지금 추천한다면 accept보다 큰 수정이 필요하다는 쪽이다. 가장 효과적인 개선은 최신 결과·그림의 일치, 강한 형상 baseline 하나, 그리고 ‘한 scalar barrier’보다 ‘인증 가능한 저비용 로컬 field evaluation’에 초점을 맞춘 novelty 설명이다.
