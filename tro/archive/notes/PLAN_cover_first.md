# 플랜: cover를 주 방법으로, tabulated를 compile 변형으로 (옵션 A)

작성일 2026-09-28. 대상 논문: `tro/main.tex` (현재 11페이지, 한도 12페이지).
실행자는 이 문서를 위에서부터 순서대로 따른다. 각 단계에 "완료 조건"이 있다.

## 0. 확정된 방향과 규칙

- **방향 (옵션 A):**
  - 주 방법은 surface cover다. 2D와 3D, 강체와 다관절에 모두 쓴다.
  - tabulated는 "평면 강체 쌍을 매끄러운 barrier 하나로 compile한 변형"으로 줄인다.
    - 장점: 제약 수 고정, 최악 step 시간 보장, 정적 장애물 여러 개를 field 하나로 처리.
    - 단점: 보수성(training target이 정함).
  - 5대 로봇과 docking은 두 방법을 나란히 보여준다.
  - 미로는 cover로 먼저 돌려보고 결과를 보고 정한다.
- **사용자 규칙 (반드시 지킬 것):**
  - 답변은 항상 한국어로 한다.
  - 요청하지 않은 변경은 하지 않는다. 그림의 색, 굵기, 배치를 "개선"하지 않는다. 바꿔야 하면 먼저 묻는다.
  - 보장은 연속 표면 전체에서 성립해야 한다. 샘플 기반 인증 주장은 금지.
  - 12페이지 이하.
  - LaTeX 수정은 Edit 도구로 한다. heredoc은 `\r`, `\C` 같은 백슬래시를 망가뜨린다.
- **컴파일:**
  - 명령: `"$TEMP/claude/c--Users-siwon-Desktop-swarm-bp-sdf-cbf/15d6a89a-2317-450e-94da-60e031c4f7fd/scratchpad/tectonic/tectonic.exe" -o <scratchpad>/build tro/main.tex`
  - 페이지 수는 `pypdf`로 센다.
- **사용자 전달물:** Overleaf에는 `tro/main.tex`, `tro/references.bib`, `tro/figs/*.png`를 복사한다.

## 1. 결정 (2026-09-28 사용자: "나머지 너가 알아서 해" → 아래 기본값으로 확정)

- Fig 1 (a)의 접촉 별표는 crimson 빨강으로 유지한다(사용자가 좋아함).
- (a)의 화살표는 B_i, B_j에 절대 닿지 않는다. 코드가 두 몸체에서 4 cm 이상 떨어진 경로를 계산해서 그린다. 새 그림에서도 이 규칙을 유지한다.

| # | 결정 | 기본값 (추천) |
|---|---|---|
| D1 | 제목 | "Certified Surface-Cover Control Barrier Functions for Non-Convex Rigid and Articulated Robots". 부제 없이 간결하게. 현재 제목 유지도 가능. |
| D2 | Fig 1 구성 | (a) first contact는 유지. (b)는 새 cover 그림: B_i의 인증 곡선 위 상자와 모서리 점, B_j의 level set. (c)는 현재 (b)의 C-space slice(tabulated)로 옮긴다. 현재 (c) M-set 산점도는 뺀다. 색은 지금 Fig 1 규칙(회색=기준/실제, 주황=인증 level set, 파랑=B_i)을 유지한다. |
| D3 | Fig 2 (docking) | (a)를 cover 결과로 바꾸고 (b) circle은 유지한다. (c),(d) 곡선은 cover, tabulated, closest, circle 4개. C⁰ Bézier는 그림에서 빼고 표와 본문에만 남긴다. 선 스타일과 색 배정은 사용자에게 보여주고 확정한다. |
| D4 | Fig 3 (5대 로봇 swap) | cover 실행(`results/five_cover_swap_cover_5.npz`)으로 바꾼다. 그림 스타일(fading, 레이아웃, (b) 색 차콜/청록)은 그대로 둔다. |
| D5 | Fig 4 (미로) | 2-3단계의 cover 실행이 성공하면 cover로 바꾼다. 실패하면 tabulated를 유지하고 한 문장으로 이유를 적는다. |

## 2. 실험·코드 (본문 숫자가 여기에 의존하므로 먼저 한다)

### 2-1. 2D cover에 certified cluster pruning 추가

- **파일:** `prototype_3d/dock_cover.py`의 `make_cover_fn`
- **배경:** 3D에는 cluster pruning이 있지만 2D에는 상자별 pruning(`ACT = .05`)만 있다. 그래서 모든 상자에서 field를 평가한다. swap 최악 66.7 ms의 주원인으로 추정한다.
- **방법:**
  - 본문 Sec. cover의 "Pruning" 규칙을 그대로 쓴다.
  - 상자를 16–32개씩 묶은 cluster에 중심 c̄와 반지름 R(상자 ball들을 포함)을 준다.
  - φ_B(c̄) − G·R − ½·M̄·r² − l_B ≥ η이면 cluster 전체를 건너뛴다.
  - G와 M̄는 cluster(와 majorant support)가 닿는 cell들의 인증된 상한이다. 3D 구현(`prototype_3d/franka3d.py` 또는 `interactive/server.py`의 `body_rows`)에 있는 방식을 재사용한다.
- **완료 조건:**
  - pruning 전후 궤적이 같다(최소 거리 차이 0.1 mm 이하).
  - 최악 step 시간이 줄어든다.

### 2-2. 평면 시간 측정을 Franka와 같은 방식으로

- **새 파일:** `prototype_3d/planar_timing.py`. `prototype_3d/franka_timing.py`를 본뜬다.
- **대상:** 5대 로봇 20 instance, swap, docking. 두 방법 모두 잰다.
- **방식:** 단일 스레드, 다른 작업 없이 격리해서 잰다. median과 95th percentile을 보고한다(Franka 표와 같은 기준).
- **완료 조건:**
  - `results/planar_timing.json`이 생긴다.
  - cover p95가 10 ms 제어 주기 안에 드는지 확인한다. 안 들면 그 사실을 그대로 기록한다(본문 Limitations에 쓴다).

### 2-3. 미로를 cover로 실행 (가장 큰 위험 요소)

- **새 파일:** `maze_cover.py`. `reference_maze_escape.py`를 재사용한다.
- **벽:**
  - `designed_corridor_walls()`의 벽 두 영역(실제 형상 polygon)을 쓴다.
  - frame과 여유를 덮는 2D B-spline field(cell 1–1.25 cm)를 벽 polygon까지의 정확한 signed distance에 맞춘다.
  - level은 `tight_level_2d`(Prop. level)로 벽 polygon 위에서 인증한다.
- **로켓:**
  - `Rocket.setup()`의 실제 형상 polygon에 `body_field`, `tight_level_2d`, `cover_2d`(5 mm)를 적용한다.
- **barrier:** 로켓 cover를 벽 field에 대고 건다. 벽은 정지해 있으므로 로켓 속도 열만 있다.
- **그대로 쓸 것:** waypoint(`few_waypoints`), nominal, γ, dt, step 수, 구간 검증(`G.physical_gaps`, interval audit).
- **출력:**
  - `results/reference_maze_cover/results.json`과 `trajectory.npz`
  - D5가 cover로 정해지면 `tro/figs/rocket_maze.png`. 그림 스타일은 현재 `draw()`를 그대로 쓴다(G 구, 여백, 선 굵기 포함).
- **기록할 값:**
  - 두 level, 상자 수, 최대 barrier 수
  - 시간 med/p95
  - 검증된 최소 간격, min h
  - 탈출 시간, 경로 길이, 개입 비율, slack 여부
- **완료 조건:**
  - 출구에 도달한다.
  - 검증된 최소 간격이 0보다 크다.
  - slack이 없다.

### 2-4. docking cover 궤적을 그림용으로

- `prototype_3d/dock_cover.py`는 이미 `dock_cover_traj_5.npz`를 저장한다.
- `make_paper_figs.py fig_dock2`가 이 파일을 읽도록 한다. D3 기본값 기준이다.

## 3. 그림

작업마다 그리고 나서 이미지를 직접 열어 확인한다.

- **Fig 1** (`make_paper_figs.py fig_overview`, D2):
  - (b) cover 그림: 상자는 보이도록 확대(예: 6–8 cm)하고 caption에 "enlarged"라고 적는다.
  - 3등분 레이아웃과 제목 기준선은 그대로 둔다.
- **Fig 2** (`fig_dock2`, D3): 색과 선 스타일 표를 사용자에게 보여준 뒤 확정한다.
- **Fig 3** (`fig_five`, D4):
  - 입력을 `five_cover_swap_cover_5.npz`(traj, h, d)로 바꾼다.
  - 로봇 외곽선은 cover가 쓰는 인증 곡선이 아니라 필드 level curve다. 무엇을 그릴지 한 줄로 확인한다. 기본값은 지금처럼 실제 형상(회색) + 인증 B-spline 경계(색)다.
- **Fig 4** (D5): 2-3단계 결과를 보고 정한다.
- **Fig 5, 6:** 변경 없음.

## 4. 본문 (`tro/main.tex`)

### 4-1. 앞부분

- **Abstract:**
  - 첫 문장은 first-contact reduction이다.
  - 다음으로 주 방법 cover를 설명한다. 인증된 level set을 Bernstein cover로 덮고, 곡률 majorant로 모서리 barrier를 만든다. 2D와 3D, 강체와 다관절에 쓴다.
  - 그다음 tabulated를 소개한다. 평면 쌍을 barrier 하나로 compile하고 Bernstein B&B로 인증한다.
  - 결과 문장은 다음 순서로 쓴다: 5대 로봇에서 두 방법 모두 무충돌, cover는 7배 가깝게(2.5 vs 19.7 mm) → docking → Franka → 두 팔.
- **Intro (현재 135–144행):**
  - "certify in two ways"는 유지한다.
  - 순서를 cover(일반) → tabulated(평면 특수화, 제약 1개)로 바꾼다.
- **Contributions:**
  1. first contact와 인증 영역
  2. surface-cover barrier (n = 2, 3; 강체 SE(2)/SE(3), 다관절)
  3. tabulated compile 변형
  4. 평가 (평면에서 두 방법 비교 포함)
- **Relation to Prior Work:** 거의 그대로 둔다.

### 4-2. 방법 섹션 순서 교체

- **Sec. Regions:** 변경 없음.
- **Sec. cover (tabulated 앞으로 이동):**
  - 도입의 "For spatial rigid bodies…"를 "For bodies in ℝⁿ, n ∈ {2, 3}"로 일반화한다.
  - Prop. taylor의 "corners v₁…v₈"을 "2ⁿ corners"로 바꾼다.
  - 평면 강체용 미분을 한 줄 추가한다. 두 몸체가 움직일 때 c와 v_k의 속도는 v_i + ω_i J(·−p_i)이고, 이를 B의 frame에서 표현한다. `make_cover_fn`의 R6와 일치시킨다.
  - Pruning 문단은 2D에도 적용된다고 명시한다.
- **Sec. tabulated (뒤로, 축약):**
  - 제목은 "Compiling a Planar Pair into One Barrier"로 한다.
  - 남길 것: field 식, h, Prop. hdot, M 집합, Def, Thm, Alg. 1, Prop. bb(증명 축약), band와 regularity(3줄), Representations(3줄).
  - 정적 장애물 union을 field 하나로 처리한다는 문단은 이 변형의 장점으로 강조한다.
- **Safety filter:** "one per active planar pair, one per active box corner otherwise"의 순서를 바꾸고 문구를 맞춘다.

### 4-3. 실험

- **Planar setup:** cover 파라미터를 추가한다.
  - 몸체당 workspace field 1 cm cell
  - Prop. level로 실제 형상 polygon 위에서 인증한 가장 촘촘한 level(0.27–0.36 mm)
  - 5 mm 상자, ACT/η
- **Table I (representations):** 유지하되 본문을 줄인다.
  - tabulated field에 관한 표임을 밝힌다.
  - C² B-spline이 cover의 κ majorant에도 필요하다는 연결 문장을 하나 넣는다.
- **Docking:**
  - 본문 순서는 cover 먼저, tabulated 다음으로 한다.
  - Table II 열 순서: cover | tabulated (B-spline) | tabulated (C⁰) | closest | circle
  - 2-2의 시간을 반영한다.
- **5대 로봇 (Table III에 cover 열 추가):**
  - 숫자 출처: `results/five_cover_stats_5.json`. 같은 실행의 tabulated 재측정값을 쓴다.

    | | tabulated | cover |
    |---|---|---|
    | 충돌 | 0 | 0 |
    | 최소 실제 거리 | 19.7 mm | 2.5 mm |
    | min h | −2.4 mm | 0.0 mm |
    | 목표 도달 | 13 | 12 |
    | 시간 med/max | 0.74/5.1 ms | 1.14/11.7 ms |
    | barrier 최대 | 1개/쌍 | 927 |

  - 시간은 2-2의 격리 측정값(med/p95)으로 교체한다.
  - 기존 표의 tabulated 0.66/16.3 ms와 섞지 않는다.
  - 공통 12개 instance: 평균 6.86 s / 6.79 s, 경로 16.43 / 16.44 m.
  - 오프라인 비용: cover는 5개 몸체 합계 9.3 s(N에 비례), tabulated는 10쌍 × 11–16 s(N²에 비례).
  - swap 문단: `results/five_cover_swap_5.json`. 12.16 s, 3.0 mm, h +0.03 mm, 최대 barrier 4160.
- **Static:**
  - m개 장애물(0.37–0.45 ms, m과 무관)은 tabulated의 강점으로 남긴다.
  - 미로는 D5에 따라 쓴다.
  - sampling/ball-world 비교는 페이지가 모자라면 한 문장으로 줄인다.
- **Franka, 두 팔:** 변경 없음.

### 4-4. 끝부분

- **Limitations:**
  - Conservatism: tabulated는 training data가 정하고, cover는 l_A + l_B + r + κ항이 정한다.
  - 비용: cover의 제약 수(수천 개)와 최악 step 시간(2-2 결과).
  - Discrete time: 음수 h는 tabulated에서만 나타났다(cover는 min h 0.0 mm).
  - 공간적으로 변하는 level은 future work로 한 문장만 둔다.
- **Conclusion:** 순서를 cover → tabulated로 바꾼다.
- **파일 머리 주석:** OPEN ITEMS와 스크립트 목록을 갱신한다.
- **표현 주의:**
  - "patch마다 오차가 bound된다"고 쓰지 않는다. 인증되는 것은 φ < l 한쪽 조건뿐이다.
  - 필요하면 한 문장만 쓴다: "all certificates are cellwise, so their intervals shrink with the cell size".

## 5. 검증과 전달

1. 컴파일한다. 12페이지 이하, `??`와 undefined reference가 없어야 한다.
2. 본문의 모든 새 숫자를 json과 대조한다. 숫자 → 파일 목록을 답변에 붙인다.
3. 바뀐 그림을 모두 열어서 확인한다.
4. 사용자에게 한국어로 보고한다: 바뀐 것, 복사할 파일 목록, 남은 결정.
5. 사용자가 구조를 확정하면 메모리 `tro_unified_paper.md`를 갱신한다. 내용: "cover가 주 방법, tabulated는 compile 변형".

## 6. 권장 실행 순서와 위험

1. 2-1 → 2-2: 싸고, 시간 주장의 근거가 된다.
2. 2-3 미로: 가장 큰 위험이다. 벽 field 해상도와 상자 수에 따라 barrier 수와 시간이 커질 수 있다.
3. 사용자에게 D1–D5를 확인한다. 2-3 결과를 같이 보여준다.
4. 그림(3장) → 본문(4장) → 검증(5장).

- **페이지 위험:** 방법 순서 교체와 tabulated 축약으로 약 0.3쪽이 줄고, 평면 미분식과 표 열 추가로 약 0.2쪽이 는다. 12쪽을 넘으면 sampling/ball-world 문단과 Table I 본문부터 줄인다.

## 7. 진행 상황 (2026-09-28, Opus 5.5 실행)

모든 단계 완료. 컴파일 결과 12페이지, 미정의 참조 0.

- 2-1 cluster pruning: `prototype_3d/dock_cover.py make_cover_fn`. docking 궤적이 동일하다(nose 0.6932 m, 6.80 mm).
- 2-2 시간: 단일 스레드로 전부 재측정했다.
  - 파일: `results/dock2.json`, `stats.json`, `five_cover_stats_5.json`, `five_cover_swap_5.json`
  - swap은 1 cm cell, 5 mm 상자로 다시 돌렸다.
- 2-3 미로 cover: `maze_cover.py` → `results/reference_maze_cover/`
  - 32.5 s에 탈출, 검증 간격 ≥ 1.5 mm, 최대 barrier 470개, 1.2/2.3 ms
- 그림: Fig 1–4를 모두 cover 기준으로 바꿨다. 이전 그림은 `tro/figs_unused/`에 있다.
- 본문: 방법 섹션 순서 교체, "compiled field" 용어 통일, Table II와 III를 재측정값으로 갱신.
- 이전 본문: `tro/main_before_cover_first_2026-09-28.tex`
