# 재작업 진행 상황: safe set 이산화 없는 방식 (2026-09-29 밤, 사용자 부재 중 자동 실행)

## 사용자 지시 (그대로)
"smooth SDF를 safe set에서 이산화하지 않는 방식(지금 방식)대로 계속 가고, 논문 모든 실험 재실험하고, 바꿀 부분 바꿔.
비주얼(컬러 테마 등)은 이미 완성했으니 그런 건 냅둬. 크레딧 다 쓸 때까지 자동화해서 다 끝내 놔."
- 5분마다 cron(73011f74)이 이 파일을 보고 재개한다.

## 방법 (새 cover = "summed-field barrier + Bernstein 계수 조건")
- safe set: H = { q : phi_B(T x) >= l_B for all x in S_A } (S_A = {phi_A = l_A}, 연속, box 없음). = 기존 Thm cover의 연속 조건.
- 강제: box Q마다 lifted barrier  h~(x) = phi_A(x) - l_A + phi_B(T x) - l_B  (S_A 위에서 h~ = h).
  Psi~(x,u) = d/dt phi_B(T x) + gamma h~(x) >= 0 for all x in Q  <=  Bernstein 계수 조건 (box당 4^n개, u에 affine)
  - phi_A: box 위 정확한 tricubic/bicubic Bernstein 계수 (오프라인)
  - phi_B(T x): 중심 Taylor 2차 모델 + 3계 도함수 majorant M3 나머지
    - 값 나머지: gamma M3 r^3/6 (상수), gradient 나머지: (M3 r^2/2) sum_j L_j |u_j| -> epigraph t_j >= |u_j| (QP 유지, 정확한 동치)
  - phi_A - l_A 항이 box 밖으로 튀어나온 점을 보상 -> 1차 보수성 상쇄 (standoff O(r^2))
- 보수성은 전부 u 쪽 (허용 입력 집합). Prop feas(u=0 feasible)는 "모든 활성 계수 c_k >= 0일 때"로 약해짐 -> 본문에 명시.
- pruning: box의 h~ 하한 >= eta 이면 생략 (점별 연속성 논증). cluster pruning도 h~ 기준 (phi_A 하한 bmin 포함).
- 코드: prototype_3d/star_bern_timing.py (3D 정적 장애물 버전, 검증됨: star tube slack 0, 16 ms vs cover 21 ms)

## TODO (순서대로; 완료 시 [x])
- [x] 0. 공용 모듈 prototype_3d/summed.py: 2D(SE(2) 두 몸체) / 3D(링크-정적, 링크-링크) 행 생성 + QP(epigraph) 헬퍼
- [x] 1. star tube: startube_setup.py에 method 'summed' 추가 -> 결과 json, star_paper.png 재생성(스타일 그대로)
- [x] 2. Franka 30 trials: franka3d.py -> franka_results.json; franka_timing.py -> \timeOurs 등 매크로
- [x] 3. dual arm 20 trials: dual3d.py -> dual_results.json; dual_fig.py -> dual_paper.png; 타이밍 \timeDual*
- [x] 4. planar docking: dock_cover.py (summed) -> Table II, Fig 2 (make_paper_figs dock2)
- [x] 5. five robots: five_cover.py stats/swap (summed) -> Table III, Fig 3 (make_paper_figs five)
- [x] 6. maze: maze_cover.py (summed) -> 본문 숫자, rocket_maze.png
- [x] 7. Fig 1(b) overview: 활성 corner 점 -> 활성 box 표시 (스타일 유지)
- [x] 8. 본문: Abstract/Intro/Contrib/Sec cover 재작성(정리·증명)/Filter(Prop feas)/실험 숫자/Limitations/Conclusion/머리 주석
- [x] 9. 컴파일(12쪽 이하, ?? 없음), 숫자-json 대조표, 그림 확인
- [x] 10. 최종 보고서 (이 파일 하단) + 메모리 갱신

## 로그
- 02:58 star tube 비교(최초 구현, cubic 64행/box): corner 20.8 ms med, 7076 rows, gap 15.06 mm / summed gap 5.8 mm, slack 0, 45 ms
- 03:05 공용 모듈 완성: prototype_3d/summed.py (2차 Taylor 모델 + 2차 Bernstein, 3D 27행/box, 2D 9행/box), summed3d.py
  - 수치 검증: summed_check.py (3D, 78,880 샘플, 최소 여유 +0.019), summed_check2d.py (2D, 92,340 샘플, +0.0069) -> 하한 유효
- 연결 완료: startube_setup.py('ours'=summed, 'corner'=이전), franka3d.py, dual3d.py, franka_timing.py, dock_cover.py,
  cspace_experiments.simulate_team('summed'), five_cover.py, maze_escape.simulate(summed_fn), maze_cover.py
  - 이전 결과 백업: *_corner.json, tro/figs_corner_backup/
  - QP 주의: 고정 몸체(v_max=0)의 입력은 제거해야 함 (NNLS 발산) -> simulate_team에서 처리
- 03:25 docking summed: 5 mm box gap 1.97 mm (corner 6.8), 2.5 mm box 1.60 mm (corner 3.9); 두 level 합 1.67 mm
  -> 이산화 보수성 5 mm에서 ~0.3 mm (2차). slack 0, 6.0 ms med, 최대 384 box / 3456 행. dock_summed_results.json
- 03:3x 문제 발견: star tube에서 팔이 튜브 가장자리에서 멈춤 (gradient 나머지 항이 너무 큼). 원인 2개, 둘 다 수정:
  1. M3 bound가 너무 넓은 창(C^2 majorant, 축당 6 cell) -> A쪽: box ball이 닿는 cell 최대(오프라인 상수, box_bound),
     B쪽: 반 cell 격자 trilinear majorant (LinearMajorant, 축당 ~2 cell, Lipschitz)
  2. 관절별 |u_j| epigraph가 실제 점 속도를 3-10배 과대평가 -> 쌍별 상대 twist epigraph a >= |E u|
     (3D: V,Omega 6개/쌍, 2D: 3개/쌍), |ydot| <= sum|V_a| + (|c-p|+r) sum|Omega_a|
  - 재검증 통과 (3D 129,380 샘플 +0.0066, 2D 92,340 샘플 +0.0006)
  - 인터페이스: rows(..., Sb, eta, umax, auxmax) -> (A, T, C, kept, h_low); solve(u_nom, A, T, C, E, G_in, h_in, z)
- docking v2: 5 mm gap 1.66 mm (= 두 level 합 1.67 -> 이산화 보수성 ~0), 2.5 mm 1.58, tv 5.9, 6.0 ms, slack 0
- maze v1(구 bound): 32.5 s, 구간 검증 0.85 mm, h_low >= 2.7 mm, 1.9/3.5 ms, 126 box -> v2 재실행 중
- maze v2 완료: 32.5 s, 23.53 m, 구간 검증 >= 0.23 mm, h_low >= 2.09 mm, slack 0, 41% 개입, 2.05/3.61 ms,
  최대 127 box / 1143 행 -> tro/figs/rocket_maze.png 교체 (results/reference_maze_summed/)
- 본문: Sec. IV 전면 재작성(연속 safe set Thm, lifted barrier, Prop coef/tight, Thm invariance), 필터 절(Prop feas 약화),
  초록/서론/기여/Fig1 캡션, Sec. III kappa 정의 이동. 컴파일 OK 12쪽. (주의: bash heredoc은 \\ 를 \ 로 바꿈 -> LaTeX는 Edit/Write로)
- star v2 (bound 개선 후): gap 5.4 mm지만 여전히 튜브 모서리에서 정지. 진단: 손이 모서리 근처(h_low 0.6 mm),
  1/2 M3 r^2 > |grad| -> box 안에서 gradient 방향이 불확실해 물러나는 것도 인증 불가 -> 얼어붙음(안전, 교착)
- 해결: 온라인 적응 세분화 (summed.rows): 활성 box 중 1/2 M3_B r^2 > 0.1 이면 표면을 포함하는 2^n 자식으로 (최대 2단계).
  rows() 반환: (A, T, C, n_box, h_low, active_idx). 검증: summed_check_children.py (자식 데이터 오차 1e-16,
  부모 표면점 3555개 모두 자식에 포함), 계수 하한 검증 재통과.
- five v2(세분화 전) 참고: 20 inst 충돌 0, 최소 1.0 mm, 12/20 도달, pooled 4.8/9.4 ms, 최대 2214 행;
  swap 12.41 s, 0.76 mm, p95 58 ms (무거움) -> 최종 코드로 재실행 중
- star v3 (세분화 전면, near=inf): **목표 도달 4.08 s**, gap 2.98 mm, slack 0 — 그러나 210 ms/step, 5509 box
- 세분화를 h 하한 < 1 cm인 box로 제한 (REFINE_NEAR=.01): 궤적 자세에서 build 16 ms (세분화 없음과 동일)
- 최종 파라미터: REFINE_THETA=.1, DEPTH=2, NEAR=.01 (논문 Sec. IV Refinement 문단에 기재)
- star v4 (near 1cm, 27행): 미도달, 1.26 mm, 141 ms, 16.5만 행 -> 행 수가 병목
- 행 모드 비교 (summed.ROW_MODE): 'bernstein' 3^n행/box, 'single' 1행/box(속도항 전부 norm -> 접선 미끄러짐 막힘,
  Franka trial 1 미도달), 'vertex' 2^n행/box (속도의 d-선형항은 꼭짓점에서 정확, 2차항만 penalty) -> **vertex 채택**
  - 검증 통과 (3D 241,420 / 2D 92,340 샘플)
- 근본 trade-off (보고서에 명시): 경계 근처 허용 접선속도 ~ gamma h / sigma -> 속도비례 standoff, 진행 느림.
  theta 0.03이면 Franka trial 1 도달하지만 p95 200 ms. theta=0.1 유지.
- star (vertex, near 1cm): 미도달(손이 별 안쪽까지 들어가지만 멈춤), 1.53 mm, 101 ms -> star만 near=inf로 재실행 중
  (SUMMED_REFINE_NEAR 환경변수). 백업: startube_setup_results_vertex_near1cm.json, star_paper_vertex_near1cm.png
- 04:4x 실행 중(vertex 최종): Franka 30, dual 20, five stats+swap+dock, maze, star(near=inf)
- 최종 결과(vertex): maze 32.51 s / 0.48 mm 구간검증 / 2.2,3.6 ms / 126 box 504행 (본문 반영)
  five stats: 충돌0, 0.60 mm, h_low>=0, 12 도달, 최대 255 box/1020 제약, 공통12 평균 6.98 s 16.48 m
  swap: 12.46 s, 0.84 mm, h_low +0.17, 2716 제약, setup levels .29-.35 mm, boxes 1240-1426
  dock: 0.698 m, 1.68 mm (=레벨합), tv 6.4, 401 box/1604, 2.5mm: 1.57 mm (본문/표 반영, 시간 제외)
- TIME TODO (격리 측정 후 기입): Table II cover time(현재 4.9 그대로), Table III 시간 행, swap 문단 시간,
  \timeOurs/\timeFranka*/\timeDual* (franka_timing.py ours dual), star 문단 시간
- 06:5x 3D 결과: Franka 충돌0, 인증 최소 1.1 mm (corner 8.4), 중앙 5.7, 도달 17 (corner 23), slack 0, 4374 box/34,992
  dual 충돌0, 1.6 mm (corner 10.3), 도달 15 (corner 17), 활성 8, 3cm이내 16, 6730 box/53,840
  star (vertex, near=inf): 4.13 s 도달, 3.03 mm, slack 0, 5410 box/43,280, 114/167 ms(동시실행) -> 그림 star_paper.png OK(통과)
- 본문 반영 완료: docking/five/maze/Franka/star/dual/초록. placeholder: TIMEFIVE, TIMESWAPMED, TIMESWAPP95, TIMESTARMED,
  TIMESTARP95, Table II cover time(4.9), Table III time 행, 5대 setup 시간(34 s), 오프라인(9.5 s), \time* 매크로
- 격리 타이밍 결과: Franka 13.3/46.5 ms (corner 9.4/18.7), dual 7.6/78 ms (corner 3.8/28.5) -> 매크로 반영
- dual 그림: trial 7, 7.62 s 도달 -> 캡션 "slide past each other ... 7.6 s" 반영. 페이지 12 맞춤(결론 압축 등)
- 07:0x 격리 타이밍 체인 실행 중 (bd944xms1): dual_fig -> franka_timing ours dual -> five stats -> swap -> dock -> star
- 본문: Prop coef를 vertex 형태로 재작성 (값: Bernstein min, 속도: 꼭짓점 + sigma_V, sigma_Omega). 현재 13쪽 -> 줄여야 함
- 03:35 (이전) 실행: star tube v2, five stats+swap v2, Franka 30, dual 20 -> 코드 변경으로 중단/폐기. 남은 것: timing(franka_timing.py ours dual),
  그림(dock2/five/overview via make_paper_figs, dual_fig), 실험 절 숫자, Limitations/Conclusion

## 최종 보고 (2026-09-29 07:35 완료)
논문: tro/main.tex (12쪽, 미정의 참조/placeholder 없음), PDF 사본 tro/main_summed_draft.pdf
이전 본문: tro/main_before_summed_2026-09-29.tex, 이전 그림: tro/figs_corner_backup/

### 방법 (Sec. IV 전면 재작성)
- safe set: S_A 위 모든 점에서 phi_B >= l_B (연속, box 없음) + Thm (충돌 없음)
- 강제: box마다 lifted barrier h~ = phi_A - l_A + phi_B - l_B; 값 = 2차 Bernstein 계수 최소(상수),
  속도 = box 꼭짓점 2^n개에서 정확 + 2차 나머지; Taylor 나머지(M_A box별, M_B Lipschitz majorant)는
  상대 twist 보조변수 a >= |(V,Omega)|로 입력 쪽에만 (QP 유지, 정확한 동치); 접촉 근처 적응 세분화
- 필터: Prop feas 약화(경계 근처 2차 band 제외), Cor/Thm invariance

### 결과 (모두 충돌 0, slack 0, 인증 하한 >= 0)
| 실험 | 이전(corner) | 새 방식 | 출처 |
|---|---|---|---|
| docking gap | 6.8 mm | 1.7 mm (= 두 level 합), 6.2 ms | prototype_3d/dock_summed_results.json, dock_summed_iso.log |
| 5대 20 inst 최소거리 | 2.5 mm | 0.6 mm, 12/20 도달, 2.0/10.2 ms | results/five_cspace_summed_stats_5.json |
| swap | 3.8 mm | 0.84 mm, 12.46 s, 2.8/22 ms | results/five_cspace_summed_swap_5.json |
| 미로 | 1.5 mm 구간검증 | 0.48 mm, 32.5 s, 2.2/3.6 ms | results/reference_maze_summed/ |
| Franka 30 | 8.4 mm, 23 도달, 9.4/18.7 ms | 1.1 mm, 17 도달, 13.3/46.5 ms | prototype_3d/franka_results.json, timing_results.json |
| star tube | 15.1 mm | 3.0 mm, 4.1 s 도달, 96/109 ms (전 box 세분화) | prototype_3d/startube_setup_results.json |
| 두 팔 20 | 10.3 mm, 17 도달, 3.8/28.5 ms | 1.6 mm, 15 도달, 7.6/78 ms | prototype_3d/dual_results.json |

### 사용자가 판단해야 할 것
1. trade-off: 훨씬 가까이 가지만(safe set 정확) 접촉 근처 접선 속도가 제한되어 도달 수 감소, 비용 증가.
   조절 손잡이: REFINE_THETA(0.03이면 도달↑, 비용↑↑), box 크기, horizon.
2. star tube만 SUMMED_REFINE_NEAR=inf (모든 활성 box 세분화). 1 cm 제한이면 모서리에서 멈춤. 본문에 명시함.
3. 수치 검증 스크립트: summed_check.py, summed_check2d.py, summed_check_children.py (모두 통과)
4. Overleaf 복사: tro/main.tex, tro/figs/*.png (dock, five, overview, rocket_maze, star, dual 재생성)
