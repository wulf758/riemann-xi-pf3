# H804 Unit-Deficit Moment-Curvature Barrier Audit

Classification: `unit_deficit_barrier_supported_with_finite_base`

## Derived Barrier

H803 needs

```text
T_{n+1}/T_n >= theta_n = C_{n-1}C_{n+1}/C_n^2.
```

The exact deficit is

```text
1 - theta_n = (8n^2 - 4n - 3)/(n^2(2n - 1)^2).
```

Since this is larger than `1/n^2` for `n>=2`, the stronger unit barrier

```text
T_{n+1}/T_n >= 1 - 1/n^2
```

implies the H803 inequality.

## Checks

| check | value |
| --- | --- |
| rows_tested | `125` |
| n_range | `[2, 126]` |
| all_unit_implies_threshold_exact | `True` |
| all_threshold_formula_agrees_reported | `True` |
| all_h803_conditions_available | `True` |
| unit_barrier_tail_start_midpoint | `8` |
| unit_barrier_failures_before_tail | `[2, 3, 4, 5, 6, 7]` |

## Tightest Rows

Minimum H803 threshold slack:

```json
{
  "n": 126,
  "observed_T_next_over_T_n": "9.999656805356967230530425041972033265E-1",
  "factorial_threshold": "9.998735247852608801527979681614431176E-1",
  "unit_floor_1_minus_1_over_n2": "9.999370118417737465356512975560594608E-1",
  "observed_deficit": "3.431946430327694695749580279667350000E-5",
  "observed_deficit_constant_n2": "5.448558152788248098972033651999884860E-1",
  "threshold_deficit": "1.264752147391198472020318385568824000E-4",
  "threshold_deficit_formula": "1.264752147391198472020318385568823770E-4",
  "threshold_formula_agreement_abs": "2.295150313934596270250806346605279502E-38",
  "threshold_deficit_constant_n2": "2.007920509198266694179457468929064982E+0",
  "unit_slack": "2.866869392297651739120664114386567863E-5",
  "threshold_slack": "9.215575043584290024453603576020890000E-5",
  "unit_implies_threshold_exact": true,
  "threshold_formula_agrees_reported": true,
  "unit_barrier_holds_midpoint": true,
  "h803_condition_holds_midpoint": true,
  "h801_certified": true
}
```

Minimum unit-barrier slack from tail start:

```json
{
  "n": 126,
  "observed_T_next_over_T_n": "9.999656805356967230530425041972033265E-1",
  "factorial_threshold": "9.998735247852608801527979681614431176E-1",
  "unit_floor_1_minus_1_over_n2": "9.999370118417737465356512975560594608E-1",
  "observed_deficit": "3.431946430327694695749580279667350000E-5",
  "observed_deficit_constant_n2": "5.448558152788248098972033651999884860E-1",
  "threshold_deficit": "1.264752147391198472020318385568824000E-4",
  "threshold_deficit_formula": "1.264752147391198472020318385568823770E-4",
  "threshold_formula_agreement_abs": "2.295150313934596270250806346605279502E-38",
  "threshold_deficit_constant_n2": "2.007920509198266694179457468929064982E+0",
  "unit_slack": "2.866869392297651739120664114386567863E-5",
  "threshold_slack": "9.215575043584290024453603576020890000E-5",
  "unit_implies_threshold_exact": true,
  "threshold_formula_agrees_reported": true,
  "unit_barrier_holds_midpoint": true,
  "h803_condition_holds_midpoint": true,
  "h801_certified": true
}
```

## Tail Constants

Last `20` rows average `n^2(1-T_{n+1}/T_n)` = `0.5528625562566443`.

Last `20` rows average `n^2(1-theta_n)` = `2.0085859329770863`.

## Selected Rows

| n | observed deficit constant | threshold constant | unit slack | unit barrier |
| ---: | ---: | ---: | ---: | --- |
| 2 | `1.752475813076464430147441338290571960E+0` | `2.333333333333333333333333333333333333E+0` | `-1.881189532691161075368603345726429899E-1` | `False` |
| 3 | `1.379477632527589497335609739399206173E+0` | `2.280000000000000000000000000000000000E+0` | `-4.216418139195438859284552659991179699E-2` | `False` |
| 4 | `1.218456208865490355293459692106159891E+0` | `2.224489795918367346938775510204081632E+0` | `-1.365351305409314720584123075663499320E-2` | `False` |
| 5 | `1.124354002253719182739948547524771768E+0` | `2.185185185185185185185185185185185185E+0` | `-4.974160090148767309597941900990870700E-3` | `False` |
| 6 | `1.060554959879069168085726268181070617E+0` | `2.157024793388429752066115702479338844E+0` | `-1.682082218863032446825729671696406022E-3` | `False` |
| 7 | `1.013366291175659086659505084969225591E+0` | `2.136094674556213017751479289940828403E+0` | `-2.727814525644711563164303054943998224E-4` | `False` |
| 8 | `9.764316953207575570425991464923292352E-1` | `2.120000000000000000000000000000000000E+0` | `3.682547606131631712093883360573557000E-4` | `True` |
| 10 | `9.211624953664446780931759201580066900E-1` | `2.096952908587257617728531855955678670E+0` | `7.883750463355532190682407984199331000E-4` | `True` |
| 12 | `8.807678847318382528068497484990103104E-1` | `2.081285444234404536862003780718336486E+0` | `8.280008004733454666190989687568728444E-4` | `True` |
| 16 | `8.237120195587143784549290571966240256E-1` | `2.061394380853277835587929240374609792E+0` | `6.886249235987719591604333703256874000E-4` | `True` |
| 22 | `7.679098169730149260154770845515605056E-1` | `2.044889129259058950784207679826933466E+0` | `4.795251715433575908771134616703295339E-4` | `True` |
| 32 | `7.098051923419337327748734476616592384E-1` | `2.030990173847316704459561602418745242E+0` | `2.833933668535803390870376487679109000E-4` | `True` |
| 48 | `6.541698882954169833933329510376453120E-1` | `2.020720221606648199445983379501385062E+0` | `1.500998748717808231799770177787997778E-4` | `True` |
| 64 | `6.184915642263413213888178612662157312E-1` | `2.015562031124062248124496248992498074E+0` | `9.314170795255338833280813933930280000E-5` | `True` |
| 80 | `5.927233263621991463990793026172006400E-1` | `2.012459950160199359202563189747240960E+0` | `6.363698025590638337514385896606240000E-5` | `True` |
| 100 | `5.684779166403196342675157108228800000E-1` | `2.009974495593545617534910734577410000E+0` | `4.315220833596803657324842891771200000E-5` | `True` |
| 126 | `5.448558152788248098972033651999884860E-1` | `2.007920509198266694179457468929064982E+0` | `2.866869392297651739120664114386567863E-5` | `True` |

## Decision

H804 is a useful simplification: H803 reduces to a finite base plus the cleaner eventual inequality T_{n+1}/T_n >= 1-1/n^2. The finite H803 data supports this from n=8 onward, and the tail deficit constant is far below the factorial threshold constant. The remaining work is analytic: prove the unit-deficit barrier from the Xi kernel or from coefficient asymptotics with explicit error.

## Limitations

- The unit barrier is stronger than H803 and fails for small n, so finite base cases remain necessary.
- The report uses H803 midpoint ratios plus H801 positivity certificates; it is not a new Arb proof.
- No Xi-kernel analytic estimate is proved here.

## Next Target

`H805 saddle-point unit-deficit lemma`: Try to prove eventually T_{n+1}/T_n >= 1-1/n^2, for example from Xi-kernel log-saddle concentration such as T_n = 1 + O(1/(n log n)) with enough regularity to control the quotient.
