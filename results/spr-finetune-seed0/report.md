# SPR adaptation benchmark

Retrospective evaluation on the previously examined cohort with reportable KD. The target is −log10(KD in nM). Cohort reuse provides no independent replication or basis for candidate promotion.

Brackets give 95% paired sequence group bootstrap intervals on fixed OOF predictions. A positive Spearman difference or negative RMSE difference favors the first model.

## Pooled (158 variants, 115 groups)

| Model or paired difference | Spearman | RMSE |
| --- | ---: | ---: |
| position_ridge | 0.6127 [0.4479, 0.7037] | 0.7406 [0.6267, 0.8205] |
| position_pairwise_ridge | 0.6485 [0.4648, 0.7473] | 0.6869 [0.6024, 0.7690] |
| frozen_head | 0.3181 [0.1830, 0.4677] | 0.9465 [0.7719, 1.0734] |
| final_block | 0.3321 [0.1985, 0.4787] | 0.9431 [0.7688, 1.0699] |
| final_block minus frozen_head | 0.0140 [-0.0021, 0.0371] | -0.0033 [-0.0110, 0.0027] |
| frozen_head minus position_ridge | -0.2947 [-0.4463, -0.0593] | 0.2059 [0.1058, 0.2913] |
| frozen_head minus position_pairwise_ridge | -0.3304 [-0.4937, -0.0633] | 0.2595 [0.1137, 0.3567] |
| final_block minus position_ridge | -0.2807 [-0.4315, -0.0443] | 0.2025 [0.1039, 0.2866] |
| final_block minus position_pairwise_ridge | -0.3164 [-0.4792, -0.0507] | 0.2562 [0.1125, 0.3514] |

## Distant subset (75 variants, 67 groups)

| Model or paired difference | Spearman | RMSE |
| --- | ---: | ---: |
| position_ridge | 0.4248 [0.2031, 0.6078] | 0.6249 [0.5187, 0.7303] |
| position_pairwise_ridge | 0.4718 [0.2668, 0.6355] | 0.5787 [0.4825, 0.6816] |
| frozen_head | 0.2866 [0.0697, 0.4876] | 0.7464 [0.6029, 0.8897] |
| final_block | 0.3062 [0.0873, 0.5066] | 0.7439 [0.6023, 0.8880] |
| final_block minus frozen_head | 0.0196 [-0.0122, 0.0521] | -0.0025 [-0.0142, 0.0083] |
| frozen_head minus position_ridge | -0.1381 [-0.3649, 0.1089] | 0.1215 [-0.0180, 0.2433] |
| frozen_head minus position_pairwise_ridge | -0.1852 [-0.4044, 0.0544] | 0.1677 [0.0343, 0.2792] |
| final_block minus position_ridge | -0.1186 [-0.3415, 0.1250] | 0.1190 [-0.0187, 0.2432] |
| final_block minus position_pairwise_ridge | -0.1656 [-0.3778, 0.0706] | 0.1652 [0.0369, 0.2770] |
