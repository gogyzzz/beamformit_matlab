# A MATLAB implementation of CHiME4 baseline Beamformit


## References

- "Acoustic beamforming for speaker diarization of meetings", Xavier Anguera, Chuck Wooters and Javier Hernando, IEEE Transactions on Audio, Speech and Language Processing, September 2007, volume 15, number 7, pp.2011-2023.
- [official beamformit github](https://github.com/xanguera/BeamformIt)

## Requirements

| script | requirement |
|---|---|
| beamformit.m | MATLAB supporting audioread (also you can use on OCTAVE by installing signal package) |
| beamformit_step_by_step.mlx | MATLAB supporting mlx format |

## Implementation detail
See beamformit_step_by_step.{html,pdf} (and beamformit_step_by_step.mlx)

Not all parts are the same. 
The reference microphone acquisition and viterbi decoding are different from the original.

## How to run

See beamformit.m 

## Result

my beamformit is slightly worse, especially in real_et05 task.

```sh
## original version

local/chime4_calc_wers.sh exp/tri3b_tr05_multi_noisy beamformit_5mics exp/tri3b_tr05_multi_noisy/graph_tgpr_5k

-------------------
best overall dt05 WER 13.66% (language model weight = 11)
-------------------
dt05_simu WER: 14.34% (Average), 12.82% (BUS), 17.09% (CAFE), 11.90% (PEDESTRIAN), 15.56% (STREET)
-------------------
dt05_real WER: 12.98% (Average), 15.96% (BUS), 12.67% (CAFE), 10.02% (PEDESTRIAN), 13.26% (STREET)
-------------------
et05_simu WER: 21.33% (Average), 15.75% (BUS), 22.97% (CAFE), 22.54% (PEDESTRIAN), 24.06% (STREET)
-------------------
et05_real WER: 21.80% (Average), 30.08% (BUS), 20.62% (CAFE), 19.90% (PEDESTRIAN), 16.62% (STREET)
-------------------

## my version

local/chime4_calc_wers.sh exp/tri3b_tr05_multi_noisy bfit5ch_1008 exp/tri3b_tr05_multi_noisy/graph_tgpr_5k

-------------------
best overall dt05 WER 14.16% (language model weight = 11)
-------------------
dt05_simu WER: 14.70% (Average), 13.05% (BUS), 17.12% (CAFE), 12.39% (PEDESTRIAN), 16.24% (STREET)
-------------------
dt05_real WER: 13.61% (Average), 16.82% (BUS), 13.35% (CAFE), 11.03% (PEDESTRIAN), 13.24% (STREET)
-------------------
et05_simu WER: 22.19% (Average), 16.51% (BUS), 23.66% (CAFE), 23.81% (PEDESTRIAN), 24.77% (STREET)
-------------------
et05_real WER: 22.83% (Average), 31.30% (BUS), 21.67% (CAFE), 21.49% (PEDESTRIAN), 16.88% (STREET)
-------------------

```
