mkdir small_ecori_msei_1M_error_output

python3 ../readsynth/readsynth.py \
  -g abundances.csv \
  -o small_ecori_msei_1M_error_output/ \
  -n 1_000_000 \
  -u 400 \
  -sd 100 \
  -l 150 \
  -m1 ecori \
  -m2 msei \
  -q1 ../../resources/q_scores/resources/q_scores/novaseq_6000_SP_251bp_R1_sampled_scores.csv \
  -q2 ../../resources/q_scores/resources/q_scores/novaseq_6000_SP_251bp_R2_sampled_scores.csv
