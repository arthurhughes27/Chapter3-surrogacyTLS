library(SurrogateRank)

fig_dir <- fs::path("output", "figures", "tutorial")
fs::dir_create(fig_dir)

# Simulate multi-study, high-dimensional individual participant data:
# 5 studies, 25 treated / 25 untreated per study, 100 candidate
# surrogates, 10% of which are truly valid
ipd_data <- generate.example.data.highdim.multistudy.ipd(
  M = 5,             # number of studies
  n1 = 25,           # treated individuals per study
  n0 = 25,           # untreated individuals per study
  p = 100,           # number of candidate surrogate markers
  prop_valid = 0.10  # proportion of truly valid surrogates
)

# Screening stage: RISE is applied within each study individually,
# then the resulting study-level effects are combined via
# random-effects meta-analysis to identify markers with consistent
# evidence of surrogacy across studies
screen_meta_result <- rise.screen.meta(
  yone      = ipd_data$y1,
  yzero     = ipd_data$y0,
  sone      = ipd_data$s1,
  szero     = ipd_data$s0,
  studyone  = ipd_data$study1,
  studyzero = ipd_data$study0,
  epsilon.study = 0.2, # non-inferiority margin for within-study screening
  epsilon.meta  = 0.2  # non-inferiority margin for the meta-analysis stage
)

sig_markers    <- screen_meta_result$significant.markers
screen_weights <- screen_meta_result$screening.weights

study_level_res = screen_meta_result$screening.metrics.study %>%
  arrange(marker) %>%
  head(n = 5)

study_level_res

meta_level_res = screen_meta_result$screening.metrics.meta %>%
  arrange(p.adjusted) %>%
  head(n = 5)

meta_level_res  %>% as.data.frame()

cat(length(sig_markers), "markers retained after meta-analytic screening\n")

p1 = screen_meta_result$gamma.s.plot$screen.plot

p1

ggsave(plot = p1, path = fig_dir,
       filename = "tutorial_screen_meta.pdf",
       height = 15, width = 30, units = "cm")

# Evaluation stage: combine the retained markers into a single
# composite surrogate signature and re-evaluate its surrogacy across
# studies via meta-analysis
eval_meta_result <- rise.evaluate.meta(
  yone      = ipd_data$y1,
  yzero     = ipd_data$y0,
  sone      = ipd_data$s1,
  szero     = ipd_data$s0,
  studyone  = ipd_data$study1,
  studyzero = ipd_data$study0,
  epsilon.study = 0.2,
  epsilon.meta  = 0.2,
  markers = screen_meta_result$significant.markers,
  screening.weights = screen_meta_result$screening.weights
)

# Meta-analytic evaluation results for the composite marker
print(eval_meta_result$evaluation.metrics.meta %>% as.data.frame() )

p2 = eval_meta_result$gamma.s.plot$forest.plot

p2

ggsave(plot = p2, path = fig_dir,
       filename = "tutorial_evaluate_meta.pdf",
       height = 15, width = 30, units = "cm")
