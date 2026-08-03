library(SurrogateRank)

# 1) Simulate a single high-dimensional dataset:
#    200 candidate surrogates, 10% of them truly valid
full_data <- generate.example.data.highdim(
  n1 = 50,  # Number of treated
  n0 = 50, # Number
  p = 200,
  prop_valid = 0.10,
  corr = 0,
  mode = "simple",
  seed = 1
)
