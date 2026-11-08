settings <- list(
  manuscript = "HBOP2FE_Design__Yuexi__v1 (29).pdf",
  alpha = .10, B = 100000L, seed = 20261107L,
  looks = c(20L, 40L, 60L), prior_t = c(.5, .5), prior_c = c(.5, .5),
  grid = expand.grid(lambda = seq(.50, .99, .01), gamma = seq(.10, 1, .10),
                     eta = seq(.50, 3, .10)),
  sensitivity = expand.grid(cap = c(10, 20, 30, 50), quality = c(.50, .75, 1)),
  case_theta = list(lambda = .73, gamma = .50, eta = 2.80)
)

primary_scenarios <- data.frame(
  scenario = paste0("S", 1:6),
  purpose = c("Compatible null", "Compatible alternative", "Mild favorable ECA: null",
              "Strong favorable ECA: null", "Mild unfavorable ECA: alternative",
              "Strong unfavorable ECA: alternative"),
  p_c_target = c(.30, .30, .35, .40, .25, .20),
  p_t = c(.30, .45, .35, .40, .40, .35),
  estimand = c("Type I error", "Power", "Type I error", "Type I error", "Power", "Power")
)
case_scenarios <- data.frame(
  scenario = c("K1", "K2", "K3", "K4", "K5"),
  purpose = c("Compatible null", "Compatible alternative", "Mild favorable ECA: null",
              "Strong favorable ECA: null", "Mild unfavorable ECA: alternative"),
  p_c_target = c(.27, .27, .32, .37, .22), p_t = c(.27, .47, .32, .37, .42),
  estimand = c("Type I error", "Power", "Type I error", "Type I error", "Power")
)
