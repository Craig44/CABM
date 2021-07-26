#'
#'
#' Validate the multinomial simulator
#'
library(ibm)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
age_comp_expect = extract.run(file.path("Multinomial_Sim", "HG_DS_age.001"))
obs_age = unique(age_comp_expect$HG_DS_age$Values$age)

sim_dir = file.path("Multinomial_Sim", "sim")
sim_file_names = unique(sapply(strsplit(list.files(sim_dir), split = "\\."), "[", 1))
extensions = unique(sapply(strsplit(list.files(sim_dir), split = "\\."), "[", 2))
n_runs = length(extensions)

unique(age_comp_expect$HG_DS_age$Values$year)
year_of_interest = 1982 
expect = age_comp_expect$HG_DS_age$Values[age_comp_expect$HG_DS_age$Values$year == year_of_interest, ]
sim_df = NULL
for(sim in 1:n_runs) {
  ## Fishery AF
  sim_abm = extract.ibm.file(file = file.path(sim_dir, paste0(sim_file_names,".",extensions[sim])), quiet = T)
  abm_years = as.numeric(sim_abm[[1]]$years$value)
  year_ndx = abm_years == year_of_interest
  mat_sim_abm = Reduce(rbind, sim_abm[[1]]$Table$obs)
  class(mat_sim_abm) <- "numeric"
  this_sim = mat_sim_abm[year_ndx, ]
  sim_df = rbind(sim_df, data.frame(age = obs_age, prop = this_sim, sim = sim, year = year_of_interest))
}

abm_plot = ggplot(sim_df, aes(x = age, y = prop, group = sim)) + 
  geom_line(size = 1.4, alpha = 0.4) +
  ggtitle("ABM") +
  geom_line(data = expect, aes(x = age, y = expected), inherit.aes = F, col = "red", size = 1.3)
abm_plot
## now R's version
Rs_sim = rmultinom(n = n_runs, size = 500, prob = expect$expected)
props = sweep(Rs_sim, MARGIN = 2, STATS = colSums(Rs_sim), FUN = "/")
dimnames(props) = list(expect$age, 1:ncol(props))
R_sim_melt = melt(props)
colnames(R_sim_melt) = c("age", "sim", "prop")

R_plot = ggplot(R_sim_melt, aes(x = age, y = prop, group = sim)) + 
  geom_line(size = 1.4, alpha = 0.4) +
  ggtitle("R") +
  geom_line(data = expect, aes(x = age, y = expected), inherit.aes = F, col = "red", size = 1.3)
R_plot
  

grid.arrange(abm_plot, R_plot, ncol = 2, nrow = 1)

my_quantile <- function(x, probs) {
  tibble(x = quantile(x, probs), probs = probs)
}

summarised_R = R_sim_melt %>% 
  group_by(age) %>% 
  summarise(qs = quantile(prop,  probs = c(0.025, 0.5, 0.975)), prob = c(0.025, 0.5, 0.975))

summarised_abm = sim_df %>% 
  group_by(age) %>% 
  summarise(qs = quantile(prop,  probs = c(0.025, 0.5, 0.975)), prob = c(0.025, 0.5, 0.975))

summarised_R$model ="R"
summarised_abm$model ="ABM"

full_df = rbind(summarised_R, summarised_abm)

ggplot(full_df, aes(x = age, y = qs , col = model, linetype = model)) +
  geom_line(size = 1.5) +
  facet_wrap(~prob)


####################################
## Compare simulating from full age-structure vs truncated age-structure
######################################
set.seed(123)
full_age =  expect$expected[1:30]
trunc_age = full_age[5:15]
trunc_age[length(trunc_age)] = sum(full_age[15:(length(full_age))])

Rs_sim_full = rmultinom(n = 1000, size = 500, prob = full_age)
Rs_sim_trunc = rmultinom(n = 1000, size = 500, prob = trunc_age)

## manually truncate full age after simulation
Rs_sim_full_trunc = Rs_sim_full[5:15, ]
Rs_sim_full_trunc[nrow(Rs_sim_full_trunc), ] = colSums(Rs_sim_full[15:nrow(Rs_sim_full), ])
Rs_sim_full_trunc <- sweep(Rs_sim_full_trunc, 2, apply(Rs_sim_full_trunc, 2, sum), "/")
Rs_sim_trunc <- sweep(Rs_sim_trunc, 2, apply(Rs_sim_trunc, 2, sum), "/")

dim(Rs_sim_full_trunc)
dim(Rs_sim_trunc)

trunc_quants = apply(Rs_sim_trunc, 1, FUN = function(x){quantile(x, probs = c(0.025, 0.5, 0.975))})
full_trunc_quants = apply(Rs_sim_full_trunc, 1, FUN = function(x){quantile(x, probs = c(0.025, 0.5, 0.975))})

plot(5:15, trunc_quants["50%",], lty = 2, lwd = 3, type = "o")
lines(5:15, full_trunc_quants["50%",], lty = 2, type = "o", col = "red", lwd = 4)
lines(5:15, trunc_quants["2.5%",], lty = 3, type = "o", col = "black", lwd = 4)
lines(5:15, full_trunc_quants["2.5%",], lty = 3, type = "o", col = "red", lwd = 4)
lines(5:15, trunc_quants["97.5%",], lty = 4, type = "o", col = "black", lwd = 4)
lines(5:15, full_trunc_quants["97.5%",], lty = 4, type = "o", col = "red", lwd = 4)


