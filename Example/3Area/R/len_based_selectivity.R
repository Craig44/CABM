# length based selectivity in an age based model
##  TODO: To complete this code you might want to illustrate deistributions assumed around mean length

library(ggplot2)
library(cyrils)
ages = 1:15
k = 0.24
t0 = -0.3
l_inf = 60
cv = 0.2
n_quants = 10
a50 = 43
ato95 = 12
quants = (1:n_quants-0.5)/n_quants

length_at_age = VB(ages, k, l_inf,t0)
plot(ages, length_at_age)

## compare qnorm from standard normal (probit()) with qnorm with a specified mean and sigma
mean = length_at_age[4]
sigma = mean* cv
## prefect we will use the probit with stan
mean + qnorm(p = quants) * sigma
lengths = qnorm(p = quants, mean = mean, sd = sigma)

sum(logis(lengths, a50,ato95 )) / n_quants
age_based_select = vector()
for(i in ages) {
  lengths = qnorm(p = quants, mean = length_at_age[i], sd = length_at_age[i] * cv)
  age_based_select[i] = sum(logis(lengths, a50,ato95 )) / n_quants
}
age_based_select_5 = vector()
quant_5 = (1:5-0.5)/5
for(i in ages) {
  lengths = qnorm(p = quant_5, mean = length_at_age[i], sd = length_at_age[i] * cv)
  age_based_select_5[i] = sum(logis(lengths, a50,ato95 )) / 5
}
age_based_select_1 = vector()
quant_1 = 0.5
for(i in ages) {
  lengths = qnorm(p = quant_1, mean = length_at_age[i], sd = length_at_age[i] * cv)
  age_based_select_1[i] = sum(logis(lengths, a50,ato95 )) / 1
}
jpeg(filename = "../Figures/length_based_sel.jpeg", units = "in", res = 300, width = 14, height = 6)
par(mfrow = c(1,3))
plot(ages, length_at_age, type = "l", lwd = 2, col = "black", xlab = "Age", ylab = "Length", main = "age length relationship")
plot(1:100, logis(1:100, a50,ato95 ), type = "l", lwd = 2, col = "black", xlab = "Length", ylab = "Selectivity", main = "Length based selectivty", ylim = c(0,1))
plot(ages, age_based_select, type = "l", lwd = 2, col = "black", xlab = "Age", ylab = "", main = "Length based selectivty for age", ylim = c(0,1))
lines(ages, age_based_select_1, type = "l", lwd = 2, col = "red", lty = 2)
lines(ages, age_based_select_5, type = "l", lwd = 2, col = "blue", lty = 2)

legend('bottomright', legend = c("n_quant = 1","n_quant = 5", "n_quant = 10"), lty = c(2,2,1), col = c("red","blue","black"), lwd = 2)
dev.off()




