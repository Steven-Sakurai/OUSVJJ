par = c(0.1, 3,  0.2,  0.1,  -0.6,  3,  0,  -0.6,  -0.03, 0.025, 0.02, 0.01)

converge.list = as.numeric(unlist(read.csv('conv.csv', header = F)))
success.list = as.numeric(unlist(read.csv('success.csv', header = F)))

print("Total trials: ")
print(length(success.list))
print("Successful trails: ")
print(sum(success.list, na.rm = T))
print("Successful Rate: ")
print(sum(success.list, na.rm = T) / length(success.list))

npar = 12
est = as.numeric(unlist(read.csv('est.csv', header = F)))
est = matrix(est, nrow = length(est)/npar, ncol = npar, byrow = T)

est.par.list = data.frame(est)
names(est.par.list) = c("mu", "kappa", "theta", "alpha", "rho", "l0", "l1", "rho_j", "mu_s", "sigma_s", "mu_v", "sigma_v")

print("True: ")
print(par)
print("Mean:")
print(apply(est.par.list, 2, mean))
print("Standard Deviation:")
print(apply(est.par.list, 2, sd))

library(ggplot2)
p1 = ggplot(est.par.list) + geom_histogram(aes(mu))
p2 = ggplot(est.par.list) + geom_histogram(aes(kappa))
p3 = ggplot(est.par.list) + geom_histogram(aes(theta))
p4 = ggplot(est.par.list) + geom_histogram(aes(alpha))
p5 = ggplot(est.par.list) + geom_histogram(aes(rho))
p6 = ggplot(est.par.list) + geom_histogram(aes(l0))
p7 = ggplot(est.par.list) + geom_histogram(aes(l1))
p8 = ggplot(est.par.list) + geom_histogram(aes(rho_j))
p9 = ggplot(est.par.list) + geom_histogram(aes(mu_s))
p10 = ggplot(est.par.list) + geom_histogram(aes(sigma_s))
p11 = ggplot(est.par.list) + geom_histogram(aes(mu_v))
p12 = ggplot(est.par.list) + geom_histogram(aes(sigma_v))

library(grid)
library(gridExtra)
p = grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ncol = 3)
ggsave('distribution.eps', p, device = 'eps', width = 10, height = 8)
