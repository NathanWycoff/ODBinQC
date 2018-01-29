
rho <- mean(X/n)
sig_pi <- sqrt(rho * (1 - rho) / n)
sd((X / n - rho) / sig_pi )

sqrt(bb.var(alpha, beta, n)) / (n * sig_pi)

x <- rnorm(1e5, 0, 2.3)
sd(x)
mean(abs(diff(x))) / 1.128
