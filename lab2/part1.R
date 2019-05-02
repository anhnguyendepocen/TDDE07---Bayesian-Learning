library(mvtnorm)
library(readr)
library(matlib)

TempLinkoping <- read_delim("~/Desktop/LIU VT 2019/TDDE07/TDDE07---Bayesian-Learning/lab2/TempLinkoping.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
#View(TempLinkoping)

x = TempLinkoping[1]
y = TempLinkoping[2]


mu0 = t(c(-10,100,-100))
omega0 = c(0.01, 0.01, 0.01)
v0 = 4
sigmasq0 = 1

e = rnorm(1, mean = 0, sd = sigmasq0)

plot(TempLinkoping)

betahat = inv((t(x)*x))*(t(x)*y) 