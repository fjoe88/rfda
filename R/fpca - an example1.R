library(fda)
library(tidyverse)
library(plotly)

#https://rviews.rstudio.com/2021/05/04/functional-data-analysis-in-r/

# Function to simulate data
set.seed(999)
n_obs <- 80
time_span <- 100
time <- sort(runif(n_obs, 0, time_span))
Wiener <- cumsum(rnorm(n_obs)) / sqrt(n_obs)
y_obs <- Wiener + rnorm(n_obs, 0, .05)


times_basis = seq(0, time_span, 1)
knots    = c(seq(0, time_span, 5)) #Location of knots
n_knots   = length(knots) #Number of knots
n_order   = 4 # order of basis functions: cubic bspline: order = 3 + 1
n_basis   = length(knots) + n_order - 2

basis = create.bspline.basis(c(min(times_basis), max(times_basis)), n_basis, n_order, knots)
n_basis

PHI = eval.basis(time, basis)
dim(PHI)

matplot(
  time,
  PHI,
  type = 'l',
  lwd = 1,
  lty = 1,
  xlab = 'time',
  ylab = 'basis',
  cex.lab = 1,
  cex.axis = 1
)
for (i in 1:n_knots)
{
  abline(v = knots[i], lty = 2, lwd = 1)
}


# Least squares estimate
# estimate basis coefficient
M = ginv(t(PHI) %*% PHI) %*% t(PHI)
c_hat = M %*% Wiener

y_hat = PHI %*% c_hat

df <- data.frame(time = time,
                 Wiener = Wiener,
                 y_hat = y_hat)

# Augment data frame for plotting
p2 <- df %>% ggplot() +
  geom_line(aes(x = time, y = Wiener), col = "grey") +
  geom_point(aes(x = time, y = y_obs)) +
  geom_line(aes(x = time, y = y_hat), col = "red")
p2 + ggtitle("Original curve and least squares estimate") +
  xlab("time") + ylab("f(time)")


# estimate the variance of noise
## SSE = (Y - Xb)'(Y - Xb)
SSE = t(y_hat - y_obs) %*% (y_hat - y_obs)
sigma2 = SSE / (n_obs - n_basis)

# estimate the variance of the fitted curve
# H is the Hat matrix H
# H = X*inv(X'X)*X``
H = PHI %*% M
varYhat = diag(H %*% H * matrix(sigma2, n_obs, n_obs))

# 95% confidence interval

y_hat025 = y_hat - 1.96 * sqrt(varYhat)
y_hat975 = y_hat + 1.96 * sqrt(varYhat)


df <- mutate(df, y_hat025 = y_hat025,
             y_hat975 = y_hat975)
#names(df) <- c("time","Wiener","y_hat", "y_hat025", "y_hat975")
p3 <- df %>% ggplot() +
  geom_line(aes(x = time, y = Wiener), col = "grey") +
  geom_point(aes(x = time, y = y_obs)) +
  geom_line(aes(x = time, y = y_hat), col = "red") +
  geom_line(aes(x = time, y_hat025), col = "green") +
  geom_line(aes(x = time, y_hat975), col = "green")
p3 + ggtitle("Estimated curve with error bars") +
  xlab("time") + ylab("f(time)")


Wiener_obj <-
  smooth.basis(argvals = time,
               y = y_obs,
               fdParobj = basis)


plot(
  time,
  Wiener,
  type = "l",
  xlab = "time",
  ylab = "f(time)",
  main = "Comparison of fda package and naive smoothing estimates",
  col = "grey"
)
lines(time, y_hat, type = "l", col = "red")
lines(Wiener_obj, lwd = 1, col = "blue")
