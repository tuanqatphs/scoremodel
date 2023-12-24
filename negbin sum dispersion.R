


mu1  <- 2.3
mu2 <- 4.1
dispersion <- 2 # Common dispersion parameter for x1 and x2.

# Simulate two independent negbin variables.
# Notie the parameterization of the dispersion parameter.
x1 <- rnbinom(n = 100000, mu = mu1, size = 1/dispersion)
x2 <- rnbinom(n = 100000, mu = mu2, size = 1/dispersion)

# Their sum.
x3 <- x1 + x2

mean(x1)
mean(x2)
mean(x3) 

# Compare the variances with the theoretical population variances
var1 <- mu1 + (mu1^2)*dispersion

var1
var(x1)

# for x2 as well.
var2 <- mu2 + (mu2^2)*dispersion

var2
var(x2)


# For x3 = x1 + x2.

# Finding the new dispersion parameter.
# var(x3) = var(x1) + var(x2)  (because of independence).

# Since sum of two indep negbin variables is also negbin, with mu3 = mu1+mu2 and a dispersion:
# var(x3) = mu1 + (mu1^2)*dispersion + mu2 + (mu2^2)*dispersion
# mu3 + (mu3^2)*dispersion3 = mu1 + (mu1^2)*dispersion + mu2 + (mu2^2)*dispersion

# solving for dispersion3 gives:
# dispersion3 = ((dispersion*(mu1^2)) + (dispersion*(mu2^2)))  / (mu3^2) 

# Check the theoretical vs simulated.

mu3 <- mu1 + mu2
dispersion3 = ((dispersion*(mu1^2)) + (dispersion*(mu2^2)))  / (mu3^2) 

var3 <- mu3 + (mu3^2)*dispersion3

var3
var(x3)



