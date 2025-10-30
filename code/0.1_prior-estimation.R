# Prior development and estimation
### adapted from code by: Dr. Thomas Reicke, University of Montana

# Simulate priors 
res <- 100
x <- seq(-3,3, length.out = res)

n <- 1000
beta0 <- rlogis(n, 0, 1)      # y-intercept
beta1 <- rnorm(n, -0.25, 0.5) # weakly informative prior


# Check the histogram of the y-intercept (should be relatively uniform, as opposed to beta1)
hist(plogis(beta0))
hist(beta1)

# For loop to see how this would interact in the model
p <- matrix(NA, n, res)

for (j in 1:res){
  p[,j] <- beta0 + beta1 * x[j]
}

boxplot(plogis(p), names = round(x,digits = 1), outline = F)
# as you can see by the boxplot, this distribution would result in a prior with a negative
# slope that allows for all possible values, but is more concentrated between 0.3 and 0.7







