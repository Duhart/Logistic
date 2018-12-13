# Libraries
library(dplyr)


# Data
df <- read.csv(paste0(getwd(),'/data/data.csv'))
head(df)
summary(df)
par(pty = 's')
plot(df$x1,df$x2,col=df$y + 1,pch = 15,
     xlab = 'x1',ylab = 'x2', main = 'Fake data', asp = 1, xlim = c(0,1), ylim = c(0,1))
legend(1.05, 0.6, c('y = 0', 'y = 1'),col = c(1,2),pch = 15, xpd = TRUE, bty = 'n')


df[c(which(df$x1 == min(df$x1)),which(df$x1 == max(df$x1))),]
# Logistic regression
log_reg <- glm(y ~ x1 + x2,family = 'binomial', data = df)
summary(log_reg)
df$fitted <- log_reg$fitted.values
boxplot(fitted~y,data=df)

grid <- expand.grid(list(x1 = seq(0,1,length = 100), x2 = seq(0,1,length = 100)))
surface <-  1/(1+exp(-predict(log_reg,newdata = grid)))

image(seq(0,1,length = 100),seq(0,1,length = 100),matrix(surface,ncol = 100),
      xlab = 'x1',ylab = 'x2', main = 'Logistic regression', asp = 1, xlim = c(0,1), ylim = c(0,1))
points(df$x1,df$x2,col = df$y + 1, pch = 15)
contour(x = seq(0,1,length = 100),
        y = seq(0,1,length = 100),
        z = matrix(surface,ncol = 100),
        levels = 0.5, 
        add = TRUE)
slope <- -log_reg$coefficients[2]/log_reg$coefficients[3]


# Logistic regression by hand

# We first set the variables as our equations.
y <- df$y

# We will add a column of ones to our covariate matrix just to have a constant in our model.
# Our equations never really prohibited this. The important issue will be that the inverse
# of the derivative of F exists.
X <- data.frame(Constant = rep(1,length(y)),
                x1 = df$x1,
                x2 = df$x2)

# Since the MLE of p is 0.5, we will start with a beta that does that.
# This means all 3 betas in 0, this way logit(p)=0 and hence p=0.5 for all observations.
beta <- rep(0,length=3)

# Define how p depends on X and beta.
# Note that we will consider X to be fixed and we will be interested in p
# just as a function of beta.
p <- function(beta){
  return(1/(1 + exp(-apply(X * beta,1,sum))))
}

# We then define the function F
Fun <- function(beta){
  res <- (X %>% as.matrix() %>% t()) %*% (y - p(beta))
  res <- as.numeric(res)
  return(res)
}


DFun <- function(beta){
  res <- - t(as.matrix(X)) %*% diag(p(beta)) %*% as.matrix(X)
  return(res)
}



# We write a function to do Newton's method
# following the previously described pseudo-code
Newton <- function(beta0, maxiter = 1000, eps = 1e-8){
  current_beta <- beta0
  convergence <- 0
  iter <- 0
  norms <- list(sum(Fun(current_beta) * Fun(current_beta)))
  if(sum(Fun(current_beta) * Fun(current_beta)) <= eps){
    res <- list(current_beta,
                convergence = 1,
                iter = 0,
                norms = unlist(norms))
    convergence <- 1
  }
  while(convergence == 0 & iter < maxiter){
    if(abs(det(DFun(current_beta))) <= eps){
      norms[[iter + 1]] <- sum(Fun(current_beta) * Fun(current_beta))
      res <- list(beta = current_beta,
                  convergence = 0,
                  iter = iter + 1,
                  norms = unlist(norms))
      convergence <- 1
    }else{
      current_beta <- current_beta - solve(DFun(current_beta)) %*% Fun(current_beta)
      norms[[iter + 1]] <- sum(Fun(current_beta) * Fun(current_beta))
      iter <- iter + 1
      if(sum(Fun(current_beta) * Fun(current_beta)) <= eps){
        res <- list(beta = current_beta,
                    convergence = 1,
                    iter = iter,
                    norms = unlist(norms))
        convergence <- 1
      }
    }
  }
  if(iter >= maxiter){
    norms[[maxiter]] <- sum(Fun(current_beta) * Fun(current_beta))
    res <- list(beta = current_beta,
                convergence = 0,
                iter = maxiter,
                norms = unlist(norm))
  }
  return(res)
}


Newton2 <- function(beta,maxiter = 1000,eps = 1e-5){
  norms <- list()
  for(k in 1:maxiter){
    beta <- beta - solve(DFun(beta)) %*% Fun(current_beta)
    norms[[k+1]] <- sum(Fun(beta) * Fun(beta))
  }
  return(list(beta = beta, norms = unlist(norms)))
}

beta_hat <- Newton(beta, maxiter = 100)
beta_hat <- Newton2(beta, maxiter = 4)

dim(X)
dim(Dp(beta))
DFun(beta)
Fun(beta)

Fun(log_reg$coefficients)
DFun(log_reg$coefficients)

dp <- Dp(beta)
dp
dim(dp)
X %>% as.matrix() %>% dim()
- dp %*% as.matrix(X)

beta0 <- log_reg$coefficients


beta0 <- beta0 - solve(DFun(beta0)) %*% Fun(beta0)
norms <- sum(Fun(beta0) * Fun(beta0))
beta0
norms
?fsolve

log_reg$coefficients


aux <- pracma::fsolve(Fun,beta,DF)
