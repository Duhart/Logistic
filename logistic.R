# Libraries
library(dplyr)
library(rstan)
library(mvtnorm)

# Configure rstan
#options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

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
df$probabilities <- log_reg$fitted.values
boxplot(probabilities~y,data=df)

grid <- expand.grid(list(x1 = seq(0,1,length = 100), x2 = seq(0,1,length = 100)))
surface <-  1/(1+exp(-predict(log_reg,newdata = grid)))

image(seq(0,1,length = 100),seq(0,1,length = 100),matrix(surface,ncol = 100),
      col = rev(heat.colors(15)),
      xlab = 'x1',ylab = 'x2', main = 'Logistic regression', asp = 1, xlim = c(0,1), ylim = c(0,1))
points(df$x1,df$x2,col = unlist(lapply(df$y,function(x){ return(if_else(x==1,'black','blue'))})),
       pch = 15)
legend(1.05, 0.6, c('y = 0', 'y = 1'),col = c('blue','black'),pch = 15, xpd = TRUE, bty = 'n')

contour(x = seq(0,1,length = 100),
        y = seq(0,1,length = 100),
        z = matrix(surface,ncol = 100),
        levels = c(0.1,0.5,0.9), 
        add = TRUE)
slope <- -log_reg$coefficients[2]/log_reg$coefficients[3]

df[df$probabilities <= 0.9 & df$probabilities >= 0.1,] %>% dim()
plot(df$x1[df$probabilities <= 0.9 & df$probabilities >= 0.1], pch =15,
     df$x2[df$probabilities <= 0.9 & df$probabilities >= 0.1],
     col = 1+df$y[df$probabilities <= 0.9 & df$probabilities >= 0.1])


# Logistic regression by hand

data.frame( our_coefficients = beta_hat$beta,
            R_function = log_reg$coefficients)

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
  return(1/(1 + exp(-colSums(apply(X,1,function(k){return(k*beta)})))))
}

# We then define the function F
Fun <- function(beta){
  res <- (X %>% as.matrix() %>% t()) %*% (y - p(beta))
  res <- as.numeric(res)
  return(res)
}


DFun <- function(beta){
  res <- - t(as.matrix(X)) %*% diag(p(beta)*(1-p(beta))) %*% as.matrix(X)
  return(res)
}



# We write a function to do Newton's method
# following the previously described pseudo-code
Newton <- function(beta0, maxiter = 1000, eps = 1e-50){
  # Initialise the results
  current_beta <- beta0
  convergence <- 0
  iter <- 0
  norms <- sqrt(sum(Fun(current_beta) * Fun(current_beta)))
  
  # Start loop to convergence
  while(convergence == 0 & iter < maxiter & norms[length(norms)] > eps){
    # Check if the derivative of F is singular
    if(abs(det(DFun(current_beta))) <= eps){
      # Return current beta with no convergence
      res <- list(beta = current_beta,
                  convergence = 0,
                  iter = iter + 1,
                  norms = norms)
      # Claim it has converged to get out of loop
      convergence <- 1
    }else{
      # We have a non-singular matrix, so we proceed
      # Update beta
      current_beta <- current_beta - solve(DFun(current_beta)) %*% Fun(current_beta)
      # Update norm vector, this is for the history of convergence
      # Recall we should find a quadratic convergence
      norms <- c(norms,sqrt(sum(Fun(current_beta) * Fun(current_beta))))
      # Update number of iterations
      iter <- iter + 1
      # Check if we have converged
      if(sqrt(sum(Fun(current_beta) * Fun(current_beta))) <= eps){
        # Update convergence
        convergence <- 1
        # This is what we want to find
        res <- list(beta = current_beta,
                    convergence = 1,
                    iter = iter,
                    norms = norms)
      }
    }
  }
  # Check if we ran out of iterations
  if(iter >= maxiter){
    # Return what we have with no convergence
    res <- list(beta = current_beta,
                convergence = 0,
                iter = iter,
                norms = norms)
  }
  # Only one return!
  return(res)
}


beta_hat <- Newton(beta, maxiter = 1000, eps = 1e-10)
beta_hat <- Newton(beta)
Fun(beta_hat$beta)
det(DFun(beta_hat$beta))
plot(beta_hat$norms)


beta0 <- c(0,0,0)


beta0 <- beta0 - solve(DFun(beta0)) %*% Fun(beta0)
beta0
Fun(beta0)

beta_hat <- Newton2(beta, maxiter = 3)
plot(beta_hat$norms,type='l')
beta_hat$beta
Fun(beta_hat$beta)
Fun(log_reg$coefficients)
sqrt(sum(Fun(beta)))



plot(df$fitted,p(log_reg$coefficients))
plot(df$fitted,p(c(0,0,0)))
Fun(beta)
Fun(log_reg$coefficients)
det(DFun(beta))
det(DFun(log_reg$coefficients))


plot(df$fitted,col = df$y + 1)
points(p(log_reg$coefficients),type='l')

df[abs(p(log_reg$coefficients) - df$fitted)>0.01,] %>% dim()


plot(1/(1+exp(-apply(cbind(1,df[,1],df[,2]) * log_reg$coefficients,1,sum))),
     1/(1+exp(-apply(X * log_reg$coefficients,1,sum))))

plot(1/(1+exp(-apply(cbind(1,df[,1],df[,2]) * log_reg$coefficients,1,sum))),
     df$fitted)


plot(log_reg$coefficients[1] + X[,2]*log_reg$coefficients[2] + X[,3]*log_reg$coefficients[3],
     log_reg$linear.predictors)

plot(1/(1 + exp(-(X[,1]*log_reg$coefficients[1] + X[,2]*log_reg$coefficients[2] + X[,3]*log_reg$coefficients[3]))),
     log_reg$fitted.values)

plot(1/(1 + exp(-colSums(apply(X,1,function(k){return(k*log_reg$coefficients)})))),
1/(1 + exp(-(X[,1]*log_reg$coefficients[1] + X[,2]*log_reg$coefficients[2] + X[,3]*log_reg$coefficients[3]))))

plot(1/(1 + exp(-colSums(apply(X,1,function(k){return(k*log_reg$coefficients)})))),
     df$fitted)


cor(log_reg$coefficients[1] + X[,2]*log_reg$coefficients[2] + X[,3]*log_reg$coefficients[3],
     log_reg$linear.predictors)

plot(log_reg$coefficients[1] + X[,2]*log_reg$coefficients[2] + X[,3]*log_reg$coefficients[3],
     ylim = c(-15,15))
lines(log_reg$linear.predictors)


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
beta0 <- beta

beta0 <- beta0 - solve(DFun(beta0)) %*% Fun(beta0)
norms <- sum(Fun(beta0) * Fun(beta0))
beta0
norms


df$fitted <- df$probabilities %>% 
  lapply(function(p){ return(if_else(p>=0.6,1,0)) }) %>% 
  unlist()

?xtabs
xtabs(y~fitted,data = df)
xtabs(~y +fitted,data = df)


sqrt(mean((df$y-df$fitted)^2))
sqrt(99/100*var(df$y-df$fitted))

Calif <- data.frame(cut_off = seq(0.01,0.99,0.001))
Calif <- lapply(1:length(Calif$cut_off),function(k){
  aux <- df$probabilities %>%
    lapply(function(p){ return(if_else(p>=Calif$cut_off[k],1,0)) }) %>% unlist()
  cont_tab <- xtabs(~df$y + aux)
  Correct <- cont_tab[1,1] + cont_tab[2,2]
  Jac <- cont_tab[2,2]/(cont_tab[2,2] + cont_tab[1,2] + cont_tab[2,1])
  res <- data.frame(Correct = Correct,
                    Jaccard = Jac)
  return(res)
}) %>% do.call(rbind,.) %>% cbind(Calif,.)
head(Calif)

Calif$cut_off[which(Calif$Correct == max(Calif$Correct))]



plot(Calif$cut_off,Calif$Correct/100,type='l')
lines(Calif$cut_off,Calif$Jaccard,col='red')

b <- c(0.5, 0,0)
### Bayesian
yhat <- function(b){
  aux = as.matrix(cbind(1,df[,1:2])) %*% matrix(b,nrow = 3)
  res <- 1/(1+exp(- aux))
  return(res)
}


ell <- function(b){
  res <- sum(df$y*log(yhat(b)) + (1 - df$y)*log(1-yhat(b)))
  return(res)
}

Lik <- function(b){
  return(exp(ell(b)))
}

ker <- 

yhat(b)
yhat(blr$coefficients)

ell(b)
ell(blr$coefficients)

Lik(b)
Lik(blr$coefficients)

##### MCMC

# Step 1. Inintialise b0
b <- blr$coefficients
b_path <- b
prior <- function(b){
  return(dmvnorm(b,rep(0,3),diag(rep(50,3))))
}

# Here comes the cycle
burn <- 5000
sample_size <- 1000
for(k in 1:(burn + sample_size)){
  # Step 2. Sample a new value b_new from distribution q(b)
  b_prop <- rmnorm(n = 1, mean = b, varcov = diag(rep(2,3)))
  
  # Step 3. Compute the acceptance radio
  R <- Lik(b_prop)*prior(b_prop) / (Lik(b)*prior(b))
  
  # Step 4. Select new value of b via the Markov Chanin
  if(runif(1) < R){
    b <- b_prop
  }
  b_path <- rbind(b_path,b)
}

par(mfrow=c(3,2))
plot(b_path[nrow(b_path) - c(999:0),1],type='l') 
hist(b_path[nrow(b_path) - c(999:0),1]) 

plot(b_path[nrow(b_path) - c(999:0),2],type='l') 
hist(b_path[nrow(b_path) - c(999:0),2]) 

plot(b_path[nrow(b_path) - c(999:0),3],type='l') 
hist(b_path[nrow(b_path) - c(999:0),3]) 

apply(b_path,2,mean)
blr$coefficients


  
blr <- glm(y ~ x1 + x2, data = df, family = binomial(link = 'logit'))
blr
blr %>% summary() %>% .$cov.unscaled %>% diag() %>% mean()
aux$c
## STAN

# data
logistic_dat <- list(N = 100,
                     y = df$y, x1 = df$x1, x2 = df$x2,
                     b0_init = blr$coefficients[1],
                     b1_init = blr$coefficients[2],
                     b2_init = blr$coefficients[3],
                     sig = 20)

# run model
MCMC_log_reg <- stan(file = 'bayes_log_reg.stan', data = logistic_dat)
log_reg_stats <- print(MCMC_log_reg)



plot(MCMC_log_reg,pars = c('b0','b1','b2'))
stan_hist(MCMC_log_reg, pars = c('b0','b1','b2'))
stan_trace(MCMC_log_reg, pars = 'b0')
stan_plot(MCMC_log_reg)

plot(density(traces$data$value[traces$data$parameter == 'b0']))



par(mfrow = c(1,3))

dens <- density(traces$data$value[traces$data$parameter == 'b0'])
plot(dens$x,dens$y,col='white', xlab = 'Constant',ylab = 'Density')
polygon(c(dens$x,rev(dens$x)),c(rep(0,length(dens$x)),rev(dens$y)),
        col = rgb(1,0,0,0.3))

dens <- density(traces$data$value[traces$data$parameter == 'b1'])
plot(dens$x,dens$y,col='white', xlab = 'x1',ylab = 'Density',
     main = 'Bayesian Logistic Regression')
polygon(c(dens$x,rev(dens$x)),c(rep(0,length(dens$x)),rev(dens$y)),
        col = rgb(0,1,0,0.3))

dens <- density(traces$data$value[traces$data$parameter == 'b2'])
plot(dens$x,dens$y,col='white', xlab = 'x2',ylab = 'Density')
polygon(c(dens$x,rev(dens$x)),c(rep(0,length(dens$x)),rev(dens$y)),
        col = rgb(0,0,1,0.3))



?polygon

stan_dens(MCMC_log_reg,pars = c('b0')) +
  ggtitle('Constant') +
  xlab('Constant') +
  ylab('Density')


stan_scat(MCMC_log_reg,pars = c('b2','b0'))
stan_ac(MCMC_log_reg,pars = c('b0','b1','b2'))

aux2 <- stan_trace(MCMC_log_reg, pars = c('b0','b1','b2'))
hist(aux2$data$value[aux2$data$parameter == 'b0'])
mean(aux2$data$value[aux2$data$parameter == 'b0'])
print(MCMC_log_reg)

str(aux2$data)

print(MCMC_log_reg)

grid <- expand.grid(list(x1 = seq(0,1,length = 100), x2 = seq(0,1,length = 100)))

blr$coefficients
aux2 <- summary(MCMC_log_reg)
str(aux2$summary)
dimnames(aux2$summary)

traces <- stan_trace(MCMC_log_reg, pars = c('b0','b1','b2'))
bayes_coeff <- c(mean(traces$data$value[traces$data$parameter == 'b0']),
                 mean(traces$data$value[traces$data$parameter == 'b1']),
                 mean(traces$data$value[traces$data$parameter == 'b2']))
data.frame( our_coefficients = beta_hat$beta,
            R_function = log_reg$coefficients,
            Bayesian_regression = bayes_coeff)


grid$surface <-  1/(1+exp(-bayes_coeff[1] - bayes_coeff[2] * grid$x1 - bayes_coeff[3] * grid$x2))


par(mfrow=c(1,1))
image(seq(0,1,length = 100),seq(0,1,length = 100),matrix(grid$surface,ncol = 100),
      col = rev(heat.colors(15)),
      xlab = 'x1',ylab = 'x2', main = 'Bayesian Logistic regression', asp = 1,
      xlim = c(0,1), ylim = c(0,1))


points(df$x1,df$x2,col = unlist(lapply(df$y,function(x){ return(if_else(x==1,'black','blue'))})),
       pch = 15)

legend(1.05, 0.6, c('y = 0', 'y = 1'),col = c('blue','black'),pch = 15, xpd = TRUE, bty = 'n')

contour(x = seq(0,1,length = 100),
        y = seq(0,1,length = 100),
        z = matrix(grid$surface,ncol = 100),
        levels = c(0.1,0.5,0.9), 
        add = TRUE)
