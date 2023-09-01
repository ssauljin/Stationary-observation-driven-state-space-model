library(MASS)
# load rawdata from URL
#load(url("https://sites.google.com/a/wisc.edu/jed-frees/home/documents/data.RData"))
#load(url("https://sites.google.com/a/wisc.edu/jed-frees/home/documents/dataout.RData"))

# train <- data[,c(1:2,4,9:14,21:22,24:25)]   # Only IM claim is used
# test <- dataout[,c(1:2,4,9:14,21:22,24:25)] # Only IM claim is used
load("data.Rdata")
#test <- test[-2,]

tClaim <- test$ClaimIM
tCount <- test$FreqIM
test$FreqIM <- 0
trainp <- subset(train, yAvgIM >0)

freqglmIM <- glm.nb(FreqIM ~ TypeCity+TypeCounty+TypeMisc+TypeSchool+TypeTown+CoverageIM+lnDeductIM,
                 data=train)
lambda1 <- exp(predict(freqglmIM))

Rprox <- as.data.frame(cbind(train$Year, train$PolicyNum, train$FreqIM, lambda1))
colnames(Rprox) <- c("Year", "PolicyNum", "FreqIM", "lambda1")

Rprox <- Rprox[order(Rprox$PolicyNum, Rprox$Year),]

dzdz <- aggregate(Year~PolicyNum,Rprox,min)
colnames(dzdz)[2] <- "minYear"


zdzd <- aggregate(Year~PolicyNum,Rprox,max)
colnames(zdzd)[2] <- "maxYear"

Rprox <- merge(Rprox, dzdz, by = "PolicyNum", all.x=TRUE)
Rprox <- merge(Rprox, zdzd, by = "PolicyNum", all.x=TRUE)

Rprox$Index <- Rprox$Year - Rprox$minYear + 1
rm(dzdz, zdzd)

h1 <- function(n, alpha, beta, lambda) {
  result <- dnbinom(n, size=alpha, prob=beta/(lambda+beta))
#  mean = alpha/beta * lambda,
#  variance = alpha/beta^2 * lambda^2 + alpha/beta * lambda = mean + mean^2/alpha
  return(result) }

joint_loglik <- function(p, q, alpha_0) {
  MM <- 1
  n <- Rprox$FreqIM
  lambda1 <- Rprox$lambda1
  I       <- length(unique(Rprox$PolicyNum))
  # alpha1_0 <- 3.8
  # alpha2_0 <- 12
  # beta1_0  <- 3.8
  # beta2_0  <- 11
  # psi      <- summary(sevglmIM)$dispersion

  alpha1_0 <- alpha_0
  beta1_0  <- alpha_0

  loglik <- 0
  
  for (i in 1:I) {
  
  beta1   <- beta1_0    
  alpha1  <- alpha1_0    

  loglik  <- loglik + log(
             h1(n[MM], alpha1, beta1 , lambda1[MM]) )
  
  alpha1  <- alpha1 + n[MM] 
  beta1   <- beta1  + lambda1[MM]
  
  MM      <- MM+1
  
      while (Rprox$Index[MM %% (nrow(Rprox)+1) + MM %/% (nrow(Rprox)+1)] >1) {
        # MM %% (nrow(Rprox)+1) + MM %/% (nrow(Rprox)+1) = MM other than MM=I
        alpha1  <- q*p  * alpha1  + beta1*q*(1-p)
        beta1   <- q    * beta1 
        loglik  <- loglik + log(
                 h1(n[MM], alpha1, beta1 , lambda1[MM]) )
        
        alpha1  <- alpha1 + n[MM] 
        beta1   <- beta1  + lambda1[MM]
        
        MM      <- MM+1
      } 
  }
  return(-loglik)
  #return(list(neglik=-loglik, M=MM)) 
  }

static_loglik <- function(parm) {
  result <- joint_loglik(p=1, q=1, alpha_0=parm[1]) 
  return(result)  }

system.time(
static_fit <- optim(freqglmIM$theta, static_loglik, method="Brent", lower=0.1, upper=100)
)

incvar_loglik <- function(parm) {
  result <- joint_loglik(p=1, q=parm[1], alpha_0=parm[2]) 
  return(result)  }

system.time(
  incvar_fit <- optim(c(0.99,freqglmIM$theta), incvar_loglik,
                      method="L-BFGS-B", upper=c(0.9999999999, Inf ))
)

decvar_loglik <- function(parm) {
  result <- joint_loglik(p=parm[1], q=1, alpha_0=parm[2]) 
  return(result)  }

system.time(
  decvar_fit <- optim(c(0.99,freqglmIM$theta), decvar_loglik,
                      method="L-BFGS-B", upper=c(0.9999999999, Inf ))
)

bddvar_loglik <- function(parm) {
  result <- joint_loglik(p=parm[1], q=parm[2], alpha_0=parm[3]) 
  return(result)  }


system.time(
  bddvar_fit <- optim(c(0.99,0.99,freqglmIM$theta), bddvar_loglik,
                      method="L-BFGS-B", upper=c(0.9999999999, 0.9999999999, Inf ))
)

contvar_loglik <- function(parm) {
  
  MM <- 1
  n <- Rprox$FreqIM
  lambda1 <- Rprox$lambda1
  I       <- length(unique(Rprox$PolicyNum))
  # alpha1_0 <- 3.8
  # alpha2_0 <- 12
  # beta1_0  <- 3.8
  # beta2_0  <- 11
  # psi      <- summary(sevglmIM)$dispersion
  
  p       <- parm[1]
  alpha_0 <- parm[2]
  
  alpha1_0 <- alpha_0
  beta1_0  <- alpha_0
  
  loglik <- 0
  
  for (i in 1:I) {
    
    beta1   <- beta1_0    
    alpha1  <- alpha1_0    
    
    loglik  <- loglik + log(
      h1(n[MM], alpha1, beta1 , lambda1[MM]) )
    
    alpha1  <- alpha1  + n[MM] 
    beta1   <- beta1   + lambda1[MM] 
 #   q       <- alpha_0 / (alpha_0*p^2+(1-p^2)*beta1)
    
    MM      <- MM+1
    
    while (Rprox$Index[MM %% (nrow(Rprox)+1) + MM %/% (nrow(Rprox)+1)] >1) {
      alpha1  <- p *(alpha_0 * alpha1) / (alpha_0*p^2+(1-p^2)*beta1) +
              (1-p)*(alpha_0 *  beta1) / (alpha_0*p^2+(1-p^2)*beta1)
      beta1   <-    (alpha_0 *  beta1) / (alpha_0*p^2+(1-p^2)*beta1)
      

      loglik  <- loglik + log(
        h1(n[MM], alpha1, beta1 , lambda1[MM]) )
      
      alpha1  <- alpha1  + n[MM] 
      beta1   <- beta1   + lambda1[MM] 

      MM      <- MM+1
    } 
  }
  return(-loglik)
  #return(list(neglik=-loglik, M=MM)) 
}


system.time(
  contvar_fit <- optim(c(1,freqglmIM$theta), contvar_loglik)
)

reparmtable <- cbind(
c(  freqglmIM$theta ,  0,                1     ), # independent
c( static_fit$par[1],  1,                1     ), # static
c( incvar_fit$par[2],  1, incvar_fit$par[1]    ), # increasing
c( decvar_fit$par[2],     decvar_fit$par[1], 1 ), # decreasing
c(contvar_fit$par[2],    contvar_fit$par[1], NA)) # constant


library(numDeriv)

static_hess  <- hessian(static_loglik, static_fit$par)
diag(sqrt(solve(static_hess)))
incvar_hess  <- hessian(incvar_loglik, incvar_fit$par)
diag(sqrt(solve(incvar_hess)))
decvar_hess  <- hessian(decvar_loglik, decvar_fit$par)
diag(sqrt(solve(decvar_hess)))
contvar_hess <- hessian(contvar_loglik, contvar_fit$par)

rownames(reparmtable) <- c("$\\alpha_0$", "$p$", "$q$")

goftable <- cbind(
c( logLik(freqglmIM), AIC(freqglmIM), BIC(freqglmIM)),  
c( -static_fit$value,
  2*static_fit$value + 2*9,
  2*static_fit$value + log(nrow(train))*9),
c( -incvar_fit$value,
  2*incvar_fit$value + 2*10,
  2*incvar_fit$value + log(nrow(train))*10),
c( -decvar_fit$value,
  2*decvar_fit$value + 2*10,
  2*decvar_fit$value + log(nrow(train))*10),
c(-contvar_fit$value,
 2*contvar_fit$value + 2*10,
 2*contvar_fit$value + log(nrow(train))*10))

rownames(goftable) <- c("Loglik", "AIC", "BIC")


contvar_parmmat <- function(parm) {
  MM <- 1
  n <- Rprox$FreqIM
  lambda1 <- Rprox$lambda1
  I       <- length(unique(Rprox$PolicyNum))

  alpha_0 <- parm[1]
  p       <- parm[2]
 
  alpha1_0 <- alpha_0
  beta1_0  <- alpha_0
  
  loglik <- 0
  
  a1       <- rep(NA, nrow(Rprox))
  b1       <- rep(NA, nrow(Rprox))
  ps       <- rep(NA, nrow(Rprox))

    for (i in 1:I) {
    
    beta1   <- beta1_0  + lambda1[MM] 
    alpha1  <- alpha1_0 +       n[MM]   
  #  q       <- alpha_0 / (alpha_0*p^2+(1-p^2)*beta1)

    a1[ MM] <- alpha1 
    b1[ MM] <- beta1  
    ps[ MM] <- p

    MM      <- MM+1
    
    while (Rprox$Index[MM %% (nrow(Rprox)+1) + MM %/% (nrow(Rprox)+1)] >1) {
      alpha1  <- p *(alpha_0 * alpha1) / (alpha_0*p^2+(1-p^2)*beta1) +
              (1-p)*(alpha_0 *  beta1) / (alpha_0*p^2+(1-p^2)*beta1)
      beta1   <-    (alpha_0 *  beta1) / (alpha_0*p^2+(1-p^2)*beta1)
      
      alpha1  <- alpha1  + n[MM] 
      beta1   <- beta1   + lambda1[MM] 

      a1[ MM] <- alpha1
      b1[ MM] <- beta1
      ps[ MM] <- p
      MM <- MM+1
      
    } 
    }
  mat <- cbind(a1, b1, ps)
  return(mat) }

joint_parmmat <- function(parm) {
  MM <- 1
  n <- Rprox$FreqIM
  lambda1 <- Rprox$lambda1
  I       <- length(unique(Rprox$PolicyNum))
  
  p        <- parm[2]
  q        <- parm[3]
  alpha_0  <- parm[1]
  alpha1_0 <- alpha_0
  beta1_0  <- alpha_0
  
  loglik <- 0
  
  
  a1       <- rep(NA, nrow(Rprox))
  b1       <- rep(NA, nrow(Rprox))
  ps       <- rep(NA, nrow(Rprox))
  qs       <- rep(NA, nrow(Rprox))
  
  for (i in 1:I) {
    
    beta1   <- beta1_0  + lambda1[MM]   
    alpha1  <- alpha1_0 +       n[MM]    
    
    a1[ MM] <- alpha1
    b1[ MM] <- beta1
    ps[ MM] <- p
    qs[ MM] <- q
    
    MM      <- MM+1
    
    while (Rprox$Index[MM %% (nrow(Rprox)+1) + MM %/% (nrow(Rprox)+1)] >1) {
      
      alpha1  <- q*p  * alpha1  + beta1*q*(1-p)   + n[MM] 
      beta1   <- q    * beta1   + lambda1[MM] 

      a1[ MM] <- alpha1
      b1[ MM] <- beta1
      ps[ MM] <- p
      qs[ MM] <- q
      
      MM <- MM+1
      }
    }
  mat <- cbind(a1, b1, ps, qs)
  return(mat) }

                  #      (p *    q * alpha+ beta * q    *(1-p)) / (q *beta)  
dz <- contvar_parmmat(reparmtable[1:2,5])
Rprox$contvar_cred <- dz[,3]*dz[,1]/dz[,2]+(1-dz[,3])

dz <- joint_parmmat(reparmtable[1:3,2])
Rprox$static_cred <- dz[,3]*dz[,1]/dz[,2]+(1-dz[,3])

dz <- joint_parmmat(reparmtable[1:3,3])
Rprox$incvar_cred <- dz[,3]*dz[,1]/dz[,2]+(1-dz[,3])

dz <- joint_parmmat(reparmtable[1:3,4])
Rprox$decvar_cred <- dz[,3]*dz[,1]/dz[,2]+(1-dz[,3])
rm(dz)

postweight <- subset(Rprox, Year==maxYear)
postweight <- postweight[,c(1, 7:11)]

test$lambda1 <- exp(predict(freqglmIM, test))
test <- merge(test, postweight, by="PolicyNum", all.x=TRUE)

test$contvar_cred[is.na(test$contvar_cred)] <- 1
test$static_cred[ is.na(test$static_cred )] <- 1
test$incvar_cred[ is.na(test$incvar_cred )] <- 1
test$decvar_cred[ is.na(test$decvar_cred )] <- 1

indep_diff   <- tCount - test$lambda1
static_diff  <- tCount - test$lambda1*test$static_cred
incvar_diff  <- tCount - test$lambda1*test$incvar_cred
decvar_diff  <- tCount - test$lambda1*test$decvar_cred
contvar_diff <- tCount - test$lambda1*test$contvar_cred

Poisson.Deviance <- function(pred, obs){
  2*(sum(pred)-sum(obs)+sum(log((obs/pred)^(obs))))/length(pred) }


indep_RMSE <- sqrt(mean(indep_diff^2))
indep_MAE  <- mean(abs( indep_diff))
indep_DEV  <- Poisson.Deviance(test$lambda1, tCount)

static_RMSE <- sqrt(mean(static_diff^2))
static_MAE  <- mean(abs( static_diff))
static_DEV  <- Poisson.Deviance(test$lambda1*test$static_cred, tCount)

incvar_RMSE <- sqrt(mean(incvar_diff^2))
incvar_MAE  <- mean(abs( incvar_diff))
incvar_DEV  <- Poisson.Deviance(test$lambda1*test$incvar_cred, tCount)

decvar_RMSE <- sqrt(mean(decvar_diff^2))
decvar_MAE  <- mean(abs( decvar_diff))
decvar_DEV  <- Poisson.Deviance(test$lambda1*test$decvar_cred, tCount)

contvar_RMSE <- sqrt(mean(contvar_diff^2))
contvar_MAE  <- mean(abs( contvar_diff))
contvar_DEV  <- Poisson.Deviance(test$lambda1*test$contvar_cred, tCount)

valtable <- rbind(
c(indep_RMSE, static_RMSE, incvar_RMSE, decvar_RMSE, contvar_RMSE),
c(indep_MAE,  static_MAE,  incvar_MAE,  decvar_MAE,  contvar_MAE ),
c(indep_DEV,  static_DEV,  incvar_DEV,  decvar_DEV,  contvar_DEV ))

colnames(valtable) <- c("Independent", "Static", "Increasing", "Decreasing", "Constant")
rownames(valtable) <- c("RMSE", "MAE", "DEV")

round(valtable, 4)

############

estable <- summary(freqglmIM)$coefficients[,c(1,4)]
                 

library(knitr)
library(kableExtra)
options(knitr.kable.NA = '')
options(knitr.table.format = "latex")

kable(reparmtable,digits=3,booktabs = T,
      linesep = c("", "", "","",  "","","", "","", "", "","\\hline"),
      bottomrule="\\hhline{======}",escape = FALSE) %>% kable_styling(latex_options = c("hold_position", "scale_down"))

kable(goftable,digits=3,booktabs = T,
      linesep = c("", "", "","",  "","","", "","", "", "","\\hline"),
      bottomrule="\\hhline{======}",escape = FALSE) %>% kable_styling(latex_options = c("hold_position", "scale_down"))

kable(valtable,digits=4,booktabs = T,
      linesep = c("", "", "","",  "","","", "","", "", "","\\hline"),
      bottomrule="\\hhline{======}",escape = FALSE) %>% kable_styling(latex_options = c("hold_position", "scale_down"))
  

#kable(estable,digits=4,booktabs = T,
#      linesep = c("", "", "","",  "","","", "","", "", "","\\hline"),
#      bottomrule="\\hhline{=======}",escape = FALSE) %>% kable_styling(latex_options = c("hold_position", "scale_down")) %>%
#  add_header_above(c(" "=1, " " = 2, "Independent" = 2, "Dependent" = 2)) %>%
#  add_header_above(c(" "=1, "Frequency" = 2,"Severity" = 4))

cred5s <- subset(postweight, Index==5)

ggcol <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n] }

plot( density(cred5s$static_cred), col="black", ylim=c(0,3.8),
      main = "Credibility factors", xlab = bquote(hat(theta)[6]))
lines(density(cred5s$incvar_cred), col=ggcol(8)[5])
lines(density(cred5s$decvar_cred), col=ggcol(8)[6])
lines(density(cred5s$contvar_cred), col=ggcol(8)[8])

legend("topright",
       legend=c("Static", "Increasing", 
                "Decreasing", "Constant"),
       col=c("black",ggcol(8)[c(5,6,8)]), lty=1, cex=0.7)

credtable <- rbind(
summary(cred5s$static_cred),
summary(cred5s$incvar_cred),
summary(cred5s$decvar_cred),
summary(cred5s$contvar_cred))

rownames(credtable) <- c("Static", "Increasing","Decreasing", "Constant")



static_loglik <- function(parm) {
  result <- joint_loglik(p=1, q=1, alpha_0=parm[1]) 
  return(result)  }

static_fit

xsq <- 1:20/10
psq <- 1:20/20
ysq <- rep(NA, 100)
zsq <- matrix(rep(NA, 20*20), nrow=20)
for (i in 1:100) {
  ysq[i] <- static_loglik(xsq[i])
}
plot(xsq, ysq, type='l')

for (i in 1:20) {
  for (j in 1:20) {
  zsq[i,j] <- incvar_loglik(c(psq[i], xsq[j])) }
}

incvar_liks <- zsq

# plot the 3D surface
inc_surface <- persp(psq, xsq, incvar_liks,
      xlab = expression(q),
      ylab = expression(alpha_0),
      zlab = expression(-Loglik),
      xlim = c(0,1),
      main = "Increasing variance model surface plot")

inc_optimum <- trans3d(incvar_fit$par[1],
                       incvar_fit$par[2],
                       incvar_fit$value,
                       pmat = inc_surface)
points(inc_optimum, pch = 16,col = 2)


