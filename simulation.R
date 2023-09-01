T = 50
R = 5000

par(mfrow=c(2,2), mar=c(4.1, 4.1, 3.1, 1.1), 
    oma=c(0,0,0,0), mgp=c(2,0.8,0))

ggcol <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n] }


#### Increasing variance ####

theta_inc <- matrix(rep(NA, R*T), nrow=R)
y_inc     <- matrix(rep(NA, R*T), nrow=R)

lambda = 1
qs     = 0.8
qss    = 0.8

for (j in 1:R) {
  set.seed(j)
  alpha = 3
  beta  = 3
  
for (t in 1:T) {
  theta_inc[j,t] <- rgamma(1, shape=alpha, rate=beta) 
      y_inc[j,t] <- rpois( 1, lambda*theta_inc[j,t])
           alpha <- qs *alpha + beta*(qss-qs) + y_inc[j,t]
           beta  <- qss*beta  + lambda
}}

plot( density(theta_inc[,50]), col=ggcol(4)[4],
      main = "Increasing variance", xlab = bquote(theta[t]))
lines(density(theta_inc[,20]), col=ggcol(4)[3])
lines(density(theta_inc[,5 ]), col=ggcol(4)[2])
lines(density(theta_inc[,1 ]), col=ggcol(4)[1])

legend("topright",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)


#### Decreasing variance ####

theta_dec <- matrix(rep(NA, R*T), nrow=R)
y_dec     <- matrix(rep(NA, R*T), nrow=R)

lambda = 1
qs     = 0.8
qss    = 1

for (j in 1:R) {
  set.seed(j)
  alpha = 3
  beta  = 3
  
  for (t in 1:T) {
    theta_dec[j,t] <- rgamma(1, shape=alpha, rate=beta) 
    y_dec[j,t] <- rpois( 1, lambda*theta_dec[j,t])
    alpha <- qs *alpha + beta*(qss-qs) + y_dec[j,t]
    beta  <- qss*beta  + lambda
  }}

plot( density(theta_dec[,1 ]), col=ggcol(4)[1], ylim = c(0,3),
      main = "Decreasing variance", xlab = bquote(theta[t]))
lines(density(theta_dec[,20]), col=ggcol(4)[2])
lines(density(theta_dec[,5 ]), col=ggcol(4)[3])
lines(density(theta_dec[,50]), col=ggcol(4)[4])



legend("topright",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)



#### Converging variance ####
theta_conv <- matrix(rep(NA, R*T), nrow=R)
y_conv     <- matrix(rep(NA, R*T), nrow=R)

lambda = 1
qs     = 0.8
qss    = 0.9

for (j in 1:R) {
  set.seed(j)
  alpha = 3
  beta  = 3
  
  for (t in 1:T) {
    theta_conv[j,t] <- rgamma(1, shape=alpha, rate=beta) 
    y_conv[j,t] <- rpois( 1, lambda*theta_conv[j,t])
    alpha <- qs *alpha + beta*(qss-qs) + y_conv[j,t]
    beta  <- qss*beta  + lambda
  }}

plot( density(theta_conv[,1 ]), col=ggcol(4)[1],
      main = "Converging variance", xlab = bquote(theta[t]),
      ylim = c(0,1.1))
lines(density(theta_conv[,5 ]), col=ggcol(4)[2])
lines(density(theta_conv[,5 ]), col=ggcol(4)[3])
lines(density(theta_conv[,50]), col=ggcol(4)[4])



legend("topright",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)



#### Constant variance ####

theta_const <- matrix(rep(NA, R*T), nrow=R)
y_const     <- matrix(rep(NA, R*T), nrow=R)

lambda = 1
qs     = 1
qss    = 1

for (j in 1:R) {
  set.seed(j)
  alpha = 3
  beta  = 3
  
  for (t in 1:T) {
    theta_const[j,t] <- rgamma(1, shape=alpha, rate=beta) 
    y_const[j,t] <- rpois( 1, lambda*theta_const[j,t])
    alpha <- qs *alpha + beta*(qss-qs) + y_const[j,t]
    beta  <- qss*beta  + lambda
    qs    <- 0.9*qss
    qss   <- 3/(3*0.81 + beta*0.19)
  }}

plot( density(theta_const[,1 ]), col=ggcol(4)[1],
      main = "Constant variance", xlab = bquote(theta[t]))
lines(density(theta_const[,20]), col=ggcol(4)[2])
lines(density(theta_const[,5 ]), col=ggcol(4)[3])
lines(density(theta_const[,50]), col=ggcol(4)[4])


T = 50
R = 5000

par(mfrow=c(2,2), mar=c(4.1, 4.1, 3.1, 1.1), 
    oma=c(0,0,0,0), mgp=c(2,0.8,0))

ggcol <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n] }


#### Increasing variance ####

theta_inc <- matrix(rep(NA, R*T), nrow=R)
y_inc     <- matrix(rep(NA, R*T), nrow=R)

lambda = 1
qs     = 0.8
qss    = 0.8

for (j in 1:R) {
  set.seed(j)
  alpha = 3
  beta  = 3
  
  for (t in 1:T) {
    theta_inc[j,t] <- rgamma(1, shape=alpha, rate=beta) 
    y_inc[j,t] <- rpois( 1, lambda*theta_inc[j,t])
    alpha <- qs *alpha + beta*(qss-qs) + y_inc[j,t]
    beta  <- qss*beta  + lambda
  }}

plot( density(theta_inc[,50]), col=ggcol(4)[4],
      main = "Increasing variance", xlab = bquote(theta[t]))
lines(density(theta_inc[,20]), col=ggcol(4)[3])
lines(density(theta_inc[,5 ]), col=ggcol(4)[2])
lines(density(theta_inc[,1 ]), col=ggcol(4)[1])

legend("topright",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)


#### Decreasing variance ####

theta_dec <- matrix(rep(NA, R*T), nrow=R)
y_dec     <- matrix(rep(NA, R*T), nrow=R)

lambda = 1
qs     = 0.8
qss    = 1

for (j in 1:R) {
  set.seed(j)
  alpha = 3
  beta  = 3
  
  for (t in 1:T) {
    theta_dec[j,t] <- rgamma(1, shape=alpha, rate=beta) 
    y_dec[j,t] <- rpois( 1, lambda*theta_dec[j,t])
    alpha <- qs *alpha + beta*(qss-qs) + y_dec[j,t]
    beta  <- qss*beta  + lambda
  }}

plot( density(theta_dec[,1 ]), col=ggcol(4)[1], ylim = c(0,3),
      main = "Decreasing variance", xlab = bquote(theta[t]))
lines(density(theta_dec[,20]), col=ggcol(4)[2])
lines(density(theta_dec[,5 ]), col=ggcol(4)[3])
lines(density(theta_dec[,50]), col=ggcol(4)[4])



legend("topright",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)

####

plot( density(log(theta_inc[,50]+1e-16)), col=ggcol(4)[4],
      main = "Increasing variance", xlab = bquote(log(theta[t])),
      xlim=c(-8, 8), ylim=c(0, 0.66))
lines(density(log(theta_inc[,20])), col=ggcol(4)[3])
lines(density(log(theta_inc[,5 ])), col=ggcol(4)[2])
lines(density(log(theta_inc[,1 ])), col=ggcol(4)[1])

legend("topleft",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)

plot( density(log(theta_dec[,50])), col=ggcol(4)[4],
      main = "Decreasing variance", xlab = bquote(log(theta[t])),
      xlim=c(-2.5,2))
lines(density(log(theta_dec[,20])), col=ggcol(4)[3])
lines(density(log(theta_dec[,5 ])), col=ggcol(4)[2])
lines(density(log(theta_dec[,1 ])), col=ggcol(4)[1])

legend("topleft",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)

plot( density(log(theta_conv[,50])), col=ggcol(4)[4],
      main = "Converging variance", xlab = bquote(log(theta[t])))
lines(density(log(theta_conv[,20])), col=ggcol(4)[3])
lines(density(log(theta_conv[,5 ])), col=ggcol(4)[2])
lines(density(log(theta_conv[,1 ])), col=ggcol(4)[1])

legend("topleft",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)


plot( density(log(theta_const[,50])), col=ggcol(4)[4],
      main = "Constant variance", xlab = bquote(log(theta[t])))
lines(density(log(theta_const[,20])), col=ggcol(4)[3])
lines(density(log(theta_const[,5 ])), col=ggcol(4)[2])
lines(density(log(theta_const[,1 ])), col=ggcol(4)[1])

legend("topleft",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)

plot( density(theta_conv[,1 ]), col=ggcol(4)[1],
      main = "Converging variance", xlab = bquote(theta[t]),
      ylim = c(0,1.1))
lines(density(theta_conv[,5 ]), col=ggcol(4)[2])
lines(density(theta_conv[,5 ]), col=ggcol(4)[3])
lines(density(theta_conv[,50]), col=ggcol(4)[4])



legend("topright",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)

plot( density(theta_const[,1 ]), col=ggcol(4)[1],
      main = "Constant variance", xlab = bquote(theta[t]))
lines(density(theta_const[,20]), col=ggcol(4)[2])
lines(density(theta_const[,5 ]), col=ggcol(4)[3])
lines(density(theta_const[,50]), col=ggcol(4)[4])



legend("topright",
       legend=c(bquote(theta[1 ]), bquote(theta[5 ]),
                bquote(theta[20]), bquote(theta[50])),
       col=ggcol(4), lty=1, cex=1)


varss <- rep(NA, T)
for (i in 1:T) {
  varss[i] <- var(theta_const[,i])
}


plot( 1:50, theta_inc[  1,1:50], col="blue", ylim=c(0, 3),
      main = "1st trajectory", ylab = bquote(theta[t]), type='l', xlab=bquote(t))
lines(1:50, theta_dec[  1,1:50], col="red")
lines(1:50, theta_conv[ 1,1:50], col="green")
lines(1:50, theta_const[1,1:50], col="purple")

legend("topleft",
       legend=c("Increasing", "Decreasing",
                "Converging", "Constant"),
       col=c("blue", "red", "green", "purple"), lty=1, cex=0.7)



plot( 1:50, theta_inc[  2,1:50], col="blue", ylim=c(0, 3),
      main = "2nd trajectory", ylab = bquote(theta[t]), type='l', xlab=bquote(t))
lines(1:50, theta_dec[  2,1:50], col="red")
lines(1:50, theta_conv[ 2,1:50], col="green")
lines(1:50, theta_const[2,1:50], col="purple")

legend("topleft",
       legend=c("Increasing", "Decreasing",
                "Converging", "Constant"),
       col=c("blue", "red", "green", "purple"), lty=1, cex=0.7)


plot( 1:50, theta_inc[  3,1:50], col="blue", ylim=c(0, 3),
      main = "3rd trajectory", ylab = bquote(theta[t]), type='l', xlab=bquote(t))
lines(1:50, theta_dec[  3,1:50], col="red")
lines(1:50, theta_conv[ 3,1:50], col="green")
lines(1:50, theta_const[3,1:50], col="purple")

legend("topleft",
       legend=c("Increasing", "Decreasing",
                "Converging", "Constant"),
           col=c("blue", "red", "green", "purple"), lty=1, cex=0.7)

plot( 1:50, theta_inc[  4,1:50], col="blue", ylim=c(0, 6),
      main = "4th trajectory", ylab = bquote(theta[t]), type='l', xlab=bquote(t))
lines(1:50, theta_dec[  4,1:50], col="red")
lines(1:50, theta_conv[ 4,1:50], col="green")
lines(1:50, theta_const[4,1:50], col="purple")

legend("topleft",
       legend=c("Increasing", "Decreasing",
                "Converging", "Constant"),
       col=c("blue", "red", "green", "purple"), lty=1, cex=0.7)


plot( 1:50, log(theta_inc[  1,1:50]), col="blue", ylim = c(-3,2),
      main = "1st trajectory", ylab = bquote(log(theta[t])), type='l', xlab=bquote(t))
lines(1:50, log(theta_dec[  1,1:50]), col="red")
lines(1:50, log(theta_conv[ 1,1:50]), col="green")
lines(1:50, log(theta_const[1,1:50]), col="purple")

legend("top",
       legend=c("Increasing", "Decreasing",
                "Converging", "Constant"),
       col=c("blue", "red", "green", "purple"), lty=1, seg.len = 2,
       inset=0.01, cex=0.7, horiz=TRUE,
       box.lty=0, x.intersp = 0.5, xjust=0.5)

plot( 1:50, log(theta_inc[  2,1:50]), col="blue", ylim = c(-3,2),
      main = "2nd trajectory", ylab = bquote(log(theta[t])), type='l', xlab=bquote(t))
lines(1:50, log(theta_dec[  2,1:50]), col="red")
lines(1:50, log(theta_conv[ 2,1:50]), col="green")
lines(1:50, log(theta_const[2,1:50]), col="purple")

legend("top",
       legend=c("Increasing", "Decreasing",
                "Converging", "Constant"),
       col=c("blue", "red", "green", "purple"), lty=1, seg.len = 2,
       inset=0.01, cex=0.7, horiz=TRUE,
       box.lty=0, x.intersp = 0.5, xjust=0.5)


plot( 1:50, log(theta_inc[  3,1:50]), col="blue", ylim = c(-3,2),
      main = "3rd trajectory", ylab = bquote(log(theta[t])), type='l', xlab=bquote(t))
lines(1:50, log(theta_dec[  3,1:50]), col="red")
lines(1:50, log(theta_conv[ 3,1:50]), col="green")
lines(1:50, log(theta_const[3,1:50]), col="purple")

legend("top",
       legend=c("Increasing", "Decreasing",
                "Converging", "Constant"),
       col=c("blue", "red", "green", "purple"), lty=1, seg.len = 2,
       inset=0.01, cex=0.7, horiz=TRUE,
       box.lty=0, x.intersp = 0.5, xjust=0.5)

plot( 1:50, log(theta_inc[  4,1:50]), col="blue", ylim = c(-2,3),
      main = "4th trajectory", ylab = bquote(log(theta[t])), type='l', xlab=bquote(t))
lines(1:50, log(theta_dec[  4,1:50]), col="red")
lines(1:50, log(theta_conv[ 4,1:50]), col="green")
lines(1:50, log(theta_const[4,1:50]), col="purple")

legend("top",
       legend=c("Increasing", "Decreasing",
                "Converging", "Constant"),
       col=c("blue", "red", "green", "purple"), lty=1, seg.len = 2,
       inset=0.01, cex=0.7, horiz=TRUE,
       box.lty=0, x.intersp = 0.5, xjust=0.5)







plot( 1:50, theta_inc[  1,1:50], col="blue", ylim=c(0, 6),
      main = "Increasing variance", ylab = bquote(theta[t]), type='l', xlab=bquote(t))
lines(1:50, theta_inc[  2,1:50], col="red")
lines(1:50, theta_inc[ 3,1:50], col="green")
lines(1:50, theta_inc[4,1:50], col="purple")

legend("topleft",
       legend=c("1st", "2nd",
                "3rd", "4th"),
       col=c("blue", "red", "green", "purple"), lty=1, cex=0.7)



plot( 1:50, theta_dec[ 1,1:50], col="blue", ylim=c(0, 6),
      main = "Decreasing variance", ylab = bquote(theta[t]), type='l', xlab=bquote(t))
lines(1:50, theta_dec[ 2,1:50], col="red")
lines(1:50, theta_dec[ 3,1:50], col="green")
lines(1:50, theta_dec[ 4,1:50], col="purple")
legend("topleft",
       legend=c("1st", "2nd",
                "3rd", "4th"),
       col=c("blue", "red", "green", "purple"), lty=1, cex=0.7)


plot( 1:50, theta_conv[ 1,1:50], col="blue", ylim=c(0, 6),
      main = "Converging variance", ylab = bquote(theta[t]), type='l', xlab=bquote(t))
lines(1:50, theta_conv[ 2,1:50], col="red")
lines(1:50, theta_conv[ 3,1:50], col="green")
lines(1:50, theta_conv[ 4,1:50], col="purple")

legend("topleft",
       legend=c("1st", "2nd",
                "3rd", "4th"),
       col=c("blue", "red", "green", "purple"), lty=1, cex=0.7)

plot( 1:50, theta_const[ 1,1:50], col="blue", ylim=c(0, 6),
      main = "Constant variance", ylab = bquote(theta[t]), type='l', xlab=bquote(t))
lines(1:50, theta_const[ 2,1:50], col="red")
lines(1:50, theta_const[ 3,1:50], col="green")
lines(1:50, theta_const[ 4,1:50], col="purple")
legend("topleft",
       legend=c("1st", "2nd",
                "3rd", "4th"),
       col=c("blue", "red", "green", "purple"), lty=1, cex=0.7)