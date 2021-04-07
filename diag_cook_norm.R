#-----------------------------------------------------------------------#
# Comandos:
#        > fit.model <- ajuste
#        > attach(dados)
#       title <- 'title'
#        > source("diag_cook_norm")     
#-----------------------------------------------------------------------#
X <- model.matrix(fit.model)
n <- nrow(X)
p <- ncol(X)
H <- X%*%solve(t(X)%*%X)%*%t(X)
h <- diag(H)
r <- resid(fit.model)
s <- sqrt(sum(r*r)/(n-p))
ts <- r/(s*sqrt(1-h))
di <- (1/p)*(h/(1-h))*(ts^2)
si <- lm.influence(fit.model)$sigma
tsi <- r/(si*sqrt(1-h))
a <- max(tsi)
b <- min(tsi)
#
plot(di,xlab="Indice", ylab="Distancia de Cook", pch=16, main=title)
cut = mean(di) + 2.5*sd(di)
abline(cut,0,lty=2)
identify(di, n=3)
#-----------------------------------------------------------------------#
