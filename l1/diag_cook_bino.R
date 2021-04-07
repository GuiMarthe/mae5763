#-----------------------------------------------------------------------#
# Comandos:
#        > fit.model <- ajuste
#        > attach(dados)
#        > source("diag_cook_bino")
#-----------------------------------------------------------------------#
X <- model.matrix(fit.model)
n <- nrow(X)
p <- ncol(X)
w <- fit.model$weights
W <- diag(w)
H <- solve(t(X)%*%W%*%X)
H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
h <- diag(H)
ts <- resid(fit.model,type="pearson")/sqrt(1-h)
td <- resid(fit.model,type="deviance")/sqrt(1-h)
di <- (h/(1-h))*(ts^2)
#
plot(di,xlab="Indice", ylab="Distancia de Cook",pch=16, main=title)
cut = mean(di) + 4*sd(di)
abline(cut,0,lty=2,lwd=2, cex=2)
identify(di, n=1)
#-----------------------------------------------------------------------#
