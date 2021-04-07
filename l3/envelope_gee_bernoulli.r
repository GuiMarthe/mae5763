#-----------------------------------------------------------------------#
# PROGRAMA PARA GERAR ENVELOPE EM REGRESSÃO LOGÍSTICA COM
# ESTRUTURA DE CORRELAÇÃO  
#-----------------------------------------------------------------------#
#										    
#  Descrição dos argumentos da função:
#						    	
#  modelo: objeto do tipo gee(Y~X, family=binomial
#  (link=logit),...).
#	    
#  dados: banco de dados utilizado no ajuste do modelo.
#			    
#  diag_binomial: objeto gerado a partir da função
#  diag_gee_binomial.
#             
#  opcaoEnvelope: caso opção seja "Normal", então será gerado
#  o envelope usual. Caso o interesse seja analisar o
#  meio-normal, então  opção deve ser "Meio-normal".		    
#										    
#-----------------------------------------------------------------------#

envelope_gee_binomial <- function(modelo, dados, diag_binomial, opcaoEnvelope="Normal")
{
  
  X <- model.matrix(as.formula(paste("~ ", modelo$call$formula[3])), dados) 
  y <- modelo$y
  beta <- coef(modelo)
  R <- modelo$work
  mi <- fitted(modelo)
  individuo <- modelo$id
  estrutura.cor <- modelo$call$corstr
  
  repet <- dim(modelo$work)[1]
  ue <- modelo$nobs/repet
  
  N <- nrow(X)
  p <- ncol(X)
  
  #-----------------------------------------------------------------------#
  # GERAR VARIÁVEIS BINÁRIAS 
  #-----------------------------------------------------------------------#
  
  repl <- 25  
  #Numero de réplicas indicado por Tan, Qu and Kutner
  
  random.y<-array(dim=c(N,repl))
  
  mi.fit <- matrix(mi,ue,byrow=T)
  
  indice<-matrix(c(1:repet),repet,repet,byrow=T)
  
  cont<-1
  for (i in 1:ue)
  {
    cat("i=",i,"\n")
    flag<-T
    k<-0
    TS<-matrix(c(1:repet),repet)
    gama<-0
    media<-mi.fit[i,]
    lambda<-log(matrix(1,repet,repet)+sqrt(diag((1-media)/media))%*%R%*%sqrt(diag((1-media)/media)))
    lambda<-round(lambda,6)
    
    while (flag)
    {
      k<-k+1
      gama.min<-min(lambda[lambda>0])
      rs <- indice[lambda==gama.min]
      S <- unique(rs)
      
      for (j in 1:length(S))
        for (l in j:length(S)) 
          if (lambda[S[j],S[l]]==gama.min) rs<-c(S[j],S[l])
      
      S <- sort(unique(rs))
      
      for (j in 1:repet)
        if (all(lambda[j,S]>0)) S <- c(S,j)
      S <- sort(unique(S))
      
      for (j in 1:length(S))
        for (l in 1:length(S))
          lambda[S[j],S[l]] <- lambda[S[j],S[l]]-gama.min
      
      ts <- rep(0,repet)
      for (j in 1:repet)
        if (is.element(j,S)) ts[j] <- 1
      
      TS <- cbind(TS,ts)
      gama <- c(gama,gama.min)
      
      if (any(lambda)==0) flag <- F
    }
    
    TS <- TS[,-1]
    gama <- gama[-1]
    m <- k
    
    for(j in 1:repl) 
    {
      x <- rpois(1,gama[1])
      for (k in 2:m) x <- cbind(x,rpois(1,gama[k]))
      U <- TS%*%t(x)  
      Z <- as.vector(rep(0,repet))
      Z[U==0] <- 1
      random.y[cont:(cont+repet-1),j] <- Z
    }
    
    cont<-cont+repet  	
  }
  
  #-----------------------------------------------------------------------#
  # ENVELOPE SIMULADO PARA REGRESSÃO LOGÍSTICA 
  #-----------------------------------------------------------------------#
  
  orig.res<-diag_binomial$ResPadronizado
  if (opcaoEnvelope=="Meio-normal") ABSorig.res <- abs(orig.res)
  if (opcaoEnvelope=="Normal") ABSorig.res <- orig.res
  SORTorig.res<-sort(ABSorig.res)
  
  dados2<-cbind(dados,random.y)
  
  random.res<-array(dim=c(N,repl))    
  cat("A simulacao acaba quando k for ",repl,"\n")
  for(k in 1:repl) 
  {
    cat("k:",k,"\n")
    temp.gee.all<-gee(as.formula(paste("random.y[,k] ~ ", modelo$call$formula[3])), family = 
                        binomial(link = logit), data = dados2, id=individuo, corstr=estrutura.cor,
                      scale.fix = T, scale.value = 1)
    summary(temp.gee.all)
    
    y <- temp.gee.all$y
    beta <- coef(temp.gee.all)
    cat("beta:",beta,"\n")
    R <- temp.gee.all$work
    cat("R:",R,"\n")
    mi <- fitted(temp.gee.all)
    
    #Matriz C <- A * Delta	
    #Ligação canônica -> Delta=Identidade
    A <- diag(mi*(1-mi),N)
    C <- A
    
    #Matriz Omega - variância e covariância de y
    Omega <- matrix(0,N,N)
    invOmega <- matrix(0,N,N)
    l <- 1
    while (l<N)
    {
      Omega[l:(l+repet-1),l:(l+repet-1)] <- sqrt(A[l:(l+repet-1),l:(l+repet-1)])%*%R%*%sqrt(A[l:(l+repet-1),l:(l+repet-1)])
      invOmega[l:(l+repet-1),l:(l+repet-1)] <-solve(Omega[l:(l+repet-1),l:(l+repet-1)])
      l <- l+repet
    }
    
    #Matrizes H e W
    W <- C%*%invOmega%*%C
    H <- solve(t(X)%*%W%*%X)
    raizW <- matrix(0,N,N)
    l <- 1
    while (l<N)
    {
      auto<-eigen(W[l:(l+repet-1),l:(l+repet-1)])
      raizW[l:(l+repet-1),l:(l+repet-1)] <- auto$vectors%*%sqrt(diag(auto$values))%*%t(auto$vectors)
      l <- l+repet
    }
    
    H <- raizW%*%X%*%H%*%t(X)%*%raizW
    h <- diag(H)
    
    #Resíduo Padronizado
    random.rsd <- as.vector(rep(0,N))
    part.rsd <- raizW%*%solve(C)%*%(y-mi)
    for (l in 1:N)
    {
      e <- as.vector(rep(0,N))
      e[l] <- 1
      random.rsd[l] <- t(e)%*%part.rsd/sqrt(1-h[l])
    }	
    
    random.res[,k]<-random.rsd
  }
  
  if (opcaoEnvelope=="Meio-normal") ABSrandom.res<-abs(random.res)
  if (opcaoEnvelope=="Normal") ABSrandom.res<-random.res
  
  SORTrandom.res<-array(dim=c(N,repl))    
  
  for(k in 1:repl) {SORTrandom.res[,k]<-sort(ABSrandom.res[,k])}
  
  descritiva<-array(dim=c(N,3))
  for(k in 1:N) {  descritiva[k,1]<-min(SORTrandom.res[k,])
  descritiva[k,2]<-median(SORTrandom.res[k,])
  descritiva[k,3]<-max(SORTrandom.res[k,])}
  
  Z<-array(dim=c(N,1))
  for(i in 1:N)   {Z[i]<-qnorm((i+N-1/8)/(2*N+1/2))}
  
  final<-cbind(Z,descritiva,SORTorig.res)
  faixa <- range(final[,5],final[,2],final[,4])
  
  if (opcaoEnvelope=="Meio-normal")
  {
    
    par(mfrow=c(1,1))
    par(pty="s")
    plot(final[,1],final[,5],xlab="Valor Esperado da Estatística de Ordem Meio-Normal",
         ylab="Valor Absoluto Ordenado do Resíduo Padronizado", ylim=faixa, pch=16)
    par(new=TRUE)
    #
    lines(final[,1],final[,2])
    lines(final[,1],final[,3],lty=2)
    lines(final[,1],final[,4])
  }
  
  if (opcaoEnvelope=="Normal")
  {
    par(pty="s")
    qqnorm(final[,5],xlab="Percentil da N(0,1)",
           ylab="Resíduo de Pearson", ylim=faixa, pch=16,main="")
    par(new=TRUE)
    #
    qqnorm(final[,2],axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1, main="")
    par(new=TRUE)
    qqnorm(final[,4],axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1, main="")
    par(new=TRUE)
    qqnorm(final[,3],axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2, main="")
  }
  
}
#-----------------------------------------------------------------------#
