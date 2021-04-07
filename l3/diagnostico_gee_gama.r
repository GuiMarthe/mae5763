#-----------------------------------------------------------------------#
#
# DIAGN?STICO GEE PARA DISTRIBUI??O GAMA - LIGA??O LOG	
#
#-----------------------------------------------------------------------#
#  Com esse ajuste, calculamos algumas medidas para an?lise de
#  diagn?stico e constru?mos os gr?ficos conforme Tan, Qu e
#  Kutner (1997) e Venezuela (2003).
#-----------------------------------------------------------------------#
#		 								                         
#  Descri??o dos argumentos da fun??o:	
#					             
#  modelo: objeto do tipo gee(Y~X, family=Gamma(link=log),
#  ...).
#                    
#  dados: banco de dados utilizado no ajuste do modelo.
#                            
#  umaJanela: argumento l?gico que d? a op??o de imprimir
#  todos os gr?ficos 	da an?lise de diagn?stico numa mesma
#  janela ou um gr?fico em 	cada janela. Por default, essa
#  op??o ? TRUE.
#   		            
#  selGraficos: vetor de uns e zeros que d? a op??o de n?o
#  imprimir todos os gr?ficos da an?lise de diagn?stico.
#  Para n?o imprimir um gr?fico dentre os quatro indicados,
#  basta colocar 0 na posi??o correspondente do gr?fico.
#  Altera??o na configura??oo default c(1,1,1,1) deve ser
#  utilizada junto com uma Janela=FALSE.
#    
#  identifica: vetor de quantidades de observa??eses a serem
#  identificadas em cada gr?fico da an?lise de diagn?stico.
#  Por exemplo, para identificar um ponto no 2o. gr?fico e 1
#  no 4o. gr?fico, basta fazer  identifica=c(0,1,0,1).
#                                                 
#-----------------------------------------------------------------------#

diag_gee_gama <- function(modelo, dados, umaJanela=T, selGraficos = c(1,1,1,1), 
				identifica=c(0,0,0,0))
{

X <- model.matrix(as.formula(paste("~ ", modelo$call$formula[3])), dados) 
y <- modelo$y
beta <- coef(modelo)
R <- modelo$work
mi <- fitted(modelo)
individuo <- modelo$id

repet <- dim(modelo$work)[1]
ue <- modelo$nobs/repet

N <- nrow(X)
p <- ncol(X)

#Calculo do Res?duo de Pearson
r<-(y-mi)*(1/sqrt(mi^2))

#Calculo de phi
invphi<-as.numeric(sum(r^2)/(N-p))
phi<-1/invphi

#Matriz C <- A * Delta  #Liga??o canonica -> Delta=Identidade
A <- diag(mi^2,N)
Delta <- diag(1/mi,N)
C <- Delta%*%A

#Matriz Omega - vari?ncia e covari?ncia de y
Omega <- matrix(0,N,N)
invOmega <- matrix(0,N,N)
l <- 1
while (l<N)
  {
  cat(l, ' | ', repet, '\n')
  Omega[l:(l+repet-1),l:(l+repet-1)] <- sqrt(A[l:(l+repet-1),l:(l+repet-1)])%*%R[1:repet,1:repet]%*%sqrt(A[l:(l+repet-1),l:(l+repet-1)])
  invOmega[l:(l+repet-1),l:(l+repet-1)] <-solve(Omega[l:(l+repet-1),l:(l+repet-1)])
  l <- l+repet
  }
Omega <- invphi*Omega
invOmega <- phi*invOmega

#Matriz H e W
W <- C%*%invOmega%*%C
H <- solve(t(X)%*%W%*%X)
raizW <- matrix(0,N,N)
l <- 1
while (l<N)
  {
  auto<-eigen(W[l:(l+repet-1),l:(l+repet-1)])
  raizW[l:(l+repet-1),l:(l+repet-1)] <- auto$vectors%*%sqrt(diag(auto$values,repet))%*%t(auto$vectors)
  l <- l+repet
  }
H <- raizW%*%X%*%H%*%t(X)%*%raizW
h <- diag(H)

#Ponto Alavanca por UE
hue<-as.vector(rep(0,ue))
haux <- matrix(h,ue,repet,byrow=T)
for (i in 1:ue)
  hue[i] <- sum(haux[i,])/repet

#Res?duo Padronizado
rsd <- as.vector(rep(0,N))
part.rsd <- raizW%*%solve(C)%*%(y-mi)
for (l in 1:N)
  {
  e <- as.vector(rep(0,N))
  e[l] <- 1
  rsd[l] <- t(e)%*%part.rsd/sqrt(1-h[l])
  }	

#Dist?ncia de Cook
cd <- as.vector(rep(0,N))
for (l in 1:N)
  {
  cd[l] <- (rsd[l])^2*h[l]/((1-h[l]))
  }	
  
#-----------------------------------------------------------------------#
#  Constru??o de Gr?ficos de Diagn?stico 
#-----------------------------------------------------------------------#
# Para  identificar os pontos que  mais  se destacam  em algum
# gr?fico, use o comando  identify(...) colocando em n o 
# n?mero de pontos que se destacaram.
#
#-----------------------------------------------------------------------#

if (umaJanela) par(mfrow=c(2,2))
if (!umaJanela) par(mfrow=c(1,1))

labelsGraf = paste("(", individuo, ",", rep(1:repet, ue), ")", sep="")

if (selGraficos[1])
{
plot(individuo,h,xlab="Unidade Experimental", ylab="H por Observação", pch=16)
cut <- 2*p/N
abline(cut,0,lty=2)
if (identifica[1]>0) identify(dados$individuo,h,labels=individuo,n=identifica[1])
#if (!umaJanela) savePlot("graf1.jpeg", type="jpeg")
}

if (selGraficos[2])
{
plot(hue,xlab="Unidade Experimental", ylab="H por Unidade Experimental", pch=16)
abline(cut,0,lty=2)
if (identifica[2]>0) identify(hue,labels=dados$individuo, n=identifica[2])
#if (!umaJanela) savePlot("graf2.jpeg", type="jpeg")
}

if (selGraficos[3])
{
plot(individuo,cd,xlab="Unidade Experimental", ylab="Distância de Cook", pch=16)
if (identifica[3]>0) identify(individuo,labels=labelsGraf,cd,n=identifica[3])
#if (!umaJanela) savePlot("graf3.jpeg", type="jpeg")
}

if (selGraficos[4])
{
plot(individuo,rsd,xlab="Unidade Experimental", ylab="Resíduo Padronizado", pch=16)
if (identifica[4]>0) identify(individuo,labels=labelsGraf,rsd,n=identifica[4])
#if (!umaJanela) savePlot("graf4.jpeg", type="jpeg")
}
saida = list(ResPadronizado=rsd, DistCook=cd, h=h, hue=hue, phi=phi)
return(invisible(saida))

}
#-----------------------------------------------------------------------#