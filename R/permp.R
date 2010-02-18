permp <- function(x, nperm, n1, n2, total.nperm=NULL)
#	 Correctly calculates permutation probabilities
#	 Belinda Phipson
#	 16 February 2010. Last modified 16 February 2010.
{
	 if(is.null(total.nperm)){
		if((n1+n2)/n1==2) total.nperm<-(choose((n1+n2),n1))/2
		else total.nperm<-choose((n1+n2),n1)
	 }
	 z<-gauss.quad.prob(128,l=0,u=0.5/total.nperm)
	 prob<-rep(z$nodes,length(x))
	 x2<-rep(x,each=128)
	 Y<-matrix(pbinom(x2,prob=prob,size=nperm),128,length(x))
	 int<-0.5/total.nperm*colSums(z$weights*Y)
	 (x+1)/(nperm+1)-int
}