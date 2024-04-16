mixt2_univ<-
function(x,alpha,mu,sd)

{
	f1<-dnorm(x,mu[1],sd[1])
    	f2<-dnorm(x,mu[2],sd[2])
    	p=alpha*f1+(1-alpha)*f2
	p
}
	
