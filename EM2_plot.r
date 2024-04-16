EM2_plot<-
function(x, title_plot)
{

  testit = function(x)
  {
    p1 = proc.time()
    Sys.sleep(x)
    proc.time() - p1 # The cpu usage should be negligible
  }
  
	n<-length(x)
	p<-numeric(n)
	f1<-numeric(n)
	f2<-numeric(n)
	w<-matrix(0,n,2)
      eps<-0.01
	itermax<-100
#	hist(x,freq=FALSE)
#
# initialize parameter
#
	alpha0<-runif(1,0,1)
	mu10<-runif(1,min(x),max(x))
	mu20<-runif(1,min(x),max(x))
	sd10<-sd(x[1:round(n/2)])
	sd20<-sd(x[round(n/2):n])

#
# first iteration
#
	iter<-0
	for (i in 1:n) {
	    f1[i]<-dnorm(x[i],mu10,sd10)
	    f2[i]<-dnorm(x[i],mu20,sd20)
	    p[i]<-alpha0*f1[i]+(1-alpha0)*f2[i]
      }
#
# initial likelihood
#
      L0<-sum(log(p))
#  	cat(iter,L0, "\t",alpha0,"\t",mu10,sd10,"\t",mu20,sd10,"\n")
	cat(file="mixtdata.txt",iter,L0,"\t",alpha0,"\t",mu10,sd10,"\t",mu20,sd10,"\n")
	mu0<-c(mu10,mu20)
	sd0<-c(sd10,sd20)
	maxf<-max(mixt2_univ(x,alpha0,mu0,sd0),max(hist(x)$density))
	hist(x, freq=F, main=paste(title_plot, "-  iteration n = ", iter, sep=""), ylim=c(0,maxf+0.005), xlab="weight", ylab="density", col="green") 
	curve(mixt2_univ(x,alpha0,mu0,sd0),min(x),max(x),add=TRUE,col="red",lwd=2)	
	
  readline(prompt="Press [enter] to continue")
	
	L1<-L0

	alpha<-alpha0
	for (iter in 1:itermax){
#
# E-step
#
	w[,1]<-alpha*f1/p
	w[,2]<-(1-alpha)*f2/p
#
# M-step
#
	alpha<-mean(w[,1])
	mu1<-sum(w[,1]*x)/sum(w[,1])
	mu2<-sum(w[,2]*x)/sum(w[,2])
	var1<-sum(w[,1]*((x-mu1)^2))/sum(w[,1])
	var2<-sum(w[,2]*((x-mu2)^2))/sum(w[,2])
	sd1<-sqrt(var1)
	sd2<-sqrt(var2)
	for (i in 1:n) {
	    f1[i]<-dnorm(x[i],mu1,sd1)
	    f2[i]<-dnorm(x[i],mu2,sd2)
	    p[i]<-alpha*f1[i]+(1-alpha)*f2[i]
      }
      L<-sum(log(p))
#	cat(iter,L, "\t",alpha ,"\t",mu1,sd1,"\t",mu2,sd2,"\n")
	cat(file="mixtdata.txt",iter,L,"\t",alpha,"\t",mu1,sd1,"\t",mu2,sd2,"\n",append=T)
	mu<-c(mu1,mu2)
	sd<-c(sd1,sd2)
	maxf<-max(mixt2_univ(x,alpha,mu,sd),max(hist(x)$density))
	hist(x, freq=F, main=paste(title_plot, "-  iteration n = ", iter, sep=""), ylim=c(0,maxf+0.005), xlab="weight", ylab="density", col="green") 
	curve(mixt2_univ(x,alpha,mu,sd),min(x),max(x),add=TRUE,col="red",lwd=2)	
	testit(0.2)
	if (abs(L1-L)<eps) break
	L1<-L
}
	cat(" initial values: \n L0 = ", L0, "\n alpha0 = ", alpha0, "\n", "mu10 = ", mu10, "\n sd10 =", sd10,"\n mu20 =", mu20, "\n sd20 =", sd20,"\n\n")
	cat(" final estimates: \n L = ", L, "\n alpha = ", alpha, "\n mu1 =", mu1, "\n sd1 =", sd1,"\n mu2 =", mu2, "\n sd2 =", sd2,"\n number of iterations =", iter, "\n")
	y<-seq(min(x),max(x),0.1)
	mu<-c(mu1,mu2)
	sd<-c(sd1,sd2)
	maxf<-max((mixt2_univ(x,alpha,mu,sd)))
	hist(x, freq=F, main=title_plot, ,ylim=c(0,maxf+0.005), xlab="weight", ylab="frequency density", col="green") 
	curve(mixt2_univ(x,alpha,mu,sd),min(x),max(x),add=TRUE,col="red",lwd=2)	
	em2<-c(alpha,mu,sd)


}
