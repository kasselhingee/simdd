bfind=function(lambda) {
  q=length(lambda)
  fb=function(b) 1-sum(1/(b+2*lambda))
  if(sum(lambda^2)==0) b0=q
  else b0=stats::uniroot(fb,interval=c(1,q))$root
  b0
}

dimq=function(A) {
  q=0
  if(is.vector(A)) q=length(A)
  if(is.matrix(A)) {if(nrow(A)==ncol(A)) q=nrow(A)}
  q
}

dimset=function(mu,A) {
  q1=length(mu); q2=dimq(A)
  if(q1>1 & q2>1 & q1!=q2) stop("dimset: mismatch of dimensions\n")
  q=max(q1,q2)
  q
}

rBingham=function(nsim,Aplus,q=dimq(Aplus),mtop=1000) {
  values=NULL; ntry=0; nleft=nsim; mloop=0; minfg=Inf; maxfg=-Inf #G
  Aminus=-Aplus # so Aminus is  minus the concentration parameters
  qA=dimq(Aminus)
  if(q<=1 & qA <=1) stop("error in rBingham: need to set a value of q>1 \n")
  if(q>1 & qA>1 & q!=qA) stop("error in rBingham: dimension mismatch for q and Aplus \n")
  q=max(q,qA)
  if(is.vector(Aminus)){switch=0; lambda=Aminus; Gamma=diag(q)}
  if(is.matrix(Aminus)) {switch=1; ee=eigen(Aminus,symmetric=TRUE); lambda=ee$values; Gamma=ee$vectors}
  lambda=lambda-min(lambda)
  b0=bfind(lambda)
  phi=1+2*lambda/b0
  while(nleft>0 & mloop<mtop) {#G  
    x=matrix(stats::rnorm(nleft*q),nleft,q)*(matrix(1,nleft,1)%*%
                 matrix(1/sqrt(phi),1,q))
    r=sqrt((x*x)%*%rep(1,q))
    x=x/(matrix(r,nleft,1)%*%matrix(1,1,q)) # so x is acg
    u=((x*x)*(matrix(1,nleft,1)%*%matrix(lambda,1,q)))%*%rep(1,q)
    v=stats::runif(nleft)
    logfg=-u+.5*(q-b0)+(q/2)*log((1+2*u/b0)*b0/q)
    fb=exp(logfg)
    fg=exp(-u+.5*(q-b0))*((1+2*u/b0)*b0/q)^(q/2)
    minfg=min(minfg,min(fg)); maxfg=max(maxfg,max(fg))#G
    ind=(v<fg) #G
    nacc=sum(ind); ntry=ntry+nleft; nleft=nleft-nacc; mloop=mloop+1#G
     if(nacc>0) values=rbind(values,x[ind,,drop=FALSE])#G
  }#G
  if(switch==1) values=values%*%t(Gamma)
  summ=c(ntry,(nsim-nleft)/ntry,(nleft==0),mloop=mloop,minfg,maxfg)#G
  names(summ)=c("ntry","efficiency","success","mloops","minfg","maxfg") #G
  attr(values,"summary")=summ
  values
}

rFisherBingham=function(nsim,mu=0,Aplus=0, q=dimset(mu,Aplus), mtop=1000) {
  qA=dimset(mu,Aplus)
  if(q<=1 & qA <=1) stop("error in rFisherBingham: need to set a value of q>1 \n")
  if(q>1 & qA>1 & q!=qA) stop("error in rFisherBingham: dimension mismatch for q and mu and Aplus \n")
  q=max(q,qA)
  Aminus=-Aplus # so Aminus is  minus the concentration parameters
  values=NULL; ntry=0; nleft=nsim; mloop=0; minfg=Inf; maxfg=-Inf;#G
  mu0=diag(q)[,q]
  if(length(mu)==1) {switchmu=0; kappa=mu}
  if(length(mu)==q) {
    kappa=sqrt(sum(mu^2))
    if(kappa>0) mu0=mu/kappa
  }
  if(is.vector(Aminus) &length(Aminus)<=1) Aminus=rep(0,q)
  if(is.vector(Aminus)) Aminus=diag(Aminus)
  A1minus=Aminus-(kappa/2)*mu0%*%t(mu0)
  ee=eigen(A1minus,symmetric=TRUE); lambda1=ee$values; Gamma1=ee$vectors
  lambda1min=min(lambda1); lambda1=lambda1-lambda1min
  b0=bfind(lambda1)
  phi=1+2*lambda1/b0
  while(nleft>0 & mloop<mtop) {
    x=matrix(stats::rnorm(nleft*q),nleft,q)*(matrix(1,nleft,1)%*%
                 matrix(1/sqrt(phi),1,q))
    r=sqrt((x*x)%*%rep(1,q))
    x=x/(matrix(r,nleft,1)%*%matrix(1,1,q)) # so x is acg
    xlhs=x%*%t(Gamma1)
    ulhs=kappa*(xlhs%*%mu0) - apply(xlhs*(xlhs%*%Aminus),1,sum)
    urhs=((x*x)*(matrix(1,nleft,1)%*%matrix(lambda1,1,q)))%*%rep(1,q)
    v=stats::runif(nleft)
    logfg=ulhs+.5*(q-b0)-kappa/2+lambda1min+(q/2)*log((1+2*urhs/b0)*b0/q)
    fg=exp(logfg)
    minfg=min(minfg,min(fg)); maxfg=max(maxfg,max(fg))
    ind=(v<fg) 
    nacc=sum(ind) 
    x1=xlhs[ind,,drop=FALSE]
    nacc=sum(ind); ntry=ntry+nleft; nleft=nleft-nacc; mloop=mloop+1
     if(nacc>0) values=rbind(values,x1)
  }
  summ=c(ntry,(nsim-nleft)/ntry,(nleft==0),mloop=mloop,minfg,maxfg)
  names(summ)=c("ntry","efficiency","success","mloops","minfg","maxfg") 
  attr(values,"summary")=summ
  if ((mloop == mtop) && (nleft > 0)){warning("Iterations reached mtop. Exiting before completing requested nsim simulations.")}
  values
}

is.s3=function(vv,eps=1e-25) {
  # check to see if vv is a  unit 4-vector
  # or a matrix whose cols are
  isgood=FALSE
  vmat=vv
  if(is.vector(vv)) {n=1; vmat=matrix(vv,1,4)}
  if(is.matrix(vmat)&nrow(vmat)>=1 &ncol(vmat)==4) {
    n=nrow(vmat)
    norm=sqrt(apply(vmat*vmat,1,sum))
    discrepency=(norm-1)^2
    if(max(discrepency)<eps) isgood=TRUE
  }
  isgood
}

is.so3=function(XX,eps=1e-25) {
  isgood=FALSE
  Xarray=XX
  if(is.matrix(XX) & nrow(XX)==3 & ncol(XX)==3) Xarray=array(XX,dim=c(1,3,3))
  if(is.array(Xarray) & length(dim(Xarray))==3 & dim(Xarray)[1]>=1 &
     dim(Xarray)[2]==3 & dim(Xarray)[3]==3) {
    n=dim(Xarray)[1]
    discrepency=rep(0,n)
    for(i in 1:n) {
      X=Xarray[i,,]
      discrepency[i]=sum((t(X)%*%X-diag(3))^2) + (det(X)-1)^2
    }
    if(max(discrepency)<eps) isgood=TRUE
  }
  isgood
}

s3toso3=function(vv) {
# assumes v is a unit vector in R^4 or matrix whose rows are
# output is a 3 by 3 rotation matrix X or array
  isgood=is.s3(vv)
  if(isgood==FALSE) return("argument must be 4-vector or matrix whose rows are")
  vmat=vv
  isvector=FALSE
  if(is.vector(vv)) {isvector=TRUE; vmat=matrix(vv,1,4)}
  n=nrow(vmat)
  Xarray=array(0,dim=c(n,3,3))
  for(i in 1:n) {
    v=vmat[i,]
    X=matrix(0,3,3)
    X[1,1]= v[1]^2+v[2]^2-v[3]^2-v[4]^2
    X[2,2]= v[1]^2-v[2]^2+v[3]^2-v[4]^2
    X[3,3]= v[1]^2-v[2]^2-v[3]^2+v[4]^2
    X[1,2]=2*( v[2]*v[3]-v[1]*v[4]); X[2,1]=2*( v[2]*v[3]+v[1]*v[4])
    X[1,3]=2*( v[1]*v[3]+v[2]*v[4]); X[3,1]=2*(-v[1]*v[3]+v[2]*v[4])
    X[2,3]=2*(-v[1]*v[2]+v[3]*v[4]); X[3,2]=2*( v[1]*v[2]+v[3]*v[4]) 
    Xarray[i,,]=X
  }
  XX=Xarray
  if(n==1) XX=Xarray[1,,]
  XX
}

so3tos3=function(XX) {
  # assumes XX is a 3 by 3 rotation matrix
  # or n by 3 by 3 array whose blocks are.
  # output is a (matrix of) unit vectors v in R^4;  let V=vv^T
  isgood=is.so3(XX)
  if(isgood==FALSE) return("argument must be 3 by 3 rotation matrix or 3 by 3 by n array whose initial blocks are")
  Xarray=XX; ismatrix=FALSE
  if(is.matrix(XX)) {ismatrix=TRUE; Xarray=array(XX,dim=c(1,3,3))}
  n=dim(Xarray)[1]; vmat=matrix(0,n,4)
  for(i in 1:n) {
    X=Xarray[i,,]
    V=matrix(0,4,4)
    V[1,2]=(X[3,2]-X[2,3]); V[3,4]=(X[3,2]+X[2,3])
    V[1,3]=(X[1,3]-X[3,1]); V[2,4]=(X[1,3]+X[3,1])
    V[1,4]=(X[2,1]-X[1,2]); V[2,3]=(X[2,1]+X[1,2])
    V=(V+t(V))/4
    V[1,1]=(1+X[1,1]+X[2,2]+X[3,3])/4
    V[2,2]=(1+X[1,1]-X[2,2]-X[3,3])/4
    V[3,3]=(1-X[1,1]+X[2,2]-X[3,3])/4
    V[4,4]=(1-X[1,1]-X[2,2]+X[3,3])/4
    ind=order(diag(V))[4]
    v=V[,ind]; v=v/sqrt(sum(v^2))
    vmat[i,]=v
    }
  vv=vmat
  if(ismatrix==TRUE) vv=vmat[1,]
  vv
}

mf2b=function(Fmat) {
  # calculate 4 by 4 symmetric parameter A matrix for Bingham
  # input is 3 by 3 parameter matrix Fmat for matrix Fisher
  A=matrix(0,4,4)
  A[1,2]=Fmat[3,2]-Fmat[2,3]; A[3,4]=Fmat[3,2]+Fmat[2,3] 
  A[1,3]=Fmat[1,3]-Fmat[3,1]; A[2,4]=Fmat[1,3]+Fmat[3,1]
  A[1,4]=Fmat[2,1]-Fmat[1,2]; A[2,3]=Fmat[2,1]+Fmat[1,2]
  A=A+t(A)
  A[2,2]=-2*(Fmat[2,2]+Fmat[3,3])
  A[3,3]=-2*(Fmat[1,1]+Fmat[3,3])
  A[4,4]=-2*(Fmat[1,1]+Fmat[2,2])
  A=A-mean(diag(A))*diag(4)
  A
}
b2mf=function(A) {
  # calculate 3 by 3 Fmat parameter matrix for matrix Fisher
  # input A is a 4 by 4 symmetric matrix 
  Fmat=matrix(0,3,3)
   A=A-A[1,1]*diag(4) # convention for Fmat calculation
  Fmat[1,1]=-(A[3,3]+A[4,4]-A[2,2])/4
  Fmat[2,2]=-(A[2,2]+A[4,4]-A[3,3])/4
  Fmat[3,3]=-(A[2,2]+A[3,3]-A[4,4])/4
  Fmat[2,3]=(A[3,4]-A[1,2])/2; Fmat[3,2]=(A[3,4]+A[1,2])/2
  Fmat[1,2]=(A[2,3]-A[1,4])/2; Fmat[2,1]=(A[2,3]+A[1,4])/2
  Fmat[3,1]=(A[2,4]-A[1,3])/2; Fmat[1,3]=(A[2,4]+A[1,3])/2
  Fmat
}

  # rmatrix.fisher3 function definition
  # simulate from the  3 by 3 matrix Fisher distribution
  # Fmat is the 3 by 3 parameter matrix.
  # output is a 3 by 3 by nsim array of rotation matrices.

rFisher.SO3=function(nsim,Fmat) s3toso3(rBingham(nsim,mf2b(Fmat)))

rBingham.Grassmann=function(nsim,Aplus=0, q=dimq(Aplus),r=1,mtop=1000) {
  ndone=0; nleft=nsim; mloop=0; ntry=0
  values=array(0,dim=c(nsim,q,r))
  minfg=Inf; maxfg=-Inf #G
  Aminus=-Aplus # so Aminus is  minus the concentration parameters
  qA=dimq(Aminus)
  if(q<=1 & qA <=1) stop("error in rBingham.Grassmann: need to set a value of q>1 \n")
  if(q>1 & qA>1 & q!=qA) stop("error in rBingham.Grassmann: dimension mismatch for q and Aplus \n")
  q=max(q,qA)
  if(is.vector(Aminus) &length(Aminus)==1) Aminus=rep(0,q)
  if(is.vector(Aminus)){switch=0; lambda=Aminus; Gamma=diag(q)}
  if(is.matrix(Aminus)) {
    switch=1
    ee=eigen(Aminus,symmetric=TRUE); lambda=ee$values; Gamma=ee$vectors
  }
  lambda=lambda-min(lambda)
  b0=bfind(lambda)
  u0=(q-b0)/2
  phi=1+2*lambda/b0
  rescale.phi=function(v) v/sqrt(phi)
  on=function(Z) qr.Q(qr(Z))
  # so on othonormalizes the columns of Z
  extract.lambda=function(X) eigen(t(X)%*%(matrix(lambda,q,r)*X))$values
  extract.phi=function(X) eigen(t(X)%*%(matrix(phi,q,r)*X))$values
  rotate=function(v) {
    if(switch==0) return(v)
    else return(Gamma%*%v)
  }
  comp.fn=function(u) (q/2)*log(1+2*u/b0)-u-(q/2)*log(1+2*u0/b0)+u0
  lf=1
  for(j in 1:r) if(lambda[j]>u0) lf=lf*exp(-comp.fn(lambda[j]))
  while(nleft>0 & mloop<mtop) {#G
    X1=array(stats::rnorm(nleft*q*r),dim=c(nleft,q,r))
    X2=X1; X3=X2; U.lambda=matrix(0,nleft,r); U.phi=U.lambda
    for(i in 1:nleft) {
        for(j in 1:r) X2[i,,j]=rescale.phi(X1[i,,j])
        X3[i,,]=on(X2[i,,])
        U.lambda[i,]=extract.lambda(X3[i,,])
        U.phi[i,]=extract.phi(X3[i,,])
    }
    v=stats::runif(nleft)
    fg=lf*exp(-apply(U.lambda,1,sum)+.5*(q-b0)*r)*(apply(U.phi*b0/q,1,prod))^(q/2)
    minfg=min(minfg,min(fg)); maxfg=max(maxfg,max(fg))#G
    ind=(v<fg) #G
    nacc=sum(ind)
    if(nacc>0)  {
      values[(ndone+1):(ndone+nacc),,]=
        aperm(apply(X3[ind,,,drop=FALSE],c(1,3),rotate),c(2,1,3))
      ndone=ndone+nacc; ntry=ntry+nleft; nleft=nleft-nacc; mloop=mloop+1
   }
  }
  summ=c(ntry,(nsim-nleft)/ntry,(nleft==0),mloop=mloop,minfg,maxfg)#G
  names(summ)=c("ntry","eff","success","mloops","minfg","maxfg") #G
  attr(values,"summary")=summ
  values
}

rBessel=function(nsim,k1,k2,alpha,mtop=1000) {
  values=NULL; ntry=0; nleft=nsim; mloop=0; minfg=Inf; maxfg=-Inf #G
  lambda=c(0,.5*(k1-alpha^2/k2)); lambdamin=0
  if(lambda[2]<0) {lambdamin=lambda[2]; lambda=lambda-lambda[2]}
  q=2
  b0=bfind(lambda)
  phi=1+2*lambda/b0
  den=besselI(k2,0)
  while(nleft>0 & mloop<mtop) {#G  
    x=matrix(stats::rnorm(nleft*q),nleft,q)*(matrix(1,nleft,1)%*%
                 matrix(1/sqrt(phi),1,q))
    r=sqrt((x*x)%*%rep(1,q))
    x=x/(matrix(r,nleft,1)%*%matrix(1,1,q)) # so x is acg
    u=((x*x)*(matrix(1,nleft,1)%*%matrix(lambda,1,q)))%*%rep(1,q)
    v=stats::runif(nleft)
    f=exp(k1*(x[,1]-1)+lambdamin)*besselI(sqrt(k2^2+alpha^2*x[,2]^2),0)/den
    gi=exp(.5*(q-b0))*((1+2*u/b0)*b0/q)^(q/2)
    fg=f*gi
    minfg=min(minfg,min(fg)); maxfg=max(maxfg,max(fg))#G
    ind=(v<fg) 
    nacc=sum(ind); ntry=ntry+nleft; nleft=nleft-nacc; mloop=mloop+1#G
     if(nacc>0) values=rbind(values,x[ind,,drop=FALSE])#G
  }
  summ=c(ntry,(nsim-nleft)/ntry,(nleft==0),mloop=mloop,minfg,maxfg)#G
  names(summ)=c("ntry","eff","success","mloops","minfg","maxfg") #G
  attr(values,"summary")=summ 
  values
}

rvMsin.torus=function(nsim,k1,k2,alpha,mtop=1000) {
  X=rBessel(nsim,k1,k2,alpha)
  summ=attr(X,"summary")
  kappa=sqrt(k2^2+alpha*X[,2]) # conditional concentrations
  Y=matrix(0,nsim,2)
  for(i in 1:nsim) Y[i,]=rFisherBingham(1,c(k2,alpha*X[i,2]))
  values=cbind(X,Y); attr(values,"summary")=summ
  values
  
}

