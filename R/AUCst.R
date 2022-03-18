#' Polynomial
#'
#' @param sk
#' @param t
#' @param theta
#' @return
#' @export
#'
#'
fp.func3=function(sk,t,theta){
  AUC.st.mat=c()
  for (s.curr in sk)
  {
    tsvec=cbind(1, s.curr, t, s.curr^2, t^2, s.curr*t, t^3, s.curr^3, (t^2)*s.curr, t*(s.curr^2))
    teta=as.vector(as.matrix(tsvec)%*%theta)
    teta[teta>=100]=100
    teta[teta<=-100]=-100
    AUC.st=exp(teta)/(1+exp(teta))
    AUC.st.mat=rbind(AUC.st.mat, AUC.st)}
  return(AUC.st.mat)}

#' calculate the log-likelihood
#'
#' @param theta
#' @param sk
#' @param d.times
#' @param n1.h.mat
#' @param n2.h.mat
#' @param eta.sk.mat
#' @param sk.le.dh
#' @param Ewt.die
#' @return
#' @export
#'
#'
logL=function(theta,sk,d.times,n1.h.mat,n2.h.mat,eta.sk.mat,sk.le.dh,Ewt.die){
  AUC.seq.mat=fp.func3(sk, d.times,theta)
  logL.seq= Ewt.die*eta.sk.mat*sk.le.dh*(n1.h.mat*log(AUC.seq.mat)+n2.h.mat*log(1-AUC.seq.mat))
  logL=-sum(logL.seq)}

#' calclualte AUCst
#' @param da
#' @param da.long
#' @param sk
#' @param par0
#' @param tseq.eval
#' @param resample
#' @param nsap
#' @param seed
#' @return
#' @export
#'
main1.sub.func=function(da, da.long, sk, par0=c(0.5,0,0,0,0,0,0,0,0,0), tseq.eval, resample=0, nsap=1,seed=12345){
  s.num=length(sk)
  t.num=length(tseq.eval)
  resap.AUC.mat=c()
  n=dim(da)[1]

  long.basic=data.frame(cbind(rep(da$obs_id,rep(s.num,n)),rep(sk,n)))
  names(long.basic)=c("obs_id","vtime")
  marker0=merge(long.basic,da.long,by=c("obs_id","vtime"),all=T)
  marker0$eta=!is.na(marker0$Zt)*1
  marker0=merge(marker0,da,by="obs_id")

  set.seed(seed)
  for( sap in 1:nsap){

    if (resample==1){ EXPwt=rexp(n, 1)}
    if (resample==0){ EXPwt=rep(1, n) }

    d.times=sort(unique(da$Y[da$delta==1]))
    md=length(d.times)
    Ymat=matrix(rep(da$Y,md),n)
    dt.mat=matrix(rep(d.times,rep(n,md)),n)
    dt.mat.tb=dt.mat+0.00001;
    n.h=apply((Ymat>dt.mat.tb),2,sum)
    death.mat=(Ymat==dt.mat)*da$delta
    marker0$eta[marker0$eta==0]=0
    RS.data=marker0

    RS.seq.list=RS.mat.list=RS.d.list=RS.dmat.list=n1.h.list=n2.h.list=list()
    Wn1.h.list=Wn2.h.list=list()
    n1n2.wt=matrix(rep(EXPwt,md),n)
    for (i in 1:s.num)
    {
      RS.seq.list[[i]]=RS.data[RS.data$vtime==sk[i],]
      RS.mat.list[[i]]<-matrix(rep(RS.seq.list[[i]]$Zt,md),n)
      RS.d.list[[i]]=apply(RS.mat.list[[i]]*death.mat,2,sum,na.rm=TRUE)
      RS.d.list[[i]][RS.d.list[[i]]==0]<-NA
      RS.dmat.list[[i]]=matrix(rep(RS.d.list[[i]],rep(n,md)),n)
      n1.h.list[[i]]=apply((RS.dmat.list[[i]]>RS.mat.list[[i]])*(Ymat>dt.mat.tb),2,sum,na.rm=TRUE)
      n2.h.list[[i]]=apply((RS.dmat.list[[i]]<=RS.mat.list[[i]])*(Ymat>dt.mat.tb),2,sum,na.rm=TRUE)

      Wn1.h.list[[i]]=apply((RS.dmat.list[[i]]>RS.mat.list[[i]])*(Ymat>dt.mat.tb)*n1n2.wt,2,sum,na.rm=TRUE)
      Wn2.h.list[[i]]=apply((RS.dmat.list[[i]]<=RS.mat.list[[i]])*(Ymat>dt.mat.tb)*n1n2.wt,2,sum,na.rm=TRUE)
    }

    n1.h.mat=n2.h.mat=c()
    Wn1.h.mat=Wn2.h.mat=c()
    eta.sk.mat=c()
    for (i in 1:s.num)
    {
      n1.h.mat=rbind(n1.h.mat,n1.h.list[[i]])
      n2.h.mat=rbind(n2.h.mat,n2.h.list[[i]])
      Wn1.h.mat=rbind(Wn1.h.mat, Wn1.h.list[[i]])
      Wn2.h.mat=rbind(Wn2.h.mat, Wn2.h.list[[i]])
      eta.sk.mat=rbind(eta.sk.mat, RS.d.list[[i]])
    }
    eta.sk.mat=(eta.sk.mat/eta.sk.mat)
    eta.sk.mat[is.na(eta.sk.mat)]=0
    sk.mat=matrix(rep(sk,each=md), ncol=md, byrow=TRUE)
    d.times.mat=matrix(rep(d.times,rep(s.num,md)),s.num)
    sk.le.dh=(sk.mat<=d.times.mat)

    da0<-cbind(da$Y,da$delta,EXPwt)
    da0.sort=da0[order(da0[,1], decreasing = F),]
    Ewt.dseq=da0.sort[da0.sort[,2]==1,3]
    Ewt.die=matrix(rep(Ewt.dseq,rep(s.num,md)),s.num)

    resap.fit=optim(par=par0, fn=logL, d.times=d.times, n1.h.mat=Wn1.h.mat, n2.h.mat=Wn2.h.mat, sk=sk,
                    eta.sk.mat=eta.sk.mat, sk.le.dh=sk.le.dh, control = list(maxit = 10000,reltol=1e-8),Ewt.die=Ewt.die)

    s.mat=matrix(rep(sk,each=t.num), ncol=t.num, byrow=TRUE)
    t.mat=matrix(rep(tseq.eval,rep(s.num,t.num)),s.num)
    slet.mat=1*(s.mat<=t.mat)
    slet.mat[slet.mat==0]=NA
    resap.AUC0=fp.func3(sk,tseq.eval,resap.fit$par)*slet.mat
    resap.AUC.vec= c(as.vector(t(resap.AUC0)), resap.fit$convergence)
    resap.AUC.mat=rbind(resap.AUC.mat, resap.AUC.vec)
    print(sap)
  }


  resap.AUC.mat=as.data.frame(resap.AUC.mat)
  names(resap.AUC.mat)[ncol(resap.AUC.mat)]<-"cov1"
  resap.AUC.mat$minAUC=apply(resap.AUC.mat[,1:(t.num*s.num)], 1, FUN=min, na.rm=T)
  resap.AUC.mat$maxAUC=apply(resap.AUC.mat[,1:(t.num*s.num)], 1, FUN=max, na.rm=T)
  resap.AUC.mat$cov2=ifelse(resap.AUC.mat$minAUC<0.00001,1,0)
  resap.AUC.mat$cov3=ifelse(resap.AUC.mat$maxAUC>0.99999,1,0)
  resap.AUC.mat.cov=as.matrix(resap.AUC.mat[resap.AUC.mat$cov1==0&resap.AUC.mat$cov2==0&resap.AUC.mat$cov3==0,])
  RESAP.sd=apply(resap.AUC.mat.cov[,1:(t.num*s.num)], 2, FUN=sd)
  RESAP.sd.mat=t(matrix(RESAP.sd, t.num, s.num))

fit.eta=optim(par=par0, fn=logL, d.times=d.times, n1.h.mat=n1.h.mat, n2.h.mat=n2.h.mat, sk=sk,
                eta.sk.mat=eta.sk.mat, sk.le.dh=sk.le.dh, control = list(maxit = 10000,reltol=1e-8),Ewt.die=1)

  theta.est=fit.eta$par
  conv.est=fit.eta$convergence
  AUC.est=fp.func3(sk,tseq.eval,theta.est)*slet.mat
  return(list(conv=conv.est, AUC=AUC.est, theta=theta.est, RS= RS.seq.list, SE=RESAP.sd.mat))
}

