rocplot=function(data.null,data.alt){
  n.null=length(data.null)
  n.alt=length(data.alt)
  n=n.null+n.alt
  res=c(data.null,data.alt)
  frac=unique(sort(res))
  fp=0
  tp=0
  fn=0
  tn=0
  for(i in 1:length(frac)){
    fp[i]=0
    tp[i]=0
    fn[i]=0
    tn[i]=0
    
    dis=res[1:n.null]>frac[i]
    tn[i]=tn[i]+sum(dis==0)
    fp[i]=fp[i]+sum(dis==1)
     
    dis=res[(n.null+1):n]>frac[i]
    fn[i]=fn[i]+sum(dis==0)
    tp[i]=tp[i]+sum(dis==1)
     
  }
  tpr=tp/(tp+fn)
  fpr=fp/(fp+tn)
  fdr=fp/(fp+tp)
  fdr[(fp+tp)==0]=0
  td=fp+tp
  return(list(tpr=tpr,fpr=fpr,fdr=fdr,td=td,tp=tp,tn=tn,fp=fp,fn=fn))
}