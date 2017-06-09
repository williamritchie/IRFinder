gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

DESeqDataSetFromIRFinder <- function(filePaths,designMatrix,designFormula){
    res=c()
    libsz=c()
    irtest=read.table(filePaths[1])
    if (irtest[1,1]=="Chr"){irtest=irtest[-1,]}
    irnames=unname(apply(as.matrix(irtest),1,FUN=function(x){return(paste0(x[4],"/",x[1],":",x[2],"-",x[3],":",x[6]))}))
    n=1
    for (i in as.vector(files$V1)){
        print(paste0("processing file ",n," at ",i))
        irtab=read.table(i)
        if (irtab[1,1]=="Chr"){irtab=irtab[-1,]}
        #rn=unname(apply(irtab,1,FUN=function(x){return(paste0(x[4],"/",x[1],":",x[2],"-",x[3],":",x[6]))}))
        #row.names(irtab)=rn
        #tmp1=round(as.numeric(as.vector(irtab[irnames,9])))
        #tmp2=as.numeric(as.vector(irtab[irnames,19]))
        tmp1=as.numeric(as.vector(irtab[,9]))
        tmp2=as.numeric(as.vector(irtab[,19]))
        tmp3=tmp1+tmp2
        res=cbind(res,tmp1)
        libsz=cbind(libsz,tmp3)
        n=n+1
    }
    res.rd=round(res)
    libsz.rd=round(libsz)
    colnames(res)=as.vector(designMatrix[,1])
    rownames(res)=irnames
    colnames(libsz)=as.vector(designMatrix[,1])
    rownames(libsz)=irnames
    colnames(libsz.rd)=as.vector(designMatrix[,1])
    rownames(libsz.rd)=irnames
    dd = DESeqDataSetFromMatrix(countData = res.rd, colData = designMatrix, design = designFormula)
    colnames(dd)=as.vector(designMatrix[,1])
    rownames(dd)=irnames
    gm=apply(libsz.rd,1,gm_mean)
    normFactors <- libsz / gm
    normFactors[normFactors==0]=1
    normalizationFactors(dd) <- normFactors
    sp=libsz-res
    final=list(dd,res,sp)
    names(final)=c("DESeq2Object","IntronDepth","SpliceDepth")
    return(final)
}

