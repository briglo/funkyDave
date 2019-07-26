#' plotSurv
#'
#' takes bioportal inputs and tests outcome relative to an elf5 split
#'
#' @param genes vector of gene symbols can be mouse or human, toupper'd to "convert" to human
#' @param df an output of cgdsr::getClinicalData with some modified fields
#' @param expdat an output of cgdsr::getProfileData with same ID as df
#' @param trimFirst logical, perform a hard trim of data first, default=F
#'
#' @return a plot and invisible list of stuff
#'
#' @examples
#' load("r_objects/181206_survivalObjects.rdata")
#' load("r_objects/180928_HALLMARKandGORSMYgeneListsforGSVA.rdata")
#' plotSurv(ngs[[sample(1:length(ngs),1)]],clindat,cdat)
#' 
#' @export
plotSurv<-function(genes,df=clindat_tcga,expdat=cdat_tcga,trimFirst=F){
require(survminer)
require(survival)
require(cowplot)
df$ELF5<-as.numeric(expdat['ELF5',rownames(df)])
if (trimFirst) {
    cuts<-quantile(df$ELF5,c(.33,.5,.67),na.rm=T)
df$binnedElf5<-as.vector(ifelse(df$ELF5<cuts[1],"low",ifelse(df$ELF5>cuts[3],"high","none")))
df<-df[!df$CLAUDIN_SUBTYPE %in% c("NC",""),]
df<-df[!is.na(df$binnedElf5),]
df<-df[!df$binnedElf5 %in% c("none","NA"),]
expdat<-expdat[,rownames(df)]
} ###up to here
 genelist<-toupper(genes[toupper(genes) %in% rownames(expdat)])
    mscore<-colSums(expdat[genelist,rownames(df)],na.rm=T)
	mq<-quantile(mscore,c(.25,.33,.5,.67,.75))
    mb<-ifelse(mscore<mq[2],'low',ifelse(mscore>mq[4],"high","mid"))
    elfExp<-t(expdat['ELF5',rownames(df)])
    eq<-quantile(elfExp,c(.25,.33,.5,.67),na.rm=T)
    eb<-ifelse(elfExp<eq[2],'low',ifelse(elfExp>eq[4],"high","mid"))
comDF<-data.frame(df,elfExp=as.vector(elfExp),binnedElf5= as.vector(eb), metascore=mscore,metabin=mb)[!is.na(as.vector(eb)),]


makeScores<-list('elf'=comDF[comDF$binnedElf5 %in% c("high",'low'),],
        'met'=comDF[comDF$metabin %in% c("high",'low'),],
		'split'=lapply(split(comDF[comDF$binnedElf5 %in% c("high",'low'),],comDF$binnedElf5[comDF$binnedElf5 %in% c("high",'low')]),function(x) {
		qs<-quantile(x$metascore,c(.25,.33,.5,.67,.75))
		binMetascore<-ifelse(x$metascore<qs[2],'low',ifelse(x$metascore>qs[4],"high","mid"))
		return(data.frame(x,binScore=binMetascore))}))


makeFits<-list('combElf'=survfit(Surv(newOS_MONTHS,newcoded_status) ~ binnedElf5, data=makeScores$elf),
'combMeta'=survfit(Surv(newOS_MONTHS,newcoded_status) ~ metabin, data=makeScores$met),

	'spl'=lapply(makeScores$split[1:2], function(X) survfit(Surv(newOS_MONTHS,newcoded_status) ~ binScore, data=X)))



 p1a<-ggsurvplot(makeFits$combElf,data=makeScores$elf, pval=T,conf.int=T) + ggtitle("elf5 bins")
	p1b<-ggsurvplot(makeFits$combMeta,data=makeScores$met, pval=T,conf.int=T) + ggtitle("metascore bins")
	p2a<-ggplot(makeScores$elf,aes(x=binnedElf5,fill=CLAUDIN_SUBTYPE)) + geom_bar(stat='count') + ggtitle("elf5 composition")
	p2b<-ggplot(makeScores$met,aes(x=metabin,fill=CLAUDIN_SUBTYPE)) + geom_bar(stat='count')+ ggtitle("metascore composition")
	p3<-ggsurvplot(makeFits$spl,data=makeScores$split[1:2], pval=T,conf.int=T)
    p4<-lapply(lapply(makeScores[[3]], function(y) data.frame(table(y$binScore,y$CLAUDIN_SUBTYPE))), function(z) ggplot(z,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat='identity'))
	p5<-ggplot(comDF,aes(x=elfExp,y=metascore,color=as.factor(newcoded_status))) + geom_point() + geom_smooth(method = "lm", se=T)
	bottom<-plot_grid(p5,ncol=2)
	top<-plot_grid(p1a$plot,p2a,p1b$plot,p2b,p3$high$plot,p4$high,p3$low$plot,p4$low,ncol=4)
	print(plot_grid(top,bottom,ncol=1))
	return(invisible(comDF))
}
