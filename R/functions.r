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



#' intgraph
#'
#' takes a summarized cellphoneDB object and plots it better than my other attempts
#'
#' @param interactionData a dataframe with source target and countAboveMean[numeric]
#' @param annotation a dataframe of fields name, type, ID, count, optional x+y coordinates for manual graph layout
#' @param scoreCut a numeric value pointing to a column of interactionData (default 0.3)
#' @param numberCut numeric value for minimum number of interactions to plot (default 20)
#' @param numberSplit numeric value for spliting numbers of interactions (i.e high vs low, default 60)
#'
#' @return a plot, invisibly the graph
#'
#' @examples
#' data("190801_newNetworkGraph.rdata")
#' IG<-intgraph(interactionData=cellphoneDB_data,annotation=cluster_anno,scoreCut=0.3, numberCut=20, numberSplit=60)
#' 
#' 
#' @export
intgraph<- function(interactionData,annotation,scoreCut=0.3, numberCut=20, numberSplit=60){
   require(ggraph)
   require(dplyr)
   require(igraph)
   require(cowplot)
   message("use like this: intgraph(scoreCut=0.3, numberCut=0, numberSplit=35)")
 tmp<-cellphoneDB_data[cellphoneDB_data[,paste0("countAboveMean",scoreCut)]>numberCut,c('source','target',paste0("countAboveMean",scoreCut))]
 colnames(tmp)[3]<-"score"
 gr<-graph_from_data_frame(tmp,directed = F,vertices=cluster_anno)
 deg=degree(gr, mode ="all")
 colname<-paste0("countAboveMean",numberCut)
   p1<-ggraph(gr, layout = 'linear', circular=T) +     geom_edge_arc(aes(width =score,alpha=score>numberSplit,colour=score>numberSplit)) + geom_node_point(aes(color=paste(name,ID),shape=type,size=count)) + geom_node_text(aes(label=paste(name,ID)), repel=F) + ggtitle(paste0("scoreCut=",scoreCut, ", numberCut=",numberCut,", numberSplit=",numberSplit))
   p2<-ggplot(tmp,aes(x=score,colour=score<numberSplit)) + geom_bar() + ggtitle("score dist")
  print( ggdraw() + draw_plot(p1 + theme(legend.justification="bottom"), 0,0,1,1) + draw_plot(p2  + theme(legend.position = "none") ,.85, 0, 0.15  , 0.15))
   return(invisible(gr))
   }


#' TSNEintgraph needs to be tested
#'
#' building on intgraph, adds network to a TSNE plor
#'
#' @param intGraph an output from intgraph
#' @param seuratObj a Seurat object
#' @param metaCol character, name of metadata column corresponding to nodes of intGraph
#' @param metaColPlot name of metadata column corresponding to nodes of intGraph if different from metaCol, default NULL
#' @param numberCut numeric value for minimum number of interactions to plot (default 20)
#' @param numberSplit numeric value for spliting numbers of interactions (i.e high vs low, default 60)
#'
#' @return a plot 
#'
#' @examples
#' data("190801_newNetworkGraph.rdata")
#' IG<-intgraph(interactionData=cellphoneDB_data,annotation=cluster_anno,scoreCut=0.3, numberCut=20, numberSplit=60)
#' load("PATH/TO/SEURAT/OBJ) 
#' TSNEintgraph(intGraph=IG,seuratObj=TA,metaCol="cellphoneDB_id",metaColPlot="cellphoneDB_id",numberCut=20,numberSplit=60)
#' 
#' 
#' @export
TSNEintgraph<-function(intGraph,seuratObj,metaCol,metaColPlot=NULL,numberCut,numberSplit){
     locs<-data.frame(do.call(rbind,lapply(split(data.frame(seuratObj@dr$tsne@cell.embeddings),seuratObj@meta.data[,metaCol]),colMeans)))
     colnames(locs)<-c('x','y')
     newAnno<-cbind(locs,as.data.frame(get.vertex.attribute(intGraph)))
p1<- ggraph(intGraph,layout = 'manual',node.position=newAnno) +     geom_edge_arc(aes(width =score,alpha=score>numberSplit,colour=score)) + geom_node_point(aes(shape=type,size=count)) + geom_node_text(aes(label=paste(name,ID)), repel=T,size=6) + theme(legend.position = "none") + xlim(-50, 50) + ylim(-50, 50) 

if(is.null(metaColPlot)) { p2<-ggplot(cbind(seuratObj@meta.data,seuratObj@dr$tsne@cell.embeddings),aes(x=tSNE_1,y=tSNE_2,colour=get(metaCol))) + geom_point() + theme(legend.position = "none") + xlim(-50, 50) + ylim(-50, 50) } else { p2<-ggplot(cbind(seuratObj@meta.data,seuratObj@dr$tsne@cell.embeddings),aes(x=tSNE_1,y=tSNE_2,colour=get(metaColPlot))) + geom_point() + theme(legend.position = "none") + xlim(-50, 50) + ylim(-50, 50) }

print(ggdraw() + draw_plot(p2 + theme(legend.position = "none"),0,0,1,1) + draw_plot(p1,0,0,1,1))

}
