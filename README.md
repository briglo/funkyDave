# funkyDave
Place where i put stuff for dave

##usage
```R
library(devtools)
github_install("https://github.com/briglo/funkyDave.git")
load('r_objects/181128_survivalObjects.rdata')
plotSurv(genes=sample(rownames(cdat),25),expdat=cdat,df=clindat)
```
having issues with below (what is the subtype column in tcga?)

## getting plotSurv objects
```R
library(cgdsr)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")
# Test the CGDS endpoint URL using a few simple API tests
test(mycgds)
# Get list of cancer studies at server
getCancerStudies(mycgds)[,1]


# Get cancer studies with breast cancer (by all means do other cancers...)
grep("brca", getCancerStudies(mycgds)[,1],value=T) 
mycancerstudy="brca_tcga"

#get datasets that are a part of your cancer study
grep('mrna',getCaseLists(mycgds,mycancerstudy)[,1],value=T) 
mycaselist="brca_tcga_rna_seq_v2_mrna"

#get specific datasets available from case list
getGeneticProfiles(mycgds,mycancerstudy)[,1] 
mygeneticprofile="brca_tcga_rna_seq_v2_mrna_median_Zscores"

#preparing clinical object
clin<-getClinicalData(mycgds,mycaselist)
clin$newOS_MONTHS<-ifelse(as.numeric(clin$OS_MONTHS)>60,60,as.numeric(clin$OS_MONTHS)) #trims to 5 years
clin$coded_status<-ifelse(clin$OS_STATUS=="DECEASED",1,0)
clin$newcoded_status<-ifelse(as.numeric(clin$OS_MONTHS)>60,0,clin$coded_status) # remakes ppl who died after 5 years alive at 5 years
#add elf5
elf<-getProfileData(mycgds,genes="ELF5",mygeneticprofile,mycaselist)
clin$ELF5<-elf[rownames(clin),"ELF5"]
saveRDS(clin,file="r_objects/190726_cgdsr_TGCA_clin.RDS") #done separate while pulling expression...

#Preparing expression object
load("r_objects/181005_epi_ano_monocle_tern_cca_POISSONgsvaHgorsmyPAM50.rdata") #or any  seurat object
gen<-toupper(rownames(epidat@data))
 hasEntrez <- clusterProfiler::bitr(gen,fromType = "SYMBOL", toType = c("ENTREZID","SYMBOL"), OrgDb = "org.Hs.eg.db")
genes<-split(hasEntrez$SYMBOL,rep(1:489,time=27))
dat<-lapply(genes, function(x) getProfileData(mycgds,genes=x,mygeneticprofile,mycaselist))
ids<-unique(unlist(lapply(dat,row.names)))
cdat_tcga<-data.frame(do.call(rbind,lapply(dat,function(x) t(x[ids,]))))
clin<-readRDS("r_objects/190726_cgdsr_TGCA_clin.RDS")
clindat_tcga<-clin[colnames(cdat_tcga),]
save(cdat_tcga,clindat_tcga,file="r_objects/190726_cgdsr_TGCA_matched.rdata")
```