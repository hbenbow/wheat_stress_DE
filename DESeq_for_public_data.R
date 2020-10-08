library(DESeq2)
library(dplyr)
library(tidyr)

setwd("~/Documents/Sobia/orphans/Expression/")
list.files("Raw_expression_data/", pattern="_count") # list all raw expression files
files<-list.files("Raw_expression_data/", pattern="_count")
metadata<-read.delim("~/Documents/Sobia/orphans/Expression/Raw_expression_data/default_metadata.txt")
metadata$Factor<-paste(metadata$Variety, metadata$Age, metadata$Stress.disease)
{
  for(i in files){
    data<-read.table(paste("Raw_expression_data/", i, sep=""), header=T)
    name<-paste(i)
    study<-substr(name, 1, nchar(name)-13)
    meta<-subset(metadata, metadata$secondary_study_accession %in%  study)
    names<-data$transcript
    data[1]<-NULL
    data<-apply( data, 2, as.integer)
    rownames(data)<-names
    order<-as.data.frame(colnames(data))
    colnames(order)<-"run_accession"
    meta<-merge(order, meta, by="run_accession", sort=F)
    dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Factor)
    dds<-DESeq(dds)
    assign(paste(name, "meta"), meta)
    assign(paste(name, "dds"), dds)
  }
  
  # DESeq2 differential expression testing
  {RESDRP000768  <-as.data.frame(results(`DRP000768_count.tsv.gz dds`, contrast=c("Factor","Chinese Spring 24 days Phosphorous starvation 10 days", "Chinese Spring 24 days none"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
    RESDRP000768$Factor<-"ChineseSpring 24 days Phosphorous starvation 10 days"
    RESDRP000768$Gene<-row.names(RESDRP000768)
    RESDRP000768<-RESDRP000768[(RESDRP000768$padj <= 0.05),]
    RESDRP000768<-na.omit(RESDRP000768)}
  {
    R1<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 12 hours",  "NIL38 from BC55F2 Remus in CM-82036 anthesis mock inoculation 12 hours")))
    R2<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 24 hours",  "NIL38 from BC55F2 Remus in CM-82036 anthesis mock inoculation 24 hours")))
    R3<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 3 hours",  "NIL38 from BC55F2 Remus in CM-82036 anthesis mock inoculation 3 hours")))
    R4<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 36 hours",  "NIL38 from BC55F2 Remus in CM-82036 anthesis mock inoculation 36 hours")))
    R5<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 48 hours",  "NIL38 from BC55F2 Remus in CM-82036 anthesis mock inoculation 48 hours")))
    R6<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 6 hours",  "NIL38 from BC55F2 Remus in CM-82036 anthesis mock inoculation 6 hours")))
    R7<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 12 hours",  "NIL51 from BC55F2 Remus in CM-82036 anthesis mock inoculation 12 hours")))
    R8<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 24 hours",  "NIL51 from BC55F2 Remus in CM-82036 anthesis mock inoculation 24 hours")))
    R9<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 3 hours",  "NIL51 from BC55F2 Remus in CM-82036 anthesis mock inoculation 3 hours")))
    R10<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 36 hours",  "NIL51 from BC55F2 Remus in CM-82036 anthesis mock inoculation 36 hours")))
    R11<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 48 hours",  "NIL51 from BC55F2 Remus in CM-82036 anthesis mock inoculation 48 hours")))
    R12<-as.data.frame(results(`ERP013829_count.tsv.gz dds`, contrast=c("Factor", "NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 6 hours",  "NIL51 from BC55F2 Remus in CM-82036 anthesis mock inoculation 6 hours")))
    
    R1$Factor<-"NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 12 hours"
    R2$Factor<-"NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 24 hours"
    R3$Factor<-"NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 3 hours"
    R4$Factor<-"NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 36 hours"
    R5$Factor<-"NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 48 hours"
    R6$Factor<-"NIL38 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 6 hours"
    R7$Factor<-"NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 12 hours"
    R8$Factor<-"NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 24 hours"
    R9$Factor<-"NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 3 hours"
    R10$Factor<-"NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 36 hours"
    R11$Factor<-"NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 48 hours"
    R12$Factor<-"NIL51 from BC55F2 Remus in CM-82036 anthesis Fusarium graminearum inoculation 6 hours" 
    R1$Gene<-row.names(R1)
    R2$Gene<-row.names(R2)
    R3$Gene<-row.names(R3)
    R4$Gene<-row.names(R4)
    R5$Gene<-row.names(R5)
    R6$Gene<-row.names(R6)
    R7$Gene<-row.names(R7)
    R8$Gene<-row.names(R8)
    R9$Gene<-row.names(R9)
    R10$Gene<-row.names(R10)
    R11$Gene<-row.names(R11)
    R12$Gene<-row.names(R12)
    RESERP013829<-rbind(R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12)
    RESERP013829<- RESERP013829[( RESERP013829$padj <= 0.05),]
    RESERP013829<-na.omit( RESERP013829)
  }
  {R1<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c("Factor", "CM-82036 anthesis Fusarium graminearum inoculation 30 hours" , "CM-82036 anthesis mock inoculation 30 hours")))
    R1$Factor<-"CM-82036 anthesis Fusarium graminearum inoculation 30 hours"
    R2<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c("Factor", "CM-82036 anthesis Fusarium graminearum inoculation 50 hours" , "CM-82036 anthesis mock inoculation 50 hours")))
    R2$Factor<-"CM-82036 anthesis Fusarium graminearum inoculation 50 hours"
    R3<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c("Factor", "NIL1 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 30 hours", "NIL1 from CM-82036 BC5 in Remus anthesis mock inoculation 30 hours")))
    R3$Factor<-"NIL1 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 30 hours"
    R4<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c("Factor", "NIL1 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 50 hours", "NIL1 from CM-82036 BC5 in Remus anthesis mock inoculation 50 hours" )))
    R4$Factor<-"NIL1 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 50 hours"
    R5<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c( "Factor","NIL2 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 30 hours", "NIL2 from CM-82036 BC5 in Remus anthesis mock inoculation 30 hours" )))
    R5$Factor<-"NIL2 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 30 hours"
    R6<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c( "Factor","NIL2 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 50 hours", "NIL2 from CM-82036 BC5 in Remus anthesis mock inoculation 50 hours")))
    R6$Factor<-"NIL2 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 50 hours"
    R7<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c( "Factor","NIL3 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 30 hours", "NIL3 from CM-82036 BC5 in Remus anthesis mock inoculation 30 hours" )))
    R7$Factor<-"NIL3 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 30 hours"
    R8<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c("Factor", "NIL3 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 50 hours", "NIL3 from CM-82036 BC5 in Remus anthesis mock inoculation 50 hours")))
    R8$Factor<-"NIL3 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 50 hours"
    R9<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c( "Factor","NIL4 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 30 hours", "NIL4 from CM-82036 BC5 in Remus anthesis mock inoculation 30 hours")))
    R9$Factor<-"NIL4 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 30 hours"
    R10<-as.data.frame(results(`ERP003465_count.tsv.gz dds`, contrast=c("Factor", "NIL4 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 50 hours", "NIL4 from CM-82036 BC5 in Remus anthesis mock inoculation 50 hours")))
    R10$Factor<-"NIL4 from CM-82036 BC5 in Remus anthesis Fusarium graminearum inoculation 50 hours"
    R1$Gene<-row.names(R1)
    R2$Gene<-row.names(R2)
    R3$Gene<-row.names(R3)
    R4$Gene<-row.names(R4)
    R5$Gene<-row.names(R5)
    R6$Gene<-row.names(R6)
    R7$Gene<-row.names(R7)
    R8$Gene<-row.names(R8)
    R9$Gene<-row.names(R9)
    R10$Gene<-row.names(R10)
    RESERP003465<-rbind(R1, R2, R3, R4, R5, R6, R7, R8, R9, R10)
    RESERP003465<-RESERP003465[(RESERP003465$padj <= 0.05),]
    RESERP003465<-na.omit(RESERP003465)}
  {R1<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Vuka three leaf stage stripe rust pathogen 87/66 1 day", "Vuka three leaf stage control")))
    R2<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Vuka three leaf stage stripe rust pathogen 87/66 2 days", "Vuka three leaf stage control")))
    R3<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Vuka three leaf stage stripe rust pathogen 87/66 3 days", "Vuka three leaf stage control")))
    R4<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Vuka three leaf stage stripe rust pathogen 87/66 5 days", "Vuka three leaf stage control")))
    R5<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Vuka three leaf stage stripe rust pathogen 87/66 7 days", "Vuka three leaf stage control")))
    R6<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Vuka three leaf stage stripe rust pathogen 87/66 9 days", "Vuka three leaf stage control")))
    R7<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Vuka three leaf stage stripe rust pathogen 87/66 11 days", "Vuka three leaf stage control")))
    R8<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Avocet+Yr5 three leaf stage stripe rust pathogen 87/66 1 day", "Avocet+Yr5 three leaf stage control")))
    R9<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Avocet+Yr5 three leaf stage stripe rust pathogen 87/66 2 days", "Avocet+Yr5 three leaf stage control")))
    R10<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Avocet+Yr5 three leaf stage stripe rust pathogen 87/66 3 days", "Avocet+Yr5 three leaf stage control")))
    R11<-as.data.frame(results(`ERP013983_count.tsv.gz dds`, contrast=c("Factor","Avocet+Yr5 three leaf stage stripe rust pathogen 87/66 5 days", "Avocet+Yr5 three leaf stage control")))
    R1$Factor<- "Vuka three leaf stage stripe rust pathogen 87/66 1 day"       
    R2$Factor<- "Vuka three leaf stage stripe rust pathogen 87/66 2 days"      
    R3$Factor<- "Vuka three leaf stage stripe rust pathogen 87/66 3 days"      
    R4$Factor<- "Vuka three leaf stage stripe rust pathogen 87/66 5 days"      
    R5$Factor<- "Vuka three leaf stage stripe rust pathogen 87/66 7 days"      
    R6$Factor<- "Vuka three leaf stage stripe rust pathogen 87/66 9 days"      
    R7$Factor<- "Vuka three leaf stage stripe rust pathogen 87/66 11 days"     
    R8$Factor<- "Avocet+Yr5 three leaf stage stripe rust pathogen 87/66 1 day" 
    R9$Factor<- "Avocet+Yr5 three leaf stage stripe rust pathogen 87/66 2 days"
    R10$Factor<- "Avocet+Yr5 three leaf stage stripe rust pathogen 87/66 3 days"
    R11$Factor<- "Avocet+Yr5 three leaf stage stripe rust pathogen 87/66 5 days"
    R1$Gene<-row.names(R1)
    R2$Gene<-row.names(R2)
    R3$Gene<-row.names(R3)
    R4$Gene<-row.names(R4)
    R5$Gene<-row.names(R5)
    R6$Gene<-row.names(R6)
    R7$Gene<-row.names(R7)
    R8$Gene<-row.names(R8)
    R9$Gene<-row.names(R9)
    R10$Gene<-row.names(R10)
    R11$Gene<-row.names(R11)
    RESERP013983<-rbind(R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11)
    RESERP013983<-RESERP013983[(RESERP013983$padj <= 0.05),]
    RESERP013983<-na.omit(RESERP013983)}
  {RESERP015130<-as.data.frame(results(`ERP015130_count.tsv.gz dds`, contrast=c("unknown (Bangladesh) grain filling Magnaporthe oryzae symptomatic", "unknown (Bangladesh) grain filling Magnaporthe oryzae asymptomatic")))
    RESERP015130$Factor<-"unknown (Bangladesh) grain filling Magnaporthe oryzae"
    RESERP015130$Gene<-row.names( RESERP015130)
    RESERP015130<-RESERP0151303[(RESERP015130$padj <= 0.05),]
    RESERP015130<-na.omit(RESERP015130)}
  {R1<-as.data.frame(results(`SRP022869_count.tsv.gz dds`, contrast=c("Factor", "Sevin 14 days Septoria tritici 4 days" , "Sevin 14 days none")))
    R1$Factor<-"Sevin 14 days Septoria tritici 4 days"
    R2<-as.data.frame(results(`SRP022869_count.tsv.gz dds`, contrast=c("Factor", "Sevin 14 days Septoria tritici 10 days" , "Sevin 14 days none")))
    R2$Factor<-"Sevin 14 days Septoria tritici 10 days"
    R3<-as.data.frame(results(`SRP022869_count.tsv.gz dds`, contrast=c( "Factor","Sevin 14 days Septoria tritici 13 days", "Sevin 14 days none")))
    R3$Factor<-"Sevin 14 days Septoria tritici 13 days"
    RESSRP022869<-rbind(R1, R2, R3)
    R1$Gene<-row.names(R1)
    R2$Gene<-row.names(R2)
    R3$Gene<-row.names(R3)
    RESSRP022869<-RESSRP022869[(RESSRP022869$padj <= 0.05),]
    RESSRP022869<-na.omit(RESSRP022869)}
  {R1<-as.data.frame(results(`SRP041017_count.tsv.gz dds`, contrast=c( "Factor","N9134 7 days Powdery mildew pathogen E09 72 hours", "N9134 7 days none" )))
    R2<-as.data.frame(results(`SRP041017_count.tsv.gz dds`, contrast=c("Factor", "N9134 7 days Powdery mildew pathogen E09 48 hours", "N9134 7 days none")))
    R3<-as.data.frame(results(`SRP041017_count.tsv.gz dds`, contrast=c( "Factor","N9134 7 days Powdery mildew pathogen E09 24 hours", "N9134 7 days none")))
    R4<-as.data.frame(results(`SRP041017_count.tsv.gz dds`, contrast=c("Factor", "N9134 7 days stripe rust pathogen CYR31 72 hours" , "N9134 7 days none")))
    R5<-as.data.frame(results(`SRP041017_count.tsv.gz dds`, contrast=c("Factor", "N9134 7 days stripe rust pathogen CYR31 48 hours" , "N9134 7 days none")))
    R6<-as.data.frame(results(`SRP041017_count.tsv.gz dds`, contrast=c("Factor", "N9134 7 days stripe rust pathogen CYR31 24 hours" , "N9134 7 days none")))
    R1$Factor<- "N9134 7 days Powdery mildew pathogen E09 72 hours"
    R2$Factor<- "N9134 7 days Powdery mildew pathogen E09 48 hours"
    R3$Factor<- "N9134 7 days Powdery mildew pathogen E09 24 hours"
    R4$Factor<- "N9134 7 days stripe rust pathogen CYR31 72 hours" 
    R5$Factor<- "N9134 7 days stripe rust pathogen CYR31 48 hours" 
    R6$Factor<- "N9134 7 days stripe rust pathogen CYR31 24 hours" 
    R1$Gene<-row.names(R1)
    R2$Gene<-row.names(R2)
    R3$Gene<-row.names(R3)
    R4$Gene<-row.names(R4)
    R5$Gene<-row.names(R5)
    R6$Gene<-row.names(R6)
    RESSRP041017<-rbind(R1, R2, R3, R4, R5, R6)
    RESSRP041017<-RESSRP041017[(RESSRP041017$padj <= 0.05),]
    RESSRP041017<-na.omit(RESSRP041017)}
  {R1<-as.data.frame(results(`SRP045409_count.tsv.gz dds`, contrast=c("Factor", "TAM107 7 days 6 hour of drought&heat combined stress", "TAM107 7 days none"  )))
    R2<-as.data.frame(results(`SRP045409_count.tsv.gz dds`, contrast=c( "Factor","TAM107 7 days 1 hour of drought&heat combined stress", "TAM107 7 days none"  )))
    R3<-as.data.frame(results(`SRP045409_count.tsv.gz dds`, contrast=c("Factor", "TAM107 7 days 6 hour of heat stress"                 , "TAM107 7 days none"  )))
    R4<-as.data.frame(results(`SRP045409_count.tsv.gz dds`, contrast=c("Factor", "TAM107 7 days 1 hour of heat stress"                 , "TAM107 7 days none"  )))
    R5<-as.data.frame(results(`SRP045409_count.tsv.gz dds`, contrast=c("Factor", "TAM107 7 days 6 hour of drought stress"              , "TAM107 7 days none"  )))
    R6<-as.data.frame(results(`SRP045409_count.tsv.gz dds`, contrast=c("Factor", "TAM107 7 days 1 hour of drought stress"              , "TAM107 7 days none"  )))
    R1$Factor<- "TAM107 7 days 6 hour of drought&heat combined stress"
    R2$Factor<- "TAM107 7 days 1 hour of drought&heat combined stress"
    R3$Factor<- "TAM107 7 days 6 hour of heat stress"                 
    R4$Factor<- "TAM107 7 days 1 hour of heat stress"                 
    R5$Factor<- "TAM107 7 days 6 hour of drought stress"              
    R6$Factor<- "TAM107 7 days 1 hour of drought stress"              
    R1$Gene<-row.names(R1)
    R2$Gene<-row.names(R2)
    R3$Gene<-row.names(R3)
    R4$Gene<-row.names(R4)
    R5$Gene<-row.names(R5)
    R6$Gene<-row.names(R6)
    RESSRP045409<-rbind(R1, R2, R3, R4, R5, R6)
    RESSRP045409<-RESSRP045409[(RESSRP045409$padj <= 0.05),]
    RESSRP045409<-na.omit(RESSRP045409)}
  {R1<-as.data.frame(results(`SRP048912_count.tsv.gz dds`, contrast=c("Factor", "Janz*2 NIL1 susceptible 2 days Fusarium pseudograminearum inoculation 5 days", "Janz*2 NIL1 susceptible 2 days mock inoculation 5 days"                      )))
    R2<-as.data.frame(results(`SRP048912_count.tsv.gz dds`, contrast=c("Factor", "Janz*2 NIL1 resistant 2 days Fusarium pseudograminearum inoculation 5 days"  , "Janz*2 NIL1 resistant 2 days mock inoculation 5 days"                        )))
    R3<-as.data.frame(results(`SRP048912_count.tsv.gz dds`, contrast=c( "Factor","Janz*2 NIL1 susceptible 2 days Fusarium pseudograminearum inoculation 3 days", "Janz*2 NIL1 susceptible 2 days mock inoculation 3 days"                      )))
    R4<-as.data.frame(results(`SRP048912_count.tsv.gz dds`, contrast=c( "Factor","Janz*2 NIL1 resistant 2 days Fusarium pseudograminearum inoculation 3 days"  , "Janz*2 NIL1 resistant 2 days mock inoculation 3 days"    )))
    R1$Factor<- "Janz*2 NIL1 susceptible 2 days Fusarium pseudograminearum inoculation 5 days"
    R2$Factor<- "Janz*2 NIL1 resistant 2 days Fusarium pseudograminearum inoculation 5 days"  
    R3$Factor<- "Janz*2 NIL1 susceptible 2 days Fusarium pseudograminearum inoculation 3 days"
    R4$Factor<- "Janz*2 NIL1 resistant 2 days Fusarium pseudograminearum inoculation 3 days"  
    R1$Gene<-row.names(R1)
    R2$Gene<-row.names(R2)
    R3$Gene<-row.names(R3)
    R4$Gene<-row.names(R4)
    RESSRP048912<-rbind(R1, R2, R3, R4)
    RESSRP048912<-RESSRP048912[(RESSRP048912$padj <= 0.05),]
    RESSRP048912<-na.omit(RESSRP048912)}
  {RESSRP060670<-as.data.frame(results(`SRP060670_count.tsv.gz dds`, contrast=c("Factor","Chinese Spring anthesis Fusarium graminearum inoculation 4 days", "Chinese Spring anthesis mock inoculation 4 days"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
    RESSRP060670$Factor<-"Chinese Spring anthesis Fusarium graminearum inoculation 4 days"
    RESSRP060670$Gene<-row.names(RESSRP060670)
    RESSRP060670<-RESSRP060670[(RESSRP060670$padj <= 0.05),]
    RESSRP060670<-na.omit(RESSRP060670)}
  {R1<-as.data.frame(results(`SRP068165_count.tsv.gz dds`, contrast=c("Factor",  "Giza 168 9 days PEG 6000 12 hours"  , "Giza 168 9 days none")))
    R2<-as.data.frame(results(`SRP068165_count.tsv.gz dds`, contrast=c("Factor",  "Giza 168 9 days PEG 6000 2 hours"   , "Giza 168 9 days none")))
    R3<-as.data.frame(results(`SRP068165_count.tsv.gz dds`, contrast=c( "Factor", "Gemmiza 10 9 days PEG 6000 12 hours", "Gemmiza 10 9 days none")))
    R4<-as.data.frame(results(`SRP068165_count.tsv.gz dds`, contrast=c( "Factor", "Gemmiza 10 9 days PEG 6000 2 hours" , "Gemmiza 10 9 days none")))
    R1$Factor<- "Giza 168 9 days PEG 6000 12 hours"  
    R2$Factor<- "Giza 168 9 days PEG 6000 2 hours"   
    R3$Factor<- "Gemmiza 10 9 days PEG 6000 12 hours"
    R4$Factor<- "Gemmiza 10 9 days PEG 6000 2 hours" 
    R1$Gene<-row.names(R1)
    R2$Gene<-row.names(R2)
    R3$Gene<-row.names(R3)
    R4$Gene<-row.names(R4)
    RESSRP068165<-rbind(R1, R2, R3, R4)
    RESSRP068165<-RESSRP068165[(RESSRP068165$padj <= 0.05),]
    RESSRP068165<-na.omit(RESSRP068165)}
  {RESSRP078208<-as.data.frame(results(`SRP078208_count.tsv.gz dds`, contrast=c("Factor","Chara 3 days  Fusarium pseudograminearum inoculation 36 hours", "Chara 3 days  mock inoculation 36 hours"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
    RESSRP078208$Factor<-"Chara 3 days  Fusarium pseudograminearum inoculation 36 hours"
    RESSRP078208$Gene<-row.names(RESSRP078208)
    RESSRP078208<-RESSRP078208[(RESSRP078208$padj <= 0.05),]
    RESSRP078208<-na.omit(RESSRP078208)}
}



dir.create("~/Documents/Sobia/orphans/Expression/DEgenes")
setwd("~/Documents/Sobia/orphans/Expression/DEgenes")
write.csv(RESERP013829,file= "RES_ERP013983.csv", row.names=F)
write.csv(RESDRP000768,file= "RES_DRP000768.csv", row.names=F)
write.csv(RESERP003465,file= "RES_ERP003465.csv", row.names=F)
write.csv(RESSRP022869,file= "RES_SRP022869.csv", row.names=F)
write.csv(RESSRP041017,file= "RES_SRP041017.csv", row.names=F)
write.csv(RESSRP045409,file= "RES_SRP045409.csv", row.names=F)
write.csv(RESSRP048912,file= "RES_SRP048912.csv", row.names=F)
write.csv(RESSRP060670,file= "RES_SRP060670.csv", row.names=F) 
write.csv(RESSRP068165,file= "RES_SRP068165.csv", row.names=F) 
write.csv(RESSRP078208,file= "RES_SRP078208.csv", row.names=F)
write.csv(all_data, file='~/Documents/Sobia/orphans/Expression/DEgenes/all_data.csv')
