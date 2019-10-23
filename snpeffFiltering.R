## v3.2.1
# good for snpEff2human 1.1

## import libraries
library(plyr)


## set cutoff. varients present in 'cutoff' or more strains will be removed
cutoff<- 3



dir <- "Z:/path/to/snpeffs/"         ## include final '/'



finaldir <- "filtered/"  


#### preprocessing #####

## find files
files <- list.files(dir, pattern = ".tsv")
nfiles <- length(files)


## read in data
allCPRC <- NULL
for(i in 1:nfiles){
  print( paste0( c("Reading... ",files[i]),collapse="") )
  Ncomments <- max(which(substr(readLines(paste0(dir,files[i]),n=100),1,1)=="#")) 
  temp <- read.delim(paste0(dir,files[i]), header=T, skip=Ncomments-1)                        ## read in
  temp$CPRC <- apply(temp,1,function(x) paste(x[1],x[2],x[4],x[5]))                           ## add CPRC col
  allCPRC <- c(allCPRC,levels(factor(temp$CPRC)))                                             ## add these CPRCs to growing list of all CPRCs
  name <- paste0(c("x",i),collapse="")                                                        ## build variable name
  temp$worm <- name
  fname<- unlist(strsplit(files[i],"\\."))[1]
  ffname <- paste0(fname,"_cut",cutoff)
  CPRCcount <- data.frame(table(allCPRC))
  
  WORMSAMPLE<-strsplit(as.character(temp$WORMSAMPLE),":")
  temp$altPercent<-sapply(WORMSAMPLE,function(x) {
    as.numeric(strsplit(x[2],",")[[1]][2])/as.numeric(x[3])
    } )
  assign(paste0(name,"name"), ffname)
  assign(name, temp)                                                                          ## assign var
  
}
print( paste0( c("Read ",nfiles," files."),collapse="") )


#### Juicy bits ####### 
#### Run any of the code blocks you want, in any order.

## find unique variants 

CPRCcount <- data.frame(table(allCPRC))                                  ## counts times each CPRC shows up.
colnames(CPRCcount)[2] <- "numWorms"                                     ## rename  column

ParSummary <- data.frame(table(CPRCcount$numWorms))                      ## aggregate into summary of parental var distribution
colnames(ParSummary)[2] <- "count"
barplot(ParSummary$count, names.arg=ParSummary$numWorms, 
        xlab="# of worms", ylab="# of variants")                         ## plot it out

parental <- parentalRemaining <- CPRCcount[which(CPRCcount$numWorm>=cutoff),]          ## assemble parentals into a DF
parentalSnpEff <- cbind(temp,"numWorms"=0)[0,]                           ## empty dataframe but with all the col names ready
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  indexofparentalCPRC <- unique(match(temp$CPRC,as.character(parentalRemaining$allCPRC)))
  indexofparentalCPRC <- indexofparentalCPRC[!is.na(indexofparentalCPRC)]
  theseParentalSnpEffs <- temp[which(!is.na(match(temp$CPRC,as.character(parentalRemaining$allCPRC)))),]
  parentalSnpEff <- rbind(parentalSnpEff,
                          cbind(theseParentalSnpEffs,
                                numWorms=  parental$numWorms[match(theseParentalSnpEffs$CPRC,parental$allCPRC)]
                          )
                          )
  parentalRemaining <- parentalRemaining[-c(indexofparentalCPRC),]
  
  temp <- temp[which(is.na(match(temp$CPRC,as.character(parental$allCPRC)))),]             ## remove parentals from worm data
  assign(varName, temp)  
}
parentalSnpEff <- parentalSnpEff[order(as.numeric(as.roman(gsub("MtDNA","L",as.character(parentalSnpEff$X.CHROM)))),parentalSnpEff$POS),-c(6:24,47,48,50,51)]



unique(x1$CPRC)[which(unique(x1$CPRC) %in% unique(x2$CPRC))]

## Extract homozygous mutations
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(temp$AC=="2"),]
  
  assign(varName, temp)  
}

# ## Remove upstream and downstream mutations
# for( i in 1:nfiles){
#   varName <- paste0("x",i)
#   temp <- get(varName)
#   
#   temp <- temp[grep("STREAM",temp$Effect,invert=TRUE),]
#   
#   assign(varName, temp)  
# }



## remove Introns
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[grep("intron",temp$ANN_Annotation,invert=TRUE),]
  
  assign(varName, temp)  
}

## remove synonymous
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(temp$ANN_Annotation!="synonymous_variant"),]
  
  assign(varName, temp)  
}


## remove intergenic
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(temp$ANN_Annotation!="intergenic_region"),]
  
  assign(varName, temp)  
}

## remove pseudogenes
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(temp$ANN_Transcript_BioType!="pseudogene"),]
  
  assign(varName, temp)  
}

## remove "."
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(temp$ANN_Gene_Name!="."),]
  
  assign(varName, temp)  
}



## select chomosome
# x1 <- x1[which(x1$Chromo="I"),]
# x2 <- x2[which(x2$Chromo=="I"),]
# x3 <- x3[which(x3$Chromo="I"),]
# x4 <- x4[which(x4$Chromo=="II"),]
# x5 <- x5[which(x5$Chromo=="III"),]

## find genes common between unique mutations
combs <- combn(paste0("x", 1:nfiles),2)
ncombs <- length(combs)/2
noncomp <- data.frame("gene"=NULL,"strains"=NULL, stringsAsFactors = FALSE)
for(i in 1:ncombs){
  temp1 <- get(combs[1,i])
  temp2 <- get(combs[2,i])
  noncomplements <- levels(factor(temp1$EFF_Gene_Name[which(is.na(match(temp1$EFF_Gene_Name,temp2$EFF_Gene_Name))==0)])) 
  print(paste(combs[1,i],combs[2,i], length(noncomplements)))
  if(length(noncomplements>0)){
    for(j in 1:length(noncomplements)){
        if(length(which(noncomp$gene==noncomplements[j]))==1) {
        atrow <- which(noncomp$gene==noncomplements[j])
        knownstrains <- as.character(noncomp[atrow,2])
        nowstrains <- paste(knownstrains,combs[1,i],combs[2,i])
        uniquenow <- paste(unique(unlist(strsplit(nowstrains, " "))), collapse = " ")
        noncomp$strains[atrow] <- uniquenow
      }else if(length(which(noncomp$gene==noncomplements[j]))==0){
        noncomp <- rbind(noncomp, data.frame("gene" = noncomplements[j], "strains"=paste(combs[1,i], combs[2,i]), stringsAsFactors = FALSE))
      }else {
        print("WAIT... IM CONFUSED")
      }
        
    }
    name <- paste0(c("noncomp",combs[1,i],combs[2,i]), collapse="_")
    assign(name, noncomplements)
  }
}



## displays number of variants
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  print(paste0(varName," ",length(levels(factor(temp$CPRC)))))
  
  assign(varName, temp)  
}


## write out to a file
dir.create(paste0(dir,finaldir))
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  tempname <- get(paste0("x",i,"name"))
    
  write.table(temp, paste0(dir,finaldir,tempname, ".tsv") ,quote=FALSE, sep="\t", row.names=FALSE)
  
}
write.table(parentalSnpEff, paste0(dir,finaldir, "parentals.tsv") ,quote=FALSE, sep="\t", row.names=FALSE)


write.table(noncomp, paste0(dir,finaldir, "noncomp.csv") ,quote=FALSE, sep="\t", row.names=FALSE)

