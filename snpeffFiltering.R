## v3.1.7
## Fixed for new snpEff output style

## import libraries
library(plyr)
library(readr)



## set cutoff. varients present in 'cutoff' or more strains will be removed
cutoff<- 3


dir <- "Z:/WholeGenomeSequenceData/ningAnalysis/"         ## include final '/'



finaldir <- "common_to_1_and_5_only/"
#  chromosomes <- c("I","II","III","IV","V","X")  ## only necesary for the var/chomosome output.

#### preprocessing #####

## find files
files <- list.files(dir, pattern = "tsv")
nfiles <- length(files)


## read in data
allCPRC <- NULL
nRead <- 0
for(i in 1:nfiles){
  print( paste0( c("Reading... ",files[i]),collapse="") )
  pathToFile <- paste0(dir,files[i])
  Ncomments <- max(which(substr(readLines(pathToFile,n=100),1,1)=="#"))
  comment <- readLines(pathToFile, n=Ncomments)
  temp <- read_delim(pathToFile, delim="\t", comment="", skip = Ncomments-1)              ## read in
  if(exists("temp")){if(length(temp)>0){nRead<-nRead+1}}
  
  temp$CPRC <- apply(temp,1,function(x) paste(x[1],x[2],x[4],x[5]))                           ## add CPRC col
  allCPRC <- c(allCPRC,levels(factor(temp$CPRC)))                                             ## add these CPRCs to growing list of all CPRCs
  name <- paste0(c("x",i),collapse="")                                                   ## build variable name
  temp$worm <- name
  fname<- unlist(strsplit(files[i],"\\."))[1]
  CPRCcount <- data.frame(table(allCPRC))
  
  assign(paste0(name,"name"), fname)
  assign(name, temp)                                                                          ## assign var
  
}
print(paste0("Read ",nRead," of ",nfiles," files."))


#### Juicy bits ####### 
#### Run any of the code blocks you want, in any order.

## identify parental variants (>= cutoff)

CPRCcount <- data.frame(table(allCPRC))                          ## counts times each CPRC shows up.
colnames(CPRCcount)[2] <- "numWorms"                                     ## rename  column

ParSummary <- data.frame(table(CPRCcount$numWorms))                ## aggregate into summary of parental var distribution
colnames(ParSummary)[2] <- "count"
barplot(ParSummary$count,names.arg=ParSummary$numWorms,xlab="# of worms", ylab="# of variants") ## plot it out

parental <- CPRCcount[which(CPRCcount$numWorm>=cutoff),]          ## assemble parentals into a DF

### These lines for pulling out vars common to a subgroup but not all. Requires code be run twice: once for all, once for subgroup.
# write.table(x1[match(parental$allCPRC,x1$CPRC),], paste0(dir,finaldir, "varsCommonToAll.csv") ,quote=FALSE, sep=",", row.names=FALSE)
# write.table(x1[match(parental$allCPRC[-match(parental16$allCPRC,parental$allCPRC)],x1$CPRC),], paste0(dir,finaldir, "varsCommonToGroup4.csv") ,quote=FALSE, sep=",", row.names=FALSE)


for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(is.na(match(temp$CPRC,parental$allCPRC))),]             ## remove parentals from worm data
  
  assign(varName, temp)  
}

### Variants present in a group but not others
VAR_MUST_BE_ONLY_IN <- c(1, 5)
CPRCtemp <- CPRCcount[which(CPRCcount$numWorms == length(VAR_MUST_BE_ONLY_IN)),]
for( i in VAR_MUST_BE_ONLY_IN){
  varName <- paste0("x",i)
  temp <- get(varName)
  CPRCtemp <- CPRCtemp[which(CPRCtemp$allCPRC %in% intersect(temp$CPRC, CPRCtemp$allCPRC)),]
}
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  temp <- temp[which(temp$CPRC %in% CPRCtemp$allCPRC),]
  assign(varName, temp)  
}


## Extract homozygous mutations
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(temp$AC=="2"),]
  
  assign(varName, temp)  
}

## Remove upstream and downstream mutations
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[grep("STREAM",temp$EFF_Effect,invert=TRUE),]
  
  assign(varName, temp)  
}



## remove Introns
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(temp$EFF_Effect!="intron_variant"),]
  
  
  assign(varName, temp)  
}

## remove synonymous
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(temp$EFF_Effect!="synonymous_variant"),]
  
  assign(varName, temp)  
}


## remove intergenic
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(temp$EFF_Effect!="intergenic_region"),]
  
  assign(varName, temp)  
}

## Remove indels
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  
  temp <- temp[which(    (apply(as.matrix(as.character(temp[,"REF"])),1,nchar) == 1) & (apply(as.matrix(as.character(temp[,"ALT"])),1,nchar) == 1) ),]
  
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


## displays number of variants per chromosome
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  print(varName)
  for (j in 1:length(chromosomes)){
  print(paste0("  ",chromosomes[j]," ",length(levels(factor(temp$CPRC[which(temp$X.CHROM==chromosomes[j])])))))
  }
  assign(varName, temp)  
}


## write out to a file
append<-"_1n5"

dir.create(paste0(dir,finaldir))
for( i in 1:nfiles){
  varName <- paste0("x",i)
  temp <- get(varName)
  tempname <- get(paste0("x",i,"name"))
  writename <- paste0(tempname,append)
    
  writeTo <- paste0(dir,finaldir,writename, ".csv")
  write.table(temp, writeTo ,quote=FALSE, sep=",", row.names=FALSE)
  print( paste0( c("Wrote to ",writeTo),collapse="") )
}

write.table(noncomp, paste0(dir,finaldir, "noncomp.csv") ,quote=FALSE, sep=",", row.names=FALSE)

