library(sommer)
library(dplyr)
library(gaston)

test_year <- 2019

impVCF <- read.vcf("F:\\WAC_Carter_2022_postimp_filt.vcf.gz")

all.pheno = read.csv("allBLUEs.csv") %>% 
  rename(line = Name1) %>% rename(LOC = ENV)

all.pheno$line <- toupper(all.pheno$line)
all.pheno$line <- chartr(" ","-", all.pheno$line)
all.pheno <- all.pheno[order(all.pheno$line,all.pheno$LOC),]

allM <- as.matrix(impVCF)
allM = allM[sort(row.names(allM)),]
head(row.names(allM))
row.names(allM) <- chartr(" ","-", row.names(allM))
row.names(allM) <- toupper(row.names(allM))

all.pheno <- all.pheno[all.pheno$line %in% row.names(allM), ]
allM <- allM[row.names(allM) %in% all.pheno$line, ]

names <- transform(all.pheno, LOC = reshape::colsplit(LOC, split = "\\_", names = c('Name', 'Year', 'Trial')))
all.pheno$ENV<- names$LOC$Name
all.pheno$Year <- names$LOC$Year

#Fix MAF
swap.scores = function(col){
  if (mean(col) > 1) {return(abs(col -2))} else {
    return(col)}
}

M = apply(allM, MARGIN = 2, FUN = swap.scores)
maf = apply(M, MARGIN = 2, FUN = function(x) mean(x)/2)

#Van Raiden Method G Matrix

M.cent = scale(M, center = TRUE, scale = FALSE)
cov.mat = M.cent%*%t(M.cent)
cov.mat[c(1:3, 71:73), c(1:3, 71:73)]

denom = as.numeric(2*(t(maf)%*%(1-maf)))
G = cov.mat/denom
G[c(1:3, 71:73), c(1:3, 71:73)]

#Subset Test Set
tester <- all.pheno[all.pheno$Year == test_year,]

#Remove Test Set in Train
#Remove Test Set in Train
all.pheno$BUAC[all.pheno$Year == test_year] <- NA
all.pheno$Can_Cover[all.pheno$Year == test_year] <- NA
all.pheno$NDVI[all.pheno$Year == test_year] <- NA
all.pheno$NDRE2[all.pheno$Year == test_year] <- NA
all.pheno$NDRE1[all.pheno$Year == test_year] <- NA
all.pheno$NWI2[all.pheno$Year == test_year] <- NA
all.pheno$NWI1[all.pheno$Year == test_year] <- NA
all.pheno$MTVI[all.pheno$Year == test_year] <- NA


#Index
itt <- list("NDVI","NDRE1","NDRE2","NWI1","NWI2","MTVI")

index <- Map(function(x) {
  
  i.form <- as.formula(paste("BUAC~",x))
  model <- mmer(fixed = i.form,
                random=~vsr(line,Gu=G),
                rcov=~units,nIters=30,
                data=all.pheno, verbose = FALSE)
  model
  
}, itt)

save(index,file = paste0(test_year,"_indices_fixed.Rdata"))

#Multi

itt <- list("NDVI+NWI2","NDVI+NDRE2","NDVI+MTVI")

multi2 <- Map(function(x) {
  
  i.form <- as.formula(paste("BUAC~",x))
  model <- mmer(fixed = i.form,
                random=~vsr(line,Gu=G),
                rcov=~units,nIters=30,
                data=all.pheno, verbose = FALSE)
  model
  
}, itt)

save(multi2,file = paste0(test_year,"_multi2_fixed.Rdata"))


itt <- list("NDVI+NWI2+NDRE2","NDVI+NDRE2+MTVI","NDVI+MTVI+NWI2")

multi3 <- Map(function(x) {
  
  i.form <- as.formula(paste("BUAC~",x))
  model <- mmer(fixed = i.form,
                random=~vsr(line,Gu=G),
                rcov=~units,nIters=30,
                data=all.pheno, verbose = FALSE)
  model
  
}, itt)

save(multi3,file = paste0(test_year,"_multi3_fixed.Rdata"))


multiblup <- mmer(fixed = BUAC~NDVI+NWI2+NDRE2+MTVI,
                  random=~vsr(line,Gu=G),
                  rcov=~units,nIters=30,
                  data=all.pheno, verbose = FALSE)

save(multiblup,file = paste0(test_year,"_multi4_fixed.Rdata"))