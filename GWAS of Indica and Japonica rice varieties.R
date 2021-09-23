library(tidyverse)

#Random sample of 10000 SNPs from the full data set
data.geno = read_csv("Rice_44K_genotypes.csv.gz",
                      na=c("NA","00"))
data.geno = data.geno %>% select(-`6_17160794_1`)
names(data.geno)[1] = "ID"
data.geno.10000 = data.geno[,c(1,sample(2:ncol(data.geno),10000))]

geno.numeric = data.geno.10000[,-1] %>%
  lapply(factor) %>%
  as.data.frame() %>% 
  data.matrix()
geno.numeric.fill =
  apply(geno.numeric, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm=T)
    x})

#PCA of first 10 PCs
geno.pca = prcomp(geno.numeric.fill, 
                   rank.=10)
pcvar = geno.pca$sdev^2
pcvar.pct = tibble(pctvar=pcvar/sum(pcvar) * 100,
                    PC=1:length(pcvar))
ggplot(data = pcvar.pct[1:10,], aes(x = PC, y = pctvar)) +
  geom_bar(stat = "identity")

#Population structure visible from PC1 and PC2 & PC2 and PC3
PCs = as_tibble(geno.pca$x) %>%
  mutate(ID=data.geno.10000$ID) %>%
  select(ID, everything())
ggplot(data = PCs, aes(x = PC1, y = PC2)) +
  geom_point()
ggplot(data = PCs, aes(x = PC3, y = PC2)) +
  geom_point()

data.pheno = read_csv("RiceDiversity.44K.MSU6.Phenotypes.csv")
data.pheno.pca = full_join(data.pheno, PCs, by = c("NSFTVID"="ID"))

# PCA plots to explore if subgroups vary by Amylose content, Pericarp color, or Region.  
ggplot(data = data.pheno.pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = `Amylose content`))
ggplot(data = data.pheno.pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = `Pericarp color`))
ggplot(data = data.pheno.pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = `Region`)) +
  scale_color_brewer(type="qual", palette = "Set1")
ggplot(data = data.pheno.pca, aes(x = PC3, y = PC2)) +
  geom_point(aes(color = `Region`)) +
  scale_color_brewer(type="qual", palette = "Set1")

data.geno.10000.fs = matrix("",nrow=nrow(data.geno.10000)*2,ncol=ncol(data.geno.10000)-1+6)
for (i in 1:nrow(data.geno.10000)) {
  data.geno.10000.fs[(i-1)*2+1,1:6] <- data.geno.10000[[i,1]]
  data.geno.10000.fs[(i-1)*2+2,1:6] <- data.geno.10000[[i,1]]
  data.geno.10000.fs[(i-1)*2+1,-1:-6] <- substr(data.geno.10000[i,-1],1,1)
  data.geno.10000.fs[(i-1)*2+2,-1:-6] <- substr(data.geno.10000[i,-1],2,2)
}
data.geno.10000.fs[is.na(data.geno.10000.fs)] = -9
write.table(data.geno.10000.fs,file="rice.data.fastStructure.input.str", col.names = FALSE, row.names = FALSE, quote = FALSE)

fam = tibble(
  FID=data.geno.10000$ID,
  IID=data.geno.10000$ID,
  PID=0,
  MID=0,
  Sex=0,
  Ptype=-9)
write.table(fam,file="rice.data.fastStructure.input.fam",col.names = FALSE, row.names = FALSE, quote = FALSE)
bim = data.geno.10000.fs[,-1:-6]
colnames(bim) = colnames(data.geno.10000)[-1]
bim[bim=="-9"] = NA
bim = apply(bim,2,function(x) unique(na.omit(x)))
bim = t(bim) %>% 
  as_tibble() %>%
  mutate(SNP_ID=colnames(bim), cM=0)
bim = bim %>% 
  separate(SNP_ID,into = c("chromosome","position"),sep="_",remove=FALSE) %>%
  select(chromosome, SNP_ID, cM, position, allele1=V1, allele2=V2)
write.table(bim,file="rice.data.fastStructure.input.bim",col.names = FALSE, row.names = FALSE, quote = FALSE)

#fastStructure population assignment
system("python /usr/local/src/fastStructure/structure.py -K 4 --input=rice.data.fastStructure.input --output=rice.fastStructure.out --format=str")
fs_results = read_delim("rice.fastStructure.out.4.meanQ", delim="  ", col_names = FALSE, col_types = 'nnnn')
fs_results = fs_results %>% 
  mutate(ID=data.geno.10000$ID) %>% 
  select(ID, pop1=X1, pop2=X2, pop3=X3, pop4=X4)
fs_results$assignedPop = apply(fs_results[,-1], 1, which.max)
fs_results$maxPr = apply(fs_results[,2:5],1,max) 
fs_results = fs_results %>% 
  arrange(assignedPop,desc(maxPr)) %>%
  mutate(plot.order=row_number())
fs_results_long = fs_results %>% pivot_longer(pop1:pop4, 
                                               names_to="population",
                                               values_to="proportion")
fs_results_long %>%
  ggplot(aes(x=plot.order, y=proportion, color=population, fill=population)) + 
  geom_col()  +
  ylab("genome proportion") + 
  scale_color_brewer(type="div") + scale_fill_brewer(type="div")

fs_results = fs_results %>% mutate(assignedPop=as.character(assignedPop))
geno.pca.pop = full_join(fs_results, PCs, by = "ID")
ggplot(data = geno.pca.pop, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = assignedPop))
ggplot(data = geno.pca.pop, aes(x = PC3, y = PC2)) +
  geom_point(aes(color = assignedPop))

# GWAs of Flowering time variation
library(statgenGWAS)
pheno.geno.pca.pop = left_join(geno.pca.pop, data.pheno, by=c("ID" = "NSFTVID"))
colnames(pheno.geno.pca.pop) = make.names(colnames(pheno.geno.pca.pop))

# Histogram and Boxplot of 4 populations made by fastStructure
ggplot(data=pheno.geno.pca.pop, aes(x=Flowering.time.at.Aberdeen)) +
  geom_histogram(binwidth = 10) + 
  facet_wrap(facets= ~ assignedPop) +
  ggtitle("Flowering time at Aberdeen by Population")
ggplot(data=pheno.geno.pca.pop, aes(x=Flowering.time.at.Aberdeen)) +
  geom_boxplot() + 
  facet_wrap(facets= ~ assignedPop) +
  ggtitle("Flowering time at Aberdeen by Population")

# ANOVA of Flowering time between populations
aov2 = aov(Flowering.time.at.Aberdeen ~ assignedPop,data=pheno.geno.pca.pop)
summary(aov2) # p-value of 0.0523

# GWAS
data.geno = read_csv("Rice_44K_genotypes.csv",
                      na=c("NA","00"))  %>%
  rename(ID=X1) %>% 
  dplyr::select(-`6_17160794_1`)
data.geno = data.geno %>% as.data.frame() %>% column_to_rownames("ID")
data.map = data.frame(SNP=colnames(data.geno))
data.map = data.map %>%
  separate(SNP, c("chr","pos"), sep = "_", remove = FALSE, convert = TRUE) %>%
  column_to_rownames("SNP")
data.pheno.small = data.pheno %>%
  set_names(make.names(colnames(.))) %>% # fixes the names
  dplyr::rename(genotype=NSFTVID) %>%
  select(genotype, where(is.numeric)) %>%
  as.data.frame()
data.cv = geno.pca.pop %>%
  as.data.frame() %>%
  column_to_rownames("ID")
gData.rice = createGData(geno=data.geno, map = data.map, pheno = data.pheno.small, covar = data.cv)
gData.rice.recode = gData.rice %>% codeMarkers(verbose = TRUE)
data.kinship = kinship(gData.rice.recode$markers)
nullmat = matrix(0, ncol=413,nrow=413, dimnames = dimnames(data.kinship))

gwas.noCorrection = runSingleTraitGwas(gData = gData.rice.recode,
                                       traits = "Flowering.time.at.Aberdeen",
                                       kin = nullmat)
summary(gwas.noCorrection)
plot(gwas.noCorrection, plotType = "qq")
plot(gwas.noCorrection, plotType = "manhattan")

gwas.PCA = runSingleTraitGwas(gData = gData.rice.recode,
                               traits = "Flowering.time.at.Aberdeen",
                               kin = nullmat,
                               covar = c("PC1", "PC2", "PC3", "PC4"))
summary(gwas.PCA)
plot(gwas.PCA, plotType = "qq")
plot(gwas.PCA, plotType = "manhattan")

gwas.K = runSingleTraitGwas(gData = gData.rice.recode,
                             traits = "Flowering.time.at.Aberdeen",
                             kin = data.kinship)
summary(gwas.K)
plot(gwas.K, plotType = "qq")
plot(gwas.K, plotType = "manhattan")

sigSnps = gwas.K$signSnp[[1]]
head(arrange(sigSnps, pValue), 10)
# SNP: 6_9531374 has the highest significance, located within the LOC_Os06g16590 gene.