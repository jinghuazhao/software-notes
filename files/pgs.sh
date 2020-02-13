# 13-2-2020 JHZ

echo "SUMSTATS (Base data)"
# Standard GWAS QC

gunzip -c GIANT.height.gz |\
awk 'NR==1 || ($6 > 0.01) && ($10 > 0.8) {print}' |\
gzip  > Height.gz

# Ambiguous SNPs

gunzip -c Height.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
gzip > Height.noambig.gz

# Duplicate SNPs

gunzip -c Height.noambig.gz |\
awk '{ print $1}' |\
sort |\
uniq -d > duplicated.snp

gunzip -c Height.noambig.gz  |\
grep -vf duplicated.snp |\
gzip - > Height.QC.gz

# Update Effect Size

R --no-save -q <<END
  dat <- read.table(gzfile("Height.QC.gz"), header=T)
  data <- within(dat, {OR <- log(OR)})
  write.table(dat, "Height.QC.Transformed", quote=F, row.names=F)
END

echo "Reference (Target data)"
# Standard GWAS QC

plink \
    --bfile EUR \
    --maf 0.05 \
    --hwe 1e-6 \
    --geno 0.01 \
    --make-bed \
    --mind 0.01 \
    --write-snplist \
    --out EUR.QC

# Clumping

plink \
    --bfile EUR.QC \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump Height.QC.Transformed \
    --clump-snp-field SNP \
    --clump-field P \
    --out EUR

awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp

## PRS
# Generate PRS

awk '{print $1,$8}' Height.QC.Transformed > SNP.pvalue

echo "0.001 0 0.001" > range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

plink \
    --bfile EUR.QC \
    --score Height.QC.Transformed 1 4 11 header \
    --q-score-range range_list SNP.pvalue \
    --extract EUR.valid.snp \
    --out EUR

plink \
    --bfile EUR.QC \
    --indep-pairwise 200 50 0.25 \
    --out EUR

# The first 6 PCs
plink \
    --bfile EUR.QC \
    --extract EUR.prune.in \
    --pca 6 \
    --out EUR

R --no-save -q <<END
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
# Read in the phenotype file 
phenotype <- read.table("EUR.height", header=T)
# Read in the PCs
pcs <- read.table("EUR.eigenvec", header=F)
# The default output from plink does not include a header
# To make things simple, we will add the appropriate headers
# (1:6 because there are 6 PCs)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
# Read in the covariates (here, it is sex)
covariate <- read.table("EUR.covariate", header=T)
# Now merge the files
pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))
# We can then calculate the null model (model with PRS) using a linear regression 
# (as height is quantitative)
null.model <- lm(Height~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
# And the R2 of the null model is 
null.r2 <- with(summary(null.model),r.squared)
prs.result <- NULL
for(i in p.threshold){
    # Go through each p-value threshold
    prs <- read.table(paste0("EUR.",i,".profile"), header=T)
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
    # relevant columns
    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    # Now perform a linear regression on Height with PRS and the covariates
    # ignoring the FID and IID from our model
    model <- lm(Height~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
    # model R2 is obtained as 
    model.r2 <- summary(model)$r.squared
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    prs.coef <- summary(model)$coeff["SCORE",]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    # We can then store the results
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}
# Best result is:
prs.result[which.max(prs.result$R2),]
END
