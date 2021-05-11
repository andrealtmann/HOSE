### description

#HOSE: HOmogeneity of SNP Effect sizes

#This script can be used to test for the
#homogeneity of odds ratios (ORs) between GWAS summary statistics
#The work is based on the Woolf test for hmogeneity of ORs. A similar method
#using the raw data has been used by Liu et al. (2012) Hum Genet. https://pubmed.ncbi.nlm.nih.gov/21858542/


# required packages
### required packages
#for quickly loading the summary statitics
library(data.table)
#for plotting
library(qqman)

### internal helper function ###

#convert chr and BP into a label (chr-BP)
mk_labels <- function(x){
  res <- apply(x, 1, function(y){
    paste(as.integer(y[1]), as.integer(y[2]), sep="-")
  })
  return(res)
}

#compute Odds ratio from MAF in cases and MAF in controls
getOR <- function(cases_maf, controls_maf){
  myOR = (cases_maf * (1-controls_maf)) / ((1-cases_maf)*controls_maf)
  return(myOR)
}

#wrapper function for data.table's fread (so .gz files can be used)
fread_wrap <- function(fname){
  f_str <- fname
  fend <- rev(unlist(strsplit(fname, split="\\.")[[1]]))[1]
  if (fend == "gz"){
     f_str <- paste("gzcat ", fname)
  }
  return(fread(f_str, data.table=F))
}

### plotting functions ###

#wrapper for qqman's qq function
woolf_qq <- function(x){
  qq(x[,1])
}

#wrapper for qqman's manhattan function
woolf_man <- function(inp, pv, ...){
  dummy <- data.frame(inp[["GWAS_1"]], P=pv)
  manhattan(dummy, chr="V1", bp="V2", snp="V3", ...)
}

### core functions ###

#GWAS sample size for case and controls
mk_sample_size_info <- function(N_cases, N_controls){
  res <- c(2*N_cases, 2*N_controls)
}

#Harmonize the data so that the same rows are available for GWAS1, GWAS2 and the MAF info
#gwas_a: filename to the first GWAS
#gwas_b: filename to the second GWAS summary stats
#maf: filename to the maf in the population (assumed healthy!)
harmonize_input <- function(gwas_a, gwas_b, maf){
  message("reading GWAS A: ", gwas_a)
  gwa1 <- fread_wrap(gwas_a)
  rownames(gwa1) <- mk_labels(gwa1)


  message("reading GWAS B: ", gwas_b)
  gwa2 <- fread_wrap(gwas_b)
  rownames(gwa2) <- mk_labels(gwa2)

  message("reading minor allele frequencies: ", maf)
  mafs <- fread_wrap(maf)
  rownames(mafs) <- mk_labels(mafs)


  keep_snp <- intersect(intersect(rownames(gwa1), rownames(gwa2)), rownames(mafs))
  message(length(keep_snp))

  res <- list()
  res[["GWAS_1"]] <- gwa1[keep_snp,]
  res[["GWAS_2"]] <- gwa2[keep_snp,]
  res[["MAF"]]    <- mafs[keep_snp,]

  return(res)
}

#one row input:
#BETA1       SE1       P1 ACA1  ACO1   BETA2    SE2     P2  ACA2  ACO2      MAF
sumstat_woolf_test <- function(x, add_n=0.00001){

  x_maf <- as.double(x["MAF"])

  Rh = x_maf/(1.0 - x_maf)
  mafs <- c()
  OR <- c()

  #For GWAS 1 (convert sample size and effect size and maf to maf for cases):
  ca1_NCHROBS <- as.integer(x["ACA1"])
  co1_NCHROBS <- as.integer(x["ACO1"])
  #NOTE: beta converted to odds_ratio
  Dn1 = ca1_NCHROBS / (exp(as.double(x["BETA1"])) * Rh + 1)
  De1 = ca1_NCHROBS - Dn1
  co1_MAF = x_maf
  ca1_MAF = De1 / ca1_NCHROBS
  OR <- c(OR, getOR(ca1_MAF + add_n, co1_MAF + add_n))
  mafs <- c(mafs, ca1_MAF, co1_MAF)

  #For GWAS 2 (same for gwas 2):
  ca2_NCHROBS <- as.integer(x["ACA2"])
  co2_NCHROBS <- as.integer(x["ACO2"])
  #NOTE: beta converted to odds_ratio
  Dn2 = ca2_NCHROBS / (exp(as.double(x["BETA2"])) * Rh + 1)
  De2 = ca2_NCHROBS - Dn2
  co2_MAF = x_maf
  ca2_MAF = De2 / ca2_NCHROBS
  OR <- c(OR, getOR(ca2_MAF + add_n, co2_MAF + add_n))
  mafs <- c(mafs, ca2_MAF, co2_MAF)

  #add a little count to avoid maf of 0
  mafs <- mafs + add_n

  #non-carriers cases, #carrier cases, #non-carrier controls, #carrier controls
  w <- c()
  tmp <- c( (1.0 - ca1_MAF) * ca1_NCHROBS,     ca1_MAF  * ca1_NCHROBS,     (1.0 - co1_MAF) * co1_NCHROBS,     co1_MAF  * co1_NCHROBS)
  w <- c(w, 1.0/sum(1.0 / (tmp + add_n)))
  tmp <- c( (1.0 - ca2_MAF) * ca2_NCHROBS,     ca2_MAF  * ca2_NCHROBS,     (1.0 - co2_MAF) * co2_NCHROBS,     co2_MAF  * co2_NCHROBS)
  w <- c(w, 1.0/sum(1.0 / (tmp + add_n)))

  logOR <- 0.0
  OR <- log(OR)
  logOR = sum(OR * w)/sum(w)

  chisq_sta = sum(w * (OR - logOR)^2)
  pv = pchisq(chisq_sta, 1, lower.tail=F)
  res <- c(pv, logOR, OR[1], OR[2])

}


##woolf_test for two summary stats
#input:
woolf_test <- function(inp, sample1=NA, sample2=NA){

  allele_not_eq <- sum(inp[["GWAS_1"]][,4] != inp[["GWAS_2"]][,4])
  if (allele_not_eq > 0){
    stop("GWAS 1 and GWAS 2 have mismatching alleles")
  }
  allele_not_eq <- sum(inp[["GWAS_1"]][,4] != inp[["MAF"]][,4])
  if (allele_not_eq > 0){
    stop("GWAS 1 and MAF have mismatching alleles")
  }

  #for gwas 1/2:
  #V6=beta     V7=se     V8=p    V9=N_ca    V10=N_Co

  #get information for GWAS1
  if (ncol(inp[["GWAS_1"]]) < 10){
    p1 <- cbind(inp[["GWAS_1"]][,6:8], NA, NA)
  } else {
    p1 <- inp[["GWAS_1"]][,6:10]
  }
  #add sample size information
  if (is.na(sample1)){
    p1[,4:5] <- p1[,4:5] * 2
  } else {
    p1[,4] <- sample1[1]
    p1[,5] <- sample1[2]
  }
  colnames(p1) <- paste(c("BETA","SE","P","ACA","ACO"), 1, sep="")

  if (is.na(p1[1,"ACA"]) | is.na(p1[1,"ACO"])){
    stop("GWAS 1 sample size missing; either add sample size columns to the summary statistics or provide number of cases and controls.")
  }

  #get information for GWAS1
  if (ncol(inp[["GWAS_2"]]) < 10){
    p2 <- cbind(inp[["GWAS_2"]][,6:8], NA, NA)
  } else {
    p2 <- inp[["GWAS_2"]][,6:10]
  }
  #add sample size information
  if (is.na(sample2[1])){
    p2[,4:5] <- p2[,4:5] * 2
  } else {
    p2[,4] <- sample2[1]
    p2[,5] <- sample2[2]
  }
  colnames(p2) <- paste(c("BETA","SE","P","ACA","ACO"), 2, sep="")

  if (is.na(p2[1,"ACA"]) | is.na(p2[1,"ACO"])){
    stop("GWAS 2 sample size missing; either add sample size columns to the summary statistics or provide number of cases and controls.")
  }


  MAF <- as.double(inp[["MAF"]][,6])

  mdat <- cbind(p1, p2, MAF)
  print(head(mdat))

  return(t(apply(mdat,1, sumstat_woolf_test)))

}

#correct p-values for genomic inflation
gc_correct <- function(pv){

  chisq_stat = qchisq(pv,1, lower.tail=F)
  lambda_gc = median(chisq_stat)/qchisq(0.5,1)

  message("genomic inflation: ", lambda_gc)

  new_chisq_stat = chisq_stat/lambda_gc
  new_pv = pchisq(new_chisq_stat, df=1, lower.tail=F)

  return(new_pv)
}
