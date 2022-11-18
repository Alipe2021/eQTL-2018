library(data.table)
library(tidyverse)
library(CMplot)
library(cpgen)
library(fdrtool)
set_num_threads(8)

rm(list = ls())
plink.exe <- "/opt/bio/plink-1.9.0/plink --threads 8"
emmax.exe <- "/opt/bio/Emmax-20120210/emmax-intel64"
ekinf.exe <- "/opt/bio/Emmax-20120210/emmax-kin-intel64"

workdir = "~/LabDisk2T/03-Project/eQTL.20221112/08-GeneE-Pheno/TWAS/EMMAX/"
setwd(workdir)
RawGeno = "../../../01-Genotype/Genotype.vcf.gz"
RawPhen = "../../../02-Phenotype/Phenotype.txt"
RawExpr = "../../../03-GeneExpr/GeneTPM.txt"
GenInfo = "~/LabDisk2T/00-PubData/11-DataBase/zma_gene_miRNA.loc.txt"

exp_data <- as.data.frame(fread(RawExpr, header = T))
phe_data <- as.data.frame(fread(RawPhen, header = T))
gen_date <- as.data.frame(fread(GenInfo, header = F, sep = "\t"))
colnames(gen_date) <- c("Gene", "Chr", "Start", 
                        "End", "Strand", "Type", "Annotation")
exp_data[1:5, 1:5]
phe_data[1:5, 1:5]
gen_date[1:5, 1:5]

# for loop
for (i in 4:ncol(phe_data)) {
  
  # i <- 2
  trait_name <- colnames(phe_data)[i]
  out_dir <- paste(workdir, trait_name, sep = "/")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  out_pfx <- paste(out_dir, trait_name, sep = "/")
  # 
  list_file <- paste0(out_pfx, "_List.txt")
  phen_file <- paste0(out_pfx, "_Pheno.txt")
  gene_file <- paste0(out_pfx, "_Expr.txt")
  geno_file <- paste0(out_pfx, "_Geno")
  kin_file <- paste0(out_pfx, ".aIBS.kinf")
  res_file <- paste0(out_pfx, "_LMM.xls")
  
  # Step 01: Output Phenotype data
  trait_dat <- na.omit(phe_data[, c(1, i)])
  head(trait_dat); dim(trait_dat)
  
  # Step 01: Fetch sample list
  sample_list <- intersect(phe_data$Taxa, colnames(exp_data))
  length(sample_list)
  write.table(x = sample_list, file = list_file, quote = F, 
              row.names = F, col.names = F )
  
  trait_dat <- subset(trait_dat, subset = trait_dat$Taxa %in% sample_list)
  trait_dat <- arrange(trait_dat, Taxa)
  fwrite(trait_dat, file = phen_file, sep = "\t")

  # Step 02: Calculate Kinship
  cmd1 <- paste(plink.exe, "--vcf", RawGeno, "--keep-fam", list_file,
               "--maf 0.05 --recode 12 transpose", "--out", geno_file)

  cat("\n", paste0(rep("=",80),collapse = ""),
      "\nCommand1:\n", cmd1, 
      "\n", paste0(rep("=",80),collapse = ""),
      "\n")
  system(cmd1, intern = FALSE)
  
  cmd2 <- paste(ekinf.exe, "-s -d 6 -M 10 -o", kin_file, geno_file)
  cat("\n", paste0(rep("=", 80), collapse = ""),
      "\nCommand2:\n", cmd2,
      "\n", paste0(rep("=", 80), collapse = ""),
      "\n")
  system(cmd2, intern = FALSE)

  # Step 03: Format Gene expression 
  kinf_data <- as.matrix(fread(kin_file))
  expr_data <- column_to_rownames(exp_data, "GeneID")
  expr_data <- subset(expr_data, select = colnames(expr_data) %in% sample_list)
  expr_data <- as.matrix(t(expr_data))
  expr_data[1:5, 1:5]; dim(expr_data)
  kinf_data[1:5, 1:5]  

  # y <- trait_dat[, 2]
  out.res <- cGWAS.emmax(y = as.vector(trait_dat[, 2]), M = expr_data, 
                         A = kinf_data, verbose = TRUE, seed = 777)  
  
  out.res <- as.data.frame(out.res)
  rownames(out.res) <- colnames(expr_data)
  head(out.res)
  
  out_dat <- arrange(out.res, p_value)
  out_dat['fdr'] <- as.numeric(p.adjust(out_dat$p_value, method = "fdr"))
  head(out_dat) 
  out_dat <- rownames_to_column(out_dat, "Gene")
  out_dat <- merge(out_dat, gen_date, by = "Gene")
  fwrite(out_dat, file = res_file, quote = F, row.names = F)

  gc()
}

  