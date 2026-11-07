# CellMol_SexChrom_Ovies_aries

**Name:** Jil M. Duhaylongsod

**Study question:**

Can the sex-chromosome complement of an animal be identified through whole-genome sequencing by comparing the read depth of its X chromosome with its autosomes? 

## Part A. Selection and of domestic sheep (Ovis aries) as the study species for sex-chromosome analysis. 


| **Item** | **Value** |
|---|---|
| **Animal Species** | Domestic sheep (*Ovis aries*) |
| **Assembly Accession** | GCA_022416915.1 - https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022416915.1/ |
| **Assembly BioSample** | SAMN23392243 - https://www.ncbi.nlm.nih.gov/biosample/SAMN23392243/ |
| **Sex System** | XY/XX |
| **WGS Run Accession** | SRR17023436 - https://www.ncbi.nlm.nih.gov/sra/?term=SAMN23392243 |
| **WGS BioSample** | SAMN23392243 - https://www.ncbi.nlm.nih.gov/biosample/SAMN23392243 |
| **Same Donor?** | Yes |

## Part B. Retrieve data and map reads in Galaxy

<img width="509" height="302" alt="Screenshot 2026-11-04 205940" src="https://github.com/user-attachments/assets/38811cd4-8996-424e-8a9b-e76703f616c6" />

**Figure 1.** Successfully imported the genomic FASTA file and corresponding sequence/assembly report from the selected NCBI genome assembly into Galaxy for subsequent analysis.

<img width="785" height="441" alt="image" src="https://github.com/user-attachments/assets/9d0f6728-61fc-462c-add8-0207fe7443ad" />

**Figure 2.** Selection of the SRR accession input option in Galaxy's Faster Download and Extract Reads tool to retrieve whole-genome sequencing data from NCBI SRA.

<img width="954" height="470" alt="Screenshot 2026-11-07 185055" src="https://github.com/user-attachments/assets/af96a359-f222-43c3-bd09-2f6804510a64" />

**Figure 3.** Successful completion of the Unzip Collection tool in Galaxy, generating forward and reverse FASTQ datasets for the paired-end whole-genome sequencing reads.

<img width="959" height="448" alt="Screenshot 2026-11-08 014620" src="https://github.com/user-attachments/assets/5a105859-1ab4-4c0b-87ae-4da80bb44001" />

**Figure 4.** Successful FastQC analysis of the forward paired-end WGS reads (SRR17023436), showing the sequencing quality metrics and basic statistics of the imported FASTQ data in Galaxy.

<img width="959" height="437" alt="image" src="https://github.com/user-attachments/assets/34d54d2d-4b3e-4285-9146-64e457d98d37" />

**Figure 5.** Successful FastQC analysis of the reverse paired-end WGS reads (SRR17023436), displaying the sequencing quality metrics and basic statistics of the reverse FASTQ data in Galaxy.

<img width="958" height="446" alt="Screenshot 2026-11-07 201328" src="https://github.com/user-attachments/assets/af36e6b4-08bf-42ea-b3db-1e979eb111a9" />

<img width="959" height="444" alt="Screenshot 2026-11-08 015328" src="https://github.com/user-attachments/assets/ab257096-0553-4713-85f8-9af147d2dea8" />

**Figure 6.** Configuration and successful execution of the BWA-MEM2 tool in Galaxy, using the domestic sheep reference genome (GCA_022416915.1) and paired-end forward and reverse WGS reads (SRR17023436) to generate the aligned BAM file for downstream analysis.

<img width="959" height="449" alt="Screenshot 2026-11-07 235136" src="https://github.com/user-attachments/assets/df2688d9-54a4-4115-a6bd-85555843f800" />

<img width="959" height="440" alt="image" src="https://github.com/user-attachments/assets/bcf73c64-3048-4036-8029-3647fd59da4e" />

**Figure 7.** Successful alignment of paired-end whole-genome sequencing reads (SRR17023436) to the domestic sheep reference genome (GCA_022416915.1) using BWA-MEM2, followed by Samtools coverage analysis to determine chromosome-level sequencing depth for subsequent sex-chromosome dosage analysis.

## Part C. Label chromosomes for the R analysis

| **Accession Number** | **Chromosome Group** |
|---|---|
| CM039926.1 | Autosomes (A) |
| CM039927.1 | Autosomes (A) |
| CM039928.1 | Autosomes (A) |
| CM039929.1 | Autosomes (A) |
| CM039930.1 | Autosomes (A) |
| CM039932.1 | Autosomes (A) |
| CM039952.1 | X Chromosome (X) |

The six autosomal accessions, CM039926.1, CM039927.1, CM039928.1, CM039929.1, CM039930.1, and CM039932.1, were grouped as autosomes because they represent non-sex chromosomes and were used to establish the baseline sequencing depth of the sample. Their coverage values were used to calculate the median autosomal depth, which served as the reference for comparing sex-chromosome coverage. Meanwhile, CM039952.1 was identified as the X chromosome and classified separately from the autosomes. Its sequencing depth was compared with the median autosomal depth to calculate the X:A ratio. This comparison helped determine the sample’s likely sex-chromosome complement based on the relative coverage of the X chromosome.

## Part D Run the analysis

<img width="947" height="454" alt="image" src="https://github.com/user-attachments/assets/51112786-ab47-4085-be5e-618bbf7d1ec1" />

**Figure 8.** Successful creation of a dedicated folder containing the chromosome_coverage.tsv and chromosome_groups.csv data files, along with the sex_depth_analysis.R script, to organize the chromosome coverage data, chromosome classifications, and R code required for sex-chromosome dosage analysis.

**Rscript**

# Sex-chromosome dosage from whole-genome DNA read depth
# Input: Samtools coverage TSV + student chromosome labels
cov <- read.delim("chromosome_coverage.tsv", comment.char = "",
check.names = FALSE, stringsAsFactors = FALSE)
names(cov)[1] <- "rname" # Accepts the '#rname' header
if (!"meandepth" %in% names(cov)) stop("Missing meandepth column")
key <- read.csv("chromosome_groups.csv", stringsAsFactors = FALSE)
if (!all(c("rname", "group") %in% names(key))) stop("Bad key columns")
key$rname <- trimws(key$rname)
key$group <- toupper(trimws(key$group))
if (anyDuplicated(key$rname)) stop("Repeated rname in key file")
if (!all(key$group %in% c("A", "X", "Y", "Z", "W")))
stop("Use only A, X, Y, Z or W in the group column")
x <- merge(key, cov[, c("rname", "meandepth")],
by = "rname", all.x = TRUE, sort = FALSE)
x$meandepth <- suppressWarnings(as.numeric(x$meandepth))
if (anyNA(x$meandepth)) stop("Missing depths: check rname matches")
if (any(x$meandepth < 0)) stop("Negative depth is not valid")
if (sum(x$group == "A") < 5) stop("Label at least 5 autosomes")
if ("X" %in% x$group && "Z" %in% x$group)
stop("Choose either an XY or a ZW species, not both")
if (sum(x$group == "X") == 1 && !any(x$group == "Z")) {
target <- "X"
one_copy <- "XY-like"
two_copies <- "XX-like"
} else if (sum(x$group == "Z") == 1 && !any(x$group == "X")) {
target <- "Z"
one_copy <- "ZW-like"
two_copies <- "ZZ-like"
} else stop("Label exactly one X or exactly one Z chromosome")
A <- median(x$meandepth[x$group == "A"])
if (!is.finite(A) || A <= 0) stop("No usable autosomal coverage")
x$relative_depth <- x$meandepth / A
ratio <- x$relative_depth[x$group == target]
call <- if (ratio >= 0.80) two_copies else if (ratio <= 0.65)
one_copy else "Ambiguous"
if (A < 1) call <- "Inconclusive: very low autosomal depth"

# Additional Y or W evidence can support, but not prove, a call
partner <- if (target == "X") "Y" else "W"
partner_ratio <- if (partner %in% x$group)
x$relative_depth[x$group == partner][1] else NA_real_

# Export results and a simple, reproducible figure
dir.create("results", showWarnings = FALSE)
write.csv(x, "results/chromosome_relative_depth.csv", row.names = FALSE)
png("results/chromosome_depth.png", width = 1000, height = 580)
par(mar = c(8, 4, 3, 1))
barplot(x$relative_depth, names.arg = x$rname, las = 2,
ylab = "Depth relative to median autosome",
main = "Autosomal and sex-chromosome depth", cex.names = 0.75)
abline(h = c(0.5, 1.0), lty = 2)
dev.off()

summary <- c(paste("Median autosome depth:", round(A, 3)),
paste(target, "/ A ratio:", round(ratio, 3)),
paste(partner, "/ A ratio:", round(partner_ratio, 3)),
paste("Dosage interpretation:", call),
"Compare with BioSample sex and discuss limitations.")
writeLines(summary, "results/sex_inference.txt")
cat(paste(summary, collapse = "\n"), "\n")

<img width="1365" height="726" alt="42f72fa3-945d-4ccd-856d-4b1db70fd01d" src="https://github.com/user-attachments/assets/e5917641-5bb0-4dae-9f0f-5396ab508683" />

**Figure 9.** Setting up and opening the sex chromosome depth analysis script in RStudio.


## Part E. Interpret your results

<img width="1000" height="580" alt="chromosome_depth" src="https://github.com/user-attachments/assets/5579359a-9e10-4d80-8993-ced19e237b99" />

**Figure 10.** Bar plot of per-chromosome sequencing depth relative to the median autosome for your friend's dataset. The dashed horizontal lines indicate the expected coverage thresholds for two-copy (1.0) and one-copy (0.5) chromosomal states across the analyzed chromosomes (CM039926.1 through CM039952.1).

The relative depth bar plot displays uniform per-chromosome sequencing coverage across all analyzed chromosomes (CM039926.1 through CM039952.1), with every sequence maintaining a relative depth clustering tightly around the $1.0$ baseline. Because all examined chromosomes align near the two-copy threshold without dropping toward the $0.5$ half-coverage line, the sample exhibits consistent chromosomal dosage representation with no evidence of hemizygous deletions or whole-chromosome copy number variations.

**1. What species did you choose, and what evidence supports a diploid XX/XY system?**

I selected the domestic sheep (Ovis aries) as my study species. Domestic sheep normally have a diploid chromosome number of 54 (2n = 54), consisting of 26 pairs of autosomes and one pair of sex chromosomes. They follow the XX/XY sex-determination system, in which females typically have XX chromosomes and males typically have XY chromosomes.

**2. What are your assembly accession, WGS run, and both BioSample accessions? Are the donors the same?**

The selected genome assembly has the accession number GCA_022416915.1 and the assembly name ASM2242691v1. Its corresponding BioSample accession is SAMN23392243. I also confirmed that the genome assembly and whole-genome sequencing (WGS) data came from the same donor, ensuring that the sequencing results can be appropriately compared with the donor’s recorded biological sex. 

**3. Which sequences were used as autosomes and as X or Z? How did you verify their identities?**

For my analysis, I used seven autosomal sequences and one X chromosome sequence. The autosomes were labeled A, while the X chromosome was labeled X in the chromosome_groups.csv file.
The chromosome identifiers used were:
Autosomes: CM039926.1, CM039927.1, CM039928.1, CM039929.1, CM039930.1, CM039931.1, and CM039932.1.
X chromosome: CM039952.1.
These identifiers were used in the chromosome coverage analysis and assigned to their respective groups based on the chromosome report from the selected genome assembly.

**4. What are the median autosomal depth, X:A or Z:A ratio, and inferred sex-chromosome complement?**

The median autosomal sequencing depth was 22.793, while the calculated X:A ratio was 0.973. The Y:A ratio was recorded as NA because no Y-chromosome depth value was available in the analysis. Since the X:A ratio is close to the expected value of 1.0 for an XX individual, it indicates that the X chromosome has approximately the same sequencing coverage as the autosomes. Based on these results, the sample is inferred to have an XX-like sex-chromosome complement.

**5. Does your result agree with the recorded sex of the WGS donor, if reported? Explain any discrepancy or uncertainty.**

Yes, my results agree with the recorded sex of the WGS donor, which is female. The R analysis showed an X:A ratio of 0.973, indicating an XX-like chromosome complement. This is consistent with the typical XX chromosome arrangement of female domestic sheep. However, the genome assembly and WGS reads have not yet been confirmed to come from the same individual. Therefore, the comparison should be interpreted with caution. The result supports the recorded sex of the WGS donor but does not independently confirm the identity of the assembly donor.

**6. Why is a chromosome-level assembly alone insufficient for a reliable read-depth sex inference?**

A chromosome-level genome assembly alone is insufficient because it only provides the reference DNA sequences and does not show the actual sequencing coverage of the individual being studied. To determine the sex-chromosome complement, whole-genome sequencing reads must be mapped to the reference genome and their coverage compared between the sex chromosomes and autosomes. This allows researchers to estimate whether the individual has one or two copies of the X chromosome. Therefore, read-depth analysis is necessary to make a more reliable inference of the individual's sex-chromosome complement.

## Limitations

This activity had some limitations that may have affected the accuracy of the results. One limitation was the sequencing coverage of the WGS data, since low or uneven coverage can affect the calculation of chromosome depth and the X:A ratio. Another possible limitation was the presence of repetitive DNA sequences, which may cause errors when aligning the reads to the reference genome. Also, the Y:A ratio was recorded as NA, meaning that no Y-chromosome depth value was available for the analysis. Because of this, we could not use the Y chromosome to further confirm the sex-chromosome complement. The analysis also depended on the quality and completeness of the reference genome assembly. Despite these limitations, the X:A ratio of 0.973 was close to the expected value for an XX individual, suggesting that the sample has an XX-like sex-chromosome complement.

## Methods and resources

**NCBI assembly record:** https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022416915.1/

**NCBI SRA record:** https://www.ncbi.nlm.nih.gov/sra/?term=SAMN23392243

**NCBI BioSample record:** https://www.ncbi.nlm.nih.gov/biosample/SAMN23392243

**Galaxy Link:** https://usegalaxy.org/u/jilduhaylongsod/h/sex-dosage-ovis-aries
