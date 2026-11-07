# BiocManager::install("Rsubread")
library(Rsubread)

fc <- featureCounts(files = "./data/SRR23652846.Aligned.sortedByCoord.out.bam",
                    annot.ext = "./data/ref/Mus_musculus.GRCm39.113.gtf",
                    isGTFAnnotationFile = TRUE,
                    GTF.attrType = "gene_id",
                    useMetaFeatures = TRUE)

fc <- featureCounts(files = "./data/SRR23652846.Aligned.sortedByCoord.out.bam",
                    annot.ext = "./data/ref/Mus_musculus.GRCm39.113.gtf",
                    # annot.inbuilt     = "mm10",               # built-in: "mm10", "mm39", "hg19", "hg38", "R64-1-1" (yeast), etc.
                    # annot.ext       = "genes.gtf.gz",       # alternative: your GTF/SAF file
                    isGTFAnnotationFile = TRUE,
                    GTF.featureType = "exon",
                    GTF.attrType    = "gene_id",
                    
                    isPairedEnd       = TRUE,                 # TRUE for paired-end, FALSE for single-end
                    # strandSpecific    = 1,                    # 0=unstranded, 1=forward (fr-secondstrand), 2=reverse (fr-firststrand)
                    # nthreads          = 8,                    # use your cores (M1 has 8–10)
                    # countMultiMappingReads = TRUE,            # often TRUE for RNA-seq
                    # allowMultiOverlap = FALSE,                # usually FALSE for genes
                    verbose           = TRUE
)

# View the count matrix
head(fc$counts)

# Save counts to CSV
write.csv(fc$counts, file = "./data/SRR23652846_counts.csv", row.names = TRUE)

print("FeatureCounts completed successfully.")

