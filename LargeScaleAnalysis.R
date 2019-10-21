## Oct 21 2019 ## 

## Initial Analyses of the large scale CRISPR Screen with Brunello

# We will make sure that the barcode to sample associations are correct
# Look at the fastq files which should make have the barcodes in the readname
# CCGAGTTA -> Bru-3

setwd('~/Documents/UT Austin Postdoc/Analysis/Bru_hyg_day0_sequencing/')
control_1 = read.table('./Bru_hyg_1_1.txt', header= T)
control_2 = read.table('./Bru_hyg_2_1.txt', header= T)

increased_1 = read.table('./Bru_hyg_3_1.txt', header= T)

reduced_1 = read.table('./Bru_hyg_4_1.txt', header= T)
reduced_2 = read.table('./Bru_hyg_5_1.txt', header= T)
reduced_3 = read.table('./Bru_hyg_6_1.txt', header= T)

plasmid = read.table('./Bru_hyg_parent.txt', header = T)

all_sgRNA_counts = cbind(plasmid$Count, control_1$Count, control_2$Count,
                         increased_1$Count, 
                         reduced_1$Count, reduced_2$Count, reduced_3$Count )
row.names( all_sgRNA_counts )  = plasmid$sgRNA
colnames(all_sgRNA_counts) = c( "Plasmid", "Control_1", "Control_2", 
                                "Increased_1", 
                                "Reduced_1", "Reduced_2", "Reduced_3" )
# write.csv(all_sgRNA_counts, file = "All_sgRNA_Counts_Merged.csv")

# Total Number of reads per condition
colSums(all_sgRNA_counts)

# We will define a new matrix that is normalized to sequencing depth
# In future work, we can think of better normalization methods
normalized_all_sgRNA_counts = apply(all_sgRNA_counts , 1, function (x) {1000000 * x / colSums(all_sgRNA_counts) })
normalized_all_sgRNA_counts  = t(normalized_all_sgRNA_counts)

# Number of sgRNAs with counts less than threshold
threshold = 1
round ( apply (all_sgRNA_counts, 2, function (x){sum (x < threshold) }) / nrow(all_sgRNA_counts), 2 ) 

plot(all_sgRNA_counts[,2], all_sgRNA_counts[,3], pch = 19, cex = .2, xlab = "Control_1", ylab = "Control_2")
plot(all_sgRNA_counts[,1], all_sgRNA_counts[,3], pch = 19, cex = .2, xlab = "Plasmid", ylab = "Control_2", 
     xlim = c(0, 300))

# Overall correlation between the sgRNA counts
round ( cor (all_sgRNA_counts, method = "spearman") , 3 ) 


# Is there a depletion for essential genes in Control population compared to Plasmid? 
essentials = read.csv('./K562_essential_genes_depmap.csv')

gene_names = unlist(lapply ( strsplit(row.names(all_sgRNA_counts), split = "_"), "[[" , 1 ) )
color_vector = rep ("gray", nrow (all_sgRNA_counts))
color_vector[gene_names %in% essentials$Genes ]  = "red"
plot(all_sgRNA_counts[,1], all_sgRNA_counts[,3], pch = 19, cex = .2, xlab = "Plasmid", ylab = "Control_2", 
     xlim = c(0, 300), col = color_vector)

boxplot (  ( all_sgRNA_counts[,2] + all_sgRNA_counts[,3] ) / all_sgRNA_counts[,1] ~ color_vector, 
           varwidth = T, ylab = "Control Sum to Plasmid Ratio", names = c("Nonessential", "Essential")) 

boxplot (  ( all_sgRNA_counts[,2] + all_sgRNA_counts[,3] ) / all_sgRNA_counts[,1] ~ color_vector, 
           varwidth = T, ylab = "Control to Plasmid Ratio", names = c("Nonessential", "Essential"), 
           ylim = c (0, 3 ) ) 

# What is the ratio of the Control sgRNAs in Control to Plasmid, or Control to Sorted? 
control_sgRNAs = normalized_all_sgRNA_counts[gene_names == "Control", ]
boxplot(control_sgRNAs)

# Read counts for known translation factors 
essential_translation_regulators = read.csv('./essential_translation_genes_GO.csv')
sgRNAs_essential_translation =   gene_names %in% essential_translation_regulators$x 
nonessential_translation_regulators = read.csv('./nonessential_translation_genes_GO.csv')
sgRNAs_nonessential_translation =   gene_names %in% nonessential_translation_regulators$x  

# One strategy is to take the sgRNAs that enriched compared to control and look for enrichment
control_mean = (normalized_all_sgRNA_counts[,2] + normalized_all_sgRNA_counts[,3])  / 2
reduced_mean = (normalized_all_sgRNA_counts[,6] + normalized_all_sgRNA_counts[,7])  / 2

ratio_reduced =  reduced_mean / (control_mean + 1 )  

enriched_in_reduced = ratio_reduced > 1.25
sum ( gene_names == "Control" ) 
write.csv ( sort ( table ( gene_names[enriched_in_reduced] ), decreasing = T ) , file = "Enriched_Genes_Reduced.csv") 

ratio_increased =  normalized_all_sgRNA_counts[,4] / (control_mean + 1 )  
enriched_in_increased = ratio_increased > 1.6
sort ( table ( gene_names[enriched_in_increased] ), decreasing = T ) 

# Ovelap between increased and reduced
multiple_sgRNAs_reduced = table ( gene_names[enriched_in_reduced] ) > 1
list_of_genes_reduced = unlist ( dimnames(table ( gene_names[enriched_in_reduced] )) ) 
length ( list_of_genes_reduced[multiple_sgRNAs_reduced] )

multiple_sgRNAs_increased = table ( gene_names[enriched_in_increased] ) > 1
list_of_genes_increased = unlist ( dimnames(table ( gene_names[enriched_in_increased] )) ) 
length( list_of_genes_increased[multiple_sgRNAs_increased] )

sum ( list_of_genes_reduced[multiple_sgRNAs_reduced] %in% list_of_genes_increased[multiple_sgRNAs_increased] ) 

length(unique(gene_names) ) 

# Two guides
fmat = matrix(nrow =2, ncol = 2)
fmat[1,] = c(212, 1958-212)
fmat[2,] = c(1484- 212, 19115 - (1958 +1484 - 212)  )
fisher.test(fmat)

# Three guides
fmat[1,] = c(2, 212)
fmat[2,] = c(121, 19115 - 335  )
fisher.test(fmat)

# Is there an enrichment of these lists among the GO term translation regulators. 
length ( unique ( gene_names[ sgRNAs_nonessential_translation ]  )  ) 
sum ( list_of_genes_reduced[multiple_sgRNAs_reduced] %in% unique ( gene_names[ sgRNAs_nonessential_translation ]  )  ) 

fmat[1,] = c(71, 1958 - 71)
fmat[2,] = c(670, 19115 - ( 670 + 1958 - 71) )
fisher.test(fmat)


length ( unique ( gene_names[ sgRNAs_essential_translation ]  )  ) 
sum ( list_of_genes_reduced[multiple_sgRNAs_reduced] %in% unique ( gene_names[ sgRNAs_essential_translation ]  )  ) 

fmat[1,] = c(28, 1958 - 28)
fmat[2,] = c( 249 - 28, 19115 - ( 249 + 1958 - 28) )
fisher.test(fmat)

# Another alternative is to look at the read counts of all sgRNAs that map to the translation regulators
boxplot( normalized_all_sgRNA_counts[sgRNAs_essential_translation,], 
         main = "Essential_Translation All sgRNAs" )
write.csv(normalized_all_sgRNA_counts[sgRNAs_essential_translation,], file = "Normalized_sgRNA_Essential_Translation.csv")

boxplot( normalized_all_sgRNA_counts[sgRNAs_nonessential_translation,], ylim = c(0, 80),
         main = "Nonessential_Translation All sgRNAs" )

write.csv(normalized_all_sgRNA_counts[sgRNAs_nonessential_translation,], file = "Normalized_sgRNA_NonEssential_Translation.csv")
