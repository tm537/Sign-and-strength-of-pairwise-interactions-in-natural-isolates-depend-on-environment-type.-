getwd()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sangeranalyseR")
library(sangeranalyseR)

my_aligned_contigs2 <- SangerAlignment(ABIF_Directory      = "abseq",
                                      processMethod       = "REGEX",
                                      REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                                      REGEX_SuffixReverse = "_[0-9]*_R.ab1$")
writeFasta(my_aligned_contigs2, outputDir = "data")


#### reading in trimmed

# Read the FASTA file
fasta_file2 <- "data/Sanger_contigs_alignment.fa"
seqs <- readDNAStringSet(fasta_file2)

# Convert to data frame
df <- data.frame(
  isolate = names(seqs),         # header names (e.g., >natural_isolates/isolate1)
  seq = as.character(seqs),      # the actual sequences
  stringsAsFactors = FALSE
)

df<-df %>%
  filter(!isolate %in% c("isolate21", "isolate22", "isolate35", "isolate36"))

count(df$isolate)

# do a phylogenetic tree of all consensus sequences
all_SS2 <- DNAStringSet(df$seq)
names(all_SS2) <- df$isolate

library(msa)

tree2 <- msa(all_SS2, method = 'ClustalW') %>%
  as.DNAbin(unmasked(.)) %>%
  ape::dist.dna(., model = 'JC69') %>%
  ape::njs()


options(ignore.negative.edge = TRUE)


library(dplyr)
library(ggtree)


pt_t<-read.csv('site_phylo_t.csv',header=T)

df <- merge(df, pt_t, by = "isolate")


d_labs <- tibble(
  label = tree2$tip.label,  
  isolate_name = sub(".*?(\\d+)$", "\\1", df$isolate),  
  genus = taxa$Genus[match(tree2$tip.label, df$isolate)]  
)



p_all <- ggtree(tree2) %<+% d_labs +
  geom_tiplab(aes(label = isolate_name), size = 4, offset = 0.01, family = 'Helvetica') +
  geom_tiplab(aes(label = paste0('italic(', genus, ')~sp.')),
              parse = TRUE, size = 4, offset = 0.03, family = 'Helvetica') +
  theme(legend.position = 'none') 

p_all




######  


different <- 'green'
general<- 'black'
local <- 'skyblue'


d_labs1 <- tibble(
  label = tree2$tip.label,
  isolate_name = sub("isolate", "", d_labs$label),
  genus = taxa$Genus[match(tree2$tip.label, df$isolate)],
  site = df$site[match(tree2$tip.label, df$isolate)]
)


p_all1 <- ggtree(tree2) %<+% d_labs1 +
  geom_tiplab(aes(label = isolate_name), size = 5, offset = 0.01) +
  geom_tiplab(aes(label = paste0('italic(', genus, ')~sp.')),
              parse = TRUE, size = 5, offset = 0.03) +
  geom_tippoint(aes(color = site), size = 3) +   # color tips by site
  theme(legend.position = 'right') +
  scale_color_manual(values = c(
    "different" = "green",
    "general"   = "black",
    "local"     = "skyblue"
  )) +
  theme(legend.position = "bottom")


p_all1

ggsave("tree.png", plot=last_plot(), path = "data",width=15 ,height=8)



tree2$edge.length <- log(tree2$edge.length + 1)

p_all2 <- ggtree(tree2) %<+% d_labs1 +
  geom_tiplab(aes(label = isolate_name), size = 5, offset = 0.01, family = 'Helvetica') +
  geom_tiplab(aes(label = paste0('italic(', genus, ')~sp.')), parse = TRUE, size = 5, offset = 0.063,       
              hjust = 0.9, family = 'Helvetica') +
  geom_tippoint(aes(color = site), size = 3) +
  scale_color_manual(values = c(
    "different" = "green",
    "general"   = "black",
    "local"     = "skyblue"
  )) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 15),   #
    legend.title = element_text(size = 15),  # 
    legend.key.size = unit(0.5, "cm")        #
  )


p_all2 <- p_all2 + xlim_tree(0.00015)  # smaller value = narrower tree


p_all2


ggsave("tree2.png", plot=last_plot(),width=15 ,height=8)








# taxonomy 

library(tidyverse)

taxa <- dada2::assignTaxonomy(df$seq, 'data/rdp_train_set_16.fa.gz', minBoot = 40) %>%
  data.frame(row.names = NULL, stringsAsFactors = FALSE) %>%
  bind_cols(., df)

## extracting values 

library(Biostrings)
library(msa)
library(ape)


all_SS1 <- DNAStringSet(df$seq)
names(all_SS1) <- df$isolate

# Run multiple sequence alignment
alignment <- msa(all_SS1, method = 'ClustalW')

# Convert to DNAbin and calculate pairwise distances
dist_matrix <- alignment %>%
  unmasked() %>%
  as.DNAbin() %>%
  ape::dist.dna(model = 'JC69')

###

# Convert 'dist' object to a full matrix
dist_df <- as.data.frame(as.matrix(dist_matrix))

# Save to CSV
write.csv(dist_df, "pairwise_distances.csv", row.names = TRUE)






