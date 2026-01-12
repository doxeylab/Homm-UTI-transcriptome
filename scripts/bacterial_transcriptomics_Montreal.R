# Create heatmap with bacterial transcriptomics results without background noise

library(tximport)
library(Biostrings)
library(stringr)
library(tidyheatmaps)
library(tidyr)
library(dplyr)
library(tidyverse)

create.tx2gene <- function(filepath) {
  #'
  #' 'filepath' points to the cds file for a given pathogen
  #'
  id_file <- readDNAStringSet(filepath)
  seq_name = names(id_file)
  df <- data.frame(entry = seq_name, stringsAsFactors = FALSE)
  df$transcript_id <- str_extract(df$entry, "(lcl\\|[^ ]+)")
  df$locus_id <- str_extract(df$entry, "(?<=\\[locus_tag=)[^\\]]+")
  df$gene_id
  df$protein_id <- str_extract(df$entry, "(?<=\\[protein_id=)[^\\]]+")
  df$protein_name <- str_extract(df$entry, "(?<=\\[protein=)[^\\]]+")
  df <- df[, c("transcript_id", "locus_id", "protein_id", "protein_name")]
  return(df)
}

generate.txi <- function(filepath, tx2gene) {
  #'
  #' 'filepath' points to the location of the salmon output files
  #' 'tx2gene' is the object produced after running create.tx2gene function
  #'
  quant.files <- list.files(filepath, pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)
  sample.names <- basename(dirname(quant.files))
  names(quant.files) <- sample.names
  txi <- tximport(quant.files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE)
  return(txi)
}

all_cds <- "/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/all/all.fna"
all_tx2gene <- create.tx2gene(all_cds)
all_txi <- generate.txi("/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/all/salmon_output", all_tx2gene)

# Add venn info
reads.meta.counts.239$venn_group <- pred.df[match(reads.meta.counts.239$ptid, rownames(pred.df)) , ]$group

# Actually plot heatmap

bt_heatmap <- function(txi, 
                       metadata, 
                       annotation_column = c("ptid", "venn_group", "uti", "ecoli", "klebsiella", "proteus", "enterococcus", "otheruropathogen", "pyelo_sbj", "mainpathogen"), # nolint: line_length_linter.
                       patients) {
  
  # Ensure patient IDs are characters
  patients <- as.character(patients)
  
  # Identify missing patients
  missing_patients <- setdiff(patients, colnames(txi$abundance))
  
  # Create a zero-filled matrix for the missing patients
  if (length(missing_patients) > 0) {
    zeros_mat <- matrix(0, nrow = nrow(txi$abundance), ncol = length(missing_patients))
    colnames(zeros_mat) <- missing_patients
    rownames(zeros_mat) <- rownames(txi$abundance)
    
    # Bind to original matrix
    txi$abundance <- cbind(txi$abundance, zeros_mat)
  }
  
  # Reorder to match 'patients' vector
  txi$abundance <- txi$abundance[, patients]

  # Log-transform
  abundance <- apply(txi$abundance, MARGIN = c(1, 2), FUN = function(x) log10(x + 1))
  
  # Long-format abundance
  txi_long <- as.data.frame(abundance) %>%
    rownames_to_column(var = "Transcript.ID") %>%
    pivot_longer(-Transcript.ID, names_to = "ptid", values_to = "Expression")
  
  # Ensure matching types
  metadata$ptid <- as.character(metadata$ptid)
  txi_long$ptid <- as.character(txi_long$ptid)
  
  # Merge metadata
  txi_long <- left_join(txi_long, metadata[, annotation_column], by = "ptid")
  txi_long <- txi_long[order(txi_long$mainpathogen) , ]
  
  txi_long <- txi_long %>%
    mutate(
      prefix = sub("_.*", "", Transcript.ID),  # grab everything before '_' (or entire string if '_' not present)
      transcriptome = case_when(
        startsWith(prefix, "HGG68") ~ "E. coli",
        startsWith(prefix, "KPHS")  ~ "Klebsiella pneumoniae",
        startsWith(prefix, "WMS")   ~ "Enterococcus faecalis",
        startsWith(prefix, "PMI")   ~ "Proteus mirabilis",
        startsWith(prefix, "PA")    ~ "Pseudomonas aeruginosa",
        TRUE ~ "Unknown"
      )
    )
  
  # assign mainpathogen colours
  mainpathogen_colors <- c(
    "E. coli" = "#E41A1C",            # Red
    "E. coli / Proteus" = "#377EB8",  # Blue
    "Enterococcus" = "#4DAF4A",       # Green
    "Klebsiella" = "#984EA3",         # Purple
    "Proteus" = "#FF7F00",            # Orange
    "Pseudomonas Aeruginosa" = "#A65628"  # Brown
  )
  
  # assign venn colours
  venn_group_colors <- c(
    "All_Negative" = "#1B9E77",         
    "All_Positive" = "#D95F02",
    "Clin_HR_Negative" = "#7570B3",
    "MT_HR_Clin_HR_Negative" = "#E7298A",
    "MT_HR_Clin_Tax_Negative" = "#66A61E",
    "MT_HR_Negative" = "#E6AB02",
    "MT_Tax_Negative" = "#A6761D",
    "MT_Tax_Clin_Tax_Negative" = "#666666",
    "Only_Clin_HR_Positive" = "#1F78B4",
    "Only_Clin_Tax_Positive" = "#33A02C",
    "Only_MT_HR_Positive" = "#FB9A99",
    "Only_MT_Tax_Positive" = "#B15928"
  )
  
  
  # Generate heatmap
  ann_colors <- list(
    mainpathogen = mainpathogen_colors,
    venn_group = venn_group_colors
  )
  
  heatmap <- tidyheatmap(
    df = txi_long, 
    rows = Transcript.ID,
    columns = ptid,
    values = Expression,
    annotation_col = annotation_column[1:length(annotation_column)],
    annotation_row = transcriptome,
    annotation_color = ann_colors,
    show_rownames = FALSE,
    colors = c("white", "red"),
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  
  return(list(heatmap = heatmap, txi_long = txi_long, abundance = abundance))
  
}

# Testing

all_heatmap <- bt_heatmap(all_txi, reads.meta.counts.239, patients = rownames(pred.df))
# View(all_heatmap)




