library(rentrez)

# Function to retrieve NCBI gene data
retrieve_ncbi_data <- function(ids) {
  # Initialize data frame
  results <- data.frame(
    gene_symbol = character(),
    full_name = character(),
    description = character(),
    length = numeric(),
    chromosome = character(),
    start = numeric(),
    end = numeric(),
    strand = character(),
    exon = numeric(),
    num_transcripts = numeric(),
    dna_seq = character(),
    mRNA_seq = character(),
    protein_seq = character(),
    stringsAsFactors = FALSE
  )
  
  for (id in ids) {
    cat("Retrieving data for:", id, "\n")
    
    # Search for the gene ID
    search_info <- entrez_search(db = "gene", term = paste0(id, "[Gene Name] AND Homo sapiens[Organism]"))
    
    if (length(search_info$ids) == 0) {
      cat("No result for", id, "\n")
      next
    }
    
    search_id <- search_info$ids[1]
    gene_summary <- entrez_summary(db = "gene", id = search_id)
    
    # Extract gene information
    gene_full_name <- ifelse(!is.null(gene_summary$nomenclaturename), gene_summary$nomenclaturename, NA)
    gene_description <- ifelse(!is.null(gene_summary$summary), gene_summary$summary, NA)
    gene_chromosome <- ifelse(!is.null(gene_summary$chromosome), gene_summary$chromosome, NA)
    
    gene_start <- NA
    gene_end <- NA
    gene_exon <- NA
    if (!is.null(gene_summary$genomicinfo)) {
      gene_start <- gene_summary$genomicinfo$chrstart
      gene_end <- gene_summary$genomicinfo$chrstop
      gene_exon <- gene_summary$genomicinfo$exoncount
    }
    
    gene_length <- ifelse(!is.na(gene_start) & !is.na(gene_end), abs(gene_end - gene_start) + 1, NA)
    gene_strand <- ifelse(!is.na(gene_end) & !is.na(gene_start) & gene_end > gene_start, "forward strand", "reverse strand")
    print(search_id)
    # Get linked sequence data (nuccore and protein)
    nu_link <- entrez_link(dbfrom = "gene", id = search_id, db = "nuccore")
    protein_link <- entrez_link(dbfrom = "gene", id = search_id, db = "protein")
    
    dna_seq <- entrez_fetch(db="nuccore", id=nu_link$links$gene_nuccore_refseqgene, rettype="fasta")
    mRNA_seq <- entrez_fetch(db="nuccore", id=nu_link$links$gene_nuccore_refseqrna, rettype="fasta")
    protein_seq <- entrez_fetch(db="protein",id=protein_link$links$gene_protein_refseq,rettype = "fasta") 
    
    
    
    # Count the Number of Transcripts
    num_transcripts <- length(nu_link$links$gene_nuccore_refseqrna)
    
    # Store results
    results <- rbind(results, data.frame(
      gene_symbol = id,
      full_name = gene_full_name,
      description = gene_description,
      length = gene_length,
      chromosome = gene_chromosome,
      start = gene_start,
      end = gene_end,
      strand = gene_strand,
      exon = gene_exon,
      num_transcripts = num_transcripts,
      dna_seq = dna_seq,
      mRNA_seq = mRNA_seq,
      protein_seq = protein_seq,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

# Define gene list
gene_list <- c("CHEK2") 

# Retrieve data
ncbi_data <- retrieve_ncbi_data(gene_list)

# Save results
print(ncbi_data)
write.csv(ncbi_data, "ncbi_results_final.csv", row.names = FALSE)
cat("Results saved to ncbi_results_final.csv\n")
