library(rentrez)

# Function to retrieve NICBI data for given gene list
retrieve_ncbi_data <- function(ids) {
  # Initialize data frame
  results <- data.frame(gene_symbol=character(),
                        full_name = character(),
                        description = character(),
                        length = numeric(),
                        chromosome = character(),
                        map_location = character(),
                        start = numeric(),
                        end = numeric(),
                        strand = character(),
                        exon = numeric(),
                        stringsAsFactors = FALSE)
  
  for (id in ids) {
    cat("Retrieving data for:", id, "\n")
    
    # Get NCBI data 
    search_info <- entrez_search(db="gene",term=paste0(id,"[Gene Name] AND Homo sapiens[Organism]"))
    
    if(length(search_info$ids) ==0 ){
      cat("No result for",id,"\n")
      next
    }
    
    search_id <- search_info$ids[1]
    gene_summary <- entrez_summary(db = "gene", id = search_id)
    
    # Extract information
    gene_full_name <- ifelse(!is.null(gene_summary$nomenclaturename),gene_summary$nomenclaturename,NA)
    gene_description <- ifelse(!is.null(gene_summary$summary),gene_summary$summary,NA)
    gene_chromosome <-ifelse(!is.null(gene_summary$chromosome),gene_summary$chromosome,NA)
    gene_location <- ifelse(!is.null(gene_summary$maplocation),gene_summary$maplocation,NA)
    
    gene_start <- 0
    gene_stop <- 0
    gene_exon <- 0
    # Check if genomic information is available
    if (!is.null(gene_summary$genomicinfo)) {
      gene_start <- gene_summary$genomicinfo$chrstart
      gene_end <- gene_summary$genomicinfo$chrstop
      gene_exon <- gene_summary$genomicinfo$exoncount
    } else {
      gene_start <- NA
      gene_end <- NA
      gene_exon <- NA
    }
    
    gene_length <- ifelse(!is.na(gene_start) & !is.na(gene_end), abs(gene_end - gene_start) + 1, NA) 
    gene_strand <- ifelse(gene_end > gene_start, "forward strand", "reverse strand")
    
    results <- rbind(results,data.frame(gene_symbol=id,
                                        full_name = gene_full_name,
                                        description = gene_description,
                                        length = gene_length,
                                        chromosome = gene_chromosome,
                                        map_location = gene_location,
                                        start = gene_start,
                                        end = gene_end,
                                        strand = gene_strand,
                                        exon = gene_exon,
                                        stringsAsFactors = FALSE))
    
    
    
  }     
  
  return(results)
}

# Define gene list
gene_list <- c("MSH2", "CDKN2A", "HOXB13", "AIP", "CEBPA")

# Retrieve data
ncbi_data <- retrieve_ncbi_data(gene_list)

# Save results
print(ncbi_data)
write.csv(ncbi_data, "ncbi_results.csv", row.names = FALSE)
cat("Results saved to ncbi_results.csv\n")