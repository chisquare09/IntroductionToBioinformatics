library(UniprotR)

# Function to retrieve UniProt data for given protein IDs
retrieve_uniprot_data <- function(ids) {
  # Initialize data frame
  results <- data.frame(Acc_No = character(), Entry.Name = character(),
                        Sequence = character(), Function = character(),
                        Gene.Ontology = character(), Length = numeric(),
                        stringsAsFactors = FALSE)
  
  for (id in ids) {
    cat("Retrieving data for:", id, "\n")
    
    # Get UniProt data with error handling
    name_info <- tryCatch(GetNamesTaxa(id), error = function(e) data.frame(Entry.Name = NA))
    seq_info <- tryCatch(GetSequences(id), error = function(e) data.frame(Sequence = NA))
    function_info <- tryCatch(GetProteinFunction(id), error = function(e) data.frame(Function = NA))
    go_info <- tryCatch(GetProteinGOInfo(id), error = function(e) data.frame(Gene.Ontology = NA))
    length_info <- tryCatch(GetSeqLength(id), error = function(e) data.frame(Length = NA))
    
    # Combine results
    entry_name <- ifelse(nrow(name_info) > 0, name_info$Entry.Name[1], NA)
    sequence <- ifelse(nrow(seq_info) > 0, seq_info$Sequence[1], NA)
    function_text <- ifelse(nrow(function_info) > 0, function_info$Function..CC.[1], NA)
    go_text <- ifelse(nrow(go_info) > 0, go_info$Gene.Ontology..GO.[1], NA)
    length_val <- ifelse(nrow(length_info) > 0, length_info$Length[1], NA)
    
    # Append to results
    results <- rbind(results, data.frame(Acc_No = id, Entry.Name = entry_name,
                                         Sequence = sequence, Function = function_text,
                                         Gene.Ontology = go_text, Length = length_val,
                                         stringsAsFactors = FALSE))
  }
  return(results)
}

# Define protein IDs
id_list <- c("P35568", "O15503", "P01308", "P01344", "P61371")

# Retrieve data
uniprot_data <- retrieve_uniprot_data(id_list)

# Save results
print(uniprot_data)
write.csv(uniprot_data, "uniprot_results.csv", row.names = FALSE)
cat("Results saved to uniprot_results.csv\n")