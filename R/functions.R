#' Load BAM file
#' 
#' This function uses Rsamtools::scanBam to read a user-specified BAM file and
#' converts it to a DataFrame object that is compatible with the next function.
#' 
#' @usage read_bam(bam_file_path)
#' @param bam_file_path Path with forward slashes enclosed in quotation marks
#' @return BAM DataFrame
#' @author Isidora Bijelovic\cr Politecnico di Milano\cr Maintainer: Isidora 
#' Bijelovic\cr E-Mail: <isidora.bijelovic@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Binary_Alignment_Map}\cr
#' @importFrom Rsamtools scanBam
#' @export
read_bam <- function(bam_file_path) {
  bam <- Rsamtools::scanBam(bam_file_path)
  bam_df <- S4Vectors::DataFrame(bam)
  return(bam_df)
}

#' Translate CIGAR string
#' 
#' This function reads the CIGAR string and other information of a specified alignment from the loaded
#' BAM file to visualize its alignment to the reference. It uses BSgenome.Hsapiens
#' libraries to load reference sequences. 
#' 
#' @usage alignment_from_cigar(bam_df, qname_input, genome)
#' @param bam_df BAM DataFrame object
#' @param qname_input Read ID
#' @param genome Genome assembly in quotation marks
#' @return NULL
#' @author Isidora Bijelovic\cr Politecnico di Milano\cr Maintainer: Isidora 
#' Bijelovic\cr E-Mail: <isidora.bijelovic@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Binary_Alignment_Map}\cr
#' @examples
#' \dontrun{
#' data("loaded_Bam", package = "CIGARTranslate") # loads myBam
#' alignment_from_cigar(myBam, "A00491:61:HFCF7DMXX:1:1222:1389:1705", "hg38")
#' }
#' @importFrom Biostrings getSeq
#' @importFrom Biostrings reverseComplement
#' @export
alignment_from_cigar <- function(bam_df, qname_input, genome) {
  if (class(bam_df) != "DFrame") {
    stop("Please import your BAM file using read_bam() first.")
  }
  
  bam_df <- as.data.frame(bam_df)
  bam_row <- subset(bam_df, qname == qname_input)
  
  if (nrow(bam_row) == 0) {
    stop("The specified entry does not exist in the BAM file.")
  }
  
  rname <- bam_row$rname
  pos <- bam_row$pos
  cigar <- bam_row$cigar
  seq <- bam_row$seq
  mapq <- bam_row$mapq
  rlen <- bam_row$qwidth
  strand <- bam_row$strand
  
  genome_package <- switch(as.character(genome),
                           "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
                           "hg38" = "BSgenome.Hsapiens.UCSC.hg38",
                           stop("Sorry, the package only works for hg19 and hg38 genome assemblies for now."))
  
  library(genome_package, character.only = TRUE)
  
  reference_seq <- Biostrings::getSeq(Hsapiens, rname, pos, pos + sum(as.numeric(as.vector(strsplit(gsub("\\D", " ", cigar), " "))[[1]]))-1)
  if (strand == "-") {
    reference_seq <- Biostrings::reverseComplement(reference_seq)
    seq <- Biostrings::reverseComplement(Biostrings::DNAString(seq))
  }
  
  aligned_read <- character()
  aligned_reference <- character()
  
  cigar_parts <- regmatches(cigar, gregexpr("[0-9]+[A-Za-z]", cigar))[[1]]
  
  for (part in cigar_parts) {
    length <- as.integer(gsub("[A-Za-z]", "", part))
    operation <- gsub("[0-9]", "", part)
    
    if (operation %in% c("M", "=", "X")) {
      aligned_read <- paste0(aligned_read, substr(seq, 1, length))
      aligned_reference <- paste0(aligned_reference, substr(reference_seq, 1, length))
      seq <- substr(seq, length + 1, nchar(seq))
      reference_seq <- as.character(substr(reference_seq, length + 1, length(reference_seq)))
    } else if (operation == "D") {
      aligned_read <- paste0(aligned_read, paste(rep("-", length), collapse = ""))
      aligned_reference <- paste0(aligned_reference, substr(reference_seq, 1, length))
      reference_seq <- as.character(substr(reference_seq, length + 1, length(reference_seq)))
    } else if (operation == "I") {
      aligned_read <- paste0(aligned_read, substr(seq, 1, length))
      aligned_reference <- paste0(aligned_reference, paste(rep("-", length), collapse = ""))
      seq <- substr(seq, length + 1, nchar(seq))
    } else if (operation == "S" || operation == "H") {
      seq <- substr(seq, length + 1, nchar(seq))
    }
    
  }
  
  print(paste("CIGAR:", cigar))
  print(paste("Ref: ", aligned_reference))
  print(paste("Read:", aligned_read))
  print(paste("The mapping quality is", mapq))
  print(paste("Strand:", strand))
}