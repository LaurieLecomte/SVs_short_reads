# Based on previous work of Marc-Andre Lemay (https://github.com/malemay/soybean_sv_paper/blob/master/scripts/fix_sniffles.R)
# Main function -----------------------------------------------------------
add_ALT <- function(input_vcf, output_vcf, refgenome = NULL) {
  
  # Disable scientific notation
  options(scipen=999)
  
  # Checking if the reference genome has been supplied
  if(is.null(refgenome) || !file.exists(refgenome)) {
    stop("Reference genome file is not supplied or does not exist")
  }
  
  # Checking if the reference genome is indexed
  if(!file.exists(paste0(refgenome, ".fai"))) {
    stop("Reference genome .fai index does not exist")
  }
  
  # Opening a connection to read the vcf file
  input_con <- file(input_vcf, open = "rt")
  on.exit(close(input_con), add = TRUE)
  
  # Opening another connection to the output file
  output_con <- file(output_vcf, open = "wt")
  on.exit(close(output_con), add = TRUE)
  
  # Lines are then used one by one to output the header
  while(grepl("^#", 
              cur_line <- scan(input_con, what = character(), sep = "\n", n = 1, quiet = TRUE))
  ) {
    if (grepl('^#CHROM', cur_line)) {
      cat(paste0(cur_line, "\n"), file = output_con)
    }
  }
 
  # Reading the vcf from file wnad ignoring header lines
  vcf <- read.table(input_vcf, comment.char = "#", stringsAsFactors = FALSE)
  
  # Creating indices vectors for each SV type
  dels <- which(grepl("SVTYPE=DEL", vcf[[8]]))
  ins  <- which(grepl("SVTYPE=INS", vcf[[8]]))
  dups <- which(grepl("SVTYPE=DUP", vcf[[8]]))
  invs <- which(grepl("SVTYPE=INV", vcf[[8]]))
  
  # Extracting some useful information for each variant
  chrs   <- vcf[[1]]
  starts <- vcf[[2]]
  altseq <- vcf[[5]]
  svlen  <- as.integer(sub(".*SVLEN=(-?[0-9]+).*", "\\1", vcf[[8]]))
  #ends <- as.integer(sub(".*;END=(-?[0-9]+).*", "\\1", vcf[[8]]))      ######### addition by me
  # END field can be at the beginning or in the middle of INFO fields, so we need to extract accordingly
  ends <- ifelse(test = grepl("^END=", x = vcf[[8]]),
                 yes = as.integer(sub("^END=(-?[0-9]+).*", "\\1", vcf[[8]])),
                 no = as.integer(sub(".*;END=(-?[0-9]+).*", "\\1", vcf[[8]]))
         )
  svtype <- sub(".*SVTYPE=([A-Z]+);.*", "\\1", vcf[[8]])
  
  #ins_len <- as.integer(sub(".*;INSLEN=(-?[0-9]+).*", "\\1", vcf[[8]]))  ######### addition by me
  imprecise <- ifelse(grepl("IMPRECISE", vcf[[8]]), "IMPRECISE", "PRECISE")
  
  consensus <- ifelse(grepl("CONSENSUS", vcf[[8]]), 
                      sub(".*;CONSENSUS=([A-Z]+);.*", "\\1", vcf[[8]]),
                      NA) ######### addition by me, absent if imprecise
  
  #avg_starts <- floor(as.integer(sub(".*;AVG_START=(-?[0-9.]+).*", "\\1", vcf[[8]])))
  #avg_ends <- floor(as.integer(sub(".*;AVG_END=(-?[0-9.]+).*", "\\1", vcf[[8]])))
  #avg_lens <- floor(as.integer(sub(".*;AVG_LEN=(-?[0-9.]+).*", "\\1", vcf[[8]])))
  widths <- ends - starts ######### addition by me
    
  
  supp_vecs <- sub(".*;SUPP_VEC=([0-1]+);.*", "\\1", vcf[[8]])
  supps <- sub(".*;SUPP=([0-9]+);.*", "\\1", vcf[[8]])
    
  callers <- sapply(X = supp_vecs, FUN = function(x) {
    switch(x, 
                    "100" = 'delly',
                    "110" = 'delly + manta',
                    "101" = 'delly + smoove',
                    "010" = 'manta',
                    "011" = 'manta + smoove',
                    "001" = 'smoove',
                    "111" = 'delly + manta + smoove') })
  
  # Computing the replacement information for each variant
  del_info <- del_process(chrs[dels], starts[dels], widths[dels], refgenome = refgenome)
  ins_info <- ins_process(chrs[ins],  starts[ins],  consensus[ins], refgenome = refgenome, 
                          svlen[ins], caller = callers[ins], altseq[ins])
  dup_info <- dup_process(chrs[dups], starts[dups], widths[dups], refgenome = refgenome)
  inv_info <- inv_process(chrs[invs], starts[invs], widths[invs], refgenome = refgenome)
  
  # Assigning the results to the right columns of the vcf file
  vcf[[2]][dels] <- del_info$pos
  vcf[[4]][dels] <- del_info$ref
  vcf[[5]][dels] <- del_info$alt
  vcf[[8]][dels] <- paste0('SVTYPE=', del_info$svtype, ';SVLEN=', del_info$svlen, ';END=', del_info$end,
                          ';SUPP_VEC=', supp_vecs[dels], ';SUPP=', supps[dels])

  vcf[[2]][ins]  <- ins_info$pos
  vcf[[4]][ins]  <- ins_info$ref
  vcf[[5]][ins]  <- ins_info$alt
  vcf[[8]][ins] <- paste0('SVTYPE=', ins_info$svtype, ';SVLEN=', ins_info$svlen, ';END=', ins_info$end,
                          ';SUPP_VEC=', supp_vecs[ins], ';SUPP=', supps[ins])
  

  vcf[[2]][dups] <- dup_info$pos
  vcf[[4]][dups] <- dup_info$ref
  vcf[[5]][dups] <- dup_info$alt
  vcf[[8]][dups] <- paste0('SVTYPE=', dup_info$svtype, ';SVLEN=', dup_info$svlen, ';END=', dup_info$end,
                           ';SUPP_VEC=', supp_vecs[dups], ';SUPP=', supps[dups])

  vcf[[4]][invs] <- inv_info$ref
  vcf[[5]][invs] <- inv_info$alt
  vcf[[8]][invs] <- paste0('SVTYPE=', inv_info$svtype, ';SVLEN=', inv_info$svlen, ';END=', inv_info$end,
                           ';SUPP_VEC=', supp_vecs[invs], ';SUPP=', supps[invs])

  
  # Format genotypes
  for (i in 10:ncol(vcf)){
    vcf[[i]] <- sub("([0-1/.]+):.*", "\\1", vcf[[i]])
  }
  
  # Format FORMAT field
  vcf[[9]] <- 'GT'
  
  
  # Writing the data.frame to the output file
  write.table(vcf, file = output_con, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # Writing
  return(invisible(NULL))
}

# DELs --------------------------------------------------------------------
del_process <- function(chr, start, width, refgenome) {
  
  # The deletion start needs to be before the actual deletion
  pos <- ifelse(start == 1, 1, start - 1)
  
  # The last deleted nucleotide is width units further
  end <- pos + width # do we need to add 1 extra nt ?
  
  # Creating a GRanges object for querying the ref sequence
  ref_range <- GenomicRanges::GRanges(seqnames = chr, 
                                      ranges = IRanges::IRanges(start = pos, end = end))
  
  # Querying the reference sequence
  ref <- Rsamtools::scanFa(refgenome, ref_range)
  
  # Creating a GRanges object for the alt sequence
  alt_range <- GenomicRanges::GRanges(seqnames = chr,
                                      ranges = IRanges::IRanges(start = pos, end = pos))
  
  # Querying the alt sequence
  alt <- Rsamtools::scanFa(refgenome, alt_range)
  
  # Returning the formatted information
  list(pos = as.integer(pos), 
       ref = unname(as.character(ref)), 
       alt = unname(as.character(alt)),
       svlen = as.integer(0 - width),
       end = as.integer(end), 
       svtype = 'DEL'
  )
}


# INSs --------------------------------------------------------------------
ins_process <- function(chr, start, cons, refgenome, ins_length, caller, altseq) {
  
  # The start of the insertion must be offset by one
  pos <- ifelse(start == 1, 1, start - 1)
  
  # Creating a GRanges object for querying the reference sequence
  ref_range <- GenomicRanges::GRanges(seqnames = chr,
                                      ranges = IRanges::IRanges(start = pos, end = pos))
  
  # Querying the reference nucleotide
  ref <- Rsamtools::scanFa(refgenome, ref_range)
  #ref <- unname(as.character(ref))
  
  
  # The alternate sequence is the reference sequence to which the alt_seq is added
  #alt <- paste0(cons, alt_seq)
  # Add consensus seq as alt if consensus is not NA or if callers set includes delly
  #alt <- if (grepl('delly', caller)) {
  #  alt = cons
  #  } 
  alt <- ifelse(test = grepl('delly', caller),
                yes = cons,
                no = altseq)
  
  # Returning the formatted information
  list(pos = as.integer(pos), 
       ref = unname(as.character(ref)), 
       alt = alt,
       svlen = as.integer(ins_length),
       end = as.integer(pos),
       svtype = 'INS'
  )
}



# DUPs --------------------------------------------------------------------
dup_process <- function(chr, start, width, refgenome) {
  
  # The start of the duplication must be offset by one
  pos <- ifelse(start == 1, 1, start - 1)
  
  # The last duplicated nucleotide is width units further
  end <- pos + width
  
  # Creating a GRanges object for querying the alt sequence
  alt_range <- GenomicRanges::GRanges(seqnames = chr,
                                      ranges = IRanges::IRanges(start = pos, end = end))
  
  # Querying the alt sequence
  alt <- Rsamtools::scanFa(refgenome, alt_range)
  
  
  # Creating a GRanges object for querying the reference sequence
  ref_range <- GenomicRanges::GRanges(seqnames = chr,
                                      ranges = IRanges::IRanges(start = pos, end = pos))
  
  # Querying the reference sequence
  ref <- Rsamtools::scanFa(refgenome, ref_range)
  
  # Returning the formatted information
  list(pos = as.integer(pos), 
       ref = as.character(unname(ref)), 
       alt = as.character(unname(alt)),
       svlen = as.integer(width),
       end = as.integer(end),
       svtype = 'DUP')
}



# INVs --------------------------------------------------------------------

inv_process <- function(chr, start, width, refgenome) {
  
  # In the inversion case, the variation truly starts at "start" so we need not update the position
  # However, the end position will be start + width -1
  end <- start + width - 1
  
  # Creating a GRanges object to query the reference sequence
  ref_range <- GenomicRanges::GRanges(seqnames = chr,
                                      ranges = IRanges::IRanges(start = start, end = end))
  
  # Querying the reference sequence
  ref <- Rsamtools::scanFa(refgenome, ref_range)
  
  # The alt sequence is simply the reverse complement
  alt <- sapply(unname(as.character(ref)), revcomp, USE.NAMES = FALSE)
  
  list(pos = as.integer(start), 
       ref = as.character(unname(ref)), 
       alt = alt,
       svlen = as.integer(width), 
       end = as.integer(end),
       svtype = 'INV')
  
}

revcomp <- function(sequence) {
  
  # Checking that only one sequence is provided
  stopifnot(length(sequence) == 1)
  
  # The lookup table that will be used for replacement
  rep_table <- c("A" = "T",
                 "T" = "A",
                 "G" = "C",
                 "C" = "G",
                 "N" = "N")
  
  # Splitting the sequence into its constituent nucleotides
  sequence <- strsplit(sequence, "")[[1]]
  
  # Replacing the nucleotides by their complement
  sequence <- rep_table[sequence]
  
  # Returning the inverted sequence
  paste0(rev(sequence), collapse = "")
}
