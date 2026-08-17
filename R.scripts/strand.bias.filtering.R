strand.bias.filtering <- function(df.subtype, classifier.all) {

  rc_complement <- function(k) {
    chartr("ACGT", "TGCA", paste(rev(strsplit(k, "")[[1]]), collapse = ""))
  }

  # 1. Identify the target subtype from df.subtype
  st <- unique(df.subtype$subtype)
  if (length(st) != 1) stop("df.subtype must contain exactly one subtype.")

  # 2. Extract subtype label from sample_name prefix (format: "{subtype}_...")
  classifier.all$subtype_label <- sub("_.*", "", classifier.all$sample_name)

  # 3. Verify the label matches — print available labels if subtype not found
  available <- sort(unique(classifier.all$subtype_label))
  if (!st %in% available) {
    stop(paste0(
      "Subtype '", st, "' not found in classifier sample_name prefixes.\n",
      "Available prefixes: ", paste(available, collapse = ", "), "\n",
      "Adjust df.subtype$subtype values or the extraction pattern in sub() to match."
    ))
  }

  # 4. Subset rows belonging to the target subtype
  classifier.sub <- classifier.all[classifier.all$subtype_label == st, ]
  
  # 5. k-mers are column names — exclude non-k-mer columns
  non_kmer_cols <- c("sample_name", "subtype_label")
  kmer_cols     <- setdiff(colnames(classifier.all), non_kmer_cols)

    # 4b. Ensure k-mer columns are numeric (guards against factor/character import)
  classifier.sub[, kmer_cols] <- lapply(
    classifier.sub[, kmer_cols, drop = FALSE],
    function(x) as.numeric(as.character(x))
  )

  # 6. Compute mean presence/absence frequency per k-mer across subtype samples
  #    classifier values are 0/1, so mean = proportion of isolates containing k-mer
  #    This is equivalent to the avg frequency used in PORT-EK enrichment output
  kmer_means <- colMeans(classifier.sub[, kmer_cols, drop = FALSE])
  # kmer_means is a named numeric vector: names = k-mer strings, values = mean frequency

    # 7. Diagnostic: confirm that the enriched k-mers themselves are in the classifier
  k_found  <- sum(df.subtype$kmer %in% names(kmer_means))
  rc_check <- sapply(df.subtype$kmer, rc_complement)
  rc_found <- sum(rc_check %in% names(kmer_means))
  message(sprintf(
    "Subtype %s: %d/%d enriched k-mers found in classifier; %d/%d RC(K) found.",
    st, k_found, nrow(df.subtype), rc_found, nrow(df.subtype)
  ))
  
  # 8. Look up RC(K) frequency using direct named-vector indexing
  df.subtype <- df.subtype %>%
    mutate(
      kmer_rc     = sapply(kmer, rc_complement),
      
      # kmer_means[kmer_rc] does a named lookup: returns NA if RC(K) not in k-mer universe
      avg_rc      = unname(kmer_means[kmer_rc]),

      # NA → RC(K) was never observed in any sample (not just below threshold)
      avg_rc_safe = coalesce(avg_rc, 0),

      bias_ratio  = avg / pmax(avg_rc_safe, 1e-6),

      strand = case_when(
        avg_rc_safe < 1e-6 ~ "sense-specific (RC frequency ~ 0)",
        bias_ratio  > 2    ~ "sense-specific",
        bias_ratio  < 0.5  ~ "antisense-specific",
        TRUE               ~ "symmetric"
      )
    )

  return(df.subtype)
}
