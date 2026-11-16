# The simplest case is to check within the AST classifier, not the ST classifier
# In the AST classifier (built from RC sequences):
#   Non-canonical run → RC(K) is a column, K may also be a column
#   Canonical run     → only one of K or RC(K) is a column (the canonical form)
validate.noncanonical.mode <- function(df, classifier.AST) {
  
k_test  <- df$kmer[1]
rc_test <- chartr("ACGT","TGCA", paste(rev(strsplit(k_test,"")[[1]]),collapse=""))

# Check in the AST classifier, not the ST classifier
cat("K     in AST classifier:", k_test  %in% colnames(classifier.AST), "\n")
cat("RC(K) in AST classifier:", rc_test %in% colnames(classifier.AST), "\n")

# Non-canonical: K = FALSE (or low freq), RC(K) = TRUE  → correct
# Canonical:     K = TRUE,               RC(K) = FALSE  → canonical collapsed them
}
