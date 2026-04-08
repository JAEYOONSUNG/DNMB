# Internal helper that mirrors the small subset of qdap::beg2char used in DNMB.
# It returns the substring before the nth literal occurrence of `char`.
.dnmb_beg2char <- function(x, char, noc = 1L) {
  if (length(x) == 0L) {
    return(character())
  }

  stopifnot(length(char) == 1L)
  noc <- as.integer(noc)
  if (is.na(noc) || noc < 1L) {
    noc <- 1L
  }

  sep <- gsub("\\\\", "", char)

  vapply(x, function(value) {
    hits <- gregexpr(sep, value, fixed = TRUE)[[1]]
    if (identical(hits[1], -1L) || length(hits) < noc) {
      return(value)
    }
    substr(value, 1L, hits[noc] - 1L)
  }, FUN.VALUE = character(1))
}
