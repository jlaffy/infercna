TPM <- function(m, dividebyten=TRUE) {
  if(dividebyten) {
    m <- 10*(2^(m)-1)}
  else if(!dividebyten) {
    m <- 2^(m)-1}
  m
}

logTPM <- function(m, dividebyten=TRUE) {
  # same as TPM() applies.
  # logTPM(tpm) == TPM(logtpm)

  if(dividebyten) {
    m <- log2(m / 10 + 1)}
  else if(!dividebyten) {
    m <- log(m + 1, 2)}
  m
}
