TPM <- function(logtpm, dividebyten=TRUE) {
  if(dividebyten) {
    tpm <- 10*(2^(logtpm)-1)}
  else if(!dividebyten) {
    tpm <- 2^(logtpm)-1}
  return(tpm)
}

logTPM <- function(tpm, dividebyten=TRUE) {
  # same as TPM() applies.
  # logTPM(tpm) == TPM(logtpm)

  if(dividebyten) {
    logtpm <- log2(tpm / 10 + 1)}
  else if(!dividebyten) {
    logtpm <- log(tpm+1, 2)}
  return(logtpm)
}
