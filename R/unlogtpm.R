#' @title Convert logTPM to TPM 
#' @description Convert logTPM to TPM, i.e. using 10*(2^(TPM)-1).
#' @param m matrix of logTPM values (gene rows; cell columns)
#' @param bulk if bulk then instead uses 2^(TPM)-1. i.e. no scaling. Default: F
#' @return TPM matrix 
#' @details TPM/10 is used for single cells since 100,000 is a more reasonable estimate than 1,000,000 for the number of RNA transcripts in a cell. 1,000,000 is reasonable estimate for bulk samples that contain multiple cells.
#' @rdname unlogtpm
#' @export 
unlogtpm = function(m, bulk = F) {
    # wrapper around scalop::tpm since scalop::tpm is confusing..
    # in that it does not generate tpm from counts, but rather removes log, scaling and pseudocount
    if (has_dim(m)) m = as.matrix(m)
    if (bulk) x = 1
    else x = 10
    (2^m) * x - 1 
    #x * (2^(m) - 1)
}

