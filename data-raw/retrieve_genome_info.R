
quickmart = function(genome = 'hg19',
                     host = 'ensembl.org') {
    if (genome == 'hg19') {
        dataset = 'hsapiens_gene_ensembl'
        host = "grch37.ensembl.org"
    }
    else if (genome == 'hg38') {
        dataset = 'hsapiens_gene_ensembl'
        host = "useast.ensembl.org"
    }
    else if (genome == 'mm10') {
        dataset = "mmusculus_gene_ensembl"
    }
    else {
        stop('genome "', genome, '" not found')
    }
    # Choose which species to use and server to download from
    biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
             dataset = dataset,
             host = host,
             port = 80)
}

getGenome = function(mart, add_attribute = 'hgnc_symbol') {
    filters = c('chromosome_name', 'hgnc_symbol')
    values = list(c(1:22, 'X', 'Y'), load('genes.rda'))
    attributes = c('chromosome_name',
                   'band',
                   'strand',
                   'start_position',
                   'end_position')
    attributes = unique(c(attributes, add_attribute))
    result = biomaRt::getBM(attributes = attributes,
                            mart = mart,
                            filters = filters,          
                            values = values,
                            uniqueRows = T)
    if ('hgnc_symbol' %in% add_attribute) {
        result = dplyr::filter(result, hgnc_symbol != "")
    }
    if ('mgi_symbol' %in% add_attribute) {
        result = dplyr::filter(result, mgi_symbol != "")
    }
    result
}
