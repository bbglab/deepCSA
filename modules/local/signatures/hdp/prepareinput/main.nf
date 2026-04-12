
process PREPARE_INPUT {

    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/hdp_wrapper'

    input:
    tuple val(meta) , val(type), path(matrix)
    path (features_groups_file)

    output:
    tuple val(meta), path("*.hdp.rds"), path("*.hdp.treelayer.rds") , emit: input_data
    tuple val(meta), path("*.csv")                                  , emit: csv_matrices
    path "versions.yml"                                             , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def features = task.ext.features_groups ? "\'TRUE\'" : "\'FALSE\'"
    def features_groups_list = task.ext.features_groups ?: ""
    """
    # First, create an R script
    cat <<EOF > process_data.R
    data = read.table("${matrix}", header = FALSE)
    rownames(data) <- data[,c(1)]
    colnames(data) <- data[c(1),]
    data <- data[-c(1),]
    data[,-c(1)] <- sapply(data[,-c(1)], as.numeric)
    data <- data[,-c(1)]
    write.table(data, file = "${prefix}.before_round.hdp.csv")
    data <- round(data)
    saveRDS(data, file = "${prefix}.hdp.rds")
    write.table(data, file = "${prefix}.hdp.csv")
    EOF

    # Run the R script
    Rscript process_data.R

    cat <<EOFF > process_metadata.R
    data = read.table("${matrix}", header = FALSE)
    data = data[,c(1, 1)]
    data <- data[-c(1),]
    colnames(data) <- c("sample", "individual")
    data\\\$group = "L"
    if (${features} != 'FALSE') {
        features_table = read.table("${features_groups_file}", header = TRUE, sep = "\\t", check.names = FALSE)
        features_cols = trimws(unlist(strsplit("${features_groups_list}", ",")))
        features_cols = features_cols[features_cols %in% colnames(features_table)]
        if (length(features_cols) > 0) {
            categorical_features = features_table[, features_cols, drop = FALSE]
            categorical_features[] = lapply(categorical_features, as.character)
            categorical_features\\\$sample = features_table[, 1]
            data = merge(data, categorical_features, by.x = "sample", by.y = "sample", all.x = TRUE)
        }
    }
    saveRDS(data, file = "${prefix}.hdp.treelayer.rds")
    write.table(data, file = "${prefix}.hdp.treelayer.csv")
    EOFF

    # Run the R script
    Rscript process_metadata.R
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        HDP     : original
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        HDP     : original
    END_VERSIONS
    """

}
