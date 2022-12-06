#' @export
RunModularityClustering <- function(
    ## Code from Seurat 
    SNN = matrix(),
    modularity = 1,
    resolution = 0.8,
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = 0,
    print.output = TRUE,
    temp.file.location = NULL,
    edge.file.name = NULL,
    do_par = FALSE
) {
    edge_file <- edge.file.name %||% ''
    adj_size <- pryr::object_size(SNN)
    if (adj_size > 5e8) {
        options(future.globals.maxSize=1.5*adj_size)
    }
    
    clusters_df <- future_map(resolution, function(.resolution) {
        RunModularityClusteringCpp(
            SNN,
            modularity,
            .resolution,
            algorithm,
            n.start,
            n.iter,
            random.seed,
            print.output,
            edge_file
        ) %>% as.matrix(ncol=1)
    }) %>% 
        bind_cols() %>% 
        apply(2, as.character) %>% 
        data.frame()
    
    colnames(clusters_df) <- paste0('Clust', resolution)
    
    return(clusters_df)
}

