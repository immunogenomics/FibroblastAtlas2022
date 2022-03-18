## Adapted from Symphony codebase
library(Rcpp)

# source('/data/srlab/ik936/Roche/referencemapping/symphonycode')
source('/data/srlab2/ik936/Roche/referencemapping/symphonycode/R/buildReferenceFromHarmonyObj.R')
source('/data/srlab2/ik936/Roche/referencemapping/symphonycode/R/utils.R')
sourceCpp('/data/srlab2/ik936/Roche/referencemapping/symphonycode/src/utils.cpp')


.buildReferenceFromHarmonyObj <- function(
    R, Z_orig, Z_corr, betas, #output object from HarmonyMatrix()
    metadata,
    vargenes_means_sds,
    pca_loadings,           # genes `x PCs
    save_uwot_path = NULL,
    verbose = TRUE
) {
    set.seed(111) # for reproducibility
    
    res <- list(meta_data = metadata)
    res$vargenes = vargenes_means_sds
    res$loadings = pca_loadings
    
    res$R <- R
    res$Z_orig <- Z_orig
    res$Z_corr <- Z_corr
    res$betas <- betas
#     res$centroids <- t(singlecellmethods:::cosine_normalize_cpp(R %*% t(Z_corr), 1))
    
    ## centroids already stored in betas
    res$centroids <- singlecellmethods:::cosine_normalize_cpp(betas[1, , ], 1)
    res$cache <- compute_ref_cache(res$R, res$Z_corr)
    res$save_uwot_path <- save_uwot_path
    return(res)
}


.mapQuery <- function(exp_query, 
                     metadata_query, 
                     ref_obj, 
                     vars = NULL, # query batch variables to harmonize over
                     do_normalize = TRUE) { 
    
    if (do_normalize) {
        exp_query <- singlecellmethods::normalizeData(exp_query, 1e4, 'log')
    }
    
    ## Synchronize and scale query genes
    idx_shared_genes = which(ref_obj$vargenes$symbol %in% rownames(exp_query))
    shared_genes = ref_obj$vargenes$symbol[idx_shared_genes]
    exp_query_scaled = singlecellmethods::scaleDataWithStats(
        exp_query[shared_genes, ], ref_obj$vargenes$mean[idx_shared_genes], ref_obj$vargenes$stddev[idx_shared_genes], 1)
    
    # To add rows of zeroes for missing genes, start with full matrix of zeroes
    exp_query_scaled_sync = matrix(0, nrow = length(ref_obj$vargenes$symbol), ncol = ncol(exp_query))  

    # Rows getting filled with exp_query_scaled values
    exp_query_scaled_sync[idx_shared_genes,] = exp_query_scaled
    rownames(exp_query_scaled_sync) = ref_obj$vargenes$symbol
    colnames(exp_query_scaled_sync) = colnames(exp_query)
    
    # 1. Project
    Z_pca_query = t(ref_obj$loadings) %*% exp_query_scaled_sync
    
    # 2. Cluster
    Z_pca_query_cos <- singlecellmethods:::cosine_normalize_cpp(Z_pca_query, 2)
    R_query = soft_cluster(ref_obj$centroids, Z_pca_query_cos, 0.1)
    
    # 3. Correct
    # Make query design matrix
    if (!is.null(vars)) {
        design = droplevels(metadata_query)[,vars] %>% as.data.frame()
        onehot = design %>% 
            map(function(.x) model.matrix(~0 + .x)) %>% 
            reduce(cbind)
        Xq = cbind(1, intercept = onehot) %>% t()
    } else { #Treat query as 1 batch
        Xq = Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))), sparse = TRUE)
    }
    
    Zq_corr = moe_correct_ref(as.matrix(Z_pca_query), 
                              as.matrix(Xq), 
                              as.matrix(R_query), 
                              as.matrix(ref_obj$cache[[1]]), 
                              as.matrix(ref_obj$cache[[2]])) # calls cpp code
        
    return(list(Z = Zq_corr, Zq_pca = Z_pca_query, R = R_query, Xq = Xq, 
                meta_data = metadata_query))
}

                
                
get_ortho <- function(genes_domain, genes_target, from, to) {
    tryCatch({
        library(biomaRt)
    }, error = function(e) {
        stop('Must either provide orthologs_table or install biomaRt')
    })

    attr_domain <- ifelse(from == 'human', 'hgnc_symbol', 'mgi_symbol')
    attr_target <- ifelse(to == 'human', 'hgnc_symbol', 'mgi_symbol')

    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    
    mart_domain <- list(human, mouse)[[as.integer(factor(from, c('human', 'mouse')))]]
    mart_target <- list(human, mouse)[[as.integer(factor(to, c('human', 'mouse')))]]

    orthologs <- getLDS(
        attributesL = attr_domain, 
        filtersL = attr_domain, 
        valuesL = genes_domain, 
        martL = mart_domain, 
        attributes = attr_target, 
        filters = attr_target,
        values = genes_target,
        mart = mart_target, 
        uniqueRows = FALSE
    )
    return(orthologs)
}


dict_to_map <- function(orthologs_dict) {
    colnames(orthologs_dict) <- c('to', 'from') 
    X <- data.table(orthologs_dict)[
        , w := 1 / .N, by = from
    ][]
    i <- factor(X$to)
    j <- factor(X$from)

    map_mat <- sparseMatrix(
        i = as.integer(i), 
        j = as.integer(j), 
        x = X$w, 
        dimnames = list(to = levels(i), from = levels(j))
    )  
    return(map_mat)
}


map_species <- function(exprs_matrix, genes_target, from='mouse', to='human', orthologs_table=NULL, round_fxn=identity) {
    genes_domain <- rownames(exprs_matrix)
    if (is.null(orthologs_table)) {
        message('No orthologues DF provided, pulling data from Biomart')
        orthologs_table <- get_ortho(genes_domain, genes_target, from, to)
    } else {
        message('CAUTION: orthologs_table should have columns MGI.symbol and HGNC.symbol')
    }
    
    ## convert ortholog pairs into sparse linear map 
    linear_map <- dict_to_map(orthologs_table)
    
    ## apply linear map 
    genes_use <- intersect(genes_domain, colnames(linear_map))
    exprs_mapped <- linear_map[, genes_use] %*% exprs_matrix[genes_use, ]
    exprs_mapped <- round_fxn(exprs_mapped)
    return(exprs_mapped)
}

                
## X: feature matrix
## R: cluster probability matrix
## labels: discrete labels to learn 
make_models <- function(X, R, labels, min_cells=8) {
    labels <- as.character(labels)
    label_levels <- unique(labels)

    models <- map(seq_len(ncol(R)), function(k) {
        Rk <- R[, k]
        idx <- which(Rk>0.05) 
    #     types <- factor(ref$meta_data$cell_type, levels = types_all)
        ## Minimum number of cells per class (from glmnet)
        labels_use <- names(which(table(labels[idx]) >= min_cells))
        idx <- intersect(idx, which(labels %in% labels_use))

        if (length(unique(labels[idx]))==1) {
            .label <- unique(labels[idx])
            res <- function(.xnew) {
                res <- matrix(0, nrow=nrow(.xnew), ncol=length(label_levels))
                colnames(res) <- label_levels
                res[, .label] <- 1
                return(res)
            }      
            environment(res) <- rlang::env(.label=.label, label_levels=label_levels)
        } else {
            ## Experimental feature. Poor testing performance 
#             if (use_weights) {
#                 y <- factor(labels[idx])
#                 level_weights <- 1 / prop.table(map_dbl(split(Rk[idx], y), sum)) 
#                 w <- level_weights[labels[idx]]
#             } else {
#                 w <- Rk[idx]
#             }
            
            model <- glmnet(
                x = X[idx, ], 
                y = labels[idx], 
                family = 'multinomial', 
                alpha = 0, 
                lambda = 1, ## TODO: estimate this with CV
                weights = Rk[idx]
            )      
            res <- function(.xnew) {
                pred0 <- glmnet:::predict.multnet(object = model, newx = .xnew, type = 'response')[, , 1]
                res <- matrix(0, nrow=nrow(.xnew), ncol=length(label_levels))
                colnames(res) <- label_levels
                res[, colnames(pred0)] <- pred0
                return(res)
            }
            environment(res) <- rlang::env(model=model, label_levels=label_levels)
        }
        return(res)
    })    
    return(models)    
}


make_predictions <- function(ref_obj, labels_to_predict, query_obj) {
    if (!labels_to_predict %in% names(ref_obj$models)) {
        stop(as.character(glue('{labels_to_predict} not available to predict in this Atlas.')))
    }
    ## per-cluster predictions
    pred <- map(ref_obj$models[[labels_to_predict]], do.call, list(.xnew=t(query_obj$Z)))

    ## dot product with cluster probabilities
    map(seq_len(nrow(query_obj$R)), function(k) {
        Diagonal(x = t(query_obj$R[k, ])) %*% pred[[k]]
    }) %>% 
        purrr::reduce(`+`)

}     
                
make_models_umap <- function(X, R, U, do.cv=FALSE) {
    models <- map(seq_len(ncol(R)), function(k) {
        Rk <- R[, k]
        idx <- which(Rk>0.05) 
        
        if (do.cv) {
            lambda_opt <- cv.glmnet(
                x = X[idx, ], 
                y = U[idx, ], 
                family = 'mgaussian', 
                alpha = 0, 
                weights = Rk[idx]
            )$lambda.1se            
        } else {
            lambda_opt <- 10
        }
        
        model <- glmnet(
            x = X[idx, ], 
            y = U[idx, ], 
            family = 'mgaussian', 
            alpha = 0, 
            lambda = lambda_opt, 
            weights = Rk[idx]
        )      
        res <- function(.xnew) {
            predict(object = model, newx = .xnew)[, , 1]
        }
        environment(res) <- rlang::env(model=model)
        return(res)
    })    
    return(models)    
}

                