suppressPackageStartupMessages({
    # library(sf)
    source('../libs.R')
    # library(gagarin)
    ## For Louvain clustering 
    source('~/ModularityClustering/R/modularity_clustering.R')
    sourceCpp('~/ModularityClustering/src/RModularityOptimizer.cpp')
    source('../../R/utils_plotting.R')
    source('../../R/utils.R')
})

colors_types <- list(DP='lightgrey', Immuno=muted('red'), Perivascular=muted('blue'), Other='lightgrey')

plot_niche_expression <- function(niches, gene) {
    spots$metadata %>% 
        ggplot(aes(y, -x)) + 
            geom_point(data = . %>% subset(Niche %in% niches), size = .2, alpha = .4, shape = 21, fill='white') + 
            geom_point(data = cbind(cells$metadata, as.matrix(spots$z)), aes_string(color = gene), shape = '.', alpha = .4) + 
            scale_color_gradient2(low = 'white', high = muted('red')) + 
            labs(title = paste(niches, collapse = ' ')) + 
            NULL    
}


pal <- list(
    Other = 'lightgrey',
    Perivascular = muted('red'),
    Vascular = muted('blue'),
    Lymphoid = 'yellow',
    Lympho_vascular = 'green'
)

make_spots <- function(cells, pow=1) {
    ## Do Delaunay triangulation to find spatial neighbors 
    snn <- gagarin::getSpatialNeighbors(dplyr::select(cells$metadata, x, y))
    snn <- snn + Diagonal(n=nrow(snn)) ## add self as neighbor     
    for (i in seq_len(pow-1)) {
        snn <- crossprod(snn)
    }
    
    ## Pool information from neighboring cells to define spot 
    spots <- list()
    spots$metadata <- tibble(
        dplyr::select(cells$metadata, SpotID=CellID, LibraryID, x, y, x_img, y_img),
        ncells = colSums(snn), 
        area = as.numeric(snn %*% cells$metadata$area)
    )
    spots$intensity <- snn %*% as.matrix(cells$intensity)    
    return(spots)
}


plotFeatures_split <- function(
    data_mat, dim_df, features, split_by, nrow = 1, qlo = 0.05, qhi = 1, 
    order_by_expression = FALSE, pt_shape = 16, pt_size = 0.5, 
    no_guide = FALSE, .xlim = c(NA, NA), .ylim = c(NA, NA), color_high = muted("blue")
) {
    split(seq_len(length(split_by)), split_by) %>% map(function(idx) {
        plotFeatures(data_mat[, idx], dim_df[idx, ], features, nrow = 1, no_guide = TRUE)
    }) 
}



do_norm <- function(obj) {
    obj[['intensity_norm']] <- sweep(obj$intensity, 1, obj$metadata$area, '/')
    return(obj)
}


do_scale <- function(obj, z_thresh, within_batch=FALSE) {
    message("SCALE")
    if (within_batch) {
        .res <- split(seq_len(nrow(obj$metadata)), obj$metadata$LibraryID) %>% map(function(idx) {
            list(
                idx=idx,
                z=obj$intensity_norm[idx, ] %>% log1p() %>% scale() %>% pmax(-z_thresh) %>% pmin(z_thresh) 
            )
        }) 
        obj$z <- matrix(
            0, 
            nrow=nrow(obj$intensity_norm), 
            ncol=ncol(obj$intensity_norm),
            dimnames=dimnames(obj$intensity_norm)
        )
        
        for (x in .res) {
            obj$z[x$idx, ] <- x$z
        }        
    } else {
        obj$z <- obj$intensity_norm %>% log1p() %>% scale() %>% pmax(-z_thresh) %>% 
            pmin(z_thresh) %>% data.frame()
    }
    return(obj)
}

do_pca <- function(obj) {
    message('PCA')
    svd_res <- svd(t(as.matrix(obj$z)))
    V <- with(svd_res, sweep(v, 2, d, '*'))
    V <- data.frame(V)
    rownames(V) <- rownames(obj$z)
    colnames(V) <- paste0('PC', 1:ncol(V))
    U <- data.frame(svd_res$u)
    colnames(U) <- paste0('PC', 1:ncol(U))
    rownames(U) <- colnames(obj$intensity)

    obj$V <- V
    obj$loadings <- U    
    return(obj)
}


do_umap <- function(obj, embedding, resname) {
    message('UMAP')
    obj[[resname]] <- uwot::umap(obj[[embedding]], ret_extra = 'fgraph')
    colnames(obj[[resname]]$embedding) <- paste0('UMAP', 1:2)
    return(obj)
}


harmonize <- function(obj, var='LibraryID', theta=1) {
    fig.size(4, 6)
    obj$H <- HarmonyMatrix(
        obj$V, obj$metadata, var, theta=theta, 
        do_pca=FALSE, plot_convergence=TRUE,
        epsilon.cluster = -Inf, epsilon.harmony = -Inf, max.iter.cluster = 10, max.iter.harmony = 10
    )
    return(obj)
}


do_louvain <- function(obj, embedding, resolutions, pow=1) {
    adjmat <- obj[[embedding]]$fgraph
    diag(adjmat) <- 1
    for (i in seq_len(pow-1)) {
        adjmat <- adjmat %*% adjmat
    }
    obj$Clusters <- RunModularityClustering(adjmat, 1, resolutions)
    # obj$Clusters <- paste0('C', RunModularityClustering(adjmat, 1, resolutions))
    return(obj)
}


plot_heatmap <- function(markers, features, .scale=FALSE) {
    X <- markers %>% 
        # subset(feature %in% c('ASMA', 'CD146', 'CD3', 'CD31', 'CD45', 'CD90', 'PDPN', 'CD68')) %>% 
        # subset(!group %in% c('C22')) %>% 
        dplyr::select(SCORE=auc, feature, group) %>% 
        tidyr::spread(group, SCORE) %>% 
        tibble::column_to_rownames('feature') %>% 
        as('matrix') 
    
    if (.scale) X <- t(scale(t(X)))
    
    X <- X[features, ]
    Heatmap(X, column_names_rot = 45)    
}

do_markers <- function(obj) {
    obj$Markers <- apply(obj$Clusters, 2, function(y) {
        wilcoxauc(t(obj$z), y)
    })
    
    names(obj$Markers) <- colnames(obj$Clusters)
    return(obj)
}


do_subcluster <- function(obj, embedding, resolution, clusters, clusters_zoom, pow=1) {
    idx <- which(clusters %in% clusters_zoom)
    adjmat <- obj[[embedding]]$fgraph
    adjmat <- adjmat[idx, idx]
    diag(adjmat) <- 1
    for (i in seq_len(pow - 1)) {
        adjmat <- adjmat %*% adjmat
    }
    # obj$Clusters <- RunModularityClustering(adjmat, 1, resolutions)
    new_clusters <- RunModularityClustering(adjmat, 1, resolution)
    new_clusters <- unlist(new_clusters)
    new_clusters <- paste0(clusters[idx], '_', new_clusters)
    
    res <- clusters
    res[idx] <- new_clusters
    return(res)
}

nice_plot_coloc <- function(spots, obj, niches, types, color='red', alpha=.1, show_libname=TRUE) {
    plt <- spots$metadata %>% 
        ggplot(aes(y, -x)) + 
            geom_point(
                data = . %>% arrange(Niche_nice %in% niches), 
                shape = '.', aes(color = Niche_nice %in% niches), fill = NA
            ) + 
            theme_void() + 
            scale_fill_manual(values = c(color)) + 
            scale_color_manual(values = c('lightgrey', 'black')) + 
            guides(
                color = 'none', 
                fill = 'none', 
                alpha = 'none'
            ) + 
            facet_wrap(~LibraryID, scales='free') + 
            theme(strip.text = element_text(size = 12)) + 
            NULL

    if (!show_libname) {
        plt <- plt + theme(strip.text = element_blank())
    }
    
    plt + 
        geom_point(
            data = subset(obj$metadata, Subtype %in% types), 
            shape = 21, 
            size = 1, 
            aes(fill = Subtype),
            alpha = alpha
        ) + 
        NULL
    
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plot_biaxial <- function(dat, x, y, x0, y0, facet_by=NULL) {
    library(viridis)    
    dat <- data.frame(dat) 
    dat$density <- get_density(dat[, x], dat[, y], n = 100)     
    
    if (!is.null(facet_by))
        dat$VAR <- facet_by
    
    plt <- ggplot(dat, aes_string(x, y, color = 'density')) + 
        geom_point(shape = '.', alpha = .1, position = position_jitter(width = .02, height = .02)) + 
        geom_hline(yintercept = y0, linetype = 2, color = 'red') + 
        geom_vline(xintercept = x0, linetype = 2, color = 'red') + 
        # scale_color_gradient2_tableau() + 
        labs(x = x, y = y) + 
        scale_color_viridis() + 
        guides(color = 'none') + 
        NULL
    
    if (!is.null(facet_by))
        plt <- plt + facet_wrap(~VAR, scales='free')
    
    return(plt)
    
}


coloc_heatmap <- function(stats) {
    stats %>% 
        dplyr::mutate(estimate = exp(estimate)) %>% 
        dplyr::select(Niche_nice, Subtype, estimate) %>% 
        tidyr::spread(Niche_nice, estimate) %>% 
        tibble::column_to_rownames('Subtype') %>% 
        as.matrix() %>% 
        Heatmap(column_names_rot = 45)
}


get_coloc_stats <- function (obj) {    
    
    split(seq_len(nrow(obj$metadata)), obj$metadata$LibraryID) %>% 
        imap(function(.idx, .libname) {
            X <- with(obj$metadata[.idx, ], table(Niche_nice, Subtype)) %>% data.frame()

            map(unique(X$Subtype), function(Subtype_level) {
                model <- lme4::glmer(
                    Subtype_test ~ 1 + (1|Niche_nice), 
                    X %>% dplyr::mutate(Subtype_test = case_when(
                        Subtype == Subtype_level ~ 1, 
                        TRUE ~ 0
                    )), 
                    binomial,
                    weights = Freq
                )
                res <- data.frame(ranef(model)$Niche_nice) 
                colnames(res) <- c('beta')
                res <- tibble::rownames_to_column(res, 'Niche_nice')
                res$Subtype <- Subtype_level

                ## compute p values
                .sim <- arm::sim(model, 100)
                res$sigma <- apply(data.frame(.sim@ranef$Niche_nice), 2, sd)
                res$pval <- with(res, exp(pnorm(-beta/sigma, log.p = TRUE, lower.tail = TRUE))) ## 1-tailed
                return(res)
            }) %>% 
                bind_rows() %>% 
                cbind(LibraryID = .libname) %>% 
                dplyr::select(LibraryID, Subtype, Niche_nice, beta, sigma, pval) %>% 
                arrange(-beta/sigma)

        }) %>% 
        bind_rows() %>% 
        arrange(pval)
    
    # split(seq_len(nrow(obj$metadata)), obj$metadata$LibraryID) %>% 
    # imap(function(.idx, .libname) {
    #     X <- with(obj$metadata[.idx, ], table(Niche_nice, Subtype)) %>% data.frame()
    #     expand.grid(Niche = unique(X$Niche_nice), Type = unique(X$Subtype)) %>% 
    #         apply(1, function(vals) {
    #             glm(y ~ 1 + x, family = poisson, X %>% dplyr::mutate(y = Subtype == 
    #                 vals[["Type"]], x = Niche_nice == vals[["Niche"]]), 
    #                 weights = Freq) %>% broom::tidy() %>% subset(term == 
    #                 "xTRUE") %>% dplyr::mutate(Niche_nice = vals[["Niche"]], 
    #                 Subtype = vals[["Type"]]) %>% dplyr::select(-term) %>% 
    #                 dplyr::select(Niche_nice, Subtype, everything())
    #         }) %>%     
    #     bind_rows() %>% 
    #     dplyr::mutate(LibraryID = .libname)
    # }) %>% 
    #     bind_rows()
}

