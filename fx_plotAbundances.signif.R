plotAbundances.signif <- function (x, k = "meta20", 
                                   group_by = "condition", 
                                   shape_by = NULL, 
                                   col_clust = TRUE, 
                                   k_pal = CATALYST:::.cluster_cols, facets_ncol=4) {
  require(ggplot2)
  require(ggsignif)
  k <- CATALYST:::.check_k(x, k)
  CATALYST:::.check_cd_factor(x, group_by)
  CATALYST:::.check_cd_factor(x, shape_by)
  CATALYST:::.check_pal(k_pal)
  stopifnot(is.logical(col_clust), length(col_clust) == 1)
  # shapes <- CATALYST:::.get_shapes(x, shape_by)
  # 
  # if (is.null(shapes)) {
  #   shape_by <- NULL
  # }
  
  ns <- table(cluster_id = cluster_ids(x, k), sample_id = sample_ids(x))
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(df$sample_id, x$sample_id)
  
  for (i in c(shape_by, group_by)) {
    df[[i]] <- x[[i]][m]
  }
  
  ###
  df %<>% droplevels() %>% 
    dplyr::filter(! is.na(.[[group_by]]))
  ###
  
  comp.levels <- levels(df[[group_by]])
  comps <- combn(comp.levels, 2) %>%
    as.data.frame() %>%
    sapply(function(x) as.character(x), simplify = FALSE)
  
  
  p <- ggplot(df, aes(y =.data[["Freq"]],
                      x =.data[[group_by]],
                      color=.data[[group_by]],
                      fill=.data[[group_by]])
              ) + 
    geom_boxplot(alpha=0.5, outlier.shape=NA)
  
  if(!is.null(shape_by)){
    p <- p+geom_jitter(width=0.2,
                       mapping=aes(shape=.data[[shape_by]]))
  }else{
    p <- p+geom_jitter(width=0.2)
  }
    
    p <- p + labs(x = NULL, y = "Proportion [%]") + 
    theme_classic() + 
    geom_signif(comparisons=comps,
                color='black',
                size=0.5, textsize=5,
                step_increase=0.25,
                tip_length=0)+
    scale_color_manual(values=k_pal)+
    scale_fill_manual(values=k_pal)+
    facet_wrap('cluster_id', scales='free', ncol=facets_ncol)+
    scale_y_continuous(expand=expansion(mult=c(0, 0.1*length(comps))))+
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  return(p)
}
