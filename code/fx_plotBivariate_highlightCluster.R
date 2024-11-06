plotBivariates.highlightClusters <- function(x, dimpair, meta, koi, bg.evts=1e6){
  require(tidyverse)
  require(ggrepel)
  
  if(class(x)=='SingleCellExperiment'){
    umap.df <- DRdf(x, dr='UMAP', 
                    md=c(colnames(ei(sce)), names(cluster_codes(sce)), rownames(sce)))
  }else if(class(x)=='data.frame'){
    umap.df <- x
  }else{stop('Input is not SCE or data frame')}
  
  my.colors <- c('firebrick','seagreen','orange','darkorchid4','navy','hotpink',
                 'indianred','darkseagreen','khaki','lavender','steelblue','pink',
                 'darkred','yellowgreen','yellow','violet','skyblue','coral')
  
  
  if(length(my.colors) < length(koi)){
    my.colors <- colorRampPalette(my.colors)(n.colors)
  }
  
  h.val.x <- umap.df %>% pull(dimpair[1]) %>% sort(decreasing=TRUE) %>% .[1:200] %>% mean(.)/10
  h.val.y <- umap.df %>% pull(dimpair[2]) %>% sort(decreasing=TRUE) %>% .[1:200] %>% mean(.)/10
  
  centroids <- umap.df %>%
    dplyr::select(all_of(c(meta, dimpair))) %>%
    dplyr::filter(.[[meta]] %in% koi) %>%
    dplyr::group_by(!!sym(meta)) %>%
    summarize(across(all_of(dimpair), median))
  
  plot <- umap.df %>%
    ggplot(aes(x=.[[ dimpair[1] ]], 
               y=.[[ dimpair[2] ]]))+
    geom_point(data= .%>% slice_sample(n=bg.evts),
               mapping=aes(x=.data[[ dimpair [1] ]], 
                           y=.data[[ dimpair[2] ]]),
               size=0.1, color='gray60', alpha=0.2)+
    geom_point(data= umap.df %>% dplyr::filter(!!sym(meta) %in% koi),
               mapping=aes(color=!!sym(meta), 
                           x=.data[[ dimpair[1] ]], 
                           y=.data[[ dimpair[2] ]]),
               size=0.3, alpha=0.1)+
    geom_point(data=centroids, 
               mapping=aes(x=.data[[ dimpair[1] ]],
                           y=.data[[ dimpair[2] ]]),
               size=4, color='white', shape=16)+
    geom_point(data=centroids, 
               size=4,  shape=1,
               mapping=aes(color=!!sym(meta),
                           x=.data[[ dimpair[1] ]],
                           y=.data[[ dimpair[2] ]]))+
    geom_density2d(data=.%>% slice_sample(n=bg.evts), color='gray90',
                   h=c(h.val.x, h.val.y), n=200, bins=20,
                   mapping=aes(x=.data[[ dimpair[1] ]],
                               y=.data[[ dimpair[2] ]]))+
    geom_label_repel(data=centroids, 
                     size=4,
                     mapping=aes(label=!!sym(meta),
                                 color=!!sym(meta),
                                 x=.data[[ dimpair[1] ]],
                                 y=.data[[ dimpair[2] ]]))+
    scale_color_manual(values=my.colors)+
    theme_classic()+
    labs(x=dimpair[1], y=dimpair[2])+
    theme(aspect.ratio=1, legend.position='none')
  
  return(plot)
}