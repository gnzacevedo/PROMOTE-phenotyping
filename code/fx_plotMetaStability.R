plotStabilityRibbons <- function(sce, labels.repel=FALSE){
  require(tidygraph)
  require(ggraph)
  
  x <- clustree(sce@metadata$cluster_codes, prefix='meta', return='graph')
  
  nodes <- as.data.frame(activate(x, nodes)) %>%
    rename(name='node', 
           x='meta',
           y='sc3_stability')
  
  edges <- as.data.frame(activate(x, edges)) %>%
    mutate(from=paste0('meta',from_meta,'C',from_clust),
           to=paste0('meta', to_meta,'C',to_clust))
  
  graph <- tbl_graph(nodes=nodes, edges=edges, directed=TRUE) 
  
  lo <- create_layout(graph, layout=nodes %>% mutate(x=as.integer(as.character(x))))
  
  my.colors <- c('firebrick','seagreen','orange','darkorchid4','navy','hotpink',
                 'indianred','darkseagreen','khaki','lavender','steelblue','pink',
                 'darkred','yellowgreen','yellow','violet','skyblue','coral')
  
  plot <- ggraph(lo)+
    geom_edge_link(mapping=aes(width=in_prop, color=from_clust), alpha=0.5, show.legend=FALSE)+
    geom_node_point(mapping=aes(size=size), 
                    shape=16,
                    color='black', 
                    alpha=1)+
    geom_node_label(mapping=aes(label=cluster),
                   color='black', repel=labels.repel, 
                   box.padding=0.1, label.r=unit(0, 'mm'), 
                   fill=alpha('white'), alpha=0.7,
                   max.overlaps=10, force=0.5, force_pull=10)+
    geom_line(data=nodes %>% group_by(x) %>%
                summarize(y=mean(y)) %>%
                mutate(x=as.integer(x)),
              color='navy', linewidth=0.5, linetype='dashed',
              mapping=aes(x=x, y=y))+
    labs(x='Metaclusters', y='SC3 stability')+
    scale_size_continuous(range=c(3,16), name='Clusters in metacluster',
                          breaks=c(50, 150, 250))+
    scale_edge_width_continuous(range=c(0.2, 4))+
    scale_edge_color_manual(values=colorRampPalette(my.colors)(max(as.numeric(edges$from_clust))))+
    scale_color_manual(values=colorRampPalette(my.colors)(max(as.numeric(nodes$cluster))))+
    theme_classic()+theme(legend.position='bottom')
  
  return(plot)
}

plotStabilityBoxes <- function(sce){
  require(tidygraph)
  require(ggraph)
  
  x <- clustree(sce@metadata$cluster_codes, prefix='meta', return='graph') %>%
    activate(nodes) %>% as.data.frame()
  
  x %<>%
    group_by(meta) %>%
    mutate(total.clusters=sum(size)) %>%
    mutate(relative.size=size/total.clusters) %>%
    mutate(weighted.stability=sc3_stability*relative.size) %>%
    ungroup()
  
  plot <- x %>%
    ggplot(aes(x=meta, y=sc3_stability))+
    geom_line(data=.%>%
                group_by(meta) %>% 
                summarize(sc3_stability=min(sc3_stability)),
              inherit.aes = FALSE, mapping=aes(x=meta, y=sc3_stability, group=NA, color='Min-max'),
              linewidth=0.5)+
    geom_line(data=.%>%
                group_by(meta) %>% 
                summarize(sc3_stability=max(sc3_stability)),
              inherit.aes = FALSE, mapping=aes(x=meta, y=sc3_stability, group=NA, color='Min-max'),
              linewidth=0.5)+
    geom_boxplot(outlier.shape=NA)+
    geom_line(data=.%>%
                group_by(meta) %>% 
                summarize(sc3_stability=mean(sc3_stability)),
              inherit.aes = FALSE, mapping=aes(x=meta, y=sc3_stability, group=NA,
                                               color='Mean stability'),
              linewidth=1)+
    geom_line(data=.%>%
                group_by(meta) %>% 
                summarize(w.stability=sum(weighted.stability)/n_distinct(cluster)) %>%
                ungroup() %>%
                mutate(w.stability=(w.stability-min(w.stability))/(max(w.stability)-min(w.stability))),
              inherit.aes = FALSE, mapping=aes(x=meta, y=w.stability, group=NA,
                                               color='Weighted mean stability'),
              linewidth=0.5)+
    scale_color_manual(values=c('Mean stability'='firebrick', 
                                'Min-max'='indianred',
                                'Weighted mean stability'='seagreen'), guide='legend', name='')+
    scale_size_continuous(name='Clusters in metacluster')+
    geom_quasirandom(alpha=0.7, shape=16, color='steelblue', 
                     mapping=aes(size=size))+
    theme_classic()+theme(legend.position='bottom')
  
  return(plot)
}

highlight=c('meta22C16','meta22C18','meta22C19','meta22C20')
