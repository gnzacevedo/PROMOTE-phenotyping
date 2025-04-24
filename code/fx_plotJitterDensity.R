
plotJitterDensity <- function(x, fts, color.option='turbo', n.points=1e5){
  require(ggdist)
  require(tidyverse)
  require(magrittr)
  
  if('SingleCellExperiment' %in% class(x)){
    d <- assay(sce, 'exprs') %>% t() %>%
      as.data.frame() %>%
      select(all_of(fts)) %>%
      slice_sample(n=n.points)
  }
  
  if('Matrix' %in% class(x) | 'matrix' %in% class(x)){
    d <- as.data.frame(x) %>%
    select(all_of(fts)) %>%
    slice_sample(n=n.points)
  }
  
  d.l <- d %>%
    pivot_longer(cols = all_of(fts), names_to='Marker', values_to='Expression')
  
  
  # Function to calculate density for y values
  calculateDensity <- function(m) {
    dens <- density(m)
    point.dens <- sapply(m, FUN=function(n) dens$y[which.min(abs(dens$x-n))])
    # Return the density value at the mean of y_values
    return(point.dens)
  }
  
  # Calculate density for each group
  density.d.l <- apply(d, MARGIN=2, FUN=calculateDensity) %>%
    as.data.frame() %>%
    mutate(across(all_of(fts), function(x) (x-min(x))/(max(x)-min(x)))) %>%
    pivot_longer(cols=all_of(fts), values_to='point.density')
  
  plot.data <- d.l %>% mutate(point.density=density.d.l$point.density)
  summary.data <- plot.data %>%
    group_by(Marker) %>%
    summarize(median=median(Expression),
              q75=quantile(Expression, 0.75),
              q25=quantile(Expression, 0.25))
  
  plot <- plot.data %>% ggplot(aes(x=Marker, y=Expression))+
    geom_jitter(mapping=aes(color=point.density), size=0.5, alpha=0.7)+
    stat_pointinterval(color='black')+
    geom_point(data=.%>% group_by(Marker) %>% summarize(Expression=median(Expression)),
               shape=1, size=5)+
    #scale_color_gradientn(colors=heat.colors)+
    scale_color_viridis_c(option=color.option)+
    theme_classic()
  
  return(plot)
  
}
