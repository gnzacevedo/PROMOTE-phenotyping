sink('_log.txt')
cat(as.character(Sys.time()), ':: Install and load libraries \n')

chooseCRANmirror(ind=1)
pckgs.cran <- c('tidyverse','clustree','magrittr','tidygraph','ggraph', 'ggbeeswarm', 'parallelly')

pckgs.bioc <- c('CATALYST','FlowSOM','ConsensusClusterPlus','flowCore')

invisible(
  lapply(c(pckgs.cran, 'foreach','doParallel','BiocManager'), FUN=function(p) {
  
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies=TRUE, quiet=TRUE)
  }
  
  library(p, character.only=TRUE)
  return(NULL)
  })
)

invisible(
  lapply(pckgs.bioc, FUN=function(p){
    
    if(!requireNamespace(p, quietly=TRUE)){
      BiocManager::install(p, ask=FALSE)
    }
    
    library(p, character.only = TRUE)
    return(NULL)
  })
)

pckgs <- c(pckgs.cran, pckgs.bioc)

cat(as.character(Sys.time()), ':: Load up functions \n')
source('code/fx_clusterSaveSOM.R')
source('code/fx_plotMetaStability.R')


cat(as.character(Sys.time()), ':: Loading data \n')
sink()

sce <- readRDS('CD4_sce_scaled_UMAP.RDS')

sink('_log.txt', append=TRUE)
cat(as.character(Sys.time()), ':: Setting up parallelization \n')
sink()

ncores <- parallelly::availableCores()

cl <- makeCluster(ncores)
registerDoParallel(cl)

sink('_log.txt', append=TRUE)
cat(as.character(Sys.time()), ':: :: Running clustering \n')
    
sizes <- c(8,9)

fts <- c("KLRG1","CD25","S1PR1","CD95", "CD127",
         "CD69","PD1","CXCR6","CD161",
         "CCR4", "CCR6", "CXCR3","CD45RA","LAG3", "CX3CR1",
         "CXCR5","CCR7")

foreach(g=sizes, .packages=pckgs,
        .export=c('plotStabilityRibbons','clusterSaveSOM','plotStabilityBoxes')) %dopar% {
          
          #for(g in 10:20){          
          sink(file='_log.txt', append = TRUE)
          cat(as.character(Sys.time()), ':: :: Clustering with gridsize ',g, 'x', g, ' started \n')
          sink()
          
          x <- clusterSaveSOM(sce, features=fts, xdim=g, ydim=g, maxK=40,
                              som.rlen=25, seed=1)
          
          sink(file='_log.txt', append = TRUE)
          cat(as.character(Sys.time()), ':: :: :: plotting for gridsize ',g, 'x', g, '\n')
          sink()
          
          png(filename=paste0('ribbons_gridsize_', g,'.png'), width=1500, height=1200)
          rp <- plotStabilityRibbons(x, labels.repel=TRUE)+
            scale_x_continuous(breaks=seq.int(0, 40, 2))+
            theme(panel.grid.major.x=element_line())
          print(rp)
          dev.off()
          
          png(filename=paste0('tree_gridsize_', g,'.png'), width=1200, height=1200)
          tp <- clustree(x@metadata$cluster_codes, prefix='meta', 
                         node_size_range=c(2,10), node_alpha=0.5,
                         edge_width=0.5,  node_colour='sc3_stability')+
            scale_color_viridis_c(option='magma', direction = -1)
          print(tp)
          dev.off()
          
          png(filename=paste0('boxes_gridsize_', g,'.png'), width=1500, height=1200)
          bp <- plotStabilityBoxes(x)+
            theme(panel.grid.major.x=element_line())
          print(bp)
          dev.off()
          
          sink(file='_log.txt', append = TRUE)
          cat(as.character(Sys.time()), ':: :: :: :: Done with gridsize ',g, 'x', g, '\n')
          sink()
          
        }

sink('_log.txt', append=TRUE)
cat(as.character(Sys.time()), ':: Stopping parallel workers \n')
stopCluster(cl)

cat(as.character(Sys.time()), ':: Cleaning up \n')
rm(sce)
gc()
cat(as.character(Sys.time()), ':: Finished.')
sink()