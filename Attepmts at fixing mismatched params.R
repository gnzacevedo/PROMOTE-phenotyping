fcs.files <- dir(data.dir, pattern='.fcs', ignore.case=TRUE, full.names=TRUE, recursive=TRUE)

fcs.files.paths <- data.frame('File'=basename(fcs.files),
                              'Path'=fcs.files)

h <- read.FCSheader(fcs.files)

names(h) <- basename(names(h))

# get the parameter names
get.Pnames <- function(x){
  keys <- grep('P[[:digit:]]+N', names(x))
  vals <- x[keys]
  return(vals)
}

p.names <- lapply(h, get.Pnames)

param.df <- data.frame(p.names) %>%
  t() %>%
  as.data.frame()

panel.opts <- param.df %>%
  distinct()


sample.grps <- list()

for(i in 1:nrow(panel.opts)){
  sample.grps[[i]] <- panel.opts[i,] %>%
    left_join(param.df %>% 
                rownames_to_column(var='Sample')) %>%
    pull(Sample)
}

sample.grps


fcs.files1 <- fcs.files.paths %>% 
  dplyr::filter((File %>% gsub(' ', '\\.', .)) %in% sample.grps[[1]]) %>%
  pull(Path)

fs1 <- read.flowSet(fcs.files1, transformation=FALSE, alter.names=TRUE, which.lines=2e3)

fcs.files2 <- fcs.files.paths %>% 
  dplyr::filter((File %>% gsub(' ', '\\.', .)) %in% sample.grps[[2]]) %>%
  pull(Path)

fs2 <- read.flowSet(fcs.files2, transformation=FALSE, alter.names=TRUE, which.lines=2e3)

fs <- rbind2(fs1, fs2)