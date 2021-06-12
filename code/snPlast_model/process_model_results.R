
rm(list = ls())

library(tidyverse)

setwd("../../data/simulations/raw")
files = list.files(pattern="*.csv")

m <- NULL

for (f in files){
   mat <- read.csv(f, header=FALSE)
   colnames(mat) <- seq(1, ncol(mat))
   m <- as_tibble(mat) %>%
     mutate(iter=seq(1, nrow(.))) %>%
     gather(key="trial", value="well", -iter) %>%
     mutate(ACh=stringr::str_extract(f, "(?<=ACh)\\d\\.\\d+"),
            DA=stringr::str_extract(f, "(?<=DA)\\d\\.\\d+"),
            wmax=stringr::str_extract(f, "(?<=wmax)\\d+"),
           trial=as.numeric(trial)) %>%
     arrange(iter, trial) %>%
     bind_rows(., m)
 }
 
saveRDS(m, "../simdat.rds")