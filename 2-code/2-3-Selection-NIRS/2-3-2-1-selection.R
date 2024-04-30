rm(list = ls())
ini = 6
set.seed(ini)

#### > 0. PACKAGES AND FUNCTIONS ####
library(tidyverse)
library(magrittr)
library(stringr)
library(purrr)
library(FactoMineR)
library(factoextra)
np.split.Multiblock.DX <- function (file = NULL, wd = NULL, export = FALSE){
  # intro ####
  if(is.null(wd) == TRUE) wd <- getwd()
  if((("samples" %in% dir(wd))==FALSE)&(export == TRUE)) dir.create("samples")
  if (is.null(file))  stop("You must provide a file name")
  
  # start import ####
  lines <- readLines(file)
  block_pat <- "##BLOCKS=(.*)"
  bc <- grep(block_pat, lines)
  bc <- sub(block_pat, "\\1", lines[bc])
  bc <- as.integer(bc)
  st_pat <- "##TITLE="
  st <- grep(st_pat, lines)
  st <- st[-1]
  end_pat <- "##END="
  end <- grep(end_pat, lines)
  end <- end[-length(end)]
  if (length(st) != length(end)) stop("Block starts and stops did not match")
  nb <- length(st)
  if (nb != bc) stop("Block count in link block did not match the number of blocks found")
  
  # block decomposition ####
  blocks <- vector("list", nb)
  for (i in 1:nb) blocks[[i]] <- (lines[st[i]:end[i]])
  
  fnames <- sapply(st,function(X,vec,sub_ch) gsub(sub_ch, "\\1", vec[X]), vec = lines, sub_ch = st_pat)
  fnames <- str_squish(fnames) %>% str_replace_all(" ","_") %>% str_replace_all("\\+","p") %>% str_replace_all("\\-","m")
  
  if (anyDuplicated(fnames))warning("Duplicated sample names found.\n\t\tLater samples will overwrite earlier samples\n\t\tunless you edit the original multiblock file.")
  names(blocks) <- paste0("S_",fnames)
  if(export == TRUE) for (i in 1:nb) writeLines(blocks[[i]], paste0(wd,"/samples/",fnames[i],".dx"))
  invisible(blocks)
}
np.read.block <- function (jdx){
  sstt <- grep( "^\\s*##XYDATA\\s*=\\s*\\(X\\+\\+\\(Y\\.\\.Y\\)\\)$", jdx)
  send <- grep("^\\s*##END\\s*=", jdx)
  
  metadata <- jdx[1:(sstt[1]-1)]
  
  Format <- c("metadata", "XYY")
  FirstLine <- c(1, sstt)
  LastLine <- c(sstt[1] - 1, send)
  
  DF <- data.frame(Format, FirstLine, LastLine, stringsAsFactors = FALSE)
  keep_lines <- (1+DF[2,2]):(DF[2,3]-2)
  fmr <- sapply(jdx[keep_lines],str_split,pattern = " ", simplify = TRUE) %>% sapply(as.numeric) %>% t()
  rownames(fmr) <- NULL
  
  VL <- list()
  VL$Dataguide <- DF
  VL$Metadata <- metadata
  VL$spectra <-  cbind(fmr[,1],rowMeans(fmr[,-1])) %>% as.data.frame()
  colnames(VL$spectra) <- c("wl","int")
  
  fmr <- metadata[grep("YFACTOR",metadata)] %>% str_split(" ")
  VL$spectra$int <- as.numeric(fmr[[1]][2])*VL$spectra$int
  
  fmr <- (grep("CONCENTRATIONS",metadata)+1):(grep("DELTAX",metadata)-1)
  VL$concentration <- str_remove_all(metadata[fmr],"\\(") %>%
    str_remove("\\)") %>% str_split(",",simplify = TRUE)
  fmr <- paste0(VL$concentration[,1],"_",VL$concentration[,3])
  VL$concentration <- as.numeric(VL$concentration[,2])
  names(VL$concentration) <- fmr
  return(VL)
}
np.export.date <- function(spl){
  spl$name <- str_split(spl$Metadata[grep("##TITLE",spl$Metadata)],"= ",simplify = TRUE)[,2] %>% str_trim() %>% str_squish()
  spl$date <- str_split(spl$Metadata[grep("##DATE",spl$Metadata)],"= ",simplify = TRUE)[,2] %>% str_trim() %>% str_squish()
  spl$time <- str_split(spl$Metadata[grep("##TIME",spl$Metadata)],"= ",simplify = TRUE)[,2] %>% str_trim() %>% str_squish()
  return(spl)
}
np.import <- function(name = "sample.dx"){
  np.split.Multiblock.DX(file = name) %>% lapply(np.read.block) %>% lapply(np.export.date)
}

#### > 1. Data ####
df = np.import('1-data/1-2-NIRS/MyDiv_subset1.dx')
df.1 = 
  df |> 
  purrr::map_df(
    .f = ~{
      .x$spectra |> 
        dplyr::mutate(sample = paste0(.x$name, "-",.x$time))
    }
  )

df.large = 
  pivot_wider(df.1, 
              names_from = wl,
              values_from = int)

# Plot all check 
p.measures = 
  ggplot(df.1 , 
       aes(x = wl, y = int, color = sample,
           alpha = .001)) + 
  geom_line() + 
  labs(x = "Length wave", y = "Measurement") +
  scale_x_reverse() + 
  theme_bw() + 
  theme(legend.position = 'none')

p.measures

ggsave(p.measures, filename = '3-plots/2-3-2-NIRS-measurements.png')

#### > 2. Sample selection #####
#### >> 2.1 PCA ####
n.tot.1 = 100
n.tot.2 = 252

# PCA measure
pca.nirs = PCA(df.large[,-1], 
               scale.unit = T,
               ncp = 10, 
               graph = F)

# Eigen values
p.eigenvalues = 
  fviz_eig(pca.nirs, addlabels = TRUE, ylim = c(0, 50))

p.eigenvalues

ggsave(p.eigenvalues, 
       filename = '3-plots/2-3-2-eigenvalues.png')

# Variables effects
pca.var <- get_pca_var(pca.nirs)

# fviz_pca_var(pca.nirs, col.var = "black")
# fviz_contrib(pca.nirs, choice = "var", axes = 1, top = 100)
# fviz_contrib(pca.nirs, choice = "var", axes = 2, top = 100)

pca.var = get_pca_var(pca.nirs)
pca.var.c = pca.var$contrib |> 
  data.frame()

pca.var.c$id = 
  unlist(rownames(pca.var.c)) 

pca.var.c = pivot_longer(pca.var.c, 
                         1:10, 
                         names_to = 'dim', 
                  values_to = 'loading')

p.lw.contribution = 
  ggplot(data = pca.var.c, 
       aes(x = as.numeric(id), 
           y = loading, color = dim)) + 
  geom_line() + 
  labs(x = 'Lengthwave', y = 'Loadings',
       color = 'Dimention') +
  theme_bw()

p.lw.contribution

ggsave(p.lw.contribution, 
       filename = '3-plots/2-3-2-loadings.png')

# Distribution across axes
p.point.distribution = 
  ggpubr::ggarrange(
    fviz_pca_ind(pca.nirs, geom="point"),
    fviz_pca_ind(pca.nirs, geom="point", axes = c(1,3)),
    fviz_pca_ind(pca.nirs, geom="point", axes = c(2,3)),
    ncol = 3
  )

p.point.distribution

ggsave(p.point.distribution, 
       filename = '3-plots/2-3-2-PCA-individuals.png')

# Species effect
fviz_pca_ind(pca.nirs,
             geom.ind = "point", 
             col.ind = 
               df.large$sample |> 
               str_split('-') |> 
               map_chr(.f = ~{.x[3]}), # color by groups
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Species"
)

# Month effect
fviz_pca_ind(pca.nirs,
             geom.ind = "point", 
             col.ind = 
               df.large$sample |> 
               str_split('-') |> 
               map_chr(.f = ~{.x[4]}),
             addEllipses = TRUE,
             legend.title = "Species")

#### >> 2.2 kmeans ####
library(ClusterR)
pca.ind = get_pca_ind(pca.nirs)$coord[,1:3]
wss <- function(k) {
  kmeans(pca.ind, k, nstart = 10)$tot.withinss
}
k.values <- 1:15
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


nb.k.selected = 4

k.means.nirs = kmeans(pca.ind, nb.k.selected, nstart = 10)

cluster.ind = 
  pca.ind |>
  data.frame() |> 
  mutate(cluster = k.means.nirs$cluster |> 
           unname() |> 
           factor())

p.cluster = 
  ggpubr::ggarrange(
    ggplot(data = cluster.ind, aes(x = Dim.1, y = Dim.2, color = cluster)) +
      geom_point() + 
      theme_bw(),
    ggplot(data = cluster.ind, aes(x = Dim.1, y = Dim.3, color = cluster)) +
      geom_point() + 
      theme_bw(),
    ggplot(data = cluster.ind, aes(x = Dim.2, y = Dim.3, color = cluster)) +
      geom_point() + 
      theme_bw(),
    ncol = 3, common.legend = T
  )

p.cluster


ggsave(p.cluster, filename = '3-plots/2-3-2-cluster.png')

library(plotly)

plot_ly(x=cluster.ind$Dim.1, 
        y=cluster.ind$Dim.2, 
        z=cluster.ind$Dim.3, 
        type="scatter3d", 
        mode="markers", 
        color=cluster.ind$cluster)

#### >> 2.3 Second kmeans for 100 samples ####
cluster.2.ind = 
  map_df(.x = 1: nb.k.selected, 
       .f = ~{
         d = cluster.ind |> 
           filter(cluster == .x)
         
         k.sub = kmeans(d[1:2], 
                        round(n.tot.1/nb.k.selected), 
                        nstart = 10)
         c.ind = 
           d |>
           data.frame() |> 
           mutate(c.2 = k.sub$cluster |> unname())
       })

cluster.2.ind = 
  cluster.2.ind |>
  mutate(cc = paste0(cluster, '-', c.2))

# Random selection in cluster
df.selection.100 = 
  map_df(
    .x = unique(cluster.2.ind$cc),
    .f = ~{
      cluster.2.ind |> 
        filter(cc == .x) |> 
        sample_n(1)
    })

p.selection.100 = 
  ggpubr::ggarrange(
    ggplot(data = cluster.2.ind, aes(x = Dim.1, y = Dim.2)) +
      geom_point(data = cluster.2.ind, aes(x = Dim.1, y = Dim.2)) +
      geom_point(data = df.selection.100, aes(x = Dim.1, y = Dim.2), color = 'red') +
      theme_bw(),
    ggplot(data = cluster.2.ind, aes(x = Dim.1, y = Dim.3)) +
      geom_point(data = cluster.2.ind, aes(x = Dim.1, y = Dim.3)) +
      geom_point(data = df.selection.100, aes(x = Dim.1, y = Dim.3), color = 'red') +
      theme_bw(),
    ggplot(data = cluster.2.ind, aes(x = Dim.2, y = Dim.3)) +
      geom_point(data = cluster.2.ind, aes(x = Dim.2, y = Dim.3)) +
      geom_point(data = df.selection.100, aes(x = Dim.2, y = Dim.3), color = 'red') +
      theme_bw(),
    ncol = 3, common.legend = T
  )

p.selection.100

ggsave(p.selection.100, filename = '3-plots/2-3-2-selected-100.png')

# Joint data 
df.100 = 
  df.large |>
  select(sample) |> 
  add_column(as.data.frame(pca.ind)) |>
  right_join(df.selection.100 |> 
               select(-c.2), 
             by = c('Dim.1', 'Dim.2', 'Dim.3'))

write_csv(df.100, file = '1-data/2-3-2-selection-100.csv')

#### >> 2.4 selections 250 ####
cluster.3.ind = 
  map_df(.x = 1: nb.k.selected, 
         .f = ~{
           d = cluster.ind |> 
             filter(cluster == .x)
           
           k.sub = kmeans(d[1:2], 
                          round(n.tot.2/nb.k.selected), 
                          nstart = 10)
           c.ind = 
             d |>
             data.frame() |> 
             mutate(c.3 = k.sub$cluster |> unname())
         })

cluster.3.ind = 
  cluster.3.ind |>
  mutate(cc = paste0(cluster, '-', c.3))

# Random selection in cluster
df.selection.250 = 
  map_df(
    .x = unique(cluster.3.ind$cc),
    .f = ~{
      d = cluster.3.ind |> 
        filter(cc == .x)
      d.j = inner_join(d, df.100 |>
                         select('Dim.1', 'Dim.2', 'Dim.3'), 
                       by = c('Dim.1', 'Dim.2', 'Dim.3'))
      if(nrow(d.j)>0){d.j}else{sample_n(d,1)}
    })

p.selection.250 = 
  ggpubr::ggarrange(
    ggplot(data = cluster.3.ind, aes(x = Dim.1, y = Dim.2)) +
      geom_point(data = cluster.3.ind, aes(x = Dim.1, y = Dim.2)) +
      geom_point(data = df.selection.250, aes(x = Dim.1, y = Dim.2), color = 'red') +
      theme_bw(),
    ggplot(data = cluster.3.ind, aes(x = Dim.1, y = Dim.3)) +
      geom_point(data = cluster.3.ind, aes(x = Dim.1, y = Dim.3)) +
      geom_point(data = df.selection.250, aes(x = Dim.1, y = Dim.3), color = 'red') +
      theme_bw(),
    ggplot(data = cluster.3.ind, aes(x = Dim.2, y = Dim.3)) +
      geom_point(data = cluster.3.ind, aes(x = Dim.2, y = Dim.3)) +
      geom_point(data = df.selection.250, aes(x = Dim.2, y = Dim.3), color = 'red') +
      theme_bw(),
    ncol = 3, common.legend = T
  )

p.selection.250

ggsave(p.selection.250, filename = '3-plots/2-3-2-selected-250.png')

# Joint data 
df.250 = 
  df.large |>
  select(sample) |> 
  add_column(as.data.frame(pca.ind)) |>
  right_join(df.selection.250 |> 
               select(-c.3), 
             by = c('Dim.1', 'Dim.2', 'Dim.3'))

write_csv(df.250, file = '1-data/2-3-2-selection-250.csv')

#### > 3. Selection control ####
#### >> 3.1 Control match selection 100 and 250 ####
sum(df.100$sample %in% df.250$sample)


#### >> 3.2 Check groups balance ####

# Should not be selected: '15-AE2-Ca-Mar'
a = unique(df.large$sample)
a = a[grepl('15-AE2-Ca-Mar', a)]
a %in% df.100$sample
a %in% df.250$sample

df.large$species = 
  df.large$sample |> 
  str_split('-') |> 
  map_chr(.f = ~{.x[3]})

df.large$plot = 
  df.large$sample |> 
  str_split('-') |> 
  map_chr(.f = ~{paste(.x[1],'-',.x[2])})

df.large$month = 
  df.large$sample |> 
  str_split('-') |> 
  map_chr(.f = ~{.x[4]})

df.100$species = 
  df.100$sample |> 
  str_split('-') |> 
  map_chr(.f = ~{.x[3]})

df.100$plot = 
  df.100$sample |> 
  str_split('-') |> 
  map_chr(.f = ~{paste(.x[1],'-',.x[2])})

df.250$plot = 
  df.250$sample |> 
  str_split('-') |> 
  map_chr(.f = ~{paste(.x[1],'-',.x[2])})

df.100$month = 
  df.100$sample |> 
  str_split('-') |> 
  map_chr(.f = ~{.x[4]})

df.design = readxl::read_xlsx(path ='1-data/1-0-design.xlsx') |> 
  select(plot = plotName, TSR = Tree_sp_rich, myc = Myc_type) |> 
  mutate(treatment = paste0(TSR, '-', myc)) |> 
  mutate(plot = str_replace_all(plot, '\\.', ' - ')) |> 
  group_by(plot, treatment) |> 
  summarise()

df.100 = df.100 |> 
  left_join(df.design, by = 'plot')

df.large = df.large |> 
  left_join(df.design, by = 'plot')

dist.species.tot = 
  df.large |> 
  group_by(species) |> 
  summarise(n = n())

dist.plot.tot = 
  df.large |> 
  group_by(treatment) |> 
  summarise(n = n())

dist.month.tot = 
  df.large |> 
  group_by(month) |> 
  summarise(n = n())

dist.species = 
  df.100 |> 
  group_by(species) |> 
  summarise(n = n())

dist.plot = 
  df.100 |> 
  group_by(treatment) |> 
  summarise(n = n())

dist.month = 
  df.100 |> 
  group_by(month) |> 
  summarise(n = n())


p.selection = 
  ggpubr::ggarrange(
  ggplot(data = dist.species, aes(x = species, y = n)) + 
    geom_bar(data = dist.species, aes(x = species, y = n/nrow(df.100)),
             stat = 'identity', fill = 'darkred') + 
    geom_bar(data = dist.species.tot, aes(x = species, y = n/nrow(df.large)),
             stat = 'identity', fill = NA, color = 'black') + 
    labs(y = "Proportion") + 
    theme_bw(),
  
  ggplot(data = dist.month, aes(x = month, y = n)) + 
    geom_bar(data = dist.month, aes(x = month, y = n/nrow(df.100)),
             stat = 'identity', fill = 'darkred') + 
    geom_bar(data = dist.month.tot, aes(x = month, y = n/nrow(df.large)),
             stat = 'identity', fill = NA, color = 'black') + 
    labs(y = "Proportion") + 
    theme_bw(),
  
  ggplot(data = dist.plot, aes(x = treatment, y = n)) + 
    geom_bar(data = dist.plot, aes(x = treatment, y = n/nrow(df.100)),
             stat = 'identity', fill = 'darkred') + 
    geom_bar(data = dist.plot.tot, aes(x = treatment, y = n/nrow(df.large)),
             stat = 'identity', fill = NA, color = 'black') + 
    labs(y = "Proportion") + 
    theme_bw(),
  nrow = 3
)

p.selection

ggsave(p.selection, filename = '3-plots/2-3-2-selection-control.png')

L = list(
  dist.m = left_join(dist.month, dist.month.tot, by = 'month') |> mutate(d = n.x - n.y),
  dist.s = left_join(dist.species, dist.species.tot, by = 'species') |> mutate(d = n.x - n.y),
  dist.p = left_join(dist.plot, dist.plot.tot, by = 'treatment') |> mutate(d = n.x - n.y)
)

l = c(ini, mean(L$dist.m$d),mean(L$dist.s$d),mean(L$dist.p$d))

df.selected = 
  df.100 |>
  select(sample) |> 
  mutate(select.100 = 'Yes') |> 
  full_join(
    df.250 |> 
      select(sample) |>
      mutate(select.250 = 'Yes'),
    by = 'sample', 
  )

writexl::write_xlsx(df.selected, '1-data/1-2-NIRS/df-selected.xlsx')
