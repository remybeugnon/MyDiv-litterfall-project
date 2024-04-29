rm(list = ls())
set.seed(123)
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

df = np.import('1-data/1-2-NIRS/MyDiv.dx')
# df$S_48mE4mBemMay
  
df.1 = 
  df |> 
  purrr::map_df(
    .f = ~{
      .x$spectra |> 
        dplyr::mutate(sample = paste0(.x$name, "-",.x$time))
    }
  )

df.large = 
  pivot_wider(df.1, names_from = wl, values_from = int)
# Plot all check 
ggplot(df.1 , 
       aes(x = wl, y = int, color = sample,
           alpha = .001)) + 
  geom_line() + 
  theme_bw() + 
  theme(legend.position = 'none')

#### > 2. Sample selection #####
#### >> 2.1 PCA ####
n.tot = 100
# PCA measure
pca.nirs = PCA(df.large[,-1], 
               scale.unit = T,ncp = 8, 
               graph = F)

# Eigen values
fviz_eig(pca.nirs, addlabels = TRUE, ylim = c(0, 50))

# Variables effects
pca.var <- get_pca_var(pca.nirs)
fviz_pca_var(pca.nirs, col.var = "black")

fviz_contrib(pca.nirs, choice = "var", axes = 1, top = 100)
fviz_contrib(pca.nirs, choice = "var", axes = 2, top = 100)
pca.var = get_pca_var(pca.nirs)
df = pca.var$contrib |> 
  data.frame()
df$id = unlist(rownames(df)) 

df = pivot_longer(df, 1:8, names_to = 'dim', 
                  values_to = 'loading')
ggplot(data = df, 
       aes(x = as.numeric(id), y = loading, color = dim)) + 
  geom_line()

# Check 4500 - 5500 N 
# 
fviz_pca_ind(pca.nirs, 
             # pointsize = "cos2", 
             # pointshape = 21, 
             geom="point")

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

pca.ind = get_pca_ind(pca.nirs)$coord[,1:2]
wss <- function(k) {
  kmeans(pca.ind, k, nstart = 10 )$tot.withinss
}
k.values <- 1:15
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


nb.k.selected = 4
k.means.nirs = kmeans(pca.ind, nb.k.selected, nstart = 10)

fviz_cluster(k.means.nirs, geom = "point", data = pca.ind) +
  theme_bw()

cluster.ind = 
  pca.ind |>
  data.frame() |> 
  mutate(cluster = k.means.nirs$cluster |> unname())

#### >> 2.3 Second kmeans ####
cluster.2.ind = 
  map_df(.x = 1: nb.k.selected, 
       .f = ~{
         d = cluster.ind |> 
           filter(cluster == .x)
         
         k.sub = kmeans(d[1:2], 
                        round(n.tot/nb.k.selected), 
                        nstart = 10)
         c.ind = 
           d |>
           data.frame() |> 
           mutate(c.2 = k.sub$cluster |> unname())
       })

ggplot(data = cluster.2.ind,
       aes(x = Dim.1, y = Dim.2,
           color = factor(c.2))) +
  geom_point() +
  facet_grid(cols = vars(cluster))

cluster.2.ind = 
  cluster.2.ind |>
  mutate(cc = paste0(cluster, '-', c.2))

#### >> 2.4 Final selection ####
df.selection = 
  map_df(
    .x = unique(cluster.2.ind$cc),
    .f = ~{
      cluster.2.ind |> 
        filter(cc == .x) |> 
        sample_n(1)
    })

# Final selection
ggplot(data = cluster.2.ind, 
       aes(x = Dim.1, y = Dim.2)) + 
  geom_point(data = cluster.2.ind, 
             aes(x = Dim.1, y = Dim.2)) + 
  geom_point(data = df.selection,
             aes(x = Dim.1, y = Dim.2,
                 color = factor(cluster))) +
  theme_bw() + 
  theme(legend.position = 'none')
