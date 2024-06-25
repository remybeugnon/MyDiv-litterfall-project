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

# Functions to read .dx files
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
# First subset
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

# Second subset
# Need to do the same for subset 2 
df = np.import('1-data/1-2-NIRS/MyDiv_subset2_all.dx')
df.2 = 
  df |> 
  purrr::map_df(
    .f = ~{
      .x$spectra |> 
        dplyr::mutate(sample = paste0(.x$name, "-",.x$time))
    }
  )

df.large.2 = 
  pivot_wider(df.2, 
              names_from = wl,
              values_from = int) 

df.all = bind_rows(df.large |> mutate(sub = 1),
                   df.large.2|> mutate(sub = 2))
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

p.measures.2 = 
  ggplot(df.2 , 
         aes(x = wl, y = int, color = sample,
             alpha = .001)) + 
  geom_line() + 
  labs(x = "Length wave", y = "Measurement") +
  scale_x_reverse() + 
  theme_bw() + 
  theme(legend.position = 'none')

p.measures.2

# PCA measure
pca.nirs = PCA(df.all[,-c(1,189)], 
               scale.unit = T,
               ncp = 10, 
               graph = F)

p.point.distribution = 
  ggpubr::ggarrange(
    fviz_pca_ind(pca.nirs, geom="point"),
    fviz_pca_ind(pca.nirs, geom="point", axes = c(1,3)),
    fviz_pca_ind(pca.nirs, geom="point", axes = c(2,3)),
    ncol = 3
  )

p.point.distribution

# predict subset 2 
df.pred.sub.2 = 
  predict.PCA(pca.nirs, df.large.2[,-1])$coord



# distribution 
p.point.distribution = 
  ggpubr::ggarrange(
    # Dimention 1 - 2 
    fviz_pca_ind(pca.nirs, geom="point") + 
      geom_point(data = df.pred.sub.2,
                 aes(x = Dim.1, y = Dim.2),
                 color = 'red'),
    # Do the same for 1-3 and 2-3
    fviz_pca_ind(pca.nirs, geom="point", axes = c(1,3)) +
      geom_point(data = df.pred.sub.2,
                 aes(x = Dim.1, y = Dim.3),
                 color = 'red'), 
    fviz_pca_ind(pca.nirs, geom="point", axes = c(2,3)) +
      geom_point(data = df.pred.sub.2,
                 aes(x = Dim.2, y = Dim.3),
                 color = 'red'),
    ncol = 3
  )

p.point.distribution

d = df.pred.sub.2 |>
  data.frame() |>
  
  filter(Dim.1 < (-20) & Dim.2 < (-25))

# 3D plot 
library(plotly)

#convert matrix to data.frame
df.pred.sub.2.tab <- 
  predict.PCA(pca.nirs, df.all[,-1])$coord |> 
  as.data.frame() |> 
  mutate(sub = df.all$sub) |> 
  mutate(name = df.all$sample)|>
  
  filter(Dim.1 < (-20) & Dim.2 < (-25))

#3D plot subset 2
plot_ly(x=df.pred.sub.2.tab$Dim.1, 
        y=df.pred.sub.2.tab$Dim.2, 
        z=df.pred.sub.2.tab$Dim.3, 
        type="scatter3d", 
        mode="text", 
        color=df.pred.sub.2.tab$sub,
        alpha = .5, size = .1, 
        text = df.pred.sub.2.tab$name)%>% add_text(textfont = t, textposition = "top right")

# Need to merge subset 1 and 2 datasets in one

#convert matrix to data.frame
pca.ind = get_pca_ind(pca.nirs)$coord[,1:3]
pca.ind.tab <- as.data.frame(pca.ind)

#3D plot subset 1
plot_ly(x=pca.ind.tab$Dim.1, 
        y=pca.ind.tab$Dim.2, 
        z=pca.ind.tab$Dim.3, 
        type="scatter3d", 
        mode="markers", 
        color=pca.ind.tab$subset)


#3D plot subset 1 and 2
pca.ind.tab$Subset <- 1
df.pred.sub.2.tab$Subset <- 2

Subset.1.2 <- bind_rows(pca.ind.tab, df.pred.sub.2.tab)

plot_ly(x=Subset.1.2$Dim.1, 
        y=Subset.1.2$Dim.2, 
        z=Subset.1.2$Dim.3, 
        type="scatter3d", 
        mode="markers", 
        color=as.factor(Subset.1.2$Subset))


#Version 2 - 3D plot subset 1 and 2
plot_subsets <- plot_ly() %>%
  add_trace(
    x = pca.ind.tab$Dim.1,
    y = pca.ind.tab$Dim.2,
    z = pca.ind.tab$Dim.3,
    type = "scatter3d",
    mode = "markers",
    name = "Subset 1",
    marker = list(color = 'red')
  ) %>%
add_trace(
    x = df.pred.sub.2.tab$Dim.1,
    y = df.pred.sub.2.tab$Dim.2,
    z = df.pred.sub.2.tab$Dim.3,
    type = "scatter3d",
    mode = "markers",
    name = "Subset 2",
    marker = list(color = 'blue')
  ) 

plot_subsets


