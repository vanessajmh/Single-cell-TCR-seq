#Function to add a in the data with the clone frequency interval
intv <- function(x){
  if(x == 0) "None = 0"
  else if(x <= 0.0001) "Rare = 0.0001"
  else if(x <= 0.001) "Small = 0.001"
  else if(x <= 0.01) "Medium = 0.01"
  else if(x <= 0.1) "Large = 0.1"
  else "Hyperexpanded = 1"
}

#Function to process rds files from the utility dataset
process_rds<- function(x){
  
  Rare = 0.0001
  Small = 0.001
  Medium = 0.01
  Large = 0.1
  Hyperexpanded = 1
  
  filename=file_path_sans_ext(basename(x))
  sample<- readRDS(x)
  tcr<- sample@meta.data %>%
    drop_na() %>%
    filter(db.class =="singlet") %>% #remove estimated doublets by scDblFinder
    mutate(cell_subset=ifelse(grepl(paste(c("CD4","Tfh","Th1","Treg"), collapse="|"), PT.annot), "CD4+", 
                              ifelse(grepl("CD8",PT.annot ), "CD8+", "other"))) %>%
    mutate(Tissue=ifelse(grepl("LT",filename), "TME", 
                         ifelse(grepl("LB",filename), "Blood", "other"))) %>%
    group_by(orig.ident, cell_subset,Tissue, CTaa) %>%
    summarize(count=n()) %>%
    group_by(orig.ident, cell_subset,Tissue) %>%
    mutate(frequency_separated=count/sum(count)) %>% # for cd4 and cd8 separately
    select(orig.ident,CTaa,cell_subset,Tissue, count, frequency_separated )  %>%
    mutate(cloneSize=unlist(lapply(frequency_separated, intv))) %>%
    distinct()
  return(tcr)
}

#Function to calculate and plot the Jaccard and Morisita scores
plotDissimilarity <- function(x, 
                              cell_subset=c("CD4+","CD8+"),
                              method = c( "jaccard",  "morisita" ),
                              clustering=c("ward.D", "ward.D2", "single", "complete", "average", 
                                           "mcquitty" , "median", "centroid" ),
                              colorBy=NULL,
                              label_colors=NULL) {
  
  if (missing(x)) stop("x is missing.")
  if (is.null(colorBy)) stop("at least one group column name is expected.")
  if (is.null(method)) stop("a distance method is expected.")
  if (is.null(clustering)) stop("a clustering method is expected.")
  
  cell <- match.arg(cell_subset)
  methodChoice <- match.arg(method)
  clust<- match.arg(clustering)
 
  tmp <- x  %>% 
    filter(cell_subset==cell) %>%
    setDT()
  
  sNames <- as.character(unique(tmp$orig.ident))
  
  dat <- data.table::dcast(data = tmp, CTaa~orig.ident, value.var = "count", fun.aggregate = sum)
  
  simmat <- dat[, vegan::vegdist(t(.SD), method = methodChoice, diag = TRUE, upper = TRUE, binary = FALSE), 
                .SDcols=sNames] 
  
  groups <- x %>% 
    select(all_of(colorBy), orig.ident) %>% 
    filter(cell_subset==cell) %>%
    distinct() %>% 
    column_to_rownames("orig.ident")
  
  p <- ComplexHeatmap::pheatmap(as.matrix(simmat),
                                cluster_rows = TRUE, cluster_cols = TRUE, name = " ",  
                                treeheight_row = 0L, clustering_distance_rows = simmat, 
                                clustering_distance_cols = simmat,
                                annotation_col=groups, show_colnames=FALSE, labels_col = NA, 
                                annotation_colors = label_colors,column_title =cell,
                                show_rownames=FALSE, clustering_method = clust, silent = FALSE, fontsize =4)
  
  return(p)
}

#Function to concatenate the clustering outputs into one dataframe
concat_res<- function(list){
  # py_run_string("
  # my_set = list[[edges]]
  # my_list = list(my_set)
  # ")
  # my_list <- py$my_list
  
  tcr_results<- as.data.frame(list[["df"]])
  tcr_results$cluster<- as.character(tcr_results$cluster)
  
  motifs<- list[["motifs"]] %>% 
    rownames_to_column("cluster") %>%
    mutate(cluster=as.numeric(cluster)-1) 
  
  tcr_results<- merge(tcr_results, motifs, by="cluster")
  
}

#Function to identify TME-enriched clusters
enrich_cl<- function(clr, dataset, packages){
 
  library(dplyr)
  
  mat_metai<- dataset %>%
    dplyr::filter(cluster == clr) %>%
    dplyr::group_by(Type) %>%
    dplyr::summarize(in_cluster=n())
  
  other <- dataset %>%
    dplyr::filter(cluster != clr) %>%
    dplyr::select(junction_aa, Type) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Type) %>%
    dplyr::summarize(other=n())
  
  if(nrow(mat_metai)!=2){
    lev<- c( "TME", "Blood")
    name<-lev[!lev %in% as.character(mat_metai$Type) ]
    if(name != "TME")   mat_metai <- mat_metai %>% add_row("Type"="Blood", "in_cluster"=0, .before = 1)
  }
  res_stats<- data.frame()
  if(nrow(mat_metai)==2){
    mat<-merge(mat_metai, other)
    mat<- setNames(data.frame(t(mat[,-1])), mat[,1])
    mat<- mat[,c(2,1)]
    
    odds_ratio <- (mat[1,1] * mat[2,2]) / (mat[1,2] * mat[2,1])
    
    res_stats= mat %>%
      replace(., is.na(.), 0) %>%
      summarise(data = list(rstatix::fisher_test( 
        ., alternative="greater"
      ) ))%>%
      tidyr::unnest_wider(data) %>%   
      select(p) %>%
      mutate(cluster=clr) %>%
      mutate(odds=odds_ratio) %>%
      cbind(., mat_metai %>% mutate(count=1) %>% reshape2::dcast(., count~Type, value.var = 
                                                                   "in_cluster") %>% 
              select(-count))
  }
  return(res_stats)
}


#Function to parallelize the enrichment calculation
prepa_enrich<-function(res, packages){
  
  filter_cl<- res %>% group_by(cluster) %>% summarize(count=n()) %>% filter(count>2)
  chain=unique(res$chain)
  
  meta <- reshape2::melt(data_all[[chain]][,1:3]) %>%
    filter(value==1) %>%
    select(-value) %>%
    rename(Type=variable) %>%
    merge(res %>% filter(size>2), ., by.x='junction_aa', by.y=paste0("TR",chain))
  
  list_cl <- as.list(unique(meta$cluster))
  
  cl <- makeCluster(4, type = "SOCK", rscript_args = "--vanilla", useXDR = TRUE)

  clusterExport(cl , varlist = c("list_cl"), envir=environment())
 
  repList <- pblapply(cl = cl,
                      list_cl,
                      enrich_cl,
                      dataset=meta,
                      packages=packages)
  
  parallel::stopCluster(cl)
  
  results <- do.call(rbind, repList)
}

#Function to generate the network graphs
gen_graph <- function(df){
  
  chain=unique(df[["motifs"]]$chain)
  
  my_df <- data.frame(value = unlist(df[["edges"]]))
  my_df2 <- my_df %>% separate(value, c("A", "B"))
  my_df2<-merge(my_df2, res_dfs[[chain]][,c(1:2)], by.x="A", by.y="junction_aa")
  my_df2<-merge(my_df2, res_dfs[[chain]][,c(1:2)], by.x="B", by.y="junction_aa")
  my_df2<- my_df2 %>% mutate(from=paste0(A, "_", cluster.x))
  my_df2<- my_df2 %>% mutate(to=paste0(B, "_", cluster.y))
  
  graph<- my_df2 %>% 
    filter(cluster.x==cluster.y) %>% 
    filter(cluster.x %in% tcr_results[[chain]]$cluster | cluster.y %in% tcr_results[[chain]]$cluster) %>%
    merge(., res_dfs[[chain]] %>% select(-junction_aa) %>% distinct(), 
          by.x="cluster.x", by.y="cluster") %>%
    select(-cluster.x, -cluster.y, -A, -B) 
  
  keep <-tcr_results[[chain]]$cluster   
  metadata <- res_dfs[[chain]] %>%
    filter(cluster %in% keep ) %>%
    merge(., data[c(paste0("TR",chain), "Tissue")], by.x="junction_aa", by.y=paste0("TR",chain)) %>%
    mutate(size=scales::rescale(size, c(15,25)),
           color= ifelse(Tissue==  "Blood", mypal$Tissue[["Blood"]], mypal$Tissue[["TME"]] ),
           title=junction_aa,
           junction_aa=paste0(junction_aa, "_", cluster)) %>%
    distinct() 
  
  #make graph
  g <- igraph::graph_from_data_frame(graph, directed=F, vertices = metadata)
}

