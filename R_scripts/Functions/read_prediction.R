read_prediction <- function(x, bio_ID){
  d_predictors <- read_csv(x) %>%
    dplyr::select(-1) %>%
    dplyr::rename(Index = `Shape Index`) %>%
    gather(cluster, probability, 1:5) %>%
    mutate(cluster = as.numeric(cluster) + 1) %>%
    right_join(d_classes) %>%
    mutate(median_weighted = probability * median) %>%
    dplyr::select(Index, median_weighted, Protein) %>%
    group_by(Index, Protein) %>%
    mutate(int_weighted = sum(median_weighted)) %>%
    mutate(bio_ID = bio_ID) %>%
    left_join(img_meta) %>%
    drop_na(cell_ID)
  
  return(d_predictors)
}