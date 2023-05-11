normalize_core_mean <- function(data, proteome_core, SA_incl){
  
  data %>%
    filter(label != "Dimethyl-n-0") %>%
    filter(cell_ID %in% SA_incl) %>%
    drop_na(int) %>%
    filter(Protein %in% proteome_core) %>%
    group_by(cell_ID) %>%
    summarise(median = median(int)) %>%
    ungroup() %>%
    summarise(median_total = mean(median)) %>%
    pull(median_total) -> proteome_core_median
  
  data  %>%
    filter(label != "Dimethyl-n-0") %>%
    filter(cell_ID %in% SA_incl) %>%
    drop_na(int) %>%
    filter(Protein %in% proteome_core) %>%
    group_by(cell_ID) %>%
    summarise(median_sample = median(int)) %>%
    ungroup() %>%
    mutate(median_total = proteome_core_median) %>%
    mutate(norm_factor = median_total / median_sample) %>%
    dplyr::select(cell_ID, norm_factor) -> proteome_core_norms
  
  data %>%
    filter(label != "Dimethyl-n-0") %>%
    filter(cell_ID %in% SA_incl) %>%
    drop_na(int) %>%
    left_join(proteome_core_norms) %>%
    mutate(int_core = int * norm_factor) %>%
    dplyr::select(-int, - int) %>%
    mutate(int_core = log2(int_core)) -> d_norm
  
  return(d_norm)
}
