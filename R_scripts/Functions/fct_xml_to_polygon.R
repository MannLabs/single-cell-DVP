########### BIAS xml file reader function  to sf-compatible data ############
if(!"tidyverse" %in% (.packages())){
  library(tidyverse) 
}

if(!"XML" %in% (.packages())){
  library(XML) 
}

if(!"sf" %in% (.packages())){
  library(sf) 
}

xml_to_polygon <- function(x){
  data <- xmlParse(x)
  
  xml_data <- xmlToList(data) # Read xml data
  
  n_elements <- as.numeric(xml_data$ShapeCount) # Extract number of shapes
  n_shapes <- n_elements
  
  ## XML format correct
  paste("Shape names correct:", n_elements == sum(grepl("Shape_", names(xml_data))))
  paste("XML length correct:", n_elements == length(names(xml_data)) - 8)
  
  ##  Put shapes into st_sf format one by one
  
  for(i in seq(from = 9, to = length(names(xml_data)), by = 1)){
    
    if(i == 9){polygons_list <- c()} # Empty variable toappend data as it comes
    
    # Wrangling
    tibble::enframe(unlist(xml_data[[i]])) %>%
      mutate(dim = str_replace(name, "_.*", "")) %>%
      dplyr::filter(dim %in% c("X", "Y")) %>%
      mutate(value = as.numeric(value)) %>%
      mutate(name = as.numeric(str_replace(name, ".*_", ""))) %>%
      spread(dim, value) %>%
      column_to_rownames("name") -> pol_i
    
    
    pol_i_mat <- as.matrix((pol_i))
    pol_i_list <- st_polygon(list(rbind(pol_i_mat, pol_i_mat[1,]))) # Close the polygon, and put into correct data format
    
    polygons_list <- append(polygons_list, pol_i_list)
  }
  
  # Put into geom_sf readable format
  polygons_sfc = lapply(polygons_list, FUN = st_polygon)
  
  return(polygons_sfc)
}
  