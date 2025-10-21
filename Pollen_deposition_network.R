library(dplyr)


Total_deposition_22_24 <- read.csv("data/Total_deposition_22_24.csv")
Total_deposition_22_24$year <- substr(Total_deposition_22_24$ID, start = 1, stop = 2)


pol_matrix <- Total_deposition_22_24[, -c(1, ncol(Total_deposition_22_24))]  # Remove ID and year columns

pol_matrix_22 <- pol_matrix[Total_deposition_22_24$year == "22", ]
pol_matrix_23 <- pol_matrix[Total_deposition_22_24$year == "23", ]
pol_matrix_24 <- pol_matrix[Total_deposition_22_24$year == "24", ]
pol_matrix_24[is.na(pol_matrix_24)] <- 0


for (i in 22:24) {
  
  pol_matrix_year <- pol_matrix[Total_deposition_22_24$year == i, ]
  pol_matrix_year[is.na(pol_matrix_year)] <- 0
  
  pol_matrix_year <- pol_matrix_year %>% mutate_at(1:length(pol_matrix_year), as.numeric)
  
  # Create a bipartite graph from the incidence matrix
  bg <- graph_from_biadjacency_matrix(pol_matrix_year)
  
  # Generate bipartite projections
  proj <- bipartite_projection(bg, multiplicity = TRUE)
  
  # Detect communities using a clustering algorithm
  communities <- cluster_louvain(proj$proj1)
  
  # Assign cluster membership as a node color
  V(proj$proj1)$color <- membership(communities)
  
  # Define the layout for the circle
  l <- layout_in_circle(proj$proj1)
  
  # Define vertex size by degree
  deg <- centr_degree(proj$proj1, mode="all")
  V(proj$proj1)$size <- 5 * sqrt(deg$res)
  
  # Shorten species names to the desired format (e.g., nam_nam)
  V(proj$proj1)$name <- sapply(
    V(proj$proj1)$name, 
    function(x) {
      parts <- strsplit(x, "_")[[1]]
      paste0(substr(parts[1], 1, 3), "_", substr(parts[2], 1, 3))
    }
  )
  
  # Calculate label positions based on vertex coordinates
  label_pos <- apply(l, 1, function(coords) {
    angle <- atan2(coords[2], coords[1])  # Calculate angle
    list(
      x = coords[1] + 0.2 * cos(angle),
      y = coords[2] + 0.2 * sin(angle)
    )
  })
  
  # Extract label x and y positions
  label_x <- sapply(label_pos, function(pos) pos$x)
  label_y <- sapply(label_pos, function(pos) pos$y)
  
  # Plot the graph with adjusted labels
  plot(
    proj$proj1,
    layout = l,
    edge.width = E(proj$proj1)$weight*0.4,
    vertex.label = NA,  # Suppress default labels
    vertex.size = V(proj$proj1)$size*0.5,
    vertex.color = zelena,
    edge.color = modra
  )
  
}
