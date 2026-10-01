# Simulations using Gaussian noise
library(cassowaryr)
library(igraph)

# Compute diameter, length, and ratio

compute_stringy_metrics <- function(x, y) {
  
  sc <- scree(x, y, out.rm = FALSE, binner = NULL)
  mst <- cassowaryr:::gen_mst(sc$del, sc$weights)
  
  diameter_length <- igraph::diameter(mst)
  
  total_length <- sum(igraph::E(mst)$weight)
  
  ratio <- diameter_length / total_length
  
  data.frame(
    diameter = diameter_length,
    length = total_length,
    ratio = ratio
  )
}




run_stringy_simulation <- function(n_start,
                                   n_end,
                                   step,
                                   B = 1000,
                                   seed = 403,
                                   file = "stringy05_results.csv") {
  
  
  n_values <- seq(n_start, n_end, by = step)
  
  all_results <- list()
  counter <- 1
  
  for (n in n_values) {
    
    message("Running n = ", n)
    
    for (sim in seq_len(B)) {
      
      current_seed <- seed + counter
      set.seed(current_seed)
      
      x <- rnorm(n, mean = 0, sd = 1)
      y <- rnorm(n, mean = 0, sd = 1)
      
      metrics <- compute_stringy_metrics(x, y)
      
      all_results[[counter]] <- data.frame(
        n = n,
        simulation_run = sim,
        diameter = metrics$diameter,
        length = metrics$length,
        ratio = metrics$ratio
      )
      
      counter <- counter + 1
    }
  }
  
  results <- do.call(rbind, all_results)
  
  write.csv(results, file = file, row.names = FALSE)
  
  return(results)
}


results <- run_stringy_simulation(
  n_start = 8000,
  n_end = 10000,
  step = 1000,
  B = 1000,
  seed = 403,
  file = "blog/stringy05_growth_rate/stringy05_7000_10000.csv"
)












# Simulations using uniform distribution
library(cassowaryr)
library(igraph)

compute_stringy_metrics <- function(x, y) {
  
  sc <- scree(x, y, out.rm = FALSE, binner = NULL)
  
  mst <- cassowaryr:::gen_mst(
    del = sc$del,
    weights = sc$weights
  )
  
  mst_weights <- igraph::E(mst)$weight
  
  # Total MST length: denominator of stringy05
  total_length <- sum(mst_weights)
  
  # Stringy05 numerator: weighted diameter of the MST
  stringy_path <- igraph::get_diameter(
    mst,
    directed = FALSE,
    weights = mst_weights
  )
  
  stringy_diameter <- igraph::diameter(
    mst,
    directed = FALSE,
    weights = mst_weights
  )
  
  stringy_num_vertices <- length(stringy_path)
  stringy_num_edges <- stringy_num_vertices - 1L
  
  # Penrose quantity: longest single MST edge
  penrose_edge_id <- which.max(mst_weights)
  penrose_diameter <- mst_weights[penrose_edge_id]
  
  # A single edge always contains two vertices
  penrose_num_edges <- 1L
  penrose_num_vertices <- 2L
  
  # Original stringy05 index
  stringy05 <- stringy_diameter / total_length
  
  data.frame(
    stringy_diameter = stringy_diameter,
    stringy_num_edges = stringy_num_edges,
    stringy_num_vertices = stringy_num_vertices,
    penrose_diameter = penrose_diameter,
    penrose_num_edges = penrose_num_edges,
    penrose_num_vertices = penrose_num_vertices,
    total_mst_length = total_length,
    stringy05 = stringy05
  )
}


# Run simulations for different sample sizes
run_stringy_simulation <- function(
    n_start,
    n_end,
    step,
    B = 1000,
    seed = 403,
    file = "stringy05_results.csv") {
  
  n_values <- seq(
    from = n_start,
    to = n_end,
    by = step
  )
  
  total_runs <- length(n_values) * B
  all_results <- vector("list", total_runs)
  
  counter <- 1L
  
  for (n in n_values) {
    
    message("Running n = ", n)
    
    for (sim in seq_len(B)) {
      
      current_seed <- seed + counter
      set.seed(current_seed)
      
      # Independent uniform points
      x <- runif(n, min = 0, max = 1)
      y <- runif(n, min = 0, max = 1)
      
      metrics <- compute_stringy_metrics(x, y)
      
      all_results[[counter]] <- data.frame(
        n = n,
        simulation_run = sim,
        seed = current_seed,
        metrics
      )
      
      counter <- counter + 1L
    }
  }
  
  results <- do.call(rbind, all_results)
  
  write.csv(
    results,
    file = file,
    row.names = FALSE
  )
  
  message("Results saved to: ", file)
  results
}




results_1 <- run_stringy_simulation(
  n_start = 50,
  n_end = 175,
  step = 25,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "stringy05_uniform_50_175.csv"
  )
)


results_2 <- run_stringy_simulation(
  n_start = 200,
  n_end = 500,
  step = 50,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "stringy05_uniform_200_500.csv"
  )
)


results_3 <- run_stringy_simulation(
  n_start = 550,
  n_end = 800,
  step = 50,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "stringy05_uniform_550_800.csv"
  )
)

results_4 <- run_stringy_simulation(
  n_start = 850,
  n_end = 1000,
  step = 50,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "stringy05_uniform_850_1000.csv"
  )
)


results_5 <- run_stringy_simulation(
  n_start = 1100,
  n_end = 2000,
  step = 100,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "stringy05_uniform_1100_2000.csv"
  )
)



results_6 <- run_stringy_simulation(
  n_start = 3000,
  n_end = 5000,
  step = 1000,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "stringy05_uniform_3000_5000.csv"
  )
)

