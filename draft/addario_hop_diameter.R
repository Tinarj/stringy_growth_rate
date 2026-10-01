# Simulations for Addario-Berry hop diameter
library(cassowaryr)
library(igraph)

# Compute Addario-Berry quantities
compute_addario_metrics <- function(x, y) {
  
  sc <- scree(
    x,
    y,
    out.rm = FALSE,
    binner = NULL
  )
  
  mst <- cassowaryr:::gen_mst(
    del = sc$del,
    weights = sc$weights
  )
  
  mst_weights <- igraph::E(mst)$weight
  
  # Addario-Berry diameter:
  # longest path based on number of edges
  # Ignore Euclidean edge weights
  
  addario_path <- igraph::get_diameter(
    mst,
    directed = FALSE,
    weights = NA
  )
  
  addario_vertices <- as.integer(addario_path)
  
  addario_num_vertices <- length(addario_vertices)
  # Hop diameter = number of edges
  addario_hop_diameter <- addario_num_vertices - 1L
  
 
  # Euclidean length of the same hop-diameter path
  
  vertex_pairs <- as.vector(
    rbind(
      addario_vertices[-length(addario_vertices)],
      addario_vertices[-1]
    )
  )
  
  addario_edge_ids <- igraph::get_edge_ids(
    mst,
    vp = vertex_pairs,
    directed = FALSE
  )
  
  addario_path_length <- sum(
    mst_weights[addario_edge_ids]
  )
  
  
  data.frame(
    addario_hop_diameter = addario_hop_diameter,
    addario_num_vertices = addario_num_vertices,
    addario_path_length = addario_path_length
  )
}


# Run Addario-Berry simulations for different sample sizes
run_addario_simulation <- function(
    n_start,
    n_end,
    step,
    B = 100,
    seed = 403,
    file = "addario_results.csv") {
  
  n_values <- seq(
    from = n_start,
    to = n_end,
    by = step
  )
  
  total_runs <- length(n_values) * B
  
  all_results <- vector(
    "list",
    total_runs
  )
  
  counter <- 1L
  
  
  for (n in n_values) {
    
    message("Running n = ", n)
    
    for (sim in seq_len(B)) {
      
      # Same seed construction as my original simulations
      current_seed <- seed + counter
      
      set.seed(current_seed)
      
      # Reproduce the same independent uniform points
      x <- runif(
        n,
        min = 0,
        max = 1
      )
      
      y <- runif(
        n,
        min = 0,
        max = 1
      )
      
      metrics <- compute_addario_metrics(
        x,
        y
      )
      
      
      all_results[[counter]] <- data.frame(
        n = n,
        simulation_run = sim,
        seed = current_seed,
        metrics
      )
      
      
      counter <- counter + 1L
    }
  }
  
  
  results <- do.call(
    rbind,
    all_results
  )
  
  
  write.csv(
    results,
    file = file,
    row.names = FALSE
  )
  
  
  message(
    "Results saved to: ",
    file
  )
  
  
  results
}



results_addario_1 <- run_addario_simulation(
  n_start = 50,
  n_end = 175,
  step = 25,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "addario_uniform_50_175.csv"
  )
)


results_addario_2 <- run_addario_simulation(
  n_start = 200,
  n_end = 500,
  step = 50,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "addario_uniform_200_500.csv"
  )
)


results_addario_3 <- run_addario_simulation(
  n_start = 550,
  n_end = 800,
  step = 50,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "addario_uniform_550_800.csv"
  )
)


results_addario_4 <- run_addario_simulation(
  n_start = 850,
  n_end = 1000,
  step = 50,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "addario_uniform_850_1000.csv"
  )
)


results_addario_5 <- run_addario_simulation(
  n_start = 1100,
  n_end = 2000,
  step = 100,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "addario_uniform_1100_2000.csv"
  )
)


results_addario_5 <- run_addario_simulation(
  n_start = 3000,
  n_end = 5000,
  step = 1000,
  B = 100,
  seed = 403,
  file = paste0(
    "blog/stringy05_growth_rate/",
    "addario_uniform_3000_5000.csv"
  )
)
