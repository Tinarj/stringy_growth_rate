library(cassowaryr)
library(tidyverse)
library(knitr)

# Comparing weighted diameter version and hop diameter version of stringy05
compare_stringy <- function(x, y) {
  
  sc <- scree(x, y, binner = NULL, out.rm = FALSE)
  mst <- cassowaryr:::gen_mst(sc$del, sc$weights)
  
  w <- igraph::E(mst)$weight
  L <- sum(w)
  
  dw <- igraph::diameter(mst, weights = w)
  dh <- igraph::diameter(mst, weights = NA)
  
  sw <- dw / L
  sh <- dh / (length(x) - 1)
  
  data.frame(
    weighted_diameter = dw,
    hop_diameter = dh,
    mst_length = L,
    weighted_stringy = sw,
    hop_stringy = sh,
    difference = sw - sh
  )
}


# 2. Compare structured patterns ----------------------------------------

set.seed(1050)

n <- 500
t <- seq(-1, 1, length.out = n)
noise_sd <- 0.005

patterns <- list(
  "Degree 1 vs 2" = list(x = t, y = t^2),
  "Degree 2 vs 3" = list(x = t^2, y = t^3),
  "Sine wave" = list(x = t, y = sin(2 * pi * t))
)

pattern_data <- lapply(names(patterns), function(name) {
  p <- patterns[[name]]
  data.frame(
    structure = name,
    x = p$x + rnorm(n, sd = noise_sd),
    y = p$y + rnorm(n, sd = noise_sd)
  )
})

pattern_results <- do.call(rbind, lapply(pattern_data, function(d) {
  data.frame(structure = d$structure[1], compare_stringy(d$x, d$y))
}))

kable(pattern_results, digits = 4)

ggplot(do.call(rbind, pattern_data), aes(x, y)) +
  geom_point(size = 0.5) +
  facet_wrap(~ structure, scales = "free") +
  theme_bw()


# 3. Gaussian-noise simulation ------------------------------------------

simulate_stringy <- function(
    n_start = 25, n_end = 150, steps = 25,
    n_sim = 10, seed_start = 1050,
    save_csv = TRUE,
    csv_file = "stringy_noise_comparison.csv",
    append = TRUE
) {
  
  ns <- seq(n_start, n_end, by = steps)
  
  existing <- save_csv && append && file.exists(csv_file)
  
  if (existing) {
    previous <- read.csv(csv_file)
    seed <- max(seed_start, max(previous$seed) + 1)
    run_id <- max(previous$run_id) + 1
  } else {
    seed <- seed_start
    run_id <- 1
  }
  
  results <- vector("list", length(ns) * n_sim)
  i <- 1
  
  for (n in ns) {
    for (b in seq_len(n_sim)) {
      
      set.seed(seed)
      
      row <- data.frame(
        run_id = run_id,
        n = n,
        replicate = b,
        seed = seed,
        compare_stringy(rnorm(n), rnorm(n))
      )
      
      results[[i]] <- row
      
      if (save_csv) {
        write.table(
          row, csv_file, sep = ",",
          row.names = FALSE,
          col.names = !existing,
          append = existing
        )
        existing <- TRUE
      }
      
      seed <- seed + 1
      run_id <- run_id + 1
      i <- i + 1
    }
    
    message("Finished n = ", n, ": ", n_sim, " simulations")
  }
  
  if (save_csv) read.csv(csv_file) else do.call(rbind, results)
}


# 4. Run simulation ------------------------------------------------------

results <- simulate_stringy(
  n_start = 50,
  n_end = 200,
  steps = 25,
  n_sim = 1000,
  seed_start = 1050,
  csv_file = "blog/stringy05_growth_rate/stringy05_weighted_hop_50_200.csv",
  append = TRUE
)

results2 <- simulate_stringy(
  n_start = 250,
  n_end = 500,
  steps = 50,
  n_sim = 1000,
  seed_start = 1050,
  csv_file = "blog/stringy05_growth_rate/stringy05_weighted_hop_250_500.csv",
  append = TRUE
)

results3 <- simulate_stringy(
  n_start = 550,
  n_end = 750,
  steps = 50,
  n_sim = 1000,
  seed_start = 1050,
  csv_file = "blog/stringy05_growth_rate/stringy05_weighted_hop_550_750.csv",
  append = TRUE
)

results4 <- simulate_stringy(
  n_start = 800,
  n_end = 1000,
  steps = 100,
  n_sim = 1000,
  seed_start = 1050,
  csv_file = "blog/stringy05_growth_rate/stringy05_weighted_hop_800_1000.csv",
  append = TRUE
)




# Running structured version multiple times.
n <- 500
n_sim <- 100
seed_start <- 1050
noise_sd <- 0.005
t <- seq(-1, 1, length.out = n)
patterns <- list(
  "Degree 1 vs 2" = list(x = t, y = t^2),
  "Degree 2 vs 3" = list(x = t^2, y = t^3),
  "Sine wave" = list(x = t, y = sin(2 * pi * t))
)

results <- bind_rows(lapply(names(patterns), function(name) {
  p <- patterns[[name]]
  bind_rows(lapply(seq_len(n_sim), function(b) {
    seed <- seed_start + b - 1
    set.seed(seed)
    x <- p$x + rnorm(n, sd = noise_sd)
    y <- p$y + rnorm(n, sd = noise_sd)
    data.frame(structure = name, run = b, seed = seed,
               compare_stringy(x, y))
  }))
}))

summary_data <- results |>
  group_by(structure) |>
  summarise(
    weighted_min = min(weighted_stringy),
    weighted_mean = mean(weighted_stringy),
    weighted_max = max(weighted_stringy),
    hop_min = min(hop_stringy),
    hop_mean = mean(hop_stringy),
    hop_max = max(hop_stringy),
    .groups = "drop"
  )

knitr::kable(summary_data, digits = 4,
             caption = "100 simulations per structure: weighted vs hop Stringy")

results |>
  pivot_longer(c(weighted_stringy, hop_stringy),
               names_to = "definition", values_to = "stringy") |>
  ggplot(aes(stringy, colour = definition)) +
  geom_density(linewidth = 0.9) +
  facet_wrap(~ structure, scales = "free_y") +
  labs(x = "Raw Stringy", y = "Density", colour = "Definition") +
  theme_bw()


write.csv(results, "blog/stringy05_growth_rate/stringy_structures_100_results.csv", row.names = FALSE)
write.csv(summary_data, "blog/stringy05_growth_rate/stringy_structures_100_summary.csv", row.names = FALSE)
