library(tidyverse)
library(sbmR)


set.seed(42)
n_blocks <- 3    # Total number of blocks
block_size <- 40 # How many nodes in each block

network <- sim_basic_block_network(
  n_blocks = n_blocks,     
  n_nodes_per_block = block_size,  
  return_edge_propensities = TRUE
)

network$edge_propensities %>% 
  ggplot(aes(x = block_1, y = block_2)) +
  geom_tile(aes(fill = propensity)) +
  theme_minimal()

visualize_network(network, width = '100%')

my_sbm <- create_sbm(network)

collapse_results <- my_sbm %>%
  collapse_blocks(report_all_steps = TRUE, sigma = 1.01)

library(tidyverse)
visualize_collapse_results(collapse_results, heuristic = 'dev_from_rolling_mean') +
  geom_vline(xintercept = 3, color = 'orangered') +
  xlim(0,20)

my_sbm <- choose_best_collapse_state(my_sbm, collapse_results, heuristic = "dev_from_rolling_mean", verbose = TRUE)

merged_state <- my_sbm %>% 
  get_state() %>% 
  select(id, parent)

nodes_w_inferred_block <- network$nodes  %>%
  left_join(merged_state, by = 'id') %>% 
  rename(inferred = parent)

table(nodes_w_inferred_block$inferred,nodes_w_inferred_block$block) %>% knitr::kable()


num_sweeps <- 150

sweep_results <- mcmc_sweep(my_sbm, num_sweeps = num_sweeps, eps = 0.3, track_pairs = TRUE)

sweep_df <- sweep_results$sweep_info %>% 
  mutate(sweep = 1:n(),
         label = 'eps = 0.1') %>% 
  pivot_longer(entropy_delta:num_nodes_moved)

sweep_df %>% 
  ggplot(aes(x = sweep, y = value)) +
  geom_line() +
  facet_grid(name~., scales = "free_y") +
  labs(
    title = glue::glue('Result of {num_sweeps} MCMC sweeps'),
    subtitle = "Entropy Delta of sweep and number of nodes moved for sweep"
  )

sweep_results$pairing_counts %>% 
  ggplot(aes(x = node_a, y = node_b, fill = proportion_connected)) +
  geom_raster() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = "node a", 
       y = "node b",
       fill = "Proportion of sweeps\nin same block",
       title = "Node block consensus matrix",
       subtitle = glue::glue("True number of blocks: {n_blocks}"))

visualize_propensity_network(sweep_results, proportion_threshold = 0.6)
