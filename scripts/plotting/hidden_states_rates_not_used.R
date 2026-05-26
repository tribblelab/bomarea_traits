library(tidyverse)
setwd("~/Desktop/bomarea_traits/output/hisse/unconstrained/")

post <- bind_rows(
  read.table("hisse_rj_run_1.log", header = TRUE, check.names = FALSE),
  read.table("hisse_rj_run_2.log", header = TRUE, check.names = FALSE)
)

hidden_rates <- post %>%
  select(
    `speciation_hidden[1]`,
    `speciation_hidden[2]`,
    `extinction_hidden[1]`,
    `extinction_hidden[2]`
  ) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "multiplier") %>%
  mutate(
    process = ifelse(grepl("speciation", parameter), "Speciation", "Extinction"),
    hidden_state = ifelse(grepl("\\[1\\]", parameter), "Hidden Umbel", "Hidden Compound")
  )

observed_rates <- post %>%
  select(
    `speciation_observed[1]`,
    `speciation_observed[2]`,
    `extinction_observed[1]`,
    `extinction_observed[2]`
  ) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "rate") %>%
  mutate(
    process = ifelse(grepl("speciation", parameter), "Speciation", "Extinction"),
    observed_state = ifelse(grepl("\\[1\\]", parameter), "Observed Umbel", "Observed Compound")
  )

ggplot(hidden_rates, aes(x = hidden_state, y = multiplier, fill = process)) +
  geom_violin(color = "black", linewidth = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ process, scales = "free_x") +
  scale_fill_manual(values = c(
    "Speciation" = "#ff960d",
    "Extinction" = "#00b3ff"
  )) +
  labs(x = "Hidden state", y = "Multiplier") +
  theme_classic(base_size = 14)

ggplot(observed_rates, aes(x = observed_state, y = rate, fill = process)) +
  geom_violin(color = "black", linewidth = 0.1) +
  facet_wrap(~ process, scales = "free_y") +
  scale_fill_manual(values = c(
    "Speciation" = "#ff960d",
    "Extinction" = "#00b3ff"
  )) +
  labs(x = "Observed state", y = "Rate") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")
