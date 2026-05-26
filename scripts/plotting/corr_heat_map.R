cor_mat <- matrix(
  c(
    1,      0.12,     0.508, -0.073,
    NA,      1,      NA,     -0.584,
    NA,   NA,      1,         0.494,
    NA,   NA,     NA,             1
  ),
  nrow = 4,
  byrow = TRUE
)

colnames(cor_mat) <- c("size", "density", "flowers", "type")
rownames(cor_mat) <- c("size", "density", "flowers", "type")

cor_long <- as.data.frame(as.table(cor_mat))

ggplot(cor_long, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(fontface = "bold", aes(label = round(Freq, 2)), size = 8) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "white"
  ) +
  theme_minimal(base_size = 25) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Correlation"
  ) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  coord_fixed()

ggsave("~/Desktop/bomarea_traits/figures/corr_heat_map.pdf", width = 10, height = 10)
