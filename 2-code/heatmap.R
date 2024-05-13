df = 
  tibble(
    C = runif(120),
    N = runif(120),
    month = factor(month.name) |> rep(10)
  ) |>
  pivot_longer(cols = 1:2)


ggplot(data = df, 
       aes(x = month, y = name, fill = value)) + 
  geom_raster()
