
fct_map_chreplot <- function(dta, 
                             var_fill,
                             title_glue = "",
                             subtitle_glue = "",
                             legend_glue = "",
                             fltr_year = NULL,
                             force_breaks = NULL, 
                             label_form = label_number(0.01),
                             probs = seq(from = 0, to = 1, by = 0.16),
                             ...) {
  
  if (is.null(dta)) return(plot_spacer()) #guide_area()
  
  dta <- dta |> mutate(across(any_of(var_fill), ~ ., .names = "fill_var"))
  
  if (all(is.null(force_breaks))) {
    force_breaks <- quantile(dta$fill_var, na.rm = T, probs = probs)
  } 
  
  force_breaks <-  c(
    min(force_breaks, dta$fill_var, na.rm = T),
    force_breaks,
    max(force_breaks, dta$fill_var, na.rm = T)
  ) |>
    unique() |>
    sort()
  
  labels <- str_c(
    force_breaks[-length(force_breaks)] |> label_form(),
    "-",
    force_breaks[-1]|> label_form()
  ) |> 
    unique()
  
  if (!is.null(fltr_year)) {
    dta <- dta |> filter(year %in% fltr_year)
  }
  
  dta <-
    dta |>
    mutate(
      fill_var_cat = cut(
        fill_var,
        include.lowest = T,
        breaks = force_breaks,
        labels = labels
      ) |> fct_relevel(labels),
      title = glue(title_glue),
      subtitle = glue(subtitle_glue),
      legendtitle = glue(legend_glue)
    )
  
  dta_titles <- 
    dta |>
    st_drop_geometry() |> 
    summarise(across(contains("title"), ~ first(.))) |> 
    as.list() |> 
    keep(~.x != "" && !is.na(.x))
  
  dta |> 
    ggplot() +
    geom_sf(
      mapping = aes(fill = fill_var_cat),
      color = "white",
      size = 0.1,
      show.legend = TRUE
    ) +
    scale_fill_viridis_d(
      dta_titles$legendtitle,
      direction = -1, drop = F) +
    labs(title = dta_titles$title, subtitle = dta_titles$subtitle) +
    geom_sf(data = geom_nuts2, fill = NA, colour = "darkred") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )
}

fct_plot_scatter <- function(dta,
                             x_var,
                             y_var,
                             colour_var = NULL,
                             group_var = colour_var,
                             x_title_glue = "",
                             y_title_glue = "",
                             title_glue = "",
                             subtitile_glue = "",
                             legendtitle_glue = "",
                             scale_x_label = label_number(0.001),
                             scale_y_label = scale_x_label,
                             add_point = T,
                             add_smooth = TRUE,
                             add_line = TRUE,
                             expand_x = 0,
                             expand_y = 0,
                             
                             ...
                             
) {
  if (is.null(dta)) return(plot_spacer())
  
  dta <- 
    dta |>
    mutate(
      across(all_of(x_var), ~ ., .names = "x_var"),
      across(all_of(y_var), ~ ., .names = "y_var"),
      x_title = glue(x_title_glue),
      y_title = glue(y_title_glue),
      title = glue(title_glue),
      subtitle = glue(subtitile_glue),
      legendtitle = glue(legendtitle_glue)
    )
  
  if (!is.null(colour_var)) {
    dta <- dta |> mutate(across(all_of(colour_var), ~ ., .names = "colour_var"))
  }
  
  if (!is.null(group_var)) {
    dta <- dta |> mutate(across(all_of(group_var), ~ ., .names = "group_var"))
  }
  
  dta_titles <-
    dta |>
    st_drop_geometry() |>
    summarise(across(contains("title"), ~ first(.))) |>
    as.list() |>
    keep(~ .x != "" && !is.na(.x))
  
  plt <- dta |> ggplot()  + aes(x = x_var, y = y_var)
  
  if (!is.null(colour_var)) {
    plt <- dta |> ggplot()  +
      aes(x = x_var, y = y_var, colour = colour_var)
  }
  
  if (!is.null(group_var)) {
    plt <- dta |> ggplot() + aes(x = x_var, y = y_var, group = group_var)
  }
  
  if (!is.null(group_var) && !is.null(colour_var)) {
    plt <- dta |> ggplot() + aes(x = x_var, y = y_var, group = group_var, colour = colour_var)
  }
  
  plt <- 
    plt +
    scale_x_continuous(label = scale_x_label) +
    scale_y_continuous(label = scale_y_label) +
    expand_limits(x = expand_x, y = expand_y) +
    xlab(dta_titles$x_title) +
    ylab(dta_titles$y_title) +
    labs(title = dta_titles$title,
         subtitle = dta_titles$subtitle) +
    scale_color_brewer(dta_titles$legendtitle, palette = "Set1") +
    theme_bw()
  
  if (add_point) plt <- plt + geom_point()
  if (add_smooth) plt <- plt + geom_smooth(method = 'loess', formula = 'y ~ x')
  if (add_line) plt <- plt + geom_abline(slope = 1, intercept = 0, linetype = 3)
  
  plt
}

