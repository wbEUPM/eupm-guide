#' PL: decode polygon IDs by their types
#' @export
pl_decode_id <- function(dta) {
  dta |> 
    mutate(id1_NUTS1 = str_sub(id, 1, 2),
           id2_NUTS2 = str_sub(id, 3, 5),
           id3_NUTS3 = str_sub(id, 6, 7),
           id4_LAU1 = str_sub(id, 8, 9),
           id5_LAU2 = str_sub(id, 10, 11),
           type_LAU = str_sub(id, 12)
    ) |> 
    pl_rename_levels_dta()
}

pl_rename_levels_dta <- function(dta) {
  dta |>
    mutate(
      level = case_when(
        id1_NUTS1 == "00" ~ "level0",
        id2_NUTS2 == "00" ~ "level1",
        id3_NUTS3 == "00" ~ "level3",
        id4_LAU1 == "00" ~ "level4",
        id5_LAU2 == "00" ~ "level5",
        .default = "level6",
      )
    )
}

#' PL: Collect data from poland API
#' @export
pl_collect_api <- function(var_ids, bd_raw, pin_prefix = "", ...) {
  var_ids |>
    walk( ~ {
      empl <-
        get_data_by_variable(
          varId = .x,
          unitParentId = NULL,
          unitLevel = NULL,
          aggregateId = NULL,
          year = 1990:2030,
          lang = "en"
        )
      bd_raw |> pin_write(empl,
                          name = str_c(pin_prefix, .x),
                          type = "parquet")
    }, .progress = T)
}

#' PL: reload and clean API variables
#' @export
pl_load_recode_raw <-
  function(var_ids, bd_raw, pin_prefix, meta_vars = NULL, ...) {
    dta <- 
      var_ids |>
      map(~ {
        bd_raw |>
          pin_read(name = str_c(pin_prefix, "_", .x)) |>
          mutate(var = .x)
      }) |>
      bind_rows() |>
      dplyr::select(id, var, year, val) |> pl_decode_id() |>
      dplyr::filter(type_LAU %in% c("0", "1", "2", "3")) |>
      dplyr::select(id, year, var, val)
    
    if (!is.null(meta_vars)) {
      var_match <- 
        meta_vars |>
        dplyr::filter(var_id_api %in% var_ids) |> 
        dplyr::select(var_id_api, var = var_id) |> 
        mutate(var = ifelse(is.na(var), var_id_api, var))
      dta <- dta |> 
        rename(var_id_api = var) |> 
        left_join(var_match, by = join_by(var_id_api)) |> 
        dplyr::select(-var_id_api)
    }
    
    dta |>
      dplyr::select(id, year, var, val)
  }


#' Rename admin levels 
#' @export
pl_rename_levels <- function(dta) {
  dta |> 
    mutate(
      level = case_when(
        level == "level0" ~ "Country",
        level == "level1" ~ "NUTS 1",
        level == "level2" ~ "Regions",
        level == "level3" ~ "NUTS 2",
        level == "level4" ~ "NUTS 3",
        level == "level5" ~ "LAU 1",
        level == "level6" ~ "LAU 2"
      )
    )
  
}


# Generic functions --------------------------------------------------------

#' Return name of a geometries level from the pre-defined geometery list
#' 
get_layer_name <- function(nm) {
  nm <- nm |> str_split("_|\\.", simplify = T) |> as.character()
  if (length(nm) == 3) {
    nnm_out <- str_c(nm[[3]], " (", nm[[2]], ")")
  } else if (length(nm) == 2) {
    nnm_out <- nm[[2]]
  } else {
    nnm_out <- nm
  }
  nnm_out
  
}


#' Plot admin boundaries
#' @export
plot_admin_map <- function(dta) {
  dta <- dta |>
    imap( ~ {
      out <- list(dta = .x, name = get_layer_name(.y))
      out$name_legend <- str_c(out$name, " (", nrow(.x), ")")
      out
    }) |>
    unname(force = F)
  
  list(ggplot()) |>
    append(rev(dta)) |>
    reduce( ~ {
      .x +
        geom_sf(aes(colour = .y$name_legend, linewidth = .y$name_legend),
                data = .y$dta,
                fill = NA)
    }) +
    scale_linewidth_manual(
      NULL, 
      breaks = dta |> map_chr("name_legend"),
      values = rev(seq(from = 0.25, by = 0.25, length.out = length(dta)))
      ) +
    scale_color_brewer(NULL, 
                       breaks = dta |> map_chr("name_legend"),
                       palette = "Set1") +
    theme_bw() +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.01, 0.01),
      legend.justification = c("left", "bottom"),
      legend.box.just = "center",
      legend.margin = margin(0, 0, 0, 0)
    )
}

#' Get N of missing data points
#' @export
get_missing <- function(dta, geoms, range = c(2011, 2019:2030), focus_levels = 4:5) {
  geoms |>
    imap( ~ .x |> mutate(level = .y, nn = n())) |>
    bind_rows() |>
    st_drop_geometry() |>
    left_join(dta, by = join_by(id)) |>
    as_tibble() |>
    group_by(level) |> 
    mutate(nn = max(nn, na.rm = TRUE)) |> 
    group_by(level, year, var, nn) |>
    summarise(n = n(),
              missing = str_c(sum(is.na(val) | val == 0)),# "/", n()),
              .groups = "drop") |>
    dplyr::filter(year %in% range)  |>
    dplyr::select(-n) |>
    pivot_wider(values_from = missing, names_from = year) |> 
    mutate(across(matches("\\d{4}"), ~ ifelse(is.na(.), as.character(nn), .))) |> 
    dplyr::filter(level %in% str_c("level", focus_levels)) |> 
    pl_rename_levels()
}

add_varnames <- function(dta, meta) {
  dta |>
    left_join(meta |>
                mutate(var_name = str_c(var_name, " (", var_units, ")")) |>
                dplyr::select(var = var_id, var_name),
              by = join_by(var))
  
}

get_format_missing <- function(dta, geoms, meta, range = c(2019:2030)) {
  in_range <- c(min(dta$year), range) |> unique()
  dta |>
    get_missing(geoms, in_range) |>
    add_varnames(meta) |>
    dplyr::select(var_name, level, matches("\\d{3,}")) |>
    arrange(var_name) |>
    group_by(var_name) |>
    as_flextable(hide_grouplabel = TRUE)
}


#' Get variables-specific coverage of by years
#' @export
get_coverage <- function(dta, keep_all = FALSE, years = NULL) {
  
  y_rel <- c(2019:2100, years) |> unique() |> sort()
  
  dta <- dta |> 
    group_by(var, year) |>
    summarise(n = n()) |> 
    ungroup() |>
    pivot_wider(names_from = var, values_from = n) |> 
    arrange(year)
  
  if (!keep_all) dta <- dta |> dplyr::filter(year %in% y_rel)
  
  dta
}


# Plotting functions --------------------------------------------------------

get_map <- function(dta, legend_name = NULL, title, ...) {
  if (nrow(dta) == 0) {
    out_plt <- dta |> ggplot() + theme_bw() +
      ggtitle(str_c(title, ". NO DATA"))
    return(out_plt)
  }
  dta |> 
    ggplot() + 
    geom_sf() + 
    aes(fill = val) +
    scale_fill_viridis_c(
      trans = "log", 
      option = "B", 
      direction = -1,
      # breaks = c(1, 5, 10, 20, 50, 100),
      n.breaks = 7,
      label = scales::number_format(0.1),
      name = legend_name ,
      guide = guide_legend()
    ) + 
    theme_bw() + 
    ggtitle(title)
}

plot_var <- 
  function(dta, var_choice, geoms, meta, focus = 4:5, years = 2023) {
    var_nm <- tibble(var = var_choice) |> add_varnames(meta) |> pull(var_name)
    if (is.na(var_nm)) var_nm <- var_choice
    dta <- dta |> dplyr::filter(var %in% var_choice)
    
    geoms |>
      imap(~ tibble(level = .y, data = list(.x))) |>
      bind_rows() |>
      dplyr::filter(level  %in% str_c("level", focus)) |> 
      pl_rename_levels() |>
      mutate(data = map2(level, data, ~ {
        .y |>
          left_join(dta, by = join_by(id)) |>
          dplyr::filter(year == years) |>
          get_map(legend_name = years, title = .x) #+
          # theme(plot.margin = unit(c(0.0, 0.0, 0.0, 0.0), "inches")) 
      })) |>
      pull(data) |> 
      reduce(~ .x + .y) + 
      plot_annotation(title = var_nm)
  }


plot_var_by_year <- 
  function(dta, var_choice, geoms, meta, focus = 4:5, years = c(2011, 2018, 2023)) {
    var_nm <- tibble(var = var_choice) |> add_varnames(meta) |> 
      mutate(var_name = str_c(var_name, " '", var, "'")) |> 
      pull(var_name)
    if (is.na(var_nm)) var_nm <- var_choice
    dta <- dta |> dplyr::filter(var %in% var_choice)
    
    all_plts <- 
      geoms |>
      imap(~ tibble(level = .y, data = list(.x))) |>
      bind_rows() |>
      dplyr::filter(level  %in% str_c("level", focus)) |> 
      pl_rename_levels() |>
      mutate(data = map2(level, data, ~ {
        adm_level_name <- .x
        adm_level <- .y
        all_dta <- adm_level |> left_join(dta, by = join_by(id)) |> dplyr::filter(!is.na(val))
        all_years <- all_dta |> pull(year) |> unique()
        year_min <- ifelse(min(years) %in% all_years, min(years), min(all_years))
        year_max <- ifelse(max(years) %in% all_years, max(years), max(all_years))
        year_mid <- ifelse(any(years[!years %in% c(min(years), max(years))] %in% all_years), 
                           years[!years %in% c(min(years), max(years))], 
                           all_years[!years %in% c(year_min, year_max)])
        years_focus <- c(year_min, year_mid, year_max) |> unique() |> sort()
        years_focus |>
          map(~{
            all_dta |> dplyr::filter(year == .x) |>
              get_map(
                legend_name = NULL,
                title = str_c(.x, ": ", adm_level_name))
          })
      }))  |>
      pull(data) |> 
      unlist(recursive = F) 
    
    if (length(all_plts) == 1) {
      return(all_plts[[1]])
    } else {
      all_plts |> reduce(~ .x + .y) + plot_annotation(title = var_nm)
    }
  }


# get_codes <- function(dta) {
#   dta |> 
#     mutate(ID2_NUTS2 = str_sub(Code, 1, 2),
#            ID4_LAU1 = str_sub(Code, 3, 4),
#            ID5_LAU2 = str_sub(Code, 5, 6),
#            LAU2_type = str_sub(Code, 7, 7),
#     ) |>
#     
#     # Removing regional aggregates
#     filter(ID2_NUTS2 != "00", ID4_LAU1 != "00", ID5_LAU2 != "00", LAU2_type != "0") |> 
#     
#     # Filtering only 1 - urban, 2 - rural, 3 - peri-urban, and 8 - capital gminas.
#     filter(LAU2_type %in% c(1, 2, 3, 8)) |> 
#     select(matches("ID\\d"), everything()) |> 
#     select(-Code)
# }

get_plot_dta <- function(dta, geoms = NULL, var_ids = NULL, years = NULL, ...) {
  dta <- 
    dta |> 
    dplyr::filter(var_id %in% var_ids, year %in% years, ...) #|> 
  # pivot_wider(names_from = var_id, values_from = value)
  
  if (!is.null(geoms)) {
    reg_all <- names(geoms)
    reg_all_dta <- dta |> dplyr::select(matches("ID\\d")) |> mutate(ID6_LAU3 = NA) |> 
      dplyr::select(-where(~ all(is.na(.)))) |> names()
    reg_max <- reg_all_dta |> max() |> str_remove("ID\\d_")
    
    return(left_join(geoms[[reg_max]], dta))
    
    
  } else {
    return(dta)
  }
}


get_stats <- function(dta, keep_all = FALSE, years = NULL) {
  
  y_rel <- c(min(dta$year), 2011, 2019, max(dta$year), years) |> 
    unique() |> sort()
  
  dta <- dta  |> 
    group_by(var_id, year) |> 
    summarise(
      n = n(), 
      mean = mean(val, na.rm = TRUE),
      sd = sd(val, na.rm = TRUE),
      min = min(val, na.rm = TRUE),
      max = max(val, na.rm = TRUE), 
      .groups = "drop"
    )
  
  if (!keep_all) dta <- dta |> dplyr::filter(year %in% y_rel)
  
  dta
}

format_stats <- function(dta) {
  dta |> 
    mutate(
      across(c(year, n), ~ number(., accuracy = 1, big.mark = "")),
      across(c(mean, sd, min, max), ~ number(., accuracy = .1, big.mark = ""))
    ) |>
    group_by(var_id) |> 
    as_flextable(hide_grouplabel = TRUE)
}


FitFlextableToPage <- function(ft, pgwidth = 6.5, size = 9) {
  ft_out <-
    ft |> fontsize(size = size, part = "all") |> 
    font(
      fontname = "Times New Roman",
      part = "all") |> 
    autofit()
  
  ft_out <- width(ft_out,
                  width = dim(ft_out)$widths * pgwidth / (flextable_dim(ft_out)$widths))
  
  ft_out <- 
    ft_out |> 
    bold(part = "header") |> 
    valign(valign  = "center", part = "all") |> 
    align(j = -1, align = "center", part = "header") |> 
    align(j = -1, part = "body")
  
  return(ft_out)
}

