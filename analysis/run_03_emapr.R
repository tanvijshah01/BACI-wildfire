library(terra); library(sf); library(here); library(dplyr)
library(ggplot2); library(scales); library(glue); library(tigris); library(png)
sf_use_s2(FALSE); options(tigris_use_cache = TRUE)
here::i_am("analysis/03_emapr_biomass_exploration.qmd")

STUDY_YEARS <- c(1990, 1992, 1993); AGG_FACTOR <- 10; SAMPLE_N <- 200000
RAW_DIR <- here("data", "raw", "emapr_biomass")
rast_paths <- setNames(file.path(RAW_DIR, glue("composite_{STUDY_YEARS}_median.tif")), STUDY_YEARS)
for (yr in STUDY_YEARS) if (!file.exists(rast_paths[as.character(yr)])) stop(glue("Missing: {yr}"))
cat("Files found\n")

rast_list <- lapply(STUDY_YEARS, function(yr) terra::rast(rast_paths[as.character(yr)]))
names(rast_list) <- STUDY_YEARS

ca_wgs84    <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |> dplyr::filter(STUSPS == "CA")
stopifnot("CA not found" = nrow(ca_wgs84) == 1)
ca_boundary <- sf::st_transform(ca_wgs84, crs = sf::st_crs(rast_list[[1]]))
ca_bbox     <- sf::st_bbox(ca_boundary)
ca_ext      <- terra::ext(ca_bbox["xmin"], ca_bbox["xmax"], ca_bbox["ymin"], ca_bbox["ymax"])
ca_vect     <- terra::vect(ca_boundary)
cat("CA extent:", paste(round(as.vector(ca_ext), 0), collapse = ", "), "\n")

# ── Summary statistics ─────────────────────────────────────────────────────────
stats_list <- lapply(STUDY_YEARS, function(yr) {
  cat(glue("Stats {yr}..."))
  r       <- terra::crop(rast_list[[as.character(yr)]], ca_ext)
  r       <- terra::mask(r, ca_vect)
  gs      <- terra::global(r, fun = c("min", "max", "mean", "sd"), na.rm = TRUE)
  n_valid <- terra::global(r, fun = "notNA")[[1]]
  vals    <- terra::spatSample(r, size = SAMPLE_N, na.rm = TRUE)[[1]]
  vals    <- vals[vals > 0]
  qs      <- quantile(vals, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
  cat(glue(" mean={round(gs$mean,1)} median={round(qs[3],1)}\n"))
  list(year = yr, mean = gs$mean, sd = gs$sd, n_valid = n_valid,
       n_total = terra::ncell(r), median = qs[3], q95 = qs[5], max = gs$max, vals = vals)
})
names(stats_list) <- STUDY_YEARS

# ── Histogram ─────────────────────────────────────────────────────────────────
hist_df <- do.call(rbind, lapply(stats_list, function(s) {
  data.frame(biomass = s$vals, year = as.factor(s$year))
}))
pHist <- ggplot(hist_df, aes(x = biomass, fill = year)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity", color = NA) +
  scale_fill_manual(values = c("1990" = "#1b7837", "1992" = "#762a83", "1993" = "#e08214")) +
  scale_x_continuous(labels = number_format(accuracy = 1)) +
  scale_y_continuous(labels = comma_format()) +
  labs(title    = "Biomass Distribution — California (1990, 1992, 1993)",
       subtitle = glue("Positive-value pixels only | {scales::comma(SAMPLE_N)} pixel sample per year"),
       x = "Biomass (Mg ha⁻¹)", y = "Pixel count (sample)", fill = "Year",
       caption = "Source: eMapR / LandTrendR, islay.ceoas.oregonstate.edu.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey40", size = 11),
        panel.grid.minor = element_blank())
ggsave(here("figures", "emapr_composite_biomass_histogram.png"),
       plot = pHist, width = 10, height = 6, dpi = 300, bg = "white")
cat("Histogram saved:", dim(png::readPNG(here("figures", "emapr_composite_biomass_histogram.png")))[2], "x",
    dim(png::readPNG(here("figures", "emapr_composite_biomass_histogram.png")))[1], "px\n")

# ── Plot A ────────────────────────────────────────────────────────────────────
trend_df <- do.call(rbind, lapply(stats_list, function(s) {
  se <- s$sd / sqrt(s$n_valid)
  data.frame(year = s$year, mean = s$mean, ymin = s$mean - se, ymax = s$mean + se)
}))
pA <- ggplot(trend_df, aes(x = factor(year), y = mean)) +
  geom_col(fill = "#2c7bb6", alpha = 0.85, width = 0.5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.12, linewidth = 0.8, color = "grey20") +
  geom_text(aes(label = round(mean, 1)), vjust = -0.8, size = 3.5, color = "grey20") +
  scale_y_continuous(labels = number_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Mean Aboveground Biomass — California (1990, 1992, 1993)",
       subtitle = "eMapR CONUS composites | Error bars = ± 1 SE across all valid California pixels",
       x = "Year", y = "Mean AGB (Mg ha⁻¹)",
       caption = "Source: eMapR / LandTrendR, islay.ceoas.oregonstate.edu.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey40", size = 11),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 11))
ggsave(here("figures", "emapr_composite_biomass_trend.png"),
       plot = pA, width = 10, height = 6, dpi = 300, bg = "white")
cat("Plot A saved\n")

# ── Maps ──────────────────────────────────────────────────────────────────────
all_forested <- unlist(lapply(stats_list, "[[", "vals"))
agb_upper    <- quantile(all_forested, 0.99, na.rm = TRUE)
cat("Color scale cap:", round(agb_upper, 1), "Mg/ha\n")

for (yr in STUDY_YEARS) {
  cat(glue("Map {yr}..."))
  lyr     <- terra::crop(rast_list[[as.character(yr)]], ca_ext)
  lyr     <- terra::aggregate(lyr, fact = AGG_FACTOR, fun = "mean", na.rm = TRUE)
  agg_res <- round(terra::res(lyr)[1])
  ca_mask <- terra::rasterize(ca_vect, lyr, field = 1, background = NA)
  lyr     <- terra::mask(lyr, ca_mask)
  lyr     <- terra::clamp(lyr, lower = 0, upper = agb_upper, values = TRUE)

  fpath <- here("figures", glue("emapr_composite_biomass_map_{yr}.png"))
  png(fpath, width = 8, height = 9, units = "in", res = 300, bg = "white")
  terra::plot(lyr,
    col  = hcl.colors(100, "viridis"), range = c(0, agb_upper),
    main = glue("Aboveground Biomass — California ({yr})"),
    plg  = list(title = glue("AGB (Mg/ha)\n(capped at {round(agb_upper)} — 99th pct)"),
                title.cex = 0.75, cex = 0.75),
    mar  = c(3, 3, 4, 7), cex.main = 1.1)
  terra::plot(ca_vect, add = TRUE, col = NA, border = "grey20", lwd = 0.6)
  mtext(glue("Source: eMapR / LandTrendR, islay.ceoas.oregonstate.edu. Aggregated to {agg_res} m."),
        side = 1, line = 1.8, cex = 0.65, col = "grey50")
  dev.off()
  dims <- dim(png::readPNG(fpath))
  cat(glue(" {dims[2]}x{dims[1]}px saved\n"))
}

cat("\n=== All figures generated ===\n")
