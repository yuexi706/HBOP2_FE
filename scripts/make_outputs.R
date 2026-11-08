# Standalone: Rscript --vanilla scripts/make_outputs.R
required <- c("ggplot2", "patchwork", "ragg", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install output dependencies: ", paste(missing, collapse = ", "))
library(ggplot2)
library(patchwork)
if (!exists("exact_oc")) source("R/design.R")
a <- readRDS("results/analysis.rds")
ex <- a$exact; mc <- a$mc; curves <- a$curves; cal <- a$calibration
theme_set(theme_classic(base_size = 10, base_family = "sans") +
  theme(axis.title = element_text(colour = "#202932"), axis.text = element_text(colour = "#303942"),
        axis.line = element_line(linewidth = .35), axis.ticks = element_line(linewidth = .3),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 9), legend.key.width = grid::unit(1.3, "lines"),
        strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
        panel.spacing = grid::unit(1, "lines"), plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(7, 10, 7, 7)))
cols <- c("BOP2-FE" = "#244B73", "HBOP2-FE" = "#B24E29")
lin <- c("BOP2-FE" = "22", "HBOP2-FE" = "solid")
export_plot <- function(p, name, width, height) {
  path <- file.path("results", "figures", name)
  ggsave(paste0(path, ".pdf"), p, width = width, height = height,
         device = grDevices::cairo_pdf, bg = "white")
  ggsave(paste0(path, ".svg"), p, width = width, height = height,
         device = svglite::svglite, bg = "white")
  ggsave(paste0(path, ".png"), p, width = width, height = height,
         device = ragg::agg_png, dpi = 300, bg = "white")
  ggsave(paste0(path, ".tiff"), p, width = width, height = height,
         device = ragg::agg_tiff, dpi = 600, compression = "lzw", bg = "white")
}
method_scale <- function(p) p + scale_color_manual(values = cols) + scale_linetype_manual(values = lin)
base_curve <- subset(curves, design %in% c("primary_hb", "primary_fixed"))
pc <- ggplot(base_curve, aes(p_t, go, color = method, linetype = method)) +
  geom_hline(yintercept = .10, linewidth = .3, color = "#9AA3AB", linetype = "dotted") +
  geom_vline(xintercept = c(.30, .45), linewidth = .3, color = "#CED3D8", linetype = "dotted") +
  geom_line(linewidth = .85) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = expression("Treatment response probability " * p[T]), y = "Probability of go")
pn <- ggplot(base_curve, aes(p_t, en, color = method, linetype = method)) +
  geom_vline(xintercept = c(.30, .45), linewidth = .3, color = "#CED3D8", linetype = "dotted") +
  geom_line(linewidth = .85) + scale_y_continuous(limits = c(20, 60), breaks = seq(20, 60, 10)) +
  labs(x = expression("Treatment response probability " * p[T]), y = "Expected current-trial N")
f1 <- (method_scale(pc) | method_scale(pn)) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
export_plot(f1 & theme(legend.position = "bottom"), "figure_1_primary_operating_curves", 7.2, 3.55)

sens_ids <- sprintf("sensitivity_%02d", which(a$settings$sensitivity$quality == 1 &
                                                       a$settings$sensitivity$cap != 30))
rob <- subset(curves, design %in% c("primary_hb", "primary_fixed", sens_ids))
rob$label <- ifelse(rob$design == "primary_fixed", "BOP2-FE", paste0("HBOP2-FE: ESS = ", rob$ess))
rob_cols <- c("BOP2-FE" = "#244B73", "HBOP2-FE: ESS = 10" = "#008678",
              "HBOP2-FE: ESS = 20" = "#9C7915", "HBOP2-FE: ESS = 30" = "#B24E29",
              "HBOP2-FE: ESS = 50" = "#805A9A")
rob_null <- subset(rob, p_t >= .20 & p_t <= .45)
rob_alt <- subset(rob, p_t >= .35 & p_t <= .60); rob_alt$p_c_target <- rob_alt$p_t - .15
rp <- function(d, x, ylab) ggplot(d, aes(.data[[x]], go, color = label, linetype = label)) +
  geom_vline(xintercept = .30, color = "#B1B8BF", linewidth = .35, linetype = "dotted") +
  geom_line(linewidth = .7) + scale_color_manual(values = rob_cols) +
  scale_linetype_manual(values = c("BOP2-FE" = "22", "HBOP2-FE: ESS = 10" = "solid",
              "HBOP2-FE: ESS = 20" = "44", "HBOP2-FE: ESS = 30" = "solid", "HBOP2-FE: ESS = 50" = "42")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = expression("Target comparator truth " * p[C]^tar), y = ylab) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
f2 <- (rp(rob_null, "p_t", "Type I error under the target null") +
         geom_hline(yintercept = .1, linewidth = .3, linetype = "dashed", color = "#7D8790") |
       rp(rob_alt, "p_c_target", "Power at target truth + 0.15")) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
export_plot(f2 & theme(legend.position = "bottom"), "figure_2_mismatch_robustness", 7.2, 4.0)

sen <- subset(ex, startsWith(design, "sensitivity_"))
sen <- unique(sen[, c("ess", "scenario", "go", "en")])
sen$scenario <- factor(sen$scenario, levels = paste0("S", 1:6))
scols <- c("S1" = "#244B73", "S2" = "#B24E29", "S3" = "#008678", "S4" = "#805A9A",
           "S5" = "#9C7915", "S6" = "#57656F")
sp <- function(d, y, title) ggplot(d, aes(ess, .data[[y]], color = scenario, shape = scenario,
                                        linetype = scenario)) +
  geom_line(linewidth = .6, show.legend = TRUE) + geom_point(size = 1.8, show.legend = TRUE) +
  scale_color_manual(values = scols, drop = FALSE) +
  scale_shape_manual(values = c(16, 17, 15, 18, 1, 2), drop = FALSE) +
  scale_linetype_manual(values = c("solid", "22", "44", "42", "13", "longdash"), drop = FALSE) +
  scale_x_continuous(breaks = c(5, 15, 30, 50)) + labs(x = expression("Effective ECA size " * ESS[C]), y = title)
f3 <- (sp(subset(sen, scenario %in% c("S1", "S3", "S4")), "go", "Type I error") +
         geom_hline(yintercept = .1, linetype = "dotted", color = "#9AA3AB", linewidth = .3) |
       sp(subset(sen, scenario %in% c("S2", "S5", "S6")), "go", "Power") |
       sp(subset(sen, scenario %in% c("S1", "S2")), "en", "Expected current-trial N")) +
       plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
export_plot(f3 & theme(legend.position = "bottom"), "figure_3_borrowing_sensitivity", 7.2, 3.7)

gd <- rbind(cbind(method = "HBOP2-FE", a$designs$primary_hb$grid),
            cbind(method = "BOP2-FE", a$designs$primary_fixed$grid))
sel <- subset(cal, design %in% c("primary_hb", "primary_fixed"))
f4 <- ggplot(gd, aes(type1, power)) +
  geom_point(aes(color = feasible), size = .45, alpha = .35) +
  geom_vline(xintercept = .10, linetype = "dashed", linewidth = .4, color = "#4F5962") +
  geom_point(data = sel, aes(exact_type1, exact_power), inherit.aes = FALSE,
             shape = 23, size = 2.9, fill = "#E7AF3B", color = "#242D35", stroke = .55) +
  facet_wrap(~method, nrow = 1) + scale_color_manual(values = c("FALSE" = "#A8ADB2", "TRUE" = "#244B73"),
    labels = c("Above alpha ceiling", "Feasible")) +
  scale_x_continuous(limits = c(0, .6), breaks = seq(0, .6, .1), oob = scales::oob_keep) +
  coord_cartesian(xlim = c(0, .6), ylim = c(0, 1)) +
  labs(x = "Exact type I error", y = "Exact compatible power")
export_plot(f4, "figure_4_exact_calibration", 7.2, 3.55)

heat <- do.call(rbind, lapply(c("case_published", "case_fixed"), function(id) {
  z <- a$designs[[id]]; d <- z$posterior
  d$decision <- vapply(seq_len(nrow(d)), function(i) {
    b <- z$boundaries[z$boundaries$n == d$n[i], ]
    if (d$x[i] <= b$futility_max) "Futility / no-go" else if (d$x[i] >= b$efficacy_min) "Efficacy / go" else "Continue"
  }, character(1))
  d$design <- if (id == "case_published") "HBOP2-FE: manuscript = grid optimum" else
    "BOP2-FE: separately calibrated"
  d
}))
hb <- do.call(rbind, lapply(split(heat, heat$design), function(d) {
  id <- if (grepl("HBOP2", d$design[1])) "case_published" else "case_fixed"
  b <- a$designs[[id]]$boundaries; b$design <- d$design[1]
  b$label <- paste0("F <= ", b$futility_max, "\nE >= ", b$efficacy_min); b
}))
heat$n <- factor(heat$n, levels = c(60, 40, 20)); hb$n <- factor(hb$n, levels = c(60, 40, 20))
f5 <- ggplot(heat, aes(x, n, fill = decision)) + geom_tile(height = .65, width = 1) +
  geom_text(data = hb, aes(x = 66, y = n, label = label), inherit.aes = FALSE,
             size = 2.65, lineheight = 1.1, hjust = .5) +
  facet_wrap(~design, ncol = 1) +
  scale_fill_manual(values = c("Futility / no-go" = "#D89E89", "Continue" = "#E5E9EB", "Efficacy / go" = "#497897"),
                    breaks = c("Futility / no-go", "Continue", "Efficacy / go")) +
  scale_x_continuous(breaks = seq(0, 60, 10), limits = c(-.5, 71), expand = c(0, 0)) +
  labs(x = "Cumulative number of treatment responses", y = "Current-trial N") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.spacing = grid::unit(1.3, "lines"))
export_plot(f5, "figure_5_case_decision_regions", 7.2, 4.5)

kc <- subset(curves, design %in% c("case_published", "case_fixed"))
kc$label <- factor(kc$design, levels = c("case_published", "case_fixed"),
                  labels = c("HBOP2-FE: manuscript = optimum", "BOP2-FE: optimized"))
kcols <- c("HBOP2-FE: manuscript = optimum" = "#B24E29", "BOP2-FE: optimized" = "#244B73")
kplot <- function(y, yl) ggplot(kc, aes(p_t, .data[[y]], color = label, linetype = label)) +
  geom_vline(xintercept = c(.27, .47), linewidth = .3, linetype = "dotted", color = "#BDC4CA") +
  geom_line(linewidth = .8) + scale_color_manual(values = kcols) +
  scale_linetype_manual(values = c("solid", "22")) +
  labs(x = expression("Treatment response probability " * p[T]), y = yl)
f6 <- (kplot("go", "Probability of go") + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) |
       kplot("en", "Expected current-trial N") + scale_y_continuous(limits = c(20, 60))) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
export_plot(f6 & theme(legend.position = "bottom"), "figure_6_case_operating_curves", 7.2, 3.7)

vm <- subset(mc, metric %in% c("go", "en"))
vm$metric <- factor(vm$metric, levels = c("go", "en"), labels = c("Go probability", "Expected current-trial N"))
f7 <- ggplot(vm, aes(z)) +
  geom_histogram(binwidth = .4, boundary = 0, fill = "#497897", color = "white", linewidth = .25) +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "#B24E29", linewidth = .45) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = -4:4) + coord_cartesian(xlim = c(-4.5, 4.5)) +
  labs(x = "(Monte Carlo estimate - exact value) / reference SE", y = "Number of scenario-design comparisons")
export_plot(f7, "figure_7_monte_carlo_validation", 7.2, 3.45)

source("scripts/make_report.R")
message("Seven figures exported in PDF, SVG, 300 dpi PNG and 600 dpi TIFF.")
