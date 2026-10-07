#### Chapter 2 analysis ####
# Soil legacy and shade effects on wilding conifer seedling biomass and
# ectomycorrhization (Pinus contorta, Pinus radiata, Pseudotsuga menziesii).
#
# Input : potdat_sup.rds (built from the raw datasheets)
# Output: fitted models, tables, Figures 2.4 and 2.5
#
# Models are fitted with brms; fits are saved via `file =` and reloaded on
# later runs, so the full script only needs to sample once.

remove(list = ls())
library(brms)      # Bürkner 2017, 2018, 2021; loo comes with brms

# shared sampler settings (4 chains x 7,500 post-warmup draws = 30,000)
SEED  <- 67
ITER  <- 10000
WARM  <- 2500
CTRL  <- list(adapt_delta = 0.999, max_treedepth = 15)

fmt <- function(x, d = 2)
  sprintf("%.*f [%.*f, %.*f]", d, median(x), d,
          quantile(x, 0.025), d, quantile(x, 0.975))

species_keys  <- c("PICO", "PIRA", "PSME")
species_names <- c(PICO = "P. contorta", PIRA = "P. radiata", PSME = "Ps. menziesii")


#### 1. Data ####
dat <- readRDS("potdat_sup.rds")
dat <- dat[dat$analysed, ]                       # 506 seedlings
dat$bulkP <- as.numeric(dat$bulk == "P")         # 1 = pine legacy bulk soil

#### 2. Biomass: log-scale Ricker function (Equation 2.2) ####
#this is the main non-linear mixed effects model
# log(mass) = log a_k + log(PPFD transmitted) - b_k * PPFD transmitted + u
# k = species x inoculum; u = random intercepts for region and region x tent
# PPFD transmitted is on a 0-100 scale.

formula_mass <- bf(
  final_mass ~ loga + log(ppfd_trans) - exp(logb) * ppfd_trans + re,
  loga ~ 0 + species:inoculate,
  logb ~ 0 + species:inoculate,
  re   ~ 0 + (1 | region/tent),
  nl   = TRUE
)

# priors: peak near 27% transmission (log b ~ -3.3);
# peak mass ~0.5-2 g implies log a ~ -3 to -1.5
priors_mass <- c(
  prior(normal(-2, 1),     class = "b",  nlpar = "loga"),
  prior(normal(-3.3, 0.5), class = "b",  nlpar = "logb"),
  prior(exponential(1),    class = "sd", nlpar = "re")
)

mod_mass <- brm(formula_mass, data = dat, prior = priors_mass,
                family = lognormal(), chains = 4, cores = 4,
                iter = ITER, warmup = WARM, seed = SEED, control = CTRL,
                file = "ch2_models/mod_mass")
summary(mod_mass)


#### 3. Bulk soil check ####
# adds a single effect of pine legacy bulk (sterilised) soil on amplitude
formula_bulk <- bf(
  final_mass ~ loga + log(ppfd_trans) - exp(logb) * ppfd_trans + re,
  loga ~ 0 + species:inoculate + bulkP,
  logb ~ 0 + species:inoculate,
  re   ~ 0 + (1 | region/tent),
  nl   = TRUE
)
priors_bulk <- c(priors_mass,
                 prior(normal(0, 1), class = "b", nlpar = "loga", coef = "bulkP"))

mod_bulk <- brm(formula_bulk, data = dat, prior = priors_bulk,
                family = lognormal(), chains = 4, cores = 4,
                iter = ITER, warmup = WARM, seed = SEED, control = list(adapt_delta = 0.9999, max_treedepth = 15),
                file = "ch2_models/mod_bulk")
summary(mod_bulk)

loo_mass <- loo(mod_mass, reloo = TRUE)
loo_bulk <- loo(mod_bulk, reloo = TRUE)

bulk_cmp <- loo_compare(loo_mass, loo_bulk)
bulk_cmp
saveRDS(bulk_cmp, "bulk_soil_loo.rds")


#### 4. Hypothesis tests (Equations 2.3-2.6; Tables 2.2-2.3) ####
draws   <- as.data.frame(as_draws_df(mod_mass))
mean_re <- rowMeans(draws[, grep("^r_region__re\\[", names(draws))])

min_trans <- min(dat$ppfd_trans)
max_trans <- max(dat$ppfd_trans)

# derived quantities for one species x inoculum group
ricker_group <- function(sp, inoc) {
  loga <- draws[[paste0("b_loga_species", sp, ":inoculate", inoc)]]
  b    <- exp(draws[[paste0("b_logb_species", sp, ":inoculate", inoc)]])
  peak_trans <- 1 / b                                         # Eq 2.2
  peak_mass  <- exp(loga + log(peak_trans) - 1 + mean_re)     # Eq 2.3
  mass_min   <- exp(loga + log(min_trans) - b * min_trans + mean_re)
  list(loga = loga, b = b,
       peak_shade = 100 - peak_trans,
       peak_mass  = peak_mass,
       slope      = (mass_min - peak_mass) / (peak_trans - min_trans))  # Eq 2.4, g per % reduced
}

hyp <- do.call(rbind, lapply(species_keys, function(sp) {
  G <- ricker_group(sp, "G")
  P <- ricker_group(sp, "P")
  cross      <- (G$loga - P$loga) / (G$b - P$b)                # Eq 2.5
  does_cross <- is.finite(cross) & cross > 0 & cross < max_trans
  data.frame(
    species      = sp,
    peak_shade_G = fmt(G$peak_shade, 1), peak_shade_P = fmt(P$peak_shade, 1),
    shade_diff   = fmt(P$peak_shade - G$peak_shade, 1),
    peak_mass_G  = fmt(G$peak_mass),     peak_mass_P  = fmt(P$peak_mass),
    mass_diff    = fmt(P$peak_mass - G$peak_mass),
    slope_G      = fmt(G$slope, 4),      slope_P      = fmt(P$slope, 4),
    slope_diff   = fmt(P$slope - G$slope, 4),
    pr_mass      = round(mean(P$peak_mass > G$peak_mass), 3),        # H1
    pr_steeper   = round(mean(abs(P$slope) > abs(G$slope)), 3),      # H2
    pr_cross     = round(mean(does_cross), 3)                        # H3
  )
}))

write.csv(hyp, "table_2_2_2_3_hypotheses.csv", row.names = FALSE)


#### 5. Ectomycorrhization ####

## 5a. Table 2.4: raw ECM colonization by species and inoculum
tab24 <- aggregate(ecm_colonization ~ species + inoculate, data = dat,
                   FUN = function(x) c(mean = mean(x), sd = sd(x),
                                       median = median(x), n = length(x)))
tab24 <- data.frame(tab24[, 1:2], round(tab24$ecm_colonization, 1))
tab24
write.csv(tab24, "table_2_4_ecm_raw.csv", row.names = FALSE)

# P. contorta in grassland soil: bimodal colonization
pg <- dat$ecm_colonization[dat$species == "PICO" & dat$inoculate == "G"]
round(100 * c(below_20 = mean(pg < 20), above_80 = mean(pg > 80)), 1)

## 5b. Zero-inflated beta model
# colonization (proportion) ~ light x inoculum x species;
# probability of zero colonization ~ species x inoculum
priors_ecm <- c(
  prior(normal(0, 1),   class = "b"),
  prior(normal(0, 1),   class = "b", dpar = "zi"),
  prior(exponential(1), class = "sd")
)

mod_ecm <- brm(
  bf(ecm_prop ~ ppfd_prop * inoculate * species + (1 | region/tent),
     zi       ~ species * inoculate),
  data = dat, prior = priors_ecm, family = zero_inflated_beta(),
  chains = 4, cores = 4, iter = ITER, warmup = WARM, seed = SEED,
  control = list(adapt_delta = 0.999, max_treedepth = 17),
  file = "mod_ecm"
)
summary(mod_ecm)

## 5c. Predicted colonization at mean light, P - G difference, Pr(P > G)
mean_prop <- mean(dat$ppfd_prop)

ecm_pred <- do.call(rbind, lapply(species_keys, function(sp) {
  nd   <- data.frame(ppfd_prop = mean_prop, inoculate = c("G", "P"), species = sp)
  ep   <- posterior_epred(mod_ecm, newdata = nd, re_formula = NA) * 100
  diff <- ep[, 2] - ep[, 1]
  data.frame(species = sp,
             G = fmt(ep[, 1], 1), P = fmt(ep[, 2], 1),
             P_minus_G = fmt(diff, 1),
             pr_P_gt_G = round(mean(diff > 0), 3))
}))

write.csv(ecm_pred, "ecm_predicted.csv", row.names = FALSE)

## 5d. Effect of light on colonization (logit scale, per % PPFD transmitted)
de <- as.data.frame(as_draws_df(mod_ecm))
coef_or_0 <- function(nm) if (nm %in% names(de)) de[[nm]] else 0

ecm_light <- do.call(rbind, lapply(species_keys, function(sp) {
  slope_G <- de$b_ppfd_prop + coef_or_0(paste0("b_ppfd_prop:species", sp))
  slope_P <- slope_G + de$`b_ppfd_prop:inoculateP` +
    coef_or_0(paste0("b_ppfd_prop:inoculateP:species", sp))
  data.frame(species = sp,
             slope_G = fmt(slope_G / 100, 4),
             slope_P = fmt(slope_P / 100, 4),
             P_minus_G = fmt((slope_P - slope_G) / 100, 4))
}))

ecm_light

## 5c. Effect of light on colonization for all soil treatments)
ecm_light_avg <- do.call(rbind, lapply(species_keys, function(sp) {
  slope_G <- de$b_ppfd_prop + coef_or_0(paste0("b_ppfd_prop:species", sp))
  slope_P <- slope_G + de$`b_ppfd_prop:inoculateP` +
    coef_or_0(paste0("b_ppfd_prop:inoculateP:species", sp))
  data.frame(species = sp,
             slope_avg = fmt((slope_G + slope_P) / 2 / 100, 3),   # averaged over soils
             P_minus_G = fmt((slope_P - slope_G) / 100, 3))       # soil difference
}))
ecm_light_avg

write.csv(ecm_light, "ecm_light_slopes.csv", row.names = FALSE)


#### 6. Ectomycorrhization and biomass (Table 2.5) ####
# ECM was measured at harvest, so these models test association,
# not a causal effect of ECM on growth.

formula_ecm_only <- bf(
  final_mass ~ loga + log(ppfd_trans) - exp(logb) * ppfd_trans + re,
  loga ~ 0 + species + ecm_prop,
  logb ~ 0 + species,
  re   ~ 0 + (1 | region/tent),
  nl   = TRUE
)
formula_soil_ecm <- bf(
  final_mass ~ loga + log(ppfd_trans) - exp(logb) * ppfd_trans + re,
  loga ~ 0 + species:inoculate + ecm_prop,
  logb ~ 0 + species:inoculate,
  re   ~ 0 + (1 | region/tent),
  nl   = TRUE
)
priors_ecm_mass <- c(priors_mass,
                     prior(normal(0, 1), class = "b", nlpar = "loga", coef = "ecm_prop"))

mod_ecm_only <- brm(formula_ecm_only, data = dat, prior = priors_ecm_mass,
                    family = lognormal(), chains = 4, cores = 4,
                    iter = ITER, warmup = WARM, seed = SEED, control = list(adapt_delta = 0.9999, max_treedepth = 15),
                    file = "ch2_models/mod_ecm_only")
mod_soil_ecm <- brm(formula_soil_ecm, data = dat, prior = priors_ecm_mass,
                    family = lognormal(), chains = 4, cores = 4,
                    iter = ITER, warmup = WARM, seed = SEED, control = list(adapt_delta = 0.9999, max_treedepth = 15),
                    file = "ch2_models/mod_soil_ecm")

mod_soil_ecm <- readRDS("mod_soil_ecm.rds")
summary(mod_soil_ecm)          # check Rhat = 1.00


summary(mod_ecm_only)
summary(mod_soil_ecm)

loo_ecm_only <- loo(mod_ecm_only, reloo = TRUE)
loo_soil_ecm <- loo(mod_soil_ecm, reloo = TRUE)
tab25 <- loo_compare(loo_ecm_only, loo_soil_ecm, loo_mass)   # loo_mass from section 3
tab25
saveRDS(tab25, "table_2_5_loo.rds")

#mediation check (uses mod_soil_ecm)
## 6b. How much of the soil effect does ECM account for?
# ratio of pine legacy to grassland peak mass, without ECM (total)
# and with ECM held constant (direct)
peak_log_ratio <- function(mod, sp) {
  dr <- as.data.frame(as_draws_df(mod))
  la <- function(i) dr[[paste0("b_loga_species", sp, ":inoculate", i)]]
  lb <- function(i) dr[[paste0("b_logb_species", sp, ":inoculate", i)]]
  (la("P") - la("G")) - (lb("P") - lb("G"))
}

mediation <- do.call(rbind, lapply(species_keys, function(sp) {
  tot <- peak_log_ratio(mod_mass,     sp)
  dir <- peak_log_ratio(mod_soil_ecm, sp)
  data.frame(species       = sp,
             total_ratio   = fmt(exp(tot)),
             direct_ratio  = fmt(exp(dir)),
             share_via_ecm = round(1 - median(dir) / median(tot), 2),
             pr_direct_gt1 = round(mean(dir > 0), 3))
}))

write.csv(mediation, "ecm_mediation.csv", row.names = FALSE) #this became table 2.6


#### 7. Figure 2.4: biomass across the shade gradient ####
# curves: posterior median biomass (lognormal median, matching Eq 2.3)
#         at the average region; ribbons: 95% credible intervals
# points: seedlings, sized by ECM colonization
observed_only <- FALSE    # TRUE: draw curves only across the observed light range

pdf("figure_2_4_biomass.pdf", width = 5.1, height = 6)
par(mfrow = c(3, 1),
    omi = c(0.4, 0.2, 0.2, 0.9),
    mai = c(0.4, 0.4, 0.2, 0.2),
    xpd = NA)

for (sp in species_keys) {
  d_sp <- dat[dat$species == sp, ]
  
  trans_seq <- if (observed_only) seq(min_trans, max_trans, length.out = 200) else
    seq(1, 100, length.out = 200)
  shade_seq <- 100 - trans_seq
  
  pred <- lapply(c(G = "G", P = "P"), function(inoc) {
    g <- ricker_group(sp, inoc)
    m <- sapply(trans_seq, function(x) exp(g$loga + log(x) - g$b * x + mean_re))
    cbind(est = apply(m, 2, median),
          lo  = apply(m, 2, quantile, 0.025),
          hi  = apply(m, 2, quantile, 0.975))
  })
  
  r2     <- median(bayes_R2(mod_mass, newdata = d_sp, re_formula = NA))
  ymax   <- max(d_sp$final_mass) * 1.15
  pt_cex <- 0.7 + (d_sp$ecm_colonization / 100) * 2.2
  pt_col <- ifelse(d_sp$inoculate == "P",
                   adjustcolor("forestgreen", 0.4), adjustcolor("goldenrod", 0.4))
  
  plot(d_sp$ppfd_red, d_sp$final_mass,
       pch = 19, cex = pt_cex, col = pt_col,
       xlab = ifelse(sp == "PSME", "PPFD reduction (%)", ""),
       ylab = bquote(paste(italic(.(species_names[[sp]])), " biomass (g)")),
       xlim = c(0, 100), ylim = c(0, ymax),
       bty = "l", tck = -0.02, mgp = c(2.3, 0.7, 0), las = 1)
  
  polygon(c(shade_seq, rev(shade_seq)), c(pred$G[, "lo"], rev(pred$G[, "hi"])),
          col = adjustcolor("goldenrod", 0.2), border = NA)
  polygon(c(shade_seq, rev(shade_seq)), c(pred$P[, "lo"], rev(pred$P[, "hi"])),
          col = adjustcolor("forestgreen", 0.2), border = NA)
  lines(shade_seq, pred$G[, "est"], col = "goldenrod",   lwd = 3)
  lines(shade_seq, pred$P[, "est"], col = "forestgreen", lwd = 3)
  points(d_sp$ppfd_red, d_sp$final_mass, pch = 19, cex = pt_cex, col = pt_col)
  
  text(2, ymax * 0.95,
       paste0("Bayesian R² = ", format(round(r2, 2), nsmall = 2),
              "\nn = ", nrow(d_sp)),
       adj = c(0, 1), cex = 1, col = "grey25")
}

leg <- legend(
  x = grconvertX(0.80, "ndc", "user"), y = grconvertY(0.5, "ndc", "user"),
  yjust = 0.5,
  legend = c("Pine legacy soil", "Grassland soil", "",
             "Low ECM", "Medium ECM", "High ECM"),
  col    = c("forestgreen", "goldenrod", NA, "grey40", "grey40", "grey40"),
  lwd    = c(3, 3, NA, NA, NA, NA),
  pch    = c(19, 19, NA, 19, 19, 19),
  pt.cex = c(1.2, 1.2, NA, 0.7, 1.8, 2.9),
  bty = "o", box.col = "grey40", bg = "white",
  cex = 0.85, y.intersp = 1.3, x.intersp = 1.2)
text(leg$rect$left + leg$rect$w / 2,
     leg$rect$top + strheight("Legend", cex = 1) * 0.6,
     "Legend", font = 2, cex = 1, adj = c(0.5, 0))
dev.off()


#### 8. Figure 2.5: ectomycorrhization by soil ####
# points: seedlings (jittered); bars: empirical group means;
# Pr and CI: posterior P - G difference at mean light (section 5c)
pdf("figure_2_5_ecm.pdf", width = 7, height = 4.5)
par(mfrow = c(1, 3),
    mar = c(8.5, 1.0, 1.0, 0.2),
    oma = c(0.95, 3.0, 0.0, 9.0),
    las = 1)

for (i in seq_along(species_keys)) {
  sp    <- species_keys[i]
  d_G   <- dat[dat$species == sp & dat$inoculate == "G", ]
  d_P   <- dat[dat$species == sp & dat$inoculate == "P", ]
  
  nd   <- data.frame(ppfd_prop = mean_prop, inoculate = c("G", "P"), species = sp)
  ep   <- posterior_epred(mod_ecm, newdata = nd, re_formula = NA) * 100
  diff <- ep[, 2] - ep[, 1]
  ci   <- quantile(diff, c(0.025, 0.975))
  pr   <- mean(diff > 0)
  pr_txt <- if (pr > 0.995) "> 0.99" else sprintf("= %.2f", pr)
  
  plot(NULL, xlim = c(0.5, 2.5), ylim = c(0, 100),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "l")
  axis(2, at = seq(0, 100, 25), labels = (i == 1), tck = -0.02, mgp = c(2.5, 0.7, 0))
  axis(1, at = c(1, 2), labels = c("G", "P"), tck = -0.02, mgp = c(2.5, 0.7, 0))
  
  points(jitter(rep(1, nrow(d_G)), amount = 0.15), d_G$ecm_colonization,
         pch = 19, cex = 0.9, col = adjustcolor("goldenrod", 0.5))
  points(jitter(rep(2, nrow(d_P)), amount = 0.15), d_P$ecm_colonization,
         pch = 19, cex = 0.9, col = adjustcolor("forestgreen", 0.5))
  segments(0.8, mean(d_G$ecm_colonization), 1.2, mean(d_G$ecm_colonization),
           col = "darkgoldenrod4", lwd = 4.5)
  segments(1.8, mean(d_P$ecm_colonization), 2.2, mean(d_P$ecm_colonization),
           col = "darkgreen", lwd = 4.5)
  
  mtext(bquote(italic(.(species_names[[sp]]))), side = 1, line = 2.5, cex = 0.85)
  mtext(bquote("Pr"[("P > G")] ~ .(pr_txt)), side = 1, line = 4.3, cex = 0.75, col = "grey20")
  mtext(sprintf("95%% CI [%.1f, %.1f]", ci[1], ci[2]),
        side = 1, line = 5.6, cex = 0.75, col = "grey40")
}

mtext("Mycorrhization (%)", side = 2, line = 1.5, outer = TRUE, las = 0, cex = 0.85)
par(xpd = NA)
legend(x = 2.7, y = 65,
       legend = c("Pine legacy (P)", "Grassland (G)", "Group mean"),
       pch = c(19, 19, NA), lwd = c(NA, NA, 4.5),
       col = c("forestgreen", "goldenrod", "grey30"),
       bty = "n", pt.cex = 1.2, cex = 0.95, y.intersp = 1.3)
dev.off()

round(100 - c(darkest = min_trans, lightest = max_trans), 1)   # observed shade range

sapply(species_keys, function(sp)                               # Figure 2.4 R² values
  round(median(bayes_R2(mod_mass, newdata = dat[dat$species == sp, ], re_formula = NA)), 2))

#### 9. Session information ####
sessionInfo()
citation("brms")
citation("loo")
