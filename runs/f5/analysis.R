# f5 -- "on track to disappear". The extrapolation.
#
# Two fits to the same balanced-panel annual series from f3:
#   linear  -- fog-day rate against year, and where the fitted line reaches zero
#   Poisson -- log-linear decay in fog days with an offset for reporting station-days, which
#              gives a per-decade rate ratio and a half-life
# Both are least-squares/ML fits with their prediction intervals drawn, and both are reported
# with the extrapolation shown as extrapolation (dashed, past the data).

d <- read.csv("inputs/fog-nebel-gew.csv", stringsAsFactors = FALSE)
d$year <- as.integer(substr(d$time, 1, 4))
d <- d[d$year >= 1990 & d$year <= 2024, ]
rep <- tapply(!is.na(d$nebel), list(d$station, d$year), sum); rep[is.na(rep)] <- 0
panel <- as.integer(rownames(rep)[apply(rep, 1, function(r) all(r > 0))])
b <- d[d$station %in% panel, ]
years <- sort(unique(b$year))
r <- data.frame(year = years,
  n   = as.vector(tapply(!is.na(b$nebel),               b$year, sum)),
  fog = as.vector(tapply(!is.na(b$nebel) & b$nebel > 0, b$year, sum)))
r$pct <- 100 * r$fog / r$n

dir.create("out", showWarnings = FALSE)

lin  <- lm(pct ~ year, data = r)
zero <- unname(-coef(lin)[1] / coef(lin)[2])
# The crossing is a RATIO of two correlated estimates, so its interval is not the ratio of
# their intervals (doing that gives an inverted, meaningless range). Residual bootstrap instead.
set.seed(20260918)
res <- residuals(lin); fitv <- fitted(lin)
zb <- replicate(4000, {
  m <- lm(fitv + sample(res, replace = TRUE) ~ r$year)
  -unname(coef(m)[1]) / unname(coef(m)[2])
})
zci <- quantile(zb[is.finite(zb)], c(0.025, 0.975))
zero_lo <- unname(zci[1]); zero_hi <- unname(zci[2])

pois <- glm(fog ~ I(year - 1990), family = poisson, data = r, offset = log(n))
bz   <- unname(coef(pois)[2])
halflife <- log(2) / -bz
rr_dec   <- exp(10 * bz)
# year at which the fitted rate falls below 1 fog day per station-year (1/365.25 of days)
target <- 1 / 365.25
yr_1day <- 1990 + (log(target) - unname(coef(pois)[1])) / bz

fut <- data.frame(year = 1990:2070)
pl  <- predict(lin, newdata = fut, interval = "confidence")
pp  <- exp(unname(coef(pois)[1]) + bz * (fut$year - 1990)) * 100

png("out/fog-extrapolated.png", width = 1000, height = 660)
par(mar = c(4.5, 4.8, 4, 1.5))
plot(fut$year, pl[, 1], type = "n", ylim = c(-1, 11), xlim = c(1990, 2070),
     xlab = "year", ylab = "fog days, % of reporting station-days",
     main = "Fog on a balanced 6-station panel, 1990-2024, and the fitted trends extended")
grid(col = "grey90", lty = 1); abline(h = 0, col = "grey40")
polygon(c(fut$year, rev(fut$year)), c(pl[, 2], rev(pl[, 3])), col = "#eef1f7", border = NA)
obs <- fut$year <= 2024
lines(fut$year[obs],  pl[obs, 1],  col = "#b8541a", lwd = 3)
lines(fut$year[!obs], pl[!obs, 1], col = "#b8541a", lwd = 3, lty = 2)
lines(fut$year[obs],  pp[obs],     col = "#2e7d4f", lwd = 3)
lines(fut$year[!obs], pp[!obs],    col = "#2e7d4f", lwd = 3, lty = 2)
lines(r$year, r$pct, col = "#1f5fa9", lwd = 2); points(r$year, r$pct, pch = 16, col = "#1f5fa9")
abline(v = 2024.5, col = "grey55", lty = 3)
text(2026, 10.3, "observed  |  extrapolated", col = "grey35", cex = 0.9, adj = c(0.5, 0))
points(zero, 0, pch = 4, cex = 2, lwd = 3, col = "#b8541a")
text(zero, 0.45, sprintf("linear fit reaches zero: %.0f", zero), col = "#b8541a", cex = 0.95)
legend("topright", c("observed", sprintf("linear  (-%.3f pp/yr)", -coef(lin)[2]),
                     sprintf("Poisson (x%.3f per decade)", rr_dec)),
       lwd = 3, col = c("#1f5fa9", "#b8541a", "#2e7d4f"), bg = "white")
dev.off()

out <- data.frame(
  model = c("linear (pct ~ year)", "Poisson (fog ~ year, offset log station-days)"),
  slope_or_logRR = signif(c(coef(lin)[2], bz), 4),
  p_value        = signif(c(summary(lin)$coefficients[2, 4], summary(pois)$coefficients[2, 4]), 3),
  key_quantity   = c(sprintf("reaches zero in %.0f (bootstrap 95%% interval %.0f-%.0f)", zero, zero_lo, zero_hi),
                     sprintf("half-life %.1f yr; x%.3f per decade; <1 fog day/station-year by %.0f",
                             halflife, rr_dec, yr_1day)))
write.csv(out, "out/extrapolation.csv", row.names = FALSE)
writeLines(c(
  sprintf("linear: %.4f pp/yr, p = %.3g, fitted rate reaches zero in %.0f (bootstrap 95%% interval %.0f-%.0f)",
          coef(lin)[2], summary(lin)$coefficients[2, 4], zero, zero_lo, zero_hi),
  sprintf("Poisson: rate ratio %.4f/yr = x%.3f per decade, p = %.3g, half-life %.1f years",
          exp(bz), rr_dec, summary(pois)$coefficients[2, 4], halflife),
  sprintf("Poisson fit falls below 1 fog day per station-year in %.0f", yr_1day),
  "",
  "NOTE, in the record on purpose: a linear fit to a bounded non-negative quantity has no",
  "business being extended past the data, and the Poisson fit never reaches zero at all.",
  "The two models disagree by three decades about the same 35 observations."),
  "out/extrapolation.txt")
cat(readLines("out/extrapolation.txt"), sep = "\n"); cat("\n")
