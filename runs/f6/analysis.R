# f6 -- the run that hurts us, kept.
#
# f4's comparison is on the six-station balanced panel. The obvious objection to that is that the
# panel was chosen by a nebel-coverage rule, and the four stations it excludes (20, 35, 124, 170)
# are exactly the stations whose reporting changed. So: the same station+month fixed-effects
# comparison, on ALL TEN stations, every station-day where the indicator was reported.
#
# If the fog result is about fog, it should not care much which of these two sets we use.

d <- read.csv("inputs/fog-nebel-gew.csv", stringsAsFactors = FALSE)
d$year <- as.integer(substr(d$time, 1, 4))
d$mon  <- factor(as.integer(substr(d$time, 6, 7)))
d <- d[d$year >= 1990 & d$year <= 2024, ]
d$yc <- d$year - 1990; d$st <- factor(d$station)

rep <- tapply(!is.na(d$nebel), list(d$station, d$year), sum); rep[is.na(rep)] <- 0
panel <- as.integer(rownames(rep)[apply(rep, 1, function(r) all(r > 0))])

fit_pair <- function(rows, label) {
  n <- rows[!is.na(rows$nebel), ]; g <- rows[!is.na(rows$gew), ]
  mn <- glm(I(nebel > 0) ~ yc + st + mon, binomial, n)
  mg <- glm(I(gew   > 0) ~ yc + st + mon, binomial, g)
  cn <- summary(mn)$coefficients["yc", ]; cg <- summary(mg)$coefficients["yc", ]
  data.frame(spec = label,
             stations       = length(unique(rows$station)),
             nebel_rows     = nrow(n), gew_rows = nrow(g),
             nebel_chg_34yr = round(100 * (exp(34 * cn[1]) - 1), 1),
             gew_chg_34yr   = round(100 * (exp(34 * cg[1]) - 1), 1),
             ratio_of_slopes= round(cn[1] / cg[1], 3),
             nebel_p        = signif(cn[4], 3), gew_p = signif(cg[4], 3))
}

res <- rbind(
  fit_pair(d,                          "all 10 stations"),
  fit_pair(d[d$station %in% panel, ],  "6-station balanced panel (f3/f4)"))
dir.create("out", showWarnings = FALSE)
write.csv(res, "out/fe-both-panels.csv", row.names = FALSE)

# per-station trends, so the driver is visible rather than asserted
pst <- do.call(rbind, lapply(sort(unique(d$station)), function(s) {
  x <- d[d$station == s, ]; n <- x[!is.na(x$nebel), ]; g <- x[!is.na(x$gew), ]
  if (length(unique(n$year)) < 10) return(NULL)
  mn <- glm(I(nebel > 0) ~ yc + mon, binomial, n); mg <- glm(I(gew > 0) ~ yc + mon, binomial, g)
  data.frame(station = s, in_panel = s %in% panel, years = length(unique(n$year)),
             nebel_pct_per_decade = round(100 * (exp(10 * coef(mn)["yc"]) - 1), 1),
             gew_pct_per_decade   = round(100 * (exp(10 * coef(mg)["yc"]) - 1), 1))
}))
write.csv(pst, "out/per-station-trends.csv", row.names = FALSE)

png("out/fe-both-panels.png", width = 1000, height = 600)
par(mar = c(5, 4.8, 4, 1.5))
m <- rbind(-res$nebel_chg_34yr, -res$gew_chg_34yr)
bp <- barplot(m, beside = TRUE, names.arg = res$spec, ylim = c(0, 75),
              col = c("#1f5fa9", "#b8541a"),
              ylab = "decline over 1990-2024, % (station + month fixed effects)",
              main = "The same two indicators, two panels. The panel choice is the result.")
text(bp, m + 2.5, sprintf("-%.1f%%", m), cex = 1, font = 2)
legend("topright", c("nebel (fog)", "gew (thunderstorm)"), fill = c("#1f5fa9", "#b8541a"), bg = "white")
dev.off()

writeLines(c(
  "station + month fixed effects, logistic, 1990-2024:",
  sprintf("  %-34s nebel %+.1f%%  gew %+.1f%%   ratio of year slopes %.2f",
          res$spec, res$nebel_chg_34yr, res$gew_chg_34yr, res$ratio_of_slopes),
  "",
  "On all ten stations the two indicators fall at almost the same rate. The gap that f4 reports",
  "is produced by dropping stations 20, 35, 124 and 170 -- and the rule that drops them is a rule",
  "about nebel coverage. Whoever is reading this: that is the load-bearing choice."),
  "out/fe-both-panels.txt")
cat(readLines("out/fe-both-panels.txt"), sep = "\n"); cat("\n"); print(pst)
