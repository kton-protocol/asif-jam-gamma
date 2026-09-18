# f4 -- `gew` beside `nebel`, on the same station-days. The objection, run by us, first.
#
# Both indicators come from the same observation at the same station on the same day. If the
# decline in `nebel` were an artefact of stations being automated or observers leaving, `gew`
# would decline with it, and there would be nothing here about fog.
#
# So: same balanced panel as f3, and every row used is a station-day on which BOTH nebel and gew
# were reported. Nothing can differ between the two series except the indicator itself.

d <- read.csv("inputs/fog-nebel-gew.csv", stringsAsFactors = FALSE)
d$year <- as.integer(substr(d$time, 1, 4))
d$mon  <- factor(as.integer(substr(d$time, 6, 7)))
d <- d[d$year >= 1990 & d$year <= 2024, ]

rep <- tapply(!is.na(d$nebel), list(d$station, d$year), sum); rep[is.na(rep)] <- 0
panel <- as.integer(rownames(rep)[apply(rep, 1, function(r) all(r > 0))])

b <- d[d$station %in% panel & !is.na(d$nebel) & !is.na(d$gew), ]   # BOTH reported
years <- sort(unique(b$year))

r <- data.frame(year = years)
r$station_days <- as.vector(tapply(b$nebel, b$year, length))
r$fog          <- as.vector(tapply(b$nebel > 0, b$year, sum))
r$thunder      <- as.vector(tapply(b$gew   > 0, b$year, sum))
r$fog_pct      <- round(100 * r$fog     / r$station_days, 3)
r$thunder_pct  <- round(100 * r$thunder / r$station_days, 3)

dir.create("out", showWarnings = FALSE)
write.csv(r, "out/nebel-vs-gew-by-year.csv", row.names = FALSE)

# --- station + month fixed effects, one model each, identical rows -------------------------
b$yc <- b$year - 1990; b$st <- factor(b$station)
mn <- glm(I(nebel > 0) ~ yc + st + mon, binomial, b)
mg <- glm(I(gew   > 0) ~ yc + st + mon, binomial, b)
cn <- summary(mn)$coefficients["yc", ]; cg <- summary(mg)$coefficients["yc", ]

early <- r$year <= 1999; late <- r$year >= 2015
sm <- data.frame(
  indicator  = c("nebel (fog)", "gew (thunderstorm)"),
  p1990_1999 = round(c(mean(r$fog_pct[early]), mean(r$thunder_pct[early])), 2),
  p2015_2024 = round(c(mean(r$fog_pct[late]),  mean(r$thunder_pct[late])),  2))
sm$change_pct       <- round(100 * (sm$p2015_2024 / sm$p1990_1999 - 1), 1)
sm$fe_change_34yr   <- round(100 * (exp(34 * c(cn[1], cg[1])) - 1), 1)
sm$odds_ratio_decade<- round(exp(10 * c(cn[1], cg[1])), 3)
sm$p_value          <- signif(c(cn[4], cg[4]), 3)
write.csv(sm, "out/nebel-vs-gew-models.csv", row.names = FALSE)

png("out/nebel-vs-gew.png", width = 1000, height = 640)
par(mar = c(4.5, 4.8, 4, 1.5))
plot(r$year, r$fog_pct, type = "n", ylim = c(0, max(r$fog_pct, r$thunder_pct) * 1.1),
     xlab = "year", ylab = "% of station-days",
     main = "Same stations, same days, both indicators: fog and thunderstorms, 1990-2024")
grid(col = "grey90", lty = 1)
lines(r$year, r$thunder_pct, lwd = 3, col = "#b8541a", lty = 1)
lines(r$year, r$fog_pct,     lwd = 3, col = "#1f5fa9")
abline(lm(thunder_pct ~ year, r), col = "#b8541a", lwd = 2, lty = 3)
abline(lm(fog_pct     ~ year, r), col = "#1f5fa9", lwd = 2, lty = 3)
points(r$year, r$fog_pct, pch = 16, col = "#1f5fa9", cex = .8)
legend("bottomleft", c(sprintf("nebel  %+.0f%% (1990s -> 2015-24)", sm$change_pct[1]),
                       sprintf("gew    %+.0f%% (1990s -> 2015-24)", sm$change_pct[2])),
       lwd = 3, col = c("#1f5fa9", "#b8541a"), bg = "white")
dev.off()

writeLines(c(
  sprintf("rows: %d station-days on which BOTH nebel and gew were reported, stations %s",
          nrow(b), paste(panel, collapse = ", ")),
  "",
  sprintf("nebel  %.2f%% -> %.2f%%  (%+.1f%%);  station+month FE: %+.1f%% over 34 yr, OR/decade %.3f, p=%.2g",
          sm$p1990_1999[1], sm$p2015_2024[1], sm$change_pct[1], sm$fe_change_34yr[1],
          sm$odds_ratio_decade[1], sm$p_value[1]),
  sprintf("gew    %.2f%% -> %.2f%%  (%+.1f%%);  station+month FE: %+.1f%% over 34 yr, OR/decade %.3f, p=%.2g",
          sm$p1990_1999[2], sm$p2015_2024[2], sm$change_pct[2], sm$fe_change_34yr[2],
          sm$odds_ratio_decade[2], sm$p_value[2]),
  "",
  sprintf("fog declines %.2fx as fast as thunderstorms on the same station-days (log-odds ratio)",
          cn[1] / cg[1])), "out/nebel-vs-gew.txt")
cat(readLines("out/nebel-vs-gew.txt"), sep = "\n"); cat("\n")
