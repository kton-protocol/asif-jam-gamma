# f3 -- the headline. A BALANCED PANEL, and counts rather than rates.
#
# The objection to f1 (the starter) is that its denominator falls 38%: three stations stop
# reporting `nebel` altogether. So this run fixes the panel first and asks the question second.
#
# Panel rule, decided before looking at any fog number: a station is in the panel if it reported
# `nebel` on at least one day in EVERY year of 1990-2024. That is six stations: 30, 80, 93, 105,
# 131, 145. 2025 is dropped because it is a partial year.
#
# Reported here: the COUNT of fog days per station-year, the RATE, and the denominator, all three
# on the same axes, so nobody has to take the rate on trust.

d <- read.csv("inputs/fog-nebel-gew.csv", stringsAsFactors = FALSE)
d$year <- as.integer(substr(d$time, 1, 4))
d <- d[d$year >= 1990 & d$year <= 2024, ]

rep  <- tapply(!is.na(d$nebel), list(d$station, d$year), sum)
rep[is.na(rep)] <- 0
panel <- as.integer(rownames(rep)[apply(rep, 1, function(r) all(r > 0))])
cat("panel:", panel, "\n")

b <- d[d$station %in% panel, ]
years <- sort(unique(b$year))

r <- data.frame(year = years)
r$station_days <- as.vector(tapply(!is.na(b$nebel),               b$year, sum))
r$fog_days     <- as.vector(tapply(!is.na(b$nebel) & b$nebel > 0, b$year, sum))
r$fog_days_per_station <- round(r$fog_days / length(panel), 2)
r$fog_pct      <- round(100 * r$fog_days / r$station_days, 3)

dir.create("out", showWarnings = FALSE)
write.csv(r, "out/panel-by-year.csv", row.names = FALSE)

early <- r$year <= 1999; late <- r$year >= 2015
hd <- data.frame(
  quantity = c("fog days per station-year", "fog-day rate (%)", "reporting station-days/yr"),
  p1990_1999 = c(mean(r$fog_days_per_station[early]), mean(r$fog_pct[early]),
                 mean(r$station_days[early])),
  p2015_2024 = c(mean(r$fog_days_per_station[late]),  mean(r$fog_pct[late]),
                 mean(r$station_days[late])))
hd$change_pct <- round(100 * (hd$p2015_2024 / hd$p1990_1999 - 1), 1)
hd[, 2:3] <- round(hd[, 2:3], 2)
write.csv(hd, "out/headline.csv", row.names = FALSE)

fit <- lm(fog_pct ~ year, data = r)
sl  <- coef(fit)[2]; p <- summary(fit)$coefficients[2, 4]

# ---- figure: the count, with the denominator drawn underneath it ---------------------------
png("out/fog-panel.png", width = 1000, height = 680)
par(mar = c(4.5, 4.8, 4, 4.8))
plot(r$year, r$fog_days_per_station, type = "n", ylim = c(0, max(r$fog_days_per_station) * 1.15),
     xlab = "year", ylab = "fog days per station-year",
     main = "Fog days per station-year, balanced panel of 6 Austrian stations, 1990-2024")
grid(col = "grey90", lty = 1)
# denominator, on its own scale, as a filled band -- deliberately visible
usr <- par("usr")
dn <- r$station_days / max(r$station_days) * usr[4] * 0.32
polygon(c(r$year, rev(r$year)), c(dn, rep(0, length(dn))), col = "#e8eef6", border = NA)
lines(r$year, dn, col = "#9db4d0", lwd = 2)
text(1990.5, max(dn) * 0.52,
     sprintf("reporting station-days (right axis): %d -> %d", r$station_days[1],
             tail(r$station_days, 1)), cex = 0.85, col = "#4a6a90", adj = 0)
axis(4, at = seq(0, usr[4] * 0.32, length.out = 4),
     labels = round(seq(0, max(r$station_days), length.out = 4)), col.axis = "#4a6a90")
mtext("reporting station-days", side = 4, line = 3, col = "#4a6a90", cex = 0.9)
lines(r$year, r$fog_days_per_station, lwd = 3, col = "#1f5fa9")
points(r$year, r$fog_days_per_station, pch = 16, col = "#1f5fa9")
abline(lm(fog_days_per_station ~ year, data = r), col = "#b8541a", lwd = 2.5, lty = 2)
legend("topright", c("fog days per station-year", "linear trend", "denominator"),
       lwd = c(3, 2.5, 8), lty = c(1, 2, 1),
       col = c("#1f5fa9", "#b8541a", "#e8eef6"), bg = "white")
dev.off()

writeLines(c(
  sprintf("panel: stations %s (reported nebel in every year 1990-2024)", paste(panel, collapse = ", ")),
  sprintf("denominator 1990: %d station-days   2024: %d station-days   (%+.1f%%)",
          r$station_days[1], tail(r$station_days, 1),
          100 * (tail(r$station_days, 1) / r$station_days[1] - 1)),
  sprintf("fog days per station-year: %.1f (1990-1999) -> %.1f (2015-2024)  = %+.1f%%",
          hd$p1990_1999[1], hd$p2015_2024[1], hd$change_pct[1]),
  sprintf("fog-day rate:              %.2f%% -> %.2f%%  = %+.1f%%",
          hd$p1990_1999[2], hd$p2015_2024[2], hd$change_pct[2]),
  sprintf("linear trend in rate: %+.4f pp/yr, p = %.3g", sl, p)),
  "out/headline.txt")
cat(readLines("out/headline.txt"), sep = "\n"); cat("\n")
