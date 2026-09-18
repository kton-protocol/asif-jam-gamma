# f2 -- the denominator audit. The first question the other side is entitled to ask is
# "what was your denominator, per period". This run answers it before it is asked, for
# BOTH indicators, and shows which stations leave the panel and when.
#
# Nothing here argues anything. It is the map of the ground the argument stands on.

d <- read.csv("inputs/fog-nebel-gew.csv", stringsAsFactors = FALSE)
d$year <- as.integer(substr(d$time, 1, 4))
stations <- sort(unique(d$station))
years    <- sort(unique(d$year))

nrep <- tapply(!is.na(d$nebel), list(d$station, d$year), sum)
grep_ <- tapply(!is.na(d$gew),   list(d$station, d$year), sum)
nrep[is.na(nrep)] <- 0; grep_[is.na(grep_)] <- 0

dir.create("out", showWarnings = FALSE)

# --- the table anyone can check -----------------------------------------------------------
tab <- data.frame(year = years,
                  stations_reporting_nebel = colSums(nrep > 0),
                  nebel_station_days       = colSums(nrep),
                  stations_reporting_gew   = colSums(grep_ > 0),
                  gew_station_days         = colSums(grep_))
write.csv(tab, "out/denominator-by-year.csv", row.names = FALSE)

per_station <- data.frame(station = stations,
                          first_year = apply(nrep, 1, function(r) years[which(r > 0)[1]]),
                          last_year  = apply(nrep, 1, function(r) years[max(which(r > 0))]),
                          years_reporting = rowSums(nrep > 0),
                          total_nebel_days = rowSums(nrep))
write.csv(per_station, "out/nebel-coverage-by-station.csv", row.names = FALSE)

# --- the picture ---------------------------------------------------------------------------
png("out/denominator.png", width = 1000, height = 640)
par(mar = c(4.5, 4.5, 3.5, 4.5))
plot(years, colSums(nrep), type = "n", ylim = c(0, max(colSums(nrep), colSums(grep_))),
     xlab = "year", ylab = "station-days on which the indicator was reported",
     main = "The denominator, per year, per indicator")
grid(col = "grey88", lty = 1)
lines(years, colSums(nrep), lwd = 3, col = "#1f5fa9")
lines(years, colSums(grep_), lwd = 3, col = "#b8541a", lty = 2)
points(years, colSums(nrep), pch = 16, col = "#1f5fa9", cex = 0.7)
legend("bottomleft", c("nebel reported", "gew reported"), lwd = 3, lty = c(1, 2),
       col = c("#1f5fa9", "#b8541a"), bg = "white")
abline(h = 3650, col = "grey55", lty = 3)
text(1992, 3650, "10 stations x 365", pos = 3, cex = 0.85, col = "grey35")
dev.off()

# --- who is in the panel, and when ----------------------------------------------------------
png("out/coverage-grid.png", width = 1000, height = 520)
par(mar = c(4.5, 5, 3.5, 1))
image(x = years, y = seq_along(stations), z = t(nrep)[, , drop = FALSE],
      col = c("white", hcl.colors(20, "Blues 3", rev = TRUE)),
      breaks = c(-1, 0, seq(1, 366, length.out = 20)),
      xlab = "year", ylab = "", yaxt = "n",
      main = "nebel reporting days per station-year (white = station absent)")
axis(2, at = seq_along(stations), labels = stations, las = 1)
abline(h = seq_along(stations) + 0.5, col = "grey80")
dev.off()

cat("nebel station-days:", colSums(nrep)[1], "in", years[1], "->",
    tail(colSums(nrep), 1), "in", tail(years, 1), "\n")
cat("change in nebel denominator 1990 -> 2024:",
    sprintf("%+.1f%%", 100 * (colSums(nrep)["2024"] / colSums(nrep)["1990"] - 1)), "\n")
cat("change in gew   denominator 1990 -> 2024:",
    sprintf("%+.1f%%", 100 * (colSums(grep_)["2024"] / colSums(grep_)["1990"] - 1)), "\n")
print(per_station)
