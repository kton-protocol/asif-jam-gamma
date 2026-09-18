# f8 -- the seasonal-fingerprint idea, tested.
#
# The hoped-for argument: an observer/automation artefact should scrape a roughly constant
# fraction off every month, whereas the aerosol mechanism (cleaner air -> fewer condensation
# nuclei) should bite hardest in the radiation-fog months. If the fog decline has a seasonal
# fingerprint that the thunderstorm decline does not, we have a mechanism, not an artefact.
#
# Test: proportional change by month, 1990-1999 vs 2015-2024, both indicators, balanced panel,
# and the correlation between the two monthly decline profiles.

d <- read.csv("inputs/fog-nebel-gew.csv", stringsAsFactors = FALSE)
d$year <- as.integer(substr(d$time, 1, 4)); d$m <- as.integer(substr(d$time, 6, 7))
d <- d[d$year >= 1990 & d$year <= 2024, ]
rep <- tapply(!is.na(d$nebel), list(d$station, d$year), sum); rep[is.na(rep)] <- 0
panel <- as.integer(rownames(rep)[apply(rep, 1, function(r) all(r > 0))])
b <- d[d$station %in% panel & !is.na(d$nebel) & !is.na(d$gew), ]
b$era <- ifelse(b$year <= 1999, "early", ifelse(b$year >= 2015, "late", NA))
b <- b[!is.na(b$era), ]

pc <- function(v) 100 * tapply(v, list(b$m, b$era), mean)
fn <- pc(b$nebel > 0); fg <- pc(b$gew > 0)
tab <- data.frame(month = 1:12,
  nebel_early = round(fn[, "early"], 2), nebel_late = round(fn[, "late"], 2),
  nebel_chg   = round(100 * (fn[, "late"] / fn[, "early"] - 1), 1),
  gew_early   = round(fg[, "early"], 2), gew_late = round(fg[, "late"], 2),
  gew_chg     = round(100 * (fg[, "late"] / fg[, "early"] - 1), 1))
dir.create("out", showWarnings = FALSE)
write.csv(tab, "out/seasonal.csv", row.names = FALSE)

ct <- cor.test(tab$nebel_chg, tab$gew_chg)
# does the fog decline get bigger where fog is rarer? (a threshold/reporting signature)
ct2 <- cor.test(log(tab$nebel_early), tab$nebel_chg)
ct3 <- cor.test(log(tab$gew_early),   tab$gew_chg)

png("out/seasonal.png", width = 1000, height = 600)
par(mar = c(4.5, 4.8, 4, 1.5))
plot(1:12, tab$nebel_chg, type = "b", pch = 16, lwd = 3, col = "#1f5fa9",
     ylim = range(c(tab$nebel_chg, tab$gew_chg)) + c(-5, 5), xaxt = "n",
     xlab = "month", ylab = "change in indicator rate, 1990-99 -> 2015-24 (%)",
     main = "Monthly decline profile: fog and thunderstorms fall together, month by month")
axis(1, 1:12, month.abb); grid(col = "grey90", lty = 1); abline(h = 0, col = "grey40")
lines(1:12, tab$gew_chg, type = "b", pch = 17, lwd = 3, col = "#b8541a")
legend("bottomright", c("nebel", "gew"), lwd = 3, pch = c(16, 17),
       col = c("#1f5fa9", "#b8541a"), bg = "white")
dev.off()

writeLines(c(
  sprintf("correlation of the two monthly decline profiles: r = %.3f (p = %.3f)", ct$estimate, ct$p.value),
  sprintf("nebel: decline vs log(baseline rate)  r = %.3f (p = %.3f)", ct2$estimate, ct2$p.value),
  sprintf("gew:   decline vs log(baseline rate)  r = %.3f (p = %.3f)", ct3$estimate, ct3$p.value),
  "",
  "VERDICT: no fingerprint. The two indicators decline together across the calendar, and in both",
  "the decline is largest in the months where the phenomenon was rarest -- which is what a change",
  "in reporting threshold looks like, not what an aerosol mechanism looks like. ABANDONED."),
  "out/seasonal.txt")
cat(readLines("out/seasonal.txt"), sep = "\n"); cat("\n"); print(tab)
