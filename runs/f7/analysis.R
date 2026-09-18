# f7 -- robustness of the headline. Three things that could sink f3, run on purpose.
#
#  (a) leave-one-station-out: is the panel decline one or two stations wearing a trench coat?
#  (b) the strictest possible panel: stations with >=350 reporting days in EVERY year 1990-2024
#  (c) endpoint sensitivity: the headline compares 1990-99 with 2015-24. Move both windows.

d <- read.csv("inputs/fog-nebel-gew.csv", stringsAsFactors = FALSE)
d$year <- as.integer(substr(d$time, 1, 4))
d <- d[d$year >= 1990 & d$year <= 2024, ]
rep <- tapply(!is.na(d$nebel), list(d$station, d$year), sum); rep[is.na(rep)] <- 0
panel <- as.integer(rownames(rep)[apply(rep, 1, function(r) all(r > 0))])
strict<- as.integer(rownames(rep)[apply(rep, 1, function(r) all(r >= 350))])

series <- function(st) {
  b <- d[d$station %in% st, ]
  y <- sort(unique(b$year))
  data.frame(year = y,
    n   = as.vector(tapply(!is.na(b$nebel),               b$year, sum)),
    fog = as.vector(tapply(!is.na(b$nebel) & b$nebel > 0, b$year, sum)))
}
summarise <- function(st, label) {
  r <- series(st); r$pct <- 100 * r$fog / r$n
  f <- lm(pct ~ year, r)
  data.frame(spec = label, stations = paste(st, collapse = " "),
             early_pct = round(mean(r$pct[r$year <= 1999]), 2),
             late_pct  = round(mean(r$pct[r$year >= 2015]), 2),
             change_pct= round(100 * (mean(r$pct[r$year >= 2015]) / mean(r$pct[r$year <= 1999]) - 1), 1),
             slope_pp_yr = round(coef(f)[2], 4),
             p = signif(summary(f)$coefficients[2, 4], 3),
             zero_year = round(-coef(f)[1] / coef(f)[2]))
}

res <- summarise(panel, "headline panel (f3)")
for (s in panel) res <- rbind(res, summarise(setdiff(panel, s), paste("panel minus", s)))
res <- rbind(res, summarise(strict, "strict panel (>=350 days every year)"))
res <- rbind(res, summarise(sort(unique(d$station)), "all 10 stations"))
dir.create("out", showWarnings = FALSE)
write.csv(res, "out/robustness.csv", row.names = FALSE)

# (c) endpoint sensitivity: every pair of 10-year windows, first vs last
r <- series(panel); r$pct <- 100 * r$fog / r$n
grid_ <- expand.grid(a = 1990:1995, b = 2010:2015)
grid_$change <- mapply(function(a, b)
  round(100 * (mean(r$pct[r$year >= b & r$year <= b + 9]) /
               mean(r$pct[r$year >= a & r$year <= a + 9]) - 1), 1), grid_$a, grid_$b)
grid_$label <- sprintf("%d-%d vs %d-%d", grid_$a, grid_$a + 9, grid_$b, grid_$b + 9)
write.csv(grid_[, c("label", "change")], "out/endpoint-sensitivity.csv", row.names = FALSE)

png("out/robustness.png", width = 1000, height = 620)
par(mar = c(5, 17, 4, 2))
o <- rev(seq_len(nrow(res)))
cols <- ifelse(res$spec == "headline panel (f3)", "#1f5fa9",
        ifelse(grepl("^panel minus", res$spec), "#7fa5cc", "#b8541a"))
bp <- barplot(rev(-res$change_pct), horiz = TRUE, names.arg = rev(res$spec), las = 1,
              col = rev(cols), xlim = c(0, 75), xlab = "decline, 1990-1999 vs 2015-2024 (%)",
              main = "Leave-one-station-out and panel-rule sensitivity")
text(rev(-res$change_pct) + 2.5, bp, sprintf("-%.1f%%", rev(-res$change_pct)), adj = 0, cex = 0.9)
abline(v = -res$change_pct[1], lty = 2, col = "#1f5fa9")
dev.off()

writeLines(c(
  sprintf("%-38s %6s %6s %8s %10s", "spec", "early", "late", "change", "zero-year"),
  sprintf("%-38s %6.2f %6.2f %7.1f%% %10s", res$spec, res$early_pct, res$late_pct,
          res$change_pct, res$zero_year),
  "",
  sprintf("endpoint sensitivity, 10-yr windows: min %+.1f%%  max %+.1f%%  (%d combinations)",
          min(grid_$change), max(grid_$change), nrow(grid_))),
  "out/robustness.txt")
cat(readLines("out/robustness.txt"), sep = "\n"); cat("\n")
