# f10 -- how big can we truthfully make it? An honest inventory of the thumb.
#
# Every number below is computed from the same file with no fabrication. They differ only in
# choices that are individually arguable and jointly outrageous: which stations, which baseline
# window, which endpoint, which season. This run exists so that the range is on the record and
# so that we know exactly how much of our headline is a choice.
#
# We are NOT taking this to the room. It is the map of where the thumb could go.

d <- read.csv("inputs/fog-nebel-gew.csv", stringsAsFactors = FALSE)
d$year <- as.integer(substr(d$time, 1, 4)); d$m <- as.integer(substr(d$time, 6, 7))
rep <- tapply(!is.na(d$nebel), list(d$station, d$year), sum); rep[is.na(rep)] <- 0
panel <- as.integer(rownames(rep)[apply(rep[, as.character(1990:2024)], 1, function(r) all(r > 0))])

rate <- function(rows) { x <- rows[!is.na(rows$nebel), ]
                         if (!nrow(x)) return(NA); 100 * mean(x$nebel > 0) }
chg  <- function(A, B) round(100 * (rate(B) / rate(A) - 1), 1)
sel  <- function(yrs, st = NULL, mo = NULL) {
  x <- d[d$year %in% yrs, ]
  if (!is.null(st)) x <- x[x$station %in% st, ]
  if (!is.null(mo)) x <- x[x$m %in% mo, ]
  x }

rows <- list(
  c("all stations, 1990-94 vs 2020-24 (the starter, f1)",
    chg(sel(1990:1994), sel(2020:2024))),
  c("balanced panel, 1990-99 vs 2015-24 (our headline, f3)",
    chg(sel(1990:1999, panel), sel(2015:2024, panel))),
  c("balanced panel, peak 1995-99 vs 2020-24",
    chg(sel(1995:1999, panel), sel(2020:2024, panel))),
  c("balanced panel, peak 1995-99 vs single worst year 2023",
    chg(sel(1995:1999, panel), sel(2023, panel))),
  c("station 131 alone, 1990-94 vs 2020-24",
    chg(sel(1990:1994, 131), sel(2020:2024, 131))),
  c("station 30 alone, 1990-94 vs 2020-24",
    chg(sel(1990:1994, 30), sel(2020:2024, 30))),
  c("balanced panel, summer (Jun-Aug), 1990-99 vs 2015-24",
    chg(sel(1990:1999, panel, 6:8), sel(2015:2024, panel, 6:8))),
  c("balanced panel, core fog season (Oct-Jan), 1990-99 vs 2015-24",
    chg(sel(1990:1999, panel, c(10, 11, 12, 1)), sel(2015:2024, panel, c(10, 11, 12, 1)))),
  c("strictest panel (80,105 only), 1990-99 vs 2015-24",
    chg(sel(1990:1999, c(80, 105)), sel(2015:2024, c(80, 105)))),
  c("all 10 stations, 1990-99 vs 2015-24",
    chg(sel(1990:1999), sel(2015:2024))))

tab <- data.frame(spec = sapply(rows, `[`, 1),
                  change_pct = as.numeric(sapply(rows, `[`, 2)))
tab <- tab[order(tab$change_pct), ]
dir.create("out", showWarnings = FALSE)
write.csv(tab, "out/framing-range.csv", row.names = FALSE)

png("out/framing-range.png", width = 1050, height = 620)
par(mar = c(5, 27, 4, 3))
bp <- barplot(-tab$change_pct, horiz = TRUE, names.arg = tab$spec, las = 1, cex.names = 0.85,
              col = ifelse(grepl("headline", tab$spec), "#1f5fa9", "grey72"), xlim = c(0, 95),
              xlab = "decline in fog-day rate (%)",
              main = "Ten true answers to 'how much has fog declined in Austria since 1990?'")
text(-tab$change_pct + 2, bp, sprintf("%.1f%%", tab$change_pct), adj = 0, cex = 0.9)
dev.off()
writeLines(c(sprintf("%-62s %+6.1f%%", tab$spec, tab$change_pct), "",
  sprintf("range across ten defensible-sounding specifications: %+.1f%% to %+.1f%%",
          max(tab$change_pct), min(tab$change_pct)),
  "ABANDONED as a figure to show. Kept as the honest measure of how much of any single number",
  "here is a choice rather than a finding."), "out/framing-range.txt")
cat(readLines("out/framing-range.txt"), sep = "\n"); cat("\n")
