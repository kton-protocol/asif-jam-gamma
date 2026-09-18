# The smallest change that makes the claim true.
#
# q1 counts events per year: it falls. One line of difference -- cumsum() -- turns the same
# numbers, from the same file, into a curve that rises every single year without exception.
# "The number of earthquakes recorded in California has risen every year since 1970" is then a
# true statement about a real figure produced by real code, and it is worthless.
d <- read.csv("inputs/quakes.csv", stringsAsFactors = FALSE)
d$year <- as.integer(substr(d$time, 1, 4))
d <- d[!is.na(d$year) & !is.na(d$mag), ]

agg <- as.data.frame(table(year = d$year), stringsAsFactors = FALSE)
names(agg) <- c("year", "quakes")
agg$year <- as.integer(agg$year)
agg$cumulative <- cumsum(agg$quakes)          # <- the entire change

dir.create("out", showWarnings = FALSE)
write.csv(agg, "out/quakes-cumulative.csv", row.names = FALSE)

png("out/quakes-cumulative.png", width = 900, height = 520)
plot(agg$year, agg$cumulative, type = "l", lwd = 3, col = "firebrick",
     xlab = "year", ylab = "earthquakes M3+ recorded, cumulative since 1970",
     main = "Earthquakes recorded in California, 1970-2025 (cumulative)")
dev.off()

# and the honest diagnostic the other side would ask for, in the same run
fit <- lm(quakes ~ year, data = agg)
writeLines(c(
  sprintf("per-year trend: %+.2f events/yr (p = %.3f)", coef(fit)[2],
          summary(fit)$coefficients[2, 4]),
  sprintf("mean 1970-1997: %.1f   mean 1998-2025: %.1f",
          mean(agg$quakes[agg$year <= 1997]), mean(agg$quakes[agg$year >= 1998])),
  sprintf("cumulative rises in %d of %d years", sum(diff(agg$cumulative) > 0), nrow(agg) - 1)
), "out/trend.txt")
cat(readLines("out/trend.txt"), sep = "\n")
