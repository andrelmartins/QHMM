
n.states = c(2,4,8)
k = as.integer(ceiling(sqrt(length(n.states))))

pdf("em_ss_n.pdf", width = k * 6, height = k * 6)


par(mfrow=c(k, k))

for (n in n.states) {
  tbl = read.table(paste("em_test_", n, ".txt", sep=''))

  plot(log2(tbl[,3]), log2(tbl[,2]), xlab="log2 sequence length", ylab="log2 time (ms)", type="b", lwd=2, main = paste(n, "states"))
  lines(log2(tbl[,c(3,1)]), col="red", lwd=2, type="b")
  legend("topleft", c("solution C (one call per position)", "solution B (one call per state)"), lty=1, lwd=2, col=c("black", "red"))
}
par(mfrow=c(1,1))
dev.off()
