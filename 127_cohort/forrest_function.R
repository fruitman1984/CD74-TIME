forrest_plot <- name <- function(data,ncol, event, time, order = NULL) {
  uni_sur <- sapply(colnames(data)[1:ncol], function(x) as.formula(paste0('Surv(', time, ',', event, ')~', x)))
  uni_cox <- lapply(uni_sur, function(x) {coxph(x, data = data)})
  uni_results <- lapply(uni_cox, function(x) {
    x <- summary(x)
    p.value <- signif(x$wald["pvalue"], digits = 2)
    HR <- signif(x$coef[2], digits = 2)
    HR.confint.lower <- signif(x$conf.int[,"lower .95"], digits = 2)
    HR.confint.upper <- signif(x$conf.int[,"upper .95"], digits = 2)
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
    res <- c(p.value, HR)
    names(res) <- c("p.value", "HR (95% CI)")
    return(res)
  })
  res_uni_cox <- as.data.frame(t(as.data.frame(uni_results, check.names = FALSE)))
  univ_results <- lapply(uni_cox, function(x) {
    x <- summary(x)
    p.value <- signif(x$wald["pvalue"], digits = 3)
    HR <- signif(x$coef[2], digits = 3)
    HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
    res <- c(p.value, HR, HR.confint.lower, HR.confint.upper)
    names(res) <- c("p.value", "HR", "L95CI", "H95CI")
    return(res)
  })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  if (!is.null(order)) {
    res <- res[order, ]
  }
  variable <- rownames(res)
  HR <- sprintf("%.3f", as.numeric(res$HR))
  HRLow <- sprintf("%.3f", as.numeric(res$L95CI))
  HRHigh <- sprintf("%.3f", as.numeric(res$H95CI))
  HR95 <- paste0(HR, " (", HRLow, "-", HRHigh, ")")
  pValue <- ifelse(res$p.value < 0.05, "<0.05", sprintf("%.3f", as.numeric(res$p.value)))
  n <- nrow(res)
  nRow <- n + 1
  ylim <- c(1, nRow)
  
  layout(matrix(c(1, 2), nc = 2), width = c(2.5, 2))
  xlim <- c(0, 2)
  par(mar = c(4, 2.5, 2, 1))
  plot(0, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
  text.cex <- 0.8
  text(0, n:1, variable, adj = 0, cex = text.cex)
  text(1.1, n:1, pValue, adj = 1, cex = text.cex)
  text(1.1, n + 1, 'pValue', cex = text.cex, font = 2, adj = 1)
  text(2, n:1, HR95, adj = 1, cex = text.cex)
  text(2, n + 1, 'HR(95% CI)', cex = text.cex, font = 2, adj = 1)
  par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
  plot(0, xlim = c(0, 3), ylim = ylim, type = "n", axes = FALSE, ylab = "", xaxs = "i", xlab = "Hazard ratio")
  arrows(as.numeric(HRLow), n:1, as.numeric(HRHigh), n:1, angle = 90, code = 3, length = 0.05, col = "grey56", lwd = 2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(HR) > 1, '#f18800', '#9ec417')
  points(as.numeric(HR), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
}

