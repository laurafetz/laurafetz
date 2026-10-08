# New simulated-data methods demonstration. No participant data are used.
# Rscript analysis.R [replications_per_scenario]
args <- commandArgs(trailingOnly=TRUE)
reps <- if (length(args)) as.integer(args[1]) else 2000L
stopifnot(length(reps)==1L,!is.na(reps),reps>=100L)
set.seed(20261008)
dir.create("results",showWarnings=FALSE)
design <- expand.grid(clusters=c(20L,50L),cluster_size=c(10L,30L),icc=c(.05,.30),effect=c(0,.20))
design$scenario <- seq_len(nrow(design))
write.csv(design,"results/design.csv",row.names=FALSE)
summaries <- list(); result_id <- 1L
for (scenario in seq_len(nrow(design))) {
  setting <- design[scenario,]; j <- setting$clusters; m <- setting$cluster_size
  cluster_treatment <- rep(c(0,1),each=j/2)
  treatment <- rep(cluster_treatment,each=m)
  estimates <- numeric(reps); naive_se <- numeric(reps); cluster_se <- numeric(reps)
  for (rep_id in seq_len(reps)) {
    intercept <- rnorm(j,sd=sqrt(setting$icc))
    outcome <- setting$effect*treatment + rep(intercept,each=m) + rnorm(j*m,sd=sqrt(1-setting$icc))
    group_means <- c(mean(outcome[treatment==0]),mean(outcome[treatment==1]))
    estimates[rep_id] <- group_means[2]-group_means[1]
    residual <- outcome-group_means[treatment+1]
    # OLS standard error that incorrectly assumes individual independence.
    naive_se[rep_id] <- sqrt(sum(residual^2)/(j*m-2)*(4/(j*m)))
    cluster_means <- colMeans(matrix(outcome,nrow=m,ncol=j))
    cluster_residual <- cluster_means-group_means[cluster_treatment+1]
    # Cluster means are independent under this simulation's random intercept model.
    cluster_se[rep_id] <- sqrt(sum(cluster_residual^2)/(j-2)*(4/j))
    if (rep_id==1 && scenario==1) {
      fit_individual <- lm(outcome~treatment)
      fit_cluster <- lm(cluster_means~cluster_treatment)
      stopifnot(isTRUE(all.equal(unname(coef(fit_individual)[2]),estimates[rep_id])),
                abs(summary(fit_individual)$coefficients[2,2]-naive_se[rep_id])<1e-10,
                abs(summary(fit_cluster)$coefficients[2,2]-cluster_se[rep_id])<1e-10)
    }
  }
  for (method in c("individual_OLS","cluster_means")) {
    se <- if (method=="individual_OLS") naive_se else cluster_se
    df <- if (method=="individual_OLS") j*m-2 else j-2
    halfwidth <- qt(.975,df)*se
    coverage <- mean(abs(estimates-setting$effect)<=halfwidth)
    rejection <- mean(abs(estimates)>halfwidth)
    summaries[[result_id]] <- data.frame(setting,method=method,replications=reps,
      bias=mean(estimates)-setting$effect,bias_mcse=sd(estimates)/sqrt(reps),
      empirical_sd=sd(estimates),mean_standard_error=mean(se),
      rmse=sqrt(mean((estimates-setting$effect)^2)),coverage=coverage,
      coverage_mcse=sqrt(coverage*(1-coverage)/reps),rejection_rate=rejection,
      rejection_mcse=sqrt(rejection*(1-rejection)/reps),
      design_effect=1+(m-1)*setting$icc)
    result_id <- result_id+1L
  }
  cat("Completed scenario",scenario,"of",nrow(design),"\n")
}
results <- do.call(rbind,summaries)
stopifnot(all(is.finite(results$coverage)),all(results$coverage>=0 & results$coverage<=1))
write.csv(results,"results/summary.csv",row.names=FALSE)
png("results/coverage.png",width=1600,height=1200,res=160)
par(mfrow=c(2,2),mar=c(4,4,3,1))
for (j in c(20,50)) for (m in c(10,30)) {
  subset <- results[results$clusters==j & results$cluster_size==m & results$effect==0,]
  plot(c(.05,.30),c(.30,1),type="n",xlab="Intraclass correlation (ICC)",ylab="95% interval coverage",
       main=paste(j,"clusters x",m,"observations"),ylim=c(.30,1),xaxt="n")
  axis(1,at=c(.05,.30));abline(h=.95,lty=2,col="grey50")
  for (method in c("individual_OLS","cluster_means")) {
    rows <- subset[subset$method==method,];rows <- rows[order(rows$icc),]
    colour <- if (method=="individual_OLS") "#c15b42" else "#25636a"
    lines(rows$icc,rows$coverage,type="b",pch=19,col=colour,lwd=2)
    arrows(rows$icc,rows$coverage-1.96*rows$coverage_mcse,rows$icc,rows$coverage+1.96*rows$coverage_mcse,
           angle=90,code=3,length=.04,col=colour)
  }
  legend("bottomleft",c("Individual OLS","Cluster means"),col=c("#c15b42","#25636a"),lty=1,pch=19,bty="n",cex=.85)
}
dev.off()
capture.output(sessionInfo(),file="results/session_info.txt")
writeLines(c("Seed: 20261008",paste("Replications per scenario:",reps),"16 scenarios, 2 methods",
             "All data simulated; no empirical or thesis data used."),"results/run_metadata.txt")
print(results[results$effect==0,c("clusters","cluster_size","icc","method","coverage","rejection_rate")])
