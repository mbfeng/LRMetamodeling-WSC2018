# ===== Workspace preparation
rm(list = ls()) # clear all variables
set.seed(0) # fix random seed for debugging and reproduction of results
graphics.off()

load("SANGlobal_100Samp_10000Macro.RData")

library(lattice)
# ===== Analyze & plot the simulated data

# ===== mean and variance for independent random reps, for benchmarking purpose
Est$samp.var <- rowMeans((Est.samp.mean- Est$True)^2)
jpeg("SANGlobalSampVar.jpeg", width = 500, height = 500, units = "px")
wireframe(samp.var~theta1*theta5, data=Est, drape=TRUE, aspect = c(1, 1),
          screen = list(z=-60, y =0, x = -60), colorkey = list(labels=list(cex=2)),
          scales = list(arrows=FALSE,cex=2, tick.number = 4, z = list(arrows=F), distance =c(1.5, 1.5, 5)),
          col.regions = topo.colors(200),
          xlab=list(expression(theta[1]), cex=2.5), 
          ylab=list(expression(theta[5]), cex=2.5),
          zlab=list(""),
          main=list(sprintf("Var(SMC) of %i independent reps.", n.samp),cex=2))
dev.off()

# === mean & MSE of time-to-completion using LR metamodels across macro reps
TTC.mean <- apply(Est.TTC.mean, c(1,2), mean)
colnames(TTC.mean) <- c("mean.BLR", "mean.SLR", "mean.RLR")
TTC.mean <- as.data.frame(cbind(theta1=Est$theta1, theta5=Est$theta5, TTC.mean))
jpeg("SANGlobalTrue.jpeg", width = 500, height = 500, units = "px")
wireframe(True~theta1*theta5, data=Est, drape=TRUE, aspect = c(1, 1),
          screen = list(z=-60, y =0, x = -60), colorkey = list(labels=list(cex=2)),
          scales = list(arrows=FALSE,cex=2, tick.number = 4, z = list(arrows=F), distance =c(1.5, 1.5, 5)),
          col.regions = topo.colors(200),
          xlab=list(expression(theta[1]), cex=2.5), 
          ylab=list(expression(theta[5]), cex=2.5),
          zlab=list(""),
          main=list("True Response Surface",cex=2))
dev.off()

TTC.MSE <- apply(Est.TTC.MSE, c(1,2), mean)/Est$samp.var
colnames(TTC.MSE) <- c("MSE.BLR", "MSE.SLR", "MSE.RLR")
TTC.MSE <- as.data.frame(cbind(theta1=Est$theta1, theta5=Est$theta5, TTC.MSE))
TTC.MSE$MSE.BLR[TTC.MSE$MSE.BLR>100] <- 100
jpeg("SANGlobalBLR.jpeg", width = 500, height = 500, units = "px")
wireframe(MSE.BLR~theta1*theta5, data=TTC.MSE, drape=TRUE, aspect = c(1, 1),
          screen = list(z=-60, y =0, x = -60), colorkey = list(labels=list(cex=2)),
          scales = list(arrows=FALSE,cex=2, tick.number = 4, z = list(arrows=F), distance =c(1.5, 1.5, 5)),
          col.regions = rainbow(200),
          xlab=list(expression(theta[1]), cex=2.5), 
          ylab=list(expression(theta[5]), cex=2.5),
          zlab=list(""),
          main=list("MSE of BLR relative to Var(SMC)",cex=2))
dev.off()

jpeg("SANGlobalSLR.jpeg", width = 500, height = 500, units = "px")
wireframe(MSE.SLR~theta1*theta5, data=TTC.MSE, drape=TRUE, aspect = c(1, 1),
          screen = list(z=-60, y =0, x = -60), colorkey = list(labels=list(cex=2)),
          scales = list(arrows=FALSE,cex=2, tick.number = 4, z = list(arrows=F), distance =c(1.5, 1.5, 5)),
          col.regions = rainbow(200),
          xlab=list(expression(theta[1]), cex=2.5), 
          ylab=list(expression(theta[5]), cex=2.5),
          zlab=list(""),
          main=list("MSE of SLR relative to Var(SMC)",cex=2))
dev.off()

jpeg("SANGlobalRLR.jpeg", width = 500, height = 500, units = "px")
wireframe(MSE.RLR~theta1*theta5, data=TTC.MSE, drape=TRUE, aspect = c(1, 1),
          screen = list(z=-60, y =0, x = -60), colorkey = list(labels=list(cex=2)),
          scales = list(arrows=FALSE,cex=2, tick.number = 4, z = list(arrows=F), distance =c(1.5, 1.5, 5)),
          col.regions = rainbow(200),
          xlab=list(expression(theta[1]), cex=2.5), 
          ylab=list(expression(theta[5]), cex=2.5),
          zlab=list(""),
          main=list("MSE of RLR relative to Var(SMC)",cex=2))
dev.off()

# ===== Average of ?LR and SMC at each point
BLR.two <- rowMeans((0.5 * Est.samp.mean + 0.5 * Est.TTC.mean[,1,] - Est$True)^2)/Est$samp.var
SLR.two <- rowMeans((0.5 * Est.samp.mean + 0.5 * Est.TTC.mean[,2,] - Est$True)^2)/Est$samp.var
RLR.two <- rowMeans((0.5 * Est.samp.mean + 0.5 * Est.TTC.mean[,3,] - Est$True)^2)/Est$samp.var
Two.MSE <- data.frame(theta1=Est$theta1, theta5=Est$theta5, BLR=BLR.two, SLR=SLR.two, RLR=RLR.two)
Two.MSE$BLR[Two.MSE$BLR>10] <- 10

jpeg("SANGlobalBLRTwo.jpeg", width = 500, height = 500, units = "px")
wireframe(BLR~theta1*theta5, data=Two.MSE, drape=TRUE, aspect = c(1, 1),
          screen = list(z=-60, y =0, x = -60), colorkey = list(labels=list(cex=2)),
          scales = list(arrows=FALSE,cex=2, tick.number = 4, z = list(arrows=F), distance =c(1.5, 1.5, 4)),
          col.regions = rainbow(200),
          xlab=list(expression(theta[1]), cex=2.5),
          ylab=list(expression(theta[5]), cex=2.5),
          zlab=list(""),
          main=list("MSE of SMC + BLR relative to Var(SMC)",cex=2))
dev.off()

jpeg("SANGlobalSLRTwo.jpeg", width = 500, height = 500, units = "px")
wireframe(SLR~theta1*theta5, data=Two.MSE, drape=TRUE, aspect = c(1, 1),
          screen = list(z=-60, y =0, x = -60), colorkey = list(labels=list(cex=2)),
          scales = list(arrows=FALSE,cex=2, tick.number = 4, z = list(arrows=F), distance =c(1.5, 1.5, 4)),
          col.regions = rainbow(200),
          xlab=list(expression(theta[1]), cex=2.5),
          ylab=list(expression(theta[5]), cex=2.5),
          zlab=list(""),
          main=list("MSE of SMC + SLR relative to Var(SMC)",cex=2))
dev.off()

jpeg("SANGlobalRLRTwo.jpeg", width = 500, height = 500, units = "px")
wireframe(RLR~theta1*theta5, data=Two.MSE, drape=TRUE, aspect = c(1, 1),
          screen = list(z=-60, y =0, x = -60), colorkey = list(labels=list(cex=2)),
          scales = list(arrows=FALSE,cex=2, tick.number = 4, z = list(arrows=F), distance =c(1.5, 1.5, 4)),
          col.regions = rainbow(200),
          xlab=list(expression(theta[1]), cex=2.5),
          ylab=list(expression(theta[5]), cex=2.5),
          zlab=list(""),
          main=list("MSE of SMC + RLR relative to Var(SMC)",cex=2))
dev.off()



Est$LR.mean <- rowMeans(Est.LR.mean)
Est$LR.var <- rowMeans(Est.LR.var)
