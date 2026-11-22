# ===== Workspace preparation
rm(list = ls()) # clear all variables
set.seed(1) # fix random seed for debugging and reproduction of results

load("SAN_10K_Macro.RData")

library(reshape2)
library(ggplot2)
data <- as.data.frame(MSE)
data$id <- 1:nrow(data)

data.all <- melt(data,id.var="id")
ggplot(data.all, aes(x = id,y = value, colour = variable, lty = variable)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    scale_x_continuous(name = "Number of Iterations", trans = "log10",
                       breaks = c(1,10,100),
                       labels = c("1", "10", "100")) +
    scale_y_continuous(name ="Mean Squared Error", trans = "log10")+
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1),
          legend.position = c(0,0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(colour ="black"),
          legend.key.width = unit(2, "cm")) +
    labs(colour="Estimators",linetype="Estimators")

data.select <- data.all[data.all$variable %in% c("Crude","Brute", "Avg", "SM.RLR"), ]
ggplot(data.select, aes(x = id,y = value, colour = variable, lty = variable, shape=variable)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    scale_x_continuous(name = "Number of Iterations", trans = "log10",
                       breaks = c(1,10,100),
                       labels = c("1", "10", "100")) +
    scale_y_continuous(trans = "log10")+
    ggtitle("Mean Squared Errors")+
    theme(plot.title = element_text(size = 16,face="bold",hjust=0.5),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1),
          legend.position = c(0,0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(colour ="black"),
          legend.key.width = unit(2, "cm")) +
    labs(colour="Estimators",linetype="Estimators",shape="Estimators")

ggsave("SANEst_CrudeBruteSMRLR.jpeg", device = "jpeg", width = 8, height = 6, units = "in")
ggsave("SANEst_CrudeBruteSMRLR.eps", device = "eps", width = 8, height = 6, units = "in")

##########################################################################
data.BLR <- data.all[data.all$variable %in% c("PW.BLR","OW.BLR", "SM.BLR"), ]
ggplot(data.BLR, aes(x = id,y = value, colour = variable, lty = variable, shape=variable)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    scale_x_continuous(name = "Number of Iterations", trans = "log10",
                       breaks = c(1,10,100),
                       labels = c("1", "10", "100")) +
    scale_y_continuous(trans = "log10")+
    ggtitle("Mean Squared Errors")+
    theme(plot.title = element_text(size = 16,face="bold",hjust=0.5),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1),
          legend.position = c(0,0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(colour ="black"),
          legend.key.width = unit(2, "cm")) +
    labs(colour="Estimators",linetype="Estimators",shape="Estimators")

ggsave("SANEst_BLR.jpeg", device = "jpeg", width = 8, height = 6, units = "in")
ggsave("SANEst_BLR.eps", device = "eps", width = 8, height = 6, units = "in")


data.SLR <- data.all[data.all$variable %in% c("PW.SLR","OW.SLR", "SM.SLR"), ]
ggplot(data.SLR, aes(x = id,y = value, colour = variable, lty = variable, shape=variable)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    scale_x_continuous(name = "Number of Iterations", trans = "log10",
                       breaks = c(1,10,100),
                       labels = c("1", "10", "100")) +
    scale_y_continuous(trans = "log10")+
    ggtitle("Mean Squared Errors")+
    theme(plot.title = element_text(size = 16,face="bold",hjust=0.5),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1),
          legend.position = c(0,0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(colour ="black"),
          legend.key.width = unit(2, "cm")) +
    labs(colour="Estimators",linetype="Estimators",shape="Estimators")

ggsave("SANEst_SLR.jpeg", device = "jpeg", width = 8, height = 6, units = "in")
ggsave("SANEst_SLR.eps", device = "eps", width = 8, height = 6, units = "in")

data.RLR <- data.all[data.all$variable %in% c("PW.RLR","OW.RLR", "SM.RLR","Crude","Brute"), ]
ggplot(data.RLR, aes(x = id,y = value, colour = variable, lty = variable, shape=variable)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    scale_x_continuous(name = "Number of Iterations", trans = "log10",
                       breaks = c(1,10,100),
                       labels = c("1", "10", "100")) +
    scale_y_continuous(trans = "log10")+
    ggtitle("Mean Squared Errors")+
    theme(plot.title = element_text(size = 16,face="bold",hjust=0.5),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1),
          legend.position = c(0,0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(colour ="black"),
          legend.key.width = unit(2, "cm")) +
    labs(colour="Estimators",linetype="Estimators",shape="Estimators")+ 
    scale_shape_discrete(labels=c("Crude", "Brute", "RLR-PW", "RLR-OW", "RLR-M")) + 
    scale_colour_discrete(labels=c("Crude", "Brute", "RLR-PW", "RLR-OW", "RLR-M")) + 
    scale_linetype_discrete(labels=c("Crude", "Brute", "RLR-PW", "RLR-OW", "RLR-M"))
ggsave("SANEst_RLR.jpeg", device = "jpeg", width = 6, height = 5, units = "in")
ggsave("SANEst_RLR.eps", device = "eps", width = 6, height = 5, units = "in")

##########################################################################
data.PW <- data.all[data.all$variable %in% c("PW.BLR","PW.SLR", "PW.RLR","Crude","Brute"), ]
ggplot(data.PW, aes(x = id,y = value, colour = variable, lty = variable, shape=variable)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    scale_x_continuous(name = "Number of Iterations", trans = "log10",
                       breaks = c(1,10,100),
                       labels = c("1", "10", "100")) +
    scale_y_continuous(trans = "log10", limits = c(0.002,0.5))+
    ggtitle("Mean Squared Errors")+
    theme(plot.title = element_text(size = 16,face="bold",hjust=0.5),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1),
          legend.position = c(0,0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(colour ="black"),
          legend.key.width = unit(2, "cm")) +
    labs(colour="Estimators",linetype="Estimators",shape="Estimators")+ 
    scale_shape_discrete(labels=c("Crude", "Brute", "BLR-PW", "SLR-PW", "RLR-PW")) + 
    scale_colour_discrete(labels=c("Crude", "Brute", "BLR-PW", "SLR-PW", "RLR-PW")) +
    scale_linetype_discrete(labels=c("Crude", "Brute", "BLR-PW", "SLR-PW", "RLR-PW"))
ggsave("SANEst_PW.jpeg", device = "jpeg", width = 6, height = 5, units = "in")
ggsave("SANEst_PW.eps", device = "eps", width = 6, height = 5, units = "in")


data.OW <- data.all[data.all$variable %in% c("OW.BLR","OW.SLR", "OW.RLR","Crude","Brute"), ]
ggplot(data.OW, aes(x = id,y = value, colour = variable, lty = variable, shape=variable)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    scale_x_continuous(name = "Number of Iterations", trans = "log10",
                       breaks = c(1,10,100),
                       labels = c("1", "10", "100")) +
    scale_y_continuous(trans = "log10", limits = c(0.002,0.5))+
    ggtitle("Mean Squared Errors")+
    theme(plot.title = element_text(size = 16,face="bold",hjust=0.5),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1),
          legend.position = c(0,0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(colour ="black"),
          legend.key.width = unit(2, "cm")) +
    labs(colour="Estimators",linetype="Estimators",shape="Estimators")+ 
    scale_shape_discrete(labels=c("Crude", "Brute", "BLR-OW", "SLR-OW", "RLR-OW")) + 
    scale_colour_discrete(labels=c("Crude", "Brute", "BLR-OW", "SLR-OW", "RLR-OW")) +
    scale_linetype_discrete(labels=c("Crude", "Brute", "BLR-OW", "SLR-OW", "RLR-OW"))
ggsave("SANEst_OW.jpeg", device = "jpeg", width = 6, height = 5, units = "in")
ggsave("SANEst_OW.eps", device = "eps", width = 6, height = 5, units = "in")

data.SM <- data.all[data.all$variable %in% c("SM.BLR","SM.SLR", "SM.RLR","Crude","Brute"), ]
ggplot(data.SM, aes(x = id,y = value, colour = variable, lty = variable, shape=variable)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    scale_x_continuous(name = "Number of Iterations", trans = "log10",
                       breaks = c(1,10,100),
                       labels = c("1", "10", "100")) +
    scale_y_continuous(trans = "log10", limits = c(0.002,0.5))+
    ggtitle("Mean Squared Errors")+
    theme(plot.title = element_text(size = 16,face="bold",hjust=0.5),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1),
          legend.position = c(0,0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(colour ="black"),
          legend.key.width = unit(2, "cm")) +
    labs(colour="Estimators",linetype="Estimators",shape="Estimators") + 
    scale_shape_discrete(labels=c("Crude", "Brute", "BLR-M", "SLR-M", "RLR-M")) + 
    scale_colour_discrete(labels=c("Crude", "Brute", "BLR-M", "SLR-M", "RLR-M")) +
    scale_linetype_discrete(labels=c("Crude", "Brute", "BLR-M", "SLR-M", "RLR-M"))
ggsave("SANEst_SM.jpeg", device = "jpeg", width = 6, height = 5, units = "in")
ggsave("SANEst_SM.eps", device = "eps", width = 6, height = 5, units = "in")
