rm(list =ls ())

setwd("/Users/lix233/Desktop/JRSS-C/simulation_res_plot")


boxdat = rbind(read.csv('simu_2_s10.csv', header = T), 
               read.csv('simu_2_s12.csv', header = T),
               read.csv('simu_2_s15.csv', header = T),
               read.csv('simu_2_s16.csv', header = T),
               read.csv('simu_2_s17.csv', header = T))



library(plyr)
library(reshape2)
library(lattice)
library(ggplot2)


png('/Users/lix233/Desktop/JRSS-C/boxplot_error_rate_newerr.png',width=1500,height=800)
sp <- ggplot(boxdat, aes(x=method, y=error.rate, fill=method)) + 
  geom_boxplot(coef=4, outlier.shape = NA, position = "dodge", alpha = 0.6, 
               lwd = 1, fatten = 0.75) +
  ylim(0, 1.0) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1,size=20),
        axis.text.y = element_text(size=20),
        axis.ticks.length=unit(.25, "cm"),
        legend.text= element_text(size=20),
        legend.key.size =  unit(0.5, "in"),
        axis.title.y = element_text(size = 30, angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0))) + 
  labs(x = "",y = "Mis-assignment rate") +
  theme(strip.text.x = element_text(size=30, colour = "black", face = "bold"))

sp + facet_wrap(~s, ncol=6, scales="free")

dev.off()



png('/Users/lix233/Desktop/JRSS-C/boxplot_pearson_chi_sq_newerr.png',width=1500,height=800)
sp <- ggplot(boxdat, aes(x=method, y=pearson.chi.sq, fill=method)) + 
  geom_boxplot(coef=4, outlier.shape = NA, position = "dodge", alpha = 0.6, 
               lwd = 1, fatten = 0.75) +
  ylim(0, 400.0) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1,size=20),
        axis.text.y = element_text(size=20),
        axis.ticks.length=unit(.25, "cm"),
        legend.text= element_text(size=20),
        legend.key.size =  unit(0.5, "in"),
        axis.title.y = element_text(size = 30, angle = 90,margin = margin(t = 0, r = 20, b = 0, l = 0))) + 
  labs(x = "", y = "Pearson's chi-square") +
  theme(strip.text.x = element_text(size=30, colour = "black", face = "bold"))

sp + facet_wrap(~s, ncol=6, scales="free")

dev.off()







