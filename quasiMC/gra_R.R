 library("ggpubr")
 my_data<-read.delim("L3.txt",header = FALSE,sep="",dec=".")
 pdf(file='L3_WUS synthesis rate.pdf')
 ggscatter(my_data, x = "V1", y = "V6", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "WUS synthesis rate", ylab = "nuclear WUS in L3")
 dev.off()
 pdf(file='L3_WUS diffusion rate.pdf')
 ggscatter(my_data, x = "V2", y = "V7", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "WUS diffusion rate", ylab = "nuclear WUS in L3")
 dev.off()
 pdf(file='L3_nuclear WUS degradation rate.pdf')
 ggscatter(my_data, x = "V3", y = "V8", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "nuclear WUS degradation rate", ylab = "nuclear WUS in L3")
 dev.off()
 pdf(file='L3_cytoplasm WUS degradation rate.pdf')
 ggscatter(my_data, x = "V4", y = "V9", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "cytoplasmic WUS degradation rate", ylab = "nuclear WUS in L3")
 dev.off()
 pdf(file='L3_nuclear export rate.pdf')
 ggscatter(my_data, x = "V5", y = "V10", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "nuclear export rate", ylab = "nuclear WUS in L3")
 dev.off()
 
 
 
 library("ggpubr")
 my_data<-read.delim("L2.txt",header = FALSE,sep="",dec=".")
 pdf(file='L2_WUS synthesis rate.pdf')
 ggscatter(my_data, x = "V1", y = "V6", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "WUS synthesis rate", ylab = "nuclear WUS in L2")
 dev.off()
 pdf(file='L2_WUS diffusion rate.pdf')
 ggscatter(my_data, x = "V2", y = "V7", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "WUS diffusion rate", ylab = "nuclear WUS in L2")
 dev.off()
 pdf(file='L2_nuclear WUS degradation rate.pdf')
 ggscatter(my_data, x = "V3", y = "V8", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "nuclear WUS degradation rate", ylab = "nuclear WUS in L2")
 dev.off()
 pdf(file='L2_cytoplasm WUS degradation rate.pdf')
 ggscatter(my_data, x = "V4", y = "V9", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "cytoplasmic WUS degradation rate", ylab = "nuclear WUS in L2")
 dev.off()
 pdf(file='L2_nuclear export rate.pdf')
 ggscatter(my_data, x = "V5", y = "V10", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "nuclear export rate", ylab = "nuclear WUS in L2")
 dev.off()
 
 
 library("ggpubr")
 my_data<-read.delim("L1.txt",header = FALSE,sep="",dec=".")
 pdf(file='L1_WUS synthesis rate.pdf')
 ggscatter(my_data, x = "V1", y = "V6", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "WUS synthesis rate", ylab = "nuclear WUS in L1")
 dev.off()
 pdf(file='L1_WUS diffusion rate.pdf')
 ggscatter(my_data, x = "V2", y = "V7", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "WUS diffusion rate", ylab = "nuclear WUS in L1")
 dev.off()
 pdf(file='L1_nuclear WUS degradation rate.pdf')
 ggscatter(my_data, x = "V3", y = "V8", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "nuclear WUS degradation rate", ylab = "nuclear WUS in L1")
 dev.off()
 pdf(file='L1_cytoplasm WUS degradation rate.pdf')
 ggscatter(my_data, x = "V4", y = "V9", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "cytoplasmic WUS degradation rate", ylab = "nuclear WUS in L1")
 dev.off()
 pdf(file='L1_nuclear export rate.pdf')
 ggscatter(my_data, x = "V5", y = "V10", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "nuclear export rate", ylab = "nuclear WUS in L1")
 dev.off()