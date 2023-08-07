# Correlate Jackson2020 scRNAseq dataset with YTK GFP data (streamlined version)
# Using overall sds and means

require(dplyr)
require(tidyverse)
require(xlsx)
require(ggplot2)
require(reshape2)
require(diptest)
require(viridis)

# Set working directory to current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Functions
getDipPval <- function(df) {
  d.t <- dip.test(df)
  return(d.t$p.value)
}

# Load YTK GFP master data
GFP.data = read.xlsx(file = file.path("../data","GFP_overall_stats.xlsx") ,1, header=TRUE)

# Load scRNAseq counts
scRNAseq.data=readRDS(file = file.path("../data","Jackson2020_lognorm_WT.Rds"))

# Get genes included in YTK GFP data
selected.genes = GFP.data$gene


#* Pre-process data --------------------------------------------------------

# 1. Calculate metrics for all scRNAseq entries and store in scRNAseq.master

# counts
scRNAseq.YPD = scRNAseq.data$YPD
scRNAseq.MMD = scRNAseq.data$MinimalGlucose

# sums
scRNAseq.YPD.sum = apply(scRNAseq.data$YPD[,c(-2:-1)],1,sum) 
scRNAseq.MMD.sum = apply(scRNAseq.data$MinimalGlucose[,c(-2:-1)],1,sum) 

# SDs
#scRNAseq.YPD.sd = apply(scRNAseq.data$YPD[,c(-2:-1)],1,sd)
#scRNAseq.MMD.sd = apply(scRNAseq.data$MinimalGlucose[,c(-2:-1)],1,sd)

# CVs
scRNAseq.YPD.CV = apply(scRNAseq.data$YPD[,c(-2:-1)],1,sd)/apply(scRNAseq.data$YPD[,c(-2:-1)],1,mean)
scRNAseq.MMD.CV = apply(scRNAseq.data$MinimalGlucose[,c(-2:-1)],1,sd)/apply(scRNAseq.data$MinimalGlucose[,c(-2:-1)],1,mean)

# means
scRNAseq.YPD.mean = apply(scRNAseq.data$YPD[,c(-2:-1)],1,mean)
scRNAseq.MMD.mean = apply(scRNAseq.data$MinimalGlucose[,c(-2:-1)],1,mean)

# dip tests
scRNAseq.YPD.dipt = apply(scRNAseq.data$YPD[,c(-2:-1)],1,getDipPval)
scRNAseq.MMD.dipt = apply(scRNAseq.data$MinimalGlucose[,c(-2:-1)],1,getDipPval)

# merge metrics into scRNAseq.master
scRNAseq.master=data.frame(scRNAseq.YPD$SystematicName,scRNAseq.YPD$Name,scRNAseq.YPD.sum,scRNAseq.MMD.sum,scRNAseq.YPD.CV,scRNAseq.MMD.CV,scRNAseq.YPD.mean,scRNAseq.MMD.mean,scRNAseq.YPD.dipt,scRNAseq.MMD.dipt) %>%
  dplyr::rename(SystematicName = scRNAseq.YPD.SystematicName, Name = scRNAseq.YPD.Name)

#saveRDS(scRNAseq.master,file.path("/Users/kshab/Documents/ICL/YTK_heterogeneity","results/Robject/scRNAseq_master.Rds"))


# 2. Master table joining GFP.master and scRNAseq.master for selected genes

selected.genes = c('TDH3','CCW12','PGK1','HHF2','TEF1','TEF2','HHF1','HTB2','RPL18B','ALD6','PAB1','RET2','RNR1','SAC6','RNR2','POP6','RAD27','PSP2','REV1')

master <- scRNAseq.master %>% filter(Name %in% selected.genes) %>%
  inner_join(GFP.data,by = c('Name' = 'gene'))




#* Plot correlations grouped by conditions ---------------------------------
#   Melt and regroup by condition ####

# MEAN scRNAseq
master.by.mean.scRNAseq <- master %>% select(SystematicName,
                                         Name,
                                         scRNAseq.YPD.mean,
                                         scRNAseq.YPD.CV,
                                         scRNAseq.MMD.mean,
                                         scRNAseq.MMD.CV) %>%
  melt(id= c("SystematicName","Name","scRNAseq.YPD.CV","scRNAseq.MMD.CV")) %>%
  dplyr::rename("scRNAseq.mean"="value") %>%
  mutate(
    condition = if_else(grepl("YPD", variable),"YPD","YNB"),
  )

# CV scRNAseq
master.by.CV.scRNAseq <- master %>% select(SystematicName,
                                   Name,
                                   scRNAseq.YPD.mean,
                                   scRNAseq.YPD.CV,
                                   scRNAseq.MMD.mean,
                                   scRNAseq.MMD.CV) %>%
  melt(id= c("SystematicName","Name","scRNAseq.YPD.mean","scRNAseq.MMD.mean")) %>%
  dplyr::rename("scRNAseq.CV"="value") %>%
  mutate(
    condition = if_else(grepl("YPD", variable),"YPD","YNB"),
  )

# MEAN GFP
master.by.mean.GFP <- master %>% select(SystematicName,
                                   Name,
                                   YNB.mu,
                                   YNB.CV,
                                   YPD.mu,
                                   YPD.CV) %>%
  melt(id= c("SystematicName","Name","YNB.CV","YPD.CV")) %>%
  dplyr::rename("GFP.mean"="value") %>%
  mutate(
    condition = if_else(grepl("YPD", variable),"YPD","YNB"),
  )

# CV GFP
master.by.CV.GFP <- master %>% select(SystematicName,
                                       Name,
                                       YNB.mu,
                                       YNB.CV,
                                       YPD.mu,
                                       YPD.CV) %>%
  melt(id= c("SystematicName","Name","YNB.mu","YPD.mu")) %>%
  dplyr::rename("GFP.CV"="value") %>%
  mutate(
    condition = if_else(grepl("YPD", variable),"YPD","YNB"),
  )



# Merge all
master.by.condition <- master.by.mean.scRNAseq %>% inner_join(master.by.mean.GFP, by = c("SystematicName"="SystematicName",
                                                           "Name"="Name",
                                                           "condition"="condition")) %>%
  inner_join(master.by.CV.scRNAseq, by = c("SystematicName"="SystematicName",
                                  "Name"="Name",
                                  "condition"="condition")) %>%
  inner_join(master.by.CV.GFP, by = c("SystematicName"="SystematicName",
                                           "Name"="Name",
                                           "condition"="condition")) %>%
  select(Name,SystematicName,condition,GFP.mean,GFP.CV,scRNAseq.mean,scRNAseq.CV)
  
#   Plot correlation between GFP and scRNAseq ####

colour_map=c('#d4818e','#b5c6cb') # YNB/YPD

# Mean plot
ggplot(master.by.condition, aes(x=scRNAseq.mean, y=GFP.mean, color=condition))+
  geom_point(size=2)+
  scale_colour_manual(values=colour_map)+
  xlab('scRNAseq reads \n mean reads per gene (log scaled)')+
  ylab('GFP intensities \n mean intensity per gene (log scaled)')+
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(file.path("../figures","Mean_YPD_YNB_overall.svg"), width = 4, height = 2.5, bg='transparent')


# CV plot
ggplot(master.by.condition, aes(x=scRNAseq.CV, y=GFP.CV, color=condition))+
  geom_point(size=2)+
  scale_colour_manual(values=colour_map)+
  xlab('scRNAseq reads \n CV')+
  ylab('GFP intensities \n CV')+
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(file.path("../figures","CV_YPD_YNB_overall.svg"), width = 4, height = 2.5, bg='transparent')


#* Compute regression ------------------------------------------------------

require(GauPro)

norm_minmax <- function(x){
  (x- min(x)) /(max(x)-min(x))
}

#   Obtain testing datasets and reformating ######

GFP.test = read.xlsx(file = file.path("../data","GFP_overall_stats.xlsx") ,2, header=TRUE)

tested.genes = c('HCS1','HHT1','TDH2','HHT2','CDC19','RPL28')

master.test <- scRNAseq.master %>% filter(Name %in% tested.genes) %>%
  inner_join(GFP.test,by = c('Name' = 'gene'))

master.test.by.sum.scRNAseq <- master.test %>% select(SystematicName,
                                            Name,
                                            scRNAseq.YPD.sum,
                                            scRNAseq.YPD.CV,
                                            scRNAseq.MMD.sum,
                                            scRNAseq.MMD.CV) %>%
  melt(id= c("SystematicName","Name","scRNAseq.YPD.CV","scRNAseq.MMD.CV")) %>%
  dplyr::rename("scRNAseq.sum"="value") %>%
  mutate(
    condition = if_else(grepl("YPD", variable),"YPD","YNB"),
  )

# MEAN scRNAseq
master.test.by.mean.scRNAseq <- master.test %>% select(SystematicName,
                                             Name,
                                             scRNAseq.YPD.mean,
                                             scRNAseq.YPD.CV,
                                             scRNAseq.MMD.mean,
                                             scRNAseq.MMD.CV) %>%
  melt(id= c("SystematicName","Name","scRNAseq.YPD.CV","scRNAseq.MMD.CV")) %>%
  dplyr::rename("scRNAseq.mean"="value") %>%
  mutate(
    condition = if_else(grepl("YPD", variable),"YPD","YNB"),
  )

# CV scRNAseq
master.test.by.CV.scRNAseq <- master.test %>% select(SystematicName,
                                           Name,
                                           scRNAseq.YPD.mean,
                                           scRNAseq.YPD.CV,
                                           scRNAseq.MMD.mean,
                                           scRNAseq.MMD.CV) %>%
  melt(id= c("SystematicName","Name","scRNAseq.YPD.mean","scRNAseq.MMD.mean")) %>%
  dplyr::rename("scRNAseq.CV"="value") %>%
  mutate(
    condition = if_else(grepl("YPD", variable),"YPD","YNB"),
  )

# MEAN GFP
master.test.by.mean.GFP <- master.test %>% select(SystematicName,
                                        Name,
                                        YNB.mu,
                                        YNB.CV,
                                        YPD.mu,
                                        YPD.CV) %>%
  melt(id= c("SystematicName","Name","YNB.CV","YPD.CV")) %>%
  dplyr::rename("GFP.mean"="value") %>%
  mutate(
    condition = if_else(grepl("YPD", variable),"YPD","YNB"),
  )

# CV GFP
master.test.by.CV.GFP <- master.test %>% select(SystematicName,
                                      Name,
                                      YNB.mu,
                                      YNB.CV,
                                      YPD.mu,
                                      YPD.CV) %>%
  melt(id= c("SystematicName","Name","YNB.mu","YPD.mu")) %>%
  dplyr::rename("GFP.CV"="value") %>%
  mutate(
    condition = if_else(grepl("YPD", variable),"YPD","YNB"),
  )



# Merge all
master.test.by.condition <- master.test.by.mean.scRNAseq %>% inner_join(master.test.by.mean.GFP, by = c("SystematicName"="SystematicName",
                                                                                         "Name"="Name",
                                                                                         "condition"="condition")) %>%
  inner_join(master.test.by.CV.scRNAseq, by = c("SystematicName"="SystematicName",
                                           "Name"="Name",
                                           "condition"="condition")) %>%
  inner_join(master.test.by.CV.GFP, by = c("SystematicName"="SystematicName",
                                      "Name"="Name",
                                      "condition"="condition")) %>%
  select(Name,SystematicName,condition,GFP.mean,GFP.CV,scRNAseq.mean,scRNAseq.CV)



#   Plot regression with testing data added ######

# Mean YPD

x<-master.by.condition %>% filter(condition=='YPD') %>% 
  select(scRNAseq.mean) %>% pull(scRNAseq.mean)
y<-master.by.condition %>% filter(condition=='YPD') %>% 
  select(GFP.mean) %>% pull(GFP.mean)

x.test<-master.test.by.condition %>% filter(condition=='YPD') %>% 
  select(scRNAseq.mean) %>% pull(scRNAseq.mean)
y.test<-master.test.by.condition %>% filter(condition=='YPD') %>% 
  select(GFP.mean) %>% pull(GFP.mean)

# Merge training and testing for common normalization
x.all.norm <- c(x,x.test) %>% norm_minmax() %>% as.numeric()
y.all.norm <- c(y,y.test) %>% norm_minmax() %>% as.numeric()

training <- c(rep(1,length(x)),rep(0,length(x.test)))

data <- cbind(x.all.norm,y.all.norm,training) %>% as.data.frame()

x.norm <- data %>% filter(training==1) %>% pull(x.all.norm) 
y.norm <- data %>% filter(training==1) %>% pull(y.all.norm)

# Perform regression on normalized training data
gp <- GauPro(x.norm, y.norm, parallel=FALSE, kernel=kern, gpk=gpk)

colour_map= c('#000000','#b5c6cb')
ggplot(data, aes(x=x.all.norm, y=y.all.norm, colour=as.factor(training)))+
  geom_point(size=2,fill=as.factor(training))+
  stat_function(fun=function(x) gp$predict(x), col='#b5c6cb')+
  stat_function(fun=function(x) gp$predict(x)+1*gp$predict(x, se=T)$se, col='#b5c6cb', linetype="dashed")+
  stat_function(fun=function(x) gp$predict(x)-1*gp$predict(x, se=T)$se, col='#b5c6cb', linetype="dashed")+
  scale_colour_manual(values=colour_map)+
  xlab('Normalized mean scRNAseq reads (a.u.)')+
  ylab('Normalized GFP intensities (a.u.)')+
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(file.path("../figures","Mean_YPD_GP_overall_testing.svg"), width = 4, height = 2.5, bg='transparent')


# CV YPD

# Removing genes which high error bars from YNB dataset which leads to strong overfitting
#genes.to.keep = c('TDH3','PGK1','TEF1','TEF2','RPL18B','ALD6','PAB1','RET2','SAC6','RNR2','POP6','RAD27','PSP2','REV1')

#x<-master.by.condition %>% filter(condition=='YPD') %>% 
#  filter(Name %in% genes.to.keep) %>% pull(scRNAseq.CV)
#y<-master.by.condition %>% filter(condition=='YPD') %>% 
#  filter(Name %in% genes.to.keep) %>% pull(GFP.CV)

x<-master.by.condition %>% filter(condition=='YPD') %>% 
  select(scRNAseq.CV) %>% pull(scRNAseq.CV)
y<-master.by.condition %>% filter(condition=='YPD') %>% 
  select(GFP.CV) %>% pull(GFP.CV)

x.test<-master.test.by.condition %>% filter(condition=='YPD') %>% 
  select(scRNAseq.CV) %>% pull(scRNAseq.CV)
y.test<-master.test.by.condition %>% filter(condition=='YPD') %>% 
  select(GFP.CV) %>% pull(GFP.CV)

# Merge training and testing for common normalization
x.all.norm <- c(x,x.test) %>% norm_minmax() %>% as.numeric()
y.all.norm <- c(y,y.test) %>% norm_minmax() %>% as.numeric()

training <- c(rep(1,length(x)),rep(0,length(x.test)))

data <- cbind(x.all.norm,y.all.norm,training) %>% as.data.frame()

x.norm <- data %>% filter(training==1) %>% pull(x.all.norm) 
y.norm <- data %>% filter(training==1) %>% pull(y.all.norm)

# Perform regression on normalized training data
gp <- GauPro(x.norm, y.norm, parallel=FALSE, kernel=kern, gpk=gpk)

colour_map= c('#000000','#b5c6cb')
ggplot(data, aes(x=x.all.norm, y=y.all.norm, colour=as.factor(training)))+
  geom_point(size=2,fill=as.factor(training))+
  stat_function(fun=function(x) gp$predict(x), col='#b5c6cb')+
  stat_function(fun=function(x) gp$predict(x)+1*gp$predict(x, se=T)$se, col='#b5c6cb', linetype="dashed")+
  stat_function(fun=function(x) gp$predict(x)-1*gp$predict(x, se=T)$se, col='#b5c6cb', linetype="dashed")+
  scale_colour_manual(values=colour_map)+
  xlab('Normalized CV scRNAseq reads (a.u.)')+
  ylab('Normalized GFP CV (a.u.)')+
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(file.path("../figures","CV_YPD_GP_overall_testing.svg"), width = 4, height = 2.5, bg='transparent')


# Mean YNB

x<-master.by.condition %>% filter(condition=='YNB') %>% 
  select(scRNAseq.mean) %>% pull(scRNAseq.mean)
y<-master.by.condition %>% filter(condition=='YNB') %>% 
  select(GFP.mean) %>% pull(GFP.mean)

x.test<-master.test.by.condition %>% filter(condition=='YNB') %>% 
  select(scRNAseq.mean) %>% pull(scRNAseq.mean)
y.test<-master.test.by.condition %>% filter(condition=='YNB') %>% 
  select(GFP.mean) %>% pull(GFP.mean)

# Merge training and testing for common normalization
x.all.norm <- c(x,x.test) %>% norm_minmax() %>% as.numeric()
y.all.norm <- c(y,y.test) %>% norm_minmax() %>% as.numeric()

training <- c(rep(1,length(x)),rep(0,length(x.test)))

data <- cbind(x.all.norm,y.all.norm,training) %>% as.data.frame()

x.norm <- data %>% filter(training==1) %>% pull(x.all.norm) 
y.norm <- data %>% filter(training==1) %>% pull(y.all.norm)

# Perform regression on normalized training data
gp <- GauPro(x.norm, y.norm, parallel=FALSE, kernel=kern, gpk=gpk)

colour_map= c('#000000','#d4818e')
ggplot(data, aes(x=x.all.norm, y=y.all.norm, colour=as.factor(training)))+
  geom_point(size=2,fill=as.factor(training))+
  stat_function(fun=function(x) gp$predict(x), col='#d4818e')+
  stat_function(fun=function(x) gp$predict(x)+1*gp$predict(x, se=T)$se, col='#d4818e', linetype="dashed")+
  stat_function(fun=function(x) gp$predict(x)-1*gp$predict(x, se=T)$se, col='#d4818e', linetype="dashed")+
  scale_colour_manual(values=colour_map)+
  xlab('Normalized mean scRNAseq reads (a.u.)')+
  ylab('Normalized GFP intensities (a.u.)')+
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(file.path("../figures","Mean_YNB_GP_overall_testing.svg"), width = 4, height = 2.5, bg='transparent')


# CV YNB

# Removing genes which high error bars from YNB dataset which leads to strong overfitting
genes.to.keep = c('TDH3','PGK1','TEF1','TEF2','RPL18B','ALD6','PAB1','RET2','SAC6','RNR2','POP6','RAD27','PSP2','REV1')

#x<-master.by.condition %>% filter(condition=='YNB') %>% 
#  filter(Name %in% genes.to.keep) %>% pull(scRNAseq.CV)
#y<-master.by.condition %>% filter(condition=='YNB') %>% 
#  filter(Name %in% genes.to.keep) %>% pull(GFP.CV)

x<-master.by.condition %>% filter(condition=='YNB') %>% 
  pull(scRNAseq.CV)
y<-master.by.condition %>% filter(condition=='YNB') %>% 
  pull(GFP.CV)

x.test<-master.test.by.condition %>% filter(condition=='YNB') %>% 
  select(scRNAseq.CV) %>% pull(scRNAseq.CV)
y.test<-master.test.by.condition %>% filter(condition=='YNB') %>% 
  select(GFP.CV) %>% pull(GFP.CV)

# Merge training and testing for common normalization
x.all.norm <- c(x,x.test) %>% norm_minmax() %>% as.numeric()
y.all.norm <- c(y,y.test) %>% norm_minmax() %>% as.numeric()

training <- c(rep(1,length(x)),rep(0,length(x.test)))

data <- cbind(x.all.norm,y.all.norm,training) %>% as.data.frame()

x.norm <- data %>% filter(training==1) %>% pull(x.all.norm) 
y.norm <- data %>% filter(training==1) %>% pull(y.all.norm)

# Perform regression on normalized training data
gp <- GauPro(x.norm, y.norm, parallel=FALSE, kernel=kern, gpk=gpk)

colour_map= c('#000000','#d4818e')
ggplot(data, aes(x=x.all.norm, y=y.all.norm, colour=as.factor(training)))+
  geom_point(size=2,fill=as.factor(training))+
  stat_function(fun=function(x) gp$predict(x), col='#d4818e')+
  stat_function(fun=function(x) gp$predict(x)+1*gp$predict(x, se=T)$se, col='#d4818e', linetype="dashed")+
  stat_function(fun=function(x) gp$predict(x)-1*gp$predict(x, se=T)$se, col='#d4818e', linetype="dashed")+
  scale_colour_manual(values=colour_map)+
  xlab('Normalized CV scRNAseq reads (a.u.)')+
  ylab('Normalized GFP CV (a.u.)')+
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(file.path("../figures","CV_YNB_GP_overall_testing.svg"), width = 4, height = 2.5, bg='transparent')



# UMAP plotting -----------------------------------------------------------

# Load UMAP parameters and gene counts from scRNAseq data (stored in sce format)
sce.conditions = readRDS(file = file.path("../data","Jackson2020_sce_conditions_UMAP.rds"))
sce.params = readRDS(file = file.path("../data","Jackson2020_sce_params_UMAP.rds"))

# Merge counts and UMAP parameters for all genotypes and selected conditions

counts.YPD <- t(sce.conditions$YPD@assays@data@listData[["logcounts"]]) %>% 
  as.data.frame() %>% rownames_to_column('cellID')
selection.YPD <- sce.params$YPD %>% rownames_to_column('cellID') %>% inner_join(counts.YPD)

counts.YNB <- t(sce.conditions$MinimalGlucose@assays@data@listData[["logcounts"]]) %>% 
  as.data.frame() %>% rownames_to_column('cellID')
selection.YNB <- sce.params$MinimalGlucose%>% rownames_to_column('cellID') %>% inner_join(counts.YNB)

## Set condition, genotype (default is WT) and gene name for UMAP plot 

selection <- selection.YPD %>% filter(genotypeGroup=='WT(ho)')

display.gene <- dplyr::pull(selection, 'YNL301C') # Can visualize other genes here
selection$target <- display.gene

ggplot(selection, aes(x=UMAP1, y=UMAP2, colour=target))+
  geom_point(size=1.5, alpha=0.6)+
  scale_colour_viridis(option='B',limits = c(0, 5))+ #discrete=TRUE/limits = c(0.2, 1)
  #scale_colour_brewer(palette="Set3") +
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme_classic()+
  theme(panel.background =  element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.ticks = element_line(color = "black"))

#ggsave("../figures/YPD_RPL18B.svg", width = 4.2, height = 3.3, bg='transparent')

# Facet wrap for all promoters tested in this study

targets_comName <- c('TDH3','CCW12','PGK1','HHF2','TEF1','TEF2','HHF1','HTB2','RPL18B','ALD6','PAB1','RET2',
                     'RNR1','SAC6','RNR2','POP6','RAD27','PSP2','REV1','HCS1','HHT1',
                     'TDH2','HHT2','CDC19','RPL28')

targets_name <- c('YGR192C','YLR110C','YCR012W','YNL030W','YPR080W','YBR118W','YBR009C','YBL002W','YNL301C','YPL061W',
                  'YER165W','YFR051C','YER070W','YDR129C','YJL026W','YGR030C','YKL113C',
                  'YML017W','YOR346W','YKL017C','YBR010W','YJR009C','YNL031C',
                  'YAL038W','YGL103W')

targets <- data.frame(targets_name,targets_comName)

selection_fw.YPD <- selection.YPD %>% filter(genotypeGroup=='WT(ho)') %>%
  select(one_of(c("cellID","UMAP1","UMAP2",targets_name))) %>% 
  pivot_longer(cols=any_of(targets_name), names_to = "target", values_to = "intensity") %>%
  inner_join(targets, by= c("target"="targets_name"))

selection_fw.YNB <- selection.YNB %>% filter(genotypeGroup=='WT(ho)') %>%
  select(one_of(c("cellID","UMAP1","UMAP2",targets_name))) %>% 
  pivot_longer(cols=any_of(targets_name), names_to = "target", values_to = "intensity") %>%
  inner_join(targets, by= c("target"="targets_name"))

p <- ggplot(selection_fw.YPD, aes(x=UMAP1, y=UMAP2, colour=intensity))+
  geom_point(size=1.5, alpha=0.3)+
  scale_colour_viridis(option='B',limits = c(0, 7.3))+ #discrete=TRUE/limits = c(0.2, 1)
  #scale_colour_brewer(palette="Set3") +
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme_classic()+
  theme(panel.background =  element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.ticks = element_line(color = "black"))

# Use vars() to supply faceting variables:
p + facet_wrap(vars(targets_comName))

ggsave("../figures/YPD_facet.svg", width = 5, height = 5, bg='transparent')

p <- ggplot(selection_fw.YNB, aes(x=UMAP1, y=UMAP2, colour=intensity))+
  geom_point(size=1.5, alpha=0.3)+
  scale_colour_viridis(option='B',limits = c(0, 7.3))+ #discrete=TRUE/limits = c(0.2, 1)
  #scale_colour_brewer(palette="Set3") +
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme_classic()+
  theme(panel.background =  element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.ticks = element_line(color = "black"))

# Use vars() to supply faceting variables:
p + facet_wrap(vars(targets_comName))

ggsave("../figures/YNB_facet.svg", width = 5, height = 5, bg='transparent')


# Histogram reads ---------------------------------------------------------

p <- ggplot(selection_fw.YPD, aes(x=intensity))+
  geom_histogram(position='identity', bins=30, alpha=0.4)+
  theme_classic()+
  theme(panel.background =  element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.ticks = element_line(color = "black"))

# Use vars() to supply faceting variables:
p + facet_wrap(vars(targets_comName))

ggsave("../figures/YPD_hist_facet.svg", width = 5, height = 5, bg='transparent')


p <- ggplot(selection_fw.YNB, aes(x=intensity))+
  geom_histogram(position='identity', bins=30, alpha=0.4)+
  theme_classic()+
  theme(panel.background =  element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.ticks = element_line(color = "black"))

# Use vars() to supply faceting variables:
p + facet_wrap(vars(targets_comName))

ggsave("../figures/YNB_hist_facet.svg", width = 5, height = 5, bg='transparent')
