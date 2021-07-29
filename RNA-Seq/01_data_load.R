################
### function ###
################

library("fields")
library('tidyverse')
library("ggsci")
library("ggrepel")
library("RColorBrewer")
library("car")
library("multcomp")
library("VennDiagram")
library("cowplot")
library("gridExtra")
library("colorhcplot")
library("gplots")
library("TCC")
library("ggbeeswarm")
library("dendextend")
library("coin")
library("ggsignif")
library("lemon")
source("script/cb.R")
source("script/ng.Colors.R")

### GO analysis 
# load the GO analysis functions
source("GO/GOanalysis_170821/GOanalysis_functions.R")

# load the table of Gene ID and GO (as ulg)
load("GO/GOanalysis_170821/ulg.RAP_190219") # for rice

### KEGG analysis 
# load the GO analysis functions
source("GO/KEGGanalysis_200107/200107_KEGG_function.R")

# load kegg_rice
load("GO/KEGGanalysis_200107/kegg_rice") 
load("GO/KEGGanalysis_200107/kegg_description") 

### color
color5 <- c("#ff4b00","#4dc4ff","#03af7a","#fff100","#005aff")
color4 <- color5[2:5]
color3 <- color5[1:3]
color_5 <- c("#00ffff","#ff00ff","#00ff00","#ff8000","#0080ff")

### execution date
#exec.date <- str_sub(Sys.time(),3,10) %>%
#  str_split(., pattern="-") %>% 
#  unlist() %>% 
#  str_c(.,collapse = "")
exec.date <- "210703"

####################
### Experiment_1 ###
####################

# set parameters ----
dir.input <- "input"
dir.output <- "output"
#dir.create(dir.output)

# load sample attribute --------------
at <- read.delim("input/OsaYH_2-5_SampleAttribute.txt", header=T)

# load Osa transcript description 
load(file="input/180604_des_simple_organelle")

# set data read row 
drr <- des[,"NormalizationGroup"]=="data"
names(drr) <- rownames(des)
des2 <- des[drr,]

# load expression data table-------------------------------------
load(file="input/190521-2_rawcnt_Osativa_1")
load(file="input/190521-2_rawcnt_Osativa_2")
rawcnt <- cbind(rawcnt_1, rawcnt_2)
save(rawcnt, file = "input/190521-2_rawcnt_Osativa")

load(file="input/190521-3_rpm_Osativa_1")
load(file="input/190521-3_rpm_Osativa_2")
rpm <- cbind(rpm_1, rpm_2)
save(rpm, file = "input/190521-3_rpm_Osativa")

load(file="input/190521-4_cvrd1_Osativa")
load(file="input/190521-5_cvrd3_Osativa")

# label ----
condition <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- factor(unique(at$hour), levels=unique(at$hour))
hour <- factor(at$hour, levels = time)

# gene criteria ----
raw_10 <- rowSums(rawcnt[drr,])/ncol(rawcnt[drr,])>10
rpm2 <- rpm[drr,]
log2rpm <- log2(rpm2+1)

### mean and sd (rpm) ----
d <- as_tibble(t(rpm2)) %>% 
  mutate(condition, hour) %>% 
  gather(gene,expression,-condition, -hour) 

d_mean_sd_condition <- d %>% 
  group_by(condition, hour, gene) %>% 
  summarize(mean = mean(expression, na.rm=T), sd=sd(expression, na.rm=T)) %>% 
  unite(., condition, hour, col="condition", sep="_")

d_mean_sd_condition$condition <- factor(d_mean_sd_condition$condition,
                                        levels = unique(d_mean_sd_condition$condition))

rpm.ave <- d_mean_sd_condition %>% 
  dplyr::select(., -sd) %>% 
  spread(., key=condition, value = mean) %>% 
  as.data.frame()

rownames(rpm.ave) <- rpm.ave[,1]
rpm.ave <- rpm.ave[,2:66]

rpm.sd <- d_mean_sd_condition %>% 
  dplyr::select(., -mean) %>% 
  spread(., key=condition, value = sd) %>% 
  as.data.frame()

rownames(rpm.sd) <- rpm.sd[,1]
rpm.sd <- rpm.sd[,2:66]

### mean and sd (log2rpm) ----
d <- as_tibble(t(log2rpm)) %>% 
  mutate(condition, hour) %>% 
  gather(gene,expression,-condition, -hour) 

d_mean_sd_condition <-   d %>% 
  group_by(condition, hour, gene) %>% 
  summarize(mean = mean(expression, na.rm=T), sd=sd(expression, na.rm=T)) %>% 
  unite(., condition, hour, col="condition", sep="_")

d_mean_sd_condition$condition <- factor(d_mean_sd_condition$condition,
                                        levels = unique(d_mean_sd_condition$condition))

log2rpm.ave <- d_mean_sd_condition %>% 
  dplyr::select(., -sd) %>% 
  spread(., key=condition, value = mean) %>% 
  as.data.frame()

rownames(log2rpm.ave) <- log2rpm.ave[,1]
log2rpm.ave <- log2rpm.ave[,2:66]

log2rpm.sd <- d_mean_sd_condition %>% 
  dplyr::select(., -mean) %>% 
  spread(., key=condition, value = sd) %>% 
  as.data.frame()

rownames(log2rpm.sd) <- log2rpm.sd[,1]
log2rpm.sd <- log2rpm.sd[,2:66]

# save data
save(rpm2, file = "input/210703_rpm2")
save(log2rpm, file = "input/210703_log2rpm")
save(rpm.ave, file = "input/210703_rpmave")
save(log2rpm.ave, file = "input/210703_log2rpmave")

####################
### Experiment_2 ###
####################

### load data ###
# load sample attribute --------------
at_204 <- read.delim("input/OsaYH_6_7_SampleAttribute.txt", header=T, as.is=T)

# load expression data table-------------------------------------
load(file="input/190515-2_rawcnt_Osativa_1")
load(file="input/190515-2_rawcnt_Osativa_2")
rawcnt <- cbind(rawcnt_1, rawcnt_2)
save(rawcnt, file = "input/190515-2_rawcnt_Osativa")

load(file="input/190515-3_rpm_Osativa_1")
load(file="input/190515-3_rpm_Osativa_2")
rpm <- cbind(rpm_1, rpm_2)
save(rpm, file = "input/190515-3_rpm_Osativa")

rpm2 <- rpm[drr,]
log2rpm <- log2(rpm2+1)

# label ----
condition_204 <- factor(at_204$condition, levels = c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time_204 <- factor(unique(at_204$hour),levels =unique(at_204$hour))
hour_204 <- factor(at_204$hour, levels = time_204)

### mean and sd (rpm) ----
condition <- condition_204
time <- time_204
hour <- hour_204

d <- as_tibble(t(rpm2)) %>% 
  mutate(condition, hour) %>% 
  gather(gene,expression,-condition, -hour) 

d_mean_sd_condition <- d %>% 
  group_by(condition, hour, gene) %>% 
  summarize(mean = mean(expression, na.rm=T), sd=sd(expression, na.rm=T)) %>% 
  unite(., condition, hour, col="condition", sep="_")

d_mean_sd_condition$condition <- factor(d_mean_sd_condition$condition,
                                        levels = unique(d_mean_sd_condition$condition))

rpm.ave <- d_mean_sd_condition %>% 
  dplyr::select(., -sd) %>% 
  spread(., key=condition, value = mean) %>% 
  as.data.frame()

rownames(rpm.ave) <- rpm.ave[,1]
rpm.ave <- rpm.ave[,2:69]

rpm.sd <- d_mean_sd_condition %>% 
  dplyr::select(., -mean) %>% 
  spread(., key=condition, value = sd) %>% 
  as.data.frame()

rownames(rpm.sd) <- rpm.sd[,1]
rpm.sd <- rpm.sd[,2:69]

### mean and sd (log2rpm) ----
d <- as_tibble(t(log2rpm)) %>% 
  mutate(condition, hour) %>% 
  gather(gene,expression,-condition, -hour) 

d_mean_sd_condition <-   d %>% 
  group_by(condition, hour, gene) %>% 
  summarize(mean = mean(expression, na.rm=T), sd=sd(expression, na.rm=T)) %>% 
  unite(., condition, hour, col="condition", sep="_")

d_mean_sd_condition$condition <- factor(d_mean_sd_condition$condition,
                                        levels = unique(d_mean_sd_condition$condition))

log2rpm.ave <- d_mean_sd_condition %>% 
  dplyr::select(., -sd) %>% 
  spread(., key=condition, value = mean) %>% 
  as.data.frame()

rownames(log2rpm.ave) <- log2rpm.ave[,1]
log2rpm.ave <- log2rpm.ave[,2:69]

log2rpm.sd <- d_mean_sd_condition %>% 
  dplyr::select(., -mean) %>% 
  spread(., key=condition, value = sd) %>% 
  as.data.frame()

rownames(log2rpm.sd) <- log2rpm.sd[,1]
log2rpm.sd <- log2rpm.sd[,2:69]

### name ----
rawcnt_204 <- rawcnt
rpm_204 <- rpm
rpm2_204 <- rpm2
log2rpm_204 <- log2rpm
rpm.ave_204 <- rpm.ave
log2rpm.ave_204 <- log2rpm.ave

### load data in Experiment 1 ----
load(file="input/190521-2_rawcnt_Osativa")
load(file="input/190521-3_rpm_Osativa")
load(file="input/190521-4_cvrd1_Osativa")
load(file="input/190521-5_cvrd3_Osativa")
load(file="input/210703_rpm2")
load(file="input/210703_log2rpm")
load(file="input/210703_rpmave")
load(file="input/210703_log2rpmave")

### Fig. S5a ----
a <- rawcnt[drr,]
b <- apply(a, 2, sum)

c <- rawcnt_204[drr,]
d <- apply(c, 2, sum)

exp <- factor(c(rep("Exp.1", 260),rep("Exp.2",204)))
reads <- log10(c(b,d))
d <- data.frame(exp, reads)

g <- ggplot(d, aes(x = reads, fill=exp))+
  geom_histogram(color="black", size=0.25)+
  scale_fill_manual(values = c("white","gray"))+
  xlim(5,7)+
  labs(x="Total read count (log10)", y="Frequency")+
  theme_cowplot(font_size = 7, line_size=0.25)+
  theme(plot.title = element_blank(),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "top",
        legend.justification = "center")
p <- plot_grid(g)

ggsave(p, file=sprintf("%s/%s_FigS5a.pdf",dir.output, exec.date),
       width = 60, height = 40, units="mm")

### Fig. S5b ----
a <- rawcnt[drr,]
e <- apply(a, 1, mean)
e <- e[e>0]

d <- data.frame(exp=log10(e))

g <- ggplot(d, aes(x = exp))+
  geom_histogram(bins=50, color="black", size=0.25)+
  geom_vline(xintercept = 1, color=color5[1], size =0.5)+
  scale_fill_manual(values = "gray")+
  xlim(-5,5)+
  ylim(0,3000)+
  labs(x="Mean read count (log10)", y="Frequency")+
  theme_cowplot(font_size = 7, line_size=0.25)+
  theme(plot.title = element_blank(),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "top",
        legend.justification = "center")
p <- plot_grid(g)

ggsave(p, file=sprintf("%s/%s_FigS5b.pdf",dir.output, exec.date),
       width = 60, height = 40, units="mm")
