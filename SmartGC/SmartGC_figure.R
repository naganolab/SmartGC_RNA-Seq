### library----
library("tidyverse")
library("cowplot")

# color
color5 <- c("#ff4b00","#4dc4ff","#03af7a","#fff100","#005aff")
color3 <- color5[1:3]

### execution date
#exec.date <- str_sub(Sys.time(),3,10) %>%
#  str_split(., pattern="-") %>% 
#  unlist() %>% 
#  str_c(.,collapse = "")
exec.date <- "210703"

### load file ----
load("input/Temperature_comparison")
load("input/Humidity_comparison")
SmartGC_260_log <- read.csv("input/SmartGC_260_log.csv")
GC_260_log <- read.csv("input/GC_260_log.csv")
load("input/comparison_all")
load("input/comparison_light_log")

# Fig. 1a ----
# Field Light
time <- as.numeric(1:841)
rownames(comparison_all) <- c(1:4321)
value <- comparison_all[3481:4321,3]
condition <- factor(rep("Field",841))

d <- data.frame(condition, time, value)

g1 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3[1])+
  ylim(0,2500)+
  labs(x="", y="" )+
  theme_void(base_size = 7)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = 'none',
        legend.justification = "center")

print(g1)

# SmartGC Light
time <- as.numeric(1:841)
rownames(comparison_all) <- c(1:4321)
value <- comparison_all[3481:4321,5]
condition <- factor(rep("Field",841))

d <- data.frame(condition, time, value)

g2 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3[2])+
  ylim(0,2500)+
  labs(x="", y="" )+
  theme_void(base_size = 7)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = 'none',
        legend.justification = "center")

print(g2)

# GC Light
time <- as.numeric(1:841)
rownames(comparison_all) <- c(1:4321)
value <- comparison_all[3481:4321,4]
condition <- factor(rep("Field",841))

d <- data.frame(condition, time, value)

g3 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3[3])+
  ylim(0,2500)+
  labs(x="", y="" )+
  theme_void(base_size = 7)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = 'none',
        legend.justification = "center")


print(g3)

# temperature Field
time <- as.numeric(1:1441)
rownames(comparison_all) <- c(1:4321)
value <- Temperature_comparison[2881:4321,2]
condition <- factor(rep("Field",1441))

d <- data.frame(condition, time, value)

g4 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3[1])+
  ylim(15,30)+
  labs(x="", y="" )+
  theme_void(base_size = 7)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = 'none',
        legend.justification = "center")

print(g4)

# temperature SmartGC
time <- as.numeric(1:1441)
rownames(comparison_all) <- c(1:4321)
value <- Temperature_comparison[2881:4321,3]
condition <- factor(rep("Field",1441))

d <- data.frame(condition, time, value)

g5 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3[2])+
  ylim(15,30)+
  labs(x="", y="" )+
  theme_void(base_size = 22)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = 'none',
        legend.justification = "center")

print(g5)

# temperature GC
time <- as.numeric(1:1441)
rownames(comparison_all) <- c(1:4321)
value <- Temperature_comparison[2881:4321,4]
condition <- factor(rep("Field",1441))

d <- data.frame(condition, time, value)

g6 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3[3])+
  ylim(15,30)+
  labs(x="", y="" )+
  theme_void(base_size = 7)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = 'none',
        legend.justification = "center")

print(g6)

# plot
g <- plot_grid(g1,g2,g3,g2,g3, g4,g5,g6,g6,g5, nrow = 2, scale = 0.7)

ggsave(g, file=sprintf("%s_Fig1a.pdf", exec.date),
       width = 55, height = 20, units="mm")

# Fig. 1a2 ----
time <- rep(as.numeric(1:1441),3)
value <- as.matrix(comparison_all[1:1441,c(3,5,4)])
value <- as.vector(value)
condition <- factor(c(rep("FIELD",1441),rep("FL/FTH",1441),rep("CL/CTH",1441)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))

d <- data.frame(condition, time, value)

g1 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  #geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(
    breaks=seq(1,1441,120),
    labels=c("19","21","23","1","3","5","7","9","11","13","15","17","19"))+
  ylim(0,2500)+
  #ggtitle("Irradiance")+
  labs(x="", y=expression(paste("(Î¼mol/m"^2, "/s)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=6, colour = "black"),
        legend.position = 'none')

ggsave(g1, file=sprintf("%s_Fig1a2.pdf", exec.date),
       width = 30, height = 30, units="mm")

# Fig. 1b ----
# GC condition
pre_light <- c(rep(0, 600), rep(867, 840),rep(NA,4321))
pre_temp <- c(rep(22, 600), rep(27, 840),rep(NA,4321))
pre_hum <- c(rep(65, 1440),rep(NA,4321))

#light
time <- rep(as.numeric(1:5761),4)
value <- as.matrix(comparison_all[,c(3,5,4)])
value <- as.vector(value)
value <- c(pre_light, rep(NA,1440),value[1:4321],rep(NA,1440),value[4322:8642],rep(NA,1440),value[8643:12963])

condition <- factor(c(rep("Pre",5761), rep("FIELD",5761),rep("FL/FTH",5761),rep("CL/CTH",5761)),
                    levels = c("Pre","FIELD", "FL/FTH","CL/CTH"))

d <- data.frame(condition, time, value)

g1 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = c("black", color3))+
  geom_vline(xintercept=c(1, 1441, 2881, 4321,5761), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,5761,1440),
                     labels=c("9/17\n19:00","9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  ylim(0,2500)+
  ggtitle("Irradiance")+
  labs(x="", y=expression(paste("(ƒÊmol/m"^2, "/s)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = 'top',
        legend.justification = "center")

print(g1)

### Temperature
value <- c(pre_temp, rep(NA, 1440), Temperature_comparison[,2], rep(NA, 1440), Temperature_comparison[,3],
           rep(NA, 1440), Temperature_comparison[,4])
condition <- factor(c(rep("Pre",5761), rep("FIELD",5761),rep("FL/FTH",5761),rep("CL/CTH",5761)),
                    levels = c("Pre","FIELD", "FL/FTH","CL/CTH"))
time <- rep(1:5761,4)
d <- data.frame(condition, time, value)

g2 <- ggplot(d, aes(x = time, y = value, group = condition, color=condition))+
  geom_line(size=0.25)+
  geom_vline(xintercept=c(1, 1441, 2881, 4321, 5761), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,5761,1440),
                     labels=c("9/17\n19:00","9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  scale_color_manual(values = c("black", color3))+
  ylim(15,30)+
  ggtitle("Temperature")+
  labs(x = "", y = expression(paste("("^o, "C)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

print(g2)

### Humidity
value <- c(pre_hum, rep(NA, 1440), Humidity_comparison[,2], rep(NA, 1440), Humidity_comparison[,3],
           rep(NA, 1440), Humidity_comparison[,4]) 
condition <- factor(c(rep("Pre",5761), rep("FIELD",5761),rep("FL/FTH",5761),rep("CL/CTH",5761)),
                    levels = c("Pre","FIELD", "FL/FTH","CL/CTH"))
time <- rep(1:5761,4)
d <- data.frame(condition, time, value)

g3 <- ggplot(d, aes(x = time, y = value, group = condition, color=condition))+
  geom_line(size=0.25)+
  geom_vline(xintercept=c(1, 1441, 2881, 4321, 5761), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,5761,1440),
                     labels=c("9/17\n19:00","9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  scale_color_manual(values = c("black",color3))+
  ylim(40,100)+
  ggtitle("Relative Humidity")+
  labs(x = "", y = "(%)")+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

print(g3)

mylegend <- get_legend(g1)
g1 <- g1 + theme(legend.position="none")

g <- plot_grid(mylegend, plot_grid(g1,g2,g3, ncol=1, align = "h"),
               ncol=1, rel_heights = c(1,14))

ggsave(g, file=sprintf("%s_Fig1b.pdf", exec.date),
       width = 45, height = 80, units="mm")

# Fig. 1c ----
#light
value <- c(comparison_light_log[,3],comparison_light_log[,4],comparison_light_log[,5])
condition <- factor(c(rep("FIELD",4321),rep("FL/FTH",4321),rep("CL/CTH",4321)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))
time <- rep(1:4321,3)
d <- data.frame(condition, time, value)

g1 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,4321,1440),
                     labels=c("9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  ylim(0,2500)+
  ggtitle("Irradiance")+
  labs(x="", y=expression(paste("(ƒÊmol/m"^2, "/s)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = 'top',
        legend.justification = "center")

print(g1)

### Temperature
value <- c(Temperature_comparison[,2], SmartGC_260_log$temp_log,GC_260_log$temp_log)
condition <- factor(c(rep("FIELD",4321),rep("FL/FTH",4321),rep("CL/CTH",4321)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))
time <- rep(1:4321,3)
d <- data.frame(condition, time, value)

g2 <- ggplot(d, aes(x = time, y = value, group = condition, color=condition))+
  geom_line(size=0.25)+
  geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,4321,1440),
                     labels=c("9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  scale_color_manual(values = color3)+
  ylim(15,30)+
  ggtitle("Temperature")+
  labs(x = "", y = expression(paste("("^o, "C)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

print(g2)

### Humidity
value <- c(Humidity_comparison[,2], SmartGC_260_log$hum_log,GC_260_log$hum_log)
condition <- factor(c(rep("FIELD",4321),rep("FL/FTH",4321),rep("CL/CTH",4321)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))
time <- rep(1:4321,3)
d <- data.frame(condition, time, value)

g3 <- ggplot(d, aes(x = time, y = value, group = condition, color=condition))+
  geom_line(size=0.25)+
  geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,4321,1440),
                     labels=c("9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  scale_color_manual(values = color3)+
  ylim(40,100)+
  ggtitle("Relative Humidity")+
  labs(x = "", y = "(%)")+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

print(g3)

mylegend <- get_legend(g1)
g1 <- g1 + theme(legend.position="none")

g <- plot_grid(mylegend, plot_grid(g1,g2,g3, ncol=1, align = "h"),
               ncol=1, rel_heights = c(1,14))

ggsave(g, file=sprintf("%s_Fig1c.pdf", exec.date),
       width = 45, height = 80, units="mm")


# Fig. S24 ----
# load standard data
standard <- read.csv("input/standard_left.csv")

# calculation
s400 <- lm(standard$X400~standard$setting)
s420 <- lm(standard$X420~standard$setting)
s450 <- lm(standard$X450~standard$setting)
s530 <- lm(standard$X530~standard$setting)
s630 <- lm(standard$X630~standard$setting)
s660 <- lm(standard$X660~standard$setting)
s735 <- lm(standard$X735~standard$setting)

a <- s400$coefficients[2]+
  s420$coefficients[2]+
  s450$coefficients[2]+
  s530$coefficients[2]+
  s630$coefficients[2]+
  s660$coefficients[2]+
  s735$coefficients[2]

b <- s400$coefficients[1]+
  s420$coefficients[1]+
  s450$coefficients[1]+
  s530$coefficients[1]+
  s630$coefficients[1]+
  s660$coefficients[1]+
  s735$coefficients[1]

a1 <- a*15+b
a2 <- a*1000+b

set <- c(rep(standard$setting,7),15,1000)
value <- c(as.matrix(standard[,2:8]) %>% as.vector(),a1,a2)
set2 <- c(rep(c("400","420","450","530","630","660","735"),each=12),rep("Sum",2))
d <- data.frame(set, value, set2)

source("script/wavelength_to_rgb.R")
color_rgb <- c(wavelength_to_rgb(400),wavelength_to_rgb(420),wavelength_to_rgb(450),
               wavelength_to_rgb(530),wavelength_to_rgb(630),wavelength_to_rgb(660),
               wavelength_to_rgb(735),"black")

#plot
g1 <- ggplot(d, aes(x=set, y=value, color = set2))+
  geom_smooth(size=0.25, method = "lm", se = FALSE, aes(color = set2))+
  geom_point(size=0.25)+
  scale_color_manual(values = color_rgb)+
  xlim(0,1000)+
  ylim(-100,2000)+
  ggtitle("Calibration curve")+
  labs(x="Output of the light source", y=expression(paste("Irradiance (ƒÊmol/m"^2, "/s)")))+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = 'top',
        legend.justification = "center")+
  guides(color=guide_legend(nrow=2, byrow = T))


print(g1)

ggsave(g1, file=sprintf("%s_FigS24.pdf", exec.date),
       width = 80, height = 80, units="mm")

### Fig. S2 ----
# 170921 5:00-7:00
time <- rep(as.numeric(1:121),3)
rownames(comparison_all) <- c(1:4321)
value <- c(comparison_all[3481:3601,3],comparison_all[3481:3601,5],comparison_all[3481:3601,4])
condition <- factor(c(rep("FIELD",121),rep("FL/FTH",121),rep("CL/CTH",121)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))

d <- data.frame(condition, time, value)

g1 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=seq(1,121,30),
    labels=c("5:00","5:30","6:00","6:30","7:00"))+
  ylim(0,1000)+
  ggtitle("Irradiance")+
  labs(x="", y=expression(paste("(ƒÊmol/m"^2, "/s)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = 'top',
        legend.justification = "center")

# 170921 5:30-6:30
time <- c(rep(as.numeric(1:61),2),c(1:13, 13+ 100/995.026824))
rownames(comparison_all) <- c(1:4321)
value <- c(comparison_all[3511:3571,3],comparison_all[3511:3571,5],rep(0,13),100)
condition <- factor(c(rep("FIELD",61),rep("FL/FTH",61),rep("CL/CTH",14)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))

d <- data.frame(condition, time, value)

gg1 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=c(1,31,61),
    labels=c("5:30","6:00","6:30"))+
  ylim(0,100)+
  labs(x="", y=expression(paste("(ƒÊmol/m"^2, "/s)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

# temperature 5:00-7:00
value <- c(Temperature_comparison[3481:3601,2],Temperature_comparison[3481:3601,3],
           Temperature_comparison[3481:3601,4])
condition <- factor(c(rep("FIELD",121),rep("FL/FTH",121),rep("CL/CTH",121)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))
time <- rep(1:121,3)
d <- data.frame(condition, time, value)

g2 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=seq(1,121,30),
    labels=c("5:00","5:30","6:00","6:30","7:00"))+
  ylim(15,25)+
  ggtitle("Temperature")+
  labs(x="", y = expression(paste("("^o, "C)")) )+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')
print(g2)

# Humidity 5:00-7:00
value <- c(Humidity_comparison[3481:3601,2],Humidity_comparison[3481:3601,3],
           Humidity_comparison[3481:3601,4])
condition <- factor(c(rep("FIELD",121),rep("FL/FTH",121),rep("CL/CTH",121)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))
time <- rep(1:121,3)
d <- data.frame(condition, time, value)

g3 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=seq(1,121,30),
    labels=c("5:00","5:30","6:00","6:30","7:00"))+
  ylim(50,100)+
  ggtitle("Relative humidity")+
  labs(x="Time (hours)", y="(%)" )+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')
print(g3)

### evening
# 170921 17:00-19:00
time <- rep(as.numeric(1:121),3)
rownames(comparison_all) <- c(1:4321)
value <- c(comparison_all[4201:4321,3],comparison_all[4201:4321,5],comparison_all[4201:4321,4])
condition <- factor(c(rep("FIELD",121),rep("FL/FTH",121),rep("CL/CTH",121)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))

d <- data.frame(condition, time, value)

g4 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=seq(1,121,30),
    labels=c("17:00","17:30","18:00","18:30","19:00"))+
  ylim(0,1000)+
  ggtitle("Irradiance")+
  labs(x="", y=expression(paste("(ƒÊmol/m"^2, "/s)")) )+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')
print(g4)

# 170921 17:30-18:30
time <- c(rep(as.numeric(1:61),2),26+ 895.026824/995.026824, 27:61)
rownames(comparison_all) <- c(1:4321)
value <- c(comparison_all[4231:4291,3],comparison_all[4231:4291,5], 100, rep(0,35))
condition <- factor(c(rep("FIELD",61),rep("FL/FTH",61),rep("CL/CTH",36)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))

d <- data.frame(condition, time, value)

gg2 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=c(1,31,61),
    labels=c("17:30","18:00","18:30"))+
  ylim(0,100)+
  labs(x="", y=expression(paste("(ƒÊmol/m"^2, "/s)")) )+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

print(gg2)


# temperature 17:00-19:00
value <- c(Temperature_comparison[4201:4321,2],Temperature_comparison[4201:4321,3],
           Temperature_comparison[4201:4321,4])
condition <- factor(c(rep("FIELD",121),rep("FL/FTH",121),rep("CL/CTH",121)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))
time <- rep(1:121,3)
d <- data.frame(condition, time, value)

g5 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=seq(1,121,30),
    labels=c("17:00","17:30","18:00","18:30","19:00"))+
  ylim(15,25)+
  ggtitle("Temperature")+
  labs(x="Time(hour)", y = expression(paste("("^o, "C)")) )+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')
print(g5)

# Humidity 17:00-19:00
value <- c(Humidity_comparison[4201:4321,2],Humidity_comparison[4201:4321,3],
           Humidity_comparison[4201:4321,4])
condition <- factor(c(rep("FIELD",121),rep("FL/FTH",121),rep("CL/CTH",121)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))
time <- rep(1:121,3)
d <- data.frame(condition, time, value)

g6 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=seq(1,121,30),
    labels=c("17:00","17:30","18:00","18:30","19:00"))+
  ylim(50,100)+
  ggtitle("Relative humidity")+
  labs(x="Time(hour)", y="(%)" )+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')
print(g6)


mylegend <- get_legend(g1)
g1 <- g1 + theme(legend.position = "none")

g <- plot_grid(mylegend, plot_grid(g1,g4,g2,g5,g3,g6, ncol=2, align = "hv"),
               ncol=1, rel_heights = c(1,9))

ggsave(g, file=sprintf("%s_FigS2_1.pdf", exec.date),
       width = 120, height = 120, units = "mm")

p <- plot_grid(gg1,gg2, ncol=2, align = "hv")

ggsave(p, file=sprintf("%s_FigS2_2.pdf", exec.date),
       width = 60, height = 30, units = "mm")

### Fig. S3 ----
### 170920 19:00- 170921 19:00
time <- rep(as.numeric(1:1441),3)
rownames(comparison_all) <- c(1:4321)
value <- c(comparison_all[2881:4321,3],comparison_all[2881:4321,5],comparison_all[2881:4321,4])
condition <- factor(c(rep("FIELD",1441),rep("FL/FTH",1441),rep("CL/CTH",1441)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))

d <- data.frame(condition, time, value)

g1 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=seq(1,1441,120),
    labels=c("19","21","23","1","3","5","7","9","11","13","15","17","19"))+
  ylim(0,2500)+
  ggtitle("Irradiance")+
  labs(x="", y=expression(paste("(ƒÊmol/m"^2, "/s)")))+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = 'top',
        legend.justification = "center")


print(g1)

# temperature
time <- rep(as.numeric(1:1441),3)
rownames(Temperature_comparison) <- c(1:4321)
value <- c(Temperature_comparison[2881:4321,2],Temperature_comparison[2881:4321,3],
           Temperature_comparison[2881:4321,4])
condition <- factor(c(rep("FIELD",1441),rep("FL/FTH",1441),rep("CL/CTH",1441)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))

d <- data.frame(condition, time, value)

g2 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=seq(1,1441,120),
    labels=c("19","21","23","1","3","5","7","9","11","13","15","17","19"))+
  ylim(15,30)+
  ggtitle("Temperature")+
  labs(x="Time(hour)", y = expression(paste("("^o, "C)")) )+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7, colour = "black"),
        legend.position = 'none')

print(g2)

# Humidity
time <- rep(as.numeric(1:1441),3)
rownames(Humidity_comparison) <- c(1:4321)
value <- c(Humidity_comparison[2881:4321,2],Humidity_comparison[2881:4321,3],
           Humidity_comparison[2881:4321,4])
condition <- factor(c(rep("FIELD",1441),rep("FL/FTH",1441),rep("CL/CTH",1441)),
                    levels = c("FIELD", "FL/FTH","CL/CTH"))

d <- data.frame(condition, time, value)

g3 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25)+
  scale_color_manual(values = color3)+
  scale_x_continuous(
    breaks=seq(1,1441,120),
    labels=c("19","21","23","1","3","5","7","9","11","13","15","17","19"))+
  ylim(40,100)+
  ggtitle("Relative humidity")+
  labs(x="Time (hours)", y="(%)" )+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7, colour = "black"),
        legend.position = 'none')

print(g3)

mylegend <- get_legend(g1)
g1 <- g1 + theme(legend.position="none")

g <- plot_grid(mylegend, plot_grid(g1,g2,g3, ncol=1, align = "hv"),
               ncol=1, rel_heights = c(1,9))

ggsave(g, file=sprintf("%s_FigS3.pdf", exec.date),
       width = 90, height = 120, units="mm")


### Fig. S4 ----
# irradiance FL/FTH
value <- c(comparison_all[,5], comparison_light_log[,4])
condition <- factor(c(rep("Setting",4321),rep("Log",4321)), levels = c("Setting","Log"))
time <- rep(1:4321,2)
d <- data.frame(condition, time, value)

g1 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25, alpha=0.5)+
  scale_color_manual(values = c("black",color3[2]))+
  geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,4321,1440),
                     labels=c("9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  ylim(0,2000)+
  ggtitle("Irradiance")+
  labs(x="", y=expression(paste("(ƒÊmol/m"^2, "/s)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = 'top',
        legend.justification = "center")

print(g1)

# irradiance CL/CTH
value <- c(comparison_all[,4], comparison_light_log[,5])
condition <- factor(c(rep("Setting",4321),rep("Log",4321)), levels = c("Setting","Log"))
time <- rep(1:4321,2)
d <- data.frame(condition, time, value)

g2 <- ggplot(d, aes(x=time, y=value, group = condition, color=condition))+
  geom_line(size=0.25, alpha=0.5)+
  scale_color_manual(values = c("black",color3[3]))+
  geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,4321,1440),
                     labels=c("9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  ylim(0,2000)+
  ggtitle("Irradiance")+
  labs(x="", y=expression(paste("(ƒÊmol/m"^2, "/s)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = 'top',
        legend.justification = "center")

print(g2)

### Temperature FL/FTH
value <- c(Temperature_comparison[,3], SmartGC_260_log$temp_log)
condition <- factor(c(rep("Setting",4321),rep("Log",4321)), levels = c("Setting","Log"))
time <- rep(1:4321,2)
d <- data.frame(condition, time, value)

g3 <- ggplot(d, aes(x = time, y = value, group = condition, color=condition))+
  geom_line(size=0.25, alpha=0.5)+
  geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,4321,1440),
                     labels=c("9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  scale_color_manual(values = c("black",color3[2]))+
  ylim(15,30)+
  ggtitle("Temperature")+
  labs(x = "", y = expression(paste("("^o, "C)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

print(g3)

dif <- Temperature_comparison[,3] - SmartGC_260_log$temp_log

range(dif)


### Temperature CL/CTH
value <- c(Temperature_comparison[,4], GC_260_log$temp_log)
condition <- factor(c(rep("Setting",4321),rep("Log",4321)), levels = c("Setting","Log"))
time <- rep(1:4321,2)
d <- data.frame(condition, time, value)

g4 <- ggplot(d, aes(x = time, y = value, group = condition, color=condition))+
  geom_line(size=0.25, alpha=0.5)+
  geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,4321,1440),
                     labels=c("9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  scale_color_manual(values = c("black",color3[3]))+
  ylim(15,30)+
  ggtitle("Temperature")+
  labs(x = "", y = expression(paste("("^o, "C)")))+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

print(g4)

### Humidity FL/FTH
value <- c(Humidity_comparison[,3], SmartGC_260_log$hum_log)
condition <- factor(c(rep("Setting",4321),rep("Log",4321)), levels = c("Setting","Log"))
time <- rep(1:4321,2)
d <- data.frame(condition, time, value)

g5 <- ggplot(d, aes(x = time, y = value, group = condition, color=condition))+
  geom_line(size=0.25, alpha=0.5)+
  geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,4321,1440),
                     labels=c("9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  scale_color_manual(values=c("black",color3[2]))+
  ylim(40,100)+
  ggtitle("Relative Humidity")+
  labs(x = "", y = "(%)")+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

print(g5)

dif <- Humidity_comparison[,3] -  SmartGC_260_log$hum_log
range(dif)


### Humidity CL/CTH
value <- c(Humidity_comparison[,4], GC_260_log$hum_log)
condition <- factor(c(rep("Setting",4321),rep("Log",4321)), levels = c("Setting","Log"))
time <- rep(1:4321,2)
d <- data.frame(condition, time, value)

g6 <- ggplot(d, aes(x = time, y = value, group = condition, color=condition))+
  geom_line(size=0.25, alpha=0.5)+
  geom_vline(xintercept=c(1, 1441, 2881, 4321), colour="black", linetype=2, size=0.25)+
  scale_x_continuous(breaks=seq(1,4321,1440),
                     labels=c("9/18\n19:00","9/19\n19:00","9/20\n19:00","9/21\n19:00"))+
  scale_color_manual(values=c("black",color3[3]))+
  ylim(40,100)+
  ggtitle("Relative Humidity")+
  labs(x = "", y = "(%)")+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')

print(g6)


mylegend <- get_legend(g1)
g1 <- g1 + theme(legend.position="none")

p1 <- plot_grid(mylegend, plot_grid(g1,g3,g5, ncol=1, align = "hv"),
               ncol=1, rel_heights = c(1,14))

mylegend <- get_legend(g2)
g2 <- g2 + theme(legend.position="none")

p2 <- plot_grid(mylegend, plot_grid(g2,g4,g6, ncol=1, align = "hv"),
               ncol=1, rel_heights = c(1,14))

p <- plot_grid(p1,p2)

ggsave(p, file=sprintf("%s_FigS4.pdf", exec.date),
       width = 120, height = 120, units="mm")
