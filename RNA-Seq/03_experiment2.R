####################
### Experiment_2 ###
####################

condition <- factor(at_204$condition, levels = c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- time_204
hour <- hour_204

##############
### Figure ###
##############

# Fig. 1h ----
pca <- prcomp(t(log2rpm_204[raw_10,]), scale = T)
pc <- pca$x
#PC1 vs PC2
xpc <- pca$sdev[1]^2 / sum(pca$sdev^2) *100
xlab <- sprintf("PC1 (%2.1f%%)",xpc)
ypc <- pca$sdev[2]^2 / sum(pca$sdev^2) *100
ylab <- sprintf("PC2 (%2.1f%%)",ypc)

times <- hour_204
set <- factor(c(rep("FL/FTH",51),rep("CL/CTH",51),rep("FL/CTH",51),rep("CL/FTH",51)),
              levels=c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
d <- data.frame(times,set,pc)

# mean and sd of PC1 and PC2
pc_mean <- data.frame(pc) %>% 
  mutate(set,hour) %>% 
  group_by(set,hour) %>%
  summarize_all(mean, na.rm=T)

pc_sd <- data.frame(pc) %>% 
  mutate(set,hour) %>% 
  group_by(set,hour) %>%
  summarize_all(sd, na.rm=T)

pcs <- data.frame(mean_PC1 = pc_mean$PC1, mean_PC2 = pc_mean$PC2,
                  sd_PC1 = pc_sd$PC1, sd_PC2 = pc_sd$PC2)

times2 <- rep(time, 4)
set2 <- factor(c(rep("FL/FTH",17),rep("CL/CTH",17),rep("FL/CTH",17),rep("CL/FTH",17)),
              levels=c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
d2 <- data.frame(times2,set2,pcs)

# plot
g <- ggplot()+
  geom_point(data=d, aes(x=PC1, y=PC2, group=set, shape=set, fill=set), size=0.5, stroke=0.25, alpha=0.5)+
  geom_point(data = d2, aes(x=mean_PC1, y=mean_PC2, group=set2, shape=set2, fill=set2), size=1, stroke=0.25, )+
  geom_errorbarh(data = d2, height=0.2, show.legend = FALSE,size=0.25,
                 aes(y=mean_PC2, xmin = (mean_PC1 - sd_PC1), xmax = (mean_PC1 + sd_PC1)))+
  geom_errorbar(data = d2, width=0.2, show.legend = FALSE, size=0.25,
                aes(x=mean_PC1, ymin = (mean_PC2 - sd_PC2), ymax = (mean_PC2 + sd_PC2)))+
  geom_path(data=d2, aes(x=mean_PC1, y=mean_PC2, colour=set2,group=set2,),size=0.25, show.legend = FALSE)+
  #geom_text_repel(data=d2, aes(x=mean_PC1, y=mean_PC2, label=times2), size=2, show.legend = FALSE)+
  scale_color_manual(values = color4)+
  scale_fill_manual(values = color4)+
  scale_shape_manual(values = c(22:25))+
  labs(x=xlab, y=ylab)+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6, colour = "black"),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = "top",
        legend.justification = "center")+
  guides(fill=guide_legend(nrow=2,byrow=T))

g

ggsave(g, file=sprintf("%s/%s_Fig1h.pdf",dir.output, exec.date),
       width = 45, height = 45, units = "mm")

# Fig. 1j ----
set <- factor(c(rep("FL/FTH",17),rep("CL/CTH",17),rep("FL/CTH",17),rep("CL/FTH",17)),
              levels = c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
rpm2.ave <- log2rpm.ave_204[raw_10,]
soukan <- cor(rpm2.ave, method = "spearman")
souka <- stats::dist(1-soukan)
souk <- hclust(souka, method = "ward.D2")
souka <- NULL
dend <- as.dendrogram(souk) %>% set("branches_lwd",0.25)
df2 <- data.frame(states = factor(souk$labels, levels = souk$labels[souk$order]))
labels(dend) <- rep(NA, 68)
times <- rep(time,4)

p1 <- ggplot(as.ggdend(dend))+
  scale_x_continuous(expand = c(0, 0.5), breaks = 1:68, labels = NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(plot.margin= unit(c(1, 0, -1, 0), "lines"))

# Tiles and labels
p2 <- ggplot(df2,aes(states, y = 1, fill = times)) +
  geom_tile(color="black", size = 0.25) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)  

maxwidth <- unit.pmax(gp1$widths, gp2$widths)
gp1$widths <- as.list(maxwidth)
gp2$widths <- as.list(maxwidth)
a <- grid.arrange(gp1, gp2, ncol = 1, heights = c(3, 1))

ggsave(a, file= sprintf("%s/%s_Fig1j.pdf", dir.output, exec.date),
       width = 100, height = 60, units = "mm")

dend2 <- as.dendrogram(souk) %>% set("branches_lwd",0.25) %>% set("labels_cex",0.4)
p3 <- ggplot(as.dendrogram(dend2))+
  scale_x_discrete(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 3))+
  theme(plot.margin= unit(c(1, 1, 1, 1), "lines"))

ggsave(p3, file= sprintf("%s/%s_Fig1j2.pdf", dir.output, exec.date),
       width = 200, height = 100, units = "mm")

# Fig. S12 ----
pca <- prcomp(t(log2rpm_204[raw_10,]), scale = T)
pc <- pca$x
#PC1 vs PC2
xpc <- pca$sdev[1]^2 / sum(pca$sdev^2) *100
xlab <- sprintf("PC1 (%2.1f%%)",xpc)
ypc <- pca$sdev[2]^2 / sum(pca$sdev^2) *100
ylab <- sprintf("PC2 (%2.1f%%)",ypc)

times <- hour_204
set <- factor(c(rep("FL/FTH",51),rep("CL/CTH",51),rep("FL/CTH",51),rep("CL/FTH",51)),
              levels=c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
d <- data.frame(times,set,pc)

# mean and sd of PC1 and PC2
pc_mean <- data.frame(pc) %>% 
  mutate(set,hour) %>% 
  group_by(set,hour) %>%
  summarize_all(mean, na.rm=T)

pc_sd <- data.frame(pc) %>% 
  mutate(set,hour) %>% 
  group_by(set,hour) %>%
  summarize_all(sd, na.rm=T)

pcs <- data.frame(mean_PC1 = pc_mean$PC1, mean_PC2 = pc_mean$PC2,
                  sd_PC1 = pc_sd$PC1, sd_PC2 = pc_sd$PC2)

times2 <- rep(time, 4)
set2 <- factor(c(rep("FL/FTH",17),rep("CL/CTH",17),rep("FL/CTH",17),rep("CL/FTH",17)),
               levels=c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
d2 <- data.frame(times2,set2,pcs)

# plot
g <- ggplot()+
  #geom_point(data=d, aes(x=PC1, y=PC2, group=set, shape=set, fill=set), size=1, stroke=0.25, alpha=0.5)+
  geom_point(data = d2, aes(x=mean_PC1, y=mean_PC2, group=set2, shape=set2, fill=set2), size=2, stroke=0.25)+
  geom_errorbarh(data = d2, height=2, show.legend = FALSE,size=0.25, alpha=0.5,
                 aes(y=mean_PC2, xmin = (mean_PC1 - sd_PC1), xmax = (mean_PC1 + sd_PC1)))+
  geom_errorbar(data = d2, width=2, show.legend = FALSE, size=0.25, alpha=0.5,
                aes(x=mean_PC1, ymin = (mean_PC2 - sd_PC2), ymax = (mean_PC2 + sd_PC2)))+
  #geom_path(data=d2, aes(x=mean_PC1, y=mean_PC2, colour=set2,group=set2,),size=0.5, show.legend = FALSE)+
  geom_text_repel(data=d2, aes(x=mean_PC1, y=mean_PC2, label=times2),
                  segment.color = "#84919e", segment.size = 0.25, size=2, show.legend = FALSE)+
  scale_color_manual(values = color4)+
  scale_fill_manual(values = color4)+
  scale_shape_manual(values = c(22:25))+
  labs(x=xlab, y=ylab)+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7, colour = "black"),
        axis.text.y = element_text(size=7, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "top",
        legend.justification = "center")

g

ggsave(g, file=sprintf("%s/%s_FigS12a.pdf",dir.output, exec.date),
       width = 90, height = 90, units = "mm")

# plot
g <- ggplot()+
  geom_point(data=d, aes(x=PC1, y=PC2, group=set, shape=set, fill=set), size=1, stroke=0.25, alpha=1)+
  #geom_point(data = d2, aes(x=mean_PC1, y=mean_PC2, group=set2, shape=set2, fill=set2), size=2, stroke=0.25)+
  #geom_errorbarh(data = d2, height=2, show.legend = FALSE,size=0.25, alpha=0.5,
  #               aes(y=mean_PC2, xmin = (mean_PC1 - sd_PC1), xmax = (mean_PC1 + sd_PC1)))+
  #geom_errorbar(data = d2, width=2, show.legend = FALSE, size=0.25, alpha=0.5,
  #              aes(x=mean_PC1, ymin = (mean_PC2 - sd_PC2), ymax = (mean_PC2 + sd_PC2)))+
  #geom_path(data=d2, aes(x=mean_PC1, y=mean_PC2, colour=set2,group=set2,),size=0.5, show.legend = FALSE)+
  geom_text_repel(data=d, aes(x=PC1, y=PC2, label=times),
                  size=2, segment.color = "#84919e", segment.size = 0.25, show.legend = FALSE)+
  scale_color_manual(values = color4)+
  scale_fill_manual(values = color4)+
  scale_shape_manual(values = c(22:25))+
  labs(x=xlab, y=ylab)+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7, colour = "black"),
        axis.text.y = element_text(size=7, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "top",
        legend.justification = "center")

g

ggsave(g, file=sprintf("%s/%s_FigS12b.pdf",dir.output, exec.date),
       width = 105, height = 105, units = "mm")

# Fig. S11b ----
set <- factor(c(rep("FL/FTH",17),rep("CL/CTH",17),rep("FL/CTH",17),rep("CL/FTH",17)),
              levels = c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
rpm2.ave <- log2rpm.ave_204[raw_10,]
soukan <- cor(rpm2.ave, method = "spearman")
souka <- stats::dist(1-soukan)
souk <- hclust(souka, method = "ward.D2")
souka <- NULL
dend <- as.dendrogram(souk) %>% set("branches_lwd",0.25)
labels(dend) <- souk$labels[souk$order]
labels_cex(dend) <- 0.5
cols <- rep(color4,each=17)[souk$order]
labels_colors(dend) <- cols

p1 <- ggplot(as.ggdend(dend))+
  scale_x_continuous(expand = c(0, 0.5), breaks = 1:68, labels = souk$labels[souk$order])+
  scale_y_continuous(expand = c(0, 6))
#theme(plot.margin= unit(c(1, 1, 1, 1), "lines"))


ggsave(p1, file= sprintf("%s/%s_FigS11b.pdf", dir.output, exec.date),
       width = 160, height = 80, units = "mm")

### DEGs by TCC (Sun et al., 2013; Tang et al., 2015) ----
d <- rawcnt_204[drr,]
d <- d[raw_10,]
time <- time_204

dir.create(sprintf("%s/TCC_204", dir.output))
a <- read.delim(sprintf("input/condition_204.txt", dir.output))
for(i in 1:6){
  dir.create(sprintf("%s/TCC_204/%s_%s", dir.output, a[i,1],a[i,2]))
}

# parameter
param_G1 <- 3 #the number of samples for G2 (SmartGC)
param_G2 <- 3 #the number of samples for G3 (GC)
param_G3 <- 3 #the number of samples for G4 (L)
param_G4 <- 3 #the number of samples for G5 (TH)
param_FDR <- 0.05

# setting before analysis
#c1 <- "Field"
#c2 <- "L"
#param_contrast <- c(1, 0, 0, -1, 0)         
#groups <- c(1,4) 
a <- read.delim(sprintf("input/condition_204.txt", dir.output))

for (j in 1:6){
  c1 <- a[j,1]
  c2 <- a[j,2]
  param_contrast <- c(a[j,3], a[j,4], a[j,5], a[j,6])         
  groups <- c(a[j,7],a[j,8]) 
  
  # analysis
  pdf(sprintf("%s/TCC_204/%s_MAplot_%s_%s.pdf", dir.output, exec.date, c1, c2),
      width = 10, height = 7)
  
  for (i in 1:17){
    # data preparation
    target <- at_204$hour==time[i]
    
    # preprocessing
    data <- d[,target]
    data.cl <- c(rep(1, param_G1), rep(2, param_G2), rep(3, param_G3), rep(4,param_G4))
    tcc <- new("TCC", data, data.cl)
    
    # Normalization
    tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                           iteration=3, FDR=0.1, floorPDEG=0.05)
    
    # DEG detection
    design <- model.matrix(~ 0 + ~ as.factor(data.cl))
    tcc <- estimateDE(tcc, test.method="edger", FDR=param_FDR,
                      design=design, contrast=param_contrast)
    result <- getResult(tcc, sort=FALSE)
    #sum(tcc$stat$q.value < param_FDR)
    #design
    
    # save data
    tmp <- cbind(rownames(tcc$count), tcc$count, result)
    out_f <- sprintf("%s/TCC_204/%s_%s/result_%s.txt", 
                     dir.output, c1, c2, time[i]) #output file
    write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)
    tmp <- tmp[tmp$q.value < 0.05,]
    tmp <- tmp[order(tmp$rank),]
    out_f <- sprintf("%s/TCC_204/%s_%s/result_0.05_%s.txt",
                     dir.output, c1, c2, time[i]) #output file
    write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)
    # figure 
    plot(tcc, FDR=param_FDR, group=groups, # group setting
         main=sprintf("%s_%s_%s",c1, c2, time[i]))              
    legend("topright", c(paste("DEG(FDR<", param_FDR, ")", sep=""), "non-DEG"),
           col=c("magenta", "black"), pch=20)
  }
  dev.off()
}

### TCC_list ----
time <- time_204
rpm2.ave <- rpm.ave_204[raw_10,]
TCC_204 <- matrix(NA, ncol=17, nrow=nrow(rpm2.ave))
rownames(TCC_204) <- rownames(rpm2.ave)
colnames(TCC_204) <- time
condition <- read.delim(sprintf("input/condition_204.txt", dir.output))

for (k in 1:6){
  c1 <- condition[k,1]
  c2 <- condition[k,2]
  
  for (i in 1:17){
    time <- colnames(TCC_204)[i]
    a <- read.delim(sprintf("%s/TCC_204/%s_%s/result_%s.txt",dir.output, c1, c2, time))
    b <- a$q.value
    names(b) <- a$rownames.tcc.count.
    for (j in 1:length(b)){
      TCC_204[j,i] <- b[rownames(TCC_204)[j]]
    }
    save(TCC_204, file = sprintf("%s/TCC_204/%s_%s/TCC_204_list", dir.output, c1, c2))}
}

### Table S4 ----
load("output/TCC_204/FL_FTH_CL_CTH/TCC_204_list")
a <- TCC_204
load("output/TCC_204/FL_FTH_FL_CTH/TCC_204_list")
b <- cbind(a,TCC_204)
load("output/TCC_204/FL_FTH_CL_FTH/TCC_204_list")
c <- cbind(b,TCC_204)

colnames(c) <- paste(rep(c("CL/CTH","FL/CTH","CL/FTH"),each=17), rep(time,3), sep="_")
write.table(c, file=sprintf("output/%s_TableS5.csv",exec.date),sep=",")

### TCC table ----
time <- time_204
rpm2.ave <- rpm.ave[raw_10,]
condition <- read.delim(sprintf("input/condition_204.txt", dir.output))

for(i in 1:17){
  ID <- time[i]
  TCC_204_table <- matrix(NA, nrow=nrow(rpm2.ave), ncol=6)
  colnames(TCC_204_table) <- c("FL_FTH_CL_CTH", "FL_FTH_FL_CTH","FL_FTH_CL_FTH",
                               "CL_CTH_FL_CTH", "CL_CTH_CL_FTH", "FL_CTH_CL_FTH")
  rownames(TCC_204_table) <- rownames(rpm2.ave)
  for (j in 1:6){
    c1 <- condition[j,1]
    c2 <- condition[j,2]
    a <- read.delim(sprintf("%s/TCC_204/%s_%s/result_0.05_%s.txt", dir.output, c1, c2, ID))
    TCC_204_table[,j] <- is.element(rownames(rpm2.ave), a$rownames.tcc.count.)
  }
  save(TCC_204_table, file = sprintf("%s/TCC_204/%s_TCC_204_0.05_%s", dir.output, exec.date, ID))
}

### DEG FL/FTH vs others table ----
TCC_204_type_table <- matrix(NA, nrow=3, ncol=17)
rownames(TCC_204_type_table) <- c("CL/CTH", "FL/CTH", "CL/FTH")
colnames(TCC_204_type_table) <- as.character(time_204)
time <- time_204
test2 <- function(input){
  length(which(input==T))  
}

for (i in 1:17){
  ID <- time[i]
  load(sprintf("%s/TCC_204/210703_TCC_204_0.05_%s", dir.output, ID))
  a <- TCC_204_table[,1:3]
  b <- apply(a,2,test2)
  TCC_204_type_table[,ID] <- b
}

write.csv(TCC_204_type_table, 
          file=sprintf("%s/TCC_204/%s_TCC_204_type_table_FL_FTHvsothers.csv", dir.output, exec.date))

# Fig. 2b ----
# Field vs others
TCC_204_type_table <- read.csv(sprintf("%s/TCC_204/210703_TCC_204_type_table_FL_FTHvsothers.csv",
                                       dir.output),row.names = 1)
time <- as.character(time_204)
time[17] <- "19"
time <- factor(time, levels = time)
colnames(TCC_204_type_table) <- time

CL_CTH <- t(TCC_204_type_table[1,])
FL_CTH <- t(TCC_204_type_table[2,])
CL_FTH <- t(TCC_204_type_table[3,])

n <- c(CL_CTH, FL_CTH, CL_FTH)
con <- factor(c(rep("CL/CTH", 17),rep("FL/CTH", 17),rep("CL/FTH", 17)),
              levels = c("CL/CTH","FL/CTH","CL/FTH"))
time2 <- as.numeric(as.character(rep(time,3)))

d <- data.frame(n, con, time2)

fig2b <- ggplot(d, aes(x = time2, y = n, group=con, colour=con, shape=con, fill=con))+
  geom_line(size = 0.25)+
  geom_point(color="black", size=1,stroke=0.25)+
  scale_color_manual(values = color5[3:5])+
  scale_fill_manual(values = color5[3:5])+
  scale_shape_manual(values = c(23:25))+
  scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
  labs(x="Time (hours)", y="The number of genes")+
  ylim(0,3500)+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  panel_border(colour=1, size=0.5)+ 
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "top",
        legend.justification = "center")
plot(fig2b)

#ggsave(sprintf("%s/%s_Fig2b.pdf", dir.output, exec.date), 
#       width=55, height = 40, units = "mm")

### DEG Venn diagram ----
CL_CTH <- list()
FL_CTH <- list()
CL_FTH <- list()

time2 <- time_204

for (i in time2){
  a <- read.delim(sprintf("%s/TCC_204/FL_FTH_CL_CTH/result_0.05_%s.txt",
                          dir.output,i))
  CL_CTH[[i]] <- as.character(a$rownames.tcc.count.)
  
  a <- read.delim(sprintf("%s/TCC_204/FL_FTH_FL_CTH/result_0.05_%s.txt",
                          dir.output,i))
  FL_CTH[[i]] <- as.character(a$rownames.tcc.count.)
  
  a <- read.delim(sprintf("%s/TCC_204/FL_FTH_CL_FTH/result_0.05_%s.txt",
                          dir.output,i))
  CL_FTH[[i]] <- as.character(a$rownames.tcc.count.)
  
  #ven <- list(GC = GC[[i]], L = L[[i]], TH = TH[[i]])
  
  #venn.diagram(ven,
  #             filename = sprintf("%s/venndiagram/%s_venn_%s.tiff",
  #                                dir.output,exec.date, i),
  #             fill = c(3, 2, 7),           # background color
  #             alpha = c(0.5, 0.5 , 0.5),   # transparency
  #             lty = c(1, 2, 3)             # border line type
  #)
}

LIGHT <- list()
for(i in 1:17){
  LIGHT[[i]] <- intersect(CL_CTH[[i]],CL_FTH[[i]]) 
}

TH <- list()
for(i in 1:17){
  TH[[i]] <- intersect(CL_CTH[[i]],FL_CTH[[i]]) 
}

LTH <- list()
for(i in 1:17){
  LTH[[i]] <- intersect(intersect(CL_CTH[[i]],FL_CTH[[i]]),CL_FTH[[i]])
}

names(LIGHT) <- as.character(unique(hour))
names(TH) <- as.character(unique(hour))
names(LTH) <- as.character(unique(hour))

for(i in 1:17){
  ti <- as.character(unique(hour))[i]
  tmp <- des2[unlist(LIGHT[i]),]
  write.csv(tmp, 
            file=sprintf("output/TCC_204/%s_LIGHT_%s.csv",exec.date, ti))
  tmp <- des2[unlist(TH[i]),]
  write.csv(tmp, 
            file=sprintf("output/TCC_204/%s_TH_%s.csv",exec.date, ti))
  tmp <- des2[unlist(LTH[i]),]
  write.csv(tmp,
            file=sprintf("output/TCC_204/%s_LTH_%s.csv",exec.date, ti))
}

# Fig. 2d ---- 
time <- as.character(time_204)
time[17] <- "19"
time <- factor(time, levels = time)
colnames(TCC_204_type_table) <- time

a <- rep(NA, 17)
for(i in 1:17){
  a[i] <- length(LIGHT[[i]])
}

b <- rep(NA, 17)
for(i in 1:17){
  b[i] <- length(TH[[i]])
}

c <- rep(NA, 17)
for(i in 1:17){
  c[i] <- length(LTH[[i]])
}

n <- c(a,b,c)
con <- factor(c(rep("LIGHT", 17),rep("TH", 17), rep("LTH", 17)), levels = c("LIGHT","TH","LTH"))
time2 <- as.numeric(as.character(rep(time,3)))

d <- data.frame(n, con, time2)

fig2d <- ggplot(d, aes(x = time2, y = n, group=con, colour=con, shape=con, fill=con))+
  geom_line(size = 0.25)+
  geom_point(color="black", size=1, stroke=0.25)+
  scale_color_manual(values = color_5[1:3])+
  scale_fill_manual(values = color_5[1:3])+
  scale_shape_manual(values = c(21:23))+
  scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
  labs(x="Time (hours)", y="The number of genes")+
  ylim(0,1500)+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  panel_border(colour=1, size=0.5)+ 
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "top",
        legend.justification = "center")
plot(fig2d)

#ggsave(sprintf("%s/%s_Fig2d.pdf", dir.output, exec.date), 
#       width=55, height = 40, units = "mm")

### Fig. 2a-d ----
p <- plot_grid(fig2a,fig2c, fig2b, fig2d, nrow=2, align = "hv")

ggsave(p, file = sprintf("%s/%s_Fig2a-d.pdf", dir.output, exec.date),
       width = 110, height = 80, units="mm")

# Fig. 3C ---
length(LIGHT[[6]])
length(LIGHT[[17]])
length(intersect(LIGHT[[6]],LIGHT[[17]]))

ven <- list(LIGHT_8 = LIGHT[[6]], LIGHT_19_2 = LIGHT[[17]])

venn.plot <- venn.diagram(ven,
                          filename = NULL,
                          #sprintf("%s/%s_Fig3B.tiff",dir.output,exec.date),
                          labels = c("LIGHT_8","LIGHT_19_2"),
                          fill = color_5[1:2],
                          alpha = c(0.5, 0.5),
                          lty = c(1, 2)
)

pdf(sprintf("%s/%s_Fig3C.pdf", dir.output, exec.date))  
grid.draw(venn.plot)
dev.off()

LTH[[7]]


pdf("a.pdf")
geneplot(LTH[[7]])
dev.off()


geneplot("Os03g0169600")



### GO significant ----
rpm2.ave <- rpm.ave_204[raw_10,]
use.genelist <- rownames(rpm2.ave)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))
term_na <- term[term != "NA"]
time <- time_204
condition <- read.delim(sprintf("%s/condition_204.txt", dir.input))

dir.create(sprintf("%s/GO_204", dir.output))
a <- read.delim(sprintf("%s/condition_204.txt", dir.input))
for(i in 1:6){
  dir.create(sprintf("%s/GO_204/%s_%s", dir.output, a[i,1],a[i,2]))
}

for (j in 1:6){
  
  c1 <- condition[j,1]
  c2 <- condition[j,2]
  
  for (i in 1:17){
    ID <- time[i]
    load(sprintf("%s/TCC_204/210703_TCC_204_0.05_%s", dir.output, ID))
    a <-  TCC_204_table[,j]==T
    gl <- rownames(TCC_204_table)[a]
    result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
      
    fn <- sprintf("%s/GO_204/%s_%s/%s_GO_%s_%s_%s.csv", dir.output, c1, c2, exec.date, c1, c2, ID)
      
    adp <- p.adjust(result[,"p.value"], method = "BH")
      
    tmp.id <- rownames(result)[adp < 0.05]
    tmp.adp <- adp[adp < 0.05]
    tmp.description <- ng.GetGOTerms(tmp.id)
    tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
    if (length(tmp.id)==1){
      out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
    }else{
      out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
    }
    colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
    write.csv(out, fn, row.names = F)
    
    cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  
  # merge csv file
  tmp <- NULL
  for (i in 1:17){
    
    ID <- time[i]
    a <- read.csv(sprintf("%s/GO_204/%s_%s/%s_GO_%s_%s_%s.csv", dir.output, c1, c2, 
                          exec.date, c1, c2, ID))
    if (is.na(a[1,1])){
      next
    }else{
      hour <- rep(ID, nrow(a))
      a <- cbind(hour, a)
      tmp <- rbind(tmp, a)
    }
  }
  write.csv(tmp, file= sprintf("%s/GO_204/%s_GO_%s_%s.csv", dir.output, exec.date, c1, c2), row.names = F)
  
}

### GO heatmap ----
### GO
rpm2.ave <- rpm.ave_204[raw_10,]
use.genelist <- rownames(rpm2.ave)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))
term_na <- term[term != "NA"]
time <- unique(at_204$hour)
condition <- read.delim(sprintf("%s/condition_204.txt", dir.input))
dir.create(sprintf("%s/GO_204/heatmap", dir.output))
for(i in 1:6){
  dir.create(sprintf("%s/GO_204/heatmap/%s_%s", dir.output, condition[i,1],condition[i,2]))
}

for (j in 1:6){
  c1 <- condition[j,1]
  c2 <- condition[j,2]
  
  for (i in 1:17){
    ID <- time[i]
    load(sprintf("%s/TCC_204/210703_TCC_204_0.05_%s", dir.output, ID))
    a <-  TCC_204_table[,j]==T
    gl <- rownames(TCC_204_table)[a]
    result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]

    adp <- p.adjust(result[,"p.value"], method = "BH")
    
    fn <- sprintf("%s/GO_204/heatmap/%s_%s/%s_GO_%s_%s_%s.csv", dir.output, c1, c2, exec.date, c1, c2, ID)
    tmp.adp <- adp
    tmp.id <- rownames(result)
    tmp.description <- term_na
    tmp.xnn <- result[, c("xtt", "xtn", "xnt", "xnn")]
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
    colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
    write.csv(out, fn, row.names = F)
    cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  
  # merge csv file
  tmp <- NULL
  for (i in 1:17){
    
    ID <- time[i]
    a <- read.csv(sprintf("%s/GO_204/heatmap/%s_%s/210703_GO_%s_%s_%s.csv", dir.output, c1, c2, c1, c2, ID))
    if (is.na(a[1,1])){
      next
    }else{
      hour <- rep(ID, nrow(a))
      a <- cbind(hour, a)
      tmp <- rbind(tmp, a)
    }
  }
  write.csv(tmp, file= sprintf("%s/GO_204/heatmap/%s_GO_%s_%s.csv", dir.output, exec.date, c1, c2), row.names = F)
}

### heatmap ----
a <- read.csv("output/GO_204/heatmap/210703_GO_FL_FTH_CL_CTH.csv")
a1 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d1 <- as.matrix(b[,2:18])

a1$hour <- factor(a1$hour, levels=time)
b1 <- a1 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/GO_204/heatmap/210703_GO_FL_FTH_FL_CTH.csv")
a2 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d2 <- as.matrix(b[,2:18])

a2$hour <- factor(a2$hour, levels=time)
b2 <- a2 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/GO_204/heatmap/210703_GO_FL_FTH_CL_FTH.csv")
a3 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d3 <- as.matrix(b[,2:18])

a3$hour <- factor(a3$hour, levels=time)
b3 <- a3 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)


d <- cbind(d1,d2,d3)

pvalue <- function(dat){
  (min(dat) < log10(0.05)) == T
}

dat <- d[apply(d,1,pvalue),]
dat[dat < -5] <- -5
tmp <- heatmap(x = dat, Colv = NA)
dat <- dat[tmp$rowInd,]

ID <- as.character(b$ID)[apply(d,1,pvalue)]
GO_Term <- ng.GetGOTerms(ID)[tmp$rowInd]

# Table S10 ----
dID <- rev(ID[tmp$rowInd])
tex <- matrix(NA, nrow=length(dID), ncol=9)
rownames(tex) <- dID
colnames(tex) <- c("ID","GO term","condition","time","Adjusted p-value", "A & B", "A",	"B", "U")
tex[,1] <- dID
tex[,2] <- rev(GO_Term)

cond <- rep(c("CL_CTH","FL_CTH","CL_FTH"),each=17)
tim <- rep(time, 3)

rownames(b1) <- b1$ID
rownames(b2) <- b2$ID
rownames(b3) <- b3$ID
aaa <- cbind(b1[dID,2:18], b2[dID,2:18],b3[dID,2:18])
bb <- apply(aaa,1,min)

for(i in 1:length(dID)){
  
  cc <- cond[aaa[i,] == bb[i]]
  dd <- as.character(tim)[aaa[i,] == bb[i]]
  ee <- bb[i]
  
  ff <- read.csv(sprintf("output/GO_204/heatmap/210703_GO_FL_FTH_%s.csv",cc))
  
  gg <- as.numeric(ff[(ff$hour==dd & ff$ID == dID[i]),5:8])
  
  tex[i,3] <- cc
  tex[i,4] <- dd
  tex[i,5] <- sprintf("%.2e",as.numeric(ee))
  tex[i,6:9] <- gg
}

write.csv(tex, file= sprintf("%s/%s_TableS10.csv", dir.output, exec.date),
          row.names = F)

# Fig. S16a ----
dat[dat >= log10(0.05)] <- 0

png(sprintf("%s/%s_FigS16a.png",dir.output,exec.date),
    width = 40.107, height = 60.269, units = "mm", res = 1000)

par(mar=c(0,0,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,51),xaxt="n",
     at = seq(0, 1, length.out = 51))
#axis(side = 2, labels = GO_Term, las = 1,cex.axis=0.5,
#     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.6)
#par(oma=c(0,0,0,0))
#image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0))
#par(oma = c(0,0,3,0))
#title("")

dev.off()

#
pdf(sprintf("%s/%s_FigS16a_2.pdf",dir.output,exec.date),
    width = 7, height = 12)

par(oma=c(0,10,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,51),xaxt="n",
     at = seq(0, 1, length.out = 51))
axis(side = 2, labels = GO_Term, las = 1,cex.axis=0.5,
     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.5)
par(oma=c(0,0,0,0))
image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0),
           horizontal = T)
par(oma = c(0,0,3,0))
title("")

dev.off()

### GO DEG venn significant ----
rpm2.ave <- rpm.ave_204[raw_10,]
use.genelist <- rownames(rpm2.ave)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))
term_na <- term[term != "NA"]
time <- time_204
condition <- read.delim(sprintf("%s/condition_204.txt", dir.input))
dir.create("output/venndiagram_204")
dir.create("output/venndiagram_204/LIGHT")
dir.create("output/venndiagram_204/TH")
dir.create("output/venndiagram_204/LTH")

# LIGHT
for (i in 1:17){
  ID <- time[i]
  gl <- LIGHT[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram_204/LIGHT/%s_GO_LIGHT_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/LIGHT/210703_GO_LIGHT_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/%s_GO_LIGHT.csv", dir.output, exec.date), row.names = F)

# TH
for (i in 1:17){
  ID <- time[i]
  gl <- TH[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram_204/TH/%s_GO_TH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/TH/210703_GO_TH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/%s_GO_TH.csv", dir.output, exec.date), row.names = F)

# LTH
for (i in 1:17){
  ID <- time[i]
  gl <- LTH[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram_204/LTH/%s_GO_LTH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/LTH/210703_GO_LTH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/%s_GO_LTH.csv", dir.output, exec.date), row.names = F)

### GO DEG venn heatmap ----
rpm2.ave <- rpm.ave_204[raw_10,]
use.genelist <- rownames(rpm2.ave)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))
term_na <- term[term != "NA"]
time <- time_204
condition <- read.delim(sprintf("%s/condition_204.txt", dir.input))
dir.create("output/venndiagram_204/heatmap/")
dir.create("output/venndiagram_204/heatmap/LIGHT")
dir.create("output/venndiagram_204/heatmap/TH")
dir.create("output/venndiagram_204/heatmap/LTH")

# LIGHT
for (i in 1:17){
  ID <- time[i]
  gl <- LIGHT[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram_204/heatmap/LIGHT/%s_GO_LIGHT_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)
  tmp.adp <- adp
  tmp.description <- term_na
  tmp.xnn <- result[, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/heatmap/LIGHT/210703_GO_LIGHT_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/heatmap/%s_GO_LIGHT.csv", dir.output, exec.date), row.names = F)

# TH
for (i in 1:17){
  ID <- time[i]
  gl <- TH[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram_204/heatmap/TH/%s_GO_TH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)
  tmp.adp <- adp
  tmp.description <- term_na
  tmp.xnn <- result[, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/heatmap/TH/210703_GO_TH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/heatmap/%s_GO_TH.csv", dir.output, exec.date), row.names = F)

# LTH
for (i in 1:17){
  ID <- time[i]
  gl <- LTH[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram_204/heatmap/LTH/%s_GO_LTH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)
  tmp.adp <- adp
  tmp.description <- term_na
  tmp.xnn <- result[, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/heatmap/LTH/210703_GO_LTH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/heatmap/%s_GO_LTH.csv", dir.output, exec.date), row.names = F)

### heatmap ----
a <- read.csv("output/venndiagram_204/heatmap/210703_GO_LIGHT.csv")
a1 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d1 <- as.matrix(b[,2:18])

a1$hour <- factor(a1$hour, levels=time)
b1 <- a1 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/venndiagram_204/heatmap/210703_GO_TH.csv")
a2 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d2 <- as.matrix(b[,2:18])

a2$hour <- factor(a2$hour, levels=time)
b2 <- a2 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/venndiagram_204/heatmap/210703_GO_LTH.csv")
a3 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d3 <- as.matrix(b[,2:18])

a3$hour <- factor(a3$hour, levels=time)
b3 <- a3 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)


d <- cbind(d1,d2,d3)

pvalue <- function(dat){
  (min(dat) < log10(0.05)) == T
}

dat <- d[apply(d,1,pvalue),]
dat[dat < -5] <- -5
tmp <- heatmap(x = dat, Colv = NA)
dat <- dat[tmp$rowInd,]

ID <- as.character(b$ID)[apply(d,1,pvalue)]
GO_Term <- ng.GetGOTerms(ID)[tmp$rowInd]

# Table S11 ----
dID <- rev(ID[tmp$rowInd])
tex <- matrix(NA, nrow=length(dID), ncol=9)
rownames(tex) <- dID
colnames(tex) <- c("ID","GO term","condition","time","Adjusted p-value", "A & B", "A",	"B", "U")
tex[,1] <- dID
tex[,2] <- rev(GO_Term)

cond <- rep(c("LIGHT","TH","LTH"),each=17)
tim <- rep(time, 3)

rownames(b1) <- b1$ID
rownames(b2) <- b2$ID
rownames(b3) <- b3$ID
aaa <- cbind(b1[dID,2:18], b2[dID,2:18],b3[dID,2:18])
bb <- apply(aaa,1,min)

for(i in 1:length(dID)){
  
  cc <- cond[aaa[i,] == bb[i]]
  dd <- as.character(tim)[aaa[i,] == bb[i]]
  ee <- bb[i]
  
  ff <- read.csv(sprintf("output/venndiagram_204/heatmap/210703_GO_%s.csv",cc))
  
  gg <- as.numeric(ff[(ff$hour==dd & ff$ID == dID[i]),5:8])
  
  tex[i,3] <- cc
  tex[i,4] <- dd
  tex[i,5] <- sprintf("%.2e",as.numeric(ee))
  tex[i,6:9] <- gg
}

write.csv(tex, file= sprintf("%s/%s_TableS11.csv", dir.output, exec.date),
          row.names = F)

# Fig. S16b ----
dat[dat >= log10(0.05)] <- 0

png(sprintf("%s/%s_FigS16b.png",dir.output,exec.date),
    width = 40.107, height = 36.91, units = "mm", res = 1000)

par(mar=c(0,0,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,51),xaxt="n",
     at = seq(0, 1, length.out = 51))
#axis(side = 2, labels = GO_Term, las = 1,cex.axis=0.5,
#     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.6)
#par(oma=c(0,0,0,0))
#image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0))
#par(oma = c(0,0,3,0))
#title("")

dev.off()

#
pdf(sprintf("%s/%s_FigS16b_2.pdf",dir.output,exec.date),
    width = 7, height = 12)

par(oma=c(0,10,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,51),xaxt="n",
     at = seq(0, 1, length.out = 51))
axis(side = 2, labels = GO_Term, las = 1,cex.axis=0.5,
     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.5)
par(oma=c(0,0,0,0))
image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0),
           horizontal = T)
par(oma = c(0,0,3,0))
title("")

dev.off()

### KEGG significant ----
rpm2.ave <- rpm.ave_204[raw_10,]
use.genelist <- rownames(rpm2.ave)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]
time <- unique(at_204$hour)
condition <- read.delim(sprintf("%s/condition_204.txt", dir.input))
dir.create(sprintf("%s/KEGG_204", dir.output))
for(i in 1:6){
  dir.create(sprintf("%s/KEGG_204/%s_%s", dir.output, condition[i,1],condition[i,2]))
}

for (j in 1:6){
  # Field vs SmartGC
  c1 <- condition[j,1]
  c2 <- condition[j,2]
  
  for (i in 1:17){
    ID <- time[i]
    load(sprintf("%s/TCC_204/210703_TCC_204_0.05_%s", dir.output, ID))
    a <-  TCC_204_table[,j]==T
    gl <- rownames(TCC_204_table)[a]
    result <- kegg.mft(cgt=kegg, gn.test=gl)
    
    fn <- sprintf("%s/KEGG_204/%s_%s/%s_KEGG_%s_%s_%s.csv", dir.output, c1, c2, exec.date, c1, c2, ID)
    
    adp <- p.adjust(result[,"p.value"], method = "BH")
    
    tmp.id <- rownames(result)[adp < 0.05]
    tmp.adp <- adp[adp < 0.05]
    tmp.description <- kegg_description[tmp.id]
    tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
    if (length(tmp.id)==1){
      out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
    }else{
      out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
    }
    colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
    write.csv(out, fn, row.names = F)
    cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  
  # merge csv file
  tmp <- NULL
  for (i in 1:17){
    
    ID <- time[i]
    a <- read.csv(sprintf("%s/KEGG_204/%s_%s/210703_KEGG_%s_%s_%s.csv", dir.output, c1, c2, c1, c2, ID))
    if (is.na(a[1,1])){
      next
    }else{
      hour <- rep(ID, nrow(a))
      a <- cbind(hour, a)
      tmp <- rbind(tmp, a)
    }
    #cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  write.csv(tmp, file= sprintf("%s/KEGG_204/%s_KEGG_%s_%s.csv", dir.output, exec.date, c1, c2), row.names = F)
}


### KEGG heatmap ----
rpm2.ave <- rpm.ave_204[raw_10,]
use.genelist <- rownames(rpm2.ave)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]
time <- unique(at_204$hour)
condition <- read.delim(sprintf("%s/condition_204.txt", dir.input))
dir.create(sprintf("%s/KEGG_204/heatmap", dir.output))
for(i in 1:6){
  dir.create(sprintf("%s/KEGG_204/heatmap/%s_%s", dir.output, condition[i,1],condition[i,2]))
}

for (j in 1:6){
  # Field vs SmartGC
  c1 <- condition[j,1]
  c2 <- condition[j,2]
  
  for (i in 1:17){
    ID <- time[i]
    load(sprintf("%s/TCC_204/210703_TCC_204_0.05_%s", dir.output, ID))
    a <-  TCC_204_table[,j]==T
    gl <- rownames(TCC_204_table)[a]
    result <- kegg.mft(cgt=kegg, gn.test=gl)
    
    fn <- sprintf("%s/KEGG_204/heatmap/%s_%s/%s_KEGG_%s_%s_%s.csv", dir.output, c1, c2, exec.date, c1, c2, ID)
    
    adp <- p.adjust(result[,"p.value"], method = "BH")
    
    tmp.id <- rownames(result)
    tmp.adp <- adp
    tmp.description <- kegg_description[tmp.id]
    tmp.xnn <- result[, c("xtt", "xtn", "xnt", "xnn")]
    if (length(tmp.id)==1){
      out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
    }else{
      out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
    }
    colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
    write.csv(out, fn, row.names = F)
    cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  
  # merge csv file
  tmp <- NULL
  for (i in 1:17){
    
    ID <- time[i]
    a <- read.csv(sprintf("%s/KEGG_204/heatmap/%s_%s/210703_KEGG_%s_%s_%s.csv", dir.output, c1, c2, c1, c2, ID))
    if (is.na(a[1,1])){
      next
    }else{
      hour <- rep(ID, nrow(a))
      a <- cbind(hour, a)
      tmp <- rbind(tmp, a)
    }
    #cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  write.csv(tmp, file= sprintf("%s/KEGG_204/heatmap/%s_KEGG_%s_%s.csv", dir.output, exec.date, c1, c2), row.names = F)
}

# heatmap ----
a <- read.csv("output/KEGG_204/heatmap/210703_KEGG_FL_FTH_CL_CTH.csv")
a1 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d1 <- as.matrix(b[,2:18])

a1$hour <- factor(a1$hour, levels=time)
b1 <- a1 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/KEGG_204/heatmap/210703_KEGG_FL_FTH_FL_CTH.csv")
a2 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d2 <- as.matrix(b[,2:18])

a2$hour <- factor(a2$hour, levels=time)
b2 <- a2 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)


a <- read.csv("output/KEGG_204/heatmap/210703_KEGG_FL_FTH_CL_FTH.csv")
a3 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d3 <- as.matrix(b[,2:18])

a3$hour <- factor(a3$hour, levels=time)
b3 <- a3 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d <- cbind(d1,d2,d3)

pvalue <- function(dat){
  (min(dat) < log10(0.05)) == T
}

dat <- d[apply(d,1,pvalue),] # at least 1 time is significant
dat[dat < -5] <- -5
tmp <- heatmap(x = dat, Colv = NA)
dat <- dat[tmp$rowInd,]

ID <- as.character(b$ID)[apply(d,1,pvalue)]
KEGG_Term <- kegg_description[ID][tmp$rowInd]

# Table S12 ----
dID <- rev(ID[tmp$rowInd])
tex <- matrix(NA, nrow=length(dID), ncol=9)
rownames(tex) <- dID
colnames(tex) <- c("ID","KEGG term","condition","time","Adjusted p-value", "A & B", "A",	"B", "U")
tex[,1] <- dID
tex[,2] <- rev(KEGG_Term)

cond <- rep(c("CL_CTH","FL_CTH","CL_FTH"),each=17)
tim <- rep(time, 3)

rownames(b1) <- b1$ID
rownames(b2) <- b2$ID
rownames(b3) <- b3$ID
aaa <- cbind(b1[dID,2:18], b2[dID,2:18],b3[dID,2:18])
bb <- apply(aaa,1,min)

for(i in 1:length(dID)){
  
  cc <- cond[aaa[i,] == bb[i]]
  dd <- as.character(tim)[aaa[i,] == bb[i]]
  ee <- bb[i]
  
  ff <- read.csv(sprintf("output/KEGG_204/heatmap/210703_KEGG_FL_FTH_%s.csv",cc))
  
  gg <- as.numeric(ff[(ff$hour==dd & ff$ID == dID[i]),5:8])
  
  tex[i,3] <- cc
  tex[i,4] <- dd
  tex[i,5] <- sprintf("%.2e",as.numeric(ee))
  tex[i,6:9] <- gg
}

write.csv(tex, file= sprintf("%s/%s_TableS12.csv", dir.output, exec.date),
          row.names = F)

# Fig. S16c ----
dat[dat >= log10(0.05)] <- 0

png(sprintf("%s/%s_FigS13c.png",dir.output,exec.date),
    width = 95.227, height = 66, units = "mm", res = 1000)

par(mar=c(0,0,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,51),xaxt="n",
     at = seq(0, 1, length.out = 51))
#axis(side = 2, labels = GO_Term, las = 1,cex.axis=0.5,
#     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.6)
#par(oma=c(0,0,0,0))
#image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0))
#par(oma = c(0,0,3,0))
#title("")

dev.off()

#
pdf(sprintf("%s/%s_FigS16c_2.pdf",dir.output,exec.date),
    width = 7, height = 12)

par(oma=c(0,10,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,51),xaxt="n",
     at = seq(0, 1, length.out = 51))
axis(side = 2, labels = KEGG_Term, las = 1,cex.axis=0.5,
     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.5)
par(oma=c(0,0,0,0))
image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0),
           horizontal = T)
par(oma = c(0,0,3,0))
title("")

dev.off()

### KEGG DEG venn significant ----
rpm2.ave <- rpm.ave_204[raw_10,]
use.genelist <- rownames(rpm2.ave)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]
time <- unique(at_204$hour)
condition <- read.delim(sprintf("%s/condition_204.txt", dir.input))

# LIGHT
for (i in 1:17){ 
  ID <- time[i]
  gl <- LIGHT[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram_204/LIGHT/%s_KEGG_LIGHT_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/LIGHT/210703_KEGG_LIGHT_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/%s_KEGG_LIGHT.csv", dir.output, exec.date), row.names = F)

# TH
for (i in 1:17){ 
  ID <- time[i]
  gl <- TH[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram_204/TH/%s_KEGG_TH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/TH/210703_KEGG_TH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/%s_KEGG_TH.csv", dir.output, exec.date), row.names = F)

# LTH
for (i in 1:17){ 
  ID <- time[i]
  gl <- LTH[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram_204/LTH/%s_KEGG_LTH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/LTH/210703_KEGG_LTH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/%s_KEGG_LTH.csv", dir.output, exec.date), row.names = F)

### KEGG DEG venn heatmap ----
rpm2.ave <- rpm.ave_204[raw_10,]
use.genelist <- rownames(rpm2.ave)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]
time <- unique(at_204$hour)
condition <- read.delim(sprintf("%s/condition_204.txt", dir.input))

# LIGHT
for (i in 1:17){ 
  ID <- time[i]
  gl <- LIGHT[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram_204/heatmap/LIGHT/%s_KEGG_LIGHT_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)
  tmp.adp <- adp
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[, c("xtt", "xtn", "xnt", "xnn")]
  out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/heatmap/LIGHT/210703_KEGG_LIGHT_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/heatmap/%s_KEGG_LIGHT.csv", dir.output, exec.date), row.names = F)

# TH
for (i in 1:17){ 
  ID <- time[i]
  gl <- TH[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram_204/heatmap/TH/%s_KEGG_TH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)
  tmp.adp <- adp
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[, c("xtt", "xtn", "xnt", "xnn")]
  out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/heatmap/TH/210703_KEGG_TH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/heatmap/%s_KEGG_TH.csv", dir.output, exec.date), row.names = F)

# LTH
for (i in 1:17){ 
  ID <- time[i]
  gl <- LTH[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram_204/heatmap/LTH/%s_KEGG_LTH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)
  tmp.adp <- adp
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[, c("xtt", "xtn", "xnt", "xnn")]
  out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in 1:17){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram_204/heatmap/LTH/210703_KEGG_LTH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram_204/heatmap/%s_KEGG_LTH.csv", dir.output, exec.date), row.names = F)


### heatmap ----
a <- read.csv("output/venndiagram_204/heatmap/210703_KEGG_LIGHT.csv")
a1 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d1 <- as.matrix(b[,2:18])

a1$hour <- factor(a1$hour, levels=time)
b1 <- a1 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/venndiagram_204/heatmap/210703_KEGG_TH.csv")
a2 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d2 <- as.matrix(b[,2:18])

a2$hour <- factor(a2$hour, levels=time)
b2 <- a2 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)


a <- read.csv("output/venndiagram_204/heatmap/210703_KEGG_LTH.csv")
a3 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d3 <- as.matrix(b[,2:18])

a3$hour <- factor(a3$hour, levels=time)
b3 <- a3 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d <- cbind(d1,d2,d3)

pvalue <- function(dat){
  (min(dat) < log10(0.05)) == T
}

dat <- d[apply(d,1,pvalue),] # at least 1 time is significant
dat[dat < -5] <- -5
tmp <- heatmap(x = dat, Colv = NA)
dat <- dat[tmp$rowInd,]

ID <- as.character(b$ID)[apply(d,1,pvalue)]
KEGG_Term <- kegg_description[ID][tmp$rowInd]

# Table S13 ----
dID <- rev(ID[tmp$rowInd])
tex <- matrix(NA, nrow=length(dID), ncol=9)
rownames(tex) <- dID
colnames(tex) <- c("ID","KEGG term","condition","time","Adjusted p-value","A & B", "A",	"B", "U")
tex[,1] <- dID
tex[,2] <- rev(KEGG_Term)

cond <- rep(c("LIGHT","TH","LTH"),each=17)
tim <- rep(time, 3)

rownames(b1) <- b1$ID
rownames(b2) <- b2$ID
rownames(b3) <- b3$ID
aaa <- cbind(b1[dID,2:18], b2[dID,2:18],b3[dID,2:18])
bb <- apply(aaa,1,min)

for(i in 1:length(dID)){
  
  cc <- cond[aaa[i,] == bb[i]]
  dd <- as.character(tim)[aaa[i,] == bb[i]]
  ee <- bb[i]
  
  ff <- read.csv(sprintf("output/venndiagram_204/heatmap/210703_KEGG_%s.csv",cc))
  
  gg <- as.numeric(ff[(ff$hour==dd & ff$ID == dID[i]),5:8])
  
  tex[i,3] <- cc
  tex[i,4] <- dd
  tex[i,5] <- sprintf("%.2e",as.numeric(ee))
  tex[i,6:9] <- gg
}

write.csv(tex, file= sprintf("%s/%s_TableS13.csv", dir.output, exec.date),
          row.names = F)

# Fig. S16d ----
dat[dat >= log10(0.05)] <- 0

png(sprintf("%s/%s_FigS16d.png",dir.output,exec.date),
    width = 95.227, height = 35.538, units = "mm", res = 1000)

par(mar=c(0,0,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,51),xaxt="n",
     at = seq(0, 1, length.out = 51))
#axis(side = 2, labels = GO_Term, las = 1,cex.axis=0.5,
#     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.6)
#par(oma=c(0,0,0,0))
#image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0))
#par(oma = c(0,0,3,0))
#title("")

dev.off()

#
pdf(sprintf("%s/%s_FigS16d_2.pdf",dir.output,exec.date),
    width = 7, height = 12)

par(oma=c(0,10,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,51),xaxt="n",
     at = seq(0, 1, length.out = 51))
axis(side = 2, labels = KEGG_Term, las = 1,cex.axis=0.5,
     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.5)
par(oma=c(0,0,0,0))
image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0),
           horizontal = T)
par(oma = c(0,0,3,0))
title("")

dev.off()

# Fig. S20 ----
rpm2.ave <- rpm.ave[raw_10,]
withgo <- ulg[ulg[,"GOid"]=="GO:0000785", "locus"]
usegene <- intersect(withgo, rownames(rpm2.ave)) 
log2rpm_chromatin <- log2rpm.ave_204[usegene,]

con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- as.numeric(time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

ylab <- expression(paste(log[2],"(rpm)"))

### linegraph 
# FL/FTH 
con <- NULL
for(i in 1:35){
  con <- c(con, rep(usegene[i],17))
}

value <- NULL
for(i in 1:35){
  value <- c(value, as.numeric(log2rpm_chromatin[i,1:17]))
}

hours <- as.numeric(rep(time,35))
d <- data.frame(con, value, hours)

g2 <- ggplot(d,aes(x=hours, y=value, group=con))+
  geom_line(colour=color5[2],size=0.5)+
  scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
  labs(x="Time (hours)", y=ylab)+
  ggtitle("FL/FTH")+
  theme_cowplot(font_size = 6, line_size = 0.5)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")

# CL/CTH 
con <- NULL
for(i in 1:35){
  con <- c(con, rep(usegene[i],17))
}

value <- NULL
for(i in 1:35){
  value <- c(value, as.numeric(log2rpm_chromatin[i,18:34]))
}

hours <- as.numeric(rep(time,35))
d <- data.frame(con, value, hours)

g3 <- ggplot(d,aes(x=hours, y=value, group=con))+
  geom_line(colour=color5[3],size=0.5)+
  scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
  labs(x="Time (hours)", y=ylab)+
  ggtitle("CL/CTH")+
  theme_cowplot(font_size = 6, line_size = 0.5)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")

# FL/CTH 
con <- NULL
for(i in 1:35){
  con <- c(con, rep(usegene[i],17))
}

value <- NULL
for(i in 1:35){
  value <- c(value, as.numeric(log2rpm_chromatin[i,35:51]))
}

hours <- as.numeric(rep(time,35))
d <- data.frame(con, value, hours)

g4 <- ggplot(d,aes(x=hours, y=value, group=con))+
  geom_line(colour=color5[4],size=0.5)+
  scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
  labs(x="Time (hours)", y=ylab)+
  ggtitle("FL/CTH")+
  theme_cowplot(font_size = 6, line_size = 0.5)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")

# CL/FTH 
con <- NULL
for(i in 1:35){
  con <- c(con, rep(usegene[i],17))
}

value <- NULL
for(i in 1:35){
  value <- c(value, as.numeric(log2rpm_chromatin[i,52:68]))
}

hours <- as.numeric(rep(time,35))
d <- data.frame(con, value, hours)

g5 <- ggplot(d,aes(x=hours, y=value, group=con))+
  geom_line(colour=color5[5],size=0.5)+
  scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
  labs(x="Time (hours)", y=ylab)+
  ggtitle("CL/FTH")+
  theme_cowplot(font_size = 6, line_size = 0.5)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")


p <- plot_grid(g1_chro,g2_chro,g3_chro,g4_chro,g5_chro,NULL,g2,g3,g4,g5,
          ncol=5, align="hv")

ggsave(p, file=sprintf("%s/%s_FigS20.pdf",dir.output, exec.date),
       width = 160, height = 80, unit="mm")
