####################
### Experiment_1 ###
####################

# label ----
condition <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- factor(unique(at$hour), levels=unique(at$hour))
hour <- factor(at$hour, levels = time)

##############
### Figure ###
##############

# Fig. 1d ----
# Field vs others same time
FL_FTH <- rep(NA, 13)
for(i in 1:13){
  FL_FTH[i] <- cor(log2rpm.ave[raw_10,1:13],log2rpm.ave[raw_10,14:26])[i,i]
}

CL_CTH <- rep(NA, 13)
for(i in 1:13){
  CL_CTH[i] <- cor(log2rpm.ave[raw_10,1:13],log2rpm.ave[raw_10,27:39])[i,i]
}

FL_CTH <- rep(NA, 13)
for(i in 1:13){
  FL_CTH[i] <- cor(log2rpm.ave[raw_10,1:13],log2rpm.ave[raw_10,40:52])[i,i]
}

CL_FTH <- rep(NA, 13)
for(i in 1:13){
  CL_FTH[i] <- cor(log2rpm.ave[raw_10,1:13],log2rpm.ave[raw_10,53:65])[i,i]
}

value <- c(FL_FTH, CL_CTH, FL_CTH, CL_FTH)
condition <- c(rep("FL/FTH", 13),rep("CL/CTH", 13),rep("FL/CTH", 13),rep("CL/FTH", 13))
d <- data.frame(cbind(condition, as.numeric(value)))
d$condition <- factor(d$condition, levels = c("FL/FTH", "CL/CTH", "FL/CTH", "CL/FTH"))

# wilcoxon rank sum test
g <- ggplot(d,aes(x=condition,y=value))+
  geom_boxplot(outlier.colour = NA, size=0.25)+
  geom_beeswarm(aes(color=condition), size=0.5)+
  scale_color_manual(values = color4)+
  geom_signif(comparisons = list(c("FL/FTH", "CL/CTH")), annotations = "3.75e-3",
              y_position = 0.93, tip_length = 0.03, size=0.25, textsize = 2) +
  geom_signif(comparisons = list(c("FL/FTH", "FL/CTH")), annotations="2.11e-2",
              y_position = 0.94, tip_length = 0.03, size=0.25, textsize = 2) +
  geom_signif(comparisons = list(c("FL/FTH", "CL/FTH")), annotations="1.39e-1",
              y_position = 0.95, tip_length = 0.03, size=0.25, textsize = 2) +
  ylim(c(0.84,0.96))+
  labs(x="", y="Correlation coefficient (r)")+
  theme_cowplot(font_size = 6, line_size = 0.25)+
  #panel_border(colour=1, size=0.5)+
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6, colour = "black"),
        axis.text.y = element_text(size=6, colour = "black"),
        legend.position = 'none')
plot(g)

ggsave(g, file=sprintf("%s/%s_Fig1d.pdf",dir.output, exec.date),
       height = 40, width = 45, units = "mm")


pv <- rep(NA,3)

for(i in 1:3){
  cond <- unique(condition)[i+1]
  v <- value[c(1:(13), (i*13+1):((i+1)*13))]
  b <- as.factor(c(rep("FL/FTH",13),rep(cond,13)))
  d <- data.frame(v,b)
  
  a <- wilcox_test(v ~ b, d, distribution="exact")
  
  pv[i] <- pvalue_interval(a)[2]
  cat(sprintf("%s\n",i))
}

ad.pv <- p.adjust(pv, method = "BH") # 0.003745169 0.021065996 0.138868527

# Fig .1ef ----
tmp <- log2rpm.ave[raw_10,]
condition <- read.delim(sprintf("input/condition.txt", dir.output))
cor_tmp <- matrix(NA, nrow=nrow(tmp), ncol=10)
cor_tmp_pvalue <- matrix(NA, nrow=nrow(tmp), ncol=10)
cor_tmp_adjustedp <- matrix(NA, nrow=nrow(tmp), ncol=10)
rownames(cor_tmp) <- rownames(tmp)
rownames(cor_tmp_pvalue) <- rownames(tmp)
rownames(cor_tmp_adjustedp) <- rownames(tmp)

FIELD <- tmp[,c(1:13)]
FL_FTH <- tmp[,c(14:26)]
CL_CTH <- tmp[,c(27:39)]
FL_CTH <- tmp[,c(40:52)]
CL_FTH <- tmp[,c(53:65)]
co <- list(FIELD,FL_FTH,CL_CTH,FL_CTH,CL_FTH)

for (j in 1:10){
  b <-  co[[condition[j,8]]]
  c <-  co[[condition[j,9]]]
  for (i in 1:nrow(tmp)){
    cor_tmp[i,j] <- cor(as.numeric(b[i,]), as.numeric(c[i,]))
    if (is.na(cor_tmp[i,j])){next}
    cor_tmp_pvalue[i,j] <- cor.test(as.numeric(b[i,]), as.numeric(c[i,]))[[3]]
  } 
}

cor_tmp_adjustedp <- apply(cor_tmp_pvalue,2, p.adjust, method = "BH")
cor_sig <- cor_tmp_adjustedp < 0.05
cor_0.7 <- cor_tmp > 0.7

tcount <- function(x){
  length(which(x)==T)
}

gn <- rep(NA,4)
gn2 <- rep(NA,4)

for (i in 1:4){
  gn[i] <- length(which(cor_sig[,i]==T))
  gn2[i] <- length(which(cor_0.7[,i]==T))
}

nam <- c(rep("q-value < 0.05",4),rep("r > 0.7",4))
Condition <- rep(c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"),2)
Condition <- factor(Condition, levels = c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
v <- c(gn,gn2)
d <- data.frame(nam, v, Condition)

g <- ggplot(d, aes(x = nam, y = v, fill = Condition))+
  geom_bar(stat = "identity", position = "dodge", color="black", size=0.25)+
  scale_fill_manual(values = color4)+
  labs(x="", y="The number of genes")+
  theme_cowplot(font_size = 6, line_size=0.25)+
  theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = "top",
        legend.justification = 0)+
  guides(fill=guide_legend(ncol=2))
plot(g)

###
gg <- list()
for(i in 1:4){
a <- cor_tmp[,i]
df <- data.frame(value = a)
qvalue <- sort(cor_tmp[,i], decreasing = T)[length(which((cor_tmp_adjustedp[,i] < 0.05) ==T))]
cond <- unique(Condition)[i]
g1 <- ggplot(df, aes(x = value))+
     geom_histogram(binwidth = 0.05, fill = color5[i+1], colour = "black",size=0.25)+
     geom_vline(xintercept = c(0.7,qvalue), color = c("black",color5[1]), size = 0.25, linetype="dashed")+
     xlim(-1,1)+ 
     ylim(0,1500)+
     labs(x="Correlation coefficient (r)", y="The number\nof genes")+
     ggtitle(sprintf("FIELD vs %s",cond))+
     theme_cowplot(font_size = 6, line_size=0.25)+
     theme(plot.title = element_text(size=6, hjust = 0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),   
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")
gg <- c(gg, list(g1))
}

p2 <- plot_grid(gg[[1]],gg[[2]],gg[[3]],gg[[4]], ncol = 2)
p <- plot_grid(p2,g,ncol=2,rel_widths = c(2,1))

ggsave(p, file=sprintf("%s/%s_Fig1ef.pdf",dir.output, exec.date),
       width = 110, height = 40, units = "mm")

# Fig. 1g ----
pca <- prcomp(t(log2rpm[raw_10,]), scale = T)
pc <- pca$x
#PC1 vs PC2
xpc <- pca$sdev[1]^2 / sum(pca$sdev^2) *100
xlab <- sprintf("PC1 (%2.1f%%)",xpc)
ypc <- pca$sdev[2]^2 / sum(pca$sdev^2) *100
ylab <- sprintf("PC2 (%2.1f%%)",ypc)

times <- hour
set <- factor(c(rep("FIELD",52),rep("FL/FTH",52),rep("CL/CTH",52),rep("FL/CTH",52),rep("CL/FTH",52)),
              levels=c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
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

times2 <- rep(time, 5)
set2 <- factor(c(rep("FIELD",13),rep("FL/FTH",13),rep("CL/CTH",13),rep("FL/CTH",13),rep("CL/FTH",13)),
               levels=c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
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
  scale_color_manual(values = color5)+
  scale_fill_manual(values = color5)+
  scale_shape_manual(values = c(21:25))+
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

ggsave(g, file=sprintf("%s/%s_Fig1g.pdf",dir.output, exec.date),
       width = 45, height = 45, units = "mm")

# Fig. 1i ----
set <- factor(c(rep("FIELD",13),rep("FL/FTH",13),rep("CL/CTH",13),rep("FL/CTH",13),rep("CL/FTH",13)),
              levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
rpm2.ave <- log2rpm.ave[raw_10,]
soukan <- cor(rpm2.ave, method = "spearman")
souka <- stats::dist(1-soukan)
souk <- hclust(souka, method = "ward.D2")
souka <- NULL
dend <- as.dendrogram(souk) %>% set("branches_lwd",0.25)
df2 <- data.frame(states = factor(souk$labels, levels = souk$labels[souk$order]))
labels(dend) <- rep(NA, 65)

p1 <- ggplot(as.ggdend(dend))+
  scale_x_continuous(expand = c(0, 0.5), breaks = 1:65, labels = NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(plot.margin= unit(c(1, 0, -1, 0), "lines"))

# Tiles and labels
p2 <- ggplot(df2,aes(states, y = 1, fill = rep(time,5))) +
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

ggsave(a, file= sprintf("%s/%s_Fig1i1.pdf", dir.output, exec.date),
       width = 100, height = 60, units = "mm")

dend2 <- as.dendrogram(souk) %>% set("branches_lwd",0.25) %>% set("labels_cex",0.4)
p3 <- ggplot(as.dendrogram(dend2))+
  scale_x_discrete(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 2))+
  theme(plot.margin= unit(c(1, 1, 1, 1), "lines"))

ggsave(p3, file= sprintf("%s/%s_Fig1i2.pdf", dir.output, exec.date),
       width = 200, height = 100, units = "mm")

### DEGs by TCC (Sun et al., 2013; Tang et al., 2015) ----
d <- rawcnt[drr,]
d <- d[raw_10,]

# parameter
param_G1 <- 4 #the number of samples for G1 (FIELD)
param_G2 <- 4 #the number of samples for G2 (FL/FTH)
param_G3 <- 4 #the number of samples for G3 (CL/CTH)
param_G4 <- 4 #the number of samples for G4 (FL/CTH)
param_G5 <- 4 #the number of samples for G5 (CL/FTH)
param_FDR <- 0.05

dir.create(sprintf("%s/TCC", dir.output))
a <- read.delim("input/condition.txt")
for(i in 1:10){
  dir.create(sprintf("%s/TCC/%s_%s", dir.output, a[i,1],a[i,2]))
}

for (j in 1:10){
  c1 <- a[j,1]
  c2 <- a[j,2]
  param_contrast <- c(a[j,3], a[j,4], a[j,5], a[j,6], a[j,7])         
  groups <- c(a[j,8],a[j,9]) 
  
  # analysis
  pdf(sprintf("%s/TCC/%s_MAplot_%s_%s.pdf", dir.output, exec.date, c1, c2),
      width = 10, height = 7)
  
  for (i in 1:13){
    # data preparation
    target <- at$hour==time[i]
    
    # preprocessing
    data <- d[,target]
    data.cl <- c(rep(1, param_G1), rep(2, param_G2), rep(3, param_G3), rep(4,param_G4), rep(5,param_G5))
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
    out_f <- sprintf("%s/TCC/%s_%s/result_%s.txt", 
                     dir.output, c1, c2, time[i]) #output file
    write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)
    tmp <- tmp[tmp$q.value < 0.05,]
    tmp <- tmp[order(tmp$rank),]
    out_f <- sprintf("%s/TCC/%s_%s/result_0.05_%s.txt",
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
rpm2.ave <- rpm.ave[raw_10,]
TCC <- matrix(NA, ncol=13, nrow=nrow(rpm2.ave))
rownames(TCC) <- rownames(rpm2.ave)
colnames(TCC) <- time
condition <- read.delim("input/condition.txt")

for (k in 1:10){
  c1 <- condition[k,1]
  c2 <- condition[k,2]
  
  for (i in 1:13){
    time <- colnames(TCC)[i]
    a <- read.delim(sprintf("%s/TCC/%s_%s/result_%s.txt",dir.output, c1, c2, time))
    b <- a$q.value
    names(b) <- a$rownames.tcc.count.
    for (j in 1:length(b)){
      TCC[j,i] <- b[rownames(TCC)[j]]
    }
    save(TCC, file = sprintf("%s/TCC/%s_%s/TCC_list", dir.output, c1, c2))}
}

### TCC table ----
time <- factor(unique(at$hour), levels=unique(at$hour))
condition <- read.delim(sprintf("input/condition.txt", dir.output))

for(i in 1:13){
  ID <- time[i]
  TCC_table <- matrix(NA, nrow=nrow(rpm2.ave), ncol=10)
  colnames(TCC_table) <- c("FIELD_FL_FTH", "FIELD_CL_CTH", "FIELD_FL_CTH", "FIELD_CL_FTH", 
                           "FL_FTH_CL_CTH", "FL_FTH_FL_CTH", "FL_FTH_CL_FTH", 
                           "CL_CTH_FL_CTH", "CL_CTH_CL_FTH", "FL_CTH_CL_FTH")
  rownames(TCC_table) <- rownames(rpm2.ave)
  for (j in 1:10){
    c1 <- condition[j,1]
    c2 <- condition[j,2]
    a <- read.delim(sprintf("%s/TCC/%s_%s/result_0.05_%s.txt", dir.output, c1, c2, ID))
    TCC_table[,j] <- is.element(rownames(rpm2.ave), a$rownames.tcc.count.)
  }
  save(TCC_table, file = sprintf("%s/TCC/%s_TCC_0.05_%s", dir.output, exec.date, ID))
}

### Table S3 ----
load("output/TCC/FIELD_FL_FTH/TCC_list")
a <- TCC
load("output/TCC/FIELD_CL_CTH/TCC_list")
b <- cbind(a,TCC)
load("output/TCC/FIELD_FL_CTH/TCC_list")
c <- cbind(b,TCC)
load("output/TCC/FIELD_CL_FTH/TCC_list")
d <- cbind(c,TCC)

colnames(d) <- paste(rep(c("FL/FTH","CL/CTH","FL/CTH","CL/FTH"),each=13), rep(time,4), sep="_")

write.table(d, file=sprintf("output/%s_TableS3.csv",exec.date),sep=",")

### DEG Field vs others table ----
time <- factor(unique(at$hour), levels=unique(at$hour))
TCC_type_table <- matrix(NA, nrow=4, ncol=13)
rownames(TCC_type_table) <- c("FL/FTH","CL/CTH", "FL/CTH", "CL/FTH")
colnames(TCC_type_table) <- unique(at$hour)
test2 <- function(input){
  length(which(input==T))  
}

for (i in 1:13){
  ID <- time[i]
  load(sprintf("%s/TCC/210703_TCC_0.05_%s", dir.output, ID))
  a <- TCC_table[,1:4]
  b <- apply(a,2,test2)
  TCC_type_table[,i] <- b
}

write.csv(TCC_type_table, 
          file=sprintf("%s/TCC/%s_TCC_type_table_FIELDvsothers.csv",
                       dir.output, exec.date))

# Fig. 2a ----
TCC_type_table <- read.csv(sprintf("%s/TCC/210703_TCC_type_table_FIELDvsothers.csv", dir.output), row.names = 1)
colnames(TCC_type_table) <- factor(unique(at$hour), levels=unique(at$hour))
condition <- c("FL/FTH","CL/CTH","FL/CTH","CL/FTH")

a <- as.matrix(t(TCC_type_table[1:4,]))
a <- as.vector(a)
hours <- factor(rep(time,4),levels=unique(time))
condition <- factor(c(rep("FL/FTH",13),rep("CL/CTH",13),rep("FL/CTH",13),rep("CL/FTH",13)),
                    levels = c("FL/FTH","CL/CTH","FL/CTH", "CL/FTH"))

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

d <- data.frame(condition, hours, value = a)

fig2a <- ggplot(d,
            aes(x=hours, y=value, group=condition, color=condition, fill=condition, shape=condition))+
  geom_line(size = 0.25)+
  geom_point(color="black", size=1, stroke=0.25)+
  scale_color_manual(values = color4)+
  scale_fill_manual(values = color4)+
  scale_shape_manual(values = c(22:25))+
  scale_x_discrete(labels = times)+
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
        legend.justification = "center")+
  guides(color=guide_legend(ncol=2,byrow = T))
plot(fig2a)

#ggsave(sprintf("%s/%s_Fig2a.pdf", dir.output, exec.date), 
#    width=55, height = 40, units = "mm")

### DEG Venn diagram----
LIGHT <- list()
TH <- list()
LTH <- list()
UNREP <- list()

for (i in 1:13){
  
  ti <- as.character(unique(hour))[i]
  a <- read.delim(sprintf("output/TCC/FIELD_FL_FTH/result_0.05_%s.txt",ti))
  FL_FTH <- as.character(a$rownames.tcc.count.)
  
  a <- read.delim(sprintf("output/TCC/FIELD_CL_CTH/result_0.05_%s.txt",ti))
  CL_CTH <- as.character(a$rownames.tcc.count.)
  
  a <- read.delim(sprintf("output/TCC/FIELD_FL_CTH/result_0.05_%s.txt",ti))
  FL_CTH <- as.character(a$rownames.tcc.count.)
  
  a <- read.delim(sprintf("output/TCC/FIELD_CL_FTH/result_0.05_%s.txt",ti))
  CL_FTH <- as.character(a$rownames.tcc.count.)
  
  GT <- intersect(CL_CTH, CL_FTH)
  GL <- intersect(CL_CTH, FL_CTH)
  
  LIGHT <- c(LIGHT, list(setdiff(setdiff(GT, FL_FTH), FL_CTH)))
  TH <- c(TH, list(setdiff(setdiff(GL, FL_FTH), CL_FTH)))
  LTH <- c(LTH, list(setdiff(setdiff(GT, FL_FTH), unlist(LIGHT[i]))))
  UNREP <- c(UNREP, list(intersect(intersect(FL_FTH, CL_CTH), intersect(FL_CTH, CL_FTH))))
}

names(LIGHT) <- as.character(unique(hour))
names(TH) <- as.character(unique(hour))
names(LTH) <- as.character(unique(hour))
names(UNREP) <- as.character(unique(hour))

#for(i in 1:13){
#  ti <- as.character(unique(hour))[i]
#  tmp <- des2[unlist(LIGHT[i]),]
#  write.csv(tmp, 
#            file=sprintf("output/TCC/%s_LIGHT_%s.csv",exec.date, ti))
#  tmp <- des2[unlist(TH[i]),]
#  write.csv(tmp, 
#            file=sprintf("output/TCC/%s_TH_%s.csv",exec.date, ti))
#  tmp <- des2[unlist(LTH[i]),]
#  write.csv(tmp,
#            file=sprintf("output/TCC/%s_LTH_%s.csv",exec.date, ti))
#  tmp <- des2[unlist(UNREP[i]),]
#  write.csv(tmp,
#            file=sprintf("output/TCC/%s_UNREP_%s.csv",exec.date, ti))
#}

### linegraph 
LIGHT_n <- rep(NA, 13)
for(i in 1:13){
  LIGHT_n[i] <- length(LIGHT[[i]])
}
TH_n <- rep(NA, 13)
for(i in 1:13){
  TH_n[i] <- length(TH[[i]])
}
LTH_n <- rep(NA, 13)
for(i in 1:13){
  LTH_n[i] <- length(LTH[[i]])
}
UNREP_n <- rep(NA, 13)
for(i in 1:13){
  UNREP_n[i] <- length(UNREP[[i]])
}

# Fig. 2c ----
d <- data.frame(LIGHT = LIGHT_n, TH = TH_n, LTH = LTH_n, UNREP = UNREP_n,
                hour = unique(hour)) %>% 
  gather(., "LIGHT","TH","LTH","UNREP",key = "type", value = "n")
d$type <- factor(d$type, levels = c("LIGHT","TH", "LTH", "UNREP"))

fig2c <- ggplot(d, aes(x=hour, y= n , group=type,
                   color=type, fill=type, shape=type))+
  geom_line(size = 0.25)+
  geom_point(color="black", size=1 ,stroke=0.25)+
  scale_color_manual(values = color_5[1:4])+
  scale_fill_manual(values = color_5[1:4])+
  scale_shape_manual(values = c(21:24))+
  scale_x_discrete(labels = times)+
  labs(x="Time (hours)", y="The number of genes")+
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
        legend.justification = "center")+
  guides(color=guide_legend(ncol=2,byrow = T))
fig2c

#ggsave(sprintf("%s/%s_Fig2c.pdf", dir.output, exec.date), 
#       width=55, height = 40, units = "mm")

### GO significant ----
rpm2.ave <- rpm.ave[raw_10,]
use.genelist <- rownames(rpm2.ave)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))
term_na <- term[term != "NA"]
time <- unique(at$hour)
condition <- read.delim(sprintf("%s/condition.txt", dir.input))

dir.create(sprintf("%s/GO", dir.output))
for(i in 1:10){
  dir.create(sprintf("%s/GO/%s_%s", dir.output, condition[i,1],condition[i,2]))
}

for (j in 1:10){
  # Field vs SmartGC
  c1 <- condition[j,1]
  c2 <- condition[j,2]
  
  for (i in 1:13){
    ID <- time[i]
    load(sprintf("%s/TCC/210703_TCC_0.05_%s", dir.output, ID))
    a <-  TCC_table[,j]==T
    gl <- rownames(TCC_table)[a]
    result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
    
    fn <- sprintf("%s/GO/%s_%s/%s_GO_%s_%s_%s.csv", dir.output, c1, c2, exec.date, c1, c2, ID)
    
    adp <- p.adjust(result[,"p.value"],method = "BH")
    
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
  for (i in 1:13){
    
    ID <- time[i]
    a <- read.csv(sprintf("%s/GO/%s_%s/210703_GO_%s_%s_%s.csv", dir.output, c1, c2, c1, c2, ID))
    if (is.na(a[1,1])){
      next
    }else{
      hour <- rep(ID, nrow(a))
      a <- cbind(hour, a)
      tmp <- rbind(tmp, a)
    }
    #cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  write.csv(tmp, file= sprintf("%s/GO/%s_GO_%s_%s.csv", dir.output, exec.date, c1, c2), row.names = F)
}

### GO heatmap ----
rpm2.ave <- rpm.ave[raw_10,]
use.genelist <- rownames(rpm2.ave)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))
term_na <- term[term != "NA"]
time <- unique(at$hour)
condition <- read.delim(sprintf("%s/condition.txt", dir.input))

dir.create(sprintf("%s/GO/heatmap", dir.output))
for(i in 1:10){
  dir.create(sprintf("%s/GO/heatmap/%s_%s", dir.output, condition[i,1],condition[i,2]))
}

for (j in 1:10){
  # Field vs SmartGC
  c1 <- condition[j,1]
  c2 <- condition[j,2]
  
  for (i in 1:13){
    ID <- time[i]
    load(sprintf("%s/TCC/210703_TCC_0.05_%s", dir.output, ID))
    a <-  TCC_table[,j]==T
    gl <- rownames(TCC_table)[a]
    result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
    adp <- p.adjust(result[,"p.value"], method = "BH")
    
    fn <- sprintf("%s/GO/heatmap/%s_%s/%s_GO_%s_%s_%s.csv", dir.output, c1, c2, exec.date, c1, c2, ID)
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
  for (i in 1:13){
    
    ID <- time[i]
    a <- read.csv(sprintf("%s/GO/heatmap/%s_%s/210703_GO_%s_%s_%s.csv", dir.output, c1, c2, c1, c2, ID))
    if (is.na(a[1,1])){
      next
    }else{
      hour <- rep(ID, nrow(a))
      a <- cbind(hour, a)
      tmp <- rbind(tmp, a)
    }
    #cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  write.csv(tmp, file= sprintf("%s/GO/heatmap/%s_GO_%s_%s.csv", dir.output, exec.date, c1, c2), row.names = F)
}

### heatmap ----
a <- read.csv("output/GO/heatmap/210703_GO_FIELD_FL_FTH.csv")
a1 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d1 <- as.matrix(b[,2:14])

a1$hour <- factor(a1$hour, levels=time)
b1 <- a1 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/GO/heatmap/210703_GO_FIELD_CL_CTH.csv")
a2 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d2 <- as.matrix(b[,2:14])

a2$hour <- factor(a2$hour, levels=time)
b2 <- a2 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/GO/heatmap/210703_GO_FIELD_FL_CTH.csv")
a3 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d3 <- as.matrix(b[,2:14])

a3$hour <- factor(a3$hour, levels=time)
b3 <- a3 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/GO/heatmap/210703_GO_FIELD_CL_FTH.csv")
a4 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d4 <- as.matrix(b[,2:14])

a4$hour <- factor(a4$hour, levels=time)
b4 <- a4 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d <- cbind(d1,d2,d3,d4)

pvalue <- function(dat){
  (min(dat) < log10(0.05)) == T
}

dat <- d[apply(d,1,pvalue),]
dat[dat < -5] <- -5
tmp <- heatmap(x = dat, Colv = NA)
dat <- dat[tmp$rowInd,]

ID <- as.character(b$ID)[apply(d,1,pvalue)]
GO_Term <- ng.GetGOTerms(ID)[tmp$rowInd]

# Table S5 ----
dID <- rev(ID[tmp$rowInd])
tex <- matrix(NA, nrow=length(dID), ncol=9)
rownames(tex) <- dID
colnames(tex) <- c("ID","GO term","condition","time","Adjusted p-value","A & B","A","B","U")
tex[,1] <- dID
tex[,2] <-rev(GO_Term)

cond <- rep(c("FL_FTH","CL_CTH","FL_CTH","CL_FTH"),each=13)
tim <- rep(time, 4)

rownames(b1) <- b1$ID
rownames(b2) <- b2$ID
rownames(b3) <- b3$ID
rownames(b4) <- b4$ID

aaa <- cbind(b1[dID,2:14], b2[dID,2:14],b3[dID,2:14],b4[dID,2:14])
bb <- apply(aaa,1,min)

for(i in 1:length(dID)){
  
  cc <- cond[aaa[i,] == bb[i]]
  dd <- as.character(tim)[aaa[i,] == bb[i]]
  ee <- bb[i]
  
  ff <- read.csv(sprintf("output/GO/heatmap/210703_GO_FIELD_%s.csv",cc))
  
  gg <- as.numeric(ff[(ff$hour==dd & ff$ID == dID[i]),5:8])
  
  tex[i,3] <- cc
  tex[i,4] <- dd
  tex[i,5] <- sprintf("%.2e",as.numeric(ee))
  tex[i,6:9] <- gg
}

write.csv(tex, file= sprintf("%s/%s_TableS5.csv", dir.output, exec.date),
          row.names = F)

# Fig. 4a ----
dat[dat >= log10(0.05)] <- 0

png(sprintf("%s/%s_Fig4a.png", dir.output, exec.date),
    width = 40.75, height = 18.847, units = "mm", res=1000)

par(mar=c(0,0,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,52),xaxt="n",
     at = seq(0, 1, length.out = 52))
#axis(side = 2, labels = GO_Term, las = 1,cex.axis=0.5,
#     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.6)
#par(oma=c(0,0,0,0)6
#image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0))
#par(oma = c(0,0,3,0))
#title("")

dev.off()

# fig
pdf(sprintf("%s/%s_Fig4a_2.pdf", dir.output, exec.date),
    width = 7, height = 12)

par(oma=c(0,10,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,52),xaxt="n",
     at = seq(0, 1, length.out = 52))
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
rpm2.ave <- rpm.ave[raw_10,]
use.genelist <- rownames(rpm2.ave)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))
term_na <- term[term != "NA"]
time <- unique(at$hour)
condition <- read.delim(sprintf("%s/condition.txt", dir.input))

dir.create("output/venndiagram")
dir.create("output/venndiagram/LIGHT")
dir.create("output/venndiagram/TH")
dir.create("output/venndiagram/LTH")
dir.create("output/venndiagram/UNREP")

# LIGHT
for (i in 1:13){
  ID <- time[i]
  gl <- LIGHT[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram/LIGHT/%s_GO_LIGHT_%s.csv", dir.output, exec.date, ID)
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/LIGHT/210703_GO_LIGHT_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/%s_GO_LIGHT.csv", dir.output, exec.date), row.names = F)

# TH
for (i in 1:13){
  ID <- time[i]
  gl <- TH[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram/TH/%s_GO_TH_%s.csv", dir.output, exec.date, ID)
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/TH/210703_GO_TH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/%s_GO_TH.csv", dir.output, exec.date), row.names = F)

# LTH
for (i in 1:13){
  ID <- time[i]
  gl <- LTH[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram/LTH/%s_GO_LTH_%s.csv", dir.output, exec.date, ID)
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/LTH/210703_GO_LTH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/%s_GO_LTH.csv", dir.output, exec.date), row.names = F)

# UNREP
for (i in 1:13){
  ID <- time[i]
  gl <- UNREP[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram/UNREP/%s_GO_UNREP_%s.csv", dir.output, exec.date, ID)
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/UNREP/210703_GO_UNREP_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/%s_GO_UNREP.csv", dir.output, exec.date), row.names = F)

### GO DEG venn diagram heatmap----
rpm2.ave <- rpm.ave[raw_10,]
use.genelist <- rownames(rpm2.ave)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))
term_na <- term[term != "NA"]
time <- unique(at$hour)
condition <- read.delim(sprintf("%s/condition.txt", dir.input))

dir.create("output/venndiagram/heatmap/")
dir.create("output/venndiagram/heatmap/LIGHT")
dir.create("output/venndiagram/heatmap/TH")
dir.create("output/venndiagram/heatmap/LTH")
dir.create("output/venndiagram/heatmap/UNREP")

# LIGHT
for (i in 1:13){
  ID <- time[i]
  gl <- LIGHT[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram/heatmap/LIGHT/%s_GO_LIGHT_%s.csv", dir.output, exec.date, ID)
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/heatmap/LIGHT/210703_GO_LIGHT_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/heatmap/%s_GO_LIGHT.csv", dir.output, exec.date), row.names = F)

# TH
for (i in 1:13){
  ID <- time[i]
  gl <- TH[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram/heatmap/TH/%s_GO_TH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"],method = "BH")
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/heatmap/TH/210703_GO_TH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/heatmap/%s_GO_TH.csv", dir.output, exec.date), row.names = F)

# LTH
for (i in 1:13){
  ID <- time[i]
  gl <- LTH[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram/heatmap/LTH/%s_GO_LTH_%s.csv", dir.output, exec.date, ID)
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/heatmap/LTH/210703_GO_LTH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/heatmap/%s_GO_LTH.csv", dir.output, exec.date), row.names = F)

# UNREP
for (i in 1:13){
  ID <- time[i]
  gl <- UNREP[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("%s/venndiagram/heatmap/UNREP/%s_GO_UNREP_%s.csv", dir.output, exec.date, ID)
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/heatmap/UNREP/210703_GO_UNREP_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/heatmap/%s_GO_UNREP.csv", dir.output, exec.date), row.names = F)

### heatmap ----
a <- read.csv("output/venndiagram/heatmap/210703_GO_LIGHT.csv")
a1 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d1 <- as.matrix(b[,2:14])

a1$hour <- factor(a1$hour, levels=time)
b1 <- a1 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/venndiagram/heatmap/210703_GO_TH.csv")
a2 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d2 <- as.matrix(b[,2:14])

a2$hour <- factor(a2$hour, levels=time)
b2 <- a2 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/venndiagram/heatmap/210703_GO_LTH.csv")
a3 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d3 <- as.matrix(b[,2:14])

a3$hour <- factor(a3$hour, levels=time)
b3 <- a3 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/venndiagram/heatmap/210703_GO_UNREP.csv")
a4 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d4 <- as.matrix(b[,2:14])

a4$hour <- factor(a4$hour, levels=time)
b4 <- a4 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d <- cbind(d1,d2,d3,d4)

pvalue <- function(dat){
  (min(dat) < log10(0.05)) == T
}

dat <- d[apply(d,1,pvalue),]
dat[dat < -5] <- -5
tmp <- heatmap(x = dat, Colv = NA)
dat <- dat[tmp$rowInd,]

ID <- as.character(b$ID)[apply(d,1,pvalue)]
GO_Term <- ng.GetGOTerms(ID)[tmp$rowInd]

# Table S6 ----
dID <- rev(ID[tmp$rowInd])
tex <- matrix(NA, nrow=length(dID), ncol=9)
rownames(tex) <- dID
colnames(tex) <- c("ID","GO term","condition","time","Adjusted p-value",  "A & B", "A",	"B", "U")
tex[,1] <- dID
tex[,2] <- rev(GO_Term)

cond <- rep(c("LIGHT","TH","LTH","UNREP"),each=13)
tim <- rep(time, 4)

rownames(b1) <- b1$ID
rownames(b2) <- b2$ID
rownames(b3) <- b3$ID
rownames(b4) <- b4$ID
aaa <- cbind(b1[dID,2:14], b2[dID,2:14],b3[dID,2:14],b4[dID,2:14])
bb <- apply(aaa,1,min)

for(i in 1:length(dID)){
  
  cc <- cond[aaa[i,] == bb[i]]
  dd <- as.character(tim)[aaa[i,] == bb[i]]
  ee <- bb[i]
  
  ff <- read.csv(sprintf("output/venndiagram/heatmap/210703_GO_%s.csv",cc))
  
  gg <- as.numeric(ff[(ff$hour==dd & ff$ID == dID[i]),5:8])
  
  tex[i,3] <- cc
  tex[i,4] <- dd
  tex[i,5] <- sprintf("%.2e",as.numeric(ee))
  tex[i,6:9] <- gg
}

write.csv(tex, file= sprintf("%s/%s_TableS6.csv", dir.output, exec.date),
          row.names = F)

# Fig. 4b ----
dat[dat >= log10(0.05)] <- 0

png(sprintf("%s/%s_Fig4b.png", dir.output, exec.date),
    width = 40.75, height = 18.847, units = "mm", res=1000)

par(mar=c(0,0,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,52),xaxt="n",
     at = seq(0, 1, length.out = 52))
#axis(side = 2, labels = GO_Term, las = 1,cex.axis=0.5,
#     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.6)
#par(oma=c(0,0,0,0)6
#image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0))
#par(oma = c(0,0,3,0))
#title("")

dev.off()

# fig
pdf(sprintf("%s/%s_Fig4b_2.pdf", dir.output, exec.date),
    width = 7, height = 12)

par(oma=c(0,10,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,52),xaxt="n",
     at = seq(0, 1, length.out = 52))
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
rpm2.ave <- rpm.ave[raw_10,]
use.genelist <- rownames(rpm2.ave)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]
time <- unique(at$hour)
condition <- read.delim(sprintf("%s/condition.txt", dir.input))
dir.create(sprintf("%s/KEGG",dir.output))
for(i in 1:10){
  dir.create(sprintf("%s/KEGG/%s_%s", dir.output, condition[i,1],condition[i,2]))
}

for (j in 1:10){
  # Field vs SmartGC
  c1 <- condition[j,1]
  c2 <- condition[j,2]
  
  for (i in 1:13){
    ID <- time[i]
    load(sprintf("%s/TCC/210703_TCC_0.05_%s", dir.output, ID))
    a <-  TCC_table[,j]==T
    gl <- rownames(TCC_table)[a]
    result <- kegg.mft(cgt=kegg, gn.test=gl)
    
    fn <- sprintf("%s/KEGG/%s_%s/%s_KEGG_%s_%s_%s.csv", dir.output, c1, c2, exec.date, c1, c2, ID)
    
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
  for (i in 1:13){
    
    ID <- time[i]
    a <- read.csv(sprintf("%s/KEGG/%s_%s/210703_KEGG_%s_%s_%s.csv", dir.output, c1, c2, c1, c2, ID))
    if (is.na(a[1,1])){
      next
    }else{
      hour <- rep(ID, nrow(a))
      a <- cbind(hour, a)
      tmp <- rbind(tmp, a)
    }
    #cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  write.csv(tmp, file= sprintf("%s/KEGG/%s_KEGG_%s_%s.csv", dir.output, exec.date, c1, c2), row.names = F)
}

### KEGG heatmap ----
rpm2.ave <- rpm.ave[raw_10,]
use.genelist <- rownames(rpm2.ave)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]
time <- unique(at$hour)
condition <- read.delim(sprintf("%s/condition.txt", dir.input))
dir.create(sprintf("%s/KEGG/heatmap", dir.output))
for(i in 1:10){
  dir.create(sprintf("%s/KEGG/heatmap/%s_%s", dir.output, condition[i,1],condition[i,2]))
}

for (j in 1:10){
  # Field vs SmartGC
  c1 <- condition[j,1]
  c2 <- condition[j,2]
  
  for (i in 1:13){
    ID <- time[i]
    load(sprintf("%s/TCC/210703_TCC_0.05_%s", dir.output, ID))
    a <-  TCC_table[,j]==T
    gl <- rownames(TCC_table)[a]
    result <- kegg.mft(cgt=kegg, gn.test=gl)
    
    fn <- sprintf("%s/KEGG/heatmap/%s_%s/%s_KEGG_%s_%s_%s.csv", dir.output, c1, c2, exec.date, c1, c2, ID)
    
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
  for (i in 1:13){
    
    ID <- time[i]
    a <- read.csv(sprintf("%s/KEGG/heatmap/%s_%s/210703_KEGG_%s_%s_%s.csv", dir.output, c1, c2, c1, c2, ID))
    if (is.na(a[1,1])){
      next
    }else{
      hour <- rep(ID, nrow(a))
      a <- cbind(hour, a)
      tmp <- rbind(tmp, a)
    }
    #cat(sprintf("%s / %s\n", i, Sys.time()))
  }
  write.csv(tmp, file= sprintf("%s/KEGG/heatmap/%s_KEGG_%s_%s.csv", dir.output, exec.date, c1, c2), row.names = F)
}

### heatmap ----
a <- read.csv("output/KEGG/heatmap/210703_KEGG_FIELD_FL_FTH.csv")
a1 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d1 <- as.matrix(b[,2:14])

a1$hour <- factor(a1$hour, levels=time)
b1 <- a1 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/KEGG/heatmap/210703_KEGG_FIELD_CL_CTH.csv")
a2 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d2 <- as.matrix(b[,2:14])

a2$hour <- factor(a2$hour, levels=time)
b2 <- a2 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/KEGG/heatmap/210703_KEGG_FIELD_FL_CTH.csv")
a3 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d3 <- as.matrix(b[,2:14])

a3$hour <- factor(a3$hour, levels=time)
b3 <- a3 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/KEGG/heatmap/210703_KEGG_FIELD_CL_FTH.csv")
a4 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d4 <- as.matrix(b[,2:14])

a4$hour <- factor(a4$hour, levels=time)
b4 <- a4 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)


d <- cbind(d1,d2,d3,d4)

pvalue <- function(dat){
  (min(dat) < log10(0.05)) == T
}

dat <- d[apply(d,1,pvalue),] # at least 1 time is significant
dat[dat < -5] <- -5
tmp <- heatmap(x = dat, Colv = NA)
dat <- dat[tmp$rowInd,]

ID <- as.character(b$ID)[apply(d,1,pvalue)]
KEGG_Term <- kegg_description[ID][tmp$rowInd]

# Table S7 ----
dID <- rev(ID[tmp$rowInd])
tex <- matrix(NA, nrow=length(dID), ncol=9)
rownames(tex) <- dID
colnames(tex) <- c("ID","KEGG term","condition","time","Adjusted p-value", "A & B", "A",	"B", "U")
tex[,1] <- dID
tex[,2] <- rev(KEGG_Term)

cond <- rep(c("FL_FTH","CL_CTH","FL_CTH","CL_FTH"),each=13)
tim <- rep(time, 4)

rownames(b1) <- b1$ID
rownames(b2) <- b2$ID
rownames(b3) <- b3$ID
rownames(b4) <- b4$ID

aaa <- cbind(b1[dID,2:14], b2[dID,2:14],b3[dID,2:14],b4[dID,2:14])
bb <- apply(aaa,1,min)

for(i in 1:length(dID)){
  
  cc <- cond[aaa[i,] == bb[i]]
  dd <- as.character(tim)[aaa[i,] == bb[i]]
  ee <- bb[i]
  
  ff <- read.csv(sprintf("output/KEGG/heatmap/210703_KEGG_FIELD_%s.csv",cc))
  
  gg <- as.numeric(ff[(ff$hour==dd & ff$ID == dID[i]),5:8])
  
  tex[i,3] <- cc
  tex[i,4] <- dd
  tex[i,5] <- sprintf("%.2e",as.numeric(ee))
  tex[i,6:9] <- gg
}

write.csv(tex, file= sprintf("%s/%s_TableS7.csv", dir.output, exec.date),
          row.names = F)

# Fig. 4c ----
dat[dat >= log10(0.05)] <- 0

png(sprintf("%s/%s_Fig4c.png",dir.output,exec.date),
    width = 94.989, height = 46.698, units= "mm", res=1000)

par(mar=c(0,0,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,52),xaxt="n",
     at = seq(0, 1, length.out = 52))
#axis(side = 4, labels = NA, las = 1,cex.axis=1,
#     at = seq(0, 1, length.out = length(KEGG_Term)))
box(lwd=0.6)
#par(oma=c(0,0,0,0))
#image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0))
#par(oma = c(0,0,3,0))
#title("")

dev.off()

#pdf
pdf(sprintf("%s/%s_Fig4c_2.pdf",dir.output,exec.date),
    width = 10, height = 7)

par(oma=c(0,0,0,10))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,52),xaxt="n",
     at = seq(0, 1, length.out = 52))
axis(side = 4, labels = KEGG_Term, las = 1,cex.axis=0.5,
     at = seq(0, 1, length.out = length(KEGG_Term)))
box(lwd=0.5)
par(oma=c(0,0,0,0))
image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0),
           horizontal = T)
par(oma = c(0,0,3,0))
title("")

dev.off()

### KEGG DEG venn significant ----
rpm2.ave <- rpm.ave[raw_10,]
use.genelist <- rownames(rpm2.ave)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]
time <- unique(at$hour)
condition <- read.delim(sprintf("%s/condition.txt", dir.input))

# LIGHT
for (i in 1:13){
  ID <- time[i]
  gl <- LIGHT[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram/LIGHT/%s_KEGG_LIGHT_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"],method = "BH")
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/LIGHT/210703_KEGG_LIGHT_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/%s_KEGG_LIGHT.csv", dir.output, exec.date), row.names = F)

# TH
for (i in 1:13){
  ID <- time[i]
  gl <- TH[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram/TH/%s_KEGG_TH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"],method = "BH")
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/TH/210703_KEGG_TH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/%s_KEGG_TH.csv", dir.output, exec.date), row.names = F)

# LTH
for (i in 1:13){
  ID <- time[i]
  gl <- LTH[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram/LTH/%s_KEGG_LTH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"],method = "BH")
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/LTH/210703_KEGG_LTH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/%s_KEGG_LTH.csv", dir.output, exec.date), row.names = F)

# UNREP
for (i in 1:13){
  ID <- time[i]
  gl <- UNREP[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram/UNREP/%s_KEGG_UNREP_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"],method = "BH")
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/UNREP/210703_KEGG_UNREP_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/%s_KEGG_UNREP.csv", dir.output, exec.date), row.names = F)

### KEGG heatmap DEG venn diagram ----
rpm2.ave <- rpm.ave[raw_10,]
use.genelist <- rownames(rpm2.ave)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]
time <- unique(at$hour)
condition <- read.delim(sprintf("%s/condition.txt", dir.input))

# LIGHT
for (i in 1:13){
  ID <- time[i]
  gl <- LIGHT[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram/heatmap/LIGHT/%s_KEGG_LIGHT_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"],method = "BH")
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/heatmap/LIGHT/210703_KEGG_LIGHT_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/heatmap/%s_KEGG_LIGHT.csv", dir.output, exec.date), row.names = F)

# TH
for (i in 1:13){
  ID <- time[i]
  gl <- TH[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram/heatmap/TH/%s_KEGG_TH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"],method = "BH")
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/heatmap/TH/210703_KEGG_TH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/heatmap/%s_KEGG_TH.csv", dir.output, exec.date), row.names = F)

# LTH
for (i in 1:13){
  ID <- time[i]
  gl <- LTH[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram/heatmap/LTH/%s_KEGG_LTH_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"],method = "BH")
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/heatmap/LTH/210703_KEGG_LTH_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/heatmap/%s_KEGG_LTH.csv", dir.output, exec.date), row.names = F)

# UNREP
for (i in 1:13){
  ID <- time[i]
  gl <- UNREP[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("%s/venndiagram/heatmap/UNREP/%s_KEGG_UNREP_%s.csv", dir.output, exec.date, ID)
  
  adp <- p.adjust(result[,"p.value"],method = "BH")
  
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
for (i in 1:13){
  
  ID <- time[i]
  a <- read.csv(sprintf("%s/venndiagram/heatmap/UNREP/210703_KEGG_UNREP_%s.csv", dir.output, ID))
  if (is.na(a[1,1])){
    next
  }else{
    hour <- rep(ID, nrow(a))
    a <- cbind(hour, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}
write.csv(tmp, file= sprintf("%s/venndiagram/heatmap/%s_KEGG_UNREP.csv", dir.output, exec.date), row.names = F)

### heatmap ----
a <- read.csv("output/venndiagram/heatmap/210703_KEGG_LIGHT.csv")
a1 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d1 <- as.matrix(b[,2:14])

a1$hour <- factor(a1$hour, levels=time)
b1 <- a1 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/venndiagram/heatmap/210703_KEGG_TH.csv")
a2 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d2 <- as.matrix(b[,2:14])

a2$hour <- factor(a2$hour, levels=time)
b2 <- a2 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/venndiagram/heatmap/210703_KEGG_LTH.csv")
a3 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d3 <- as.matrix(b[,2:14])

a3$hour <- factor(a3$hour, levels=time)
b3 <- a3 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

a <- read.csv("output/venndiagram/heatmap/210703_KEGG_UNREP.csv")
a4 <- a
a$hour <- factor(a$hour, levels=time)
a$Adjusted.P.value <- log10(a$Adjusted.P.value)
b <- a %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d4 <- as.matrix(b[,2:14])

a4$hour <- factor(a4$hour, levels=time)
b4 <- a4 %>% dplyr::select(hour,Adjusted.P.value,ID) %>% 
  spread(hour, Adjusted.P.value)

d <- cbind(d1,d2,d3,d4)

pvalue <- function(dat){
  (min(dat) < log10(0.05)) == T
}

dat <- d[apply(d,1,pvalue),] # at least 1 time is significant
dat[dat < -5] <- -5
tmp <- heatmap(x = dat, Colv = NA)
dat <- dat[tmp$rowInd,]

ID <- as.character(b$ID)[apply(d,1,pvalue)]
KEGG_Term <- kegg_description[ID][tmp$rowInd]

# Table S8 ----
dID <- rev(ID[tmp$rowInd])
tex <- matrix(NA, nrow=length(dID), ncol=9)
rownames(tex) <- dID
colnames(tex) <- c("ID","KEGG term","condition","time","Adjusted p-value", "A & B", "A",	"B", "U")
tex[,1] <- dID
tex[,2] <- rev(KEGG_Term)

cond <- rep(c("LIGHT","TH","LTH","UNREP"),each=13)
tim <- rep(time, 4)

rownames(b1) <- b1$ID
rownames(b2) <- b2$ID
rownames(b3) <- b3$ID
rownames(b4) <- b4$ID
aaa <- cbind(b1[dID,2:14], b2[dID,2:14],b3[dID,2:14],b4[dID,2:14])
bb <- apply(aaa,1,min)


for(i in 1:length(dID)){
  
  cc <- cond[aaa[i,] == bb[i]]
  dd <- as.character(tim)[aaa[i,] == bb[i]]
  ee <- bb[i]
  
  ff <- read.csv(sprintf("output/venndiagram/heatmap/210703_KEGG_%s.csv",cc))
  
  gg <- as.numeric(ff[(ff$hour==dd & ff$ID == dID[i]),5:8])
  
  tex[i,3] <- cc
  tex[i,4] <- dd
  tex[i,5] <- sprintf("%.2e",as.numeric(ee))
  tex[i,6:9] <- gg
}

write.csv(tex, file= sprintf("%s/%s_TableS8.csv", dir.output, exec.date),
          row.names = F)

# Fig. 4d ----
dat[dat >= log10(0.05)] <- 0

png(sprintf("%s/%s_Fig4d.png",dir.output,exec.date),
    width = 94.989, height = 31.132, units = "mm", res=1000)

par(mar=c(0,0,0,0))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,52),xaxt="n",
     at = seq(0, 1, length.out = 52))
#axis(side = 4, labels = NA, las = 1,cex.axis=0.5,
#     at = seq(0, 1, length.out = length(ID)))
box(lwd=0.5)
#par(oma=c(0,0,0,0))
#image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0))
#par(oma = c(0,0,3,0))
#title("")

dev.off()

#pdf
pdf(sprintf("%s/%s_Fig4d_2.pdf",dir.output,exec.date),
    width = 10, height = 7)

par(oma=c(0,0,0,10))
image(t(dat), axes = FALSE, col=rev(ng.po.colors(64)), zlim=c(-5,0))
axis(side = 1, labels = rep(NA,52),xaxt="n",
     at = seq(0, 1, length.out = 52))
axis(side = 4, labels = KEGG_Term, las = 1,cex.axis=0.5,
     at = seq(0, 1, length.out = length(ID)))
box()
par(oma=c(0,0,0,0))
image.plot(dat, col=rev(ng.po.colors(64)), legend.only = T, zlim=c(-5,0),
           horizontal = T)
par(oma = c(0,0,3,0))
title("")

dev.off()

### geneplot ----
# function
#geneplot <- function(ID,lab = as.character(des2[ID,]$Description)){
#  con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
#  time <- unique(at$hour)
#  hour <- factor(at$hour, levels = time)
#  ylab <- expression(paste(log[2],"(rpm)"))
#  
#  if(length(ID)==1){
#    d <- data.frame(log2rpm[ID,]) %>% 
#     mutate(con, hour) %>% 
#      gather(gene,expression,-con, -hour) 
#    
#    d_mean_sd <- d %>%
#      group_by(con, hour) %>%
#      summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))
#    
#    names(lab) <- ID
#    errors <- aes(ymax = mean + sd, ymin = mean - sd)
#    
#    p <- ggplot(d_mean_sd,
#                aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
#      geom_line(size=0.25)+
#      geom_errorbar(errors, width = 0.2)+
#      geom_point(color="black",size = 1)+
#      scale_color_manual(values = color5)+
#      scale_fill_manual(values = color5)+
#     scale_shape_manual(values = c(21:25))+
#      ggtitle(sprintf("%s\n%s", ID, lab))+
#      labs(x="Time(hour)", y=ylab)+
#      theme_cowplot(font_size=6, line_size = 0.25)+
#      theme(plot.title = element_text(size=6, hjust = 0.5),
#            axis.title.x = element_text(size=6),
#            axis.title.y = element_text(size=6),
#            axis.text.x = element_text(size=6),
#            axis.text.y = element_text(size=6),
#            legend.title = element_blank(),
#            legend.text = element_text(size=6),
#            legend.position = "top",
#            legend.justification = "center")
#    print(p)
#  }else{
#    d <- data.frame(t(log2rpm[ID,])) %>% 
#      mutate(con, hour) %>% 
#      gather(gene,expression,-con, -hour) 
#    d$gene <- factor(d$gene, levels = ID)
#    
#    d_mean_sd <- d %>%
#      group_by(con, hour, gene) %>%
#      summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))
#    
#    names(lab) <- ID
#    errors <- aes(ymax = mean + sd, ymin = mean - sd)
#    
#    p <- d_mean_sd %>%
#      group_by(gene) %>% 
#      nest() %>% 
#      mutate(plot = map2(data, gene, 
#                         ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
#                           geom_line(size=0.25)+
#                           geom_errorbar(errors, width = 0.2)+
#                           geom_point(color="black",size = 1)+
#                           scale_color_manual(values = color5)+
#                           scale_fill_manual(values = color5)+
#                           scale_shape_manual(values = c(21:25))+
#                           ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
#                           #labs(x="Time(hour)", y=ylab)+
#                           theme_cowplot(font_size=6, line_size = 0.25)+
#                           theme(plot.title = element_text(size=6, hjust = 0.5),
#                                 axis.title.x = element_text(size=6),
#                                 axis.title.y = element_text(size=6),
#                                 axis.text.x = element_text(size=6),
#                                 axis.text.y = element_text(size=6),
#                                 legend.title = element_blank(),
#                                 legend.text = element_text(size=6),
#                                 legend.position = "top",
#                                 legend.justification = "center")))
#    p$plot
#  }
#}

# Fig. 2e-g ----
#a <- read.csv("output/TCC/200311_LIGHT_7.csv")
#b <- read.csv("output/TCC_204/210506_LIGHT_7.csv")
#c <- intersect(a$X, b$X)

#pdf("210506_LIGHT.pdf")
#geneplot(c)
#dev.off()

#a <- read.csv("output/TCC/200311_TH_7.csv")
#b <- read.csv("output/TCC_204/210506_TH_7.csv")
#c <- intersect(a$X, b$X)

#pdf("210506_TH.pdf")
#geneplot(c)
#dev.off()

#a <- read.csv("output/TCC/200311_LTH_7.csv")
#b <- read.csv("output/TCC_204/210506_LTH_7.csv")
#c <- intersect(a$X, b$X)

#pdf("210506_LTH.pdf")
#geneplot(c)
#dev.off()

ID <- c("Os03g0169600","Os10g0496900","Os04g0444900")
lab <- c("Dof12","PORB","PHD21")
names(lab) <- ID

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

ylab <- expression(paste(log[2],"(rpm)"))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.1, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time(hour)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size = 6),
                             axis.title.y =  element_text(size = 6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

# Experiment 2
con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- factor(time, levels=time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

d <- as_tibble(t(log2rpm_204[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 
d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

p2 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.1, size = 0.25)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color4)+
                       scale_fill_manual(values = color4)+
                       scale_shape_manual(values = c(22:25))+
                       scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
                       labs(x="Time(hour)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour="black", size=0.5)+
                       theme(plot.title = element_blank(),
                             axis.title.x =  element_text(size = 6),
                             axis.title.y =  element_text(size = 6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- g_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p2$plot[[1]]+ theme(legend.position="none")
g5 <- p2$plot[[2]]+ theme(legend.position="none")
g6 <- p2$plot[[3]]+ theme(legend.position="none")

p <- plot_grid(g1,g2,g3,g4,g5,g6, nrow = 2, align="hv")

ggsave(p, file=sprintf("%s/%s_Fig2efg.pdf",dir.output, exec.date),
       width = 120, height = 60, units="mm")

# Fig. 3 ----
a <- read.csv("input/191025_clockgene.csv")
ID <- as.character(a$ID)[c(6,8,35,34,48)]
lab <- as.character(a$name)[c(6,8,35,34,48)]
names(lab) <- ID
ylab <- expression(paste(log[2],"(rpm)"))

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

# Experiment 2
con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- factor(time, levels=time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

d <- as_tibble(t(log2rpm_204[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 
d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

p2 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color4)+
                       scale_fill_manual(values = color4)+
                       scale_shape_manual(values = c(22:25))+
                       scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour="black", size=0.5)+
                       theme(plot.title = element_blank(),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p1$plot[[4]]+ theme(legend.position="none")
g5 <- p1$plot[[5]]+ theme(legend.position="none")

h1 <- p2$plot[[1]]+ theme(legend.position="none")
h2 <- p2$plot[[2]]+ theme(legend.position="none")
h3 <- p2$plot[[3]]+ theme(legend.position="none")
h4 <- p2$plot[[4]]+ theme(legend.position="none")
h5 <- p2$plot[[5]]+ theme(legend.position="none")

p <- plot_grid(mylegend, 
               plot_grid(g1,g2,NULL,h1,h2,NULL,g3,g4,g5,h3,h4,h5,nrow = 4, align="hv"),
               ncol=1, rel_heights = c(1,15))

ggsave(p, file=sprintf("%s/%s_Fig3.pdf",dir.output, exec.date),
       width = 160, height = 160, units="mm")

# Fig. S8 and TableS2 ----
### molecular timetable method (Ueda et al., PNAS, 2004; Higashi et al., Frontiers in Plant Science, 2016)
# label
dir.create("output/MTB")

condition <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- factor(unique(at$hour), levels=unique(at$hour))
hour <- factor(at$hour, levels = time)

# target value
ZT <- rep(rep(c(0,2,4,6,8,10,12,14,16,18,20,22,24), each=4), 5)
names(ZT) = colnames(log2rpm)
ZT.field <- ZT[1:52]

# MTB by FIELD data
set.seed(1234)
log2rpm.field <- t(log2rpm[raw_10,1:52])
log2rpm.field.normalized <- apply(log2rpm.field,2,scale)
rownames(log2rpm.field) <- colnames(log2rpm)[1:52]
rownames(log2rpm.field.normalized) <- colnames(log2rpm)[1:52]

# use all sample
fit.gnl = list()
for (gn in colnames(log2rpm.field)){
  cat(gn)
  cat("\n")
  tmp = NULL
  ex = log2rpm.field[,gn]
  tv = ZT.field[names(ex)]
  co = rep(NA, 1440)
  names(co) = 0:1339
  for (b in 1:1440){
    predict = sqrt(2)*cos(2*pi*(tv-(b-1)/60)/24)
    co[b] = cor(ex,predict) 
  }
  tmp = list(1/60*as.numeric(names(which(co == max(co,na.rm=T)))))
  names(tmp) = gn
  fit.gnl = c(fit.gnl,tmp)
}
fit.gnl = fit.gnl[!is.na(fit.gnl)]

# Caluculate correlation between observed values and fitted values
cor.gnl = list()
for (gn in names(fit.gnl)) {
  predict = sqrt(2)*cos(2*pi*(ZT[rownames(log2rpm.field)]-fit.gnl[[gn]])/24)
  cor = list(cor(log2rpm.field[,gn], predict))
  names(cor) = gn
  cor.gnl = c(cor.gnl, cor)
}

save(fit.gnl,file="output/MTB/fit.gnl_raw_field")
save(cor.gnl,file="output/MTB/cor.gnl_raw_field")

### Summarize molecular peak time of time-indicating genes
load("output/MTB/fit.gnl_raw_field")
load("output/MTB/cor.gnl_raw_field")

CV <- function(x){
  sd(x)/mean(x)
}
cv.cut <- apply(log2rpm.field[,names(cor.gnl)],2,CV)

### parameter setting
#cv = 0.15
#pdf(sprintf("%s/%s_parameter_raw.pdf",dir.output ,exec.date))
#
#for (th in c(0.8,0.85,0.9,0.93,0.935,0.94,0.945,0.95)){
#  tmp = sum((cor.gnl>th)&(cv.cut > cv)) 
#  hist(unlist(cor.gnl),main=sprintf("r = %s, %s genes", th, tmp))
#  abline(v = th, col = "red")
#}
#
#dev.off()

# th = 0.935 cv = 0.15
MPT.all = fit.gnl # Molecular peak time

th= 0.935
cv= 0.15
gnl.tig = names(cor.gnl)[(unlist(cor.gnl)>th)&(cv.cut > cv)] # Time-indicating genes

MPT = list()
for (gn in gnl.tig){
  tmp =  MPT.all[[gn]] # Molecular peak time
  names(tmp) = gn
  MPT = c(MPT,tmp)
}

gnl.tig = gnl.tig[order(unlist(cor.gnl)[gnl.tig], decreasing = T)] # 143 genes

# TableS2
de <- des2[gnl.tig,2:3]
mptime <- as.numeric(MPT)+19
mptime <- as.numeric(sprintf("%.2f", mptime))

for(i in 1:length(mptime)){
  if(mptime[i]>24){mptime[i]<- mptime[i]-24}
}

r <- as.numeric(cor.gnl[gnl.tig])
r <- as.numeric(sprintf("%.2f", r))
a <- cv.cut[gnl.tig]
a <- as.numeric(sprintf("%.2f", a))

tb <- data.frame(de[,1],r,a,mptime,de[,2])
colnames(tb) <- c("ID","r","a","Molecular Peak Time (hour)","Description")

tb <- tb[order(tb[,4]),]

write.csv(tb, file=sprintf("%s/%s_TableS2.csv", dir.output, exec.date),row.names = F)


### time inference from field data
log2rpm.normalized.tig <- apply(t(log2rpm)[,names(MPT)],2,scale)
rownames(log2rpm.normalized.tig) <- colnames(log2rpm)

BTs = rep(NA, 260)
names(BTs) = colnames(log2rpm)
for (i in 1:260){
  cat(sprintf("%s\n",i))
  ex = log2rpm.normalized.tig[i,]
  tv = unlist(MPT)
  co = rep(NA, 1440)
  names(co) = 0:1439
  for (b in 1:1440){
    predict = sqrt(2)*cos(2*pi*(tv-(b-1)/60)/24)
    co[b] = cor(ex,predict) 
  }
  BTs[i] = 1/60*as.numeric(names(which(co == max(co,na.rm=T))))
}

BTs
BTs_gap = BTs - ZT
BTs_gap[BTs_gap > 12] <- BTs_gap[BTs_gap > 12] - 24


times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

gg <- list()
for (i in 1:5){
  cond <- unique(condition)[i]
  d <- data.frame(v = BTs_gap[condition==cond], h = hour[condition==cond])
  g <- ggplot(d,aes(x=h,y=v))+
    geom_beeswarm(fill=color5[i],shape=i+20, color="black", size=1, stroke = 0.25)+
    geom_hline(yintercept = c(-1,0,1),color="black",lty=2, size=0.25)+
    scale_x_discrete(labels = times)+
    labs(x="", y="r")+
    ggtitle(cond)+
    ylim(c(-2.3,2.3))+
    theme_cowplot(font_size = 7,line_size = 0.25)+
    theme(plot.title = element_text(size=7, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          legend.position = 'none')
  gg <- c(gg,list(g))
}

pp <- plot_grid(gg[[1]],gg[[2]],gg[[3]],gg[[4]],gg[[5]], scale=0.9)

# measurement noise
mn <- rep(NA, 143)
for(i in 1:143){
  a <- abs(log2rpm.field.normalized[,gnl.tig[i]] - 
             sqrt(2)*cos(2*pi*(ZT.field-BTs[1:52])/24))
  mn[i] <- sd(a)
}
mean(mn)
sd(mn)

# Fig. S8a-d
df <- data.frame(value = as.numeric(cor.gnl))
g1 <- ggplot(df, aes(x = value))+
  geom_histogram(fill = "white", colour = "black",
                 size=0.25, breaks=seq(0, 1, by=0.025))+
  geom_vline(xintercept = 0.935, color = color5[1], size = 0.25, linetype="dashed")+
  xlim(0,1)+ 
  ylim(0,500)+
  xlab("r")+
  ylab("The number of gene")+
  theme_cowplot(font_size = 7, line_size=0.25)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

df <- data.frame(value = cv.cut)
g2 <- ggplot(df, aes(x = value))+
  geom_histogram(fill = "white", colour = "black",
                 size=0.25, breaks=seq(0, 4, by=0.1))+
  geom_vline(xintercept = 0.15, color = color5[1], size = 0.25, linetype="dashed")+
  xlim(0,4)+ 
  ylim(0,2500)+
  xlab("a")+
  ylab("The number of gene")+
  theme_cowplot(font_size = 7, line_size=0.25)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

value = as.numeric(MPT)

df <- data.frame(value)
g3 <- ggplot(df, aes(x = value))+
  geom_histogram(fill = "white", colour = "black",
                 size=0.25, breaks=seq(0, 24, by=2))+
  scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2),
                     labels = c("19","21","23","1","3","5","7","9","11","13","15","17","19"))+ 
  scale_y_continuous(limits=c(0,40), breaks=seq(0,40,10))+
  xlab("Molecular peak time")+
  ylab("The number of gene")+
  theme_cowplot(font_size = 7, line_size=0.25)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

# FL/FTH 7:00 replicate 1
ex = log2rpm.normalized.tig[77,]
tv = unlist(MPT)
co = rep(NA, 1440)
names(co) = 0:1439
for (b in 1:1440){
    predict = sqrt(2)*cos(2*pi*(tv-(b-1)/60)/24)
    co[b] = cor(ex,predict) 
}

BT = 1/60*as.numeric(names(which(co == max(co,na.rm=T))))

df <- data.frame(ex,tv)
plot(tv,ex)
curve(sqrt(2)*cos(2*pi*(x-(682-1)/60)/24), from = 0, to = 24, add=T)

g4 <- ggplot(df, aes(x=tv, y=ex))+
     geom_point(size=0.5)+
     geom_line(aes(y=sqrt(2)*cos(2*pi*(tv-(682-1)/60)/24)),
               size=0.5, color=color5[1])+
     geom_vline(xintercept = BT, size=0.5, linetype="dashed")+
     scale_x_continuous(limits=c(0,24), breaks=seq(0,24,2),labels = times)+ 
     scale_y_continuous(limits=c(-3,3), breaks=seq(-3,3,1))+
     xlab("Molecular peak time (hours)")+
     ylab("Normalized expression of\ntime-indicating genes")+
     ggtitle("FL/FTH 7:00 replicate 1")+
     theme_cowplot(font_size = 7, line_size=0.25)+
     theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

p <- plot_grid(plot_grid(g1,g2,g3,g4, ncol=2, align="hv", scale = 0.9), pp,
               nrow=2)

ggsave(p, file=sprintf("%s/%s_FigS8.pdf",dir.output, exec.date),
       width = 150, height = 160, units = "mm")

# Fig. 5ab ----
# load data
sugar <- read.csv("input/190222_sugar.csv",row.names = 1)
lab <- c("Starch(mg/gFW)","Sucrose(mg/gFW)","Fructose(mg/gFW)","Glucose(mg/gFW)",
         "Hexose(mg/gFW)","Soluble sugar(mg/gFW)","NSC(mg/gFW)","Sucrose/Starch ratio(mg/gFW)")

con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(sugar) %>% 
  mutate(con, hour) %>% 
  gather(type,value,-con, -hour) 

d$type <- factor(d$type, levels = colnames(sugar))

d_mean_sd <- d %>%
  group_by(con, hour, type) %>%
  summarize(mean = mean(value, na.rm=T), sd = sd(value, na.rm=T))
errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

g <- d_mean_sd %>%
  group_by(type) %>% 
  nest() %>% 
  mutate(plot = map2(data, type, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, size = 0.25, width = 0.2)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       labs(x="Time (hours)", y="(mg/gFW)")+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.25)+
                       theme(plot.title = element_blank(),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6, colour = "black"),
                             axis.text.y = element_text(size=6, colour = "black"),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- g$plot[[1]] + theme(legend.position="none")
g2 <- g$plot[[2]] + theme(legend.position="none")
gg <- get_legend(g$plot[[1]])


p <- plot_grid(gg, plot_grid(g1,g2), nrow=2, rel_heights = c(1,10))

ggsave(p, file = sprintf("%s/%s_Fig5ab.pdf", dir.output, exec.date),
       width = 160, height = 50, units="mm")

### Table S14 ----
# ANOVA
dat_lm <- d %>%
  group_by(type, hour) %>% 
  nest() %>% 
  mutate(result = purrr::map(data, ~lm(value ~ con, data = .)))

dat_anova <- dat_lm %>% 
  mutate(anova = purrr::map(result, ~Anova(.))) %>% 
  mutate(mcomp = purrr::map(result, ~glht(.,linfct = mcp(con = "Tukey")))) %>% 
  mutate(alpha = purrr::map(mcomp, ~cld(., level = 0.05, decreasing =T) %>% .$mcletters %>% .$Letters))

d2 <- dat_anova %>% 
  dplyr::select(type,hour, alpha) %>% 
  unnest(.,cols=c(alpha)) %>% 
  add_column(., con=rep(c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"),104))

d3 <- inner_join(d_mean_sd, d2)
colnames(d3)[6] <- "significance"
d4 <- d3[order(d3$type,d3$hour),]

#write.csv(d4,"output/210703_multicomparison.csv",row.names = F)

d5 <- d4 %>% 
      dplyr::select(type, hour, con, significance)
d5 <- d5[d5$type=="starch"|d5$type=="sucrose",] %>% 
      spread(., hour, significance)

write.csv(d5,"output/210703_TableS14.csv",row.names = F)

# Fig. 5c-e ----
a <- read.csv("input/genelist_sugar_metabolism.csv")
rownames(a) <- a$Gene.ID..RAP.DB.
ID<- c("Os03g0735000","Os09g0298200","Os12g0641400","Os05g0518600")
lab <- as.character(a[ID,"Gene.name"])
names(lab) <- ID
ylab <- expression(paste(log[2],"(rpm)"))

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.1, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y =  element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

# Experiment 2
con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- factor(time, levels=time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

d <- as_tibble(t(log2rpm_204[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 
d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

p2 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.1, size = 0.25)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color4)+
                       scale_fill_manual(values = color4)+
                       scale_shape_manual(values = c(22:25))+
                       scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       labs(x="Time (hours)", y=ylab)+
                       panel_border(colour="black", size=0.5)+
                       theme(plot.title = element_blank(),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "none",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- g_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p1$plot[[4]]+ theme(legend.position="none")
g5 <- p2$plot[[1]]+ theme(legend.position="none")
g6 <- p2$plot[[2]]+ theme(legend.position="none")
g7 <- p2$plot[[3]]+ theme(legend.position="none")
g8 <- p2$plot[[4]]+ theme(legend.position="none")

p1 <- plot_grid(g1,g2,g3,g4,g5,g6,g7,g8, nrow = 2, align="hv")

ggsave(p1, file=sprintf("%s/%s_Fig5c-e.pdf",dir.output, exec.date),
       width = 160, height = 60, units="mm")

# Fig. 6a ----
condition <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- factor(unique(at$hour), levels=unique(at$hour))
hour <- factor(at$hour, levels = time)

# GO:0005840 (ribosome) mean value
rpm2.ave <- rpm.ave[raw_10,]
withgo <- ulg[ulg[,"GOid"]=="GO:0005840", "locus"]
tmp <- intersect(withgo, rownames(rpm2.ave))
usegene <- intersect(withgo, rownames(rpm2.ave)) 
log2rpm_ribosome <- log2rpm[usegene,hour=="13"]
log2rpm_ribosome_scaled <- apply(log2rpm_ribosome,1,scale)

co <- factor(rep(c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"),each=4),
             levels=c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))

log2rpm_ribosome_scaled_mean <- as_tibble(log2rpm_ribosome_scaled) %>% 
  mutate(co) %>% 
  group_by(co) %>%
  summarize_all(mean, na.rm=T)

value <- t(log2rpm_ribosome_scaled_mean)[2:267,] %>% as.vector() %>% as.numeric()
condition <- c(rep("FIELD", 266),rep("FL/FTH", 266),rep("CL/CTH", 266),rep("FL/CTH", 266),rep("CL/FTH", 266))
d <- data.frame(condition, as.numeric(value))
d$condition <- factor(d$condition, levels = c("FIELD","FL/FTH","CL/CTH", "FL/CTH", "CL/FTH"))

g <- ggplot(d,aes(x=condition,y=value))+
  geom_boxplot(outlier.colour = NA, size=0.25)+
  geom_beeswarm(aes(color=condition), size=0.2)+
  geom_signif(comparisons = list(c("FIELD", "FL/FTH")), annotations="3.21e-75", y_position = 2.2, tip_length = 0.03,
              size=0.25, textsize = 2) +
  geom_signif(comparisons = list(c("FIELD", "CL/CTH")), annotations="1.14e-98", y_position = 2.4, tip_length = 0.03,
              size=0.25, textsize = 2) +
  geom_signif(comparisons = list(c("FIELD", "FL/CTH")), annotations="1.39e-106", y_position = 2.6, tip_length = 0.03,
              size=0.25, textsize = 2) +
  geom_signif(comparisons = list(c("FIELD", "CL/FTH")), annotations="7.60e-51", y_position = 2.8, tip_length = 0.03,
              size=0.25, textsize = 2) +
  scale_color_manual(values = color5)+
  ylim(c(-2,3))+
  labs(x="", y="Normalized expression levels\n(Z-score)")+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  panel_border(colour=1, size=0.25)+
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7, colour = "black"),
        axis.text.y = element_text(size=7, colour = "black"),
        legend.position = 'none')
plot(g)

ggsave(g, file=sprintf("%s/%s_Fig6a.pdf",dir.output, exec.date),
       height = 50, width = 60, units = "mm")


pv <- rep(NA,4)

for(i in 1:4){
cond <- unique(condition)[i+1]
v <- value[c(1:(266), (i*266+1):((i+1)*266))]
b <- as.factor(c(rep("FIELD",266),rep(cond,266)))
d <- data.frame(v,b)
 
a <- wilcox_test(v ~ b, d, distribution="exact")

pv[i] <- pvalue_interval(a)[2]
cat(sprintf("%s\n",i))
}

ad.pv <- p.adjust(pv,method = "BH")

### DEGs between Field and FL/FTH ----
rpm3 <- log2rpm[raw_10,]

d2 <- as_tibble(t(rpm3)) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d_mean_condition_hour <- d2 %>% 
  group_by(con, gene, hour) %>% 
  summarize(mean = mean(expression, na.rm=T)) 

d_mean_condition_hour_FIELD <- d_mean_condition_hour %>% 
  filter(., con == "FIELD")

d_mean_condition_hour_FL_FTH <- d_mean_condition_hour %>% 
  filter(., con == "FL/FTH")

DEG_UP <- list()
DEG_DOWN <- list()

for (i in 1:13){
  time <- as.character(unique(hour))[i]
  a <- read.delim(sprintf("output/TCC/FIELD_FL_FTH/result_0.05_%s.txt",time))
  
  d_f <- d_mean_condition_hour_FIELD %>% 
    filter(., hour == time)
  d_s <- d_mean_condition_hour_FL_FTH %>% 
    filter(., hour == time)
  
  UP <- intersect(as.character(a$rownames.tcc.count.), d_f$gene[d_f$mean > d_s$mean])
  DOWN <- intersect(as.character(a$rownames.tcc.count.), d_f$gene[d_f$mean < d_s$mean])
  
  DEG_UP <- c(DEG_UP,list(UP))
  DEG_DOWN <- c(DEG_DOWN, list(DOWN))
}

DEG_FIELD_FL_FTH_UP <- unique(unlist(DEG_UP))
DEG_FIELD_FL_FTH_DOWN <- unique(unlist(DEG_DOWN))

### FIELD gene, FIELD vs FL/FTH ----
rpm3 <- log2rpm[raw_10,]

d <- as_tibble(t(rpm3)) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d_mean_condition <- d %>% 
  group_by(con, gene) %>% 
  summarize(mean = mean(expression, na.rm=T))

d_mean_FIELD <- d_mean_condition %>% 
  filter(., con=="FIELD") 

d_mean_FL_FTH <- d_mean_condition %>% 
  filter(., con=="FL/FTH") 

d_mean_ratio <- d_mean_FIELD$mean - d_mean_FL_FTH$mean
names(d_mean_ratio) <- d_mean_FIELD$gene

FIELD_up_FL_FTH_2 <- intersect(names(d_mean_ratio[d_mean_ratio > 1]), DEG_FIELD_FL_FTH_UP)
FIELD_down_FL_FTH_2 <- intersect(names(d_mean_ratio[d_mean_ratio < -1]), DEG_FIELD_FL_FTH_DOWN)

# Fig. 6ce ----
a <- read.delim("input/200505_phenylpropanoid_selected.txt")
ID1 <- as.character(a$id)
lab1 <- as.character(a$alias)
names(lab1) <- ID1
ylab <- expression(paste(log[2],"(rpm)"))

a <- read.delim("input/PRgenes.txt")
rownames(a) <- a$rap
b <- a[intersect(a$rap,FIELD_up_FL_FTH_2),]
ID2 <- as.character(b$rap)
lab2 <- as.character(b$alias)
names(lab2) <- ID2

PRgene <- lab2

ID <- c(ID1,ID2)
lab <- c(lab1,lab2)
names(lab) <- ID

con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.1, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=7),
                             legend.position = "top",
                             legend.justification = "center")))


g1 <- p$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p$plot[[1]]+ theme(legend.position="none")
g2 <- p$plot[[2]]+ theme(legend.position="none")
g3 <- p$plot[[3]]+ theme(legend.position="none")
g4 <- p$plot[[4]]+ theme(legend.position="none")
g5 <- p$plot[[5]]+ theme(legend.position="none")
g6 <- p$plot[[6]]+ theme(legend.position="none")
g7 <- p$plot[[7]]+ theme(legend.position="none")
g8 <- p$plot[[8]]+ theme(legend.position="none")
g9 <- p$plot[[9]]+ theme(legend.position="none")
g10 <- p$plot[[10]]+ theme(legend.position="none")
g11 <- p$plot[[11]]+ theme(legend.position="none")
g12 <- p$plot[[12]]+ theme(legend.position="none")
g13 <- p$plot[[13]]+ theme(legend.position="none")
g14 <- p$plot[[14]]+ theme(legend.position="none")

p1 <- plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,NULL,NULL,g11,g12,NULL,NULL, g13,g14,NULL,NULL, ncol=4, align="hv")
pp <- plot_grid(mylegend, p1,
                ncol=1, rel_heights = c(1,12))

ggsave(pp, file=sprintf("%s/%s_Fig6ce.pdf",dir.output, exec.date),
       width = 160, height = 130, units="mm")

# Fig. S6 ----
pca <- prcomp(t(log2rpm[raw_10,]), scale = T)
pc <- pca$x
#PC1 vs PC2
xpc <- pca$sdev[1]^2 / sum(pca$sdev^2) *100
xlab <- sprintf("PC1 (%2.1f%%)",xpc)
ypc <- pca$sdev[2]^2 / sum(pca$sdev^2) *100
ylab <- sprintf("PC2 (%2.1f%%)",ypc)

times <- hour
set <- factor(c(rep("FIELD",52),rep("FL/FTH",52),rep("CL/CTH",52),rep("FL/CTH",52),rep("CL/FTH",52)),
              levels=c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
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

times2 <- rep(time, 5)
set2 <- factor(c(rep("FIELD",13),rep("FL/FTH",13),rep("CL/CTH",13),rep("FL/CTH",13),rep("CL/FTH",13)),
               levels=c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
d2 <- data.frame(times2,set2,pcs)

# plot
g <- ggplot()+
  #geom_point(data=d, aes(x=PC1, y=PC2, group=set, shape=set, fill=set), size=1, stroke=0.25, alpha=0.5)+
  geom_point(data = d2, aes(x=mean_PC1, y=mean_PC2, group=set2, shape=set2, fill=set2), size=2, stroke=0.25, )+
  geom_errorbarh(data = d2, height=2, show.legend = FALSE,size=0.25, alpha=0.25,
                 aes(y=mean_PC2, xmin = (mean_PC1 - sd_PC1), xmax = (mean_PC1 + sd_PC1)))+
  geom_errorbar(data = d2, width=2, show.legend = FALSE, size=0.25, alpha=0.25,
                 aes(x=mean_PC1, ymin = (mean_PC2 - sd_PC2), ymax = (mean_PC2 + sd_PC2)))+
  #geom_path(data=d2, aes(x=mean_PC1, y=mean_PC2, colour=set2,group=set2,),size=0.5, show.legend = FALSE)+
  geom_text_repel(data=d2, aes(x=mean_PC1, y=mean_PC2, label=times2),
                  segment.color = "#84919e", segment.size = 0.25, size=2, show.legend = FALSE)+
  scale_color_manual(values = color5)+
  scale_fill_manual(values = color5)+
  scale_shape_manual(values = c(21:25))+
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

ggsave(g, file=sprintf("%s/%s_FigS6a.pdf",dir.output, exec.date),
       width = 90, height = 90, units = "mm")

#
# plot
g <- ggplot()+
  geom_point(data=d, aes(x=PC1, y=PC2, group=set, shape=set, fill=set), size=1, stroke=0.25, alpha=1)+
  #geom_point(data = d2, aes(x=mean_PC1, y=mean_PC2, group=set2, shape=set2, fill=set2), size=2, stroke=0.25, )+
  #geom_errorbarh(data = d2, height=2, show.legend = FALSE,size=0.25, alpha=0.25,
  #               aes(y=mean_PC2, xmin = (mean_PC1 - sd_PC1), xmax = (mean_PC1 + sd_PC1)))+
  #geom_errorbar(data = d2, width=2, show.legend = FALSE, size=0.25, alpha=0.25,
  #              aes(x=mean_PC1, ymin = (mean_PC2 - sd_PC2), ymax = (mean_PC2 + sd_PC2)))+
  #geom_path(data=d2, aes(x=mean_PC1, y=mean_PC2, colour=set2,group=set2,),size=0.5, show.legend = FALSE)+
  geom_text_repel(data=d, aes(x=PC1, y=PC2, label=times),
                  segment.color = "#84919e", segment.size = 0.25, size=2, show.legend = FALSE)+
  scale_color_manual(values = color5)+
  scale_fill_manual(values = color5)+
  scale_shape_manual(values = c(21:25))+
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

ggsave(g, file=sprintf("%s/%s_FigS6b.pdf",dir.output, exec.date),
       width = 105, height = 105, units = "mm")

# Fig. S7a ----
set <- factor(c(rep("FIELD",13),rep("FL/FTH",13),rep("CL/CTH",13),rep("FL/CTH",13),rep("CL/FTH",13)),
              levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
rpm2.ave <- log2rpm.ave[raw_10,]
soukan <- cor(rpm2.ave, method = "spearman")
souka <- stats::dist(1-soukan)
souk <- hclust(souka, method = "ward.D2")
souka <- NULL
dend <- as.dendrogram(souk) %>% set("branches_lwd",0.25)
labels(dend) <- souk$labels[souk$order]
labels_cex(dend) <- 0.5
cols <- rep(color5,each=13)[souk$order]
labels_colors(dend) <- cols

p1 <- ggplot(as.ggdend(dend))+
  scale_x_continuous(expand = c(0, 0.5), breaks = 1:65, labels = souk$labels[souk$order])+
  scale_y_continuous(expand = c(0, 6))
  #theme(plot.margin= unit(c(1, 1, 1, 1), "lines"))


ggsave(p1, file= sprintf("%s/%s_FigS7a.pdf", dir.output, exec.date),
       width = 160, height = 80, units = "mm")

# Fig. S14 ----
### KEGG
# KEGG photosynthesis
a <- kegg_rice[,"passway"]=="dosa00195"
tmp <- kegg_rice[a,"locus"]
kegg_photosynthesis_gene <- intersect(rownames(rpm2.ave),tmp)

# KEGG antenna
a <- kegg_rice[,"passway"]=="dosa00196"
tmp <- kegg_rice[a,"locus"]
kegg_anntena_gene <- intersect(rownames(rpm2.ave),tmp)

# GO:0015979 (photosynthesis)
withgo <- ulg[ulg[,"GOid"]=="GO:0015979", "locus"]
go_photosynthesis_gene <- intersect(withgo, rownames(rpm2.ave))

# GO:0009765 (photosynthesis, light-harvesting)
withgo <- ulg[ulg[,"GOid"]=="GO:0009765", "locus"]
go_lightharvesting_gene <- intersect(withgo, rownames(rpm2.ave))

# GO:0019684 (photosynthesis, light-reaction)
withgo <- ulg[ulg[,"GOid"]=="GO:0019684", "locus"]
go_lightreaction_gene <- intersect(withgo, rownames(rpm2.ave))

# kegg anntena or go light harvesting (=kegg anntena)
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- factor(unique(at$hour), levels=unique(at$hour))
hour <- factor(at$hour, levels = time)

usegene <- union(kegg_anntena_gene, go_lightharvesting_gene)
ID <- usegene
lab <- des2[ID,"gene.alias"]
names(lab) <- ID
ylab <- expression(paste(log[2],"(rpm)"))

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y =  element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=7),
                             legend.position = "top",
                             legend.justification = "center")))

# Experiment 2
con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- factor(time, levels=time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

d <- as_tibble(t(log2rpm_204[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 
d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

p2 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color4)+
                       scale_fill_manual(values = color4)+
                       scale_shape_manual(values = c(22:25))+
                       scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour="black", size=0.5)+
                       theme(plot.title = element_blank(),
                             axis.title.x =  element_text(size=6),
                             axis.title.y =  element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p1$plot[[4]]+ theme(legend.position="none")
g5 <- p1$plot[[5]]+ theme(legend.position="none")
g6 <- p1$plot[[6]]+ theme(legend.position="none")
g7 <- p1$plot[[7]]+ theme(legend.position="none")
g8 <- p1$plot[[8]]+ theme(legend.position="none")
g9 <- p1$plot[[9]]+ theme(legend.position="none")
g10 <- p1$plot[[10]]+ theme(legend.position="none")
g11 <- p1$plot[[11]]+ theme(legend.position="none")
g12 <- p1$plot[[12]]+ theme(legend.position="none")
g13 <- p1$plot[[13]]+ theme(legend.position="none")
g14 <- p1$plot[[14]]+ theme(legend.position="none")
g15 <- p1$plot[[15]]+ theme(legend.position="none")

h1 <- p2$plot[[1]]+ theme(legend.position="none")
h2 <- p2$plot[[2]]+ theme(legend.position="none")
h3 <- p2$plot[[3]]+ theme(legend.position="none")
h4 <- p2$plot[[4]]+ theme(legend.position="none")
h5 <- p2$plot[[5]]+ theme(legend.position="none")
h6 <- p2$plot[[6]]+ theme(legend.position="none")
h7 <- p2$plot[[7]]+ theme(legend.position="none")
h8 <- p2$plot[[8]]+ theme(legend.position="none")
h9 <- p2$plot[[9]]+ theme(legend.position="none")
h10 <- p2$plot[[10]]+ theme(legend.position="none")
h11 <- p2$plot[[11]]+ theme(legend.position="none")
h12 <- p2$plot[[12]]+ theme(legend.position="none")
h13 <- p2$plot[[13]]+ theme(legend.position="none")
h14 <- p2$plot[[14]]+ theme(legend.position="none")
h15 <- p2$plot[[15]]+ theme(legend.position="none")

p <- plot_grid(mylegend, 
               plot_grid(g1,g2,g3,g4,h1,h2,h3,h4,g5,g6,g7,g8,h5,h6,h7,h8,
                         g9,g10,g11,g12,h9,h10,h11,h12,g13,g14,g15,NULL,h13,h14,h15,NULL,
                         nrow = 8, align="hv"),
               ncol=1, rel_heights = c(1,20))

ggsave(p, file=sprintf("%s/%s_FigS14.pdf",dir.output, exec.date),
       width = 160, height = 200, units="mm")

# Table S13 ----
load("input/Nagano2012/mpar_150124")
mpar_kegg_antenna <- mpar[kegg_anntena_gene,]
write.csv(mpar_kegg_antenna, file=sprintf("%s/%s_TableS13.csv",
                                          dir.output, exec.date))

# Fig. S10 ----
a <- read.csv("input/191025_clockgene.csv")
ID <- as.character(a$ID[c(1:5,7,9:13)])
lab <- as.character(a$name[c(1:5,7,9:13)])
names(lab) <- ID
lab[9:11] <- "ELF4"
ylab <- expression(paste(log[2],"(rpm)"))

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=7),
                             legend.position = "top",
                             legend.justification = "center")))

# Experiment 2
con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- factor(time, levels=time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

d <- as_tibble(t(log2rpm_204[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 
d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

p2 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color4)+
                       scale_fill_manual(values = color4)+
                       scale_shape_manual(values = c(22:25))+
                       scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour="black", size=0.5)+
                       theme(plot.title = element_blank(),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p1$plot[[4]]+ theme(legend.position="none")
g5 <- p1$plot[[5]]+ theme(legend.position="none")
g6 <- p1$plot[[6]]+ theme(legend.position="none")
g7 <- p1$plot[[7]]+ theme(legend.position="none")
g8 <- p1$plot[[8]]+ theme(legend.position="none")
g9 <- p1$plot[[9]]+ theme(legend.position="none")
g10 <- p1$plot[[10]]+ theme(legend.position="none")
g11 <- p1$plot[[11]]+ theme(legend.position="none")

h1 <- p2$plot[[1]]+ theme(legend.position="none")
h2 <- p2$plot[[2]]+ theme(legend.position="none")
h3 <- p2$plot[[3]]+ theme(legend.position="none")
h4 <- p2$plot[[4]]+ theme(legend.position="none")
h5 <- p2$plot[[5]]+ theme(legend.position="none")
h6 <- p2$plot[[6]]+ theme(legend.position="none")
h7 <- p2$plot[[7]]+ theme(legend.position="none")
h8 <- p2$plot[[8]]+ theme(legend.position="none")
h9 <- p2$plot[[9]]+ theme(legend.position="none")
h10 <- p2$plot[[10]]+ theme(legend.position="none")
h11 <- p2$plot[[11]]+ theme(legend.position="none")

p <- plot_grid(mylegend, 
               plot_grid(g1,g2,g3,g4,h1,h2,h3,h4,g5,g6,g7,g8,h5,h6,h7,h8,
                         g9,g10,g11,NULL,h9,h10,h11,NULL,
                         nrow = 6, align="hv"),
               ncol=1, rel_heights = c(1,20))

ggsave(p, file=sprintf("%s/%s_FigS10.pdf",dir.output, exec.date),
       width = 160, height = 160, units="mm")

# Fig. S11 ----
a <- read.csv("input/191025_clockgene.csv")
ID <- as.character(a$ID[c(15:21,23,24,22)])
lab <- as.character(a$name[c(15:21,23,24,22)])
names(lab) <- ID
lab[8] <- "PHOT1"
lab[9] <- "PHOT2"
ylab <- expression(paste(log[2],"(rpm)"))

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=7),
                             legend.position = "top",
                             legend.justification = "center")))

# Experiment 2
con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- factor(time, levels=time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

d <- as_tibble(t(log2rpm_204[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 
d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

p2 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color4)+
                       scale_fill_manual(values = color4)+
                       scale_shape_manual(values = c(22:25))+
                       scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour="black", size=0.5)+
                       theme(plot.title = element_blank(),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p1$plot[[4]]+ theme(legend.position="none")
g5 <- p1$plot[[5]]+ theme(legend.position="none")
g6 <- p1$plot[[6]]+ theme(legend.position="none")
g7 <- p1$plot[[7]]+ theme(legend.position="none")
g8 <- p1$plot[[8]]+ theme(legend.position="none")
g9 <- p1$plot[[9]]+ theme(legend.position="none")
g10 <- p1$plot[[10]]+ theme(legend.position="none")

h1 <- p2$plot[[1]]+ theme(legend.position="none")
h2 <- p2$plot[[2]]+ theme(legend.position="none")
h3 <- p2$plot[[3]]+ theme(legend.position="none")
h4 <- p2$plot[[4]]+ theme(legend.position="none")
h5 <- p2$plot[[5]]+ theme(legend.position="none")
h6 <- p2$plot[[6]]+ theme(legend.position="none")
h7 <- p2$plot[[7]]+ theme(legend.position="none")
h8 <- p2$plot[[8]]+ theme(legend.position="none")
h9 <- p2$plot[[9]]+ theme(legend.position="none")
h10 <- p2$plot[[10]]+ theme(legend.position="none")


p <- plot_grid(mylegend, 
               plot_grid(g1,g2,g3,g4,h1,h2,h3,h4,g5,g6,g7,g8,h5,h6,h7,h8,g9,g10,NULL,NULL,h9,h10,NULL,NULL,
                         nrow = 6, align="hv"),
               ncol=1, rel_heights = c(1,12))

ggsave(p, file=sprintf("%s/%s_FigS11.pdf",dir.output, exec.date),
       width = 160, height = 180, units="mm")

# Fig. S12 ----
a <- read.csv("input/191025_clockgene.csv")
ID <- as.character(a$ID[c(32,33,36:38,26:31)])
lab <- c(rep("RVE",4),"LNK","PIL11","PIL12","PIL13","PIL14","PIL15","PIL16")
names(lab) <- ID
ylab <- expression(paste(log[2],"(rpm)"))

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=7),
                             legend.position = "top",
                             legend.justification = "center")))

# Experiment 2
con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- factor(time, levels=time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

d <- as_tibble(t(log2rpm_204[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 
d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

p2 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color4)+
                       scale_fill_manual(values = color4)+
                       scale_shape_manual(values = c(22:25))+
                       scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour="black", size=0.5)+
                       theme(plot.title = element_blank(),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p1$plot[[4]]+ theme(legend.position="none")
g5 <- p1$plot[[5]]+ theme(legend.position="none")
g6 <- p1$plot[[6]]+ theme(legend.position="none")
g7 <- p1$plot[[7]]+ theme(legend.position="none")
g8 <- p1$plot[[8]]+ theme(legend.position="none")
g9 <- p1$plot[[9]]+ theme(legend.position="none")
g10 <- p1$plot[[10]]+ theme(legend.position="none")
g11 <- p1$plot[[11]]+ theme(legend.position="none")

h1 <- p2$plot[[1]]+ theme(legend.position="none")
h2 <- p2$plot[[2]]+ theme(legend.position="none")
h3 <- p2$plot[[3]]+ theme(legend.position="none")
h4 <- p2$plot[[4]]+ theme(legend.position="none")
h5 <- p2$plot[[5]]+ theme(legend.position="none")
h6 <- p2$plot[[6]]+ theme(legend.position="none")
h7 <- p2$plot[[7]]+ theme(legend.position="none")
h8 <- p2$plot[[8]]+ theme(legend.position="none")
h9 <- p2$plot[[9]]+ theme(legend.position="none")
h10 <- p2$plot[[10]]+ theme(legend.position="none")
h11 <- p2$plot[[11]]+ theme(legend.position="none")

p <- plot_grid(mylegend, 
               plot_grid(g1,g2,g5,h1,h2,h5,g3,g4,NULL,h3,h4,NULL,
                         g6,g7,g8,h6,h7,h8,g9,g10,g11,h9,h10,h11,
                         nrow = 8, align="hv"),
               ncol=1, rel_heights = c(1,20))

ggsave(p, file=sprintf("%s/%s_FigS12.pdf",dir.output, exec.date),
       width = 120, height = 200, units="mm")

# Fig. S15 ----
a <- read.csv("input/genelist_sugar_metabolism.csv")
rownames(a) <- a$Gene.ID..RAP.DB.
ID<- c("Os06g0160700","Os02g0744700","Os04g0624600","Os07g0412100",
       "Os02g0534400","Os04g0409900","Os02g0229400","Os09g0553200",
       "Os01g0851700","Os07g0523600","Os07g0523400")
lab <- as.character(a[ID,"Gene.name"])
names(lab) <- ID
ylab <- expression(paste(log[2],"(rpm)"))

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=7),
                             legend.position = "top",
                             legend.justification = "center")))

# Experiment 2
con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- factor(time, levels=time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

d <- as_tibble(t(log2rpm_204[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 
d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

p2 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color4)+
                       scale_fill_manual(values = color4)+
                       scale_shape_manual(values = c(22:25))+
                       scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour="black", size=0.5)+
                       theme(plot.title = element_blank(),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p1$plot[[4]]+ theme(legend.position="none")
g5 <- p1$plot[[5]]+ theme(legend.position="none")
g6 <- p1$plot[[6]]+ theme(legend.position="none")
g7 <- p1$plot[[7]]+ theme(legend.position="none")
g8 <- p1$plot[[8]]+ theme(legend.position="none")
g9 <- p1$plot[[9]]+ theme(legend.position="none")
g10 <- p1$plot[[10]]+ theme(legend.position="none")
g11 <- p1$plot[[11]]+ theme(legend.position="none")

h1 <- p2$plot[[1]]+ theme(legend.position="none")
h2 <- p2$plot[[2]]+ theme(legend.position="none")
h3 <- p2$plot[[3]]+ theme(legend.position="none")
h4 <- p2$plot[[4]]+ theme(legend.position="none")
h5 <- p2$plot[[5]]+ theme(legend.position="none")
h6 <- p2$plot[[6]]+ theme(legend.position="none")
h7 <- p2$plot[[7]]+ theme(legend.position="none")
h8 <- p2$plot[[8]]+ theme(legend.position="none")
h9 <- p2$plot[[9]]+ theme(legend.position="none")
h10 <- p2$plot[[10]]+ theme(legend.position="none")
h11 <- p2$plot[[11]]+ theme(legend.position="none")

p <- plot_grid(mylegend, 
               plot_grid(g1,g2,g3,g4,h1,h2,h3,h4,g5,g6,g7,g8,h5,h6,h7,h8,
                         g9,g10,g11,NULL,h9,h10,h11,NULL,
                         nrow = 6, align="hv"),
               ncol=1, rel_heights = c(1,20))

ggsave(p, file=sprintf("%s/%s_FigS15.pdf",dir.output, exec.date),
       width = 160, height = 160, units="mm")

# Fig. S16 ----
a <- read.csv("input/genelist_sugar_metabolism.csv")
rownames(a) <- a$Gene.ID..RAP.DB.
ID <- as.character(a$Gene.ID..RAP.DB.[c(145:151,153,155,156,160,161,163)])
lab <- as.character(a[ID,"Gene.name"])
names(lab) <- ID
ylab <- expression(paste(log[2],"(rpm)"))

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=7),
                             legend.position = "top",
                             legend.justification = "center")))

# Experiment 2
con <- condition_204
time <- as.character(time_204)
time[time=="19_2"] <- "19"
time <- factor(time, levels=time)

hour <- as.character(hour_204)
hour[hour=="19_2"] <- "19"
hour <- as.numeric(hour)

d <- as_tibble(t(log2rpm_204[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 
d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

p2 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = ., aes(x=hour, y=mean, group=con, colour=con, shape=con, fill=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke=0.25)+
                       scale_color_manual(values = color4)+
                       scale_fill_manual(values = color4)+
                       scale_shape_manual(values = c(22:25))+
                       scale_x_continuous(breaks = seq(-5,19,2), limits = c(-5, 19))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour="black", size=0.5)+
                       theme(plot.title = element_blank(),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p1$plot[[4]]+ theme(legend.position="none")
g5 <- p1$plot[[5]]+ theme(legend.position="none")
g6 <- p1$plot[[6]]+ theme(legend.position="none")
g7 <- p1$plot[[7]]+ theme(legend.position="none")
g8 <- p1$plot[[8]]+ theme(legend.position="none")
g9 <- p1$plot[[9]]+ theme(legend.position="none")
g10 <- p1$plot[[10]]+ theme(legend.position="none")
g11 <- p1$plot[[11]]+ theme(legend.position="none")
g12 <- p1$plot[[12]]+ theme(legend.position="none")
g13 <- p1$plot[[13]]+ theme(legend.position="none")

h1 <- p2$plot[[1]]+ theme(legend.position="none")
h2 <- p2$plot[[2]]+ theme(legend.position="none")
h3 <- p2$plot[[3]]+ theme(legend.position="none")
h4 <- p2$plot[[4]]+ theme(legend.position="none")
h5 <- p2$plot[[5]]+ theme(legend.position="none")
h6 <- p2$plot[[6]]+ theme(legend.position="none")
h7 <- p2$plot[[7]]+ theme(legend.position="none")
h8 <- p2$plot[[8]]+ theme(legend.position="none")
h9 <- p2$plot[[9]]+ theme(legend.position="none")
h10 <- p2$plot[[10]]+ theme(legend.position="none")
h11 <- p2$plot[[11]]+ theme(legend.position="none")
h12 <- p2$plot[[12]]+ theme(legend.position="none")
h13 <- p2$plot[[13]]+ theme(legend.position="none")

p <- plot_grid(mylegend, 
               plot_grid(g1,g2,g3,g4,h1,h2,h3,h4,g5,g6,g7,g8,
                         h5,h6,h7,h8,g9,g10,g11,g12,h9,h10,h11,h12,
                         g13,NULL,NULL,NULL,h13,
                         ncol = 4, align="hv"),
               ncol=1, rel_heights = c(1,20))

ggsave(p, file=sprintf("%s/%s_FigS17.pdf",dir.output, exec.date),
       width = 160, height = 200, units="mm")

# Fig. S17 ----
# GO:0000785 (chromatin)
rpm2.ave <- rpm.ave[raw_10,]
withgo <- ulg[ulg[,"GOid"]=="GO:0000785", "locus"]
usegene <- intersect(withgo, rownames(rpm2.ave)) 
log2rpm_chromatin <- log2rpm.ave[usegene,]
ylab <- expression(paste(log[2],"(rpm)"))

### linegraph 
tim <- c("19","23","3","7","11","15","19")
# FIELD 
con <- NULL
for(i in 1:35){
  con <- c(con, rep(usegene[i],13))
}

value <- NULL
for(i in 1:35){
  value <- c(value, as.numeric(log2rpm_chromatin[i,1:13]))
}

hours <- rep(time,35)
d <- data.frame(con, value, hours)

g1_chro <- ggplot(d,aes(x=hours, y=value, group=con))+
  geom_line(colour=color5[1],size=0.5)+
  scale_x_discrete(labels = times)+
  ggtitle("FIELD")+
  labs(x="Time (hours)", y=ylab)+
  theme_cowplot(font_size = 6, line_size = 0.5)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")


# FL/FTH 
con <- NULL
for(i in 1:35){
  con <- c(con, rep(usegene[i],13))
}

value <- NULL
for(i in 1:35){
  value <- c(value, as.numeric(log2rpm_chromatin[i,14:26]))
}

hours <- rep(time,35)
d <- data.frame(con, value, hours)

g2_chro <- ggplot(d,aes(x=hours, y=value, group=con))+
  geom_line(colour=color5[2],size=0.5)+
  scale_x_discrete(labels = times)+
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
  con <- c(con, rep(usegene[i],13))
}

value <- NULL
for(i in 1:35){
  value <- c(value, as.numeric(log2rpm_chromatin[i,27:39]))
}

hours <- rep(time,35)
d <- data.frame(con, value, hours)

g3_chro <- ggplot(d,aes(x=hours, y=value, group=con))+
  geom_line(colour=color5[3],size=0.5)+
  scale_x_discrete(labels = times)+
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
  con <- c(con, rep(usegene[i],13))
}

value <- NULL
for(i in 1:35){
  value <- c(value, as.numeric(log2rpm_chromatin[i,40:52]))
}

hours <- rep(time,35)
d <- data.frame(con, value, hours)

g4_chro <- ggplot(d,aes(x=hours, y=value, group=con))+
  geom_line(colour=color5[4],size=0.5)+
  scale_x_discrete(labels = times)+
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
  con <- c(con, rep(usegene[i],13))
}

value <- NULL
for(i in 1:35){
  value <- c(value, as.numeric(log2rpm_chromatin[i,53:65]))
}

hours <- rep(time,35)
d <- data.frame(con, value, hours)

g5_chro <- ggplot(d,aes(x=hours, y=value, group=con))+
  geom_line(colour=color5[5],size=0.5)+
  scale_x_discrete(labels = times)+
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

### output is in 03_experiment2.R ###

# Fig. S18 ----
a <- read.delim("input/PRgenes.txt")
rownames(a) <- a$rap
b <- a[intersect(a$rap,FIELD_down_FL_FTH_2),]

ID1 <- c(as.character(b$rap),"Os01g0731100")
lab1 <- c("PR2","Similar to Pathogen-related protein")
names(lab1) <- ID1

a <- read.delim("input/PR_related.txt")
rownames(a) <- a$rap

ID2 <- as.character(a$rap)[1:3]
lab2 <- as.character(a$alias)[1:3]

ID <- c(ID1, ID2)
lab <- c(lab1, lab2)
names(lab) <- ID
ylab <- expression(paste(log[2],"(rpm)"))

con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.1, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))


g1 <- p$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p$plot[[1]]+ theme(legend.position="none")
g2 <- p$plot[[2]]+ theme(legend.position="none")
g3 <- p$plot[[3]]+ theme(legend.position="none")
g4 <- p$plot[[4]]+ theme(legend.position="none")
g5 <- p$plot[[5]]+ theme(legend.position="none")

# plot
g <- plot_grid(mylegend, plot_grid(g1,g2,NULL,g3,g4,g5,align = "hv", nrow=2),
               nrow=2, rel_heights = c(1,8))

ggsave(g, file=sprintf("%s/%s_FigS18.pdf",dir.output, exec.date),
       width = 150, height = 80, units="mm")

# Fig. S19 ----
ID <- c("Os03g0348200","Os08g0139700","Os08g0167800")
lab <- as.character(des2[ID3,"Description"])
names(lab) <- ID
ylab <- expression(paste(log[2],"(rpm)"))

con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.1, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=6),
                             legend.position = "top",
                             legend.justification = "center")))


g1 <- p$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p$plot[[1]]+ theme(legend.position="none")
g2 <- p$plot[[2]]+ theme(legend.position="none")
g3 <- p$plot[[3]]+ theme(legend.position="none")

p1 <- plot_grid(g1,g2,g3, ncol=3, align="hv")
pp <- plot_grid(mylegend, p1,
                ncol=1, rel_heights = c(1,5))

ggsave(pp, file=sprintf("%s/%s_FigS19.pdf",dir.output, exec.date),
       width = 150, height = 40, units="mm")

# Fig. S20 ----
a <- read.delim("input/herbivore_response(Ye_et_al_2019).txt")
ID <- as.character(a$rap)
lab <- as.character(a$alias)
names(lab) <- ID
ylab <- expression(paste(log[2],"(rpm)"))

# Experiment 1
con <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- unique(at$hour)
hour <- factor(at$hour, levels = time)

d <- as_tibble(t(log2rpm[ID,])) %>% 
  mutate(con, hour) %>% 
  gather(gene,expression,-con, -hour) 

d$gene <- factor(d$gene, levels = ID)

d_mean_sd <- d %>%
  group_by(con, hour, gene) %>%
  summarize(mean = mean(expression, na.rm=T), sd = sd(expression, na.rm=T))

errors <- aes(ymax = mean + sd, ymin = mean - sd)

times <- time
times[times=="19_2"] <- "19"
times <- as.numeric(as.character(times))

p1 <- d_mean_sd %>%
  group_by(gene) %>% 
  nest() %>% 
  mutate(plot = map2(data, gene, 
                     ~ggplot(data = .,
                             aes(x=hour, y=mean, group=con, color=con, fill=con, shape=con))+
                       geom_line(size=0.25)+
                       geom_errorbar(errors, width = 0.2, size = 0.25)+
                       geom_point(color="black",size = 1, stroke = 0.25)+
                       scale_color_manual(values = color5)+
                       scale_fill_manual(values = color5)+
                       scale_shape_manual(values = c(21:25))+
                       scale_x_discrete(labels = times)+
                       ggtitle(sprintf("%s\n%s", .y, lab[.y]))+
                       labs(x="Time (hours)", y=ylab)+
                       theme_cowplot(font_size = 6, line_size = 0.25)+
                       panel_border(colour=1, size=0.5)+
                       theme(plot.title = element_text(size=6, hjust = 0.5),
                             axis.title.x = element_text(size=6),
                             axis.title.y = element_text(size=6),
                             axis.text.x = element_text(size=6),
                             axis.text.y = element_text(size=6),
                             legend.title = element_blank(),
                             legend.text = element_text(size=7),
                             legend.position = "top",
                             legend.justification = "center")))

g1 <- p1$plot[[1]]
mylegend <- get_legend(g1)
g1 <- p1$plot[[1]]+ theme(legend.position="none")
g2 <- p1$plot[[2]]+ theme(legend.position="none")
g3 <- p1$plot[[3]]+ theme(legend.position="none")
g4 <- p1$plot[[4]]+ theme(legend.position="none")
g5 <- p1$plot[[5]]+ theme(legend.position="none")
g6 <- p1$plot[[6]]+ theme(legend.position="none")
g7 <- p1$plot[[7]]+ theme(legend.position="none")
g8 <- p1$plot[[8]]+ theme(legend.position="none")
g9 <- p1$plot[[9]]+ theme(legend.position="none")
g10 <- p1$plot[[10]]+ theme(legend.position="none")
g11 <- p1$plot[[11]]+ theme(legend.position="none")
g12 <- p1$plot[[12]]+ theme(legend.position="none")
g13 <- p1$plot[[13]]+ theme(legend.position="none")
g14 <- p1$plot[[14]]+ theme(legend.position="none")
g15 <- p1$plot[[15]]+ theme(legend.position="none")
g16 <- p1$plot[[16]]+ theme(legend.position="none")

p <- plot_grid(mylegend, 
               plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,
                         g9,g10,g11,g12,g13,g14,g15,g16,
                         nrow = 4, align="hv"),
               ncol=1, rel_heights = c(1,9))

ggsave(p, file=sprintf("%s/%s_FigS20.pdf",dir.output, exec.date),
       width = 160, height = 120, units="mm")
