# blast ----
rowname <- c("query_name","query_length","subject_name",
             "subject_length","alignment_length","query_frame",
             "query_start","query_end","subject_frame",
             "subject_start","subject_end","query_sequence",
             "subject_sequence","identity%","evalue","bitscore",
             "mismatch","gapopen","gaps","taxonomyID",
             "scientific_name","name","blast_name","kingdom",
             "title")

a <- read.delim("input/OsaYH_2-7_unmapped/blast_result.tsv",header = F)
colnames(a) <- rowname

b <- a[a$blast_name!="monocots" & a$blast_name!="eudicots",]
rownames(b) <- b$query_name

# load sample attribute --------------
at <- read.csv("input/OsaYH_2-7_unmapped/SmartGC_SampleAttribute.csv")
genelist <- read.table("input/OsaYH_2-7_unmapped/genes.list")[,1] %>% as.character()
b <- b[genelist[53680:53698],]

# load data ----
load("input/210703_rawcnt_unmapped")

# Table17 ----
condition <- factor(at$condition,
                    levels=c("exp1_FIELD","exp1_FL-FTH", "exp1_CL-CTH", "exp1_FL-CTH", "exp1_CL-FTH",
                    "exp2_FL-FTH", "exp2_CL-CTH", "exp2_FL-CTH", "exp2_CL-FTH"))

rawcnt_sum <- data.frame(t(rawcnt_extra)) %>% mutate(condition) %>% 
  group_by(condition) %>% 
  summarize_all(sum) %>% 
  dplyr::select(., -condition) %>% 
  t()
colnames(rawcnt_sum) <- unique(condition)
rawcnt_blast <- cbind(rawcnt_sum,b[,c(20:25,1:19)])

write.csv(rawcnt_blast, file=sprintf("%s/%s_TableS17.csv",dir.output, exec.date))

### Fig. S21 ----
condition <- factor(at$condition, levels = c("FIELD","FL/FTH","CL/CTH","FL/CTH","CL/FTH"))
time <- factor(unique(at$hour), levels=unique(at$hour))
hour <- factor(at$hour, levels = time)

rawcnt_virus <- rawcnt[des$NormalizationGroup=="virus",]
rawcnt_virus_max <- apply(rawcnt_virus, 1, max) %>% 
  sort(., decreasing = T)
rawcnt_virus_mean <- apply(rawcnt_virus, 1, mean) %>% 
  sort(., decreasing = T)

length(which(rawcnt_virus_max >= 1)) # 115

cvrd_virus <- cvrd3[des$NormalizationGroup=="virus",]
cvrd_virus_max <- apply(cvrd_virus, 1, max) %>% 
  sort(., decreasing = T)
cvrd_virus_mean <- apply(cvrd_virus, 1, mean) %>% 
  sort(., decreasing = T)

plot(cvrd_virus_max[rawcnt_virus_max > 1], log2(rawcnt_virus_max[rawcnt_virus_max > 1]+1))

d <- data.frame(count=log2(rawcnt_virus_max[rawcnt_virus_max > 1]+1),
                cvrd=cvrd_virus_max[rawcnt_virus_max > 1])

ylab <- expression(paste(log[2],"(Max read number + 1)"))

p <- ggplot(d,aes(x=cvrd, y=count))+
  geom_point(size=2, shape=21, stroke=0.5)+
  labs(x="Max Coverage (depth > 3) [%]", y=ylab)+
  theme_cowplot(font_size = 7, line_size = 0.25)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7, colour = "black"),
        axis.text.y = element_text(size=7, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_blank())

ggsave(p, file=sprintf("%s/%s_FigS21.pdf",dir.output, exec.date),
       width = 80, height = 80, units="mm")

# Fig. 6f ----
a <- cvrd1[des$NormalizationGroup=="virus",]
des_virus <- des[des$NormalizationGroup=="virus",]
b <- apply(a,1,max)
c <- sort(b,decreasing = T)
d <- c[c>0]
e <- des[names(d),"Short_description"]

rpm_virus <- rpm[des$NormalizationGroup=="virus",]
rpm_virus2124 <- rpm_virus[grep("Oryza sativa alphaendornavirus",des_virus$Short_description),]
rpm_virus_204 <- rpm_204[des$NormalizationGroup=="virus",]
rpm_virus2124_204 <-  rpm_virus_204[grep("Oryza sativa alphaendornavirus",des_virus$Short_description),]

log2rpm_virus2124 <- log2(c(rpm_virus2124,rpm_virus2124_204)+1)

a <- read.delim("input/PRgenes.txt")
rownames(a) <- a$rap
b <- a[intersect(a$rap,FIELD_up_FL_FTH_2),]
ID2 <- as.character(b$rap)
lab2 <- as.character(b$alias)
names(lab2) <- ID2
log2rpm_prgene <- cbind(log2rpm[ID2,],log2rpm_204[ID2,])  

col <- factor(c(rep("FIELD",52),rep("Exp1_the_others",208),rep("Exp2",204)),
              levels = c("FIELD","Exp1_the_others","Exp2"))
cols <- c(color5[1],"black","gray")
d <- data.frame(virus=log2rpm_virus2124,t(log2rpm_prgene),col)
lab <- PRgene

g <- list()
for (i in 1:6){
  gene <- rownames(log2rpm_prgene)[i]
  p <- ggplot(d,aes_(x=as.name(gene), y=as.name("virus"), color=col, group=col))+
    geom_point(size=1, shape=21, stroke=0.2)+
    scale_color_manual(values = cols)+
    ggtitle(sprintf("%s\n%s", gene, lab[gene]))+
    theme_cowplot(font_size = 7, line_size = 0.25)+
    panel_border(colour=1, size=0.5)+
    theme(plot.title = element_text(size=6, hjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=6, colour = "black"),
          axis.text.y = element_text(size=6, colour = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.position = "top",
          legend.justification = "center")
  g <- c(g,list(p))
}

mylegend <- get_legend(g[[1]])
g1 <- g[[1]] + theme(legend.position="none")
g2 <- g[[2]] + theme(legend.position="none")
g3 <- g[[3]] + theme(legend.position="none")
g4 <- g[[4]] + theme(legend.position="none")
g5 <- g[[5]] + theme(legend.position="none")
g6 <- g[[6]] + theme(legend.position="none")

g <- plot_grid(mylegend, plot_grid(g1,g2,g3,g4,g5,g6, ncol=2, align = "hv"),
               ncol=1, rel_heights = c(1,9))

ggsave(g, file=sprintf("%s/%s_Fig6f.pdf",dir.output, exec.date),
       width = 40, height = 80, units="mm")

# Fig S23 ----
rawcnt_extra_selected <- rawcnt_extra[c(2,11,12,15,16,18),]

col <- factor(c(rep("FIELD",52),rep("Exp1_the_others",208),rep("Exp2",204)),
              levels = c("FIELD","Exp1_the_others","Exp2"))
cols <- c(color5[1],"black","gray")
xlab <- expression(paste(log[2],"(rpm)"))

d <- data.frame(t(rawcnt_extra_selected),t(log2rpm_prgene),col)

gg <- list()
for (j in 1:6){
  prg <- colnames(d)[j+6]
  
  g <- list()
  for (i in 1:6){
    gene <- colnames(d)[i]
    p <- ggplot(d,aes_(x=as.name(prg), y=as.name(gene), color=col, group=col))+
      geom_point(size=1, shape=21, stroke=0.2)+
      scale_color_manual(values = cols)+
      ggtitle(sprintf("contig%s",i))+
      labs(x=xlab, y="the number\nof reads")+
      theme_cowplot(font_size = 7, line_size = 0.25)+
      panel_border(colour=1, size=0.25)+
      theme(plot.title = element_text(size=6, hjust=0.5),
            axis.title.x = element_text(size=6, colour = "black"),
            axis.title.y = element_text(size=6, colour = "black"),
            axis.text.x = element_text(size=6, colour = "black"),
            axis.text.y = element_text(size=6, colour = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size=6),
            legend.position = "top",
            legend.justification = "center")
    g <- c(g,list(p))
  }
  
  mylegend <- get_legend(g[[1]])
  title <- ggdraw() + draw_label(PRgene[j], fontface='bold', size = 7)
  g1 <- g[[1]] + theme(legend.position="none")
  g2 <- g[[2]] + theme(legend.position="none")
  g3 <- g[[3]] + theme(legend.position="none")
  g4 <- g[[4]] + theme(legend.position="none")
  g5 <- g[[5]] + theme(legend.position="none")
  g6 <- g[[6]] + theme(legend.position="none")
  
  pp <- plot_grid(title, plot_grid(g1,g2,g3,g4,g5,g6, ncol=3, align="hv"), 
                  ncol=1, rel_heights = c(0.1, 0.9))
  gg <- c(gg, list(pp))
  
}

ppp <- plot_grid(mylegend, plot_grid(gg[[1]],gg[[2]],gg[[3]],gg[[4]],gg[[5]],gg[[6]], ncol=2),
                 ncol=1, rel_heights = c(1,19))

ggsave(ppp, file=sprintf("%s/%s_FigS23.pdf",dir.output, exec.date),
       width = 160, height = 170, units="mm")
