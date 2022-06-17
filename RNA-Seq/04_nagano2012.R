### Nagano et al., Cell, 2012 ----
### load data
load("input/Nagano2012/m.ex.2009")
load("input/Nagano2012/m.ex.itoh")
load("input/Nagano2012/m.ex.ts")
load("input/Nagano2012/wt_2009_tsukuba")
at_2009 <- dget("input/Nagano2012/attribute_2009_110301")

# Field
a <- at_2009[1:108,]
att <- a[order(a$time),]
att[att$hour==24,"day"] <- att[att$hour==24,"day"]+1
att[att$hour==24,"hour"] <- 0

ra <- rep(NA, nrow(att))
for(i in 1:nrow(att)){
  w <- wt$month==att[i,]$month&
    wt$day==att[i,]$day&
    wt$hour==att[i,]$hour&
    wt$min==att[i,]$min
  ra[i] <- wt[w,"radiation"]
}

L_D <- ra > 0.2
at_2009_Field <- cbind(att,L_D)

# GC
a <- at_2009[109:152,]
at_2009_GC <- a[order(a$time),]
L_D <- c("L","D","D","D","D","D","L","L",
         "L","D","D","D","D","D","D","D","L","L",
         "L","D","D","D","D","D","D","D","L","L",
         rep(c("D","D","L","L"),4))
L_D <- L_D=="L"
at_2009_GC <- cbind(at_2009_GC,L_D)

colnames(m.ex.2009) <- 1:108
colnames(m.ex.itoh) <- 109:136
colnames(m.ex.ts) <- 137:152

### mean ----
ex_n8_2009 <- m.ex.2009[,1:96]
ex <- cbind(ex_n8_2009, m.ex.ts)
cond <- c(rep("Field_N8_2009", ncol(ex_n8_2009)), rep("GC_Temp", ncol(m.ex.ts)))

cond <- factor(cond, levels = c("Field_N8_2009","GC_Temp"))

genelist <- sort(intersect(rownames(ex_n8_2009),rownames(rpm2.ave)))

dat_mean <- as_tibble(t(ex[genelist,])) %>% 
  mutate(cond) %>% 
  gather(genes,value,-cond) %>% 
  group_by(cond, genes) %>% 
  summarize(mean = mean(value, na.rm=T))

dat_GC_Temp <- dat_mean[dat_mean$cond=="GC_Temp",]$mean
dat_Field_N8_2009 <- dat_mean[dat_mean$cond=="Field_N8_2009",]$mean

d_relative <- dat_Field_N8_2009 -dat_GC_Temp
names(d_relative) <- genelist

### comparison ----
gene_cell_2 <- genelist[d_relative > 1]
gene_cell_0.5 <- genelist[d_relative < -1]

Field_up_2 <- intersect(gene_cell_2, FIELD_up_FL_FTH_2)
Field_down_2 <- intersect(gene_cell_0.5, FIELD_down_FL_FTH_2)

### Fig. 6d ----
cell_2 <- as.numeric((d_relative > 1) & (genelist %in% FIELD_up_FL_FTH_2))
cell_0.5 <- as.numeric((d_relative < -1) & (genelist %in% FIELD_down_FL_FTH_2))
col <- as.factor(cell_2 + cell_0.5*3 + 1)
d <- data.frame(d_FL_FTH = d_mean_ratio[genelist], d_relative, col)
cols <- c("black",color5[1],color5[2])

g <- ggplot(d, aes(x = d_FL_FTH, y = d_relative, color = col, group=col, alpha=col))+
  geom_point(size=0.25)+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  scale_alpha_manual(values = c(0.2,1,1))+
  geom_hline(yintercept=0, col = "black", linetype="dashed", size=0.25)+
  geom_vline(xintercept=0, col = "black",  linetype="dashed", size=0.25)+
  xlim(-8,8)+
  ylim(-8,8)+
  labs(x="FIELD vs FL/FTH", y="Paddy field vs Growth chamber\n(Nagano2012)")+
  ggtitle("")+
  theme_cowplot(line_size=0.25)+
  panel_border(colour=1, size=0.5)+
  theme(plot.title = element_text(size=7, hjust = 0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7, colour = "black"),
        axis.text.y = element_text(size=7, colour = "black"),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.key.size = unit(3,"line"))
plot(g)

ggsave(g, file=sprintf("%s/%s_Fig6d.pdf",dir.output, exec.date),
       width = 55, height =55, units = "mm")

### Table S16,17,18,19 ----
# S16S17
up_FIELD_FL_FTH <- d_mean_ratio[FIELD_up_FL_FTH_2]
write.csv(up_FIELD_FL_FTH,file=sprintf("%s/%s_TableS16.csv",dir.output,exec.date))

down_FIELD_FL_FTH <- d_mean_ratio[FIELD_down_FL_FTH_2]
write.csv(down_FIELD_FL_FTH,file=sprintf("%s/%s_TableS17.csv",dir.output,exec.date))

# S18S19
up_FIELD_FL_FTH <- d_mean_ratio[Field_up_2]
up_FIELD_FL_FTH_cell <- d_relative[Field_up_2]
d <- data.frame(up_FIELD_FL_FTH,up_FIELD_FL_FTH_cell)
write.csv(d,file=sprintf("%s/%s_TableS18.csv",dir.output,exec.date))

down_FIELD_FL_FTH <- d_mean_ratio[Field_down_2]
down_FIELD_FL_FTH_cell <- d_relative[Field_down_2]
d <- data.frame(down_FIELD_FL_FTH,down_FIELD_FL_FTH_cell)
write.csv(d,file=sprintf("%s/%s_TableS19.csv",dir.output,exec.date))
