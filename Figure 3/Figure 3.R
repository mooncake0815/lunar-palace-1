rm(list=ls())
##setwd("")
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(digest)
library(ggpubr)
library(export)
library(devtools)
library(ggord)
library(wesanderson)
library(RColorBrewer)
library(agricolae)
library(multcomp)
# 1. read input file
alpha = read.table("alpha.txt", header=T, sep="\t", comment.char="")
row.names(alpha) <- alpha$sample
# read design
design = read.table("design.txt", header=T,  sep="\t", comment.char="")
# Extract sample group information, the default is group, can be specified
sampFile = as.data.frame(design$group,row.names = row.names(design))
colnames(sampFile)[1] = "group"
row.names(sampFile) <- row.names(alpha)
## Chao1ntropy index
# add design to alpha

index =cbind(alpha[rownames(sampFile),]$Chao1, sampFile)
colnames(index) = c("Chao1","group") # add Chao1 colname is value
# 2. Statistics and plotting

model = aov(Chao1 ~ group, data=index)
Tukey_HSD <- TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
Tukey_HSD_table <- as.data.frame(Tukey_HSD$group)
out <- LSD.test(model,"group", p.adj="holm")###来住agricolae包
stat = out$groups
stat
# The group information is added to the Index

########The first way
index$stat=stat[as.character(index$group),]$groups
##############The second way
index1=index[,c(1:2)]
for (i in 1:nrow(index1)){{
  if (as.character(index1$group)[i]=="A") index1$stat[i]= as.character(stat["A","groups"])
  else  if (as.character(index1$group)[i]=="B") index1$stat[i]= as.character(stat["B","groups"])
  else if (as.character(index1$group)[i]=="C") index1$stat[i]= as.character(stat["C","groups"]
  )}}
index1$stat=factor(index$stat)
#####The third way

stat1=data.frame(group=rownames(stat),stat=stat[,-1])
index2=merge.data.frame(index1[,-3],stat1,by="group",all = T)
####The fourth way
stat1=data.frame(group=rownames(stat),stat=stat[,-1])
index3 <- inner_join(x = index1[,-3], y = stat1, by = "group")

# Set the group position to 5% of the maximum y value + the maximum Y value of each group
max=max(index[,"Chao1"])
min=min(index[,"Chao1"])
x = index[,c("group","Chao1")]


y = x %>% group_by(group) %>% summarise(Max=max(Chao1))

###Add another line to index ———— 5% of the maximum increase
y
index
index <- inner_join(index,y,by="group")
index$Max <-index$Max+ (max-min)*0.05
colnames(index)[4]="y"
rownames(index)=rownames(sampFile)

p=ggplot(data=index, mapping=aes(x=factor(group,levels =c("G1", "G2")), y=Chao1,fill=group)) +
  geom_boxplot(alpha=1,size=0.7, width=0.6,outlier.size = 0,notch = F)+
  labs(x="Groups", y="Chao1") +
  scale_fill_manual(values=c("#DA4B35","#4CB1CA")) +
  theme_bw()+
  geom_text(data=index, aes(x=group,y=y, label=stat)) +
  geom_jitter(position=position_jitter(0.2), pch=21,size=2,alpha=1)+
  theme(legend.position="right", panel.grid.minor = element_blank())
p


###################################################
##########Replace chao1 with shannon and draw again





##############################################
##########Bray-Curtis


library(vegan)
library(ggplot2)
library(ape)
library(ggalt)

df <- read.csv('feature-table.csv', row.names = 1,header = T)
df <- t(df)
bray <- vegdist(df,method = 'bray')
bray <- as.matrix(bray)
write.table(bray,"bray-curtis.txt",sep='\t')


options(warn = -1)


site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T) 
}

if (TRUE){
  option_list = list(
    make_option(c("-t", "--type"), type="character", default="bray-curtis",
                help="Distance type; bray_curtis, bray_curtis_binary, euclidean, jaccard, jaccard_binary, manhatten, unifrac, unifrac_binary [default %default]"),   
    make_option(c("-i", "--input"), type="character", default="",
                help="Input beta distance"),
    make_option(c("-d", "--design"), type="character", default="design.txt",
                help="design file"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="name of group type"),
    make_option(c("-w", "--width"), type="numeric", default=4,
                help="Width of figure"),
    make_option(c("-e", "--height"), type="numeric", default=2.5,
                help="Height of figure"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  ###default
  if (opts$input==""){opts$input=paste("",opts$type, ".txt", sep = "")}
  if (opts$output==""){opts$output=paste("pcoa_",opts$type, sep = "")}
  
  # Display input and output to confirm that they are correct
  print(paste("The distance matrix file is ", opts$input,  sep = ""))
  print(paste("Type of distance type is ", opts$type,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}




# Read input file

###

dis = read.table(opts$input) 
dis[is.na(dis)]<-0
# read design
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char="") 

#Extract sample group information, the default is group
sampFile = as.data.frame(design[,opts$group],row.names = row.names(design))
colnames(sampFile)[1] = "group"



# Statistics and plotting

# vegan:cmdscale
pcoa = cmdscale(dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
eig = pcoa$eig
points = cbind(points, sampFile[rownames(points),])
colnames(points) = c("x", "y", "z","group") 
library(RColorBrewer)
display.brewer.all()
# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group)) + 
  geom_point(alpha=1, size=3,pointshape = 20,fill.ind = points$group) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=paste("Adonis R2: 0.05; P-value=0.001",sep=""))  + 
  theme_grey() +
  scale_color_manual(values=c("#DA4B35","#4CB1CA")) +
  stat_ellipse(level = 0.6, show.legend = FALSE) +  
  theme(
    text=element_text(size=11,face="plain",color="black"),
    axis.title=element_text(size=12,face="plain",color="black"),
    axis.text = element_text(size=10,face="plain",color="black"),
    legend.title = element_text(size=11,face="plain",color="black"),
    legend.text = element_text(size=11,face="plain",color="black"),
    legend.background = element_blank()
  )+
  geom_vline(xintercept = 0, size = 0.2, color="black",linetype=2) + 
  geom_hline(yintercept = 0, size = 0.2, color="black",linetype=2) +
  geom_encircle(aes(fill=group), alpha = 0.1, show.legend = F) +
  theme_bw()+
  theme(#No legend 
    #axis.text.x=element_text(colour="black",family="Times",size=15), #Sets the font properties of the X-axis scale label
    #axis.text.y=element_text(family="Times",size=15,face="plain"), #Sets the font properties of the Y-axis scale label
    #axis.title.y=element_text(family="Times",size = 15,face="plain"), #Sets the font properties of the Y-axis title
    #axis.title.x=element_text(family="Times",size = 15,face="plain"), #Sets the font properties of the X-axis title
    #plot.title = element_text(family="Times",size=12,face="bold",hjust = 0.5), #Sets the font properties of the general title
    panel.grid.minor = element_blank())

p
library(export)
graph2ppt(file="bray-curtis.pptx", aspectr=1.7)###powerpoint
ggsave('bray-curtis.pdf', p, width = 4, height = 3)




# Sample labeling
p1=p+geom_text_repel(label=paste(rownames(points)),colour="black",size=3)
p1
library(export)
graph2ppt(file="bray_curtis_label.pptx", aspectr=1.7)###
ggsave('bray_curtis_label.pdf', p, width = 5, height = 3)


###PERMANOVA
set.seed(1)
dune.div <- adonis2(dis ~ group, data = design, permutations = 999, method="bray")

dune.div
p2=p+ geom_text(aes(x=0.5,y=0.3,label='p=0.001'))

p2
ggsave('bray-curtis.pdf', p2, width = 4, height = 3)


