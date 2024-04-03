#setwd("##")
library(vegan)
library(ggplot2)
library(ape)

df <- read.table('feature_table_noAIT.txt', row.names = 1,header = T)
df <- df[which (rowSums (df) > 10 ),]
df_group <- read.table("design_noAIT.txt.txt",header = T)

df <- t(df)
bray <- vegdist(df,method = 'bray')
bray <- as.matrix(bray)
write.table(bray,"braycurtis.txt",sep='\t')


options(warn = -1)


site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T) 
}

if (TRUE){
  option_list = list(
    make_option(c("-t", "--type"), type="character", default="braycurtis_noAIT",
                help="Distance type; bray_curtis, bray_curtis_binary, euclidean, jaccard, jaccard_binary, manhatten, unifrac, unifrac_binary [default %default]"),   
    make_option(c("-i", "--input"), type="character", default="",
                help="Input beta distance"),
    make_option(c("-d", "--design"), type="character", default="design_noAIT.txt",
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
 # scale_color_manual(values=c("#1c85f0","#f38903","#80a738")) +
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
  #geom_encircle(aes(fill=group), alpha = 0.1, show.legend = F) +
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
graph2ppt(file="bray_curtis_label.pptx", aspectr=1.7)###导出PPT格式
ggsave('bray_curtis_label.pdf', p, width = 5, height = 3)


###PERMANOVA
set.seed(1)
dune.div <- adonis2(dis ~ group, data = design, permutations = 999, method="bray")

dune.div
p2=p+ geom_text(aes(x=0.5,y=0.3,label='p=0.001'))
    
p2
ggsave('bray-curtis.pdf', p2, width = 4, height = 3)

#####################################################################
###############venn
library(tidyverse)
library(ggvenn)

df1 <- read.table("sum_g.txt",row.names = 1)
df1$LP <- rowSums(df1[,c(1:70)])
df1$ILMAH <- rowSums(df1[,c(71:115)])
df1$ISS <- rowSums(df1[,c(116:139)])



venn <- df1[,140:142]

library(VennDiagram) #加载包

ss <- t(venn)##转置表

a <- names(ss[1,])[ss[1,]>0] #如果第一行的某一列值大于0，则将这一列的列名赋值给a     

b <- names(ss[2,])[ss[2,]>0]

c <- names(ss[3,])[ss[3,]>0]

venn.plot <- venn.diagram(x=list(LP=a,ILMAH=b,ISS=c),
                          filename=NULL, col = "black",
                          lty = "dotted", lwd = 3,
                          fill= c ("red", "blue", "green"),
                          alpha = 0.5,cex=2.0,
                          label.col = c("darkred", "white", "darkblue", "white",
                                        "white", "white", "darkgreen"),
                          fontfamily = "serif",fontface = "bold",
                          cat.col =c("darkred", "darkblue", "darkgreen"),
                          cat.cex = 1.8, cat.fontface = "bold",
                          cat.fontfamily = "serif",
                          cat.dist=0.07,margin=0.2) #A\B为分组的名字。

pdf(file="venn.pdf")
grid.draw(venn.plot)
dev.off()


