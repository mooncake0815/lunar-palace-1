##setwd("")

library(igraph)
library(dplyr)
library(Hmisc)

otu_rare <- read.csv('feature-table.csv',header = T,row.names = 1,stringsAsFactors = F)
otu_rare <- otu_rare[,-71]
otu_rare <- otu_rare [which (rowSums (otu_rare) > 10 ),]

a <- otu_rare[,1:10]
a_CC <- a[,1:3]
a_PC <- a[,4:7]
a_SC <- a[,8:10]
b <- otu_rare[,11:20]
b_CC <- b[,1:3]
b_PC <- b[,4:7]
b_SC <- b[,8:10]
c <- otu_rare[,21:30]
c_CC <- c[,1:3]
c_PC <- c[,4:7]
c_SC <- c[,8:10]
d <- otu_rare[,31:40]
d_CC <- d[,1:3]
d_PC <- d[,4:7]
d_SC <- d[,8:10]
e <- otu_rare[,41:50]
e_CC <- e[,1:3]
e_PC <- e[,4:7]
e_SC <- e[,8:10]
f <- otu_rare[,51:60]
f_CC <- f[,1:3]
f_PC <- f[,4:7]
f_SC <- f[,8:10]
g <- otu_rare[,61:70]
g_CC <- g[,1:3]
g_PC <- g[,4:7]
g_SC <- g[,8:10]

split_otu <- list(a_SC,b_SC,c_SC,d_SC,e_SC,f_SC,g_SC)
## Define some colors
col_g <- "#C1C1C1"
cols <- c("#DEB99B" ,"#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D" ,"#6E4821", "#A4B423",
          "#C094DF" ,"#DC95D8" ,"#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F", "#30455C", "#F96C72","#5ED2BF")


g <- lapply(split_otu,function(x){
  occor<-WGCNA::corAndPvalue(t(x)/colSums(x),method = 'spearman')
  # occor<-WGCNA::corAndPvalue(t(x),method = 'pearson')
  mtadj<-multtest::mt.rawp2adjp(unlist(occor$p),proc='BH')
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(t(x)/colSums(x))[2])
  ## R value
  occor.r<-occor$cor
  diag(occor.r) <- 0
  occor.r[occor.p>0.05|abs(occor.r)<0.7] = 0
  occor.r[is.na(occor.r)]=0
  g <-  graph.adjacency(occor.r, weighted = TRUE, mode = 'undirected')
  # Delete autocorrelation
  g <- simplify(g)
  # Deleting isolated Nodes
  g <- delete.vertices(g, which(degree(g)==0) )
  return(g)
})
save(g,file = 'network-SC.rda')
load('network-SC.rda')


# Single graph example

g7 <- g[[7]]
# plot(g[[1]])

## Set the weight of the network in preparation for calculating modularity
E(g7)$correlation <- E(g7)$weight
E(g7)$weight <- abs(E(g7)$weight)

## Computing network module
set.seed(007)
V(g7)$modularity <- membership(cluster_fast_greedy(g7))


## Add nodes and edge colors

#Set node color by module

#The first 18 modules containing the number of nodes are selected to give different colors, and the remaining modules are given gray

V(g7)$label <- V(g7)$name
V(g7)$label <- NA
modu_sort <- V(g7)$modularity %>% table() %>% sort(decreasing = T)
top_num <- 18
modu_name <- names(modu_sort[1:18])
modu_cols <- cols[1:length(modu_name)]
names(modu_cols) <- modu_name
V(g7)$color <- V(g7)$modularity
V(g7)$color[!(V(g7)$color %in% modu_name)] <- col_g
V(g7)$color[(V(g7)$color %in% modu_name)] <- modu_cols[match(V(g7)$color[(V(g7)$color %in% modu_name)],modu_name)]
V(g7)$frame.color <- V(g7)$color


## Set the color of the edge
#The color of the edge is consistent with the color of the module

#Since the edge connects two nodes, if both nodes belong to the same module, we give them the color of the module, 
#and if two nodes belong to different modules, we give them gray


E(g7)$color <- col_g
for ( i in modu_name){
  col_edge <- cols[which(modu_name==i)]
  otu_same_modu <-V(g7)$name[which(V(g7)$modularity==i)]
  E(g7)$color[(data.frame(as_edgelist(g7))$X1 %in% otu_same_modu)&(data.frame(as_edgelist(g7))$X2 %in% otu_same_modu)] <- col_edge
}



## Calculate the layout of the network and output it

sub_net_layout <- layout_with_fr(g7, niter=999,grid = 'nogrid')
## Visualization and output
par(font.main=4)
plot(g7,layout=sub_net_layout, edge.color = E(g7)$color,vertex.size=2)
title(main = paste0('Nodes=',length(V(g7)$name),', ','Edges=',nrow(data.frame(as_edgelist(g7)))))


pdf(paste0("Example 8.pdf"), encoding="MacRoman", width=6, height=6)
par(font.main=4)
plot(g7,layout=sub_net_layout, edge.color = E(g7)$color,vertex.size=2)
title(main = paste0('Nodes=',length(V(g7)$name),', ','Edges=',nrow(data.frame(as_edgelist(g7)))))
dev.off()


# Multi-picture batch output

pdf(paste0("Example SC.pdf"), encoding="MacRoman", width=15, height=9)
par(mfrow=c(1,7),mar=c(0,0,1,0),font.main=4)
for(i in 1:7){
  g1 <- g[[i]]
  E(g1)$correlation <- E(g1)$weight
  E(g1)$weight <- abs(E(g1)$weight)
  set.seed(007)
  V(g1)$modularity <- membership(cluster_fast_greedy(g1))
  
  V(g1)$label <- V(g1)$name
  V(g1)$label <- NA
  modu_sort <- V(g1)$modularity %>% table() %>% sort(decreasing = T)
  
  top_num <- 18
  modu_name <- names(modu_sort[1:18])
  modu_cols <- cols[1:length(modu_name)]
  names(modu_cols) <- modu_name
  V(g1)$color <- V(g1)$modularity
  V(g1)$color[!(V(g1)$color %in% modu_name)] <- col_g
  V(g1)$color[(V(g1)$color %in% modu_name)] <- modu_cols[match(V(g1)$color[(V(g1)$color %in% modu_name)],modu_name)]
  V(g1)$frame.color <- V(g1)$color
  
  E(g1)$color <- col_g
  for ( i in modu_name){
    col_edge <- cols[which(modu_name==i)]
    otu_same_modu <-V(g1)$name[which(V(g1)$modularity==i)]
    E(g1)$color[(data.frame(as_edgelist(g1))$X1 %in% otu_same_modu)&(data.frame(as_edgelist(g1))$X2 %in% otu_same_modu)] <- col_edge
  }
  
  
  sub_net_layout <- layout_with_fr(g1, niter=999,grid = 'nogrid')
  plot(g1,layout=sub_net_layout, edge.color = E(g1)$color,vertex.size=2)
  title(main = paste0('Nodes=',length(V(g1)$name),', ','Edges=',nrow(data.frame(as_edgelist(g1)))))
}
dev.off()


##############
library(ggplot2)
df <- read.table("value.txt",header = T)
p1 <- ggplot(data=df, mapping=aes(x=time,y=node))+
  geom_bar(stat="identity", fill="#A4CFD9",width=0.8)
p1
p2 <- ggplot(data=df, mapping=aes(x=time,y=edge))+
  geom_bar(stat="identity", fill="#D9745C",width=0.8)
p2

library(patchwork)
p = p1+p2
p










#####################
# ADONIS

data=read.csv("feature-table.csv", header=TRUE, row.names=1)
Group=read.table("design.txt", header=TRUE)
rownames(Group)=Group[,1]
Group=Group[,-1]

means=apply(data, 1, mean)
otu=data[names(means[means>5]),]
otu=t(otu)
#Adonis analysis and Varpart analysis 
library(vegan)
adonis=adonis2(otu~group*time*site2,Group,permutations=999,method = "bray")
adonis

RDA=varpart(otu,Group[1],Group[3],Group[4],transfo="hel",permutations=999)
RDA


#######混合效应模型
library(lme4)
alpha <- read.table("alpha.txt",header = T)
design <- read.table("design.txt",header = T)
alpha$group <- design$group
alpha$time <- design$day
alpha$site <- design$site2
avd <- read.csv("avd.csv")

model1 <- lmerTest::lmer(ACE ~ time + (1|group), data = alpha)
anova(model1)
model2 <- lmerTest::lmer(Shannon ~ time + (1|group), data = alpha)
anova(model2)
model3 <- lmerTest::lmer(avd ~ time + (1|group), data = avd)
anova(model3)

library(ggplot2)
library(Rmisc)
library(cowplot)
tgc1 <- summarySE(alpha, measurevar="Chao1", groupvars="time")
tgc1$time <- c("D58","D90","D123","D156","D216","D316","D330")
P1 <- ggplot(alpha,aes(x =factor(time,levels =c("D58", "D90","D123","D156","D216","D316","D330")),
                       y = ACE)) +   
  geom_point( alpha = 1, shape = 17, size = 3,color="#3CB371") +
  stat_summary(fun=median, geom="line", aes(group=1))+
  theme_bw()

P1

ggsave("ACE-time.pdf",P1,width=6,height = 5)
P2 <- ggplot(alpha,aes(x =factor(time,levels =c("D58", "D90","D123","D156","D216","D316","D330")),
                       y = Shannon)) +   
  geom_point(alpha = 1, shape = 17, size = 3,color="skyblue") +
  stat_summary(fun=median, geom="line", aes(group=1))+
  theme_bw()

P2
ggsave("Shannon-time.pdf",P2,width=6,height = 5)
P3 <- ggplot(avd,aes(x =factor(time,levels =c("D58", "D90","D123","D156","D216","D316","D330")),
                     y = avd)) +   
  geom_point(alpha = 1, shape = 17, size = 3,color="MediumSlateBlue") +
  stat_summary(fun=median, geom="line", aes(group=1))+
  theme_bw()

P3
ggsave("avd-time.pdf",P3,width=6,height = 5)
