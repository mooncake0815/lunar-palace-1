# setwd("##")

library(ggplot2)

pcr <- read.table("qpcr.txt")
design <- read.table("design.txt",header = T,row.names = 1)
pcr1=pcr
pcr2=pcr
pcr1$group = design$group
pcr2$site = design$site2
pcr1 <- melt(pcr1, id.vars = "group", variable.name = "Variable", value.name = "Value")
pcr2 <- melt(pcr2, id.vars = "site", variable.name = "Variable", value.name = "Value")
pcr1 <- pcr1[pcr1$Value != 0, ]
pcr2 <- pcr2[pcr2$Value != 0, ]
p1=ggplot(pcr1, aes(x = group, y = Value, fill = group, color = group)) +
  geom_point(size=2,pch=17)+
  geom_boxplot(color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("#DA4B35","#4CB1CA"),name = ' ')+
  scale_color_manual(values = c("#DA4B35","#4CB1CA"),name = 'group')+
  facet_wrap(~ Variable, nrow = 1) +
  xlab("Variable") +
  ylab("Value") +
  labs(fill = "group")
p1
p2=ggplot(pcr2, aes(x = site, y = Value, fill = site, color = site)) +
  geom_point(size=2,pch=17)+
  geom_boxplot(color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("#A4CFD9","#D9745C","#DEC033"),name = ' ')+
  scale_color_manual(values = c("#A4CFD9", "#D9745C", "#DEC033"),name='site') +
  facet_wrap(~ Variable, nrow = 1) +
  xlab("Variable") +
  ylab("Value") +
  labs(fill = "site")
p2



library(tidyverse) 
library(reshape2) 
df <- read.table("qpcr-time.txt",header = T)
head(df)

df1 <- df %>%
  melt(id=c('treat','Day'))
df1 <- df1[df1$value != 0, ]

topbar <- function(x){
  return(mean(x)+sd(x)/sqrt(length(x))) 
}
bottombar <- function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
data1 <- df1

p3 <- ggplot(data1,aes(Day,value,color=treat))+
  geom_rect(aes(xmin=2,xmax=4,ymin=(-Inf),ymax=Inf),
            fill='grey90',color='grey90')+
#   geom_rect(aes(xmin=216,xmax=330,ymin=(-Inf),ymax=Inf),
#            fill='#DA4B35',color='#DA4B35',alpha=0.1)+ geom_rect(aes(xmin=216,xmax=330,ymin=(-Inf),ymax=Inf),
#           fill='#DA4B35',color='#DA4B35',alpha=0.1)+
  geom_vline(xintercept =c(1,2,4,7),linetype=2,cex=0.5)+
  stat_summary(geom = 'line',fun='mean',cex=1)+
  stat_summary(geom = 'errorbar',
               fun.min = bottombar,fun.max = topbar,
               width=0.1,cex=0.8,aes(color=treat))+
  stat_summary(geom = 'point',fun='mean',aes(fill=treat),
               size=2,pch=21,color='black')+
  theme_classic(base_size = 15)+
  theme(legend.position = 'none')+
  ylim(2,6)+
  scale_color_manual(values = c(brewer.pal(6,"Set2")[2:6]))+
  scale_fill_manual(values = c(brewer.pal(6,"Set2")[2:6]))+
  labs(y='Gene copies')
p3
ggpubr::ggarrange(p1,p2,p3, nrow = 3, ncol = 1,widths = 6,
                  heights = 10)

ggsave("Fig8.pdf",width =6,height =8,units ="in")
