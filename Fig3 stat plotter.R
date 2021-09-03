### ProportionComparison -- Andrew R Gross -- 2021-08-13
### Compares the proportion of reprogrammed lines with and without karyotypic abnormalities across various age and passage number bins

###############################################################################################
### 1 - Header ################################################################################
library(karyoploteR)
library(ggplot2)
library(reshape2)
library(ggsignif)

###### 1.1 - Functions ########################################################################

zScoreCalculator <- function(y1, y2, n1, n2) {
  p1 = y1/n1
  p2 = y2/n2
  p = (y1 + y2)/(n1 + n2)
  z = (p1 - p2) / sqrt(p * (1-p) * (1/n1 + 1/n2))
  return(z)
}
###############################################################################################
### 2 - Input Data ############################################################################

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates/Fig3 files/')
fig3a <- read.csv('Fig3stats-a.csv', row.names = 1, fileEncoding="UTF-8-BOM")
fig3b <- read.csv('Fig3stats-b.csv', row.names = 1, fileEncoding="UTF-8-BOM")
fig3c <- read.csv('Fig3stats-c.csv', row.names = 1, fileEncoding="UTF-8-BOM")
fig3d <- read.csv('Fig3stats-d.csv', row.names = 1, fileEncoding="UTF-8-BOM")


###############################################################################################
### 3 - Format     ############################################################################
####### 3.1 - Fig3 A: Age
fig3a$expandedPercent <- fig3a$Expanded.abnormal/fig3a$expanded.total*100
fig3a$unexpandedPercent <- fig3a$unexpanded.abnormal/fig3a$unexpanded.total*100
fig3a$zScore <- zScoreCalculator(fig3a$Expanded.abnormal, fig3a$unexpanded.abnormal, fig3a$expanded.total, fig3a$unexpanded.abnormal)
zScores <- c()
for(rowNum in 1:nrow(fig3a)) {
  currentRow <- fig3a[rowNum,]
  newZscore = zScoreCalculator(currentRow[,1], currentRow[,3], currentRow[,2], currentRow[,4])
  zScores <- c(zScores, newZscore)
}
fig3a$zScore <- zScores
fig3a$pValue <- 2*(1-pnorm(fig3a$zScore))
plotA <- melt(t(fig3a[5:6]))

####### 3.1 - Fig3 B: Karyotype number
fig3b$expandedPercent <- fig3b$Expanded.abnormal/fig3b$expanded.total*100
fig3b$unexpandedPercent <- fig3b$unexpanded.abnormal/fig3b$unexpanded.total*100
fig3b$zScore <- zScoreCalculator(fig3b$Expanded.abnormal, fig3b$unexpanded.abnormal, fig3b$expanded.total, fig3b$unexpanded.abnormal)
zScores <- c()
for(rowNum in 1:nrow(fig3b)) {
  currentRow <- fig3b[rowNum,]
  newZscore = zScoreCalculator(currentRow[,1], currentRow[,3], currentRow[,2], currentRow[,4])
  zScores <- c(zScores, newZscore)
}
fig3b$zScore <- zScores
fig3b$pValue <- 2*(1-pnorm(fig3b$zScore))

plotB <- melt(t(fig3b[5:6]))

####### 3.3 - Fig3 C: Passage number by first karyotype
fig3c$expandedPercent <- fig3c$Expanded.abnormal/fig3c$expanded.total*100
fig3c$unexpandedPercent <- fig3c$unexpanded.abnormal/fig3c$unexpanded.total*100
fig3c$zScore <- zScoreCalculator(fig3c$Expanded.abnormal, fig3c$unexpanded.abnormal, fig3c$expanded.total, fig3c$unexpanded.abnormal)
zScores <- c()
for(rowNum in 1:nrow(fig3c)) {
  currentRow <- fig3c[rowNum,]
  newZscore = zScoreCalculator(currentRow[,1], currentRow[,3], currentRow[,2], currentRow[,4])
  zScores <- c(zScores, newZscore)
}
fig3c$zScore <- zScores
fig3c$pValue <- 2*(1-pnorm(fig3c$zScore))

plotC <- melt(t(fig3c[5:6]))

####### 3.4 - Fig3 D: Passage number by first karyotype
fig3d$expandedPercent <- fig3d$Expanded.abnormal/fig3d$expanded.total*100
fig3d$unexpandedPercent <- fig3d$unexpanded.abnormal/fig3d$unexpanded.total*100
fig3d$zScore <- zScoreCalculator(fig3d$Expanded.abnormal, fig3d$unexpanded.abnormal, fig3d$expanded.total, fig3d$unexpanded.abnormal)
zScores <- c()
for(rowNum in 1:nrow(fig3d)) {
  currentRow <- fig3d[rowNum,]
  newZscore = zScoreCalculator(currentRow[,1], currentRow[,3], currentRow[,2], currentRow[,4])
  zScores <- c(zScores, newZscore)
}
fig3d$zScore <- zScores
fig3d$pValue <- 2*(1-pnorm(fig3d$zScore))

plotD <- melt(t(fig3d[5:6]))

###############################################################################################
### 4 - Plot       ############################################################################
####### 4.1 - Fig3 A: 

(f3aPlot <- ggplot(data = plotA, aes(x = Var2, y = value)) + 
  geom_bar(aes(fill = Var1), stat = 'identity', position=position_dodge()) + 
   scale_fill_manual(values = c('black', 'grey50')) +
   scale_y_continuous(breaks = seq(0,50,5), limits = c(0,45.1), expand = c(0,0)) +
  labs(x = 'Donor Age', 
       y = '% Abnormal') +
   theme(axis.title.x = element_text(face="italic", size=14, margin =margin(10,0,0,0)),
         axis.title.y = element_text(face="italic", size=14, margin =margin(0,10,0,0)),
         axis.text.x = element_text(size = 14),
         axis.text.y = element_text(size = 14),
         panel.background = element_rect(fill = 'white', color = 'white', size = 1),
         panel.grid = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 2),
         legend.position = 'none') ) #+ coord_equal(ratio = 0.08)

sigData <- data.frame(x=c(0.875, 1.875, 2.875, 3.875, 4.875), xend=c(1.125, 2.125, 3.125, 4.125, 5.125),
                      y=c(16, 18, 19, 18, 42), annotation=c('***', '**', '****', ' **** ', ' ** '))

f3aPlot <- f3aPlot + geom_signif(stat="identity", 
                      data = sigData,
                      aes(x=x,xend=xend, y=y, yend=y, annotation=annotation),
                      tip_length = 0,
                      vjust = 0) 

####### 4.2 - Fig3 B: 
(f3bPlot <- ggplot(data = plotB, aes(x = Var2, y = value)) + 
    geom_bar(aes(fill = Var1), stat = 'identity', position=position_dodge()) + 
    scale_fill_manual(values = c('black', 'grey50')) +
    scale_y_continuous(breaks = seq(0,50,5), limits = c(0,30.1), expand = c(0,0)) +
    labs(x = 'Karyotype # times', 
         y = '% Abnormal') +
    scale_x_discrete(labels = c("1st","Repeat (2nd - 4th)")) +
    theme(axis.title.x = element_text(face="italic", size=14, margin =margin(10,0,0,0)),
          axis.title.y = element_text(face="italic", size=14, margin =margin(0,10,0,0)),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          panel.background = element_rect(fill = 'white', color = 'white', size = 1),
          panel.grid = element_blank(),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 2),
          legend.position = 'none') ) #+   coord_equal(ratio = 0.08)

sigData <- data.frame(x=c(0.875, 1.875), xend=c(1.125, 2.125),
                      y=c(22, 20), annotation=c('****', ' **** '))

f3bPlot <- f3bPlot + geom_signif(stat="identity", 
                                 data = sigData,
                                 aes(x=x,xend=xend, y=y, yend=y, annotation=annotation),
                                 tip_length = 0,
                                 vjust = 0) 

####### 4.3 - Fig3 C: 
(f3cPlot <- ggplot(data = plotC, aes(x = Var2, y = value)) + 
    geom_bar(aes(fill = Var1), stat = 'identity', position=position_dodge()) + 
    scale_fill_manual(values = c('black', 'grey50')) +
    scale_y_continuous(breaks = seq(0,50,5), limits = c(0,30.1), expand = c(0,0)) +
    labs(x = 'iPSC passage number (range)', 
         y = '% Abnormal') +
    theme(axis.title.x = element_text(face="italic", size=14, margin =margin(10,0,0,0)),
          axis.title.y = element_text(face="italic", size=14, margin =margin(0,10,0,0)),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          panel.background = element_rect(fill = 'white', color = 'white', size = 1),
          panel.grid = element_blank(),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 2),
          legend.position = 'none') ) #+   coord_equal(ratio = 0.08)

sigData <- data.frame(x=c(0.875, 1.875), xend=c(1.125, 2.125),
                      y=c(19, 27), annotation=c('****', ' *** '))

f3cPlot <- f3cPlot + geom_signif(stat="identity", 
                                 data = sigData,
                                 aes(x=x,xend=xend, y=y, yend=y, annotation=annotation),
                                 tip_length = 0,
                                 vjust = 0) 

####### 4.4 - Fig3 D: 
(f3dPlot <- ggplot(data = plotD, aes(x = Var2, y = value)) + 
    geom_bar(aes(fill = Var1), stat = 'identity', position=position_dodge()) + 
    scale_fill_manual(values = c('black', 'grey50')) +
    scale_y_continuous(breaks = seq(0,60,5), limits = c(0,40.1), expand = c(0,0)) +
    labs(x = 'iPSC passage number (range)', 
         y = '% Abnormal') +
    theme(axis.title.x = element_text(face="italic", size=14, margin =margin(10,0,0,0)),
          axis.title.y = element_text(face="italic", size=14, margin =margin(0,10,0,0)),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          panel.background = element_rect(fill = 'white', color = 'white', size = 1),
          panel.grid = element_blank(),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 2),
          legend.position = 'none') ) #+   coord_equal(ratio = 0.03)

sigData <- data.frame(x=c(0.875, 1.875), xend=c(1.125, 2.125),
                      y=c(36, 14), annotation=c('****', ' *** '))

f3dPlot <- f3dPlot + geom_signif(stat="identity", 
                                 data = sigData,
                                 aes(x=x,xend=xend, y=y, yend=y, annotation=annotation),
                                 tip_length = 0,
                                 vjust = 0) 

############################################################################################
### Write to folder
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates/Fig3 files")
### Save plot
tiff(filename= paste0('Fig3A.tiff'), width = 420, height = 350, units = "px", pointsize = 20, res = 100)
f3aPlot
dev.off()

tiff(filename= paste0('Fig3B.tiff'), width = 290, height = 350, units = "px", pointsize = 20, res = 100)
f3bPlot
dev.off()

tiff(filename= paste0('Fig3C.tiff'), width = 320, height = 350, units = "px", pointsize = 20, res = 100)
f3cPlot
dev.off()

tiff(filename= paste0('Fig3D.tiff'), width = 320, height = 350, units = "px", pointsize = 20, res = 100)
f3dPlot
dev.off()


svg(filename= paste0('Fig3A.svg'), width = 420, height = 350, pointsize = 20)
f3aPlot
dev.off()


### Save plot
eps(filename= paste0('Fig3A.eps'), width = 420, height = 350, units = "px", pointsize = 20, res = 100)
f3aPlot
dev.off()

tiff(filename= paste0('Fig3B.tiff'), width = 290, height = 350, units = "px", pointsize = 20, res = 100)
f3bPlot
dev.off()

tiff(filename= paste0('Fig3C.tiff'), width = 320, height = 350, units = "px", pointsize = 20, res = 100)
f3cPlot
dev.off()

tiff(filename= paste0('Fig3D.tiff'), width = 320, height = 350, units = "px", pointsize = 20, res = 100)
f3dPlot
dev.off()


svg(filename= paste0('Fig3A.svg'), width = 420, height = 350, pointsize = 20)
f3aPlot
dev.off()
