library("ggplot2")
library("grid")
library('gridExtra')
library(magrittr)


g <- ggplot(mtcars, aes(x = hp, y = mpg)) +
  geom_point() +
  theme(axis.title.x = element_text(angle = -45,
                                    hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(angle = -45,
                                    hjust = 0.5, vjust = 0.5))

h <- g + annotation_custom(
  grob = textGrob(label = "The infamous mtcars dataset", rot = -45,
                  x = unit(1.1, "npc"), y = unit(1.1, "npc")))



g

g + coord_flip()

g + scale_y_reverse()

h = h + coord_flip()
h = h + scale_y_reverse()
print(h, vp = viewport(width = unit(0.5, "npc"), height = unit(0.5, "npc"), angle = 45))


grid.draw(rectGrob())
sample_vp <- viewport(x = 0, y = 0, 
                      width = 0.2, height = 1,
                      just = c("left", "bottom"))
pushViewport(sample_vp)
grid.draw(roundrectGrob())
grid.draw(h)
popViewport()

grid.draw(rectGrob())
sample_vp <- viewport(x = 0.5, y = 0.5, 
                      width = 0.5, height = 0.5,
                      just = c("left", "bottom"))
pushViewport(sample_vp)
grid.draw(roundrectGrob())
grid.draw(g)
popViewport()



#### KP objects


kp <- plotKaryotype(plot.type = 1, chromosomes = c("chr1", "chr2", "chr3", "chr4"))


#### Trying record Plot

x <- recordPlot(load=NULL, attach=NULL)
y <- grob(replayPlot(x, reloadPkgs=FALSE))

############################################################################################
### Write to folder
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates/')

### Save plot
png(filename= 'test.png', width = 2000, height = 1600, units = "px", pointsize = 12, res = 250)
x
dev.off()

### Oh shit that worked! Now how do I make it a grob?

library("ggplotify")

p1 <- as.grob(~barplot(1:10))

p2 <- as.grob(x)

p3 <- as.grob(replayPlot(x, reloadPkgs = FALSE))

p4 <- as.grob(y)

library(cowplot)

p5 <- as_grob(x)

grid.arrange(g, h, p5, nrow = 1)

### SHIT! It worked AGAIN!!!
### Put it all together...



kp <- plotKaryotype(plot.type = 1, chromosomes = c("chr1", "chr2", "chr3", "chr4"))

### Data Panel 1 ###

stackPos = 0
width = 0.2 # 10
width = 0.1 # 20. pretty good visibility
width = 0.05 #Around 40? 1/2 pixel limit approached, even with huge distances

#Big regions
kpRect(kp, data = big.regs.up, y0=stackPos, y1=stackPos+width, col="#FFDDDD", border=NA, r0=0, r1=0.8)
stackPos = stackPos + width
kpRect(kp, data = big.regs.down, y0=stackPos, y1=stackPos+width, col="#DDFFDD", border=NA, r0=0, r1=0.8)
stackPos = stackPos + width
kpRect(kp, data = big.regs.three, y0=stackPos, y1=stackPos+width, col="#FF3366", border=NA, r0=0, r1=0.8)
stackPos = stackPos + width
kpRect(kp, data = big.regs.four, y0=stackPos, y1=stackPos+width, col="#C200FB", border=NA, r0=0, r1=0.8)
stackPos = stackPos + width
kpRect(kp, data = big.regs.up, y0=stackPos, y1=stackPos+width, col="#FFDDDD", border=NA, r0=0, r1=0.8)
stackPos = stackPos + width
kpRect(kp, data = big.regs.down, y0=stackPos, y1=stackPos+width, col="#C200FB", border=NA, r0=0, r1=0.8)
stackPos = stackPos + width
kpRect(kp, data = big.regs.four, y0=stackPos, y1=stackPos+width, col="#DDFFDD", border=NA, r0=0, r1=0.8)
stackPos = stackPos + width
kpRect(kp, data = big.regs.three, y0=stackPos, y1=stackPos+width, col="#FFDDDD", border=NA, r0=0, r1=0.8)
stackPos = stackPos + width
x <- recordPlot(load=NULL, attach=NULL)
p5 <- as_grob(x)
grid.arrange(g, p5, nrow = 1)

### Now rotate.

g + coord_flip()
grid.arrange(g, p5, nrow = 1)


g + scale_y_reverse()

print(h, vp = viewport(width = unit(0.5, "npc"), height = unit(0.5, "npc"), angle = 45))
print(p5, vp = viewport(angle = 45))
print(p5)

p5
plot(p5)


p6 <- editGrob(p5, vp=viewport(angle=-90, width=unit(0.85,"npc"), 
                                     height=unit(0.85,"npc")))    

grid.arrange(g, p6, nrow = 1)
