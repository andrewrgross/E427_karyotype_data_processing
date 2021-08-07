### Ideogram generator -- Andrew R Gross -- 2021-08-02
### Third generation ideogram generator

library(karyoploteR)
library(regioneR)
library(zoo)

set.seed(1234)

#Parameters
data.points.colors <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")

num.data.points <- 3000
num.big.regions.up <- 3
num.big.regions.down <- 30

num.mid.regions <- 6000

num.marks <- 90

#Create the random fake data  

#Big regions

big.regs.up <- joinRegions(
  createRandomRegions(nregions = num.big.regions.up, length.mean = 20000000, length.sd = 10000000, non.overlapping = TRUE, mask=NA), min.dist = 1)

big.regs.down <- joinRegions(createRandomRegions(nregions = num.big.regions.down, length.mean = 10000000, length.sd = 5000000, non.overlapping = TRUE, mask=big.regs.up), min.dist = 1)

big.regs.three <- joinRegions(createRandomRegions(nregions = num.big.regions.down, length.mean = 10000000, length.sd = 5000000, non.overlapping = TRUE, mask=big.regs.up), min.dist = 1)
big.regs.four <- joinRegions(createRandomRegions(nregions = num.big.regions.down, length.mean = 10000000, length.sd = 5000000, non.overlapping = TRUE, mask=big.regs.up), min.dist = 1)


#Data points
data.points <- createRandomRegions(nregions = num.data.points, length.mean = 1, length.sd = 0, non.overlapping = TRUE, mask=NA)
mcols(data.points) <- data.frame(y=rnorm(n = num.data.points, 0.5, sd = 0.1))
dp.colors <- sample(head(data.points.colors, 2), size = num.data.points, replace = TRUE)

#and move the data points with the big regions
data.points[overlapsAny(data.points, big.regs.up)]$y <- data.points[overlapsAny(data.points, big.regs.up)]$y + runif(n=numOverlaps(data.points, big.regs.up), min = 0.1, max=0.3)
data.points[overlapsAny(data.points, big.regs.down)]$y <- data.points[overlapsAny(data.points, big.regs.down)]$y - runif(n=numOverlaps(data.points, big.regs.down), min = 0.1, max=0.3)

#markers
marks <- createRandomRegions(nregions = num.marks, length.mean = 1, length.sd = 0)
mcols(marks) <- data.frame(labels=paste0("rs", floor(runif(num.marks, min = 10000, max=99999))))

#medium regions
mid.regs <- createRandomRegions(nregions = num.mid.regions, length.mean = 5000000, length.sd = 1000000, non.overlapping = FALSE)

######################################################################################

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

kpPlotRegions(kp, data = mid.regs, r0 = 0.2, r1=1, border=NA, data.panel=1)
kpPlotRegions(kp, data = mid.regs, r0 = 0, r1=1, border=NA, data.panel=1, layer.margin = 0)

kpPlotRegions(kp, data = big.regs.up, r0 = 0, r1=0.2, layer.margin = 0, border=NA, data.panel=1)
kpPlotRegions(kp, data = big.regs.down, r0 = 0.3, r1=1, border=NA, data.panel=1)
kpPlotRegions(kp, data = big.regs.three, r0 = 0.3, r1=1, border=NA, data.panel=1, col = '#FFDDDD')
kpPlotRegions(kp, data = big.regs.four, r0 = 0.3, r1=0.5, border=NA, data.panel=1, col = '#FFDDDD')


#Data points
kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.8, numticks = 5, col="#666666", cex=0.5)
kpPoints(kp, data=data.points, pch=16, cex=0.5, col=dp.colors, r0=0, r1=0.8)

#Mean and sd of the data points.  
for(chr in seqlevels(kp$genome)) {
  chr.dp <- sort(keepSeqlevels(x = data.points, value = chr, pruning.mode = "coarse"))
  rmean <- rollmean(chr.dp$y, k = 6, align = "center")  
  rsd <- rollapply(data = chr.dp$y, FUN=sd, width=6)
  kpLines(kp, chr = chr, x=start(chr.dp)[3:(length(chr.dp)-3)], y=rmean, col=data.points.colors[3], r0=0, r1=0.8)
  kpPlotRibbon(kp, chr=chr, data=chr.dp[3:(length(chr.dp)-3)], y0=rmean-rsd, y1=rmean+rsd, r0=0, r1=0.8, col="#FF336633", border=NA)
}

#Markers
kpPlotMarkers(kp, data=marks, label.color = "#333333", r1=1.1, cex=0.5, label.margin = 5)

### Data Panel 2 ###

#medium regions and their coverage

kpPlotRegions(kp, data = mid.regs, r0 = 0.2, r1=1, border=NA, data.panel=2)
kpPlotCoverage(kp, data=mid.regs, r0=0.2, r1=0, col=data.points.colors[2], data.panel = 2)
kpPlotCoverage(kp, data=mid.regs, r0=0.2, r1=0.12, col=data.points.colors[1], data.panel = 2)

kpText(kp, chr=seqlevels(kp$genome), y=0.4, x=0, data.panel = 2, r0=0.2, r1=0, col="#444444", label="30x", cex=0.8, pos=2)
kpAbline(kp, h=0.4, data.panel = 2, r0=0.2, r1=0, col=data.points.colors[3])


data.points.colors <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")

kp <- plotKaryotype(plot.type = 1, chromosomes = c("chr1", "chr2", "chr3", "chr4"))
getCytobandColors(color.table=data.points.colors , color.schema=c("circos", "biovizbase", "only.centromeres"))
colors = getCytobandColors(color.table=data.points.colors , color.schema="biovizbase")

kpAddCytobandsAsLine(kp, lwd = 16, lend = 0)

kp <- plotKaryotype(ideogram.plotter = NULL, plot.type=2)
kpAddCytobandsAsLine(kp)
