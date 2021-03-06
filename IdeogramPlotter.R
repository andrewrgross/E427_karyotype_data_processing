### Ideogram generator -- Andrew R Gross -- 2021-08-02
### Third generation ideogram generator

###############################################################################################
### 1 - Header ################################################################################
library(karyoploteR)
library(regioneR)
library(zoo)
library(biovizBase)
library(ggplot2)
library(ggplotify)
library(scatterpie)
library(RIdeogram)
library(IdeoViz)

###### 1.1 - Functions ########################################################################
cytogeneticTableCheck <- function(dataframe) {
  for (currentRowNumber in 1:nrow(dataframe)){
    errorCounter = 0
    currentRow = dataframe[currentRowNumber,]
    chr <- currentRow[,1]
    startPos = currentRow[,2]
    endPos = currentRow[,3]
    type = currentRow[,4]
    ### Check if type is recognized
    if (all(grepl(type, c('tri', 'dup', 'mono', 'del', 'trans', 'inv', 'blank')))) {
      print(paste('Row', currentRowNumber, '; Type ', type, 'Not found.')); errorCounter = errorCounter + 1
    }
    if (currentRow[,4] != 'blank') {
      ### Subset by chromosome
      chrGRange <- hg19[seqnames(hg19) == chr]
      ### Search for the position on the chromosome
      if (all(grepl(startPos, chrGRange$name)==FALSE)) {
        if(startPos != '-') {
          print(paste('Row', currentRowNumber, '; ', chr, ', ', startPos, 'Not found.')); errorCounter = errorCounter + 1
        } 
      }
      if (all(grepl(endPos, chrGRange$name)==FALSE)) {
        if(endPos != '-') {
          print(paste('Row', currentRowNumber, '; ', chr, ', ', startPos, 'Not found.')); errorCounter = errorCounter + 1
        }
      } # End pos check loop
    } # Blank check loop
  } # for loop
  if(errorCounter == 0) {print('Check complete. All IDs recognized.')}
}

cytoToCoord <- function(dataframe, verbose = FALSE) {
  for (currentRowNumber in 1:nrow(dataframe)){
    if(verbose){print(currentRowNumber)}
    currentRow = dataframe[currentRowNumber,]
    chr <- currentRow[,1]
    startPos = currentRow[,2]
    endPos = currentRow[,3]
    type = currentRow[,4]
    
    colorNum = 7
    if(type == 'tri'){colorNum = 1}
    if(type == 'dup'){colorNum = 2}
    if(type == 'mono'){colorNum = 3}
    if(type == 'del'){colorNum = 4}
    if(type == 'trans'){colorNum = 5}
    if(type == 'inv'){colorNum = 6}
    ### Subset by chromosome
    chrGRange <- hg19[seqnames(hg19) == chr]
    ### Find start and end coordinates of chr
    firstOfChr = chrGRange[1]
    lastOfChr  = chrGRange[length(chrGRange)]
    
    ### Check if the row contains an ID
    if(type == 'blank') {         
      startCoord = 1  ;  endCoord = 1                  # If blank, bar has zero length
      
      ### Find the start position
    } else {
      if (startPos == '-'){                            # If the position is '-', the start is the beginning of the chromosome
        startCoord = 1
      } else {
        targetRowOnChr= grep(startPos, chrGRange$name) # If a position is specified, find it
        startCoord    = start(ranges(chrGRange[targetRowOnChr])[1])
      }
      if(verbose){print(paste('Cytogenic name =', startPos, ';  StartCoord =', startCoord))}
      
      ### Find end position
      if (endPos == '-'){
        ### If the end position is a dash, check which branch the start is on
        if(startPos == '-'){
          endCoord = end(ranges(lastOfChr))
        }
        if (grepl('p', startPos)) {
          endCoord = 1
        }
        if (grepl('q', startPos)) {
          endCoord = end(ranges(lastOfChr))
        }
      }else {                                          # If the end position is specifed, find it
        targetRowOnChr= grep(endPos, chrGRange$name)
        endCoord      = end(ranges(chrGRange[targetRowOnChr])[1])
      }
    }
    if(verbose){print(paste('Cytogenic name =', endPos, ';  endCoord =', endCoord))}
    ### If the starcoordinate is greater than end, invert them.
    if (startCoord > endCoord) {temp = startCoord; startCoord = endCoord; endCoord = temp}  ### IRanges must have the start before the end, so if that's not the case, flip'em
    ### Assign the coordinates to the input data frame
    dataframe$startCoord[currentRowNumber] <- startCoord
    dataframe$endCoord[currentRowNumber]   <- endCoord
    dataframe$color[currentRowNumber]      <- colorNum
  }  # End of row loop
  return(dataframe)
}

###############################################################################################
### 2 - Import Data ###########################################################################
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates/')
expandedRaw <- read.table('Expanded_abnorms.txt', sep = '\t')
unexpandedRaw <- read.table('Unexpanded_abnorms.txt', sep = '\t')

hg19 <- getIdeo("hg19")

###############################################################################################
### 3 - Format Data ###########################################################################
###### 3.1 - Subset data as necessary #####################
unexpanded1 <- unexpandedRaw[c(1:25),]         # Subset chr1 - chr9
unexpanded2 <- unexpandedRaw[c(1,24:59),]      # Subset chr1, chr9 - chr16
unexpanded3 <- unexpandedRaw[c(1,60:117),]     # Subset chr1, chr17 = chrY

expanded1 <- expandedRaw[1:50,]                # Subset chr1 - chr9
expanded2 <- expandedRaw[c(1,50:105),]         # Subset chr1, chr9 - chr16
expanded214 <- expandedRaw[c(1,50:83, 87:105),]# Subset chr1, chr9 - chr16, but hide some of 13 so 14 is visible
expanded3 <- expandedRaw[c(1,106:158),]        # Subset chr1, chr17 - chrY
expanded3y <- expandedRaw[c(1,106:136,154:158 ),]# Subset chr1, chr17 - chrY, but hide some of chrX so Y is visible


###### 3.1 - Check if each cytogenetic position's coordinate can be found #####################
cytogeneticTableCheck(unexpanded3)
cytogeneticTableCheck(expanded3)

hg19[which(seqnames(hg19) == 'chr11')][1:10]

###### 3.2 - Convert cytogenetic notation to coordinates ######################################

(dataToPlot <- cytoToCoord(unexpanded3, verbose = FALSE))

(dataToPlot <- cytoToCoord(expanded214, verbose = FALSE))

###### 3.3 - Generate GRange object to plot ###################################################

### Without stack position Assignment
grToPlot <- GRanges(
  seqnames = dataToPlot$V1,
  ranges = IRanges(start = dataToPlot$startCoord, end = dataToPlot$endCoord),
  color = dataToPlot$color)


###############################################################################################
### 4 - Plot data #############################################################################
###### 4.1 - Alternatively, I could just reformat manually ###
colorList <- c("#006837" , "#1A9850" , "#A50026" , "#D73027" , "#4575B4" , "#313695" , "#b8b8b8" )  # Final choice

width = 0.1                  # For expanded
#width = 0.1                  # For expanded
#width = 0.01421333            # For unexpanded
width = 533/1875 * 0.1            # For unexpanded
verbose = FALSE
### Optimizing for output:
### Display bars in data.panel 2
pp <- getDefaultPlotParams(plot.type=2)
pp$ideogramheight = 25
pp$data1inmargin = 0
pp$data2inmargin = 10
pp$data1height = 0
pp$data2height = 300
pp$leftmargin = 0.3
kp <- plotKaryotype(plot.type = 2, chromosomes = unique(seqnames(grToPlot)), plot.params = pp)
for(currentChr in unique(seqnames(grToPlot))) {
  if(verbose){print(currentChr)}
  currentGR <- grToPlot[which(as.vector(seqnames(grToPlot) == currentChr))]
  ### For a given chromosome, plot each entry in the input GRange
  if(verbose){print(currentGR)}
  stackPos = 0
  for (rangeNum in 1:length(currentGR)) {
    range = currentGR[rangeNum]
    color = colorList[range$color]
    kpRect(kp, data = range, data.panel = 2, y0=stackPos, y1=stackPos+width, col= color, border=NA, r0=0, r1=0.8)
    #kpRect(kp, data = range, y0=(range$stackPosA*width), y1=range$stackPosA*width+width, col= color, border=NA, r0=0, r1=0.8)
    stackPos = stackPos + width
  }
}


###### 4.2 - Plot pie charts ###
colorList <- c("#006837" , "#1A9850" , "#A50026" , "#D73027" , "#4575B4" , "#313695" , "#b8b8b8" )  # Final choice

pieData <- data.frame(abnormal = 103, normal = 430, total = 533, percent = 19.3)
pieData <- data.frame(Group = c('Abnormal minux sc', 'Abnormal in sc', 'Normal'), value = c(92, 11, 430))

piechart <- ggplot2

# Basic piechart
ggplot(pieData, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)

sizeData <- data.frame(x = c(5, 30, 100, 150), size = c(1, 533, 1874, 20))

ggplot(sizeData, aes(x = x, y = 20, size = size)) + 
  geom_point() + scale_size_area() + scale_size(range = c(0,150)) + theme(legend.position = 'none')

sizeData2 <- data.frame(x = c(0, 3000), y = c(20, 20), A = c(15,53), B = c(80,42), C = c(427, 1782), r = c(522,1877))
sizeColor <- c('grey40', 'grey60', 'grey90')
sizeColor <- c('#54278f', '#756bb1', '#f2f0f7')
sizeColor <- c('#54278f', '#9e9ac8', '#f2f0f7')


ggplot(data = sizeData) + 
  geom_scatterpie(data = sizeData2, aes(x=x, y=y, r=sqrt(r/pi)*70), cols=LETTERS[1:3], size = 1.4) + 
  coord_equal() + scale_fill_manual(values = sizeColor, labels = c('Lines with abnormalities (somatic)', 'Lines with abnormalities (on X/Y)', 'Normal Lines')) +
  theme(panel.background = element_rect(fill = 'white'), axis.text = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(),
        legend.key.size = unit(6,'mm'), legend.title = element_blank(), legend.text = element_text(size = 23), legend.position = 'None')


###### 4.3 - Plot good chromosome images ###
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")
gene_density <- GFFex(input = "gencode.v32.annotation.gff3.gz", karyotype = "human_karyotype.txt", feature = "gene", window = 1000000)

human_karyotype <- human_karyotype[,1:3]
ideogram(karyotype = human_karyotype, overlaid = NULL, label = NULL, label_type = NULL)
convertSVG("chromosome.svg", device = "png")


ideogram(karyotype, overlaid = NULL, label = NULL, label_type = NULL, synteny = NULL, colorset1, colorset2, width, Lx, Ly, output = "chromosome.svg")
convertSVG(svg, device, width, height, dpi)
convertSVG("chromosome.svg", device = "tiff", dpi = 600)

require(IdeoViz)
data(binned_multiSeries)
hg19 <- getIdeo("hg19")

data(binned_singleSeries)
chr = c('chr1', 'chr2',  'chr3',  'chr4',  'chr5',  'chr6',  'chr7',  'chr8',  'chr9')
chr = c('chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17')
chr = c('chr1', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',  'chrY')
df <- data.frame(chr=chr, start=1, end=2, values = 0)
gro <- GRanges(df)

plotOnIdeo(chrom=seqlevels(gro), ideo=hg19, values_GR = gro, value_cols= colnames(mcols(gro)), vertical=T, addScale = F)

data(hg18_ideo) # cytoBandIdeo table downloaded previously and stored as a data.frame.
plotOnIdeo(chrom=chr, ideo=hg19, values_GR=binned_singleSeries, value_cols=colnames(mcols(binned_singleSeries)), vertical=T,  
             val_range=c(-10,10))

plotOnIdeo(chrom=chr, ideo=hg19, values_GR=binned_singleSeries, value_cols=colnames(mcols(binned_singleSeries)), vertical=T)

############################################################################################
### Write to folder
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates/Output images/')

kpOutput <- recordPlot(load=NULL, attach=NULL)

### Save plot
tiff(filename= 'expanded1.tiff', width = 400, height = 2000, units = "px", pointsize = 12, res = 250)
kpOutput
dev.off()



p <- ggplot() + geom_scatterpie(aes(x=long, y=lat, group=region, r=radius), data=d,
                                cols=LETTERS[1:4], color=NA) + coord_equal()






########################################################## ##########################################################
########################################################## ##########################################################
########################################################## ##########################################################
##### Scratchwork

############################################################################################
### Write to folder

test <- as_grob(kpOutput)
test <- as.ggplot(expression(
  for(currentChr in unique(seqnames(grToPlot))) {
    if(verbose){print(currentChr)}
    currentGR <- grToPlot[which(as.vector(seqnames(grToPlot) == currentChr))]
    
    ### For a given chromosome, plot each entry in the input GRange
    pp <- getDefaultPlotParams(plot.type=2)
    pp$rightmargin <- margins[marginNum]
    kp <- plotKaryotype(plot.type = 1, chromosomes = currentChr, plot.params = pp)
    if(verbose){print(currentGR)}
    stackPos = 0
    for (rangeNum in 1:length(currentGR)) {
      range = currentGR[rangeNum]
      color = colorList[range$color]
      kpRect(kp, data = range, y0=stackPos, y1=stackPos+width, col= color, border=NA, r0=0, r1=0.8)
      stackPos = stackPos + width
    }
    ### Save the complete plot of the current chromosome
    latestPlot <- recordPlot(load=NULL, attach=NULL)
    latestPlot <- as_grob(latestPlot)
    latestPlot <- editGrob(latestPlot, vp=viewport(angle=-90, width=unit(0.85,"npc"), height=unit(0.85,"npc")))  
    plotList[[currentChr]] <- latestPlot
    marginNum = marginNum + 1
    print(paste(currentChr, 'Complete'))
  }
))

grid.arrange(test, nrow = 1)

library(cowplot)
p1 <- as.ggplot(expression(kp <- plotKaryotype(main="Human (hg19)", plot.params = pp), kpDataBackground(kp), kpAbline(kp, h=0.5, col="red")))
p2 <- as.ggplot(expression(kp <- plotKaryotype(genome="hg38", main="Human (hg38)",  plot.params = pp), kpDataBackground(kp), kpAbline(kp, h=0.5, col="blue")))
p3 <- as.ggplot(expression(plotKaryotype(genome = "mm10", main="Mouse", plot.params = pp)))
p4 <- as.ggplot(expression(plotKaryotype(genome = "dm6", main="Fruit fly")))
save_plot("Multipanel_cowplot.png", plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:9]), base_height=10)

grid.arrange(plotList[[1]], plotList[[2]], plotList[[3]],  nrow = 1)
grid.arrange(latestPlot, nrow = 1)
grid.arrange(plotList[[1]], nrow = 1)

###### 3.2.1 - Assign stack position ##########################################################
### For 8-14
dataToPlot$StackPosA = c(0,
                         0,1,
                         0,
                         0,1,2,3,4,4,0,5,6,7,8,9,10,1,11,12,
                         0,1,2,
                         0,1,1,2,3,
                         0)

### WITH stack position assingnment
grToPlot <- GRanges(
  seqnames = dataToPlot$V1,
  ranges = IRanges(start = dataToPlot$startCoord, end = dataToPlot$endCoord),
  color = dataToPlot$color,
  stackPosA = dataToPlot$StackPosA)



###### 4.2 - Loop through GRanges object and generate separate ideogram for each chromosome ###
colorList <- c("#006837", "#A50026", "#1A9850", "#D73027", "#313695", "#4575B4", "#b8b8b8")  # Final choice
margins <- c(0,0.1,.30,.40,.50,.60,.70,.80,.90,.100)
marginNum = 1

plotList <- list()
width = 0.1
for(currentChr in unique(seqnames(grToPlot))) {
  if(verbose){print(currentChr)}
  currentGR <- grToPlot[which(as.vector(seqnames(grToPlot) == currentChr))]
  
  ### For a given chromosome, plot each entry in the input GRange
  pp <- getDefaultPlotParams(plot.type=2)
  pp$rightmargin <- margins[marginNum]
  kp <- plotKaryotype(plot.type = 1, chromosomes = currentChr, plot.params = pp)
  if(verbose){print(currentGR)}
  stackPos = 0
  for (rangeNum in 1:length(currentGR)) {
    range = currentGR[rangeNum]
    color = colorList[range$color]
    kpRect(kp, data = range, y0=stackPos, y1=stackPos+width, col= color, border=NA, r0=0, r1=0.8)
    stackPos = stackPos + width
  }
  ### Save the complete plot of the current chromosome
  latestPlot <- recordPlot(load=NULL, attach=NULL)
  latestPlot <- as_grob(latestPlot)
  latestPlot <- editGrob(latestPlot, vp=viewport(angle=-90, width=unit(0.85,"npc"), height=unit(0.85,"npc")))  
  plotList[[currentChr]] <- latestPlot
  marginNum = marginNum + 1
  print(paste(currentChr, 'Complete'))
}


