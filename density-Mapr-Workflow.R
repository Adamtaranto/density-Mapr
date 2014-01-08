#############################################
## Check environment
getwd()
ls()

#############################################
## 2.1	Load required R packages

source("http://bioconductor.org/biocLite.R")
biocLite()

require("GenomicRanges")
require("rtracklayer")
require("Rsamtools")
require("png")
require("gridExtra")
require("ggplot2")


#############################################
## 2.2	Required custom R functions

source("Resources/filled.contour3.R")
source("Resources/getFeat2b.R")


#############################################
## 2.3 Optional R packages

require("fields")
require("rgl")
require("EBImage")


#############################################
## 3.1 Calculate flanking intergenic regions (FIRs)

gff <- import.gff("Sample_datasets/Mygtf.gtf", asRangedData=FALSE)

gffgene<-getFeat2b(x=gff, format="gtf", range_types=c("gene"))

strand(gffgene)<-mcols(gffgene)$score

mcols(gffgene)$score<-NULL

gffintg<-getFeat2b(x=gff, format="gtf", range_types=c("intergenic"))

length_intg<-as.data.frame(cbind(seq(1:length(ranges(gffintg))), as.numeric(mcols(gffintg)$length)))

colnames(length_intg)<-c("index", "length")

three_intg_index<-precede(gffgene, gffintg)

five_intg_index<-follow(gffgene, gffintg)

gene_data<-   as.data.frame(cbind(as.character(mcols(gffgene)$group), 
              as.character(strand(gffgene)), as.numeric(five_intg_index), 
              as.numeric(three_intg_index)))

colnames(gene_data)<-c("geneid", "strand", "FivePrime_index", "ThreePrime_index")

tempdata<-merge(x=gene_data, y=length_intg, by.x="FivePrime_index", by.y="index", all.x=TRUE)

colnames(tempdata)<-c("delete1", "geneid", "strand", "ThreePrime_index", "fiveprime")

FIRdata<-merge(x=tempdata, y=length_intg, by.x="ThreePrime_index", by.y="index", all.x=TRUE)

FIRdata$ThreePrime_index<-NULL

FIRdata$delete1<-NULL

colnames(FIRdata)<-c("geneid", "strand", "fiveprime", "threeprime")

## Export current FIR data
write.table(FIRdata,file="MyFIRs.csv", sep=",", row.names=FALSE)


#############################################
## Alt. import own FIR data

FIRdata<-read.csv(file="Resources/MyFIRs.csv",sep=",")


#############################################
## 3.2  Bin breaks setup

NumBins=40

if ((max(FIRdata$fiveprime, na.rm=TRUE)>max(FIRdata$threeprime, na.rm=TRUE)) == TRUE) {
       FIR2Bin<-FIRdata$fiveprime
} else {
       FIR2Bin<-FIRdata$threeprime
       }


FIR2Bin=FIR2Bin[which(FIR2Bin!=0)]

FIR2Bin<-na.omit(FIR2Bin)

BinSteps<-round(length(FIR2Bin)/(NumBins-1), digits=0)

FIR2BinOrd<-sort(FIR2Bin)

TempBinLimits<-FIR2BinOrd[seq(FIR2BinOrd[2*BinSteps], length(FIR2BinOrd),BinSteps)]

TempBinLimits[length(TempBinLimits)+1]<-max(FIR2Bin, na.rm=TRUE)

x<-seq(length(TempBinLimits))

fit<-nls(log(TempBinLimits) ~ a*x + b, start=c(a=0, b=0), 
       algorithm='port', weights=((x-0.5*NumBins)^2))

pred=predict(fit, x)

BinLimits=c(1, round(exp(pred),0), max(FIR2Bin))

## Export current BinLimits to file
write.table(BinLimits,file="MyBins.txt")

#############################################
## Alt. import an external set of bin breaks:

BinLimits<-as.numeric(unlist(read.table(file="Resources/MyBins.txt", header=TRUE, row.names=1)))

#############################################
## 3.3 Data binning

xbin=cut(FIRdata$fiveprime, breaks= c(BinLimits))

ybin=cut(FIRdata$threeprime, breaks= c(BinLimits))

FIRdata<-cbind(FIRdata, xbin, ybin, 
       genevalue=rep(1, length(FIRdata$fiveprime)))

GenValMatrix<-with(FIRdata, tapply(genevalue, list(xbin, ybin), sum))

#############################################
## 3.4 Heatmap drawing
## Start from here if importing both FIRdata and BinLimits

x<-1:ncol(GenValMatrix)

y<-1:nrow(GenValMatrix)

zlim = range(as.numeric(unlist(GenValMatrix)), finite=TRUE)

mypalette<-colorRampPalette(c( "white", "darkblue", "forestgreen", "goldenrod1", "orangered", "red3", "darkred"), space="rgb")

mycol=mypalette(2*max(GenValMatrix, na.rm=TRUE))

mylabels<-paste(BinLimits[1:length(BinLimits)-1], BinLimits[2:length(BinLimits)], sep = " - ", collapse = NULL)

filled.contour(x, y, z=GenValMatrix, 
       plot.title = title(main ="Phytophthora infestans genome", 
              xlab = "five prime intergenic regions", 
              ylab = "three prime intergenic regions", 
              cex.main=0.8, cex.lab=0.5), 
       key.title = title(main ="Number of genes", cex.main=0.5, line=1), 
       col=mycol, 
       levels = pretty(zlim, 2*max(GenValMatrix, na.rm=TRUE)), 
       plot.axes={axis(1,at=x, labels=mylabels, las=2, 
              cex.axis=0.5);
       axis(2,at=y, labels=mylabels, cex.axis=0.5)})


#############################################
## 3.5 Overlaying a scatter plot over a genome architecture heatmap
image_name<-paste(as.character(format(Sys.time(), "%Y%m%d%H%M%S")), "_graph", sep="")

png(filename = paste(image_name, ".png", sep=""))

par(mar=c(0,0,0,0))

filled.contour3(x, y, z=GenValMatrix, 
       col=mycol, 
       levels = pretty(zlim, 2*max(GenValMatrix, na.rm=TRUE)), 
       frame.plot = FALSE, 
       axes = FALSE)

dev.off()

## Save current graphic as pdf
quartz.save("myHeatmap.pdf", type="pdf")


## Save Heatmap as png, import as raster for use as ggplot background
img <- readPNG(paste(image_name, ".png", sep=""))

g <- rasterGrob(img, interpolate=TRUE)

## Import data points
rxlrData<-as.data.frame(read.csv('Sample_datasets/RXLR_FIRs.csv', header=TRUE))

## Map data points onto base heatmap


#Loose 39 rxlr records: 38 with at least one FIR > 64843 (the value of the 40th bin in BinLimits), and 1 which is below the 2 value in BinLimits (20).
gg<-ggplot(data=rxlrData, aes(x=rxlrData$rxlr_five, y=rxlrData$rxlr_three, geom="blank")) + 
       annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
       coord_fixed(ratio=1) + 
       geom_point(shape=21, fill="white", colour="black", size=4, alpha=0.7, na.rm=FALSE) + 
       scale_y_log10(breaks= BinLimits[2:length(BinLimits)], limits= c(BinLimits[2], BinLimits[NumBins])) + 
       scale_x_log10(breaks= BinLimits[2:length(BinLimits)], limits= c(BinLimits[2], BinLimits[NumBins])) + 
       theme(axis.text.y=element_text(size =8, vjust=0.5)) + 
       theme(axis.text.x=element_text(size=8, vjust=0.5, angle=90)) + 
       theme(axis.title.x = element_text(face="bold",size=12)) +
       xlab("five prime intergenic region") + 
       theme(axis.title.y = element_text(face="bold",size=12)) +
       ylab("three prime intergenic region")

#Correct NumBins to equal actual number of bins in BinLimits
NumBins<-as.numeric(c(length(BinLimits)))

#Loose 1 record, (190: 39629,12)
gg<-ggplot(data=rxlrData, aes(x=rxlrData$rxlr_five, y=rxlrData$rxlr_three, geom="blank")) + 
       annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
       coord_fixed(ratio=1) + 
       geom_point(shape=21, fill="white", colour="black", size=4, alpha=0.7, na.rm=FALSE) + 
       scale_y_log10(breaks= BinLimits[2:length(BinLimits)], limits= c(BinLimits[2], BinLimits[NumBins])) + 
       scale_x_log10(breaks= BinLimits[2:length(BinLimits)], limits= c(BinLimits[2], BinLimits[NumBins])) + 
       theme(axis.text.y=element_text(size =8, vjust=0.5)) + 
       theme(axis.text.x=element_text(size=8, vjust=0.5, angle=90)) + 
       theme(axis.title.x = element_text(face="bold",size=12)) +
       xlab("five prime intergenic region") + 
       theme(axis.title.y = element_text(face="bold",size=12)) +
       ylab("three prime intergenic region")

#Loose no data points. But axes are messed up, and heatmap is incorrectly scaled.
gg<-ggplot(data=rxlrData, aes(x=rxlrData$rxlr_five, y=rxlrData$rxlr_three, geom="blank")) + 
       annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
       coord_fixed(ratio=1) + 
       geom_point(shape=21, fill="white", colour="black", size=4, alpha=0.7, na.rm=FALSE) + 
       scale_y_log10(breaks= BinLimits[2:length(BinLimits)], limits= c(BinLimits[1], BinLimits[NumBins])) + 
       scale_x_log10(breaks= BinLimits[2:length(BinLimits)], limits= c(BinLimits[1], BinLimits[NumBins])) + 
       theme(axis.text.y=element_text(size =8, vjust=0.5)) + 
       theme(axis.text.x=element_text(size=8, vjust=0.5, angle=90)) + 
       theme(axis.title.x = element_text(face="bold",size=12)) +
       xlab("five prime intergenic region") + 
       theme(axis.title.y = element_text(face="bold",size=12)) +
       ylab("three prime intergenic region")

gg

quartz.save("myHeatmap_withPoints.pdf", type="pdf")

