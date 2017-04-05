#############################################
## Check environment
getwd()
ls()

#############################################
## 1.0 Install required R packages

## Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("rtracklayer")
biocLite("Rsamtools")

## CRAN packages
install.packages("png")
install.packages("gridExtra")
install.packages("ggplot2")

## Optional
install.packages("fields")
install.packages("rgl")
install.packages("EBImage")

## Other
#require(devtools)
#install_version("ggplot2", version = "0.9.1", repos = "http://cran.us.r-project.org")

## 1.1 Import packages
library("GenomicRanges","rtracklayer","Rsamtools",
       "png","gridExtra","ggplot2")

#############################################
## 1.2 Source custom R functions
source("src/filled.contour3.R")
source("src/getFeat2b.R")

#############################################
## 1.3 Optional R packages

library("fields","rgl","EBImage")

#############################################
## 2.0 Calculate flanking intergenic regions (FIRs)

x <- import.gff("Sample_datasets/Mygtf.gtf")

## Set "group_col" to gff metadata field with unique gene identifiers.
gffgene <- getFeat2b(x=x, format="GFF", group_col="gene_id", range_types=c("gene"))

strand(gffgene) <- mcols(gffgene)$score

mcols(gffgene)$score<-NULL

gffintg<-getFeat2b(x=x, format="GFF",group_col="gene_id", range_types=c("intergenic"))

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
## 3.0  Bin breaks setup

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
## 3.1 Data binning

xbin=cut(FIRdata$fiveprime, breaks= c(BinLimits))
ybin=cut(FIRdata$threeprime, breaks= c(BinLimits))

FIRdata<-cbind(FIRdata, xbin, ybin, 
       genevalue=rep(1, length(FIRdata$fiveprime)))

GenValMatrix<-with(FIRdata, tapply(genevalue, list(xbin, ybin), sum))

#############################################
## Import pre-calculated bins and FIRS
#############################################

#4.0
## Alt. import an external set of bin breaks:
BinLimits<-as.numeric(unlist(read.table(file="Resources/MyBins.txt", header=TRUE, row.names=1)))

## Alt. import own FIR data
FIRdata<-read.csv(file="Resources/MyFIRs.csv",sep=",")
xbin=cut(FIRdata$fiveprime, breaks= c(BinLimits))
ybin=cut(FIRdata$threeprime, breaks= c(BinLimits))
FIRdata<-cbind(FIRdata, xbin, ybin, 
       genevalue=rep(1, length(FIRdata$fiveprime)))
GenValMatrix<-with(FIRdata, tapply(genevalue, list(xbin, ybin), sum))

## 4.1 Heatmap drawing
## Start from here if importing both FIRdata and BinLimits

x<-1:ncol(GenValMatrix)
y<-1:nrow(GenValMatrix)

zlim = range(as.numeric(unlist(GenValMatrix)), finite=TRUE)

mypalette <- colorRampPalette(c( "white", 
                                   "darkblue", 
                                   "forestgreen", 
                                   "goldenrod1", 
                                   "orangered", 
                                   "red3", 
                                   "darkred"), 
                                   space="rgb")

mycol=mypalette(2*max(GenValMatrix, na.rm=TRUE))

mylabels <- paste(BinLimits[1:length(BinLimits)-1], 
              BinLimits[2:length(BinLimits)], 
              sep = " - ", collapse = NULL)

filled.contour(x, y, z=GenValMatrix, 
       plot.title = title(main ="Phytophthora infestans genome", 
              xlab = "Five prime intergenic regions", 
              ylab = "Three prime intergenic regions", 
              cex.main=0.8, cex.lab=0.5), 
       key.title = title(main ="Number of genes", cex.main=0.5, line=1), 
       col=mycol, 
       levels = pretty(zlim, 2*max(GenValMatrix, na.rm=TRUE)), 
       plot.axes={axis(1,at=x, labels=mylabels, las=2, 
              cex.axis=0.5);
       axis(2,at=y, labels=mylabels, cex.axis=0.5)})


#############################################
## 4.2 Overlaying a scatter plot over a genome architecture heatmap
image_name <- paste(as.character(format(Sys.time(), "%Y%m%d%H%M%S")), "_graph", sep="")

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

## Extract genes of interest from gene2gene FIR list
## Note: Reduce all values in first bin to zero, reduce all values above 2nd-to-last bin limit to last-bin minus 100 (protect from jitter)
## Headers: rxlr_five,rxlr_three,family
rxlrData<-as.data.frame(read.csv('Sample_datasets/RXLR_FIRs.csv', header=TRUE))

## Map data points onto base heatmap
## Correct NumBins to equal actual number of bins in BinLimits
NumBins<-as.numeric(c(length(BinLimits)))

gg <-(ggplot(data=FIRdata, aes(x=rxlrData$rxlr_five, y=rxlrData$rxlr_three, geom="blank",fill=factor(family))) + 
       annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
       coord_fixed(ratio=1) + 
       labs(fill = "RXLR Genes") +
       geom_point(shape=21, colour="black", size=5, alpha=0.9, na.rm=FALSE, position=position_jitter(width=0.01, height=0.01)) +
       scale_fill_manual(values = c('rxlr'="red")) +
       scale_y_log10(breaks=BinLimits[2:length(BinLimits)], limits=c(BinLimits[2], BinLimits[NumBins-1])) + 
       scale_x_log10(breaks=BinLimits[2:length(BinLimits)], limits=c(BinLimits[2], BinLimits[NumBins-1])) + 
       ggtitle("RXLR Gene Intergenic Flanking Distance \n vs Global Intergenic Distance ") +
       xlab("Five prime intergenic distance (nt)") + 
       ylab("Three prime intergenic distance (nt)") +
       theme( plot.title = element_text(size=14, face="bold.italic",hjust = 0.5),
              axis.text.y  = element_text(size =8, vjust=0.5),
              axis.text.x  = element_text(size=8, vjust=0.5, angle=90),
              axis.title.x = element_text(face="bold",size=12),
              axis.title.y = element_text(face="bold",size=12))
              )
       
# Run
gg 
# Save
quartz.save("myHeatmap_withPoints.pdf", type="pdf")