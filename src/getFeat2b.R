## Generate Gene and Intergenic Ranges from GFF/GTF Annotations ## 
##################################################################
##  This is a modified version of the getFeat function by Thomas Girke
##	Computes intergenic, intron and overlapping gene ranges from GFF
##	or GTF files. Results can be combined with existing ranges. Input
##	and output objects are of class 'GRanges'.

getFeat2b <- function(x, gff, format="gff", group_col="group",range_types=c("intergenic", "gene_red", "intron", "gene", "exon")) {
	require(GenomicRanges)
	## Input checks
	if(!missing(gff)) { x <- gff } # For backwards compatibility (if gff instead of x is used as function argument)
	if(class(x)!="GRanges") { stop("Input x needs of class GRanges (e.g. GFF/GTF file imported with import.gff)") }                                                                                                                                                                                                                   
	if(!format %in% c("gff", "GFF", "gtf", "GTF")) { stop("Format argument needs to be GFF or GTF") }
	
	## If seqlengths (chr/space lengths) slot is empty, populate it with max value in each space (chr)
	if(any(is.na(seqlengths(x)))) {
		seqlengths(x) <- as.numeric(tapply(end(x), as.character(seqnames(x)), max))
	}
	
	# If no metadata 'group' in grange object, use ID
	if(!"group" %in% colnames(mcols(x))) {paste("Default 'group' column not found, using: ",group_col)}
	if(!group_col %in% colnames(mcols(x))) {stop("Specified group_col does not exist in granges object.")}
	if(!identical(group_col,"group")) { elementMetadata(x)[,"group"] <- elementMetadata(x)[,group_col]}

	## Exon ranges
	elementMetadata(x) <- elementMetadata(x)[,c("type", "group", "score")]
	elementMetadata(x)[,"score"]<- strand(x)
	elementMetadata(x)[,"type"] <- as.character(as.data.frame(elementMetadata(x)[,"type"])[,1])
	if(format == "gff" | format == "GFF") { 
		elementMetadata(x)[,"group"] <- gsub("^.*?=(.*?)($|,|;|-).*", "\\1", elementMetadata(x)[,"group"])
	}
	if(format == "gtf" | format == "GTF") { 
		geneid <- gsub(";.*", "", elementMetadata(x)[,"group"])
		geneid <- gsub("\"| |gene_id", "", geneid)                                                                                                                                                       
		elementMetadata(x)[,"group"] <- geneid
	}
	strand(x) <- "*" # Erase strand information to make next steps strand insensitive
	x <- x[elementMetadata(x)[,"type"] != "chromosome", ] # Removes chromosome ranges from gff
	exon <- x[elementMetadata(x)[,"type"] == "exon",] # Returns only exon ranges
	
	## Generate gene ranges when input is in GFF format
	if(format == "gff" | format == "GFF") { 
		gene <- x[elementMetadata(x)[,"type"] == "gene",] # Returns only gene ranges
		elementMetadata(gene)["length"]<- as.numeric(2 + end(gene) - start(gene))
	}
	
	## Generate gene ranges when input is in GTF format
	if(format == "gtf" | format == "GTF") { 
		getgeneRanges <- function(x) {                                                                                                                                                                                                            
			tmpdf <- as.data.frame(x)                                                                                                                                                                                                              
			tmpgene <- tmpdf[order(tmpdf$group, tmpdf$start),]                                                                                                                                                                                       
			tmpgene <- tmpgene[!duplicated(tmpgene$group),]
			tmpgeneend <- tmpdf[order(tmpdf$group, -tmpdf$end),]                                                                                                                                                                                     
			tmpgeneend <- tmpgeneend[!duplicated(tmpgeneend$group),]                                                                                                                                                                                 
			tmpgene[,"end"] <- tmpgeneend$end                                                                                                                                                                                                        
			tmpgene[,"width"] <- tmpgene$end - tmpgene$start + 1                                                                                                                                                                                     
			tmpgene[,"type"] <- "gene"
			tmpgene[,"score"]<- tmpgene$score
			tmpgene <- tmpgene[order(tmpgene$seqnames, tmpgene$start),]                                                                                                                                                                              
			gr <- GRanges(seqnames = Rle(tmpgene$seqnames), ranges = IRanges(tmpgene$start, end = tmpgene$end), strand = Rle(tmpgene$strand))                                                                                                        
			elementMetadata(gr) <- tmpgene[,6:8]                                                                                                                                                                                                     
			seqlengths(gr) <- as.numeric(seqlengths(x))
			return(gr)                                                                                                                                                                                                                               
		}
		gene <- getgeneRanges(x)
		elementMetadata(gene)["length"]<- as.numeric(2 + end(gene) - start(gene))
	}

	## Get range from chromosome start to first annotation feature                                         
	start <- tapply(start(gene), as.vector(seqnames(gene)), min)
	start <- paste(names(start), start, sep="_")
	all <- paste(as.vector(seqnames(gene)), start(gene), sep="_")
	firstinter <- gene[which(all %in% start),]
	firstinter <- unique(firstinter) 
	ranges(firstinter) <- IRanges(1, start(firstinter)-1)
	elementMetadata(firstinter)["group"] <- paste("start", as.character(as.data.frame(elementMetadata(firstinter)["group"])[,1]), sep="_")
	## Get range from chromosome end to last annotation feature                 
	chrend <-  seqlengths(x) 
	end <- tapply(end(gene), as.vector(seqnames(gene)), max)
	end <- paste(names(end), end, sep="_")
	all <- paste(as.vector(seqnames(gene)), end(gene), sep="_")
	lastinter <- gene[which(all %in% end),]
	lastinter <- unique(lastinter) 
	ranges(lastinter) <- IRanges(end(lastinter)+1, chrend)
	elementMetadata(lastinter)["group"] <- paste(as.character(as.data.frame(elementMetadata(lastinter)["group"])[,1]), "end", sep="_")

	## Reduced gene ranges
	if(any(range_types %in% c("gene_red", "intergenic"))) {
		## Obtain gene ranges and collapse overlapping genes with reduce
		gene_red <- reduce(gene, min.gapwidth=0)
		ols <- as.matrix(findOverlaps(gene, gene_red))
		red_labels <- tapply(as.character(as.data.frame(elementMetadata(gene)["group"])[,1]), ols[,2], paste, collapse="_")
		elementMetadata(gene_red) <- data.frame(type="gene_red", group=red_labels)
		elementMetadata(gene_red)["type"] <- as.character(as.data.frame(elementMetadata(gene_red)["type"])[,1])
	}

	## Intergenic ranges
	if(any(range_types %in% "intergenic")) {
		tmp_gene_red <- gene_red; start(tmp_gene_red) <- start(tmp_gene_red)+1; end(tmp_gene_red) <- end(tmp_gene_red)-1
		seqlengths(tmp_gene_red) <- rep(NA, length(seqlengths(tmp_gene_red)))
		intergenic <- gaps(tmp_gene_red)
		intergenic <- intergenic[start(intergenic) != 1]
		start(intergenic) <- start(intergenic)+1; end(intergenic) <- end(intergenic)-1
		tmp <- intergenic; start(tmp) <- start(tmp)-1; end(tmp) <- end(tmp)+1
		ols <- as.matrix(findOverlaps(gene_red, tmp))
		inter_labels <- tapply(as.character(as.data.frame(elementMetadata(gene_red)["group"])[ols[,1],1]), ols[,2], paste, collapse="__")
		elementMetadata(intergenic) <- data.frame(type="intergenic", source=NA, phase=NA, group=inter_labels)
		metcol <- intersect(intersect(names(elementMetadata(firstinter)), names(elementMetadata(gene_red))), names(elementMetadata(lastinter)))
		intergenic <- intergenic[order(as.vector(seqnames(intergenic)), start(intergenic))]
		elementMetadata(intergenic)["type"] <- "intergenic"
		elementMetadata(intergenic)["length"]<- as.numeric(1 + end(intergenic) - start(intergenic))
	}
	
	## Intron Ranges
	if(any(range_types %in% "intron")) {
		tmp_exon <- exon
		seqlengths(tmp_exon) <- rep(NA, length(seqlengths(tmp_exon)))
		exonRl <- split(tmp_exon, elementMetadata(tmp_exon)[,"group"])
		intron <- mendoapply(gaps, exonRl) 
		intron <- unlist(intron)
		intron <- intron[start(intron) != 1]
		elementMetadata(intron) <- data.frame(type="intron", group=names(intron))
		names(intron) <- NULL
		elementMetadata(intron)["type"] <- as.character(elementMetadata(intron)[,"type"])
		elementMetadata(intron)["group"] <- as.character(elementMetadata(intron)[,"group"])
	}
	
	## Organinze data components in one GRanges object
	tmpall <- eval(parse(text=paste("c(", paste(range_types, collapse=", "), ")")))
	tmpall <- tmpall[order(as.vector(seqnames(tmpall)), start(tmpall))]
	return(tmpall)
}
