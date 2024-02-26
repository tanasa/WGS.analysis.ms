
###########################################################
###########################################################
###########################################################

library("deconstructSigs")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Hsapiens.UCSC.hg19")

library("reshape2")
library("ggplot2")
library("ggthemes")

library("plyr")
library("dplyr")
library("tidyr")
library("magrittr")

library("cowplot")
library("gridExtra")

library("superheat")
library("corrplot")

library("data.table")
library("splitstackshape")

library("ComplexHeatmap")
library("EnrichedHeatmap")
library("circlize")

###########################################################
###########################################################
###########################################################
############### reading the ANNOVAR FILES : 

read.file.annovar <- function(file)
{
   NAMES <- read.table(file, nrow = 1, stringsAsFactors = FALSE, sep = "\t", quote="", fill=FALSE)
   DATA <- read.table(file, skip = 1, stringsAsFactors = FALSE, sep = "\t", quote="", fill=FALSE)
   DATA <- DATA[, 1:152]
   names(DATA) <- NAMES 
   return(DATA)
}

###########################################################
###########################################################  reading the files whose names end with "hg38_multianno.txt"

files.names <- list.files(pattern="*hg38_multianno.txt") 

list.of.files <- lapply(files.names, read.file.annovar)
# list.of.files <- lapply(files.names, read.delim, header=TRUE, stringsAsFactors=F)
# it reads in all the files and places them in a list, with each element of the list being a data.frame.

names(list.of.files) <- files.names

length.list.of.files <- length(list.of.files)

#####################################################################
### in order to decompose the name of the file into a SIMPLER NAMES :
### having the particle after vcf.... ie. vcf.NAME....
### the names of the SAMPLES will start with SPCG- ...
#####################################################################

for (i in 1:length(list.of.files))
{
     short.name <- strsplit(names(list.of.files)[i],".", fixed=TRUE)[[1]][2]
     ## print(paste("The gene is", short.name))
     names(list.of.files)[[i]] <- short.name
}

######################################################################
######################################################################

index.list.of.files <- 1 

df.sample.mutations <- data.frame( Sample=character(),
                                   Chr=character(),
                                   Start=numeric(),
                                   End=numeric(),
                                   Ref=character(),
                                   Alt=character(), 
                                   stringsAsFactors=FALSE
                                  )

################################################################### 
###################################################################

for (index.list.of.files in 1:length.list.of.files)
{
    nr.sample.mutations <- data.frame( Sample <- names(list.of.files)[[index.list.of.files]],
                                       Chr=list.of.files[[index.list.of.files]]$Chr,
                                       Start=list.of.files[[index.list.of.files]]$Start,
                                       End=list.of.files[[index.list.of.files]]$End,
                                       Ref=list.of.files[[index.list.of.files]]$Ref,
                                       Alt=list.of.files[[index.list.of.files]]$Alt, 
                                       stringsAsFactors=FALSE)

    df.sample.mutations <- rbind(df.sample.mutations, nr.sample.mutations)
} 

###############################################################################################
###############################################################################################
### we did read all the ANNOVAR files and have placed all the files in a dataframe
### that has the following name (below). 
### as we have only 1 FILE, we can also use the NAME for THAT PARTICULAR SAMPLE/FILE .. 


output <- "file.DATAFRAME.from.ANNOVAR.table"

write.table(df.sample.mutations, file=paste(output, short.name, sep="."),
                                 sep="\t", 
                                 quote = FALSE, 
                                 row.names=FALSE)

dim(df.sample.mutations)

################################

df.sample.mutations$Type <- ifelse( ((grepl("-", df.sample.mutations$Ref)=="TRUE") & 
                                     (grepl("-", df.sample.mutations$Alt)=="FALSE")), 
                                      "INS", 
                                      
                                     ifelse( ((grepl("-", df.sample.mutations$Ref)=="FALSE") & 
                                              (grepl("-", df.sample.mutations$Alt)=="TRUE")), 
                                     "DEL", 
                                     "SNV")
)

###############################################################################################
###############################################################################################
################################ selecting only the SNV : 

df.sample.mutations.SNV <- df.sample.mutations[df.sample.mutations$Type=="SNV",]

dim(df.sample.mutations.SNV)

df.sample.mutations.SNV$Sample <- df.sample.mutations.SNV$Sample....names.list.of.files...index.list.of.files..

head(df.sample.mutations.SNV)

###############################################################################################
###############################################################################################
############################### converting to an input for deconstructSigs : 

## df.sample.mutations.SNV.input <- subset(df.sample.mutations.SNV, 
##                                        select=c("Sample", "Chr", "Start", "Ref", "Alt"))

## head(df.sample.mutations.SNV.input)

## or, in order to keep the nomenclature with the rest of the dataframe :

sample.mut.ref <- subset(df.sample.mutations.SNV,
                         select=c("Sample", "Chr", "Start", "Ref", "Alt"))

head(sample.mut.ref)

################################

# library(deconstructSigs)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome.Hsapiens.UCSC.hg19)

sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, 
                                sample.id = "Sample", 
                                chr = "Chr", 
                                pos = "Start", 
                                ref = "Ref", 
                                alt = "Alt", 
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

name <- unique(sample.mut.ref$Sample)

######################################################################################
########################### NORMALIZATION -- GENOME :
########################### using mutation signatures from NATURE 2013 and from COSMIC

wchSig1 = whichSignatures(tumor.ref = sigs.input, 
                         # signatures.ref = signatures.cosmic, 
                         signatures.ref = signatures.nature2013,
                         sample.id = unique(sample.mut.ref$Sample),
                         ## sample.id = "ewing", 
                         contexts.needed = TRUE,
                         tri.counts.method = "genome")

wchSig2 = whichSignatures(tumor.ref = sigs.input, 
                         signatures.ref = signatures.cosmic, 
                         # signatures.ref = signatures.nature2013,
                         sample.id = unique(sample.mut.ref$Sample),
                         ## sample.id = "ewing", 
                         contexts.needed = TRUE,
                         tri.counts.method = "genome")

################################################################################
################################################################################

pdf(paste("SIGNATURES..", name, "..normaliz.genome..ref.nature2013.pdf", sep=""), 
                                        width = 10, height = 10)
chart1 <- plotSignatures(wchSig1)
dev.off()

pdf(paste("SIGNATURES..", name, "..normaliz.genome..ref.cosmic.pdf", sep=""), 
                                        width = 10, height = 10)
chart2 <- plotSignatures(wchSig2)
dev.off()


# pdf(paste("SIGNATURES..", name, "..normaliz.genome..ref.nature2013.pie.charts.pdf", sep=""), 
#                                                                  width = 10, height = 10)
# makePie(wchSig1, sub = 'using signatures nature 2013')
# dev.off()

# pdf(paste("SIGNATURES..", name, "..normaliz.genome..ref.cosmic.pie.charts.pdf", sep=""), 
#                                                               width = 10, height = 10)
# makePie(wchSig2, sub = 'using signatures cosmic')
# dev.off()

##########################################################################################
##########################################################################################
##########################################################################################

wchSig3 = whichSignatures(tumor.ref = sigs.input, 
                         # signatures.ref = signatures.cosmic, 
                         signatures.ref = signatures.nature2013,
                         # sample.id = "ewing",
                         sample.id = unique(sample.mut.ref$Sample), 
                         contexts.needed = TRUE,
                         tri.counts.method = "default")

wchSig4 = whichSignatures(tumor.ref = sigs.input, 
                         signatures.ref = signatures.cosmic, 
                         # signatures.ref = signatures.nature2013,
                         # sample.id = "ewing", 
                         sample.id = unique(sample.mut.ref$Sample),
                         contexts.needed = TRUE,
                         tri.counts.method = "default")

pdf(paste("SIGNATURES..", name, "..param.default..ref.nature2013.pdf", sep=""), 
                                        width = 10, height = 10)
chart3 <- plotSignatures(wchSig3)
dev.off()

pdf(paste("SIGNATURES..", name, "..param.default..ref.cosmic.pdf", sep=""), 
                                        width = 10, height = 10)
chart4 <- plotSignatures(wchSig4)
dev.off()

###############################################################################
################### making the PIE CHARTS :
###############################################################################

# pdf(paste("SIGNATURES..", name, "..param.default..ref.nature2013.pie.charts.pdf", sep=""), 
#                                                                  width = 10, height = 10)
# makePie(wchSig3, sub = 'using signatures nature 2013')
# dev.off()

# pdf(paste("SIGNATURES..", name, "..param.default..ref.cosmic.pie.charts.pdf", sep=""), 
#                                                               width = 10, height = 10)
# makePie(wchSig4, sub = 'using signatures cosmic')
# dev.off()

################### ################### ################### ################### 
################### ################### ################### ################### 
################### ################### ################### ################### 

# library("VariantAnnotation")
# vcf<-readVcf("","hg38")
# vcf.to.sigs.input

COLORS_SIGNATURES = c(
"Signature.1"  = "coral1", 	
"Signature.2"  = "chocolate1",	
"Signature.3"  = "chartreuse1",	
"Signature.4"  = "brown1",	
"Signature.5"  = "blue1",
"Signature.6"  = "aquamarine",	
"Signature.7"  = "dodgerblue1",	
"Signature.8"  = "deeppink3",	
"Signature.9"  = "darkseagreen1",	
"Signature.10" = "darkorchid2",	
"Signature.11" = "darkorange1",	
"Signature.12" = "darkcyan",	
"Signature.13" = "cyan1",	
"Signature.14" = "firebrick2",	
"Signature.15" = "indianred2",	
"Signature.16" = "hotpink1",	
"Signature.17" = "maroon1",	
"Signature.18" = "magenta1",	
"Signature.19" = "linen",	
"Signature.20" = "lightgrey",	
"Signature.21" = "lightcyan",	
"Signature.22" = "palevioletred1",	
"Signature.23" = "palegreen1",	
"Signature.24" = "orchid",	
"Signature.25" = "orangered2",	
"Signature.26" = "orange2",	
"Signature.27" = "navyblue",	
"Signature.28" = "sienna3",	
"Signature.29" = "salmon3",	
"Signature.30" = "red2"
)

################### ################### ################### ################### 
################### ################### ################### ################### 
################### ################### ################### ################### 

# in order to obtain the weights, and to write the weights, please let's do :

### here it is the NOT NORMALIZED SIGNATURE :

write.table(wchSig4$weights, 
file=paste("SIGNATURES..", name, "..param.default..ref.cosmic.TABLE.txt", sep=""), 
                           sep="\t", 
                           quote = FALSE, 
                           row.names = FALSE,
                           col.names = TRUE )

#### here collecting the WEIGHTS of the SIGNATURES

w4 <- wchSig4$weights

head(w4)

w4m <- melt(w4)

################################################################
#### and making a few TYPES of PLOTS : TYPE 1
################################################################

ggplot(w4m, aes(x=variable, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("\n COSMIC signatures") +
  ylab("signature weight \n") +
  scale_colour_manual(values = COLORS_SIGNATURES) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(size=12)) +
  theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank()
         ) + 
  ggtitle(paste(name, " : deconstructSigs param default", paste="")) 

ggsave(file=paste("SIGNATURES..", name, "..param.default..ref.cosmic.PICTURE.multiple.bars.pdf", sep=""), 
       height=14, width=30, units="cm")

################################################################
#### and making a few TYPES of PLOTS : TYPE 2
################################################################

ggplot(w4m, aes(1, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("\n COSMIC signatures") +
  ylab("signature weight \n") +
  scale_colour_manual(values = COLORS_SIGNATURES) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(size=12)) +
  theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank()
         ) + 
  ggtitle(paste(name, " : deconstructSigs param default", paste="")) + coord_flip()

ggsave(file=paste("SIGNATURES..", name, "..param.default..ref.cosmic.PICTURE.bars.horizontal.pdf", sep=""), 
       height=14, width=30, units="cm")

################################################################
#### and making a few TYPES of PLOTS : TYPE 3
################################################################

ggplot(w4m, aes(1, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("\n COSMIC signatures") +
  ylab("signature weight \n") +
  scale_colour_manual(values = COLORS_SIGNATURES) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(size=12)) +
  theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank() 
  ) + 
  ggtitle(paste(name, " : deconstructSigs param default", paste="")) 

ggsave(file=paste("SIGNATURES..", name, "..param.default..ref.cosmic.PICTURE.bars.vertical.pdf", sep=""), 
       height=30, width=10, units="cm")


#### the NORMALIZED SIGNATURE : wchSig2 
################### ################### ################### ################### 
################### ################### ################### ################### 
################### ################### ################### ################### 

# in order to obtain the weights, and to write the weights, please let's do :
# the NORMALIZED SIGNATURE :

write.table(wchSig2$weights, 
file=paste("SIGNATURES..", name, "..normaliz.genome..ref.cosmic.TABLE.txt", sep=""), 
                           sep="\t", 
                           quote = FALSE, 
                           row.names = FALSE,
                           col.names = TRUE )

### collecting the WEIGHTS of the SIGNATURES

w2 <- wchSig2$weights

head(w2)

w2m <- melt(w2)

################################################################
#### and making a few TYPES of PLOTS : TYPE 1
################################################################

ggplot(w2m, aes(x=variable, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("\n COSMIC signatures") +
  ylab("signature weight \n") +
  scale_colour_manual(values = COLORS_SIGNATURES) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(size=12)) +
  theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank()
         ) + 
  ggtitle(paste(name, " : deconstructSigs normaliz genome", paste="")) 

ggsave(file=paste("SIGNATURES..", name, "..normaliz.genome..ref.cosmic.PICTURE.multiple.bars.pdf", sep=""), 
       height=14, width=30, units="cm")

################################################################
#### and making a few TYPES of PLOTS : TYPE 2
################################################################

ggplot(w2m, aes(1, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("\n COSMIC signatures") +
  ylab("signature weight \n") +
  scale_colour_manual(values = COLORS_SIGNATURES) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(size=12)) +
  theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank()
         ) + 
  ggtitle(paste(name, " : deconstructSigs normaliz genome", paste="")) + coord_flip()

ggsave(file=paste("SIGNATURES..", name, "..normaliz.genome..ref.cosmic.PICTURE.bars.horizontal.pdf", sep=""), 
       height=14, width=30, units="cm")

################################################################
#### and making a few TYPES of PLOTS : TYPE 3
################################################################

ggplot(w2m, aes(1, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("\n COSMIC signatures") +
  ylab("signature weight \n") +
  scale_colour_manual(values = COLORS_SIGNATURES) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(size=12)) +
  theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank() 
  ) + 
  ggtitle(paste(name, " : deconstructSigs normaliz genome", paste="")) 

ggsave(file=paste("SIGNATURES..", name, "..normaliz.genome..ref.cosmic.PICTURE.bars.vertical.pdf", sep=""), 
       height=30, width=10, units="cm")

################################################################
################################################################
################################################################