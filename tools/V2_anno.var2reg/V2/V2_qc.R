# This R script will create some QC output
# 
# Author: niko@well.ox.ac.uk
###############################################################################


ttheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.7)),
  rowhead = list(fg_params=list(cex = 0.7)))

plotTable = function(d, field, title, factorOrder=NULL, labels=NULL, log=F, datatable=F) {
  dat=NULL
  if ( is.null( factorOrder) ) {
    dat=table(factor(d[[field]]), exclude = NULL)
  } else {
    dat=table(factor(d[[field]], levels=factorOrder), exclude = NULL)
  }
  xax=NULL
  if ( ! is.null(labels)) {
    xax = scale_x_discrete(labels=labels) 
  } else {
    xax=scale_x_discrete()
  }
  yax=NULL
  if ( log == T ) {
    yax = scale_y_log10(name="log(count)") 
  } 
  
  plt=(ggplot(melt(dat), aes(x=factor(Var1), y=value)) +
    geom_bar(stat="Identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title) + xlab("") + ylab("count") + xax + yax)
  #tbl=tableGrob(grid.table(t(dat), theme=ttheme))
  tbl=tableGrob(t(dat), theme=ttheme)
  if (datatable) {
    grid.arrange(plt, tbl,
                 nrow=2,
                 as.table=TRUE,
                 heights=c(5,1))
  } else {
    print(plt)
  }
}

plotFracVector = function(vec, label, title, datatable=F) {
  plt= ggplot(data.frame(vec), aes(x=seq_along(vec),y=vec, labels=names(vec))) + 
       geom_bar(stat="Identity", position = "dodge") +
       ylim(0,1) +
       theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_x_discrete(name = "", limits=seq_along(vec), labels=names(vec)) +
       ggtitle(label) + xlab("") + ylab("fraction")
  tbl=tableGrob(t(vec),theme=ttheme)
  if (datatable) {
    grid.arrange(plt, tbl,
                 nrow=2,
                 as.table=TRUE,
                 heights=c(5,1))
  } else {
    print(plt)
  }
}

plotDensity = function(d, field, title=paste0("Distribution of ", field,"\nNA=", sum(is.na(d[[field]])), "/", length(d[[field]]), "=", format(sum(is.na(d[[field]]))*100/length(d[[field]]), digits=2, nsmall=2), "%, mean=", format(mean(d[[field]], na.rm = T), digits=2, nsmall=2) )) {
  print(
    ggplot(melt(d[[field]]), aes(x=value)) + 
      geom_density() + 
      geom_vline(aes(xintercept=mean(d[[field]], na.rm = T)), color="blue", linetype="dashed", size=1) +
      ggtitle(title) 
  )
}

# Constants
IMPACTS = c("UNKNOWN", "MODIFIER", "LOW", "MODERATE", "HIGH")
CONS = c("UNKNOWN", "LOW", "HIGH")
RARENESS = c("VERY_RARE", "RARE", "COMMON", "VERY_COMMON")
CHROM = c(seq(1,22), "X", "Y")
GTs=c("0/0", "1/0", "0/1", "1/1")
CHR_SIZE_hg38 = data.frame(
  CHROM=c(paste0("chr", seq(1,22)), "chrX", "chrY", "chrM"),
  size=c(
  248956422,
  242193529,
  198295559,
  190214555,
  181538259,
  170805979,
  159345973,
  145138636,
  138394717,
  133797422,
  135086622,
  133275309,
  114364328,
  107043718,
  101991189,
  90338345,
  83257441,
  80373285,
  58617616,
  64444167,
  46709983,
  50818468,
  156040895,
  57227415,
  16569)
)
CHR_SIZE_hg38$CHROM=factor(CHR_SIZE_hg38$CHROM, levels=CHR_SIZE_hg38$CHROM)

# load/install packages
list.of.packages <- c("rjson", "data.table", "ggplot2", "reshape2", "dplyr", "stringr", "grid", "gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, FUN = function(X) { do.call("require", list(X)) })

# read cmd line args
args <- commandArgs(trailingOnly = TRUE)
# args=c(
#   "/Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/cohort_test/01_anno/cohort_platypus_joincalls/cohort_platypus_joincalls.GRCh38.final.tsv.gz",
#   "/Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/cohort_test/01_V2.config.json",
#   "GRCh38",
#   "/Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/cohort_test/01_anno/cohort_platypus_joincalls/qc/new.qc.pdf"
#  )

print(args)
if (length(args) < 4) {
  stop("usage: Rscript V2_qc.R [final.tsv] [conf.json] [genome] [out.pdf]")
} 
dataF=args[1]
confF=args[2]
genome=args[3]
outF=args[4]

id=tools::file_path_sans_ext(basename(dataF))
conf=fromJSON(file=confF)

# read data
if ( endsWith(dataF, ".gz")) {
  d=fread(paste('gunzip -c', dataF), header=T, sep="\t", na.strings=c("na","NA","."))
} else {
  d=fread(dataF, header=T, sep="\t", na.strings=c("na","NA","."))
}

if ( startsWith(as.character( d[1,1] ),"chr") ) {
  CHROM=paste0("chr", CHROM)
}

# get sample names
samples=substring( colnames(d)[grep("^GT_",colnames(d))] , 4 )
print(paste("Writing QC data for", id, " and samples ", paste(samples, collapse = ",")))

# write all output to log file
pdf(outF)
grid.newpage()
grid.table(strwrap(paste0(c(paste(id, "samples:"), samples), collapse = ", "), width=50))

# check how many rows contain a comma in the respective score:
# length(which(grepl(",",d$FIRE_FIRE)))

# Vars per chrom
plotTable(d, "CHROM", "Number of (filtered) variants per chromosome", CHROM)
tab = d %>% group_by( CHROM) %>% summarise( n = n())
tab = merge(CHR_SIZE_hg38, tab, by="CHROM") 
tab$CHROM=factor(tab$CHROM, levels=CHR_SIZE_hg38$CHROM)
ggplot( tab %>% mutate( n_norm = n/size), aes(x=CHROM, y=n_norm) ) + 
          geom_bar(stat='identity') + 
          ylab("#variants normalized by chrom size") + xlab("") +
          ggtitle( "Number of (filtered) variants per chromosome, normalized")

# SNV vs INDEL
plotTable(d, "TYPE", paste0("SNVs vs INDELs, ratio=", table(d$TYPE)[["SNV"]]/table(d$TYPE)[["INDEL"]]), c("SNV", "INDEL"), datatable=T)

# Known vs unknown
plotTable(d, "IsKnown", paste0("Variants 'known' in configured databases, frac=", format( table(factor(d$IsKnown, c("0","1")))[["1"]]/length(d$IsKnown), digits=2, nsmall=2)), NULL, c("unknown", "known"), log=T, datatable=T)

# Impact
plotTable(d, "maxImpact", "Maximum impact", IMPACTS, NULL, log=T, datatable=T)

# popAF
afs=names(conf[["af"]][[genome]])
if ( ! is.null(afs) ) {
  for ( id in afs ) {
    plotTable(d, paste0("maxPopAF_cat_",id), paste0("Maximum population allele frequency (",id,")"), RARENESS, NULL, log=T, datatable=T)
  }
}
# TsTv
plotTable(d, "TsTv", paste0("TsTv ratio: ", table(d$TsTv)[["Ts"]]/table(d$TsTv)[["Tv"]]))

# annotation values with NA 
n_all = nrow(d)
n_indel = nrow(d %>% filter( TYPE=='INDEL'))
na_table = d %>% select(-starts_with("GQ_"),-starts_with("GT_")) %>% summarise_all(funs(sum(! is.na(.))))
print(
  ggplot( reshape2::melt( na_table), aes(x=variable, y=value)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  geom_hline(yintercept = n_all, col="red", linetype="dotted") +
  geom_hline(yintercept = n_all - n_indel, col="green", linetype="dotted") +
  ggtitle("# non-NA value per annotation field", "red: n_all, green: n_all - n_indel ")
)

# ROIs
rois=names(conf[["roi"]][[genome]])
if ( ! is.null(rois) ) {
  roic=(colSums(d %>% select(paste0("ROI_",rois))))
  roiprc=as.numeric(format(roic/nrow(d), digits=2, nsmall=2))
  names(roiprc)=rois
  plotFracVector(roiprc, "Fractions of variants overlapping various\nregions of interest", datatable=T)
  }

# various distributions
for ( f in  conf[["output"]][["fields"]]) {
  if ( is.numeric(d[[f[1]]] ) ) {
    #print(f[[1]])
    plotDensity(d, f[[1]])
  }
}

# individual genotypes
for ( s in samples) {
  dt=table(d[[paste0("GT_",s)]])
  plotTable(d, paste0("GT_",s), paste0("Genotype of calls in ", s, "\nHet/hom ratio: ", (dt[["1/0"]]+dt[["0/1"]])/dt[["1/1"]]), GTs)
}

# close open connections
dev.off()
