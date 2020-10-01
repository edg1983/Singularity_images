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

plotDensity = function(d, field, title=paste0("Distribution of ", field,"\nNA=", sum(is.na(d[[field]])), "/", length(d[[field]]), "=", (sum(is.na(d[[field]]))/length(d[[field]])*100), "%, mean=", mean(d[[field]], na.rm = T) )) {
  #NOPTE: pass empty data frame to vline() to avoid drawing line for each row.
  print(
    ggplot(d, aes_string(x=field)) + 
      geom_density() + 
      geom_vline(aes(xintercept=mean(d[[field]], na.rm = T)), color="blue", linetype="dashed", size=1, data=data.frame()) +
      ggtitle(title) 
  )
}

ALLOWED_GT=c("1/0", "0/1", "1/1", "1|0", "0|1", "1|1")
filterBySamples = function(d, samples) {
  allsamples=substring( colnames(d)[grep("^GT_",colnames(d))] , 4 )
  droplist=allsamples[!allsamples %in% samples]
  droplist = c(paste0("GT_", droplist), paste0("GQ_", droplist) ) 
  d.fil = d %>% select(-droplist)
  d.fil = d.fil %>% filter_at(paste0("GT_", samples), any_vars(. %in% ALLOWED_GT))
  
  return (d.fil)
}

# Constants
IMPACTS = c("UNKNOWN", "MODIFIER", "LOW", "MODERATE", "HIGH")
CONS = c("UNKNOWN", "LOW", "HIGH")
RARENESS = c("VERY_RARE", "RARE", "COMMON", "VERY_COMMON")
GTs=c("0/0", "1/0", "0/1", "1/1")

# load/install packages
list.of.packages <- c("rjson", "data.table", "ggplot2", "reshape2", "dplyr", "stringr", "grid", "gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
for ( p in list.of.packages ) {
  try(library(p, character.only = TRUE), silent=F)
}

# read cmd line args
args <- commandArgs(trailingOnly = TRUE)
# args=c(
#   "/Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/annotest/results/HICF2/HICF2.GRCh38.final.tsv.gz",
#   "/Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/annotest/config.json",
#   "GRCh38",
#   "/Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/annotest/VCF_catalog.txt.small",
#   "/Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/annotest/results/HICF2/qc/"
#  )

print(args)
if (length(args) < 4) {
  stop("usage: Rscript V2_analyze_pedigree.R [final.tsv] [conf.json] [genome] [vcf_catalog] [outdir]")
} 
dataF=args[1]
confF=args[2]
genome=args[3]
catalogF=args[4]
outdir=args[5]

id=tools::file_path_sans_ext(basename(dataF))
conf=fromJSON(file=confF)

# read data
if ( endsWith(dataF, ".gz")) {
  all=fread(paste('gunzip -c', dataF), header=T, sep="\t", na.strings=c("na","NA","."))
} else {
  all=fread(dataF, header=T, sep="\t", na.strings=c("na","NA","."))
}


# read catalog
catalog=data.frame(read.table(catalogF, col.names=c("fid", "PI", "NumSamples", "samples", "vcf", "ped")))

for (row in 1:nrow(catalog)) {
  
  outF = paste0(outdir, "/", catalog[row,]$fid, "_", catalog[row,]$PI, ".qc.pdf")

  # filter by PED
  if(!file.exists(as.character( catalog[row,]$ped))){
    print(paste("Cannot process ", outF, " as PED file is missing:", as.character( catalog[row,]$ped)) )
    next
  }
  ped=read.table(as.character( catalog[row,]$ped), col.names=c("fid", "id", "patid", "matid", "sex", "affected"))
  if ( ! all( paste0("GT_",ped$id) %in% colnames(all) ) ) {
    print(paste("Cannot process ", outF, " as not all samples existing:", paste(c(as.character(ped$id), paste0("GT_",ped$id) %in% colnames(all)), collapse=", ") )  )
    next
  }
  d=filterBySamples(all,ped$id)

  CHROM = c(seq(1,22), "X", "Y")
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
  
  # SNV vs INDEL
  plotTable(d, "TYPE", paste0("SNVs vs INDELs, ratio=", table(d$TYPE)[["SNV"]]/table(d$TYPE)[["INDEL"]]), c("SNV", "INDEL"), datatable=T)
  
  # Known vs unknown
  plotTable(d, "IsKnown", paste0("Variants 'known' in configured databases, frac=", format( table(d$IsKnown)[["1"]]/length(d$IsKnown), digits=2, nsmall=2)), NULL, c("unknown", "known"), log=T, datatable=T)
  
  # Impact
  plotTable(d, "maxImpact", "Maximum impact", IMPACTS, NULL, log=T, datatable=T)
  
  # popAF
  plotTable(d, "maxPopAF", "Maximum population allele frequency", RARENESS, NULL, log=T, datatable=T)
  
  # TsTv
  plotTable(d, "TsTv", paste0("TsTv ratio: ", table(d$TsTv)[["Ts"]]/table(d$TsTv)[["Tv"]]))
  
  
  # ROIs
  rois=names(conf[["roi"]][[genome]])
  roic=(colSums(d %>% select(paste0("ROI_",rois))))
  roiprc=as.numeric(format(roic/nrow(d), digits=2, nsmall=2))
  names(roiprc)=rois
  plotFracVector(roiprc, "Fractions of variants overlapping various\nregions of interest", datatable=T)
  
  # genotypes
  for ( s in samples) {
    dt=table(d[[paste0("GT_",s)]])
    print(plotTable(d, paste0("GT_",s), paste0("Genotype of calls in ", s, "\nHet/hom ratio: ", (dt[["1/0"]]+dt[["0/1"]])/dt[["1/1"]]), GTs))
  }
  
  # various distributions
  for ( f in  conf[["output"]][["fields"]]) {
    if ( is.numeric(d[[f[1]]] ) ) {
      #print(f[[1]])
      plotDensity(d, f[[1]])
    }
  }
  
  # close open connections
  dev.off()

}
