#!/usr/bin/Rscript

library(data.table)

args=commandArgs(TRUE)

infile=args[1]
minfreq=args[2]

outfile=gsub(".tsv", "_recount.tsv", infile)
outfile_minfreq=gsub(".tsv", paste0("_recount_minfreq_",minfreq,".tsv"), infile)

barcode_stats=read.csv(infile, sep = "\t", header = T)

if (!file.exists("rrnDB-5.8_pantaxa_stats_NCBI.tsv")) {
  system('wget -c https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.8_pantaxa_stats_NCBI.tsv.zip')
  system('unzip rrnDB-5.8_pantaxa_stats_NCBI.tsv.zip')
}

count_stats=read.csv("rrnDB-5.8_pantaxa_stats_NCBI.tsv", sep = "\t", header = T)

count_stats$sum16slist=gsub("[][]", "", count_stats$sum16slist)

num_rows=nrow(barcode_stats)
for (row in 1:num_rows) {
  species=barcode_stats$species[row]
  df=subset(count_stats, count_stats$name == species)
  #print(length(df$median))
  
  if (length(df$median)>0) {
    barcode_stats[row,15] = df$median
  } else {
    genus=barcode_stats$genus[row]
    #print(genus)
    df = subset(count_stats, count_stats$name == genus)
    df = df[,11]
    df1 = data.table::tstrsplit(df, ", ")
    df2 = as.matrix(df1)
    class(df2) = "numeric"
    median = median(df2)
    if (!is.na(median)) {
      barcode_stats[row,15] = median
    } else {
      family=barcode_stats$family[row]
      df = subset(count_stats, count_stats$name == family)
      df = df[,11]
      df1 = data.table::tstrsplit(df, ", ")
      df2 = as.matrix(df1)
      class(df2) = "numeric"
      median = median(df2)
      if (!is.na(median)) {
        barcode_stats[row,15] = median
      } else {
        order=barcode_stats$order[row]
        df = subset(count_stats, count_stats$name == order)
        df = df[,11]
        df1 = data.table::tstrsplit(df, ", ")
        df2 = as.matrix(df1)
        class(df2) = "numeric"
        median = median(df2)
        if (!is.na(median)) {
          barcode_stats[row,15] = median
        } else {
          class=barcode_stats$class[row]
          df = subset(count_stats, count_stats$name == class)
          df = df[,11]
          df1 = data.table::tstrsplit(df, ", ")
          df2 = as.matrix(df1)
          class(df2) = "numeric"
          median = median(df2)
          if (!is.na(median)) {
            barcode_stats[row,15] = median
          } else {
            phylum=barcode_stats$phylum[row]
            df = subset(count_stats, count_stats$name == phylum)
            df = df[,11]
            df1 = data.table::tstrsplit(df, ", ")
            df2 = as.matrix(df1)
            class(df2) = "numeric"
            median = median(df2)
            if (!is.na(median)) {
              barcode_stats[row,15] = median
            } else {
              if (barcode_stats$tax_id[row] == "unassigned") {
                barcode_stats[row,15] = 1
              }
            }
          }
        }
      } 
    }
  } 
}


colnames(barcode_stats)[colnames(barcode_stats) == "V15"] = "gene.count"
barcode_stats$new.count = barcode_stats$estimated.counts/barcode_stats$gene.count
barcode_stats$new.abundance = barcode_stats$new.count/sum(barcode_stats$new.count, na.rm = T)
write.table(barcode_stats, outfile, sep = "\t", quote=F, row.names = F)

barcode_stats_minfreq=subset(barcode_stats, barcode_stats$new.abundance > minfreq)
barcode_stats_minfreq$new.abundance = barcode_stats_minfreq$new.count/sum(barcode_stats_minfreq$new.count, na.rm = T)
write.table(barcode_stats_minfreq, outfile_minfreq, sep = "\t", quote=F, row.names = F)
