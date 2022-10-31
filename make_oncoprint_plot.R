#


################################################################
### EXAMPLE CALL

#module load R/3.4.0
#module load yapsa-devel/80f748e
#R -f make_oncoprint_plot --no-save --no-restore --args --input_table <onco_print table from make_oncprint_table.pl> -sampleinfo_table <sample_info table from make_oncprint_table.pl> --min_recurrence=4 --cnas_num=6 --annotation_table=<(optional) input file containing sammples in first column, then any number of columns for annotation> --group_over=<see description below> etc.. 

################################################################
### Install ComplexHeatmap
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

################################################################
### LOAD LIBRARIES

print(paste("Loading libraries..."))
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(GetoptLong)
library(YAPSA)

################################################################
### CONFIGURE DEFUALT PARAMS

version=paste("v0.15")

set.seed(123)

title=paste("Recurrently mutated genes")
annotation_table=NA
gene_info="NA"
group_over=paste("NA")
remove_features=paste("intronic_snv,SV_TAD,SV_near,UTR_3_snv,UTR_5_snv,UTR_3_indel,UTR_5_indel,synonymous_SNV,intronic_indel")
features_to_keep="NA"
min_recurrence=1;
intogen_file=paste("NA")
intogen_pvalue_column=paste("MUTSIG_PVALUE")
min_significance=1;
cnas_num = 6
selected_gene_list=NA

################################################################
### PARSE COMMAND LINE ARGUMENTS

GetoptLong(
    "title=s",                 "title of the oncoprint",
    "input_table=s",           "input oncorpint mutation matrix from oncorpint_table script (*.oncoprint.tsv)",
    "sampleinfo_table=s",      "sample info file from oncorpint_table script (*.sample_info.tsv)",
    "annotation_table=s",      "custom annotation file, with 'Sample' colummn as indentifier",
    "gene_info=s",             "", # ?
    "group_over=s",            "feature for sample ordering taken from annotation, sampleinfo or oncoprint fields (default none, e.g 'TP53' or 'CNA sex')",
    "remove_features=s",       "comma separated feature list to remove (recommend removing UTRs, intronic, SV TAD and SV near)",
    "features_to_keep=s",      "comma separated feature list to keep (this over rides remove_features)", 
    "min_recurrence=i",        "minimum recurrence threshold (default: 1)",
    "intogen_file=s",          "path to intogen results",
    "intogen_pvalue_column=s", "intogen significance column (default: MUTSIG_PVALUE)",
    "min_significance=i",      "minimum intogen p value (default: 1)",
    "cnas_num=i",              "top CNVs to be used in heatmap annotation (default: cnas_num)",
    "selected_gene_list=s",    "file containing selected list of genes, one gene per row",
    "verbose!",                "print messages"
)


print(paste("Input mutation matrix:", input_table))
print(paste("Sample info file:", sampleinfo_table))
print(paste("Sample annotation file:", annotation_table))
print(paste("Grouping fields:", group_over))
print(paste("Feature to remove:",remove_features))
print(paste("Minimum recurrence:", min_recurrence))
print(paste("Selected_gene_list:", selected_gene_list))
print(title)

################################################################
### LOAD SORTING FUNCTION

# library(matuationalDensity)

print(paste("Loading oncoprint ordering function (Daniel Huebschmann)..."))
### Function to order the mutation matrix from Daniel Huebschmann from MutationalDensity (d.huebschmann@dkfz.de)
oncoprintOrder <- function(in_matrix,
                          in_subgroup_rank_vector = NULL){
  # create sister matrix for reordering
  order_mat <- as.matrix(in_matrix)
  order_mat[order_mat == ""] <- 0
  order_mat[order_mat != 0] <- 1
  order_mat <- apply(order_mat, 2 , as.numeric)
  # reorder features in oncoprint like fashion
  feature_ind <- order(rowSums(order_mat),decreasing = TRUE)
  ordered_feature_mat <- order_mat[feature_ind,]
  # reorder samples in oncoprint like fashion
  vector_list <- split(ordered_feature_mat,seq(nrow(ordered_feature_mat)))
  # now include subgroup information
  if(!is.null(in_subgroup_rank_vector)){
    if(!is.null(names(in_subgroup_rank_vector))){
      matching_ind <- match(colnames(in_matrix),names(in_subgroup_rank_vector))
      if(all(is.finite(matching_ind))) 
        in_subgroup_rank_vector <- in_subgroup_rank_vector[matching_ind]
    }
    vector_list <- add_as_fist_to_list(vector_list,
                                       in_subgroup_rank_vector)
  }
  sample_ind <- do.call(order,c(vector_list,decreasing = TRUE))
  return(list(matrix = in_matrix[feature_ind,sample_ind],
              feature_ind = feature_ind,
              sample_ind = sample_ind))
}

################################################################
### COLOR MATRIX

#source("/home/kleinhei/Project/mmml/mmml.git/oncoprints/colorAnnotationSpecifications.R")

col <- c(
        # normal mutations are as defined by Kortine
	"frameshift_deletion" = "#5ac4fc",
	"frameshift_insertion" = "#5ac4fc", 
	"nonframeshift_deletion" = "#a4defc", 
	"nonframeshift_insertion" = "#a4defc", 
	"nonsynonymous_SNV" = "#FB9A99", 
	"splicing_snv" = "#E31A1C", 
	"splicing_indel" = "#E31A1C",
	"stopgain_snv" = "#FDBF6F",
	"stopgain_indel" = "#FDBF6F", 
	"stoploss_snv" = "#FF7F00", 
	"stoploss_indel" = "#FF7F00",
	"synonymous_SNV" = "#CAB2D6", 

        # germline in black
	"frameshift_deletion_germline" = "#000000",
	"frameshift_insertion_germline" = "#000000", 
	"nonframeshift_deletion_germline" = "#000000", 
	"nonframeshift_insertion_germline" = "#000000", 
	"nonsynonymous_SNV_germline" = "#000000", 
	"splicing_snv_germline" = "#000000", 
	"splicing_indel_germline" = "#000000",
	"stopgain_snv_germline" = "#000000",
	"stopgain_indel_germline" = "#000000", 
	"stoploss_snv_germline" = "#000000", 
	"stoploss_indel_germline" = "#000000",
	"synonymous_SNV_germline" = "#000000", 

        # ncRNA exonic greens
	"ncRNA_exonic" = "#33FF33",
	"ncRNA_exonic_snv" = "#33FF33",
	"ncRNA_exonic_indel" = "#009900",

        # sv now follow purples with  more spread
	"SV_direct" = "#800080",
	"SV_near" = "#BA55D3",
	"SV_TAD" = "#DA70D6",
        #"inversion" = "black",

        # CNVs are now recoloured to match previous CNV plots (output from GISTIC2) 
        "high_amplification" = "#663300",
	"amplification" = "#CC6600",
        #"duplication" = "#CC6600",
        "homo_del" = "#000066",
	#"loss" = "#0066CC",
        "deletion"= "#0066CC",
	#"cn_LOH" = "#CC0000",

	# LOH events in green, as this colour was not used yet
	"highAmp_LOH" = "#333300",
	"amp_LOH" = "#FFFFCC",
	"del_LOH" = "#99CCCC",
	"cn_LOH" = "#006600",

	# chromosome level events follow the cnv colouring
        "chrAmplification" = "#CC6600",
	"chrAmpLOH" = "#FFFFCC",
        "chrDeletion" = "#0066CC",
        "chrDeletionLOH" = "#99CCCC",
        "chrCnLOH" = "#006600",
        "chrHomoDel" = "#000066",

        # intronic and kataegis in greys
	"intronic_snv"="#CCCCCC",
	"intronic_indel"="#AAAAAA", 
        "kataegis" = "#888888",

        # UTR browns
	"UTR_3_snv"="#FAE5D3",  
	"UTR_5_snv"="#F6DDCC",  
	"UTR_3_indel"="#FAE5D3",
	"UTR_5_indel"="#F6DDCC",

        # upstream in yellow
        "upstream_snv" = "#FFFF66",
        "upstream_indel" = "#FFFF66",

       	# rna fusion in pink
        "rna_fusion" = "#ff0089"
        )

################################################################
### ALTER FUNCTION

alter_fun <- function(x, y, w, h, v) {
               n = sum(v)
               h = h*0.9
               w = w*0.9
               grid.rect(x, y, w, h, gp = gpar(fill = "#CCCCCC", col = NA))
               if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h, gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
             }

################################################################
### READ ONCOPRINT TABLE

## Read table
print(paste("Reading oncorpint file..."))
mat = read.table(input_table, header=T, sep="\t", check.names=FALSE, row.names=1)
mat[is.na(mat)] = ""

## Remove counts
print(paste("Removing counts..."))
mat_filtered <- data.frame(lapply(mat, function(x) { gsub(":\\d*;", ";", x) }))
rownames(mat_filtered) <- rownames(mat)
colnames(mat_filtered) <- colnames(mat)
mat <- as.data.frame(mat_filtered)

################################################################
### READ SAMPLEINFO TABLE

# read sample info table
print(paste("Reading sampleinfo file..."))
anno<- read.table(sampleinfo_table, header=T, sep="\t", row.names=1, check.names=FALSE, na.strings=c("","NA"))
annot=transpose(anno)
colnames(annot) <- rownames(anno)
rownames(annot) <- colnames(anno)

for (index in grep("CNA ch|CNA sex",names(annot), invert=T)){
  annot[,index] <- as.numeric(annot[,index])
}

################################################################
### READ ANNOTATION TABLE

## Read table
if(!is.na(annotation_table)) {
  custom_table = read.table(annotation_table, header=T, sep="\t", row.names=1)
  print(paste("Reading annotation file..."))
  if(all(rownames(custom_table) %in% colnames(mat))) {
    print(paste("not sure what to do"))
    all_annot <- merge(annot, custom_table, by="row.names")
    rownames(all_annot) <- all_annot[[1]]
    all_annot <- all_annot[-1]
    annot<-all_annot[rownames(annot),]
  } else {
    stop("Custom annotation row names don't match the oncoprint mat column names")
  }
}


################################################################
### BUILD TOP ANNOTATION

### SNV_ anno

snv_annot <- as.numeric(annot[,"SNV - total"])
snv_col <- c("lightgrey")

### INDEL_ anno

indel_annot <- as.numeric(annot[, "INDEL - total"])
indel_col <- c("lightgrey")

### SV annotation

sv_annot <- as.numeric(annot[,"SV - total"])
sv_col <- c("lightgrey")
 
### SEX

sex_mat<-annot$"CNA sex"
sex_mat <- gsub("female", "F", sex_mat)
sex_mat <- gsub("male", "M", sex_mat)

### ploidy

ploidy_mat<-annot$"CNA ploidy"
max_ploidy <- max(as.numeric(ploidy_mat), na.rm = TRUE)

### TOP CNAs

cna_mat  <- annot[,names(annot[,grep("CNA chr", names(annot))])]
na_count <- sapply(cna_mat, function(y) sum(length(which(is.na(y)))))
top_cnas <- names(sort(na_count))
cna_mat  <- as.matrix(cna_mat[,top_cnas[1:cnas_num]])

## CNA TCC

cna_tcc <- annot["CNA TCC"]

### Make heatmap annotations

column_ha = HeatmapAnnotation(
                                     "Total SVs" = anno_barplot(sv_annot, axis = TRUE, gp = gpar(fill = sv_col)),
                                     "Total INDELs" = anno_barplot(indel_annot, axis = TRUE, gp = gpar(fill = indel_col)),
                                     "Total SNVs" = anno_barplot(snv_annot, axis = TRUE, gp = gpar(fill = snv_col)),
                                     #sex_text = anno_text(sex_mat, rot = 0, just = "centre", offset = unit(2, "mm")),
                                     sex = sex_mat,
                                     ploidy = ploidy_mat,
                                     tcc = as.matrix(cna_tcc)[,1],
                                     cna = cna_mat,
                                     col = list (sex = c("M"="lightblue", "F"="pink"),
                                                 ploidy = colorRamp2(c(0, max_ploidy), c("white", "blue")), 
                                                 tcc = colorRamp2(c(0, 1), c("white", "red")),
                                                 cna = c("amp;"="grey" , "LOH;"="grey" , "del;"="grey" , "homoDel;"="grey" ,

                                                         "amp;LOH;"="grey" , "amp;del;"="grey" , "amp;homoDel;"="grey" ,
                                                         "amp;del;LOH;"="grey" ,  "amp;del;homoDel;"="grey" ,
                                                         "amp;del;homoDel;LOH;"="grey" ,

                                                         "del;LOH;"="grey" , "LOH;homoDel;"="grey" ,
                                                         "del;homoDel;LOH;"="grey" ,

                                                         "del;homoDel;"="grey"
                                                        )
                                                 ),
#                                     height=unit(13, "cm"),
                                     annotation_name_gp = gpar(fontsize=10),
                                     show_annotation_name=TRUE,
                                     annotation_height = unit(c(15, 15, 15, 5, 5, 5, 10*(cnas_num/2)), "mm"),
                                     show_legend = c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
                            )

### Custom annotations

if(!is.na(annotation_table)) {
  custom_ha <- HeatmapAnnotation(df = annot[colnames(custom_table)], annotation_name_gp = gpar(fontsize=10), show_annotation_name=T)
  column_ha = c(column_ha, custom_ha)
}


################################################################
### FILTER ONCOPRINT

## Filter features in oncoprint
print(paste("Filtering features..."))
if (features_to_keep=="NA"){
  feature_to_remove <- gsub(",",";|\\\\b",remove_features)
  feature_to_remove <- paste0("\\b",feature_to_remove,";")
  print(paste("No features_to_keep to defined. Instead removing undesired features...'", feature_to_remove,"'", sep =""))

  mat_filtered2 <- data.frame(lapply(mat, function(x) { gsub(feature_to_remove, "", x) }))

  rownames(mat_filtered2) <- rownames(mat)
  colnames(mat_filtered2) <- colnames(mat)

  mat <- mat_filtered2

} else {
  print(paste("Features_to_keep to defined. Will keep '",features_to_keep,"'", sep =""))

  # make empty dataframe of size mat, with empty elements
  mat_selection <- matrix("", nrow=nrow(mat), ncol=ncol(mat))
  rownames(mat_selection) <- rownames(mat)
  colnames(mat_selection) <- colnames(mat)

  # foreach feature to keep, use grep to populate mat_selection
  feature_list <- strsplit(features_to_keep, ",")[[1]]
  feature=NULL

  for (feature in feature_list){
    feature_bool <- grepl(paste0("\\b",feature,";"), as.matrix(mat))
    mat_selection[feature_bool] <- paste0(mat_selection[feature_bool],feature,";")
  }

  mat_selection<- as.matrix(mat_selection)
  rownames(mat_selection) <- rownames(mat)
  colnames(mat_selection) <- colnames(mat)
  mat <- mat_selection
}

###############################################################
### select genes
cat("Selecting gene from list")
if(!is.na(selected_gene_list)){
  select_genes = scan(selected_gene_list, character())
  mat_select <- mat
  #select_list <- strsplit(select_genes, ",")[[1]]
  mat_select <- mat_select[rownames(mat_select) %in% select_genes,]
  mat <- mat_select
}

################################################################
### SUBGROUPING

grouping_vector<-rep(1, ncol(mat))
names(grouping_vector)<- colnames(mat)

## Determine subgroup ordering
if(!group_over=="NA"){

    print(paste("Determining ordering vector from '",group_over,"'", sep=""))
    feature=group_over

    if (feature %in% names(annot)) {
      print(paste0("Found feature ",feature," in sampleannotation"))
      group_temp <- as.matrix(annot[feature])
      group_temp [is.na(group_temp)] = ""
      names(group_temp) <- colnames(mat)
      group_temp_uniq <- as.matrix(unique(sort(as.matrix(group_temp), decreasing=FALSE)))
      iterator=0
      grouping_vector<-setNames(as.numeric(factor(group_temp)), rownames(group_temp))[names(grouping_vector)]
      # for(level in group_temp_uniq){
      #   print(paste0("Subgrouping sampleinfo feature ",group_over," - ", level))
      #   grouping_vector[as.matrix(group_temp)==level] <- grouping_vector[as.matrix(group_temp)==level] + iterator
      #   iterator = iterator + 1
      # }
    } else if (feature %in% rownames(mat)) {
      print(paste0("Found feature ",feature," in oncoprint matrix"))
      group_temp <- as.matrix(mat[feature,])
      group_temp[!(group_temp=="")] = "mutated"
      group_temp[(group_temp=="")] = "wildtype"
      names(group_temp) <- colnames(mat)
      group_temp_uniq <- rev(unique(sort(as.matrix(group_temp), decreasing=FALSE)))
      iterator=0
      for(level in group_temp_uniq){
        print(paste0("Subgrouping gene feature ",group_over," - ", level))
        grouping_vector[as.matrix(group_temp)==level] <- grouping_vector[as.matrix(group_temp)==level] + iterator
        iterator = iterator + 1
      }
    } else {
      print(paste0("ERROR: couldnt find '",feature,"'"))
    }
}

################################################################
### ORDER COLUMNS

sample_count <- length(cna_mat[,1])

dim(mat)
mat <- mat[ order(apply(mat, 1, function(x) sum(x!="")) , decreasing=TRUE), ]
l <- apply(mat, 1, function(x) sum(x!="")) >= min_recurrence
mat <- mat[l, ]
dim(mat)

mat_order <- oncoprintOrder(mat,grouping_vector)

################################################################
### PLOT

outfile = paste(input_table,version,"pdf", sep=".")
print(outfile)

w <- (length(colnames(mat)) * .2) + 4
h <- (length(rownames(mat)) * .25) + 8.5

pdf(file=outfile, width = w, height = h)

chrArmLegend = Legend(at = c("Deletion", "Homozygous deletion", "Gain", "High gain", "LOH"), 
                      title = "ChrArmLevelCNVs", 
                      type = "points" , 
                      pch=c(25,25, 24, 24, 20), 
                      legend_gp = gpar(col = c("#0066CC", "#000066", "#CC6600", "#663300", "#006600"), 
                                       fill=c("#0066CC", "#000066", "#CC6600", "#663300", "#006600")
                                      ))

hm <- oncoPrint(mat, 
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, 
          col = col,
          column_title = title,
          row_order=NULL,
          right_annotation=NULL, # if you are using version after Dec. 2018
          #show_row_barplot=F, # if you are using version before Dec. 2018
          column_order=mat_order$sample_ind,
          show_column_names=TRUE,
          top_annotation=column_ha,
          show_pct = TRUE,
          width = 20,
)

draw(hm, annotation_legend_list = list(chrArmLegend), newpage=F)
#draw(hm, newpage=F)
#paste(rownames(mat))

# write table of exactly the matrix that was plotted:
outfile_table = paste(input_table,version,"plot","tsv", sep=".")
write.table(mat, outfile_table, sep='\t', append=FALSE, quote=FALSE, dec=".", row.names = TRUE, col.names=NA)

################################################################
### DECORATE

rownames(cna_mat)<-colnames(mat)
cna_mat2<-cna_mat[mat_order$sample_ind,]
cna_mat2[is.na(cna_mat2)]=""

# Annotate gain;loss chromosomes
for(i in (1:nrow(cna_mat2))){
  for(j in (1:ncol(cna_mat2))){
    if(grepl("del",cna_mat2[i,j])){
      decorate_annotation("cna", {grid.points(x=(i-0.5)/nrow(cna_mat2), y=(ncol(cna_mat2)-j+0.5)/ncol(cna_mat2), gp=gpar(col="#0066CC", fill="#0066CC"), pch=25, size=unit(1, "snpc")*0.07, default.units="npc")})
    }
    if(grepl("homoDel",cna_mat2[i,j])){
      decorate_annotation("cna", {grid.points(x=(i-0.5)/nrow(cna_mat2), y=(ncol(cna_mat2)-j+0.5)/ncol(cna_mat2), gp=gpar(col="#000066", fill="#0066CC"), pch=25, size=unit(1, "snpc")*0.07, default.units="npc")})
    }
    if(grepl("amp",cna_mat2[i,j])){
      decorate_annotation("cna", {grid.points(x=(i-0.5)/nrow(cna_mat2), y=(ncol(cna_mat2)-j+0.5)/ncol(cna_mat2), gp=gpar(col="#CC6600", fill="#CC6600"), pch=24, size=unit(1, "snpc")*0.07, default.units="npc")})
    }
    if(grepl("highAmp",cna_mat2[i,j])){
      decorate_annotation("cna", {grid.points(x=(i-0.5)/nrow(cna_mat2), y=(ncol(cna_mat2)-j+0.5)/ncol(cna_mat2), gp=gpar(col="#663300", fill="#663300"), pch=24, size=unit(1, "snpc")*0.07, default.units="npc")})
    }
    if(grepl("LOH",cna_mat2[i,j])){
      decorate_annotation("cna", {grid.points(x=(i-0.5)/nrow(cna_mat2), y=(ncol(cna_mat2)-j+0.5)/ncol(cna_mat2), gp=gpar(col="#006600"), pch=20, size=unit(1, "snpc")*0.07, default.units="npc")})
    }
  }
}

percentage_cna_vector<-NULL

for (i in 1:cnas_num){
  percentage_cna_vector[i]  <- as.integer(100*((sample_count - sort(na_count)[i])/sample_count))
}

# Annotate Chr arm percentation
decorate_annotation("cna", {grid.text(paste(percentage_cna_vector, "%", sep = ""), unit(-3, "mm"), (cnas_num:1-0.5)/cnas_num, just ="right")})

# Annotate sample info headers
decorate_annotation("Total SNVs",{grid.lines(c(0, 1), unit(c(median(annot$"SNV - total"), median(annot$"SNV - total")), "native"), gp = gpar(lty = 2, col = "#000000"))})
decorate_annotation("Total INDELs",{grid.lines(c(0, 1), unit(c(median(annot$"INDEL - total"), median(annot$"INDEL - total")), "native"), gp = gpar(lty = 2, col = "#000000"))})
sv_median <- median(sv_annot)
decorate_annotation("Total SVs",{grid.lines(c(0, 1), unit(c(sv_median,sv_median), "native"), gp = gpar(lty = 2, col = "#000000"))})

dev.off()
