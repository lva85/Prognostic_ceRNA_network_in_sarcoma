### set workspace
setwd("~/SARC_proj/")

##load gene expression data
if(T){
	gene_exp.mat = read.table("TCGA_SARC_geneExp.txt",header = T,sep = "\t")
	gene_exp.df = gene_exp.mat[!duplicated(gene_exp.mat$gene),]
	row.names(gene_exp.df) = gene_exp.df$gene
	gene_exp.df = gene_exp.df[,-c(1:2)]
	#retrieve the tumor samples
	gene_exp.tumor.df = gene_exp.df[,which(substring(colnames(gene_exp.df),14,15) == "01")]
	#normalization
	gene_exp.tumor.df = gene_exp.tumor.df[which(apply(gene_exp.tumor.df,1,sum)>0),]
	gene_exp.tumor.df.normalized <- t(scale(t(gene_exp.tumor.df)))

}

##load survival info
if(T){
	survival.info = read.table("data_bcr_clinical_data_patient.txt",header = T,sep = "\t")
    length(unique(gsub("-",".",survival.info$PATIENT_ID)))
    # 261
    length(unique(substring(colnames(miRNA_exp.df),1,12)))
    # 259
    length(intersect(unique(gsub("-",".",survival.info$PATIENT_ID)),unique(substring(colnames(miRNA_exp.df),1,12))))
    # [1] 259
}

## DEG
if(T){
  library(limma)
  #cancer=1, normal=0. 
  condition_table <- result.mat
  condition_table_for_limma <- model.matrix(~0+ condition_table$result)
  colnames(condition_table_for_limma) <- gsub("condition_table\\$result","", colnames(condition_table_for_limma))
  lmfit <- lmFit(gene_exp.tumor.df.normalized, condition_table_for_limma)
  
  contrast.matrix<-makeContrasts(cluster2-cluster1,levels=condition_table_for_limma)
  fit2<-contrasts.fit(lmfit,contrast.matrix)
  fit2<-eBayes(fit2)
  # fit2$coefficients
  tT=topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1)
  DEGs_cluster2TO1 <- tT[which(abs(tT$logFC) > 1 & tT$adj.P.Val < 0.01),]
  dim(DEGs_cluster2TO1)
  # [1] 2365    6
  
  DEGs_cluster2TO1.Mat <- data.frame("gene"=row.names(DEGs_cluster2TO1),DEGs_cluster2TO1)
  
  DEGs_cluster2TO1.df <- gene_exp.tumor.df.normalized[c(row.names(DEGs_cluster2TO1.Mat)),]
  DEGs_cluster2TO1.df.z3 <- as.matrix(DEGs_cluster2TO1.df)
  range(DEGs_cluster2TO1.df.z3)
  dim(DEGs_cluster2TO1.df.z3)
  
  contrast.matrix<-makeContrasts(cluster3-cluster1,levels=condition_table_for_limma)
  fit2<-contrasts.fit(lmfit,contrast.matrix)
  fit2<-eBayes(fit2)
  # fit2$coefficients
  tT=topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1)
  DEGs_cluster3TO1 <- tT[which(abs(tT$logFC) > 1 & tT$adj.P.Val < 0.01),]
  dim(DEGs_cluster3TO1)
  # [1] 15  6
  
  DEGs_cluster3TO1.Mat <- data.frame("gene"=row.names(DEGs_cluster3TO1),DEGs_cluster3TO1)
  
  DEGs_cluster3TO1.df <- gene_exp.tumor.df.normalized[c(row.names(DEGs_cluster3TO1.Mat)),]
  DEGs_cluster3TO1.df.z3 <- as.matrix(DEGs_cluster3TO1.df)
  range(DEGs_cluster3TO1.df.z3)
  
  contrast.matrix<-makeContrasts(cluster2-cluster3,levels=condition_table_for_limma)
  fit2<-contrasts.fit(lmfit,contrast.matrix)
  fit2<-eBayes(fit2)
  tT=topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1)
  DEGs_cluster2TO3 <- tT[which(abs(tT$logFC) > 1 & tT$adj.P.Val < 0.01),]
  dim(DEGs_cluster2TO3)
  # [1] 380   6
  
  DEGs_cluster2TO3.Mat <- data.frame("gene"=row.names(DEGs_cluster2TO3),DEGs_cluster2TO3)
  
  DEGs_cluster2TO3.df <- gene_exp.tumor.df.normalized[c(row.names(DEGs_cluster2TO3.Mat)),]
  DEGs_cluster2TO3.df.z3 <- as.matrix(DEGs_cluster2TO3.df)
  range(DEGs_cluster2TO3.df.z3)
  
  
  library(ComplexHeatmap)
  library(dendextend)
  library(RColorBrewer)
  library(circlize)
  DEGs_total.z3 = rbind(DEGs_cluster2TO1.df.z3,
						DEGs_cluster3TO1.df.z3,
						DEGs_cluster2TO3.df.z3)
  split.total = c(rep("C2vsC1",nrow(DEGs_cluster2TO1.df.z3)),
				  rep("C3vsC1",nrow(DEGs_cluster3TO1.df.z3)),
				  rep("C2vsC3",nrow(DEGs_cluster2TO3.df.z3)))
  
  annot_df <- data.frame(Cluster = result.mat$result)
  annot_df$Cluster = gsub("cluster1","C1(Low)",annot_df$Cluster)
  annot_df$Cluster = gsub("cluster2","C2(High)",annot_df$Cluster)
  annot_df$Cluster = gsub("cluster3","C3(Moderate)",annot_df$Cluster)
  row.names(annot_df) = row.names(result.mat)
  col = list(Cluster = c("C1(Low)" = "green", "C2(High)" = "red", "C3(Moderate)" = "blue"))
  ha <- HeatmapAnnotation("Cluster" = annot_df$Cluster, 
						  col = col,which = "column")
  
  Heatmap(DEGs_total.z3, name = "colorKey",
		  colorRamp2(c(-2, 0, 2), c("blue", "gray", "darkred")),
		  show_row_names = FALSE,show_column_names = FALSE,show_row_dend=F,
		  column_dend_height = unit(3, "cm"),cluster_columns = dend,
		  split = split.total,cluster_rows=F,top_annotation = ha
  )
  
}

## 10-year OS
if(T){
        survival.info$PATIENT_ID = gsub("-",".",survival.info$PATIENT_ID)
        substring(colnames(gene_exp.tumor.df),1,12)
        survival.info$OS_STATUS = as.numeric(substring(survival.info$OS_STATUS,1,1))
        survival.info$OS_STATUS = substring(survival.info$OS_STATUS,1,1)
        survival.info$OS_STATUS[which(survival.info$OS_STATUS == "[")] = NA
        survival.info$OS_STATUS = as.numeric(survival.info$OS_STATUS)
        survival.info$OS_MONTHS[grep("Not",survival.info$OS_MONTHS)] = NA
        survival.info$OS_MONTHS = as.numeric(survival.info$OS_MONTHS)
        survival.info$cluster = result.mat$result[match(survival.info$PATIENT_ID,substring(result.mat$sample,1,12))]
        survival.info = survival.info[!is.na(survival.info$cluster),]
        
        library("survival")
        library("survminer")
        
        survival.info.year1 = survival.info[which(survival.info$PHARMACEUTICAL_TX_ADJUVANT == "NO"),]
        survival.info.year1 = survival.info.year1[which(survival.info.year1$OS_MONTHS<120),]
        
        fit.cluster <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ cluster, data=survival.info.year1)
        
        p= ggsurvplot(fit.cluster, data = survival.info.year1,
                     pval = TRUE,
                     pval.coord = c(0, 0.03),
                     # surv.median.line = "hv", #添加中位生存曲线
                     palette=c("green", "red","blue"),  #更改线的颜色
                     legend.labs=c(paste0("C1: Low (n=",nrow(survival.info.year1[survival.info.year1$cluster == "cluster1",]),")"),
                                   paste0("C2: High (n=",nrow(survival.info.year1[survival.info.year1$cluster == "cluster2",]),")"),
                                   paste0("C3: Moderate (n=",nrow(survival.info.year1[survival.info.year1$cluster == "cluster3",]),")")), #标签
                     # legend.title="Treatment", 
                     title=paste0("Overall survival for subtypes "), #标题
                     # ylab="Disease Free survival (%)",
                     ylab="Overall survival (%)",
                     xlab = " Time (Months)", #更改横纵坐标
                     # risk.table = TRUE,
                     censor.shape = 124,censor.size = 2,conf.int = FALSE, #删失点的形状和大小
                     font.main = c(13, "bold", "black"),
                     font.x = c(12, "bold", "black"),
                     font.y = c(12, "bold", "black"),
                     break.x.by = 40)
        print(p, newpage = FALSE)
        
        
        
      }

## expression of CD1 family genes
if(T){
  #validate the gene-family related to immune
  #CD1 gene-family
  genef1 = DEGs_total.z3[startsWith(row.names(DEGs_total.z3),"CD"),]
  genef1 <- genef1[grep("CD[0-9]",row.names(genef1)),]
  genef1 <- genef1[which(!(row.names(genef1) %in% c("CD300LF","CD200R1","CD163L1"))),]
  genef1 = genef1[!duplicated(row.names(genef1)),]
  genef1 <- data.frame("gene"=row.names(genef1),genef1)
  library(reshape2)
  genef1 <- melt(genef1,id.vars=c("gene"))
  genef1$group <- ifelse(genef1$variable %in% result.mat$sample[which(result.mat$result == "cluster2")],"C2(high)",
						 ifelse(genef1$variable %in% result.mat$sample[which(result.mat$result == "cluster3")],"C3(moderate)","C1(low)") )
  
  library(ggpubr)
  
  genef1.p = ggboxplot(genef1, "gene", "value", color = "group",outlier.size=0.2,
					   palette = c("green", "red","blue"),legend.title = "Subtypes",
					   ylab = F,xlab = FALSE,size=0.001,width=0.5)+
	coord_flip()+
	scale_y_continuous(expand = c(0,0))+
	stat_compare_means(aes(group = group), label = "p.signif")
  print(genef1.p)
  
}

## expression of IL1 family genes
if(T){
  genef2 <- DEGs_total.z3[startsWith(row.names(DEGs_total.z3),"IL"),]
  genef2 = genef2[!duplicated(row.names(genef2)),]
  genef2 <- data.frame("gene"=row.names(genef2),genef2)
  genef2 <- melt(genef2,id.vars=c("gene"))
  genef2$group <- ifelse(genef2$variable %in% result.mat$sample[which(result.mat$result == "cluster2")],"C2(high)",
						 ifelse(genef2$variable %in% result.mat$sample[which(result.mat$result == "cluster3")],"C3(moderate)","C1(low)") )
  
  
  
  library(ggpubr)
  genef2.p = ggboxplot(genef2, "gene", "value", color = "group",outlier.size=0.2,
					   palette = c("green", "red","blue"),legend.title = "Subtypes",
					   add.params=list(size=0.1),dot.size=0.1,point.size=0.1,
					   ylab=FALSE,xlab=FALSE,size=0.01,width=0.5)+
	coord_flip()+
	scale_y_continuous(expand = c(0,0))+
	stat_compare_means(aes(group = group), label = "p.signif")
  
  print(genef2.p)
 
}

##help function for boxplot
BoxplotofClusters=function(tmp){
	temp.df = data.frame("gene"=row.names(tmp),tmp)
	library(reshape2)
	temp.df = melt(temp.df,id.vars=c("gene"))
	temp.df$group <- ifelse(temp.df$variable %in% result.mat$sample[which(result.mat$result == "cluster2")],"C2(high)",
							ifelse(temp.df$variable %in% result.mat$sample[which(result.mat$result == "cluster3")],"C3(moderate)","C1(low)") )
	return(temp.df)

}

## use xCell to validate the immune cell infiltration in the classification above
if(T){
	# install.packages("remotes")
	# remotes::install_github("icbi-lab/immunedeconv")

	library(ggplot2)
	library(immunedeconv)
	library(tidyverse)


	immunedeconv_xcell = deconvolute_xcell(gene_expression_matrix = gene_exp.tumor.df.normalized ,arrays = TRUE)

	#plot
	immunedeconv_xcell.df = BoxplotofClusters(immunedeconv_xcell)
	
	library(ggpubr)
	ggboxplot(immunedeconv_xcell.df, "gene", "value", color = "group",outlier.size=0.2,
			palette = c("green", "red","blue"),
			ylab = "Estimated Fraction",xlab = FALSE,size=0.001,width=0.5)+
	rotate_x_text(angle=30,hjust=1,vjust=1)+
	stat_compare_means(aes(group = group), label = "p.signif")
	

	# remove unrelated cell subtypes, i.e., Astrocytes, Hepatocytes, HSC and scores in the result
	remove.list = c("Astrocytes","Hepatocytes","HSC","ImmuneScore","StromaScore","MicroenvironmentScore")
	immunedeconv_xcell.df2 = immunedeconv_xcell.df[which(!(immunedeconv_xcell.df$gene %in% remove.list)),]
	
	library(ggpubr)
	ggboxplot(immunedeconv_xcell.df2, "gene", "value", color = "group",outlier.size=0.2,
			palette = c("green", "red","blue"),
			ylab = "Estimated Fraction",xlab = FALSE,size=0.001,width=0.5)+
	rotate_x_text(angle=30,hjust=1,vjust=1)+
	stat_compare_means(aes(group = group), label = "p.signif")
	

	# MicroenvironmentScore
	immunedeconv_xcell.df3 = immunedeconv_xcell.df[which(immunedeconv_xcell.df$gene %in% "MicroenvironmentScore"),]
	my_comparisons <- list(c("C1(low)", "C3(moderate)"), c("C1(low)", "C2(high)"), c("C3(moderate)", "C2(high)"))
	library(ggpubr)
	MicroenvironmentScore.p = ggviolin(immunedeconv_xcell.df3,x="group",y="value",
									 select=c("C1(low)","C3(moderate)","C2(high)"),
									 order=c("C1(low)","C3(moderate)","C2(high)"),
									 fill = "group",
									 color = "group",
									 width=1,xlab = FALSE,
									 ylab="MicroenvironmentScore",
									 palette = c("green", "blue","red"),
									 add = "boxplot", 
									 add.params = list(fill = "white",shape = "group",alpha=1), 
									 shape = "group")+
	stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
	stat_compare_means(label.y = 4.5,label.x = 0.9)
	print(MicroenvironmentScore.p)
	
}

## ESTIMATE: to calculate the Stromalscore, immunescore, ESTIMATE SCORE
if(T){
  library(utils)
  # rforge <- "http://r-forge.r-project.org"
  # install.packages("estimate", repos=rforge, dependencies=TRUE)
  library(estimate)
  
  write.table(gene_exp.tumor.df.normalized,"gene_exp.tumor.df.normalized.input.txt",sep = "\t",quote = F)
  input.exp = "gene_exp.tumor.df.normalized.input.txt"
  filterCommonGenes(input.f=input.exp, 
					output.f="gene_exp.tumor.df.normalized.input.commonGenes.gct", 
					id="GeneSymbol")
  
  
  estimateScore(input.ds = "gene_exp.tumor.df.normalized.input.commonGenes.gct",
				output.ds="gene_exp.tumor.df.normalized.input_estimate_score.gct", 
				platform="illumina")
  
  estimate_scores=read.table("gene_exp.tumor.df.normalized.input_estimate_score.gct",skip = 2,header = T)
  row.names(estimate_scores) = estimate_scores$NAME
  estimate_scores2=t(estimate_scores[,3:ncol(estimate_scores)])
  
  #plot
  estimate_scores2.df = BoxplotofClusters(t(estimate_scores2))
  
  library(ggpubr)
  estimate_scores.p = ggboxplot(estimate_scores2.df, "gene", "value", color = "group",outlier.size=0.2,
								palette = c("green", "red","blue"),
								ylab = FALSE,xlab = FALSE,size=0.001,width=0.5)+
	rotate_x_text(angle=0,hjust=0.5,vjust=1)+
	stat_compare_means(aes(group = group), label = "p.signif")
  print(estimate_scores.p)
  
  
  TumorPurity = cos(0.6049872018+0.0001467884 * estimate_scores2[,3])
  TumorPurity.df = data.frame(cbind(TumorPurity = as.numeric(TumorPurity), Sample = names(TumorPurity)))
  TumorPurity.df$group <- ifelse(TumorPurity.df$Sample %in% result.mat$sample[which(result.mat$result == "cluster2")],"C2(high)",
								 ifelse(TumorPurity.df$Sample %in% result.mat$sample[which(result.mat$result == "cluster3")],"C3(moderate)","C1(low)") )
  TumorPurity.df$TumorPurity = as.numeric(TumorPurity.df$TumorPurity)
  
  
  my_comparisons <- list(c("C1(low)", "C3(moderate)"), c("C1(low)", "C2(high)"), c("C3(moderate)", "C2(high)"))
  TumorPurity.p = ggviolin(TumorPurity.df,x="group",y="TumorPurity",
						   select=c("C1(low)","C3(moderate)","C2(high)"),
						   order=c("C1(low)","C3(moderate)","C2(high)"),
						   fill = "group",
						   color = "group",
						   width=1,xlab = FALSE,
						   palette = c("green", "blue","red"),
						   add = "boxplot", 
						   add.params = list(fill = "white",shape = "group",alpha=0.8), 
						   shape = "group")+
	stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ 
	stat_compare_means(label.y = 1.5,label.x = 0.9)
  print(TumorPurity.p)
  
  
}

##investigate the tumor site of patients
if(T){
	survival.info.organized = read.table("survival.info4organizing.txt",header = T,sep = "\t")
	survival.info.organized = survival.info.organized[which(!is.na(survival.info.organized$SITE_OF_TUMOR_TISSUE)),]

	library(RColorBrewer)
	library(ggpubr)
	tumor.site = data.frame(table(survival.info.organized$SITE_OF_TUMOR_TISSUE))
	tumor.site$typeNum = paste0(tumor.site$Var1, " (", tumor.site$Freq,")")

	tumor.site.group = data.frame(table(survival.info.organized$cluster,survival.info.organized$SITE_OF_TUMOR_TISSUE))
	colnames(tumor.site.group) = c("Cluster","Tumor_Site","Num")

	tumor.site.group = survival.info.organized[,c("cluster","SITE_OF_TUMOR_TISSUE")]
	tumor.site.group$cluster = ifelse(tumor.site.group$cluster == "cluster1","C1(low)",
									ifelse(tumor.site.group$cluster == "cluster2","C2(High)","C3(Moderate)"))

	tumor.site2 = tumor.site
	colnames(tumor.site2)[1] = c("Tumor_Site")
	tumor.site.group2 = merge(tumor.site.group,tumor.site2,by="Tumor_Site")
	tumor.site.group2$subtype = ifelse(tumor.site.group2$Cluster == "cluster1","C1(Low)(69)",
									 ifelse(tumor.site.group2$Cluster == "cluster2","C2(High)(68)","C3(Moderate)(121)"))
	tumor.site.group2$patient = rep(paste0("SARC(",sum(tumor.site.group2$Num),")"),nrow(tumor.site.group2))
	colnames(tumor.site.group2)
	aggregate(tumor.site.group2$Num, by=list(tumor.site.group2$Cluster),sum) 

	library(ggalluvial)
	colors=c( "#3C5488B2","#00A087B2", 
			"#F39B7FB2","#91D1C2B2", 
			"#8491B4B2", "#DC0000B2", 
			"#7E6148B2","yellow", 
			"darkolivegreen1", "lightskyblue", 
			"darkgreen",  "brown1", "darkorange1", 
			"cyan1", "royalblue4", "darksalmon", 
			"darkgoldenrod1", "darkseagreen", "darkorchid")

	ggplot(as.data.frame(tumor.site.group2),
		 aes(y = Num, axis1 = patient, axis2 = typeNum, axis3 = subtype)) +
	geom_alluvium(aes(fill = typeNum),curve_type = "quintic", alpha = 0.3,
				  width = 0, knot.pos = 0, reverse = FALSE) +
	guides(fill = F) +
	geom_stratum(width = 1/8, reverse = FALSE,alpha = 1,lwd=0.2) +
	geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2,
			  reverse = FALSE) +
	labs(x="",y="")+
	scale_fill_manual(values = colors)+
	theme_bw()+
	theme(axis.text.x = element_blank(),
		  axis.text.y = element_blank(),
		  axis.title.x = element_text(face="bold", color ="Black", size = 6, angle = 0, hjust = .5, vjust = 0.5),
		  axis.title.y = element_text(face="bold",color ="Black", size = 6, angle = 90, hjust = 0.5, vjust = 1.5),
		  legend.position="none",
		  legend.text = element_text(color ="Black", size = 14, angle = 0, hjust = 0, vjust = 0.5),
		  legend.title = element_text(face="bold",color ="Black", size = 14, angle = 0, hjust = .5, vjust = .5),
		  strip.text = element_text(face="bold", color ="Black", size = 14),
		  axis.ticks = element_blank(),
		  panel.background = element_blank(),
		  panel.grid = element_blank(),
		  panel.border = element_blank() )

}





