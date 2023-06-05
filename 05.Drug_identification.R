## set workspace
setwd("~/SARC_proj/")

load("classification.RData")

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

### CMAP analysis
# BiocManager::install("qusage")
# BiocManager::install("DrInsight")
library("DrInsight")

cmap.ref.profiles = get.cmap.ref(cmap.data.path = '14.Drug_prediction/rankMatrix.txt', 
							   probe.to.genes = probe.to.genes, drug.info = drug.info)

data("brca.tcga")

drug.ident.res = drug.ident(query.data = data.frame("geneSymbol"=row.names(Low2High.DE.res),
												  "score" = Low2High.DE.res$B), 
						  cmap.ref.profiles = cmap.ref.profiles,
						  repurposing.unit = "treatment", connectivity = "negative")

drug.ident.res.sig = drug.ident.res$drug.pvals[which(drug.ident.res$drug.pvals$FDR < 0.05),]

drug.ident.res.sig$pval = format(drug.ident.res.sig$pval,scientific = T,digits = 4)
drug.ident.res.sig$FDR = format(drug.ident.res.sig$FDR,scientific = T,digits = 4)


##analysis the drug sensitivity
if(T){
library(data.table)
drug.response <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)

drug.response.TriA = data.frame("Sample" = drug.response$V1,"Trichostatin_A"=drug.response$`Trichostatin A_437`)
drug.response.TriA$Group = result.mat$result[match(drug.response.TriA$Sample,result.mat$sample)]
drug.response.TriA$Group = gsub("cluster1","C1(Low)",drug.response.TriA$Group)
drug.response.TriA$Group = gsub("cluster2","C2(High)",drug.response.TriA$Group)
drug.response.TriA$Group = gsub("cluster3","C3(Moderate)",drug.response.TriA$Group)

my_comparisons <- list(c("C1(Low)", "C2(High)"), c("C2(High)", "C3(Moderate)"), c("C1(Low)", "C3(Moderate)"))

drug.response.TriA.p = ggplot(data = drug.response.TriA,
							  aes(y = Trichostatin_A,
								  x = Group))+
  geom_boxplot(alpha = 0.3,
			   fill = c('#47a4e9','#e94753',"#F39B7FB2"))+
  theme_bw()+
  ylab('Predicted Trichostatin A Sensitivity') +
  xlab('') +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox.test", hide.ns = TRUE)+
  stat_compare_means(label.y = 0.7,label.x = 0.7) +
  theme(legend.position="none") + 
  theme(axis.text.x=element_text(face="bold",color="Black", size=6,hjust=0.5,vjust = 0),
		axis.text.y = element_text(color ="Black", size = 6, angle = 0, hjust = 0.5, vjust = 0.5),
		axis.title.x = element_text(face="bold", color ="Black", size = 8, angle = 0, hjust = .5, vjust = 0.5),
		axis.title.y = element_text(face="bold",color ="Black", size = 8, angle = 90, hjust = 0.5, vjust = 1.5),
		axis.line = element_line(colour = "black"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		plot.background = element_blank(),
		panel.border = element_blank(),
		strip.text = element_text(face="bold", color ="Black", size = 8, angle = 0, hjust = .5, vjust = 0.5),
		plot.title.position="plot"
  )
print(drug.response.TriA.p)
