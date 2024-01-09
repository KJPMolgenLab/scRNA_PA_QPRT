require("kableExtra")
require("tidyverse")
require("compareGroups")
require("RColorBrewer")
require("stringr")
require("pheatmap")
require("DESeq2")



Dark8 = brewer.pal(8, "Dark2")
Dark8_50 = paste0(brewer.pal(8, "Dark2"), "7D")
jetcolors = colorRampPalette(c("darkblue", "skyblue", "green",
                               "yellow", "orange", "red", "darkred"))

OBcolors = colorRampPalette(c("darkblue", "skyb lue",
                              "white",  "orange", "darkorange3"))

display_tab = function(df){
  df %>% kbl(digits = 3,align = "l") %>%
    kable_classic(full_width = F, position = "left") %>%
    scroll_box(width = "900px", height ="500px",
               fixed_thead = T)
}


table_sumstat_grp = function(DF, columns, groupfactor){
  tmpdf= DF %>% select(all_of(c(columns, groupfactor)))
  table <- compareGroups(formula(paste0(groupfactor, "~.")), data = tmpdf)
  pvals <- getResults(table, "p.overall")
  export_table <- createTable(table)
  return(export_table)
}


PCAplot=function(PCA_res, Sampledf, label){
  plot(PCA_res[,])

}


# comparison fucntion of target
comparison <- function(dds_object, samples, target, randomeffect){
  require(DESeq2)
  require(limma)

  if(length(samples)==0){
    print("Sample length is 0 all samples included")
    samples = colnames(dds_object)
  }

  designform = as.formula(paste0("~ 1+",target))
  dds_filt = dds_object[,samples]
  ## no random effect
  if(length(randomeffect)==0){
    design(dds_filt) <- designform
    dds_filt <- DESeq2::DESeq(dds_filt)
    res = results(dds_filt)
    return(res)
  }
  ## with random effects
  if(length(randomeffect)==1){
    log_cpm=log2(counts(dds_filt, normalize=T)+1)
    design = model.matrix(designform, colData(dds_filt))
    rande = colData(dds_filt)[,randomeffect]
    dupcor <- duplicateCorrelation(log_cpm, design, block=rande)
    fitDupCor <- lmFit(log_cpm, design, block=rande, correlation=dupcor$consensus)
    fit<- eBayes(fitDupCor)
    topTable(fit,n=dim(fit)[1])
  }
}


# go profiler function
getGOresults = function(geneset, genereference, evcodes=FALSE){
  require(gprofiler2)
  resgo = gost(geneset, organism = "hsapiens",
               correction_method = "gSCS",
               domain_scope = "custom",
               evcodes=evcodes,
               sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "TF", "HP", "HPA"),
               custom_bg = genereference,
               numeric_ns = "ENTREZGENE_ACC")
  if(length(resgo) != 0){
    return(resgo)
  } else {
    print("no significant results")
    return(NULL)
  }
}


GOplot = function(GOtable, N, Title="GO plot", ylabel="GO term"){
  if(nrow(GOtable)<N){N=nrow(GOtable)}
  GOtable = GOtable[GOtable$parents!="character(0)",]
  Tabtoplot=GOtable[order(GOtable$p_value, decreasing = F)[1:N],]
  Tabtoplot$log10pvalue=-log10(Tabtoplot$p_value)
  Tabtoplot$genperc=Tabtoplot$intersection_size/Tabtoplot$effective_domain_size

  wrapit = function(long_phrase, cutoff){
    if(nchar(long_phrase) > cutoff){
      cutpos=ceiling(str_count(long_phrase, pattern = " ")/2)
      modx = gsub(paste0("(([^ ]* ){",cutpos,"})([^ ]*)"), "\\1\n\\3", long_phrase)
      return(modx)
    } else {
      return(long_phrase)}

  }

  Tabtoplot$term_name = sapply(Tabtoplot$term_name, wrapit, cutoff=40)

  ggplot(Tabtoplot) + geom_point(aes(x =log10pvalue,
                                     y = N:1,
                                     size=precision,
                                     colour=genperc)) +
    scale_color_viridis_c(limits = c(0,0.05))+
    labs(fill="genomic cov", size="precision")+
    xlab("- log10(p-value)") + ylab(ylabel)+
    scale_size_area(limits=c(0,1))+
    theme_bw(base_size = 12) + ggtitle(Title)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(breaks=N:1,
                       labels=Tabtoplot$term_name)+
    guides(size = guide_legend(order = 2),
           colour = guide_colorsteps(order = 1, barheight = 4))+
    xlim(0,15)

}


geneheatmap=function(GOIsEntrez,exprobj,CellID=NULL,
                     hm.cluster_rows = F,
                     hm.cluster_cols = F,
                     show_rownames = F,
                     show_colnames = F,
                     hm.scale="none",
                     annotation_row_labels = NA,
                     annotation_row_names = NA,
                     title = "log(count+1) gene expression"){
  CellIDorig=CellID
  SampleInfo = colData(exprobj) %>% as.data.frame()

  if(is.null(CellIDorig)){
    CellID = unique(SampleInfo$CellLine)
  } else {
    CellID = CellIDorig
  }


  if(length(names(GOIsEntrez))==0)
    names(GOIsEntrez) = as.character(GOIsEntrez)

  idx = match(GOIsEntrez, rownames(exprobj))
  log_2cpm=log2(counts(exprobj, normalize=T)+1)
  tmpsmpl = rownames(SampleInfo)[SampleInfo$CellLine %in% CellID]
  tmpsmpl = intersect(colnames(exprobj), tmpsmpl)
  idxsmpl = rownames(SampleInfo) %in% tmpsmpl


  dataset= log_2cpm[idx, rownames(SampleInfo)[idxsmpl]]

  colnames(dataset)= SampleInfo[rownames(SampleInfo)[idxsmpl],"label_rep"]
  rownames(dataset) = names(GOIsEntrez)

  #colors for plotting heatmap
  colors <- scales::viridis_pal()(255)


  gRNAcol = Dark8[c(1:nlevels(SampleInfo$gRNA))+nlevels(SampleInfo$CellLine)]
  names(gRNAcol) = levels(SampleInfo$gRNA)

  diffcol = brewer.pal(3,"Set1")[1:nlevels(SampleInfo$DIFF)]
  names(diffcol) = levels(SampleInfo$DIFF)

  rapacol = brewer.pal(3,"Set2")[1:nlevels(SampleInfo$RAPA)]
  names(rapacol) = levels(SampleInfo$RAPA)

  ann_colors = list(
    DIFF = diffcol,
    RAPA = rapacol,
    gRNA = gRNAcol)

  if(is.null(CellIDorig)){
  labels = SampleInfo[match(colnames(dataset), SampleInfo$label_rep),
                      c("RAPA", "gRNA", "DIFF")]
  } else {

  labels = SampleInfo[match(colnames(dataset), SampleInfo$label_rep),
                      c("RAPA", "gRNA", "DIFF",  "CellLine")]
  }
  rownames(labels)=colnames(dataset)
  labels = labels[order(as.numeric(labels$RAPA), as.numeric(labels$gRNA), as.numeric(labels$DIFF)),]

  labels = labels %>% mutate_all(as.character) %>% as.data.frame()

  #annotation_row_labels = c(15,20,24)
  tmp_annr = annotation_row_labels

  if(any(!is.na(annotation_row_labels))){
    for(i in annotation_row_labels){
      ann_colors[[colnames(mcols(exprobj))[i]]] = c("TRUE"="#000000", "FALSE"="#CCCCCC")
    }
    annotation_row_labels = mcols(exprobj)[rownames(dataset),annotation_row_labels] %>% as.data.frame()
    annotation_row_labels=annotation_row_labels %>% mutate_all(as.factor) %>% as.data.frame()

  }

  if(any(!is.na(annotation_row_names))){
    names(ann_colors) = c("DIFF", "RAPA", "gRNA", annotation_row_names)
    colnames(annotation_row_labels)=annotation_row_names
  }



  pheatmap(dataset[,rownames(labels)],
           show_rownames = show_rownames,
           show_colnames = show_colnames,
           border_color = NA,
           cluster_rows = hm.cluster_rows,
           cluster_cols = hm.cluster_cols,
           col = colors,
           fontsize=8,
           scale = hm.scale,
           annotation_col = labels,
           annotation_row = annotation_row_labels,
           annotation_colors = ann_colors,
           main=title)

}

samesign <- function(x) {abs(sum(sign(x)))==length(x)}

multiORplot = function(datatoplot=FALSE, Pval = "Pval", Padj = "Padj", SE = "SE", beta="beta", pheno = "pheno", xlims=NULL){
  starpval=convertpvaltostars(datatoplot[[Pval]])
  starpval[datatoplot[[Padj]]<0.05]="adj.p**"
  starpval[is.na(datatoplot[[beta]])]="n.a."
  CIUpper = datatoplot[[beta]] +1.96*datatoplot[[SE]]
  CILower = datatoplot[[beta]] -1.96*datatoplot[[SE]]
  if(length(xlims)==0){
    xlims=range(c(CIUpper, CILower), na.rm=T)*1.2
  }
  par(mar=c(5,16,5,2))
  betas = datatoplot[[beta]]

  plot(x=betas, y=1:length(betas),
       type="n", panel.first = grid(ny=NA),
       yaxt = "n", ylab="",
       xlim=xlims,
       xlab=expression(paste('log(OR)'%+-%95,"%CI")),
       main=paste(pheno))
  abline(v=0,col="black", lty=3)
  axis(2, at=1:length(betas),
       labels=base::rev(rownames(datatoplot)),
       las=1)
  arrows(x0=CILower, x1=CIUpper, y0=length(betas):1, y1=length(betas):1, col=rainbow(length(betas)), length=0, lwd=2,code = 3)
  points(y=length(betas):1, x=betas, pch=18, col="black")
  betas[is.na(betas)]=0
  text(y=(length(betas):1)+0.5, x=betas, labels=starpval, cex=0.7)
}

convertpvaltostars=function(x){
  sapply(x, function(x){ifelse(x<=0.01, "**", ifelse(x<=0.05, "*", ""))})
}

# plots eigenvalues of specified samples
EigengenePlot=function(data, Sampledata, samplesincl){
  Sampledata = Sampledata[samplesincl,]
  data=data[rownames(Sampledata),]
  for (i in colnames(data)){
    nf=layout(matrix(c(1:4,rep(5,4)),ncol=2),
              heights = c(12,1,1,1),
              widths = c(10,2))
    par(mar=c(0.2,4.1,3,1))
    lim=max(abs(data[,i]))
    barplot(data[,i], col=gsub("ME", "", i),border = NA, main=i, ylim=c(-lim,lim),
            ylab="ME expression")
    par(mar=c(0.1,4.1,0,1))
    a=barplot(rep(1,length(data[,i])),border = NA,
              col=ann_colors[["gRNA"]][Sampledata[,"gRNA"]], yaxt='n')
    a=barplot(rep(1,length(data[,i])),border = NA,
              col=ann_colors[["RAPA"]][Sampledata[,"RAPA"]], yaxt='n')
    a=barplot(rep(1,length(data[,i])),border = NA,
              col=ann_colors[["DIFF"]][Sampledata[,"DIFF"]], yaxt='n')
    #a=barplot(rep(1,length(data[,i])), border = NA,
    #          col=ann_colors[["CellLine"]][Sampledata[,"CellLine"]], yaxt='n')
    par(mar=c(0,0,0,0))
    plot.new()
    legend(0,0.5, legend = names(ann_colors[["gRNA"]]),fill = ann_colors[["gRNA"]], xpd=T,bty = "n")
    legend(0,0.3, legend = names(ann_colors[["RAPA"]]),fill = ann_colors[["RAPA"]], xpd=T,bty = "n")
    legend(0,0.2, legend = names(ann_colors[["DIFF"]]),fill = ann_colors[["DIFF"]], xpd=T,bty = "n")
    #legend(0,0.125, legend = names(ann_colors[["CellLine"]]),fill = ann_colors[["CellLine"]], xpd=T,bty = "n")
  }
}


