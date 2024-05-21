
# dir.amp(ps0)
# theme.amp(theme = T)
# package.amp()
# color.amp()
package.amp <- function(){
  library(phyloseq)
  library(tidyverse)
  library(ggClusterNet)
  library(EasyStat)
  library(fs)
  library(ggthemes)
  library(RColorBrewer)
  library(magrittr)
  library(MicrobiotaProcess)
  library(ggsignif)
  library(ggtree)
  library(ggtreeExtra)
  library(ggstar)
  library(MicrobiotaProcess)
  library(ggnewscale)
  library(grid)
}

dir.amp <- function(ps0,
                    smart = FALSE){
  # # #--------
  result_path <- paste("./","/result_and_plot/",sep = "")
  fs::dir_create(result_path)
 
  if (smart) {
    #---results#---------
    tax.1 = c("Fungi",
              "fungi",
              "K__Fungi",
              "k__Fungi",
              "d__Fungi",
              "d__fungi",
              "d__Eukaryota",
              "Eukaryota",
              "K:Fungi",
              "k:Fungi",
              "d:Fungi",
              "d:fungi",
              "d:Eukaryota",
              "Eukaryota",
              "D:Eukaryota"
    )
    TFdir.f <- as.data.frame(table(
      phyloseq::tax_table(ps0)[,1]))[,2][as.data.frame(table(phyloseq::tax_table(ps0)[,1]))[,1] %in%
                                           tax.1] > 10
    
    tax.2 = c("Bacteria",
              "K__Bacteria",
              "k__Bacteria",
              "d__Bacteria",
              "k:Bacteria",
              "bacteria",
              "D:bacteria",
              "D__bacteria",
              "d__bacteria",
              "d:prokaryotes",
              "d__prokaryotes",
              "prokaryotes"
    )
    
    TFdir.b <- as.data.frame(table(phyloseq::tax_table(ps0)[,1]))[,2][as.data.frame(table(phyloseq::tax_table(ps0)[,1]))[,1] %in%
                                                                        tax.2 ] > 10
    
    if (length(TFdir.f) != 0) {
      print("ITS")
      res1path <- paste(result_path,"/Base_diversity_ITS",sep = "")
      id = tax.1
    }
    
    if (length(TFdir.b) != 0) {
      print("16s")
      res1path <- paste(result_path,"/Base_diversity_16s",sep = "")
      id = tax.2
    }
  }else {
    res1path <- paste(result_path,"/Data.mining",sep = "")
    id = "No check"
  }

  fs::dir_create(res1path)
  return(list(res1path,id))
  
}

#-

theme_my = function() {
  
  mytheme1 = ggplot2::theme_bw() + ggplot2::theme(
    panel.background=  ggplot2::element_blank(),
    panel.grid=  ggplot2::element_blank(),
    legend.position="right",
    legend.title =  ggplot2::element_blank(),
    legend.background= ggplot2::element_blank(),
    legend.key= ggplot2::element_blank(),
    plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =  ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.title.x = ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.text =  ggplot2::element_text(size = 20,face = "bold"),
    axis.text.x =  ggplot2::element_text(colour = "black",size = 14),
    axis.text.y =  ggplot2::element_text(colour = "black",size = 14),
    legend.text =  ggplot2::element_text(size = 15,face = "bold")
  )
  
  mytheme2 = ggplot2::theme_bw() + ggplot2::theme(
    panel.background=  ggplot2::element_blank(),
    panel.grid=  ggplot2::element_blank(),
    legend.position="bottom",
    
    legend.title =  ggplot2::element_blank(),
    legend.background=  ggplot2::element_blank(),
    legend.key= ggplot2::element_blank(),
    plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =  ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.title.x = ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.text =  ggplot2::element_text(size = 20,face = "bold"),
    axis.text.x =  ggplot2::element_text(colour = "black",size = 14,angle = 90,hjust = 1),
    axis.text.y =  ggplot2::element_text(colour = "black",size = 14),
    legend.text =  ggplot2::element_text(size = 15,face = "bold")
  )
  
  
  
  #--
  gnum <- unique(phyloseq::sample_data(ps)$Group) %>% length()
  
  # scales::show_col(RColorBrewer::brewer.pal(9,"Set1"))
    if (gnum < 10 ) {
      colset1 <- RColorBrewer::brewer.pal(9,"Set1")
    } else {
      colset1 <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(gnum)
    }
    ##------------

    colset2 <- RColorBrewer::brewer.pal(12,"Paired")
    colset3 <- c(RColorBrewer::brewer.pal(11,"Set1"),RColorBrewer::brewer.pal(9,"Pastel1"))
    colset4 = colset3
  

  return(list(mytheme1,mytheme2,colset1,colset2,colset3,colset4))
  
}


# mytheme1.1 = ggplot2::theme_bw() + ggplot2::theme(
#   panel.background=  ggplot2::element_blank(),
#   panel.grid=  ggplot2::element_blank(),
#   legend.position="right",
#   legend.title =  ggplot2::element_blank(),
#   legend.background= ggplot2::element_blank(),
#   legend.key= ggplot2::element_blank(),
#   plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
#   axis.title.y =  ggplot2::element_text(colour = "black"),
#   axis.title.x = ggplot2::element_text(),
#   axis.text =  ggplot2::element_text(),
#   axis.text.x =  ggplot2::element_text(),
#   axis.text.y =  ggplot2::element_text(),
#   legend.text =  ggplot2::element_text()
# )
# 
# mytheme2.1 = ggplot2::theme_bw() + ggplot2::theme(
#   panel.background=  ggplot2::element_blank(),
#   panel.grid=  ggplot2::element_blank(),
#   legend.position="right",
#   
#   legend.title =  ggplot2::element_blank(),
#   legend.background=  ggplot2::element_blank(),
#   legend.key= ggplot2::element_blank(),
#   plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
#   axis.title.y =  ggplot2::element_text(),
#   axis.title.x = ggplot2::element_text(),
#   axis.text =  ggplot2::element_text(),
#   axis.text.x =  ggplot2::element_text(angle = 90),
#   axis.text.y =  ggplot2::element_text(),
#   legend.text =  ggplot2::element_text()
# )

theme.col = function(ps = ps,method = "anhei",Top = 15,
                     gnum = 5
                     ) {
  
  # 
  mytheme1 = ggthemes::theme_base()
  mytheme1.1 = ggthemes::theme_few() 
  mytheme1.2 = ggthemes::theme_par()
  mytheme1.3 = theme_bw() 
  
  mytheme2 = ggthemes::theme_base() + ggplot2::theme(
    legend.position="bottom",
    axis.text.x =  ggplot2::element_text(angle = 90,hjust = 1),
  )
  mytheme2.1 = ggthemes::theme_few() + ggplot2::theme(
    legend.position="bottom",
    axis.text.x =  ggplot2::element_text(angle = 90,hjust = 1),
  )
  mytheme2.2 = ggthemes::theme_par()+ ggplot2::theme(
    legend.position="bottom",
    axis.text.x =  ggplot2::element_text(angle = 90,hjust = 1),
  )
  mytheme2.3 = theme_bw() + ggplot2::theme(
    legend.position="bottom",
    axis.text.x =  ggplot2::element_text(angle = 90,hjust = 1),
  )
  
  
  mytheme0.1 = ggplot2::theme_bw() + ggplot2::theme(
    panel.background=  ggplot2::element_blank(),
    panel.grid=  ggplot2::element_blank(),
    legend.position="right",
    legend.title =  ggplot2::element_blank(),
    legend.background= ggplot2::element_blank(),
    legend.key= ggplot2::element_blank(),
    plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =  ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.title.x = ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.text =  ggplot2::element_text(size = 20,face = "bold"),
    axis.text.x =  ggplot2::element_text(colour = "black",size = 14),
    axis.text.y =  ggplot2::element_text(colour = "black",size = 14),
    legend.text =  ggplot2::element_text(size = 15,face = "bold")
  )
  
  mytheme0.2 = ggplot2::theme_bw() + ggplot2::theme(
    panel.background=  ggplot2::element_blank(),
    panel.grid=  ggplot2::element_blank(),
    legend.position="right",
    
    legend.title =  ggplot2::element_blank(),
    legend.background=  ggplot2::element_blank(),
    legend.key= ggplot2::element_blank(),
    plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =  ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.title.x = ggplot2::element_text(size = 24,face = "bold",colour = "black"),
    axis.text =  ggplot2::element_text(size = 20,face = "bold"),
    axis.text.x =  ggplot2::element_text(colour = "black",size = 14,angle = 90,hjust = 1),
    axis.text.y =  ggplot2::element_text(colour = "black",size = 14),
    legend.text =  ggplot2::element_text(size = 15,face = "bold")
  )
  
  
  
  #--
  if (is.null(gnum)) {
    gnum <- unique(phyloseq::sample_data(ps)$Group) %>% length()
    
  }

  #
  if (gnum < 10 ) {
    
    
    colset1 <- RColorBrewer::brewer.pal(9,"Set1")
    colset1.1 <- RColorBrewer::brewer.pal(9,"Set3")
    colset1.2 <- RColorBrewer::brewer.pal(8,"Set2")
    colset1.3 <- RColorBrewer::brewer.pal(8,"Dark2")
    colset1.4 <- RColorBrewer::brewer.pal(8,"Accent")
    
    col.time.r = RColorBrewer::display.brewer.pal(9,"YlOrRd")
    col.time.b = RColorBrewer::display.brewer.pal(9,"YlGnBu")
    
    if (method =="anhei") {
      colset1 <- ggsci::pal_aaas(alpha = 1)(9)
      colset1.1 <-ggsci::pal_jco()(9)
    } else{
      
    }
    
    
    
  } else {
    
    colset1 <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(gnum)
    colset1.1 <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set3"))(gnum)
    colset1.2 <- colorRampPalette(RColorBrewer::brewer.pal(8,"Set2"))(gnum)
    colset1.3 <- colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(gnum)
    colset1.4 <- colorRampPalette(RColorBrewer::brewer.pal(8,"Accent"))(gnum)
    
  }
  
  col.group = list(colset1,colset1.1,colset1.2,colset1.3,colset1.4)
  
  
  # 
  col.time.r = colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(gnum)
  col.time.b = colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(gnum)
  col.time.rb = colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(gnum)
  col.time.ob = colorRampPalette(RColorBrewer::brewer.pal(11,"BrBG"))(gnum)
  col.time.rblack = colorRampPalette(RColorBrewer::brewer.pal(11,"RdGy"))(gnum)
  col.time.jd  <-ggsci::rgb_gsea(n = gnum)
  
  col.time = list(col.time.r,col.time.b,col.time.rb,col.time.ob,col.time.rblack,col.time.jd)
  
  # 
  colset2.bars1 <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(Top)
  colset2.bars2 <- colorRampPalette(c(RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(9,"Pastel1")))(Top)
  colset2.bars3 <- colorRampPalette(ggsci::pal_futurama()(11))(Top)
  colset2.bars3 <- ggsci::pal_igv()(Top)
  col.bar = list(colset2.bars1,colset2.bars2,colset2.bars3)
  
  mytheme = list(mytheme1,
                 mytheme1.1,
                 mytheme1.2,
                 mytheme1.3,
                 mytheme2,
                 mytheme2.1,
                 mytheme2.2,
                 mytheme2.3,
                 mytheme0.1,
                 mytheme0.2
  )
  
  return(list(mytheme = mytheme,
              col.group =   col.group,
              col.time = col.time,
              col.bar = col.bar))
  
}




