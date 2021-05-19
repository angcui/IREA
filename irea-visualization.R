#### Immune response enrichment score visualization functions #####

#' Compass plot
#'
#' \code{IreaCompassPlot} Visualize IREA results on a compass plot
#'
#' @param    degs              A list of differentially expresssed genes that user would like to investigate
#' @param    input_celltype    Choose from one of the listed cell types that most resemble the input
#'
#'

IreaCompassPlot = function(df_irea, p_cutoff = 0.05,
                               color_by = c("pval", "ES")) {
  library(ggplot2)
  library(plyr)
  cytokines_spreadsheet = read.csv("~/Dropbox (Personal)/Hacohen/ligands/spreadsheets/cytokine_lists.csv")

  df_irea$nlog10p = -log10(df_irea$padj + 0.00001)
  df_irea$nlog10p[df_irea$nlog10p<0]=0.000000001
  `%notin%` = Negate(`%in%`)
  df_irea = subset(df_irea, Cytokine %notin% c("IL2+IL15", "IL2+IFNg"))

  df_irea$nlog10p_signed = df_irea$nlog10p * sign(df_irea$ES)
  max_barlimit = max(df_irea$ES)
  min_barlimit = min(df_irea$ES)
  text_angle = c(seq(0,-180,length.out = 43), seq(360, 180, length.out = 43))

  df_irea$Cytokine_Label = mapvalues(df_irea$Cytokine, from = cytokines_spreadsheet$Cytokine_OriginalName,
                                     to = cytokines_spreadsheet$Cytokine_DisplayName)
  df_irea$Cytokine_Label = factor(df_irea$Cytokine_Label, levels = df_irea$Cytokine_Label[order(df_irea$ES)])

  text_angle = c(seq(0,-180,length.out = 43), seq(360, 180, length.out = 43))
  ggplot(df_irea, aes(x = Cytokine_Label, y = ES)) +
    theme_void() +
    geom_segment(x = 1:86, xend = 1:86, y = 0, yend = max_barlimit, arrow = arrow(length=unit(0.15,"cm")), color = "gray30") +
    geom_bar(stat = "identity",size = 0.2, aes(fill = nlog10p_signed)) +
    ylim(min_barlimit, max_barlimit) +
    coord_polar() +
    scale_fill_gradient2(mid = "white", high = "red", low = "blue", midpoint = 0) +
    #    scale_color_gradient(low = "white", high = "black") +
    theme(axis.text.x=element_text(color = "black", size = 8,angle = text_angle+90),
          plot.title = element_text(hjust = 0.5))
}






