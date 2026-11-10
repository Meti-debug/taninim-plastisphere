# Chromagar colonies
#
# CHROMagar colony colours and total growth, stacked bars per river.
#
# Extracted from r/figures_as_run.R, lines 1-79, which is the session code
# deposited verbatim. Paths were made relative; nothing else was changed.
#
# Run from the repository root:
#   Rscript r/01_chromagar_colonies.R

source("r/00_paths.R")

#Figure 2. CHROMagar colony morphology, total growth, and ESBL confirmation across sampled rivers. (A)  CHROMagar ESBL Stacked bar chart showing the proportional distribution of colony colors for each sampled river. (B) CHROMagar Ampicillin Stacked bar chart showing the proportional distribution of colony colors on medium. Numbers inside segments indicate raw colony counts where n > 10. Colony colors correspond to putative bacterial identities based on CHROMagar manufacturer guidelines: metallic blue and turquoise colonies are indicative of Escherichia coli and related Enterobacteriaceae; pink colonies suggest Klebsiella, Enterobacter, or Citrobacter spp.; purple colonies are associated with other Gram-negative organisms; yellow tinge colonies indicate non-fermenting Gram-negatives. (C) Total colony counts comparing growth on CHROMagar (total bacteria), CHROMagar Ampicillin (ampicillin-resistant bacteria), and CHROMagar ESBL (ESBL-producing bacteria) across sampling locations. (D) ESBL confirmation from river water isolates. Oxidase-positive colonies (teal) were excluded as non-Gram-negative Enterobacteriaceae. Of the remaining oxidase-negative candidates, confirmed ESBL-positive (red) and ESBL-negative (blue) isolates are shown. Taninim River yielded the highest proportion of confirmed ESBL-producing Gram-negative bacteria (24 of 25 oxidase-negative candidates), followed by Hadera River (6 of 9), while Yarkon and Alexander rivers yielded no confirmed ESBL producers.

library(ggplot2); library(dplyr); library(ggpubr)

rivers_order <- c("Yarkon River", "Alexander River", "Hadera River", "Taninim River")

amp_data <- data.frame(
  River = rep(rivers_order, each = 10),
  Color = rep(c("Pink","Turquoise","Metallic Blue","Clear","Purple",
                "White","Yellow Tinge","Purple & White Halo",
                "Dark Purple","Offwhite"), 4),
  Count = c(30,50,80,80,50,30,20,50,50,0,
            20,100,100,80,100,50,30,0,100,0,
            30,100,100,80,100,50,50,5,50,0,
            0,100,200,200,100,200,0,0,200,0))

esbl_data <- data.frame(
  River = rep(rivers_order, each = 8),
  Color = rep(c("Pink","Turquoise","Metallic Blue","Clear",
                "Purple","White","Yellow Tinge","Offwhite"), 4),
  Count = c(1,0,0,33,0,5,0,0,
            0,1,30,0,1,5,0,10,
            100,5,50,100,80,10,0,0,
            100,100,200,100,0,100,200,0))

colony_colors <- c(
  "Pink"="#FF69B4","Turquoise"="#00BFFF","Metallic Blue"="#2F4F7F",
  "Clear"="#E8E8D0","Purple"="#8B008B","White"="#F0F0F0",
  "Yellow Tinge"="#FFEA00","Purple & White Halo"="#BA55D3",
  "Dark Purple"="#4B0082","Offwhite"="#FAF0E6")

shared_legend <- c(
  "Pink"                = "Pink — Escherichia coli",
  "Turquoise"           = "Turquoise — Enterococcus / Serratia spp.",
  "Metallic Blue"       = "Metallic Blue — Klebsiella / Enterobacter / Citrobacter",
  "Clear"               = "Clear — Pseudomonas / Stenotrophomonas",
  "Purple"              = "Purple — Klebsiella / Enterobacter / Citrobacter",
  "White"               = "White — Acinetobacter / Others",
  "Yellow Tinge"        = "Yellow Tinge — Pseudomonas aeruginosa",
  "Purple & White Halo" = "Purple & White Halo — Citrobacter / Enterobacter",
  "Dark Purple"         = "Dark Purple — Other Gram-negative",
  "Offwhite"            = "Offwhite — Other")

make_stacked <- function(data) {
  data$Color <- factor(data$Color, levels = names(colony_colors))
  ggplot(data %>% filter(Count > 0),
         aes(x = factor(River, levels = rivers_order), y = Count, fill = Color)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.45) +
    scale_fill_manual(values = colony_colors, labels = shared_legend,
                      name = NULL, drop = FALSE) +
    labs(x = NULL, y = "Number of Colonies") +
    theme_classic(base_size = 13) +
    theme(
      text = element_text(color = "black"),
      axis.text.x = element_text(size = 11.5, face = "bold", angle = 30,
                                 hjust = 1, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title.y = element_text(face = "bold", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
      axis.line = element_line(color = "black", linewidth = 1.0),
      legend.position = "right",
      legend.text = element_text(size = 8.8, face = "bold.italic", color = "black"))
}

pESBL <- make_stacked(esbl_data)   # now panel A
pAmp  <- make_stacked(amp_data)    # now panel B

final <- ggarrange(pESBL, pAmp, ncol = 2, nrow = 1,
                   labels = c("A", "B"),
                   common.legend = TRUE, legend = "right",
                   font.label = list(size = 16, face = "bold", color = "black"))

ggsave("Figure2_ESBL_Amp_Final.png", final,
       width = 14, height = 6.5, dpi = 600, bg = "white")
ggsave("Figure2_ESBL_Amp_Final.pdf", final,
       width = 14, height = 6.5, device = cairo_pdf, bg = "white")

dim(png::readPNG("Figure2_ESBL_Amp_Final.png"))[2:1]

