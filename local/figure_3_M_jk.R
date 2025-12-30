library(latex2exp)
require(tidyverse)

theme_Publication <- function(base_size=11, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            #axis.title = element_text(face = "bold",size = rel(1)),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            #panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            #legend.title = element_text(size = rel(1), face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

f <- '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/maf_bin_stratified_by_afr_maf_corrected_enrichment_meta.csv'
enrichment_stratified_by_afr_maf <- data.table::fread(f, data.table = F)


bin_type <- 'quintile'
tex_labels <- c(
    TeX('$p_a < .005, p_e$ varies'),
    TeX('$.005 \\leq p_a < .05, p_e$ varies'),
    TeX('$p_a \\geq .05, p_e$ varies'),
    TeX('$p_e \\geq .05, p_a$ varies')
)
lineplot <- enrichment_stratified_by_afr_maf %>%
    dplyr::filter(MAF > .05) %>%
    dplyr::mutate(
        per_allele_enrichment_upper = expectation_b2_normalized_meta + 1.96 * expectation_b2_normalized_meta_se,
        per_allele_enrichment_lower = expectation_b2_normalized_meta - 1.96 * expectation_b2_normalized_meta_se) %>%
    dplyr::mutate(snp_group = word(Category, 1, 1, fixed('('))) %>%
    dplyr::mutate(maf_type = word(Category, 1, 1, fixed('_'))) %>%
    dplyr::filter(str_detect(Category, bin_type)) %>%
    dplyr::filter(str_detect(Category, '_in_')) %>%
    dplyr::mutate(label = dplyr::case_when(
        startsWith(Category, 'afr_maf_quintile_in_eur_c_') ~ 'p_a >= .05$\np_a varies',
        startsWith(Category, 'eur_maf_quintile_in_afr_c_') ~ 'p_a >= .05\np_e varies',
        startsWith(Category, 'eur_maf_quintile_in_afr_lf_') ~ '.005 <= p_a < .05\np_e varies',
        startsWith(Category, 'eur_maf_quintile_in_afr_rare_') ~ 'p_e >= .05\np_e varies'
    )) %>%
    dplyr::mutate(label = factor(label, levels = c('p_e >= .05\np_e varies', '.005 <= p_a < .05\np_e varies', 'p_a >= .05\np_e varies', 'p_a >= .05$\np_a varies'))) %>%
    ggplot(aes(
        x = MAF, 
        ymin = per_allele_enrichment_lower, 
        ymax = per_allele_enrichment_upper, 
        y = expectation_b2_normalized_meta, 
        color = label, 
        linetype = label)) + 
    geom_point(size = 2, show.legend = FALSE) + 
    geom_linerange(show.legend = F, linetype = 1) + 
    geom_line() +
    theme_Publication(base_size = 10)  + 
    scale_color_manual(values = c("#377EB8","#984EA3", "#E41A1C", "#E41A1C"), labels = tex_labels) + 
    scale_linetype_manual(values = c(2, 2, 2, 1), labels = tex_labels) + 
    guides(
        color = guide_legend(title = '', position = 'bottom', nrow = 2, theme = theme(legend.byrow = TRUE)), 
        linetype = guide_legend(title = '', position = 'bottom', nrow = 2, theme = theme(legend.byrow = TRUE))) + 
    theme(legend.text=element_text(size=9)) + 
    labs(y = 'Meta-analyzed per-allele\n effect size variance', color = '', x = 'MAF (EUR or AFR)', linetype = 'MAF') + 
    scale_y_log10() + 
    scale_x_log10() +
    theme(plot.margin=unit(c(10,40,5,5),"mm"))

f <- '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/coarse_grid_enrichment_meta.csv'
maf_grid_enrichment <- f %>%
    data.table::fread(data.table = F) %>% 
    dplyr::mutate(
        eur_maf_bin = word(Category, 5, 5, '_'),
        afr_maf_bin = word(Category, 8, 8, '_'))

box_offset <- .5 - .054
box_width <- 1.5
grid_figure <- maf_grid_enrichment %>%
    dplyr::filter(Category != 'base') %>%
    dplyr::mutate(enrichment_rounded = round(expectation_b2_normalized_meta, 2)) %>%
    dplyr::mutate(
        eur_left_bound = as.numeric(str_sub(word(eur_maf_bin, 1, 1, fixed(',')), 2, -1)),
        eur_right_bound = as.numeric(str_sub(word(eur_maf_bin, 2, 2, fixed(',')), 1, -2)),
        afr_left_bound = as.numeric(str_sub(word(afr_maf_bin, 1, 1, fixed(',')), 2, -1)),
        afr_right_bound = as.numeric(str_sub(word(afr_maf_bin, 2, 2, fixed(',')), 1, -2))) %>% 
    dplyr::mutate(
        eur_left_bound = ifelse(eur_left_bound < .01, round(eur_left_bound, 3), round(eur_left_bound, 2)),
        afr_left_bound = ifelse(afr_left_bound < .01, round(afr_left_bound, 3), round(afr_left_bound, 2)),
        eur_right_bound = ifelse(eur_right_bound < .01, round(eur_right_bound, 3), round(eur_right_bound, 2)),
        afr_right_bound = ifelse(afr_right_bound < .01, round(afr_right_bound, 3), round(afr_right_bound, 2)),
        eur_maf_bin = paste0('[', round(eur_left_bound, 3), ',', round(eur_right_bound, 3), ')'),
        afr_maf_bin = paste0('[', round(afr_left_bound, 3), ',', round(afr_right_bound, 3), ')')) %>%
  ggplot(aes(x = eur_maf_bin, y = afr_maf_bin, fill = expectation_b2_normalized_meta)) + 
  geom_tile(color = 'black') + 
  theme_bw(base_size = 10) + 
  geom_text(aes(label= enrichment_rounded), color = 'black', size = 3) + 
  labs(x = 'Eur MAF Bin', y = 'Afr MAF Bin') + 
  coord_equal() + 
  scale_fill_distiller(palette = 'RdYlBu') + 
  guides(fill = guide_legend(title = 'Meta-analyzed\nper-allele\neffect size variance', position = 'right', direction = 'vertical')) + 
  theme(axis.line=element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) #+
  #geom_rect(aes(xmin = 1 - .5, xmax = 5 + .5, ymin = 1 - box_offset, ymax = 1 + box_offset), fill = 'transparent', color = "#377EB8", linewidth = box_width) + 
  #geom_rect(aes(xmin = 1 - .5, xmax = 5 + .5, ymin = 2 - box_offset, ymax = 2 + box_offset), fill = 'transparent', color = "#984EA3", linewidth = box_width) +
  #geom_rect(aes(xmin = 1 - .5, xmax = 5 + .5, ymin = 3 - box_offset, ymax = 7 + box_offset), fill = 'transparent', color = "#E41A1C", linewidth = box_width)

cowplot::plot_grid(
    grid_figure,
    lineplot,
    labels = c('a', 'b'),
    nrow = 2,
    rel_heights = c(.7, .5)
)

ggsave(
    paste0('~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/figure3', bin_type, '.pdf'), 
    width = 115,
    height = 180,
    units = 'mm')


## get some stats

maf_grid_enrichment %>%
    dplyr::filter(Category != 'base') %>%
    dplyr::mutate(enrichment_rounded = round(expectation_b2_normalized_meta, 2)) %>%
    dplyr::mutate(
        eur_left_bound = as.numeric(str_sub(word(eur_maf_bin, 1, 1, fixed(',')), 2, -1)),
        eur_right_bound = as.numeric(str_sub(word(eur_maf_bin, 2, 2, fixed(',')), 1, -2)),
        afr_left_bound = as.numeric(str_sub(word(afr_maf_bin, 1, 1, fixed(',')), 2, -1)),
        afr_right_bound = as.numeric(str_sub(word(afr_maf_bin, 2, 2, fixed(',')), 1, -2))) %>% 
    dplyr::filter(afr_maf_bin %in% c('(-0.0001,0.005]', '(0.3773,0.5]')) %>% 
    dplyr::group_by(eur_maf_bin) %>% 
    dplyr::summarise(q = expectation_b2_normalized_meta[afr_maf_bin == '(-0.0001,0.005]'] / expectation_b2_normalized_meta[afr_maf_bin == '(0.3773,0.5]'])

maf_grid_enrichment %>%
    dplyr::filter(Category != 'base') %>%
    dplyr::mutate(enrichment_rounded = round(expectation_b2_normalized_meta, 2)) %>%
    dplyr::mutate(
        eur_left_bound = as.numeric(str_sub(word(eur_maf_bin, 1, 1, fixed(',')), 2, -1)),
        eur_right_bound = as.numeric(str_sub(word(eur_maf_bin, 2, 2, fixed(',')), 1, -2)),
        afr_left_bound = as.numeric(str_sub(word(afr_maf_bin, 1, 1, fixed(',')), 2, -1)),
        afr_right_bound = as.numeric(str_sub(word(afr_maf_bin, 2, 2, fixed(',')), 1, -2))) %>%
    dplyr::filter(eur_maf_bin %in% c('(0.0499,0.09816]', '(0.3773,0.5]')) %>% 
    dplyr::group_by(afr_maf_bin) %>% 
    dplyr::summarise(q = expectation_b2_normalized_meta[eur_maf_bin == '(0.0499,0.09816]'] / expectation_b2_normalized_meta[eur_maf_bin == '(0.3773,0.5]'])
    
maf_grid_enrichment %>%
    dplyr::filter(Category != 'base') %>%
    dplyr::mutate(enrichment_rounded = round(expectation_b2_normalized_meta, 2)) %>%
    dplyr::mutate(
        eur_left_bound = as.numeric(str_sub(word(eur_maf_bin, 1, 1, fixed(',')), 2, -1)),
        eur_right_bound = as.numeric(str_sub(word(eur_maf_bin, 2, 2, fixed(',')), 1, -2)),
        afr_left_bound = as.numeric(str_sub(word(afr_maf_bin, 1, 1, fixed(',')), 2, -1)),
        afr_right_bound = as.numeric(str_sub(word(afr_maf_bin, 2, 2, fixed(',')), 1, -2))) %>%
    dplyr::filter(afr_maf_bin %in% c('(0.05,0.09816]', '(0.3773,0.5]')) %>% 
    dplyr::group_by(eur_maf_bin) %>% 
    dplyr::summarise(q = expectation_b2_normalized_meta[afr_maf_bin == '(0.05,0.09816]'] / expectation_b2_normalized_meta[afr_maf_bin == '(0.3773,0.5]'])
    
