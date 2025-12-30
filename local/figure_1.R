require(tidyverse)
require(latex2exp)

theme_Publication <- function(base_size=11, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
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
            panel.grid.major = element_line(colour="#f0f0f0"),
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


#### main fig ----
ns_stratified_by_afr_and_eur <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/nonsynonymous_eur_afr_maf_props.csv',
    data.table = F)

bin_order <- unique(ns_stratified_by_afr_and_eur$afr_maf_bin)

data_for_plotting <- ns_stratified_by_afr_and_eur %>%
    dplyr::mutate(non_synonymous = non_synonymous * 100) %>%
    dplyr::mutate(non_synonymous_rounded = str_sub(as.character(round(non_synonymous, 2)), 2, -1)) %>%
    dplyr::filter(afr_maf_bin != '[0.0, 0.002)' | eur_maf_bin != '[0.0, 0.002)') %>%
    dplyr::mutate(eur_maf_bin = factor(eur_maf_bin, levels = bin_order)) %>%
    dplyr::mutate(afr_maf_bin = factor(afr_maf_bin, levels = bin_order))

box_offset <- .5-.06
box_width <- 1.5
grid_plot <- data_for_plotting %>%
    ggplot(aes(x = eur_maf_bin, y = afr_maf_bin, fill = non_synonymous)) + 
    geom_tile(color = 'black', data = data_for_plotting %>% dplyr::filter(!is.na(non_synonymous))) + 
    theme_bw(base_size = 11) + 
    geom_text(aes(label= non_synonymous_rounded), color = 'black', size = 3) + 
    labs(x = 'Eur MAF Bin', y = 'Afr MAF Bin') + 
    coord_equal() + 
    scale_fill_distiller(palette = 'RdYlBu', na.value="white") + 
    guides(fill = guide_legend(title = 'Percent of SNPs that\n are non-synonymous', position = 'right', direction = 'vertical')) + 
    theme(axis.line=element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    geom_rect(aes(xmin = 1 - .5, xmax = 2 + box_offset, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#377EB8", linewidth = box_width) + 
    geom_rect(aes(xmin = 3 - box_offset, xmax = 3 + box_offset, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#984EA3", linewidth = box_width) +
    geom_rect(aes(xmin = 4 - box_offset, xmax = 8 + .5, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#E41A1C", linewidth = box_width) 
    #geom_rect(aes(ymin = 1 - .5, ymax = 2 + box_offset, xmin = 2 - .5, xmax = 8 + .5), fill = 'transparent', color = "#377EB8", linewidth = box_width, linetype = "dashed") + 
    #geom_rect(aes(ymin = 3 - box_offset, ymax = 3 + box_offset, xmin = 2 - .5, xmax = 8 + .5), fill = 'transparent', color = "#984EA3", linewidth = box_width, linetype = "dashed") +
    #geom_rect(aes(ymin = 4 - box_offset, ymax = 8 + .5, xmin = 2 - .5, xmax = 8 + .5), fill = 'transparent', color = "#E41A1C", linewidth = box_width, linetype = "dashed")
    #geom_rect(aes(xmin = 1 - .5, xmax = 1 + box_offset, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#377EB8", linewidth = box_width) + 
    #geom_rect(aes(xmin = 2 - .5, xmax = 2 + box_offset, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#377EB8", linewidth = box_width) + 
    #geom_rect(aes(xmin = 3 - box_offset, xmax = 3 + box_offset, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#984EA3", linewidth = box_width) +
    #geom_rect(aes(xmin = 4 - box_offset, xmax = 4 + .5, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#E41A1C", linewidth = box_width) +
    #geom_rect(aes(xmin = 5 - box_offset, xmax = 5 + .5, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#E41A1C", linewidth = box_width) +
    #geom_rect(aes(xmin = 6 - box_offset, xmax = 6 + .5, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#E41A1C", linewidth = box_width) +
    #geom_rect(aes(xmin = 7 - box_offset, xmax = 7 + .5, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#E41A1C", linewidth = box_width) +
    #geom_rect(aes(xmin = 8 - box_offset, xmax = 8 + .5, ymin = 2 - .5, ymax = 8 + .5), fill = 'transparent', color = "#E41A1C", linewidth = box_width)

bin_color_mapping <- c(
    'MAF < 0.005' = "#377EB8",
    '0.005 <= MAF < 0.05' = "#984EA3",
    'MAF >= 0.05' = "#E41A1C")

afr_logistic_predictions <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_afr_maf_logistic_predictions.csv',
    data.table = F) %>%
    dplyr::mutate(eur_maf_bin = dplyr::case_when(
        query == 'index==index' ~ 'All',
        #query == '.005 <= MAF_eur70346 < .05' ~ 'Low frequency in Europeans (0.005 <= MAF < 0.05)',
        #query == '.05 <= MAF_eur70346' ~ 'Common in Europeans (MAF >= 0.05)',
        #query == 'MAF_eur70346 < .005' ~ 'Rare in Europeans (MAF < 0.005)',
        query == '.005 <= MAF_eur70346 < .05' ~ '0.005 <= MAF < 0.05',
        query == '.05 <= MAF_eur70346' ~ 'MAF >= 0.05',
        query == 'MAF_eur70346 < .005' ~ 'MAF < 0.005',
    )) %>% dplyr::filter(
        eur_maf_bin != 'All'
    ) %>% dplyr::mutate(
        MAF = exp(log_maf),
        eur_maf_bin = factor(eur_maf_bin, names(bin_color_mapping))
    )

non_synonymous_stratified_by_eur_maf <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/prop_ns_by_afr_maf_eur_maf_stratified.csv',
        data.table = F) %>%
    dplyr::filter(eur_maf_bin != 'all') %>%
    dplyr::mutate(eur_maf_bin = dplyr::case_when(
        eur_maf_bin == '(0.05, 0.5]' ~ 'MAF >= 0.05',
        eur_maf_bin == '(0.005, 0.05]' ~ '0.005 <= MAF < 0.05',
        eur_maf_bin == '(-0.001, 0.005]' ~ 'MAF < 0.005',
        T  ~ NA)) %>% 
    dplyr::transmute(MAF = MAF_afr_unadmixed, p_non_synonymous = non_synonymous, non_synonymous_lower, non_synonymous_upper, eur_maf_bin) %>%
    dplyr::mutate(eur_maf_bin = factor(eur_maf_bin, levels = names(bin_color_mapping)))


eur_logistic_predictions <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_eur_maf_logistic_predictions.csv',
        data.table = F) %>%
    dplyr::mutate(afr_maf_bin = dplyr::case_when(
        #query == 'index==index' ~ 'All',
        #query == '.005 <= MAF_afr_unadmixed < .05' ~ 'Low frequency in Africans (0.005 <= MAF < 0.05)',
        #query == '.05 <= MAF_afr_unadmixed' ~ 'Common in Africans (MAF >= 0.05)',
        #query == 'MAF_afr_unadmixed < .005' ~ 'Rare in Africans (MAF < 0.005)',
        query == 'index==index' ~ 'All',
        query == '.005 <= MAF_afr_unadmixed < .05' ~ '0.005 <= MAF < 0.05',
        query == '.05 <= MAF_afr_unadmixed' ~ 'MAF >= 0.05',
        query == 'MAF_afr_unadmixed < .005' ~ 'MAF < 0.005',
    )) %>% dplyr::filter(
        afr_maf_bin != 'All'
    ) %>% dplyr::mutate(
        MAF = exp(log_maf),
        afr_maf_bin = factor(afr_maf_bin, names(bin_color_mapping))
    )

non_synonymous_stratified_by_afr_maf <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/prop_ns_by_eur_maf_afr_maf_stratified.csv',
        data.table = F) %>%
    dplyr::mutate(afr_maf_bin = dplyr::case_when(
        afr_maf_bin == '(0.05, 0.5]' ~ 'MAF >= 0.05',
        afr_maf_bin == '(0.005, 0.05]' ~ '0.005 <= MAF < 0.05',
        afr_maf_bin == '(-0.001, 0.005]' ~ 'MAF < 0.005',
        T  ~ NA)) %>% 
    dplyr::filter(afr_maf_bin != 'all') %>%
    dplyr::transmute(MAF = MAF_eur70346, p_non_synonymous = non_synonymous, non_synonymous_lower, non_synonymous_upper, afr_maf_bin) %>%
    dplyr::mutate(afr_maf_bin = factor(afr_maf_bin, levels = names(bin_color_mapping)))

ns_min <- min(c(afr_logistic_predictions$p_non_synonymous, eur_logistic_predictions$p_non_synonymous, non_synonymous_stratified_by_eur_maf$non_synonymous_lower, non_synonymous_stratified_by_afr_maf$non_synonymous_lower))
ns_max <- max(c(afr_logistic_predictions$p_non_synonymous, eur_logistic_predictions$p_non_synonymous, non_synonymous_stratified_by_eur_maf$non_synonymous_upper, non_synonymous_stratified_by_afr_maf$non_synonymous_upper))



afr_logistic_plot <- afr_logistic_predictions %>%
    ggplot(aes(x = MAF, y = p_non_synonymous, color = eur_maf_bin, fill = eur_maf_bin)) + 
    geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), alpha = .25, color = NA) +
    #geom_pointrange(aes(ymin = non_synonymous_lower, ymax = non_synonymous_upper), data = non_synonymous_stratified_by_eur_maf, size = .5, alpha = .25) + 
    geom_pointrange(aes(ymin = p_non_synonymous, ymax = p_non_synonymous), data = non_synonymous_stratified_by_eur_maf, size = .5, alpha = .25) + 
    geom_line() + 
    annotate('text', label = TeX('$\\gamma=-0.13$', output = 'character'), x = .002, y = .0029 * 1.13 * 1.13, parse = TRUE, color = "#377EB8", size = 4, hjust = 0) +
    annotate('text', label = TeX('$\\gamma=-0.13$', output = 'character'), x = .002, y = .0029 * 1.13, parse = TRUE, color = "#984EA3", size = 4, hjust = 0) +
    annotate('text', label = TeX('$\\gamma=-0.1\\0$', output = 'character'), x = .002, y = .0029, parse = TRUE, color = "#E41A1C", size = 4, hjust = 0) +
    scale_x_log10() + 
    scale_color_manual(values = bin_color_mapping) +
    scale_fill_manual(values = bin_color_mapping) +
    guides(color = guide_legend(title = 'European\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE)), 
           fill = guide_legend(title = 'European\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE))) + 
    labs(x = 'African MAF', y = 'Proportion of SNPs that\nare non-synonymous') +
    theme_Publication() + 
    scale_y_log10(limits = c(ns_min, ns_max)) 

eur_logistic_plot <- eur_logistic_predictions %>%
    dplyr::group_by(afr_maf_bin) %>%
    dplyr::filter(row_number() %in% c(1, n(), sample(2:n()-1, 1000))) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x = MAF, y = p_non_synonymous, color = afr_maf_bin, fill = afr_maf_bin)) + 
    geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), alpha = .25, color = NA) +
    #geom_pointrange(aes(ymin = non_synonymous_lower, ymax = non_synonymous_upper), data = non_synonymous_stratified_by_afr_maf, size = .5, alpha = .25) + 
    geom_pointrange(aes(ymin = p_non_synonymous, ymax = p_non_synonymous), data = non_synonymous_stratified_by_afr_maf, size = .5, alpha = .25) + 
    geom_line(linetype = 'dashed') + 
    annotate('text', label = TeX('$\\gamma=-0.1\\0$', output = 'character'), x = .002, y = .0057, parse = TRUE, color = "#377EB8", size = 4, hjust = 0) +
    annotate('text', label = TeX('$\\gamma=-0.047$', output = 'character'), x = .002, y = .0039, parse = TRUE, color = "#984EA3", size = 4, hjust = 0) +
    annotate('text', label = TeX('$\\gamma=-0.028$', output = 'character'), x = .002, y = .0029, parse = TRUE, color = "#E41A1C", size = 4, hjust = 0) +
    scale_x_log10() + 
    guides(color = guide_legend(title = 'African\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE)),
           fill = guide_legend(title = 'African\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE))) + 
    scale_color_manual(values = bin_color_mapping) +
    scale_fill_manual(values = bin_color_mapping) +
    labs(x = 'European MAF', y = 'Proportion of SNPs that\nare non-synonymous') +
    theme_Publication() + 
    scale_y_log10(limits = c(ns_min, ns_max))


row2 <- cowplot::plot_grid(
    afr_logistic_plot + theme(legend.position = 'none'),
    eur_logistic_plot + theme(legend.position = 'none'),
    nrow = 1, 
    labels = c('b', 'c')
)
row3 <- cowplot::plot_grid(
    cowplot::get_plot_component(afr_logistic_plot, 'guide-box-bottom', return_all = TRUE),
    cowplot::get_plot_component(eur_logistic_plot, 'guide-box-bottom', return_all = TRUE)
)

combined <- cowplot::plot_grid(
    grid_plot,
    row2,
    row3,
    nrow = 3,
    labels = c('a'),
    rel_heights = c(.7, .43, .07)
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/figure_1.pdf',
    combined,
    width = 180,
    height = 180,
    unit = 'mm'
)

#### SUPPLEMENTARY GRID FIGURES WITH N AND SE, THESE WILL NOT RUN BY DEFAULT!!!! ----
grid_se_plot <- data_for_plotting %>%
    dplyr::mutate(non_synonymous_se = non_synonymous_se * 100) %>%
    dplyr::mutate(se_rounded = str_sub(as.character(round(non_synonymous_se, 3)), 2, -1)) %>%
    ggplot(aes(x = eur_maf_bin, y = afr_maf_bin, fill = non_synonymous_se)) + 
    geom_tile(color = 'black') + 
    theme_bw(base_size = 11) + 
    geom_text(aes(label= se_rounded), color = 'black', size = 3.5) + 
    labs(x = 'Eur MAF Bin', y = 'Afr MAF Bin') + 
    coord_equal() + 
    scale_fill_distiller(palette = 'RdYlBu', na.value="white") + 
    guides(fill = guide_legend(title = 'Standard error on \npercent of variants \nthat are non-synonymous', position = 'right', direction = 'vertical')) + 
    theme(axis.line=element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(
    '~/google drive/my drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/se_grid.pdf',
    grid_se_plot,
    width = 180,
    height = 140,
    unit = 'mm'
)

grid_n_plot <- data_for_plotting %>%
    dplyr::mutate(n_millions = n / 1e6) %>%
    dplyr::mutate(n_rounded = round(n_millions, 2)) %>%
    ggplot(aes(x = eur_maf_bin, y = afr_maf_bin, fill = n_millions)) + 
    geom_tile(color = 'black') + 
    theme_bw(base_size = 11) + 
    geom_text(aes(label= n_rounded), color = 'black', size = 3) + 
    labs(x = 'Eur MAF Bin', y = 'Afr MAF Bin') + 
    coord_equal() + 
    scale_fill_distiller(palette = 'RdYlBu', na.value="white") + 
    guides(fill = guide_legend(title = 'Number of SNPs \nin bin (millions)', position = 'right', direction = 'vertical')) + 
    theme(axis.line=element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(
    '~/google drive/my drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/n_grid.pdf',
    grid_n_plot,
    width = 180,
    height = 140,
    unit = 'mm'
)

#### SUPPLEMENTARY FIGURE FOR LOGISTIC REGRESSION MATCHES DATA ----
afr_logistic_plot <- afr_logistic_predictions %>%
    dplyr::mutate(eur_maf_bin = factor(eur_maf_bin, levels = names(bin_color_mapping))) %>%
    dplyr::arrange(eur_maf_bin) %>%
    ggplot(aes(x = MAF, y = p_non_synonymous, color = eur_maf_bin, fill = eur_maf_bin)) + 
    geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), alpha = .25, color = NA) +
    geom_line() + 
    geom_pointrange(aes(ymin = non_synonymous_lower, ymax = non_synonymous_upper), data = non_synonymous_stratified_by_eur_maf, size = .5) + 
    scale_x_log10() + 
    scale_color_manual(values = bin_color_mapping) +
    scale_fill_manual(values = bin_color_mapping) +
    guides(color = guide_legend(title = 'European\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE)), 
           fill = guide_legend(title = 'European\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE))) + 
    labs(x = 'African MAF', y = 'P(Non-synonymous)') +
    theme_Publication() + 
    facet_wrap('eur_maf_bin', nrow = 3) + 
    scale_y_log10(limits = c(ns_min, ns_max))
    
eur_logistic_plot <- eur_logistic_predictions %>%
    dplyr::group_by(afr_maf_bin) %>%
    dplyr::filter(row_number() %in% c(1, n(), sample(2:n()-1, 1000))) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x = MAF, y = p_non_synonymous, color = afr_maf_bin, fill = afr_maf_bin)) + 
    geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), alpha = .25, color = NA) +
    geom_line(linetype = 'dashed') + 
    geom_pointrange(aes(ymin = non_synonymous_lower, ymax = non_synonymous_upper), data = non_synonymous_stratified_by_afr_maf, size = .5) + 
    scale_x_log10() + 
    guides(color = guide_legend(title = 'African\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE)),
           fill = guide_legend(title = 'African\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE))) + 
    scale_color_manual(values = bin_color_mapping) +
    scale_fill_manual(values = bin_color_mapping) +
    labs(x = 'European MAF', y = 'P(Non-synonymous)') +
    theme_Publication() + 
    facet_wrap('afr_maf_bin', nrow = 3) + 
    scale_y_log10(limits = c(ns_min, ns_max))

combined_logistic_agreement_plot <- cowplot::plot_grid(
    afr_logistic_plot,
    eur_logistic_plot,
    labels = c('a', 'b'),
    nrow = 1
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/logistic_agreement_plot.pdf',
    combined_logistic_agreement_plot,
    width = 180,
    height = 180,
    units = 'mm'
)


#### SUPPLEMENTARY FIGURE WITH MARGINALS AND UNSTRATIFIED SNPS ----
ns_stratified_by_afr <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/nonsynonymous_afr_maf_props.csv',
    data.table = FALSE) %>%
    dplyr::filter(afr_maf_bin != '[0.0, 0.002)') %>%
    dplyr::select(-V1) %>%
    dplyr::mutate(
        eur_maf_bin = 'All'
    )

ns_stratified_by_eur <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/nonsynonymous_eur_maf_props.csv',
    data.table = FALSE) %>%
    dplyr::filter(eur_maf_bin != '[0.0, 0.002)') %>%
    dplyr::select(-V1) %>%
    dplyr::mutate(
        afr_maf_bin = 'All'
    )

spacer_column <- ns_stratified_by_afr %>%
    dplyr::mutate(eur_maf_bin = '', non_synonymous = NA)


spacer_row <- ns_stratified_by_eur %>%
    dplyr::mutate(afr_maf_bin = '', non_synonymous = NA)

data_for_plotting <- ns_stratified_by_afr_and_eur %>%
    dplyr::select(-n, -non_synonymous_se, -level_2) %>%
    rbind(ns_stratified_by_afr) %>%
    rbind(spacer_column) %>%
    rbind(ns_stratified_by_eur) %>%
    rbind(spacer_row) %>%
    dplyr::mutate(non_synonymous = non_synonymous * 100) %>%
    dplyr::mutate(non_synonymous_rounded = str_sub(as.character(round(non_synonymous, 2)), 2, -1)) %>%
    dplyr::filter(afr_maf_bin != '[0.0, 0.002)' | eur_maf_bin != '[0.0, 0.002)') %>%
    dplyr::mutate(eur_maf_bin = factor(eur_maf_bin, levels = c('[0.0, 0.002)', ns_stratified_by_afr$afr_maf_bin, '', 'All'))) %>%
    dplyr::mutate(afr_maf_bin = factor(afr_maf_bin, levels = c('[0.0, 0.002)', ns_stratified_by_afr$afr_maf_bin, '', 'All')))

grid_plot <- data_for_plotting %>%
    ggplot(aes(x = eur_maf_bin, y = afr_maf_bin, fill = non_synonymous)) + 
    geom_tile() + 
    geom_tile(color = 'black', data = data_for_plotting %>% dplyr::filter(!is.na(non_synonymous))) + 
    theme_bw(base_size = 10) + 
    geom_text(aes(label= non_synonymous_rounded), color = 'black', size = 3) + 
    labs(x = 'Eur MAF Bin', y = 'Afr MAF Bin') + 
    coord_equal() + 
    scale_fill_distiller(palette = 'RdYlBu', na.value="white") + 
    guides(fill = guide_legend(title = 'Percent of variants \nthat are non-synonymous', position = 'right', direction = 'vertical')) + 
    theme(axis.line=element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

non_synonymous_stratified_by_eur_maf <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/prop_ns_by_afr_maf_eur_maf_stratified.csv',
        data.table = F) %>%
    dplyr::filter(eur_maf_bin != 'all') %>%
    dplyr::mutate(eur_maf_bin = dplyr::case_when(
        eur_maf_bin == '(0.05, 0.5]' ~ 'Common in Europeans (0.05 <= MAF <= 0.5)',
        eur_maf_bin == '(0.005, 0.05]' ~ 'Low frequency in Europeans (.005 < MAF <= 0.05)',
        eur_maf_bin == '(-0.001, 0.005]' ~ 'Rare in Europeans (MAF <= 0.005)',
        eur_maf_bin == '(0.005, 0.5]' ~ 'Common/Low frequency in Europeans (0.005 < MAF <= 0.5)',
    ))
non_synonymous_stratified_by_afr_maf <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/prop_ns_by_eur_maf_afr_maf_stratified.csv',
        data.table = F) %>%
    dplyr::filter(afr_maf_bin != 'all') %>%
    dplyr::mutate(afr_maf_bin = dplyr::case_when(
        afr_maf_bin == '(0.05, 0.5]' ~ 'Common in Africans (0.05 <= MAF <= 0.5)',
        afr_maf_bin == '(0.005, 0.05]' ~ 'Low frequency in Africans (.005 < MAF <= 0.05)',
        afr_maf_bin == '(-0.001, 0.005]' ~ 'Rare in Africans (MAF <= 0.005)',
        afr_maf_bin == '(0.005, 0.5]' ~ 'Common/Low frequency in Africans (0.005 < MAF <= 0.5)',
    ))

max_non_synonymous <- max(c(non_synonymous_stratified_by_eur_maf$non_synonymous_upper, non_synonymous_stratified_by_afr_maf$non_synonymous_upper))
min_non_synonymous <- min(c(non_synonymous_stratified_by_eur_maf$non_synonymous_lower, non_synonymous_stratified_by_afr_maf$non_synonymous_lower)) 
max_maf <- max(c(non_synonymous_stratified_by_eur_maf$MAF_afr_unadmixed, non_synonymous_stratified_by_afr_maf$MAF_eur70346))
min_maf <- min(c(non_synonymous_stratified_by_eur_maf$MAF_afr_unadmixed, non_synonymous_stratified_by_afr_maf$MAF_eur70346))

eur_strat_line_plot <- non_synonymous_stratified_by_eur_maf %>% 
    ggplot(aes(x = MAF_afr_unadmixed, color = eur_maf_bin, y = non_synonymous, ymax = non_synonymous_upper, ymin = non_synonymous_lower)) + 
    geom_pointrange() + 
    scale_x_log10(limits = c(min_maf, max_maf)) + 
    #scale_y_log10(limits = c(min_non_synonymous, max_non_synonymous)) + 
    geom_line() + 
    scale_color_manual(values = c("#377EB8","#984EA3", "#E41A1C")) + 
    theme_Publication() + 
    guides(color = guide_legend(title = '', nrow = 3)) + 
    labs(x = 'African MAF', y = 'Proportion of variants \nthat are non-synonymous')
    

afr_strat_line_plot <- non_synonymous_stratified_by_afr_maf %>% 
    ggplot(aes(x = MAF_eur70346, color = afr_maf_bin, y = non_synonymous, ymax = non_synonymous_upper, ymin = non_synonymous_lower)) + 
    geom_pointrange() + 
    scale_x_log10(limits = c(min_maf, max_maf)) + 
    #scale_y_log10(limits = c(min_non_synonymous, max_non_synonymous)) + 
    geom_line() + 
    scale_color_manual(values = c("#377EB8","#984EA3", "#E41A1C")) + 
    theme_Publication() + 
    guides(color = guide_legend(title = '', nrow = 3)) + 
    labs(x = 'European MAF', y = 'Proportion of variants \nthat are non-synonymous')
    
#row2 <- cowplot::plot_grid(
#    eur_strat_line_plot,
#    afr_strat_line_plot,
#    nrow = 1, 
#    labels = c('b', 'c')
#)
#
#combined <- cowplot::plot_grid(
#    grid_plot,
#    row2,
#    nrow = 2,
#    labels = c('a')
#)
#
#ggsave(
#    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/non_synonymous_figure.pdf',
#    combined,
#    width = 180,
#    height = 170,
#    unit = 'mm'
#)

bin_color_mapping <- c(
    'MAF < 0.005' = "#377EB8",
    '0.005 <= MAF < 0.05' = "#984EA3",
    'MAF >= 0.05' = "#E41A1C",
    'All' = 'black')

afr_logistic_predictions <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_afr_maf_logistic_predictions.csv',
    data.table = F) %>%
    dplyr::mutate(eur_maf_bin = dplyr::case_when(
        query == 'index==index' ~ 'All',
        #query == '.005 <= MAF_eur70346 < .05' ~ 'Low frequency in Europeans (0.005 <= MAF < 0.05)',
        #query == '.05 <= MAF_eur70346' ~ 'Common in Europeans (MAF >= 0.05)',
        #query == 'MAF_eur70346 < .005' ~ 'Rare in Europeans (MAF < 0.005)',
        query == '.005 <= MAF_eur70346 < .05' ~ '0.005 <= MAF < 0.05',
        query == '.05 <= MAF_eur70346' ~ 'MAF >= 0.05',
        query == 'MAF_eur70346 < .005' ~ 'MAF < 0.005',
    )) %>% dplyr::mutate(
        MAF = exp(log_maf),
        eur_maf_bin = factor(eur_maf_bin, names(bin_color_mapping))
    )


eur_logistic_predictions <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_eur_maf_logistic_predictions.csv',
        data.table = F) %>%
    dplyr::mutate(afr_maf_bin = dplyr::case_when(
        #query == 'index==index' ~ 'All',
        #query == '.005 <= MAF_afr_unadmixed < .05' ~ 'Low frequency in Africans (0.005 <= MAF < 0.05)',
        #query == '.05 <= MAF_afr_unadmixed' ~ 'Common in Africans (MAF >= 0.05)',
        #query == 'MAF_afr_unadmixed < .005' ~ 'Rare in Africans (MAF < 0.005)',
        query == 'index==index' ~ 'All',
        query == '.005 <= MAF_afr_unadmixed < .05' ~ '0.005 <= MAF < 0.05',
        query == '.05 <= MAF_afr_unadmixed' ~ 'MAF >= 0.05',
        query == 'MAF_afr_unadmixed < .005' ~ 'MAF < 0.005',
    )) %>% dplyr::mutate(
        MAF = exp(log_maf),
        afr_maf_bin = factor(afr_maf_bin, names(bin_color_mapping))
    )

ns_min <- min(c(afr_logistic_predictions$p_non_synonymous, eur_logistic_predictions$p_non_synonymous))
ns_max <- max(c(afr_logistic_predictions$p_non_synonymous, eur_logistic_predictions$p_non_synonymous))

afr_logistic_plot <- afr_logistic_predictions %>%
    ggplot(aes(x = MAF, y = p_non_synonymous, color = eur_maf_bin)) + 
    geom_line() + 
    geom_line(data = afr_logistic_predictions %>% dplyr::filter(eur_maf_bin == 'All'), linewidth = 2) + 
    scale_x_log10() + 
    guides(color = guide_legend(title = 'European\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE))) + 
    scale_color_manual(values = bin_color_mapping) +
    labs(x = 'African MAF', y = 'P(Non-synonymous)') +
    theme_Publication() + 
    scale_y_log10(limits = c(ns_min, ns_max))#$ + 
    #lims(y = c(ns_min, ns_max))

eur_logistic_plot <- eur_logistic_predictions %>%
    dplyr::group_by(afr_maf_bin) %>%
    dplyr::filter(row_number() %in% c(1, n(), sample(2:n()-1, 500))) %>%
    dplyr::ungroup() %>%
    {ggplot(., aes(x = MAF, y = p_non_synonymous, color = afr_maf_bin)) + 
    geom_line(data = . %>% dplyr::filter(afr_maf_bin != 'All'), linetype = 'dashed') + 
    geom_line(data = . %>% dplyr::filter(afr_maf_bin == 'All'), linewidth = 2, linetype = 'dashed') + 
    scale_x_log10() + 
    guides(color = guide_legend(title = 'African\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE))) + 
    scale_color_manual(values = bin_color_mapping) +
    labs(x = 'European MAF', y = 'P(Non-synonymous)') +
    theme_Publication() + 
    scale_y_log10(limits = c(ns_min, ns_max))}

row2 <- cowplot::plot_grid(
    afr_logistic_plot + theme(legend.position = 'none'),
    eur_logistic_plot + theme(legend.position = 'none'),
    nrow = 1, 
    labels = c('b', 'c')
)
row3 <- cowplot::plot_grid(
    cowplot::get_plot_component(afr_logistic_plot, 'guide-box-bottom', return_all = TRUE),
    cowplot::get_plot_component(eur_logistic_plot, 'guide-box-bottom', return_all = TRUE)
)

combined <- cowplot::plot_grid(
    grid_plot,
    row2,
    row3,
    nrow = 3,
    labels = c('a'),
    rel_heights = c(.7, .43, .07)
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_1_with_marginals.pdf',
    combined,
    width = 180,
    height = 180,
    unit = 'mm'
)

#### Supplementary Figure restricting all analyses to SNPs with P_a > 0.002 and p_e > 0.002 ----
ns_stratified_by_afr_and_eur <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/nonsynonymous_eur_afr_maf_props.csv',
    data.table = F)


data_for_plotting <- ns_stratified_by_afr_and_eur %>%
    dplyr::mutate(non_synonymous = non_synonymous * 100) %>%
    dplyr::mutate(non_synonymous_rounded = str_sub(as.character(round(non_synonymous, 2)), 2, -1)) %>%
    dplyr::filter(afr_maf_bin != '[0.0, 0.002)' | eur_maf_bin != '[0.0, 0.002)') %>%
    dplyr::filter(as.numeric(str_sub(word(afr_maf_bin, 1, 1, fixed(',')), 2, -1)) >= 0.002) %>%
    dplyr::filter(as.numeric(str_sub(word(eur_maf_bin, 1, 1, fixed(',')), 2, -1)) >= 0.002)
bin_order <- unique(ns_stratified_by_afr_and_eur$afr_maf_bin)
data_for_plotting <- data_for_plotting %>%
    dplyr::mutate(eur_maf_bin = factor(eur_maf_bin, levels = bin_order)) %>%
    dplyr::mutate(afr_maf_bin = factor(afr_maf_bin, levels = bin_order))

box_offset <- .5-.06
box_width <- 1.5
grid_plot <- data_for_plotting %>%
    ggplot(aes(x = eur_maf_bin, y = afr_maf_bin, fill = non_synonymous)) + 
    geom_tile(color = 'black', data = data_for_plotting %>% dplyr::filter(!is.na(non_synonymous))) + 
    theme_bw(base_size = 11) + 
    geom_text(aes(label= non_synonymous_rounded), color = 'black', size = 3) + 
    labs(x = 'Eur MAF Bin', y = 'Afr MAF Bin') + 
    coord_equal() + 
    scale_fill_distiller(palette = 'RdYlBu', na.value="white") + 
    guides(fill = guide_legend(title = 'Percent of SNPs that\n are non-synonymous', position = 'right', direction = 'vertical')) + 
    theme(axis.line=element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    geom_rect(aes(xmin = 1 - .5, xmax = 1 + box_offset, ymin = 1 - .5, ymax = 7 + .5), fill = 'transparent', color = "#377EB8", linewidth = box_width) + 
    geom_rect(aes(xmin = 2 - box_offset, xmax = 2 + box_offset, ymin = 1 - .5, ymax = 7 + .5), fill = 'transparent', color = "#984EA3", linewidth = box_width) +
    geom_rect(aes(xmin = 3 - box_offset, xmax = 7 + .5, ymin = 1 - .5, ymax = 7 + .5), fill = 'transparent', color = "#E41A1C", linewidth = box_width)

afr_logistic_predictions <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_afr_maf_logistic_predictions_afr_eur_ge_.002.csv',
    data.table = F) %>%
    dplyr::mutate(eur_maf_bin = dplyr::case_when(
        query == 'index==index' ~ 'All',
        #query == '.005 <= MAF_eur70346 < .05' ~ 'Low frequency in Europeans (0.005 <= MAF < 0.05)',
        #query == '.05 <= MAF_eur70346' ~ 'Common in Europeans (MAF >= 0.05)',
        #query == 'MAF_eur70346 < .005' ~ 'Rare in Europeans (MAF < 0.005)',
        query == '.005 <= MAF_eur70346 < .05' ~ '0.005 <= MAF < 0.05',
        query == '.05 <= MAF_eur70346' ~ 'MAF >= 0.05',
        query == 'MAF_eur70346 < .005' ~ 'MAF < 0.005',
    )) %>% dplyr::filter(
        eur_maf_bin != 'All'
    ) %>% dplyr::mutate(
        MAF = exp(log_maf),
        eur_maf_bin = factor(eur_maf_bin, names(bin_color_mapping))
    )

non_synonymous_stratified_by_eur_maf <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/prop_ns_by_afr_maf_eur_maf_stratified_afr_eur_ge_.002.csv',
        data.table = F) %>%
    dplyr::filter(eur_maf_bin != 'all') %>%
    dplyr::mutate(eur_maf_bin = dplyr::case_when(
        eur_maf_bin == '(0.05, 0.5]' ~ 'MAF >= 0.05',
        eur_maf_bin == '(0.005, 0.05]' ~ '0.005 <= MAF < 0.05',
        eur_maf_bin == '(-0.001, 0.005]' ~ 'MAF < 0.005',
        T  ~ NA)) %>% 
    dplyr::transmute(MAF = MAF_afr_unadmixed, p_non_synonymous = non_synonymous, non_synonymous_lower, non_synonymous_upper, eur_maf_bin) %>%
    dplyr::mutate(eur_maf_bin = factor(eur_maf_bin, levels = names(bin_color_mapping)))


eur_logistic_predictions <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_eur_maf_logistic_predictions_afr_eur_ge_.002.csv',
        data.table = F) %>%
    dplyr::mutate(afr_maf_bin = dplyr::case_when(
        #query == 'index==index' ~ 'All',
        #query == '.005 <= MAF_afr_unadmixed < .05' ~ 'Low frequency in Africans (0.005 <= MAF < 0.05)',
        #query == '.05 <= MAF_afr_unadmixed' ~ 'Common in Africans (MAF >= 0.05)',
        #query == 'MAF_afr_unadmixed < .005' ~ 'Rare in Africans (MAF < 0.005)',
        query == 'index==index' ~ 'All',
        query == '.005 <= MAF_afr_unadmixed < .05' ~ '0.005 <= MAF < 0.05',
        query == '.05 <= MAF_afr_unadmixed' ~ 'MAF >= 0.05',
        query == 'MAF_afr_unadmixed < .005' ~ 'MAF < 0.005',
    )) %>% dplyr::filter(
        afr_maf_bin != 'All'
    ) %>% dplyr::mutate(
        MAF = exp(log_maf),
        afr_maf_bin = factor(afr_maf_bin, names(bin_color_mapping))
    )

non_synonymous_stratified_by_afr_maf <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/prop_ns_by_eur_maf_afr_maf_stratified_afr_eur_ge_.002.csv',
        data.table = F) %>%
    dplyr::mutate(afr_maf_bin = dplyr::case_when(
        afr_maf_bin == '(0.05, 0.5]' ~ 'MAF >= 0.05',
        afr_maf_bin == '(0.005, 0.05]' ~ '0.005 <= MAF < 0.05',
        afr_maf_bin == '(-0.001, 0.005]' ~ 'MAF < 0.005',
        T  ~ NA)) %>% 
    dplyr::filter(afr_maf_bin != 'all') %>%
    dplyr::transmute(MAF = MAF_eur70346, p_non_synonymous = non_synonymous, non_synonymous_lower, non_synonymous_upper, afr_maf_bin) %>%
    dplyr::mutate(afr_maf_bin = factor(afr_maf_bin, levels = names(bin_color_mapping)))

ns_min <- min(c(afr_logistic_predictions$p_non_synonymous, eur_logistic_predictions$p_non_synonymous, non_synonymous_stratified_by_eur_maf$non_synonymous_lower, non_synonymous_stratified_by_afr_maf$non_synonymous_lower))
ns_max <- max(c(afr_logistic_predictions$p_non_synonymous, eur_logistic_predictions$p_non_synonymous, non_synonymous_stratified_by_eur_maf$non_synonymous_upper, non_synonymous_stratified_by_afr_maf$non_synonymous_upper))

bin_color_mapping <- c(
    'MAF < 0.005' = "#377EB8",
    '0.005 <= MAF < 0.05' = "#984EA3",
    'MAF >= 0.05' = "#E41A1C")

afr_logistic_plot <- afr_logistic_predictions %>%
    ggplot(aes(x = MAF, y = p_non_synonymous, color = eur_maf_bin, fill = eur_maf_bin)) + 
    geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), alpha = .25, color = NA) +
    geom_pointrange(aes(ymin = non_synonymous_lower, ymax = non_synonymous_upper), data = non_synonymous_stratified_by_eur_maf, size = .5, alpha = .25) + 
    #geom_pointrange(aes(ymin = p_non_synonymous, ymax = p_non_synonymous), data = non_synonymous_stratified_by_eur_maf, size = .5, alpha = .25) + 
    geom_line() + 
    scale_x_log10() + 
    scale_color_manual(values = bin_color_mapping) +
    scale_fill_manual(values = bin_color_mapping) +
    guides(color = guide_legend(title = 'European\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE)), 
           fill = guide_legend(title = 'European\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE))) + 
    labs(x = 'African MAF', y = 'Proportion of SNPs that\nare non-synonymous') +
    theme_Publication() + 
    scale_y_log10(limits = c(ns_min, ns_max)) 

eur_logistic_plot <- eur_logistic_predictions %>%
    dplyr::group_by(afr_maf_bin) %>%
    dplyr::filter(row_number() %in% c(1, n(), sample(2:n()-1, 1000))) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x = MAF, y = p_non_synonymous, color = afr_maf_bin, fill = afr_maf_bin)) + 
    geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), alpha = .25, color = NA) +
    geom_pointrange(aes(ymin = non_synonymous_lower, ymax = non_synonymous_upper), data = non_synonymous_stratified_by_afr_maf, size = .5, alpha = .25) + 
    #geom_pointrange(aes(ymin = p_non_synonymous, ymax = p_non_synonymous), data = non_synonymous_stratified_by_afr_maf, size = .5, alpha = .25) + 
    geom_line(linetype = 'dashed') + 
    scale_x_log10() + 
    guides(color = guide_legend(title = 'African\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE)),
           fill = guide_legend(title = 'African\nMAF', nrow = 2, theme = theme(legend.byrow = TRUE))) + 
    scale_color_manual(values = bin_color_mapping) +
    scale_fill_manual(values = bin_color_mapping) +
    labs(x = 'European MAF', y = 'Proportion of SNPs that\nare non-synonymous') +
    theme_Publication() + 
    scale_y_log10(limits = c(ns_min, ns_max))

row2 <- cowplot::plot_grid(
    afr_logistic_plot + theme(legend.position = 'none'),
    eur_logistic_plot + theme(legend.position = 'none'),
    nrow = 1, 
    labels = c('b', 'c')
)
row3 <- cowplot::plot_grid(
    cowplot::get_plot_component(afr_logistic_plot, 'guide-box-bottom', return_all = TRUE),
    cowplot::get_plot_component(eur_logistic_plot, 'guide-box-bottom', return_all = TRUE)
)

combined <- cowplot::plot_grid(
    grid_plot,
    row2,
    row3,
    nrow = 3,
    labels = c('a'),
    rel_heights = c(.67, .53, .12)
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/fig_1_.002_both.pdf',
    combined,
    width = 180,
    height = 160,
    units = 'mm'
)

#### marginal logistic models varying SNP sets ----
pmix_either <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_p_mix_params_ge.002_in_either_afr_eur.csv', 
    data.table = F) %>% 
    #dplyr::mutate(variant_set = 'p[A]>=0.002~\"and/or\"~p[E]>=0.002')
    dplyr::mutate(variant_set = '$p_A \\geq 0.002$ and/or $p_E \\geq 0.002$')
pmix_common_both <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_p_mix_params_ge.05_afr_eur.csv', 
    data.table = F) %>% 
    #dplyr::mutate(variant_set = 'p[A]>=0.05~\"and\"~p[E]>=0.05')
    dplyr::mutate(variant_set = '$p_A \\geq 0.05$ and $p_E \\geq 0.05$')
pmix_both <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_p_mix_params_ge.002_afr_eur.csv', 
    data.table = F) %>% 
    #dplyr::mutate(variant_set = 'p[A]>=0.002~\"and\"~p[E]>=0.002')
    dplyr::mutate(variant_set = '$p_A \\geq 0.002$ and $p_E \\geq 0.002$')
pmix_afr <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_p_mix_params_ge.002_in_afr.csv', 
    data.table = F) %>% 
    #dplyr::mutate(variant_set = 'p[A]>=0.002')
    dplyr::mutate(variant_set = '$p_A \\geq 0.002$')
pmix_eur <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_p_mix_params_ge.002_in_eur.csv', 
    data.table = F) %>% 
    #dplyr::mutate(variant_set = 'p[E]>=0.002')
    dplyr::mutate(variant_set = '$p_E \\geq 0.002$')

combined_pmix <- rbind(pmix_either, pmix_both, pmix_afr, pmix_eur, pmix_common_both) %>%
    dplyr::filter(w %in% c(0, 1)) %>%
    dplyr::mutate(ancestry = dplyr::case_when(
        w == 0 ~ '$p_E$', 
        w == 1 ~ '$p_A$'))

delta_log_likelihood_plot <- combined_pmix %>% 
    dplyr::mutate(variant_set = factor(variant_set, levels = rev(unique(variant_set)))) %>%
    dplyr::group_by(variant_set) %>%
    dplyr::summarise(delta_log_likelihood = c(l[ancestry == '$p_A$'] - l[ancestry == '$p_E$'])[1]) %>%
    ggplot(aes(x = variant_set, y = delta_log_likelihood)) + 
    geom_col() +
    coord_flip() +
    scale_y_log10() +
    scale_x_discrete(labels = TeX) +
    theme_Publication() +
    labs(y = 'African MAF log-likelihood - European MAF log-likelihood (log-scale)', x = 'SNP set')

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/marginal_regression_shared_snp_sets.pdf',
    delta_log_likelihood_plot,
    width = 180,
    height = 100,
    units = 'mm'
)

#loglikelihood_plot <- combined_pmix %>%
#    ggplot(aes(x = ancestry, y = l)) + 
#    geom_point() + 
#    theme_Publication() + 
#    labs(x = 'Predictor', y = 'log-likelihood') +
#    facet_wrap(variant_set ~ ., scales = 'free', labeller = label_parsed, nrow = 5) +
#    scale_x_discrete(labels = TeX) 
#
#
#beta_plot <- combined_pmix %>%
#    dplyr::mutate(beta_upper = beta + 1.96 * beta_se, beta_lower = beta - 1.96 * beta_se) %>%
#    ggplot(aes(x = ancestry, y = beta, ymin = beta_lower, ymax = beta_upper)) + 
#    geom_point() +
#    geom_linerange(color = 'red') +
#    #geom_col() + 
#    theme_Publication() + 
#    labs(x = 'Predictor', y = TeX('$\\gamma$ (regression coefficient)')) +
#    facet_wrap(variant_set ~ . , labeller = label_parsed, nrow = 5) + 
#    #geom_errorbar(color = 'red', width = .5) +
#    scale_x_discrete(labels = TeX) 

#combined <- cowplot::plot_grid(loglikelihood_plot, beta_plot, labels = c('a', 'b'), nrow  = 1)

#ggsave(
#    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/marginal_regression_shared_snp_set.pdf',
#    combined,
#    width = 140,
#    height = 160,
#    units = 'mm'
#)

#### get some stats ----

ns_stratified_by_afr_and_eur <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/nonsynonymous_eur_afr_maf_props.csv',
    data.table = F)

bin_order <- unique(ns_stratified_by_afr_and_eur$afr_maf_bin)

data_for_plotting <- ns_stratified_by_afr_and_eur %>%
    dplyr::mutate(non_synonymous = non_synonymous * 100) %>%
    dplyr::mutate(non_synonymous_rounded = str_sub(as.character(round(non_synonymous, 2)), 2, -1)) %>%
    dplyr::filter(afr_maf_bin != '[0.0, 0.002)' | eur_maf_bin != '[0.0, 0.002)') %>%
    dplyr::mutate(eur_maf_bin = factor(eur_maf_bin, levels = bin_order)) %>%
    dplyr::mutate(afr_maf_bin = factor(afr_maf_bin, levels = bin_order))

data_for_plotting %>% 
    dplyr::group_by(eur_maf_bin) %>%
    dplyr::summarise(non_synonymous[afr_maf_bin == '[0.05, 0.1)'] / non_synonymous[afr_maf_bin == '[0.38, 0.5)'])

data_for_plotting %>% 
    dplyr::group_by(afr_maf_bin) %>%
    dplyr::summarise(non_synonymous[afr_maf_bin == '[0.05, 0.1)'] / non_synonymous[afr_maf_bin == '[0.38, 0.5)'])
