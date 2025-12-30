require(tidyverse)
require(cowplot)
require(latex2exp)

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

jordans_list <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/jordans_trait_list.txt', sep = '\t', data.table = F) %>%
    dplyr::rename(t = trait_identifier)

trait_categories <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/baselineLD2.2_plus_maf_grid_afr_maf_bins/nh2_selected_traits/thresholded/name_mapping.csv',
        data.table = F) %>% 
    dplyr::rename(t = V1) %>% 
    dplyr::left_join(jordans_list)
f <- '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/alpha_mix_trait_specific.csv'
alpha_mix_df <- data.table::fread(f, data.table = F) %>%
    dplyr::left_join(dplyr::rename(trait_categories, trait = t))

short_df <- alpha_mix_df %>%
    dplyr::mutate(nh2 = n*h2) %>% 
    dplyr::arrange(desc(nh2)) %>% 
    dplyr::top_n(20) %>%
    dplyr::mutate(w_upper = w_point_estimate + 1.96 * w_se, w_lower = w_point_estimate - 1.96 * w_se)

w_plot <- short_df %>%
    dplyr::arrange(w_point_estimate) %>%
    dplyr::mutate(pretty_name = factor(pretty_name, levels = pretty_name)) %>% 
    ggplot(aes(x = pretty_name, y = w_point_estimate, ymin = w_lower, ymax = w_upper)) + 
    geom_pointrange() + 
    geom_linerange(color = 'red') + 
    coord_flip(ylim = c(0,1), clip = 'on', expand = TRUE) + 
    theme_Publication() + 
    labs(y = 'w (mixture weight)', x = 'Trait') + 
    scale_x_discrete(breaks = c(trait_categories$pretty_name)) + 
    theme(axis.text.y = element_text(size = 10))

alpha_plot <- short_df %>%
    dplyr::arrange(w_point_estimate) %>%
    dplyr::mutate(pretty_name = factor(pretty_name, levels = pretty_name)) %>%
    dplyr::mutate(alpha_lower = alpha_mean - 1.96 * alpha_se, alpha_upper = alpha_mean + 1.96*alpha_se) %>%
    ggplot(aes(x = pretty_name, y = alpha_mean, ymin = alpha_lower, ymax = alpha_upper)) + 
        geom_pointrange() + 
        geom_linerange(color = 'red') + 
        coord_flip(ylim = c(-1,0), clip = 'on', expand = TRUE) + 
        theme_Publication() + 
        scale_x_discrete(breaks = c(trait_categories$pretty_name)) + 
        #geom_hline(yintercept = -0.356066, color = 'red') + 
        #geom_hline(yintercept = -0.4052195257997297, color = 'blue') + 
        labs(y = TeX("$\\alpha_{mix}$"), x = 'Trait') +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) 


combined <- cowplot::plot_grid(
    w_plot, alpha_plot, nrow = 1, rel_widths = c(1, .605) #.6213 too high
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/figure5_short.pdf',
    combined,
    width = 180,
    height = 190,
    units = 'mm'
)




jordans_list <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/jordans_trait_list.txt', sep = '\t', data.table = F) %>%
    dplyr::rename(t = trait_identifier)

trait_categories <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/baselineLD2.2_plus_maf_grid_afr_maf_bins/nh2_selected_traits/thresholded/name_mapping.csv',
        data.table = F) %>% 
    dplyr::rename(t = V1) %>% 
    dplyr::left_join(jordans_list)
f <- '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/alpha_mix_trait_specific.csv'
alpha_mix_df <- data.table::fread(f, data.table = F) %>%
    dplyr::left_join(dplyr::rename(trait_categories, trait = t))


#w_mix_summary_df <- data.table::fread(
#        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/baselineLD2.2_plus_maf_grid_afr_maf_bins/nh2_selected_traits/thresholded/w_mix_summary.txt',
#        data.table = F)
#alpha_summary_df <- data.table::fread(
#        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/baselineLD2.2_plus_maf_grid_afr_maf_bins/nh2_selected_traits/thresholded/alpha_mix_summary.txt',
#        data.table = F)

df <- alpha_mix_df %>%
    dplyr::left_join(trait_categories) %>%
    dplyr::filter(!is.na(pretty_name)) %>%
    dplyr::arrange(category, w_mean) %>%
    dplyr::mutate(w_upper = w_mean + 1.96 * w_se, w_lower = w_mean - 1.96 * w_se)
df_list <- list()
for (i in 1:5) {
    df_cat <- df %>% dplyr::filter(category == unique(df$category)[i])
    if (i %in% 1:4) {
        df_cat <- rbind(df_cat, dplyr::mutate(df_cat[1,], pretty_name = paste0(rep(' ', i), collapse = ''), w_mean = -100, w_lower = -101, w_upper = -99))
    }
    df_list[[i]] <- df_cat
}

w_plot <- purrr::reduce(df_list, rbind) %>%
    dplyr::mutate(pretty_name = factor(pretty_name, levels = pretty_name)) %>% 
    ggplot(aes(x = pretty_name, y = w_mean, ymin = w_lower, ymax = w_upper)) + 
        geom_pointrange() + 
        coord_flip(ylim = c(0,1), clip = 'on', expand = TRUE) + 
        theme_Publication() + 
        labs(y = 'w (mixture weight)', x = 'Trait') + 
        scale_x_discrete(breaks = c(trait_categories$pretty_name)) + 
        theme(axis.text.y = element_text(size = 8))

alpha_plot <- purrr::reduce(df_list, rbind) %>%
    dplyr::mutate(pretty_name = factor(pretty_name, levels = pretty_name)) %>%
    dplyr::mutate(alpha_lower = alpha_mean - 1.96 * alpha_se, alpha_upper = alpha_mean + 1.96*alpha_se) %>%
    dplyr::mutate(
        alpha_mean = ifelse(pretty_name %in% c(' ', '  ', '   ', '    '), 100, alpha_mean), 
        alpha_lower = ifelse(pretty_name %in% c(' ', '  ', '   ', '    '), 100, alpha_lower), 
        alpha_upper = ifelse(pretty_name %in% c(' ', '  ', '   ', '    '), 100, alpha_upper)) %>%
    ggplot(aes(x = pretty_name, y = alpha_mean, ymin = alpha_lower, ymax = alpha_upper)) + 
        geom_pointrange() + 
        coord_flip(ylim = c(-1,0), clip = 'on', expand = TRUE) + 
        theme_Publication() + 
        scale_x_discrete(breaks = c(trait_categories$pretty_name)) + 
        #geom_hline(yintercept = -0.356066, color = 'red') + 
        #geom_hline(yintercept = -0.4052195257997297, color = 'blue') + 
        labs(y = TeX("$\\alpha_{mix}$"), x = 'Trait') +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
    
combined <- cowplot::plot_grid(
    w_plot, alpha_plot, nrow = 1, rel_widths = c(1, .605) #.6213 too high
)
    

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/figure5_indiv.pdf',
    combined,
    width = 180,
    height = 190,
    units = 'mm'
)

# get some statistics
alpha_mix_df %>% 
    dplyr::mutate(w_z = abs(w_mean / w_se)) %>% 
    dplyr::arrange(w_z) %>% 
    dplyr::filter(w_z > -1 * qnorm(.05/2/50, lower.tail = TRUE)) %>% 
    dim()

alpha_mix_df %>% 
    dplyr::mutate(w_z = abs((w_mean - .5)) / w_se) %>% 
    dplyr::arrange(w_z) %>% 
    dplyr::filter(w_z > -1 * qnorm(.05/2/50, lower.tail = TRUE)) %>% 
    dim()

alpha_mix_df %>% 
    dplyr::mutate(w_z = abs((w_mean - 1)) / w_se) %>% 
    dplyr::arrange(w_z) %>% 
    dplyr::filter(w_z > -1 * qnorm(.05/2/50, lower.tail = TRUE)) %>% 
    dim()

alpha_mix_df %>% 
    dplyr::mutate(w_z = abs((w_mean - 1)) / w_se) %>% 
    dplyr::arrange(w_z) %>% 
    dplyr::filter(w_z > -1 * qnorm(.05/2/50, lower.tail = TRUE)) %>% 
    dim()

alpha_mix_df %>% dplyr::arrange(alpha_mean)