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
w_df <- data.table::fread(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/simulations_2/summaries/alpha_.38_h2_.5_thresholded_.005.csv', 
    data.table = F
)
# new version with bars
w_figure <- w_df %>% 
  dplyr::group_by(true_w) %>% 
  dplyr::summarise(
    mean_w = mean(w),
    ci_dist = 1.96 * sd(w) / sqrt(40)) %>%
  dplyr::mutate(
    w_lower = mean_w - ci_dist,
    w_upper = mean_w + ci_dist) %>%
  dplyr::mutate(true_w_label = factor(true_w)) %>%
  ggplot(aes(x = true_w_label, y = mean_w, ymax = w_upper, ymin = w_lower)) +
  #geom_pointrange(color = 'black') + 
  geom_col(color = 'black') +
  geom_errorbar(color = 'black', width = .4) + 
  theme_Publication() +
  lims(y = c(-.05,1.05)) +
  geom_segment(aes(y = true_w - .006, yend = true_w + .006, xend = true_w_label), linewidth = 7.5, color = 'red') +
  labs(y = 'Estimated w (mixture weight)', x = 'True w (mixture weight)')

alpha_figure <- w_df %>% 
  dplyr::group_by(true_w) %>% 
  dplyr::summarise(
    mean_alpha = mean(alpha),
    ci_dist = 1.96 * sd(alpha) / sqrt(40)) %>%
  dplyr::mutate(
    alpha_lower = mean_alpha - ci_dist,
    alpha_upper = mean_alpha + ci_dist) %>%
  dplyr::mutate(true_w_label = factor(true_w)) %>%
  ggplot(aes(x = true_w_label, y = mean_alpha, ymax = alpha_upper, ymin = alpha_lower)) +
  #geom_pointrange(color = 'black') + 
  geom_col(color = 'black') +
  geom_errorbar(color = 'black', width = .4) + 
  theme_Publication() +
  geom_segment(aes(y = -.38 - .002, yend = -.38+ .002, xend = true_w_label), linewidth = 7.5, color = 'red') +
  labs(y = TeX("Estimated $\\alpha_{mix}$"), x = 'True w (mixture weight)')

combined_main_figure <- cowplot::plot_grid(
    w_figure, 
    alpha_figure,
    nrow = 1,
    labels = c('a', 'b')
)
  
ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/figure_6.pdf',
    combined_main_figure,
    width = 180,
    height = 80,
    units = 'mm'
)

# old version with boxplots
  w_figure <- w_df %>% 
    dplyr::mutate(true_w_label = factor(true_w)) %>%
    {ggplot(., aes(x = true_w_label, y = w)) + 
    geom_boxplot(outliers = FALSE) + 
    geom_jitter(size = .7, alpha = .3) +
    geom_segment(aes(y = true_w - .006, yend = true_w + .006, xend = true_w_label), data = dplyr::filter(., seed_start == 0), linewidth = 7.5, color = 'red') +
    lims(y = c(-.05,1.05)) +
    facet_grid(~true_w_label, scale = 'free_x') +
    labs(y = 'Estimated w (mixture weight)', x = 'True w (mixture weight)') +
    theme_Publication() +
    theme(panel.spacing.x = unit(0, "lines"), strip.background = element_blank(), strip.text.x = element_blank())}

w_figure <- w_df %>% 
    dplyr::mutate(true_w_label = factor(true_w)) %>%
    {ggplot(., aes(x = true_w_label, y = w)) + 
    geom_boxplot(outliers = FALSE) + 
    geom_jitter(size = .7, alpha = .3) +
    geom_segment(aes(y = true_w - .006, yend = true_w + .006, xend = true_w_label), data = dplyr::filter(., seed_start == 0), linewidth = 7.5, color = 'red') +
    lims(y = c(-.05,1.05)) +
    facet_grid(~true_w_label, scale = 'free_x') +
    labs(y = 'Estimated w (mixture weight)', x = 'True w (mixture weight)') +
    theme_Publication() +
    theme(panel.spacing.x = unit(0, "lines"), strip.background = element_blank(), strip.text.x = element_blank())}

alpha_figure <- w_df %>% 
    dplyr::mutate(true_w_label = factor(true_w)) %>%
    {ggplot(., aes(x = true_w_label, y = alpha)) + 
    geom_boxplot(outliers = FALSE) + 
    geom_jitter(size = .7, alpha = .3) +
    geom_segment(aes(y = -.38 - .002, yend = -.38+ .002, xend = true_w_label), data = dplyr::filter(., seed_start == 0), linewidth = 7.5, color = 'red') +
    facet_grid(~true_w_label, scale = 'free_x') +
    labs(y = TeX("Estimated $\\alpha_{mix}$"), x = 'True w (mixture weight)') +
    theme_Publication() +
    lims(y = c(NA, 0)) +
    theme(panel.spacing.x = unit(0, "lines"), strip.background = element_blank(), strip.text.x = element_blank())}

combined_main_figure <- cowplot::plot_grid(
    w_figure, 
    alpha_figure,
    nrow = 1,
    labels = c('a', 'b')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/figure_6.pdf',
    combined_main_figure,
    width = 180,
    height = 80,
    units = 'mm'
)


vary_threshold_w_plot <- w_df %>% 
    tidyr::complete(seed_start, true_w, threshold) %>% 
    dplyr::left_join(dplyr::transmute(dplyr::filter(w_df, threshold == 0), true_w, seed_start, c_mean_T0 = c_mean, alpha_mean_T0 = alpha_mean)) %>%
    dplyr::mutate(c_mean = ifelse(is.na(c_mean), c_mean_T0, c_mean)) %>%
    dplyr::mutate(alpha_mean = ifelse(is.na(alpha_mean), alpha_mean_T0, alpha_mean)) %>%
    dplyr::mutate(true_w_label = factor(true_w)) %>%
    {ggplot(., aes(x = true_w_label, y = c_mean)) + 
    geom_boxplot() + 
    geom_jitter(size = .5, alpha = .2) +
    geom_segment(aes(y = true_w - .006, yend = true_w + .006, xend = true_w_label), data = dplyr::filter(., seed_start == 0), linewidth = 7.5, color = 'red') +
    lims(y = c(-.05,1.05)) +
    facet_grid(. ~ threshold) +
    labs(y = 'Estimated w (mixture weight)', x = 'True w (mixture weight)') +
    theme_Publication()}

vary_threshold_alpha_plot <- w_df %>% 
    tidyr::complete(seed_start, true_w, threshold) %>% 
    dplyr::left_join(dplyr::transmute(dplyr::filter(w_df, threshold == 0), true_w, seed_start, c_mean_T0 = c_mean, alpha_mean_T0 = alpha_mean)) %>%
    dplyr::mutate(c_mean = ifelse(is.na(c_mean), c_mean_T0, c_mean)) %>%
    dplyr::mutate(alpha_mean = ifelse(is.na(alpha_mean), alpha_mean_T0, alpha_mean)) %>%
    dplyr::mutate(true_w_label = factor(true_w)) %>%
    {ggplot(., aes(x = true_w_label, y = alpha_mean)) + 
    geom_boxplot() + 
    geom_jitter(size = .5, alpha = .2) +
    #geom_segment(aes(y = true_w - .004, yend = true_w + .004, xend = true_w_label), data = dplyr::filter(., seed_start == 0), linewidth = 12, color = 'red') +
    geom_segment(aes(y = -.38 - .002, yend = -.38+ .002, xend = true_w_label), data = dplyr::filter(., seed_start == 0), linewidth = 7.5, color = 'red') +
    facet_grid(. ~ threshold, scale = 'free_x') +
    labs(y = 'Estimated alpha', x = 'True w (mixture weight)') +
    theme_Publication()}

combined_vary_threshold_plot <- cowplot::plot_grid(
    vary_threshold_w_plot,
    vary_threshold_alpha_plot,
    labels = c('a', 'b'),
    nrow = 2
)
ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/vary_threshold_simulation_figure.pdf',
    combined_vary_threshold_plot,
    units = 'mm',
    width = 180, 
    height = 140
)