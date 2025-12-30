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

for (suffix in c('thresholded_.05', 'n_causal_5000', 'thresholded_MAC1', 'n_causal_20000', 'thresholded_.005_indiv', 'ukbb_british_maf')) {
  f <- paste0(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/simulations_2/summaries/alpha_.38_h2_.5_',
    suffix, 
    '.csv'
  )
  w_df <- data.table::fread(
      f, 
      data.table = F
  )
  # new version with bars
  w_figure <- w_df %>% 
    dplyr::group_by(true_w) %>% 
    dplyr::summarise(
      mean_w = mean(w),
      ci_dist = 1.96 * sd(w) / sqrt(n())) %>%
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
    
  f <- paste0(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/figure_6_',
    suffix,
    '.pdf'
  )
  ggsave(
      f,
      combined_main_figure,
      width = 180,
      height = 80,
      units = 'mm'
  )
}
