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

#### main figure ----

ns_p_mix_params_grid <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_p_mix_params_ge.002_in_either_afr_eur.csv', 
        data.table = F) %>% 
    dplyr::mutate(
        beta_upper = beta + 1.96 * beta_se,
        beta_lower = beta - 1.96 * beta_se
    )

ll_bar_df <- data.frame(
        params = c('$\\gamma=0$', '$w=0$', '$w=1$', '$w_{MLE}$'),
        ll = c(-763071.8895360784, -761976.2, -759157.9, -759070.5886847276)) %>%
    dplyr::mutate(
        params = factor(params, levels = params)
    )

ll_bar_df %>% dplyr::mutate(delta = ll - ll[params == '$\\gamma=0$'])

loglikelihood_plot <- ns_p_mix_params_grid %>%
    ggplot(aes(x = w, y = l * -1)) + 
    geom_line() + 
    labs(x = 'w (mixture weight)', y = 'Negative log-likelihood') + 
    theme_Publication() + 
    lims(y = -1 * rev(c(min(ll_bar_df$ll), max(ll_bar_df$ll) + 1000))) +
    geom_vline(xintercept = 0.9466057164823091, color = 'red')

beta_plot <- ns_p_mix_params_grid %>%
    ggplot(aes(x = w, y = beta)) + 
    geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper), fill = 'grey70') + 
    geom_line() +
    geom_vline(xintercept = 0.9466057164823091, color = 'red') +
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = 'Slope \n(effect of log MAF on log odds)') 


loss_bars <- ll_bar_df %>% 
    ggplot(aes(x = params, y = ll * -1)) + 
    geom_point() +
    theme_Publication() + 
    scale_x_discrete(labels = TeX) + 
    labs(y = 'Negative log-likelihood', x = '') + 
    theme(axis.text=element_text(size=11)) + 
    lims(y = -1 * rev(c(min(ll_bar_df$ll), max(ll_bar_df$ll) + 1000))) +
    #coord_flip() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text = element_text(size = 9))




combined <- cowplot::plot_grid(
    loglikelihood_plot, 
    loss_bars,
    beta_plot, 
    nrow = 1,
    labels = c('a', 'b', 'c')
)


ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/figure_2.pdf',
    combined,
    width = 180,
    height = 70,
    units = 'mm'
)

#### supplementary figure for pmix scaled ----
ns_p_mix_params_grid <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_p_mix_params_ge.002_in_either_afr_eur.csv', 
        data.table = F) %>% 
    dplyr::mutate(
        beta_upper = beta + 1.96 * beta_se,
        beta_lower = beta - 1.96 * beta_se) %>% 
    dplyr::mutate(
        normalization = 
            'variance unnormalized'
    )
ns_p_mix_params_grid_scaled <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_p_mix_params_ge.002_in_either_afr_eur_scaled.csv', 
        data.table = F) %>% 
    dplyr::mutate(
        beta_upper = beta + 1.96 * beta_se,
        beta_lower = beta - 1.96 * beta_se
    ) %>% dplyr::mutate(
        normalization = 'variance = 1'
    )
beta_plot <- rbind(ns_p_mix_params_grid, ns_p_mix_params_grid_scaled) %>%
    ggplot(aes(x = w, y = beta, color = normalization, fill = normalization)) + 
    geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper), color = 'transparent', alpha = .25) + 
    geom_line() +
    geom_vline(xintercept = 0.9466057164823091, color = 'red') +
    scale_color_brewer(palette = 'Set1', direction = -1) +
    scale_fill_brewer(palette = 'Set1', direction = -1) +
    theme_Publication() + 
    labs(x = 'w (Mixture weight)', y = 'Slope \n(effect of log MAF on log odds)', color = '', fill = '') 


ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_2_scaled_pmix.pdf',
    beta_plot,
    width = 90,
    height = 120,
    units = 'mm'
)


#### supplementary figure with alt threshold ----
ns_p_mix_params_grid <- data.table::fread(
        '~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_aou/ns_p_mix_params_ge.002_in_either_afr_eur_05_thresholded.csv', 
        data.table = F) %>% 
    dplyr::mutate(
        beta_upper = beta + 1.96 * beta_se,
        beta_lower = beta - 1.96 * beta_se
    )


loglikelihood_plot <- ns_p_mix_params_grid %>%
    ggplot(aes(x = w, y = l)) + 
    geom_line() + 
    labs(x = 'w (mixture weight)', y = 'Log-likelihood') + 
    theme_Publication() + 
    geom_vline(xintercept = 0.80, color = 'red')

beta_plot <- ns_p_mix_params_grid %>%
    ggplot(aes(x = w, y = beta)) + 
    geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper), fill = 'grey70') + 
    geom_line() +
    geom_vline(xintercept = 0.80, color = 'red') +
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = 'Slope \n(effect of log MAF on log odds)') 





combined <- cowplot::plot_grid(
    loglikelihood_plot, 
    beta_plot, 
    nrow = 1,
    labels = c('a', 'b')
)


ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_2_alt_threshold.pdf',
    combined,
    width = 180,
    height = 70,
    units = 'mm'
)
