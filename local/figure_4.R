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

MAIN_POINT_ESTIMATE <- 0.9608306748172475

#### main fig ----

NULL_VALUE <- 1.3717327963060135
loss_df <- data.frame(
        params = c('$\\alpha=0$', '$\\alpha_{eur}$', '$\\alpha_{afr}$', '$\\alpha_{mix}$'),
        loss = c(1.3717327963060135, 1.3194045111327573, 1.1408945613544303, 1.134858093513253) / NULL_VALUE) %>%
    dplyr::mutate(
        params = factor(params, levels = params)
    )

pmix_params_varying_w <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/nh2_selected_w_unconstrained.alpha_mix_varying_w.txt', data.table = F)
LOSS_LOWER_AXIS_LIMIT <- 0.75
loss_plot <- pmix_params_varying_w %>% 
    dplyr::mutate(loss_lower = loss_point_estimate - 1.96 * loss_se, loss_upper = loss_point_estimate + 1.96 * loss_se) %>%
    dplyr::mutate(loss_point_estimate = loss_point_estimate / NULL_VALUE) %>%
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') + 
    lims(y = c(LOSS_LOWER_AXIS_LIMIT, 1))

loss_bars <- loss_df %>% 
    ggplot(aes(x = params, y = loss)) + 
    #geom_col(color = 'black') +
    geom_point() +
    theme_Publication() + 
    scale_x_discrete(labels = TeX) + 
    labs(y = 'Weighted scaled MSE', x = 'Model') + 
    lims(y = c(LOSS_LOWER_AXIS_LIMIT, 1))

pmix_params_varying_w <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/nh2_selected_w_unconstrained.alpha_mix_varying_w.txt', data.table = F)

alpha_plot <- pmix_params_varying_w %>%
    dplyr::mutate(alpha_lower = alpha_mean - 1.96 * alpha_se, alpha_upper = alpha_mean + 1.96 * alpha_se) %>%
    ggplot(aes(x = w, y = alpha_mean)) + 
    geom_ribbon(aes(ymin = alpha_lower, ymax = alpha_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$")) +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    theme(axis.title.y =element_text(size=13))

combined <- cowplot::plot_grid(
    loss_plot, 
    loss_bars,
    alpha_plot, 
    nrow = 1,
    labels = c('a', 'b', 'c')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/figure_4.pdf',
    combined,
    width = 180,
    height = 70,
    units = 'mm'
)

# p-values

z <- (0.9608306748172475 - 0.5) / 0.103789608110797
p <- (1 - pnorm(abs(z))) * 2
print(p)
z <- (0.9608306748172475 - 0) / 0.103789608110797
p <- (1 - pnorm(abs(z))) * 2
print(p)
z <- (1 - 0.9608306748172475) / 0.103789608110797
p <- (1 - pnorm(abs(z))) * 2
print(p)


#### common domain supplementary figure ----
pmix_params_varying_w_common_domain <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/nh2_selected_w_unconstrained_common_domain.alpha_mix_varying_w.txt', data.table = F)

common_domain_point_estimate <- 0.9499067956772033

loss_plot <- pmix_params_varying_w_common_domain %>% 
    dplyr::mutate(loss_lower = loss_point_estimate - 1.96 * loss_se, loss_upper = loss_point_estimate + 1.96 * loss_se) %>%
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = common_domain_point_estimate), color = 'blue')

alpha_plot <- pmix_params_varying_w_common_domain %>%
    dplyr::mutate(snp_set = 'p_a > 0.05, p_e > 0.05') %>%
    rbind(dplyr::mutate(pmix_params_varying_w, snp_set = 'p_e > 0.05')) %>%
    dplyr::mutate(alpha_lower = alpha_mean - 1.96 * alpha_se, alpha_upper = alpha_mean + 1.96 * alpha_se) %>%
    ggplot(aes(x = w, y = alpha_mean)) + 
    geom_ribbon(aes(ymin = alpha_lower, ymax = alpha_upper, fill = snp_set), alpha = .25) + 
    geom_line(aes(color = snp_set)) + 
    theme_Publication() + 
    scale_color_brewer(palette = "Set1", labels = c(TeX('$p_a \\geq 0.05$ and $p_e \\geq 0.05$'), TeX('$p_e \\geq 0.05$')), direction = -1) + 
    scale_fill_brewer(palette = "Set1", labels = c(TeX('$p_a \\geq 0.05$ and $p_e \\geq 0.05$'), TeX('$p_e \\geq 0.05$')), direction = -1) + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$"), color = '', fill = '') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = common_domain_point_estimate), color = 'blue') +
    theme(axis.title.y =element_text(size=13)) + 
    guides(fill = guide_legend(nrow = 2, position = 'top'), color = guide_legend(nrow = 2, position = 'top'))
    
loss_df <- data.frame(
        params = c('$\\alpha=0$', '$\\alpha_{eur}$', '$\\alpha_{afr}$', '$\\alpha_{mix}$'),
        loss = c(1.1766547962905998, 1.1764970893046067, 1.1721593913719583, 1.172110801277421)) %>%
    dplyr::mutate(
      loss = loss / loss[1],
        params = factor(params, levels = params)
    )

loss_bars <- loss_df %>% 
    ggplot(aes(x = params, y = loss)) + 
    geom_col() +
    theme_Publication() + 
    scale_x_discrete(labels = TeX) + 
    labs(y = 'Weighted scaled MSE', x = 'Model') + 
    theme(axis.text=element_text(size=11)) + 
  coord_cartesian(ylim = c(0.95, 1))

combined <- cowplot::plot_grid(
    loss_bars,
    loss_plot, 
    alpha_plot, 
    nrow = 1,
    labels = c('a', 'b', 'c')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_4_common_domain.pdf',
    combined,
    width = 180,
    height = 80,
    units = 'mm'
)

#### supplementary figure using random effects meta analysis across traits ----
pmix_params_varying_w_re_meta <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/alpha_mix_varying_w_meta.csv', data.table = F)

w_hat_RE <- 0.995031
loss_plot <- pmix_params_varying_w_re_meta %>% 
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = w_hat_RE), color = 'blue') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red')

alpha_plot <- pmix_params_varying_w_re_meta %>%
    ggplot(aes(x = w, y = eff)) + 
    geom_ribbon(aes(ymin = ci_low, ymax = ci_upp), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$")) +
    geom_vline(aes(xintercept = w_hat_RE), color = 'blue') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    theme(axis.title.y =element_text(size=13))

combined <- cowplot::plot_grid(
    loss_plot, 
    alpha_plot, 
    nrow = 1,
    labels = c('a', 'b')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_4_random_effects.pdf',
    combined,
    width = 180,
    height = 80,
    units = 'mm'
)


#### supplementary figure using random effects meta analysis across traits on common domain----
pmix_params_varying_w_re_meta <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/alpha_mix_varying_w_meta_common_domain.csv', data.table = F)

w_hat_RE <- 0.525600
loss_plot <- pmix_params_varying_w_re_meta %>% 
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = w_hat_RE), color = 'blue') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red')

alpha_plot <- pmix_params_varying_w_re_meta %>%
    ggplot(aes(x = w, y = eff)) + 
    geom_ribbon(aes(ymin = ci_low, ymax = ci_upp), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$")) +
    geom_vline(aes(xintercept = w_hat_RE), color = 'blue') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    theme(axis.title.y =element_text(size=13))

pmix_params_varying_w_re_meta_no_phosphate <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/alpha_mix_varying_w_meta_common_domain_no_phosphate.csv', data.table = F)

w_hat_RE_no_phosphate <- 0.525707	
loss_plot_no_phosphate <- pmix_params_varying_w_re_meta_no_phosphate %>% 
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = w_hat_RE_no_phosphate), color = 'blue') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red')

alpha_plot_no_phosphate <- pmix_params_varying_w_re_meta_no_phosphate %>%
    ggplot(aes(x = w, y = eff)) + 
    geom_ribbon(aes(ymin = ci_low, ymax = ci_upp), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$")) +
    geom_vline(aes(xintercept = w_hat_RE_no_phosphate), color = 'blue') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    theme(axis.title.y =element_text(size=13))


combined <- cowplot::plot_grid(
    loss_plot, 
    alpha_plot, 
    loss_plot_no_phosphate,
    alpha_plot_no_phosphate,
    nrow = 2,
    labels = c('a', 'b', 'c', 'd')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_4_random_effects_common_domain.pdf',
    combined,
    width = 180,
    height = 160,
    units = 'mm'
)


#### ukbb supplementary figure ----
pmix_params_varying_w_ukbb <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/nh2_selected_w_unconstrained_ukbiobank_traits.alpha_mix_varying_w.txt', data.table = F)

ukbb_point_estimate <- 0.9657406024552085

loss_plot_ukbb <- pmix_params_varying_w_ukbb %>% 
    dplyr::mutate(loss_lower = loss_point_estimate - 1.96 * loss_se, loss_upper = loss_point_estimate + 1.96 * loss_se) %>%
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = ukbb_point_estimate), color = 'blue')

alpha_plot_ukbb <- pmix_params_varying_w_ukbb %>%
    dplyr::mutate(alpha_lower = alpha_mean - 1.96 * alpha_se, alpha_upper = alpha_mean + 1.96 * alpha_se) %>%
    ggplot(aes(x = w, y = alpha_mean)) + 
    geom_ribbon(aes(ymin = alpha_lower, ymax = alpha_upper), alpha = .25) + 
    geom_line() + 
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$")) +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = ukbb_point_estimate), color = 'blue') +
    theme(axis.title.y =element_text(size=13)) 
    

combined <- cowplot::plot_grid(
    loss_plot_ukbb, 
    alpha_plot_ukbb, 
    nrow = 1,
    labels = c('a', 'b', 'c')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_4_ukbb.pdf',
    combined,
    width = 180,
    height = 80,
    units = 'mm'
)

#### pass supplementary figure ----
pmix_params_varying_w_pass <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid/nh2_selected_w_unconstrained_ukbiobank_traits.alpha_mix_varying_w.txt', data.table = F)

pass_point_estimate <- 0.9584492732875323

loss_plot_pass <- pmix_params_varying_w_pass %>% 
    dplyr::mutate(loss_lower = loss_point_estimate - 1.96 * loss_se, loss_upper = loss_point_estimate + 1.96 * loss_se) %>%
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = pass_point_estimate), color = 'blue')

alpha_plot_pass <- pmix_params_varying_w_pass %>%
    dplyr::mutate(alpha_lower = alpha_mean - 1.96 * alpha_se, alpha_upper = alpha_mean + 1.96 * alpha_se) %>%
    ggplot(aes(x = w, y = alpha_mean)) + 
    geom_ribbon(aes(ymin = alpha_lower, ymax = alpha_upper), alpha = .25) + 
    geom_line() + 
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$")) +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = pass_point_estimate), color = 'blue') +
    theme(axis.title.y =element_text(size=13)) 
    

combined <- cowplot::plot_grid(
    loss_plot_pass, 
    alpha_plot_pass, 
    nrow = 1,
    labels = c('a', 'b', 'c')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_4_pass.pdf',
    combined,
    width = 180,
    height = 80,
    units = 'mm'
)

#### YRI MAF  ----
pmix_params_varying_w_yri <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2_plus_maf_grid_yri/nh2_selected_w_unconstrained.alpha_mix_varying_w.txt', data.table = F)

yri_point_estimate <- 0.9561042056042768

loss_plot_yri <- pmix_params_varying_w_yri %>% 
    dplyr::mutate(loss_lower = loss_point_estimate - 1.96 * loss_se, loss_upper = loss_point_estimate + 1.96 * loss_se) %>%
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = yri_point_estimate), color = 'blue')

alpha_plot_yri <- pmix_params_varying_w_yri %>%
    dplyr::mutate(alpha_lower = alpha_mean - 1.96 * alpha_se, alpha_upper = alpha_mean + 1.96 * alpha_se) %>%
    ggplot(aes(x = w, y = alpha_mean)) + 
    geom_ribbon(aes(ymin = alpha_lower, ymax = alpha_upper), alpha = .25) + 
    geom_line() + 
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$")) +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = yri_point_estimate), color = 'blue') +
    theme(axis.title.y =element_text(size=13)) 
    

combined <- cowplot::plot_grid(
    loss_plot_yri, 
    alpha_plot_yri, 
    nrow = 1,
    labels = c('a', 'b', 'c')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_4_yri.pdf',
    combined,
    width = 180,
    height = 80,
    units = 'mm'
)

##### Grid only, no baselineLD  ----
pmix_params_varying_w_grid_only <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/maf_grid/nh2_selected_w_unconstrained.alpha_mix_varying_w.txt', data.table = F)

grid_only_point_estimate <- 0.9798801535662633

loss_plot_grid_only <- pmix_params_varying_w_grid_only %>% 
    dplyr::mutate(loss_lower = loss_point_estimate - 1.96 * loss_se, loss_upper = loss_point_estimate + 1.96 * loss_se) %>%
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = grid_only_point_estimate), color = 'blue')

alpha_plot_grid_only <- pmix_params_varying_w_grid_only %>%
    dplyr::mutate(alpha_lower = alpha_mean - 1.96 * alpha_se, alpha_upper = alpha_mean + 1.96 * alpha_se) %>%
    ggplot(aes(x = w, y = alpha_mean)) + 
    geom_ribbon(aes(ymin = alpha_lower, ymax = alpha_upper), alpha = .25) + 
    geom_line() + 
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$")) +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = grid_only_point_estimate), color = 'blue') +
    theme(axis.title.y =element_text(size=13)) 
    

combined <- cowplot::plot_grid(
    loss_plot_grid_only, 
    alpha_plot_grid_only, 
    nrow = 1,
    labels = c('a', 'b', 'c')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_4_grid_only.pdf',
    combined,
    width = 180,
    height = 80,
    units = 'mm'
)

#### baselineLD only ----
pmix_params_varying_w_baseline_only <- data.table::fread('~/Google Drive/My Drive/research/alkes/h2xancestry/data_from_o2/real_traits_2/baselineLD2.2/nh2_selected_w_unconstrained.alpha_mix_varying_w.txt', data.table = F)

baseline_only_point_estimate <- 0.8543494480978265

loss_plot_baseline_only <- pmix_params_varying_w_baseline_only %>% 
    dplyr::mutate(loss_lower = loss_point_estimate - 1.96 * loss_se, loss_upper = loss_point_estimate + 1.96 * loss_se) %>%
    ggplot(aes(x = w, y = loss_point_estimate)) + 
    #geom_ribbon(aes(ymin = loss_lower, ymax = loss_upper), fill = 'grey70') + 
    geom_line() + 
    theme_Publication() + 
    labs(y = 'Weighted Scaled MSE', x = 'w (mixture weight)') +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = baseline_only_point_estimate), color = 'blue')

alpha_plot_baseline_only <- pmix_params_varying_w_baseline_only %>%
    dplyr::mutate(alpha_lower = alpha_mean - 1.96 * alpha_se, alpha_upper = alpha_mean + 1.96 * alpha_se) %>%
    ggplot(aes(x = w, y = alpha_mean)) + 
    geom_ribbon(aes(ymin = alpha_lower, ymax = alpha_upper), alpha = .25) + 
    geom_line() + 
    theme_Publication() + 
    labs(x = 'w (mixture weight)', y = TeX("$\\alpha_{mix}$")) +
    geom_vline(aes(xintercept = MAIN_POINT_ESTIMATE), color = 'red') +
    geom_vline(aes(xintercept = baseline_only_point_estimate), color = 'blue') +
    theme(axis.title.y =element_text(size=13)) 
    

combined <- cowplot::plot_grid(
    loss_plot_baseline_only, 
    alpha_plot_baseline_only, 
    nrow = 1,
    labels = c('a', 'b', 'c')
)

ggsave(
    '~/Google Drive/My Drive/research/alkes/h2xancestry/chapter_3_figures/supplementary_figures/figure_4_baseline_only.pdf',
    combined,
    width = 180,
    height = 80,
    units = 'mm'
)
