#############################################
####                                     ####
####      Winner's curse                 ####
####      Simulations - Figures, etc     ####
####                                     ####
####                                     ####
####      Ninon Mounier                  ####
####      10/01/2023                     ####
#############################################

# assume than in the current directory, there is a "results/" folder
# where results have been be saved
working_dir = "results/"
setwd(working_dir)

library(tidyverse)

library(ggplot2)
theme_set(theme_bw())
library(MetBrewer)

library(ggpubr)


my_colors = met.brewer("Archambault", 100)[c(1,10,
                                             18,36,51,
                                             68,75,86,100)]


my_breaks6 = c("IVW", "IVW_3samples",
                "RAPS_3samples", "RAPS_ZP", "RAPS_UMVCUE",
                "IVW_RC")
my_labels6 = c("Naive IVW", "IVW, 3-sample design", "MR-RAPS, 3-sample design",
                "MR-RAPS, ZP", "MR-RAPS, UMVCUE",
                "IVW, regression calibration")
my_colors6 = my_colors[c(1,2,3,4,5,6,7,7,8,8)]
my_lines6 = c(rep("solid",6))


my_breaks10 = c("IVW", "IVW_3samples",
                "RAPS_3samples", "RAPS_ZP", "RAPS_UMVCUE",
                "IVW_RC", "MedianBased_RC", "WMedianBased_RC", "ModeBased_RC","WModeBased_RC")
my_labels10 = c("Naive IVW", "IVW, 3-sample design", "MR-RAPS, 3-sample design",
                "MR-RAPS, ZP", "MR-RAPS, UMVCUE",
                "IVW, regression calibration",
                "Simple median, regression calibration", "Weighted median, regression calibration",
                "Simple mode, regression calibration", "Weighted mode, regression calibration")
my_colors10 = my_colors[c(1,2,3,4,5,6,7,7,8,8)]
my_lines10 = c(rep("solid", 6), "dashed", "dotted", "dashed", "dotted")

#### Results ####

# print summary at 5e-8 (estimate + sd for each method, as well as number of IVs / F-stats)
get_results <- function(scenario, pleiotropy=F){
  res = readRDS(paste0("res_sim_", scenario, ".RDS"))

  ## mean / sd + number of IVs and F-stats
  if(pleiotropy){
    res %>%
      dplyr::transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, MedianBased_RC, WMedianBased_RC, ModeBased_RC, WModeBased_RC) %>%
      tidyr::pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, MedianBased_RC, WMedianBased_RC, ModeBased_RC, WModeBased_RC)) %>%
      mutate(name = lvls_reorder(as_factor(name), idx = match(my_breaks10, levels(as_factor(name))))) -> res_tidy
    res_tidy %>%
      group_by(name) %>%
      arrange(name) %>%
      summarise(m = mean(value), s = sd(value), bias = mean(value) - 0.2, rmse = sqrt(mean((value - 0.2)^2))) %>% as.data.frame %>% print()

    res %>%
      summarise(mean(m_combined), mean(m), mean(F_combined), mean(F_disc), mean(F_rep),
                m_pleio_disc = mean(pleiotropic_IVs_disc), perc_pleio_disc = m_pleio_disc/mean(m)*100,
                m_pleio_combined = mean(pleiotropic_IVs_combined), perc_pleio_combined = m_pleio_combined/mean(m_combined)*100) %>%
      as.data.frame %>% print()

  } else {
    res %>%
    dplyr::transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC) %>%
      tidyr::pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) %>%
      mutate(name = lvls_reorder(as_factor(name), idx = match(my_breaks6, levels(as_factor(name))))) -> res_tidy

    res_tidy %>%
      group_by(name) %>%
      arrange(name) %>%
      summarise(m = mean(value), s = sd(value), bias = mean(value) - 0.2, rmse = sqrt(mean((value - 0.2)^2))) %>% as.data.frame %>% print()

    res %>%
      summarise(mean(m_combined), mean(m), mean(F_combined), mean(F_disc), mean(F_rep))  %>%
      as.data.frame %>% print()
  }


}

get_results("balancedpleio")
get_results("correlatedpleio_strong", TRUE)
get_results("correlatedpleio_weak", TRUE)


# copy results to clipboard for all thresholds (then directly paste in supplementary file)
get_results_all <- function(scenario, pleiotropy=F){
  all_res = readRDS(paste0("res_varyT_", scenario, ".RDS"))

  ## mean / sd + number of IVs and F-stats
  if(pleiotropy){
    all_res %>%
      group_by(threshold) %>%
      summarize(m_combined = mean(m_combined),
                perc_pleio_combined = mean(pleiotropic_IVs_combined)/m_combined*100,
                m_disc = mean(m),
                perc_pleio_disc = mean(pleiotropic_IVs_disc)/m_disc*100,
                F_combined = mean(F_combined),
                F_disc = mean(F_disc),
                F_rep = mean(F_rep),
                IVW_beta = mean(IVW),
                IVW_sd = sd(IVW),
                IVW_bias = mean(IVW) - 0.2,
                IVW_rmse = sqrt(mean((IVW - 0.2)^2)),
                IVW_3samples_beta = mean(IVW_3samples, na.rm=T),
                IVW_3samples_sd = sd(IVW_3samples, na.rm=T),
                IVW_3samples_bias = mean(IVW_3samples, na.rm=T) - 0.2,
                IVW_3samples_rmse = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
                RAPS_3samples_beta = mean(RAPS_3samples, na.rm=T),
                RAPS_3samples_sd = sd(RAPS_3samples, na.rm=T),
                RAPS_3samples_bias = mean(RAPS_3samples, na.rm=T) - 0.2,
                RAPS_3samples_rmse = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
                RAPS_ZP_beta = mean(RAPS_ZP, na.rm=T),
                RAPS_ZP_sd = sd(RAPS_ZP, na.rm=T),
                RAPS_ZP_bias = mean(RAPS_ZP, na.rm=T) - 0.2,
                RAPS_ZP_rmse = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
                RAPS_UMVCUE_beta = mean(RAPS_UMVCUE, na.rm=T),
                RAPS_UMVCUE_sd = sd(RAPS_UMVCUE, na.rm=T),
                RAPS_UMVCUE_bias = mean(RAPS_UMVCUE, na.rm=T) - 0.2,
                RAPS_UMVCUE_rmse = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
                IVW_RC_beta = mean(IVW_RC, na.rm=T),
                IVW_RC_sd = sd(IVW_RC, na.rm=T),
                IVW_RC_bias = mean(IVW_RC, na.rm=T) - 0.2,
                IVW_RC_rmse = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T)),
                MedianBased_RC_beta = mean(MedianBased_RC, na.rm=T),
                MedianBased_RC_sd = sd(MedianBased_RC, na.rm=T),
                MedianBased_RC_bias = mean(MedianBased_RC, na.rm=T) - 0.2,
                MedianBased_RC_rmse = sqrt(mean((MedianBased_RC - 0.2)^2, na.rm=T)),
                WMedianBased_RC_beta = mean(WMedianBased_RC, na.rm=T),
                WMedianBased_RC_sd = sd(WMedianBased_RC, na.rm=T),
                WMedianBased_RC_bias = mean(WMedianBased_RC, na.rm=T) - 0.2,
                WMedianBased_RC_rmse = sqrt(mean((WMedianBased_RC - 0.2)^2, na.rm=T)),
                ModeBased_RC_beta = mean(ModeBased_RC, na.rm=T),
                ModeBased_RC_sd = sd(ModeBased_RC, na.rm=T),
                ModeBased_RC_bias = mean(ModeBased_RC, na.rm=T) - 0.2,
                ModeBased_RC_rmse = sqrt(mean((ModeBased_RC - 0.2)^2, na.rm=T)),
                WModeBased_RC_beta = mean(WModeBased_RC, na.rm=T),
                WModeBased_RC_sd = sd(WModeBased_RC, na.rm=T),
                WModeBased_RC_bias = mean(WModeBased_RC, na.rm=T) - 0.2,
                WModeBased_RC_rmse = sqrt(mean((WModeBased_RC - 0.2)^2, na.rm=T))) %>%
      write.table("clipboard",sep="\t", row.names = F) # to copy to clipboard and paste in excel supplementary file
  } else {
    all_res %>%
      group_by(threshold) %>%
      summarize(m_combined = mean(m_combined),
                m_disc = mean(m),
                F_combined = mean(F_combined),
                F_disc = mean(F_disc),
                F_rep = mean(F_rep),
                IVW_beta = mean(IVW),
                IVW_sd = sd(IVW),
                IVW_bias = mean(IVW) - 0.2,
                IVW_rmse = sqrt(mean((IVW - 0.2)^2)),
                IVW_3samples_beta = mean(IVW_3samples, na.rm=T),
                IVW_3samples_sd = sd(IVW_3samples, na.rm=T),
                IVW_3samples_bias = mean(IVW_3samples, na.rm=T) - 0.2,
                IVW_3samples_rmse = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
                RAPS_3samples_beta = mean(RAPS_3samples, na.rm=T),
                RAPS_3samples_sd = sd(RAPS_3samples, na.rm=T),
                RAPS_3samples_bias = mean(RAPS_3samples, na.rm=T) - 0.2,
                RAPS_3samples_rmse = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
                RAPS_ZP_beta = mean(RAPS_ZP, na.rm=T),
                RAPS_ZP_sd = sd(RAPS_ZP, na.rm=T),
                RAPS_ZP_bias = mean(RAPS_ZP, na.rm=T) - 0.2,
                RAPS_ZP_rmse = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
                RAPS_UMVCUE_beta = mean(RAPS_UMVCUE, na.rm=T),
                RAPS_UMVCUE_sd = sd(RAPS_UMVCUE, na.rm=T),
                RAPS_UMVCUE_bias = mean(RAPS_UMVCUE, na.rm=T) - 0.2,
                RAPS_UMVCUE_rmse = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
                IVW_RC_beta = mean(IVW_RC, na.rm=T),
                IVW_RC_sd = sd(IVW_RC, na.rm=T),
                IVW_RC_bias = mean(IVW_RC, na.rm=T) - 0.2,
                IVW_RC_rmse = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T))) %>%
      write.table("clipboard",sep="\t", row.names = F) # to copy to clipboard and paste in excel supplementary file

  }

}

get_results_all("balancedpleio")
get_results_all("correlatedpleio_strong", TRUE)
get_results_all("correlatedpleio_weak", TRUE)


# copy results to clipboard for all proportions of discovery/replication (then directly paste in supplementary file)
get_results_allN <- function(scenario){
  all_res_N = readRDS(paste0("res_varyN_", scenario, ".RDS"))

  ## mean / sd + number of IVs and F-stats
  all_res_N %>%
    group_by(prop) %>%
    summarize(m_combined = mean(m_combined),
              m_disc = mean(m),
              F_combined = mean(F_combined),
              F_disc = mean(F_disc),
              F_rep = mean(F_rep),
              IVW_beta = mean(IVW),
              IVW_sd = sd(IVW),
              IVW_bias = mean(IVW) - 0.2,
              IVW_rmse = sqrt(mean((IVW - 0.2)^2)),
              IVW_3samples_beta = mean(IVW_3samples, na.rm=T),
              IVW_3samples_sd = sd(IVW_3samples, na.rm=T),
              IVW_3samples_bias = mean(IVW_3samples, na.rm=T) - 0.2,
              IVW_3samples_rmse = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
              RAPS_3samples_beta = mean(RAPS_3samples, na.rm=T),
              RAPS_3samples_sd = sd(RAPS_3samples, na.rm=T),
              RAPS_3samples_bias = mean(RAPS_3samples, na.rm=T) - 0.2,
              RAPS_3samples_rmse = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
              RAPS_ZP_beta = mean(RAPS_ZP, na.rm=T),
              RAPS_ZP_sd = sd(RAPS_ZP, na.rm=T),
              RAPS_ZP_bias = mean(RAPS_ZP, na.rm=T) - 0.2,
              RAPS_ZP_rmse = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
              RAPS_UMVCUE_beta = mean(RAPS_UMVCUE, na.rm=T),
              RAPS_UMVCUE_sd = sd(RAPS_UMVCUE, na.rm=T),
              RAPS_UMVCUE_bias = mean(RAPS_UMVCUE, na.rm=T) - 0.2,
              RAPS_UMVCUE_rmse = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
              IVW_RC_beta = mean(IVW_RC, na.rm=T),
              IVW_RC_sd = sd(IVW_RC, na.rm=T),
              IVW_RC_bias = mean(IVW_RC, na.rm=T) - 0.2,
              IVW_RC_rmse = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T))) %>%
    write.table("clipboard",sep="\t", row.names = F) # to copy to clipboard and paste in excel supplementary file


}

get_results_allN("balancedpleio")
# get_results_allN("balancedpleio_largerN")


# copy results to clipboard for all (total) exposure sample sizes
get_results_Ntotexp <- function(scenario){
  all_res_N = readRDS(paste0("res_varyNtotexp_", scenario, ".RDS"))

  ## mean / sd + number of IVs and F-stats
  all_res_N %>%
    group_by(tot_N) %>%
    summarize(m_combined = mean(m_combined),
              m_disc = mean(m),
              F_combined = mean(F_combined),
              F_disc = mean(F_disc),
              F_rep = mean(F_rep),
              IVW_beta = mean(IVW),
              IVW_sd = sd(IVW),
              IVW_bias = mean(IVW) - 0.2,
              IVW_rmse = sqrt(mean((IVW - 0.2)^2)),
              IVW_3samples_beta = mean(IVW_3samples, na.rm=T),
              IVW_3samples_sd = sd(IVW_3samples, na.rm=T),
              IVW_3samples_bias = mean(IVW_3samples, na.rm=T) - 0.2,
              IVW_3samples_rmse = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
              RAPS_3samples_beta = mean(RAPS_3samples, na.rm=T),
              RAPS_3samples_sd = sd(RAPS_3samples, na.rm=T),
              RAPS_3samples_bias = mean(RAPS_3samples, na.rm=T) - 0.2,
              RAPS_3samples_rmse = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
              RAPS_ZP_beta = mean(RAPS_ZP, na.rm=T),
              RAPS_ZP_sd = sd(RAPS_ZP, na.rm=T),
              RAPS_ZP_bias = mean(RAPS_ZP, na.rm=T) - 0.2,
              RAPS_ZP_rmse = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
              RAPS_UMVCUE_beta = mean(RAPS_UMVCUE, na.rm=T),
              RAPS_UMVCUE_sd = sd(RAPS_UMVCUE, na.rm=T),
              RAPS_UMVCUE_bias = mean(RAPS_UMVCUE, na.rm=T) - 0.2,
              RAPS_UMVCUE_rmse = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
              IVW_RC_beta = mean(IVW_RC, na.rm=T),
              IVW_RC_sd = sd(IVW_RC, na.rm=T),
              IVW_RC_bias = mean(IVW_RC, na.rm=T) - 0.2,
              IVW_RC_rmse = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T))) %>%
    write.table("clipboard",sep="\t", row.names = F) # to copy to clipboard and paste in excel supplementary file


}

get_results_Ntotexp("balancedpleio")

# copy results to clipboard for all outcome sample sizes
get_results_Ny <- function(scenario){
  all_res_N = readRDS(paste0("res_varyNy_", scenario, ".RDS"))

  ## mean / sd + number of IVs and F-stats
  all_res_N %>%
    group_by(Ny) %>%
    summarize(m_combined = mean(m_combined),
              m_disc = mean(m),
              F_combined = mean(F_combined),
              F_disc = mean(F_disc),
              F_rep = mean(F_rep),
              IVW_beta = mean(IVW),
              IVW_sd = sd(IVW),
              IVW_bias = mean(IVW) - 0.2,
              IVW_rmse = sqrt(mean((IVW - 0.2)^2)),
              IVW_3samples_beta = mean(IVW_3samples, na.rm=T),
              IVW_3samples_sd = sd(IVW_3samples, na.rm=T),
              IVW_3samples_bias = mean(IVW_3samples, na.rm=T) - 0.2,
              IVW_3samples_rmse = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
              RAPS_3samples_beta = mean(RAPS_3samples, na.rm=T),
              RAPS_3samples_sd = sd(RAPS_3samples, na.rm=T),
              RAPS_3samples_bias = mean(RAPS_3samples, na.rm=T) - 0.2,
              RAPS_3samples_rmse = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
              RAPS_ZP_beta = mean(RAPS_ZP, na.rm=T),
              RAPS_ZP_sd = sd(RAPS_ZP, na.rm=T),
              RAPS_ZP_bias = mean(RAPS_ZP, na.rm=T) - 0.2,
              RAPS_ZP_rmse = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
              RAPS_UMVCUE_beta = mean(RAPS_UMVCUE, na.rm=T),
              RAPS_UMVCUE_sd = sd(RAPS_UMVCUE, na.rm=T),
              RAPS_UMVCUE_bias = mean(RAPS_UMVCUE, na.rm=T) - 0.2,
              RAPS_UMVCUE_rmse = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
              IVW_RC_beta = mean(IVW_RC, na.rm=T),
              IVW_RC_sd = sd(IVW_RC, na.rm=T),
              IVW_RC_bias = mean(IVW_RC, na.rm=T) - 0.2,
              IVW_RC_rmse = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T))) %>%
    write.table("clipboard",sep="\t", row.names = F) # to copy to clipboard and paste in excel supplementary file


}

get_results_Ny("balancedpleio")


#### Figures ####

# print figure and create pdf + png files at 5e-8
make_figure <- function(scenario, pleiotropy = F){

  res = readRDS(paste0("res_sim_", scenario, ".RDS"))

  ## Figure
  if(pleiotropy){

    res %>%
      dplyr::transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, MedianBased_RC, WMedianBased_RC, ModeBased_RC, WModeBased_RC) %>%
      tidyr::pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, MedianBased_RC, WMedianBased_RC, ModeBased_RC, WModeBased_RC)) %>%
      mutate(name = lvls_reorder(as_factor(name), idx = match(my_breaks10, levels(as_factor(name))))) -> res_tidy
    ggplot(res_tidy, aes(x=name, y=value, color=name, fill=name)) +
      geom_boxplot(size=1, lty=my_lines10) +
      scale_fill_manual(values = alpha(my_colors10, 0.2)) +
      geom_hline(yintercept = 0.2, color="black", lty=2) +
      labs(title = "",
           x = NULL,
           y = "Estimate") +
      scale_x_discrete(breaks = my_breaks10,
                       labels = my_labels10) +
      scale_color_manual(values= my_colors10,
                         breaks = my_breaks10,
                         labels = my_labels10) +
      theme(legend.position="none",
            axis.text.y = element_text(size=12, hjust=1, color = "black"),
            axis.text.x = element_text(size=15, angle = 45, hjust=1, color = "black"),
            axis.title.y = element_text(size=15, color = "black")) -> fig

    ggsave(paste0("res_sim_", scenario, ".pdf"), fig, height = 19.8, width = 30, units = "cm", dpi = 320)
    ggsave(paste0("res_sim_", scenario, ".jpg"), fig, height = 19.8, width = 30, units = "cm", dpi = 320)

  } else {

    res %>%
      dplyr::transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC) %>%
      tidyr::pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) %>%
      mutate(name = lvls_reorder(as_factor(name), idx = match(my_breaks6, levels(as_factor(name))))) -> res_tidy
    ggplot(res_tidy, aes(x=name, y=value, color=name, fill=name)) +
      geom_boxplot(size=1) +
      scale_fill_manual(values = alpha(my_colors6, 0.2)) +
      geom_hline(yintercept = 0.2, color="black", lty=2) +
      labs(title = "",
           x = NULL,
           y = "Estimate") +
      scale_x_discrete(breaks = my_breaks6,
                       labels = my_labels6, ) +
      scale_color_manual(values= my_colors6,
                         breaks = my_breaks6,
                         labels = my_labels6) +
      theme(legend.position="none",
            axis.text.y = element_text(size=12, hjust=1, color = "black"),
            axis.text.x = element_text(size=15, angle = 45, hjust=1, color = "black"),
            axis.title.y = element_text(size=15, color = "black")) -> fig

    ggsave(paste0("res_sim_", scenario, ".pdf"), fig, height = 18, width = 20, units = "cm", dpi = 320)
    ggsave(paste0("res_sim_", scenario, ".jpg"), fig, height = 18, width = 20, units = "cm", dpi = 320)

  }


  return(fig)

}


make_figure("balancedpleio")
make_figure("correlatedpleio_strong", TRUE)
make_figure("correlatedpleio_weak", TRUE)

# group pleiotropy figures
make_figure_2pleio <- function(){
  fig1 = make_figure("correlatedpleio_weak", TRUE)
  fig2 = make_figure("correlatedpleio_strong", TRUE)

  ## common figure pleio
  ggpubr::ggarrange(fig1 + theme(axis.text.x = element_blank()), NULL,
                    fig2,labels = c("A) weak pleiotropy", "", "B) strong pleiotropy"),
                    heights = c(0.372, 0.04, 0.588), nrow=3) -> figcommon_pleio

  ggsave("res_sim_correlatedpleio.pdf", figcommon_pleio, height = 33, width = 25, units = "cm", dpi = 320)
  ggsave("res_sim_correlatedpleio.png", figcommon_pleio, height = 33, width = 25, units = "cm", dpi = 320)

}

make_figure_2pleio()


# print figure and create pdf + png files for all thresholds
make_figure_all <- function(scenario, pleiotropy = F){

  all_res = readRDS(paste0("res_varyT_", scenario, ".RDS"))

  if(pleiotropy){

    all_res %>%
      transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, MedianBased_RC, WMedianBased_RC, ModeBased_RC, WModeBased_RC,
                th = threshold) %>%
      group_by(th) %>%
      summarize(IVW = mean(IVW),
                IVW_3samples = mean(IVW_3samples),
                RAPS_3samples = mean(RAPS_3samples),
                RAPS_ZP = mean(RAPS_ZP),
                RAPS_UMVCUE = mean(RAPS_UMVCUE),
                IVW_RC = mean(IVW_RC),
                MedianBased_RC = mean(MedianBased_RC),
                WMedianBased_RC = mean(WMedianBased_RC),
                ModeBased_RC = mean(ModeBased_RC),
                WModeBased_RC = mean(WModeBased_RC)) %>%
      pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, MedianBased_RC, WMedianBased_RC, ModeBased_RC, WModeBased_RC)) %>%
      mutate(name = lvls_reorder(as_factor(name), idx = match(my_breaks10, levels(as_factor(name))))) -> all_res_tidy
    ggplot(all_res_tidy, aes(x=th, y=value, color=name, lty=name)) +
      geom_line(size=1) +
      geom_hline(yintercept = 0.2, color="black", lty=2) +
      labs(title = "",
           x = "Threshold",
           y = "Estimate",
           color = "",
           linetype =  "") +
      scale_linetype_manual(values=my_lines10,breaks = my_breaks10,
                            labels = my_labels10) +
      scale_color_manual(values= my_colors10,
                         breaks = my_breaks10,
                         labels = my_labels10) +
      theme(axis.title.x = element_text(size=15, color = "black"),
            axis.title.y = element_text(size=15, color = "black"),
            axis.text.x = element_text(size=12, hjust=1, color = "black"),
            axis.text.y = element_text(size=12, hjust=1, color = "black"),
            legend.title  = element_text(size=15, color = "black"),
            legend.text = element_text(size=13, color = "black"),
            legend.position = "right",
            legend.justification='right',
            legend.direction='vertical',
            legend.key.width = unit(1, 'cm')) +
      guides(color=guide_legend(nrow=10,byrow=TRUE),
             linetype=guide_legend(nrow=10,byrow=TRUE)) -> fig_all

    ggsave(paste0("res_varyT_", scenario, ".pdf"), fig_all, height = 14, width = 32, units = "cm", dpi = 320)
    ggsave(paste0("res_varyT_", scenario, ".jpg"), fig_all, height = 14, width = 32, units = "cm", dpi = 320)


  } else{

    all_res %>%
      transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC,
                th = threshold) %>%
      group_by(th) %>%
      summarize(IVW = mean(IVW, na.rm=T),
                IVW_3samples = mean(IVW_3samples, na.rm=T),
                RAPS_3samples = mean(RAPS_3samples, na.rm=T),
                RAPS_ZP = mean(RAPS_ZP, na.rm=T),
                RAPS_UMVCUE = mean(RAPS_UMVCUE, na.rm=T),
                IVW_RC = mean(IVW_RC, na.rm=T)) %>%
      pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) %>%
      mutate(name = lvls_reorder(as_factor(name), idx = match(my_breaks6, levels(as_factor(name))))) -> all_res_tidy
    ggplot(all_res_tidy, aes(x=th, y=value, color=name)) +
      geom_line(size=1) +
      geom_hline(yintercept = 0.2, color="black", lty=2) +
      labs(title = "",
           x = "Threshold",
           y = "Estimate",
           color = "") +
      scale_color_manual(values=my_colors6,
                         breaks = my_breaks6,
                         labels =my_labels6) +
      theme(axis.title.x = element_text(size=15, color = "black"),
            axis.title.y = element_text(size=15, color = "black"),
            axis.text.x = element_text(size=12, hjust=1, color = "black"),
            axis.text.y = element_text(size=12, hjust=1, color = "black"),
            legend.title  = element_text(size=15, color = "black"),
            legend.text = element_text(size=13, color = "black"),
            legend.position = "right",
            legend.justification='right',
            legend.direction='vertical',
            legend.key.width = unit(1, 'cm')) +
      guides(color=guide_legend(nrow=6,byrow=TRUE))-> fig_all

    ggsave(paste0("res_varyT_", scenario, ".pdf"), fig_all, height = 14, width = 29.7, units = "cm", dpi = 320)
    ggsave(paste0("res_varyT_", scenario, ".jpg"), fig_all, height = 14, width = 29.7, units = "cm", dpi = 320)

  }

  return(fig_all)
}

make_figure_all("balancedpleio")
make_figure_all("correlatedpleio_strong", TRUE)
make_figure_all("correlatedpleio_weak", TRUE)

make_figure_all_rmse <- function(scenario, pleiotropy = F){

  all_res = readRDS(paste0("res_varyT_", scenario, ".RDS"))

  if(pleiotropy){

    all_res %>%
      transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, MedianBased_RC, WMedianBased_RC, ModeBased_RC, WModeBased_RC,
                th = threshold) %>%
      group_by(th) %>%
      summarize(IVW = sqrt(mean((IVW - 0.2)^2, na.rm=T)),
                IVW_3samples = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
                RAPS_3samples = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
                RAPS_ZP = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
                RAPS_UMVCUE = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
                IVW_RC = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T)),
                MedianBased_RC = sqrt(mean((MedianBased_RC - 0.2)^2, na.rm=T)),
                WMedianBased_RC = sqrt(mean((WMedianBased_RC - 0.2)^2, na.rm=T)),
                ModeBased_RC = sqrt(mean((ModeBased_RC - 0.2)^2, na.rm=T)),
                WModeBased_RC = sqrt(mean((WModeBased_RC - 0.2)^2, na.rm=T))) %>%
      pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, MedianBased_RC, WMedianBased_RC, ModeBased_RC, WModeBased_RC)) %>%
      mutate(name = lvls_reorder(as_factor(name), idx = match(my_breaks10, levels(as_factor(name))))) -> all_res_tidy
    ggplot(all_res_tidy, aes(x=th, y=value, color=name, lty=name)) +
      geom_line(size=1) +
      labs(title = "",
           x = "Threshold",
           y = "RMSE",
           color = "",
           linetype =  "") +
      scale_linetype_manual(values=my_lines10,breaks = my_breaks10,
                            labels = my_labels10) +
      scale_color_manual(values= my_colors10,
                         breaks = my_breaks10,
                         labels = my_labels10) +
      theme(axis.title.x = element_text(size=15, color = "black"),
            axis.title.y = element_text(size=15, color = "black"),
            axis.text.x = element_text(size=12, hjust=1, color = "black"),
            axis.text.y = element_text(size=12, hjust=1, color = "black"),
            legend.title  = element_text(size=15, color = "black"),
            legend.text = element_text(size=13, color = "black"),
            legend.position = "right",
            legend.justification='right',
            legend.direction='vertical',
            legend.key.width = unit(1, 'cm')) +
      guides(color=guide_legend(nrow=10,byrow=TRUE),
             linetype=guide_legend(nrow=10,byrow=TRUE)) -> fig_all

    ggsave(paste0("res_varyT_", scenario, "_rmse.pdf"), fig_all, height = 14, width = 32, units = "cm", dpi = 320)
    ggsave(paste0("res_varyT_", scenario, "_rmse.jpg"), fig_all, height = 14, width = 32, units = "cm", dpi = 320)


  } else{

    all_res %>%
      transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC,
                th = threshold) %>%
      group_by(th) %>%
      summarize(IVW = sqrt(mean((IVW - 0.2)^2, na.rm=T)),
                IVW_3samples = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
                RAPS_3samples = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
                RAPS_ZP = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
                RAPS_UMVCUE = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
                IVW_RC = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T))) %>%
      pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) %>%
      mutate(name = lvls_reorder(as_factor(name), idx = match(my_breaks6, levels(as_factor(name))))) -> all_res_tidy
    ggplot(all_res_tidy, aes(x=th, y=value, color=name)) +
      geom_line(size=1) +
      labs(title = "",
           x = "Threshold",
           y = "RMSE",
           color = "") +
      scale_color_manual(values=my_colors6,
                         breaks = my_breaks6,
                         labels =my_labels6) +
      theme(axis.title.x = element_text(size=15, color = "black"),
            axis.title.y = element_text(size=15, color = "black"),
            axis.text.x = element_text(size=12, hjust=1, color = "black"),
            axis.text.y = element_text(size=12, hjust=1, color = "black"),
            legend.title  = element_text(size=15, color = "black"),
            legend.text = element_text(size=13, color = "black"),
            legend.position = "right",
            legend.justification='right',
            legend.direction='vertical',
            legend.key.width = unit(1, 'cm')) +
      guides(color=guide_legend(nrow=6,byrow=TRUE))-> fig_all

    ggsave(paste0("res_varyT_", scenario, "_rmse.pdf"), fig_all, height = 14, width = 29.7, units = "cm", dpi = 320)
    ggsave(paste0("res_varyT_", scenario, "_rmse.jpg"), fig_all, height = 14, width = 29.7, units = "cm", dpi = 320)

  }

  return(fig_all)
}

make_figure_all_rmse("balancedpleio")
make_figure_all_rmse("correlatedpleio_strong", TRUE)
make_figure_all_rmse("correlatedpleio_weak", TRUE)

# group pleiotropy figures
make_figure_all_2pleio <- function(){
  fig1 = make_figure_all("correlatedpleio_weak", TRUE)
  fig2 = make_figure_all("correlatedpleio_strong", TRUE)

  ggpubr::ggarrange(plotlist = list(fig1, NULL, NULL, NULL, fig2, NULL), nrow=3, ncol=2,
                    heights = c(0.48, 0.04, 0.48),
                    widths = c(0.97, 0.03),
                    labels = c("A) weak pleiotropy", "", "", "", "B) strong pleiotropy", ""),
                    common.legend = T, legend="right") +
    theme(plot.margin = margin(5.5,5.5,5.5,5.5, "points")) -> figcommon_pleio_T

  ggsave("res_varyT_correlatedpleio.pdf", figcommon_pleio_T, height = 29.7, width = 32, units = "cm", dpi = 320)
  ggsave("res_varyT_correlatedpleio.png", figcommon_pleio_T, height = 29.7, width = 32, units = "cm", dpi = 320)

}

make_figure_all_2pleio()


# print figure and create pdf + png files for all proportions of discovery/replication
make_figure_allN <- function(scenario){

  all_res_N = readRDS(paste0("res_varyN_", scenario, ".RDS"))


  all_res_N %>%
    transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, prop) %>%
    group_by(prop) %>%
    summarize(IVW = mean(IVW),
              IVW_3samples = mean(IVW_3samples),
              RAPS_3samples = mean(RAPS_3samples),
              RAPS_ZP = mean(RAPS_ZP),
              RAPS_UMVCUE = mean(RAPS_UMVCUE),
              IVW_RC = mean(IVW_RC)) %>%
    pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) -> all_res_N_tidy
  ggplot(all_res_N_tidy, aes(x=prop, y=value, color=name)) +
    geom_line(size=1) +
    geom_hline(yintercept = 0.2, color="black", lty=2) +
    labs(title = "",
         x = "Proportion of the exposure sample used as discovery",
         y = "Estimate",
         color = "") +
    scale_color_manual(values=my_colors6,
                       breaks = my_breaks6,
                       labels =my_labels6) +
    theme(axis.title.x = element_text(size=15, color = "black"),
          axis.title.y = element_text(size=15, color = "black"),
          axis.text.x = element_text(size=12, color = "black"),
          axis.text.y = element_text(size=12, color = "black"),
          legend.title  = element_text(size=15, color = "black"),
          legend.text = element_text(size=13, color = "black"),
          legend.position = "right",
          legend.justification='right',
          legend.direction='vertical',
          legend.key.width = unit(1, 'cm')) +
    guides(color=guide_legend(nrow=6,byrow=TRUE))-> fig_allN


  ggsave(paste0("res_varyN_", scenario, ".pdf"), fig_allN, height =  14, width = 29.7, units = "cm", dpi = 320)
  ggsave(paste0("res_varyN_", scenario, ".jpg"), fig_allN, height =  14, width = 29.7, units = "cm", dpi = 320)



  return(fig_allN)
}

make_figure_allN("balancedpleio")
#make_figure_allN("balancedpleio_largerN")

make_figure_allN_rmse <- function(scenario){

  all_res_N = readRDS(paste0("res_varyN_", scenario, ".RDS"))


  all_res_N %>%
    transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, prop) %>%
    group_by(prop) %>%
    summarize(IVW = sqrt(mean((IVW - 0.2)^2, na.rm=T)),
              IVW_3samples = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
              RAPS_3samples = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
              RAPS_ZP = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
              RAPS_UMVCUE = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
              IVW_RC = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T))) %>%
    pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) -> all_res_N_tidy
  ggplot(all_res_N_tidy, aes(x=prop, y=value, color=name)) +
    geom_line(size=1) +
    labs(title = "",
         x = "Proportion of the exposure sample used as discovery",
         y = "RMSE",
         color = "") +
    scale_color_manual(values=my_colors6,
                       breaks = my_breaks6,
                       labels =my_labels6) +
    theme(axis.title.x = element_text(size=15, color = "black"),
          axis.title.y = element_text(size=15, color = "black"),
          axis.text.x = element_text(size=12, color = "black"),
          axis.text.y = element_text(size=12, color = "black"),
          legend.title  = element_text(size=15, color = "black"),
          legend.text = element_text(size=13, color = "black"),
          legend.position = "right",
          legend.justification='right',
          legend.direction='vertical',
          legend.key.width = unit(1, 'cm')) +
    guides(color=guide_legend(nrow=6,byrow=TRUE))-> fig_allN


  ggsave(paste0("res_varyN_", scenario, "_rmse.pdf"), fig_allN, height =  14, width = 29.7, units = "cm", dpi = 320)
  ggsave(paste0("res_varyN_", scenario, "_rmse.jpg"), fig_allN, height =  14, width = 29.7, units = "cm", dpi = 320)



  return(fig_allN)
}

make_figure_allN_rmse("balancedpleio")
#make_figure_allN_rmse("balancedpleio_largerN")


# print figure and create pdf + png files for all (total) exposure sample sizes
make_figure_Ntotexp <- function(scenario){

  all_res_N = readRDS(paste0("res_varyNtotexp_", scenario, ".RDS"))


  all_res_N %>%
    transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, tot_N) %>%
    group_by(tot_N) %>%
    summarize(IVW = mean(IVW),
              IVW_3samples = mean(IVW_3samples),
              RAPS_3samples = mean(RAPS_3samples),
              RAPS_ZP = mean(RAPS_ZP),
              RAPS_UMVCUE = mean(RAPS_UMVCUE),
              IVW_RC = mean(IVW_RC)) %>%
    pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) -> all_res_N_tidy
  ggplot(all_res_N_tidy, aes(x=tot_N, y=value, color=name)) +
    geom_line(size=1) +
    geom_hline(yintercept = 0.2, color="black", lty=2) +
    labs(title = "",
         x = "Exposure (total: discovery+replication) sample size",
         y = "Estimate",
         color = "") +
    scale_color_manual(values=my_colors6,
                       breaks = my_breaks6,
                       labels =my_labels6) +
    theme(axis.title.x = element_text(size=15, color = "black"),
          axis.title.y = element_text(size=15, color = "black"),
          axis.text.x = element_text(size=12, color = "black"),
          axis.text.y = element_text(size=12, color = "black"),
          legend.title  = element_text(size=15, color = "black"),
          legend.text = element_text(size=13, color = "black"),
          legend.position = "right",
          legend.justification='right',
          legend.direction='vertical',
          legend.key.width = unit(1, 'cm')) +
    guides(color=guide_legend(nrow=6,byrow=TRUE))-> fig_allN


  ggsave(paste0("res_varyNtotexp_", scenario, ".pdf"), fig_allN, height =  14, width = 29.7, units = "cm", dpi = 320)
  ggsave(paste0("res_varyNtotexp_", scenario, ".jpg"), fig_allN, height = 14, width = 29.7, units = "cm", dpi = 320)



  return(fig_allN)
}

make_figure_Ntotexp("balancedpleio")

make_figure_Ntotexp_rmse <- function(scenario){

  all_res_N = readRDS(paste0("res_varyNtotexp_", scenario, ".RDS"))


  all_res_N %>%
    transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, tot_N) %>%
    group_by(tot_N) %>%
    summarize(IVW = sqrt(mean((IVW - 0.2)^2, na.rm=T)),
              IVW_3samples = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
              RAPS_3samples = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
              RAPS_ZP = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
              RAPS_UMVCUE = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
              IVW_RC = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T))) %>%
    pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) -> all_res_N_tidy
  ggplot(all_res_N_tidy, aes(x=tot_N, y=value, color=name)) +
    geom_line(size=1) +
    labs(title = "",
         x = "Exposure (total: discovery+replication) sample size",
         y = "RMSE",
         color = "") +
    scale_color_manual(values=my_colors6,
                       breaks = my_breaks6,
                       labels =my_labels6) +
    theme(axis.title.x = element_text(size=15, color = "black"),
          axis.title.y = element_text(size=15, color = "black"),
          axis.text.x = element_text(size=12, color = "black"),
          axis.text.y = element_text(size=12, color = "black"),
          legend.title  = element_text(size=15, color = "black"),
          legend.text = element_text(size=13, color = "black"),
          legend.position = "right",
          legend.justification='right',
          legend.direction='vertical',
          legend.key.width = unit(1, 'cm')) +
    guides(color=guide_legend(nrow=6,byrow=TRUE))-> fig_allN


  ggsave(paste0("res_varyNtotexp_", scenario, "_rmse.pdf"), fig_allN, height =  14, width = 29.7, units = "cm", dpi = 320)
  ggsave(paste0("res_varyNtotexp_", scenario, "_rmse.jpg"), fig_allN, height = 14, width = 29.7, units = "cm", dpi = 320)



  return(fig_allN)
}

make_figure_Ntotexp_rmse("balancedpleio")


# print figure and create pdf + png files for all outcome sample sizes
make_figure_Ny <- function(scenario){

  all_res_N = readRDS(paste0("res_varyNy_", scenario, ".RDS"))


  all_res_N %>%
    transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, Ny) %>%
    group_by(Ny) %>%
    summarize(IVW = mean(IVW),
              IVW_3samples = mean(IVW_3samples),
              RAPS_3samples = mean(RAPS_3samples),
              RAPS_ZP = mean(RAPS_ZP),
              RAPS_UMVCUE = mean(RAPS_UMVCUE),
              IVW_RC = mean(IVW_RC)) %>%
    pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) -> all_res_N_tidy
  ggplot(all_res_N_tidy, aes(x=Ny, y=value, color=name)) +
    geom_line(size=1) +
    geom_hline(yintercept = 0.2, color="black", lty=2) +
    labs(title = "",
         x = "Outcome sample size",
         y = "Estimate",
         color = "") +
    scale_color_manual(values=my_colors6,
                       breaks = my_breaks6,
                       labels =my_labels6) +
    scale_x_continuous(limits = c(20000, 260000), breaks = seq(20000, 260000, 60000)) +
    theme(axis.title.x = element_text(size=15, color = "black"),
          axis.title.y = element_text(size=15, color = "black"),
          axis.text.x = element_text(size=12, color = "black"),
          axis.text.y = element_text(size=12, color = "black"),
          legend.title  = element_text(size=15, color = "black"),
          legend.text = element_text(size=13, color = "black"),
          legend.position = "right",
          legend.justification='right',
          legend.direction='vertical',
          legend.key.width = unit(1, 'cm')) +
    guides(color=guide_legend(nrow=6,byrow=TRUE))-> fig_allN


  ggsave(paste0("res_varyNy_", scenario, ".pdf"), fig_allN, height = 14, width = 29.7, units = "cm", dpi = 320)
  ggsave(paste0("res_varyNy_", scenario, ".jpg"), fig_allN, height = 14, width = 29.7, units = "cm", dpi = 320)



  return(fig_allN)
}

make_figure_Ny("balancedpleio")

make_figure_Ny_rmse <- function(scenario){

  all_res_N = readRDS(paste0("res_varyNy_", scenario, ".RDS"))


  all_res_N %>%
    transmute(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC, Ny) %>%
    group_by(Ny) %>%
    summarize(IVW = sqrt(mean((IVW - 0.2)^2, na.rm=T)),
              IVW_3samples = sqrt(mean((IVW_3samples - 0.2)^2, na.rm=T)),
              RAPS_3samples = sqrt(mean((RAPS_3samples - 0.2)^2, na.rm=T)),
              RAPS_ZP = sqrt(mean((RAPS_ZP - 0.2)^2, na.rm=T)),
              RAPS_UMVCUE = sqrt(mean((RAPS_UMVCUE - 0.2)^2, na.rm=T)),
              IVW_RC = sqrt(mean((IVW_RC - 0.2)^2, na.rm=T))) %>%
    pivot_longer(c(IVW, IVW_3samples, RAPS_3samples, RAPS_ZP, RAPS_UMVCUE, IVW_RC)) -> all_res_N_tidy
  ggplot(all_res_N_tidy, aes(x=Ny, y=value, color=name)) +
    geom_line(size=1) +
    labs(title = "",
         x = "Outcome sample size",
         y = "RMSE",
         color = "") +
    scale_color_manual(values=my_colors6,
                       breaks = my_breaks6,
                       labels =my_labels6) +
    scale_x_continuous(limits = c(20000, 260000), breaks = seq(20000, 260000, 60000)) +
    theme(axis.title.x = element_text(size=15, color = "black"),
          axis.title.y = element_text(size=15, color = "black"),
          axis.text.x = element_text(size=12, color = "black"),
          axis.text.y = element_text(size=12, color = "black"),
          legend.title  = element_text(size=15, color = "black"),
          legend.text = element_text(size=13, color = "black"),
          legend.position = "right",
          legend.justification='right',
          legend.direction='vertical',
          legend.key.width = unit(1, 'cm')) +
    guides(color=guide_legend(nrow=6,byrow=TRUE))-> fig_allN


  ggsave(paste0("res_varyNy_", scenario, "_rmse.pdf"), fig_allN, height = 14, width = 29.7, units = "cm", dpi = 320)
  ggsave(paste0("res_varyNy_", scenario, "_rmse.jpg"), fig_allN, height = 14, width = 29.7, units = "cm", dpi = 320)



  return(fig_allN)
}

make_figure_Ny_rmse("balancedpleio")


# 1x2 figure estimate + RMSE
# param should be "varyT", or "varyN", or "varyNtotexp", or "varyNy"
make_panel_rmse  <- function(scenario, param){
  if(param == "varyT"){
    fig1 = make_figure_all(scenario) +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("") +
      theme(legend.direction='horizontal',
            legend.position = "bottom",
            legend.justification = "center",
            legend.spacing.x = unit(0.25, 'cm'),
            legend.key.width = unit(1, 'cm'),
            legend.key.size = unit(1.5, 'lines'),
            legend.text=element_text(margin=margin(r=0.5,unit="inch"))) +
      scale_y_continuous(
        labels = scales::number_format(accuracy = 0.01))
    fig2 = make_figure_all_rmse(scenario) +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("") +
      theme(legend.direction='horizontal',
            legend.position = "bottom",
            legend.justification = "center",
            legend.spacing.x = unit(0.25, 'cm'),
            legend.key.width = unit(1, 'cm'),
            legend.key.size = unit(1.5, 'lines'),
            legend.text=element_text(margin=margin(r=0.5,unit="inch"))) +
      scale_y_continuous(
              labels = scales::number_format(accuracy = 0.01))
  }   else if (param == "varyN") {
    fig1 = make_figure_allN(scenario) +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("") +
      theme(legend.direction='horizontal',
            legend.position = "bottom",
            legend.justification = "center",
            legend.spacing.x = unit(0.25, 'cm'),
            legend.key.width = unit(1, 'cm'),
            legend.key.size = unit(1.5, 'lines'),
            legend.text=element_text(margin=margin(r=0.5,unit="inch"))) +
      scale_y_continuous(
        labels = scales::number_format(accuracy = 0.01))
    fig2 = make_figure_allN_rmse(scenario) +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("") +
      theme(legend.direction='horizontal',
            legend.position = "bottom",
            legend.justification = "center",
            legend.spacing.x = unit(0.25, 'cm'),
            legend.key.width = unit(1, 'cm'),
            legend.key.size = unit(1.5, 'lines'),
            legend.text=element_text(margin=margin(r=0.5,unit="inch"))) +
      scale_y_continuous(
        labels = scales::number_format(accuracy = 0.01))
  } else if(param == "varyNtotexp") {
    fig1 = make_figure_Ntotexp(scenario) +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("") +
      theme(legend.direction='horizontal',
            legend.position = "bottom",
            legend.justification = "center",
            legend.spacing.x = unit(0.25, 'cm'),
            legend.key.width = unit(1, 'cm'),
            legend.key.size = unit(1.5, 'lines'),
            legend.text=element_text(margin=margin(r=0.5,unit="inch"))) +
      scale_y_continuous(
        labels = scales::number_format(accuracy = 0.01))
    fig2 = make_figure_Ntotexp_rmse(scenario) +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("") +
      theme(legend.direction='horizontal',
            legend.position = "bottom",
            legend.justification = "center",
            legend.spacing.x = unit(0.25, 'cm'),
            legend.key.width = unit(1, 'cm'),
            legend.key.size = unit(1.5, 'lines'),
            legend.text=element_text(margin=margin(r=0.5,unit="inch"))) +
      scale_y_continuous(
        labels = scales::number_format(accuracy = 0.01))
  } else if(param == "varyNy"){
    fig1 = make_figure_Ny(scenario) +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("") +
      theme(legend.direction='horizontal',
            legend.position = "bottom",
            legend.justification = "center",
            legend.spacing.x = unit(0.25, 'cm'),
            legend.key.width = unit(1, 'cm'),
            legend.key.size = unit(1.5, 'lines'),
            legend.text=element_text(margin=margin(r=0.5,unit="inch"))) +
      scale_y_continuous(
        labels = scales::number_format(accuracy = 0.01), breaks = scales::breaks_pretty(3))
    fig2 = make_figure_Ny_rmse(scenario) +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      ylab("") +
      theme(legend.direction='horizontal',
            legend.position = "bottom",
            legend.justification = "center",
            legend.spacing.x = unit(0.25, 'cm'),
            legend.key.width = unit(1, 'cm'),
            legend.key.size = unit(1.5, 'lines'),
            legend.text=element_text(margin=margin(r=0.5,unit="inch"))) +
      scale_y_continuous(
        labels = scales::number_format(accuracy = 0.01))
  } else {
    stop("wrong parameter")
  }



  ggpubr::ggarrange(plotlist = list(fig1, fig2), nrow=1, ncol=2,
                    widths = c(0.5, 0.5),
                    labels = c("A) Estimate", "B) RMSE"), label.x = c(0, 0.02),
                    common.legend = T, legend="bottom") +
    theme(plot.margin = margin(5.5,5.5,5.5,5.5, "points")) -> figcommon

  ggsave(paste0("res_", param, "_", scenario, "_panel.pdf"), figcommon, height = 16, width = 35, units = "cm", dpi = 320)
  ggsave(paste0("res_", param, "_", scenario, "_panel.jpg"), figcommon, height = 16, width = 35, units = "cm", dpi = 320)

}
make_panel_rmse("balancedpleio", "varyT")
make_panel_rmse("balancedpleio", "varyN")
# make_panel_rmse("balancedpleio_largerN", "varyN")
make_panel_rmse("balancedpleio", "varyNtotexp")
make_panel_rmse("balancedpleio", "varyNy")
