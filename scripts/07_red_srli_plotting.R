
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 22 June 2026


# Intro ------------------------------------------------------------------------
# This script is used to create the main figures of the manuscript based on the outputs of the other scripts. 

# Libraries --------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(sf)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(grid) 
library(gridExtra)
library(cowplot)
library(ggplotify)




plot_prop_sample_f <- function(continent_name, 
                               null = null_props, 
                               null_stats = null_prop_stats) {
  
  pal <- c(
    "#ff6347",
    "#fcc200"
   )
  
  
  null <- null %>% 
    filter(LEVEL1_NAM == continent_name) %>% 
    mutate(across(contains("prop"),  ~ .x *100))
  
  null_stats <- null_stats %>% 
    filter(LEVEL1_NAM == continent_name) %>% 
    mutate(across(where(is.numeric),  ~ .x *100))
    
  
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, 
               colour = "grey") +
    geom_jitter(data = null, 
                aes(y = (prop_diff), x = dataset, col = dataset, 
                    fill = dataset), alpha = 0.05, width = 0.2, show.legend = F) +
    geom_errorbar(data = null_stats, 
                  mapping = aes(x = dataset, group = dataset,
                                ymin= q_lower  , ymax= q_upper),col = "black",  
                  width = 0.1) +
    geom_point(data = null_stats, mapping = aes(x = dataset, group = dataset,
                                                y= mean)) +
    theme_bw(base_size = 10) +
    scale_color_manual(values =pal)  +
    scale_fill_manual(values =pal)  +
    scale_y_continuous(guide = guide_axis(angle = 0), limits = c(-10, 10)) +
    labs(x = NULL, y = "% difference to representative sample") +
    #coord_flip() + 
    theme(strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text = element_text(face = "bold", size = 10), 
          legend.position = "top", 
          legend.box = "horizontal") + 
    coord_flip()
}
plot_growth_sample_f <- function(continent_name, 
                                 growth_stats = null_growth, 
                                 set = "redlist") {
  cols <- c(nationalparkcolors::park_palette("Saguaro")[1:5], "#5b2635")
  
  
  growth_stats <- null_growth %>% 
    filter(LEVEL1_NAM == continent_name) %>% 
    filter(dataset == set)
  
  ggplot(growth_stats, aes(x = growth_form, y = dif_total)) +

    geom_bar(aes(fill = growth_form),  stat = "identity") + 
    scale_fill_manual(values = c(cols)) +
    theme_bw() +
    #scale_y_continuous(breaks = c(seq(-9000, 9000, 3000)), limits = c(-9000,9000)) +
    #scale_y_continuous(breaks = c(0,seq(-2, 2, 1)), limits = c(-2.5,2.5)) +
   scale_y_continuous(breaks = c(0,seq(0, 40000, 10000)), limits = c(0,40000)) +
    
    labs(x = NULL, y = "Number of assessments missing", fill = "Growth form") +
    coord_flip() + 
    theme(          legend.position = "top",
                    legend.box = "horizontal") +
    guides(fill = guide_legend(nrow = 1))
  
  # ggplot(growth_stats, aes(x = growth_form, y = dif)) +
  #   #geom_segment(yend = redlist_country_growths$dif, y = 0, x= redlist_country_growths$growth_form, 
  #   #             col ="grey") +
  #   #geom_point(aes(col = growth_form), size = 2.5) + 
  #   geom_bar(aes(fill = growth_form),  stat = "identity") + 
  #   
  #   geom_errorbar(aes(ymin= dif_max, ymax=dif_min), width = 0.25) +
  #   scale_fill_manual(values = c(cols)) +
  #   #paletteer::scale_fill_paletteer_d("nationalparkcolors::Badlands") +
  #   #paletteer::scale_fill_paletteer_d("nationalparkcolors::Saguaro") +
  #   
  #   #geom_hline(yintercept = 0, lty = "dashed")  +
  #   theme_bw() +
  #   scale_y_continuous(breaks = c(seq(-9000, 9000, 3000)), limits = c(-9000,9000)) +
  #   #scale_y_continuous(breaks = c(0,seq(-2, 2, 1)), limits = c(-2.5,2.5)) +
  #   
  #   labs(x = NULL, y = "Difference to random", fill = "Growth form") +
  #   coord_flip() + 
  #   theme(          legend.position = "top",
  #                   legend.box = "horizontal") +
  #   guides(fill = guide_legend(nrow = 1))
  
}

plot_threat_sample_f <- function(continent_name, 
                                 threat_stats = null_threat
) {
  
  threat_stats_c <- threat_stats %>% 
    filter(LEVEL1_NAM == continent_name) %>% 
    mutate(across(where(is.numeric),  ~ .x *100))
  
  pal <- c(
    "darkred",
    
    "goldenrod3")
  
  ggplot(threat_stats_c) +
    geom_bar(
      aes(y = (mean), x = dataset, col = dataset,
          fill = dataset), alpha = 0.6, width = 0.3, stat = "identity", 
      show.legend = F) +
    # geom_jitter(data = threat_null, 
    #             aes(y = (prop_threat_cor), x = dataset, col = dataset, 
    #                 fill = dataset), alpha = 0.05, width = 0.1) +
    geom_errorbar( mapping = aes(x = dataset, 
                                 ymin= q_lower  , ymax= q_upper), width = 0.1) +
    theme_bw(base_size = 10) +
    scale_color_manual(values =pal)  +
    scale_fill_manual(values =pal)  +
    scale_y_continuous(guide = guide_axis(angle = 0), limits = c(0,100)) +
    labs(x = NULL, y = "% estimated threatened spp.") +
    #coord_flip() + 
    theme(strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text = element_text(face = "bold", size = 6), 
          legend.position = "top", 
          legend.box = "horizontal") + 
    coord_flip()
}

plot_map_f <- function(continent_name, pal_val = col_pick) {
  
  tdwg_3_map <-  tdwg_3 %>%  
    filter(LEVEL1_NAM == continent_name) 
  
  if(continent_name == "PACIFIC") {
    
    tdwg_3_p <-  tdwg_3 %>% 
      st_drop_geometry() %>% 
      st_as_sf(coords = c("lon","lat"), crs = st_crs(tdwg_3)) %>% 
      st_transform(st_crs(eqearth_crs)) %>%  
      filter(LEVEL1_NAM == continent_name) 
      
    
    ggplot(tdwg_3_p) +
      geom_sf(data = border, fill = "white", color = "black", size = 0.2) +
     geom_sf(data = tdwg_3, fill = "white", color = "black", size = 0.1) +
      geom_sf(aes(col = LEVEL1_NAM), show.legend = F,
              size = 1) +
      # paletteer::scale_fill_paletteer_d("colorblindr::OkabeIto_black") +
      scale_color_manual(values = pal_val) +
      #facet_grid(~LEVEL1_NAM, switch = "y") +
      coord_sf(crs = eqearth_crs, expand = T) +
      theme_minimal(base_size = 12) +
      labs(
        title = continent_name
      )  +
      theme(
        plot.title.position = "plot",
        plot.title = element_text(face = "bold",size = 12, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        # strip.switch.pad.grid = unit(0.4, "cm") 
      ) 
    
  } else {
    
    ggplot(tdwg_3_map) +
      geom_sf(data = border, fill = "white", color = "black", size = 0.2) +
      geom_sf(data = tdwg_3, fill = "white", color = "black", size = 0.1) +
      geom_sf(aes(fill = LEVEL1_NAM), show.legend = F, color = "transparent", size = 0.1) +
      # paletteer::scale_fill_paletteer_d("colorblindr::OkabeIto_black") +
      scale_fill_manual(values = pal_val) +
      #facet_grid(~LEVEL1_NAM, switch = "y") +
      coord_sf(crs = eqearth_crs, expand = T) +
      theme_minimal(base_size = 12) +
      labs(
        title = continent_name
      )  +
      theme(
        plot.title.position = "plot",
        plot.title = element_text(face = "bold",size = 12, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        strip.text = element_blank(),
        # strip.switch.pad.grid = unit(0.4, "cm") 
      ) 
  }
  
}

plot_curves_f <- function(continent_name, 
                       curve_data = acc_curves_cont, 
                       pal_val = col_pick) {
  
  curves_div <-  curve_data %>%  
    filter(LEVEL1_NAM == continent_name) %>% 
    mutate(across(contains("prop"),  ~ .x *100))

  
  ggplot() +
    
    geom_line(data = curve_data %>% ungroup() %>% 
                dplyr::rename(cont = LEVEL1_NAM) %>% 
                mutate(across(contains("prop"),  ~ .x *100)), 
              aes(x = yearPublished, y = prop_sum, group = cont), lwd = 0.3, 
              col = "grey", show.legend = F) +
    geom_line(data = curves_div, aes(x = yearPublished, y = prop_sum, col = LEVEL1_NAM), 
              lwd = 1, show.legend = T) +
    #paletteer::scale_colour_paletteer_d("MoMAColors::Klein") +
    #paletteer::scale_colour_paletteer_d("ggthemes::calc") +
    #paletteer::scale_colour_paletteer_d("ggthemes::Green_Orange_Teal") +
    scale_colour_manual(values = pal_val) +
    scale_y_continuous(guide = guide_axis(angle = 0), limits = c(0,100)) +
    #  paletteer::scale_colour_paletteer_d("khroma::light") +
    #paletteer::scale_colour_paletteer_d("nationalparkcolors::Saguaro") +
    # paletteer::scale_colour_paletteer_d("jcolors::pal4")+
    #geom_point(col = "red", size = 2) +
    theme_bw() +
    labs(x = "Year", y = "% spp. assessed")

}

    


# data -------------------------------------------------------------------------

null_prop_stats <- fread("output/plotting_data/red_srli_null_prop_stats_09_2026.txt") 

null_props <- fread("output/plotting_data/red_srli_null_prop_09_2026.txt")

null_growth <- 
  fread("output/plotting_data/red_srli_null_growth_stats_cont_09_2026.txt") %>% 
  mutate(n = ifelse(is.na(n), 0, n), 
           
           dif_total = value - n, 
           dif_total_prop = (value - n) / value, 
           
           
           dif = n - mean, 
           dif_max = n - q_lower,
           dif_min = n - q_upper,
           
           dif2 =  (n - mean) / mean, 
           dif2_max = (n - q_lower) /  q_lower,
           dif2_min = (n - q_upper) / q_upper,
           
           dif3 =  (n - mean) / sd, 
           dif3_max = (n - q_lower) /  sd,
           dif3_min = (n - q_upper) / sd,
           
  
        dif_prop_r = n / mean, 
         dif_prop_r_min = n / q_lower, 
         dif_prop_r_max = n / q_upper)
 

null_threat <- 
  fread("output/plotting_data/red_srli_threat_stats_cor_fact_09_2026.txt") %>% 
  mutate(dataset = paste0(dataset, "_threat"))

acc_curves_cont <- fread("output/plotting_data/red_continent_acc_curves_09_2026.txt")
redlist_country_growths <- fread("output/plotting_data/red_growths_country_comp_09_2026.txt")


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU")




# --- plotting  ----------------------------------------------------------------
continents <- unique(tdwg_3$LEVEL1_NAM)


eqearth_crs <- "+proj=eqearth"

border <- st_graticule() |>
  st_bbox() |>
  st_as_sfc() |>
  st_transform(3857) |>
  st_segmentize(500000) |>
  st_transform(st_crs(eqearth_crs)) |>
  st_cast("POLYGON")


#conts <- unique(tdwg_3$LEVEL1_NAM)
conts <- c("NORTHERN AMERICA",  "SOUTHERN AMERICA",  
           "EUROPE",  "ASIA-TEMPERATE",  "ASIA-TROPICAL", 
           "AFRICA", 
           "AUSTRALASIA", 
           "PACIFIC", "ANTARCTICA")


pal_val <- c("#FF7F0FFF","#FFB977FF",
             "#ACD98DFF","#39734CFF","#82853BFF",
             "#7c1c05",
             "#FFD94AFF",
             "#3CB7CCFF","#98D9E4FF"
            )

plot_list <- list()

plot_index <- 1

for (i in 1:length(conts)) {
  
  # generate plots
  
  col_pick2 <- pal_val[i]

  m <-plot_map_f(continent_name = conts[i], pal_val = col_pick2) + 
    theme(legend.position = "none",
          aspect.ratio = 1/1.5)
  
  
  ac <- plot_curves_f(continent_name = conts[i], pal_val = col_pick2) +  
    theme(legend.position = "none", 
          aspect.ratio = 1/1.5)
  
  a <-plot_prop_sample_f(continent_name = conts[i]) + 
    theme(legend.position = "none", 
          aspect.ratio = 1/1.5)

  b <- plot_threat_sample_f(continent_name = conts[i])+ 
    theme(legend.position = "none", 
          aspect.ratio = 1/1.5)
  
  c <- plot_growth_sample_f(continent_name = conts[i])+ 
    theme(legend.position = "none", 
          aspect.ratio = 1/1.7)
  
  # d <- plot_pca_f(continent_name = conts[i])+ 
  #   theme(legend.position = "none", 
  #         aspect.ratio = 1/1.5)
  # }
  
  # Save each plot individually into the flat list
  plot_list[[plot_index]] <- m
  plot_index <- plot_index + 1
  
  plot_list[[plot_index]] <- ac
  plot_index <- plot_index + 1
  
  plot_list[[plot_index]] <- c
  plot_index <- plot_index + 1
  
  plot_list[[plot_index]] <- a
  plot_index <- plot_index + 1
  
  plot_list[[plot_index]] <- b
  plot_index <- plot_index + 1
  
  # plot_list[[plot_index]] <- d
  # plot_index <- plot_index + 1
}


labels <- c("(a)","", "", "", "", 
            "(b)","", "", "", "", 
            "(c)","", "", "", "", 
            "(d)", "", "", "", "", 
            "(e)","", "", "", "", 
            "(f)","", "", "", "",  
            "(g)", "", "", "", "", 
            "(h)", "", "", "", "",
            "(i)", "", "", "", ""
            
)
combined_plot <- plot_grid(plotlist = plot_list, nrow = 9, ncol = 5, labels = labels) 

ggsave("09_20206_all_continents_redsrli_pubready.png", combined_plot, width = 14, height = 16, dpi = 900, 
       units = "in")


# --- growth from maps 





get_optimal_breaks <- function(x, min_step = 100) {
  max_val <- max(x, na.rm = TRUE)
  
  # We test step sizes: 100, 200, 300, 400...
  # We want the smallest step where (max_val / step) is between 5 and 7 
  # (Because 5 intervals = 6 breaks, 7 intervals = 8 breaks)
  
  step <- min_step
  found_step <- FALSE
  
  while(!found_step) {
    n_intervals <- max_val / step
    
    # Check if this step gives us 6, 7, or 8 breaks
    if (n_intervals >= 5 && n_intervals <= 7) {
      found_step <- TRUE
    } else if (n_intervals < 5) {
      # If we have too few breaks (e.g., only 4), the step is too large.
      # We stop here anyway to avoid the "1800" problem.
      found_step <- TRUE
    } else {
      # Too many breaks, increase the step size and try again
      step <- step + min_step
    }
  }
  
  # Calculate final sequence
  n_steps <- ceiling(max_val / step)
  seq(0, by = step, length.out = n_steps + 1)
}

plot.growth.map.corr <- function(growth, 
                                 data = redlist_country_growths) {
  
  redlist_country_growths_f <- data %>% 
    filter(growth_form == growth)
  
  extreme_pts <- redlist_country_growths_f %>%
    ungroup() %>% 
    # group_by(bi_class) %>%  
    slice_max(dif, n = 5)
  
  #options(warn = 1)
  extreme_pts <- redlist_country_growths_f %>%
    ungroup() %>% 
    # group_by(bi_class) %>%  
    slice_max(dif *-1, n = 5) %>% 
    rbind(extreme_pts)
  
  extreme_pts <- redlist_country_growths_f %>%
    ungroup() %>% 
    # group_by(bi_class) %>%  
    slice_max(n, n = 5) %>% 
    rbind(extreme_pts)
  
  
  extreme_pts <- redlist_country_growths_f %>%
    ungroup() %>% 
    # group_by(bi_class) %>%  
    slice_max(mean, n = 5) %>% 
    rbind(extreme_pts) %>% 
    filter(!duplicated(LEVEL3_NAM))
  
  max_f <- max(c(redlist_country_growths_f$n, redlist_country_growths_f$q_upper), na.rm = T)
  
  breaks3 <- get_optimal_breaks(max_f)
  
  
  
  continent_dif_plot <- ggplot() +
    geom_point(data = redlist_country_growths_f, aes(x = mean, y = n, col = growth_form), 
               size = 1, 
               show.legend = F) + 
    geom_text_repel(
      data = extreme_pts,
      aes(y = n, x = mean, label = LEVEL3_NAM),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2, 
      min.segment.length = 0,
      show.legend = F
    ) +
    scale_color_manual(values = cols_re_p) +
    geom_errorbar(data = redlist_country_growths_f, 
                  aes(xmin= q_lower, xmax=q_upper, x= mean, y = n, col = growth_form), 
                  width = 0.25, 
                  show.legend = F) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", col = "grey70") +
    #facet_wrap(~growth_form, scales = "free", nrow = 2) +
    scale_y_continuous(breaks = breaks3, limits = c(0, max(breaks3) +1)) +
    scale_x_continuous(breaks = breaks3, limits = c(0, max(breaks3) +1)) +
    labs(y = "Spp. in IUCN Red List", x = "Spp. in random sample", col = "Growth form", 
         caption =  growth_vec[i]) +
    theme(aspect.ratio = 1/1.5, 
          plot.caption = element_text(size = 16, face = "bold"))
  
  
  
  x <- redlist_country_growths_f %>% 
    right_join(tdwg_3, by = "LEVEL3_NAM") %>% 
    st_as_sf()
  
  
  min<- min(redlist_country_growths_f$dif, na.rm = T)
  max <- max(redlist_country_growths_f$dif, na.rm = T)
  
  if(abs(min) > abs(max)) {
    min<- plyr::round_any(min, 100, f = floor)
    
    # width <- abs(round_any((abs(min)) / 3, 100, f = ceiling))
    
    breaks2 <- seq(min, -min, length.out = 5)
    
  } else {
    max<- plyr::round_any(max, 100, f = ceiling)
    
    
    
    breaks2 <- seq(-max, max, length.out = 5)
  }
  
  diff_plot <- ggplot(redlist_country_growths_f %>% 
                        right_join(tdwg_3, by = "LEVEL3_NAM") %>% 
                        st_as_sf()) +
    #geom_sf(data = border, fill = "azure", color = "black", size = 0.2) +
    geom_sf(aes(fill = dif), color = "grey25", size = 0.2) +
    scale_fill_gradient2(  low ="#E65100",  #"#ff6347",
                           mid = "white",
                           high = "#01579B", #"#00dfff",
                           midpoint = 0, 
                           na.value = "grey", 
                           name = "Spp. difference", 
                           breaks = c(breaks2), 
                           limits = c(min(breaks2) -1 , max(breaks2) +1)) +
    # scale_fill_manual(values = pal2, na.value = "grey") +
    #bi_scale_fill(pal = "DkBlue", dim = 2, flip_axes = T) +
    coord_sf(crs = eqearth_crs, expand = T) +
    theme_minimal(base_size = 12) +
    theme(
      # panel.spacing = unit(1, "cm"),  # or "cm", "mm", etc.
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "black"),
      strip.text = element_text(color = "white", face = "bold", size = 6),
      panel.grid.major = element_line(color = "grey85", size = 0.2),
      legend.position = "right",
      legend.box = "vertical",
      axis.text = element_blank(),
      axis.ticks = element_blank(), 
      aspect.ratio = 1/1.5
      # strip.switch.pad.grid = unit(0.4, "cm") 
    )
  diff_plot + continent_dif_plot
  
} 



growth_vec <- unique(redlist_country_growths$growth_form)[c(5,3,4,2,7,1,6)]

cols <- c(nationalparkcolors::park_palette("Saguaro")[1:5], "#5b2635")
cols_re <- c(cols[c(6,4,5,3,2,1)], "black")


eqearth_crs <- "+proj=eqearth"


plot_list <- list()
plot_index <- 1


for (i in 1:length(growth_vec)) {
  
  # generate plots
  cols_re_p <- cols_re[i]
  
  a <- plot.growth.map.corr(growth = growth_vec[i]) + 
    theme(legend.position = "none") 
  
  
  
  # d <- plot_pca_f(continent_name = conts[i])+ 
  #   theme(legend.position = "none", 
  #         aspect.ratio = 1/1.5)
  # }
  
  # Save each plot individually into the flat list
  plot_list[[plot_index]] <- a
  plot_index <- plot_index + 1
  
}

combined_plot <- plot_grid(plotlist = plot_list, nrow = 7, labels = paste0("(", letters, ")")) + 
  plot_layout( guides = "collect") 


ggsave("red_growths_country_map_corr_09_2026.png", combined_plot, width = 10, height = 25, dpi = 900, 
       units = "in") 


ggsave("total_growths_country_map_corr_04_2026.png", a, width = 10, height = 5, dpi = 600, 
       units = "in") 

