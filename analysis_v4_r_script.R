# Description----
# Author: John Everett Parkinson
# Contact: jparkinson@usf.edu
# Date: 2021-02-28
# Publication: "Zoantharian endosymbiont community dynamics during a stress event"
# By: Fujiwara J, Kawamura I, Reimer JD, Parkinson JE
# Journal: Frontiers in Microbiology
# For Updates: https://github.com/parkinsonje

# Set working directory and seed----
#setwd("Add/Your/Pathway/Here/")
set.seed(123)

# Load libraries----
suppressMessages(suppressWarnings(library('reshape2')))
suppressMessages(suppressWarnings(library('dplyr')))
suppressMessages(suppressWarnings(library('ggplot2')))
suppressMessages(suppressWarnings(library('MCMC.qpcr')))
suppressMessages(suppressWarnings(library('car')))
suppressMessages(suppressWarnings(library('grid')))
suppressMessages(suppressWarnings(library('scales')))
suppressMessages(suppressWarnings(library('lubridate')))

# Define custom functions----

# Multiple plot function from: 
# <http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/>

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Load temp/light logger data and calculate summary stats----
data_logger_deep <- read.csv("input/data_logger_deep.csv") %>%
  mutate(lubri_date = mdy_hm(date)) %>%
  mutate(daylight = hour(lubri_date) > 5 & hour(lubri_date) < 19) %>%
    #define daylight hours (from 6:00 to 18:00)
  filter(daylight == TRUE) %>%
    # restrict to hours of interest (daylight hours)
  filter(lubri_date > mdy("8/19/2015") & lubri_date < mdy("9/3/2015"))
    # restrict to dates of interest (8/19/2015 to 9/2/2015)

data_logger_shallow <- read.csv("input/data_logger_shallow.csv") %>%
  mutate(lubri_date = mdy_hm(date)) %>%
  mutate(daylight = hour(lubri_date) > 5 & hour(lubri_date) < 19) %>%
    #define daylight hours (from 6:00 to 18:00)
  filter(daylight == TRUE) %>%
    # restrict to hours of interest (daylight hours)
  filter(lubri_date > mdy("8/19/2015") & lubri_date < mdy("9/3/2015"))
    # restrict to dates of interest (8/19/2015 to 9/2/2015)

data_logger_combined <- bind_rows(data_logger_deep, data_logger_shallow) %>%
  mutate(month = month(lubri_date)) %>%
  mutate(day = day(lubri_date)) %>%
  mutate(hour = hour(lubri_date)) %>%
  select(depth, month, day, hour, temp, light) %>%
  group_by(depth, hour) %>%
    # only need depth and hour for plotting purposes
  summarize(temp_mean = mean(temp), temp_sd = sd(temp), 
            light_mean = mean(light), light_sd = sd(light))
  
# Deep Temp Max Average: 28.54887 +/- 0.1929603
# Shallow Temp Max Avearge: 31.07160 +/- 2.8313019
# Deep Light Max Average: 1064.91333	+/- 492.649191
# Shallow Light Max Avearge: 25142.453 +/- 26013.0263


# Plot temp/light data----
plot_logger_light <- ggplot(data_logger_combined,
    aes(x = hour, y = light_mean, group = depth, fill = depth, color = depth)) +
  ylim(-3500, 51500) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray") +
  scale_x_continuous(breaks = seq(6, 18, by = 1)) +
  geom_point(position = position_dodge(width = 1), size =2) +
  geom_line(position = position_dodge(width = 1), size = 1) +
  geom_errorbar(aes(ymin = light_mean - light_sd, ymax = light_mean + light_sd),
                position = position_dodge(width = 1)) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) + 
    # first purple (deep), then orange (shallow)
  labs(x = "Time (24-hour clock)", y = "Illuminance (lx)") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'))


plot_logger_temp <- ggplot(data_logger_combined,
    aes(x = hour, y = temp_mean, group = depth, fill = depth, color = depth)) +
  ylim(28,34) +
  scale_x_continuous(breaks = seq(6, 18, by = 1)) +
  geom_point(position = position_dodge(width = 1), size =2) +
  geom_line(position = position_dodge(width = 1), size = 1) +
  geom_errorbar(aes(ymin = temp_mean - temp_sd, ymax = temp_mean + temp_sd),
                position = position_dodge(width = 1)) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) + 
    # first purple (deep), then orange (shallow)
  labs(x = "Time (24-hour clock)", y = expression(paste("Temperature (", degree, "C)"))) +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'))

multiplot(plot_logger_light, plot_logger_temp, cols = 2)
  # Print plots
  # ggsave(plot_logger_light, file = "fig2a_light.eps", width = 80, units = "mm", dpi = 300)
  # ggsave(plot_logger_temp, file = "fig2b_temp.eps", width = 80, units = "mm", dpi = 300)

# Load mortality data----
data_mortality <- read.csv("input/data_mortality.csv")

data_mortality$transplant <- factor(data_mortality$transplant,
  levels = c("SS", "SD", "DD", "DS"))

# Plot mortality data----
plot_mortality <- ggplot(data_mortality,
    aes(x = month, y = count, group = month, fill = status)) +
  geom_bar(position = "fill", stat = "identity", color = "black", width = 0.9) +
  scale_fill_manual(values = c("#a1d76a", "#e9a3c9", "#ffffff")) + 
    # first green (alive), then dead (pink), then missing (white)
  labs(x = "Transplant Period", y = "Proportion of Colonies with Health Status") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid')) +
  facet_wrap( ~ transplant, nrow = 2)

plot_mortality
  # Print plot
  # ggsave(plot_mortality, file = "fig3_mortality.eps", width = 160, units = "mm", dpi = 300)
  
# Load qpcr data----
data_stepone_raw <- read.csv("input/data_stepone_raw.csv") %>%
  transmute(experiment = Experiment, # Import raw stepone data; rename columns
            plate = Plate,
            colony = Colony,
            orientation = Orientation,
            sample_type = Task,
            sample_name = Sample.Name,
            assay = substr(Target.Name, 1, 1),
            quant_mean = Quantity.Mean) %>% 
  filter(sample_type != "") %>% # remove blanks
  filter(sample_type == "UNKNOWN") # exclude standards and controls

data_stepone_raw[is.na(data_stepone_raw)] <- 0
  # Replace NAs with 0 (changes non-amplifications to 0 quantity)

data_stepone_expa <- data_stepone_raw %>%
  filter(experiment == "A") %>% # retain only experiment A samples (transplants)
  group_by(experiment, plate, sample_type, sample_name, assay) %>%
  summarize(quant_mean = mean(quant_mean)) %>%
    # Calculate average initial quantity of template DNA in each sample.
    # Although stepone already provides the mean, this step simply removes
    # the extra replicate entry for each sample.
  mutate(transplant = substr(sample_name, 1, 2),
         time = substrRight(sample_name, 2),
         sample_name = gsub(".*[_]([^_]+)[_].*", "\\1", sample_name)) %>%
  filter(sample_name != "Clear25" & # unreliable amplification @ 0W
           sample_name != "White20" & # no amplification @ 0W
           sample_name != "White21" & # no amplification @ 0W
           sample_name != "White24" & # no amplification @ 4W
           sample_name != "Yellow14" & # no amplification @ 4W
           sample_name != "Yellow25") # no amplification @ 1W)

data_stepone_expb <- data_stepone_raw %>%
  filter(experiment == "B") %>% # retain only experiment B samples (intra-colony)
  group_by(experiment, plate, sample_type, colony, orientation, 
           sample_name, assay) %>%
  summarize(quant_mean = mean(quant_mean)) %>%
    # Calculate average initial quantity of template DNA in each sample.
    # Although stepone already provides the mean, this step simply removes
    # the extra replicate entry for each sample.
  mutate(sample_name = paste(colony, orientation, sep="_"))

data_stepone_expb[is.na(data_stepone_expb)] <- 0
  # Replace NAs with 0 (changes non-amplifications to 0 quantity)


data_stepone_expb$orientation <- factor(data_stepone_expb$orientation,
  levels = c("north", "east", "south", "west", "center"))
    # Re-level orientation factor

# Plot qpcr data (experiment A: transplant)----
var_transplant <- c("DD","DS","SD","SS")
  # Initiate list of transplant abbreviations for the loop
  # DD = deep to deep
  # DS = deep to shallow
  # SD = shallow to deep
  # SS = shallow to shallow

transplant_abundance_table <- data.frame()
# Create data frame to which loop output can be added

for( i in var_transplant) {
  
  A_data <- data_stepone_expa %>% 
    filter(assay=="A") %>% 
    filter(transplant==i) %>% 
    ungroup() %>%
    select(transplant, assay, time, sample_name, quant_mean) %>% 
    rename(A = quant_mean)
      # Isolate A (shallow Symbiodinium sp., r A1z) assay data
  
  S_data <- data_stepone_expa %>% 
    filter(assay=="S") %>% 
    filter(transplant==i) %>%
    ungroup() %>%
    select(quant_mean) %>% 
    rename(S = quant_mean)
      # Isolate S (shallow Cladocopium sp., or C1z_intertidal) assay data
  
  D_data <- data_stepone_expa %>% 
    filter(assay=="D") %>% 
    filter(transplant==i) %>%
    ungroup() %>%
    select(quant_mean) %>% 
    rename(D = quant_mean)
      # Isolate D (deep Cladocopium sp., or C1z_subtidal) assay data
  
  count_ASD <- cbind(A_data, S_data, D_data)
  count_ASD <- count_ASD %>% 
    mutate(total_count = A + S + D)
      # Merge assay data to calculate abundance across all three assays per sample
  
  percent_ASD <- count_ASD %>% 
    select(A, S, D) %>% 
    rename(prop_A = A, prop_S = S,prop_D = D)
  percent_ASD <- prop.table(as.matrix(percent_ASD), margin = 1)
  percent_ASD <- bind_cols(count_ASD,as.data.frame(percent_ASD)) %>%
    select(transplant, assay, time, sample_name, prop_A, prop_S, prop_D)
  percent_ASD <- melt(percent_ASD, variable.name = "Species", 
                      value.name = "Proportional_Abundance")
  percent_ASD <- percent_ASD %>% filter(Proportional_Abundance != "NaN")
    # Convert to proportional abundance by dividing assay count by total count
  
  levels(percent_ASD$Species)[levels(percent_ASD$Species)=="prop_A"] <- "A1z-shallow"
  levels(percent_ASD$Species)[levels(percent_ASD$Species)=="prop_S"] <- "C1z-shallow"
  levels(percent_ASD$Species)[levels(percent_ASD$Species)=="prop_D"] <- "C1z-deep"
    # Re-level to rename assays by species
  
  percent_ASD$sample_name[percent_ASD$sample_name=="Clear07"] <- "SS16"
  percent_ASD$sample_name[percent_ASD$sample_name=="Clear09"] <- "SS15"
  percent_ASD$sample_name[percent_ASD$sample_name=="Clear10"] <- "SS19"
  percent_ASD$sample_name[percent_ASD$sample_name=="Clear15"] <- "SS13"
  percent_ASD$sample_name[percent_ASD$sample_name=="Clear16"] <- "SS18"
  percent_ASD$sample_name[percent_ASD$sample_name=="Clear18"] <- "SS12"
  percent_ASD$sample_name[percent_ASD$sample_name=="Clear20"] <- "SS14"
  percent_ASD$sample_name[percent_ASD$sample_name=="Clear23"] <- "SS11"
  percent_ASD$sample_name[percent_ASD$sample_name=="Clear24"] <- "SS17"
    # Rename SS transplant samples for better graphing position
  
  percent_ASD$sample_name[percent_ASD$sample_name=="Red05"] <- "DD36"
  percent_ASD$sample_name[percent_ASD$sample_name=="Red09"] <- "DD35"
  percent_ASD$sample_name[percent_ASD$sample_name=="Red15"] <- "DD32"
  percent_ASD$sample_name[percent_ASD$sample_name=="Red16"] <- "DD33"
  percent_ASD$sample_name[percent_ASD$sample_name=="Red18"] <- "DD38"
  percent_ASD$sample_name[percent_ASD$sample_name=="Red20"] <- "DD34"
  percent_ASD$sample_name[percent_ASD$sample_name=="Red22"] <- "DD31"
  percent_ASD$sample_name[percent_ASD$sample_name=="Red23"] <- "DD39"
  percent_ASD$sample_name[percent_ASD$sample_name=="Red24"] <- "DD37"
    # Rename DD transplant samples for better graphing position
  
  percent_ASD$sample_name[percent_ASD$sample_name=="White11"] <- "SD22"
  percent_ASD$sample_name[percent_ASD$sample_name=="White15"] <- "SD23"
  percent_ASD$sample_name[percent_ASD$sample_name=="White16"] <- "SD27"
  percent_ASD$sample_name[percent_ASD$sample_name=="White17"] <- "SD24"
  percent_ASD$sample_name[percent_ASD$sample_name=="White18"] <- "SD26"
  percent_ASD$sample_name[percent_ASD$sample_name=="White19"] <- "SD25"
  percent_ASD$sample_name[percent_ASD$sample_name=="White22"] <- "SD29"
  percent_ASD$sample_name[percent_ASD$sample_name=="White23"] <- "SD28"
  percent_ASD$sample_name[percent_ASD$sample_name=="White27"] <- "SD21"
    # Rename SD transplant samples for better graphing position
  
  percent_ASD$sample_name[percent_ASD$sample_name=="Yellow12"] <- "DS47"
  percent_ASD$sample_name[percent_ASD$sample_name=="Yellow18"] <- "DS42"
  percent_ASD$sample_name[percent_ASD$sample_name=="Yellow20"] <- "DS45"
  percent_ASD$sample_name[percent_ASD$sample_name=="Yellow21"] <- "DS44"
  percent_ASD$sample_name[percent_ASD$sample_name=="Yellow22"] <- "DS43"
  percent_ASD$sample_name[percent_ASD$sample_name=="Yellow23"] <- "DS46"
  percent_ASD$sample_name[percent_ASD$sample_name=="Yellow24"] <- "DS49"
  percent_ASD$sample_name[percent_ASD$sample_name=="Yellow26"] <- "DS48"
  percent_ASD$sample_name[percent_ASD$sample_name=="Yellow28"] <- "DS41"
    # Rename DS transplant samples for better graphing position
  
  table_species_prop <- percent_ASD %>%
    select(transplant, sample_name, time, Species, Proportional_Abundance) %>%
    dcast(sample_name + time ~ Species) %>%
    cbind(transplant = i)
  
  transplant_abundance_table <- rbind(transplant_abundance_table, table_species_prop)

  assign(paste(i, "ra_plot_stepone", sep = "_"), 
    ggplot(percent_ASD,
      aes(x = time, y = Proportional_Abundance, group = Species, fill = Species)) +
      geom_blank() +

      geom_bar(stat = "identity", color = "black", width = 1) +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = c("#fc8d59", "#91bfdb", "#4575b4")) + 
        # first pink (A), then light blue (S), then dark blue (D)
      labs(title = paste(i), x = "Transplant Period", y = "Proportional Abundance") +
      theme_classic() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            panel.spacing = unit (0.5, "lines"),
            axis.text.x = element_text(angle = 90, hjust = 0),
            axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
            axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid')) +
      facet_wrap( ~ sample_name, nrow = 1)
  )
}

multiplot(SS_ra_plot_stepone, DD_ra_plot_stepone, SD_ra_plot_stepone,
          DS_ra_plot_stepone, cols = 2)
  # Print plots
  # ggsave(SS_ra_plot_stepone, file = "fig4a_transplant.eps", width = 80,
  #   height = 126, units = "mm", dpi = 300)
  # ggsave(DD_ra_plot_stepone, file = "fig4b_transplant.eps", width = 80,
  #   height = 126, units = "mm", dpi = 300)
  # ggsave(SD_ra_plot_stepone, file = "fig4c_transplant.eps", width = 80,
  #   height = 126, units = "mm", dpi = 300)
  # ggsave(DS_ra_plot_stepone, file = "fig4d_transplant.eps", width = 80,
  #   height = 126, units = "mm", dpi = 300)

  # ggsave(SS_ra_plot_stepone, file = "fig4a_transplant.eps", width = 160, units = "mm", dpi = 300)
  # ggsave(DD_ra_plot_stepone, file = "fig4b_transplant.eps", width = 160, units = "mm", dpi = 300)
  # ggsave(SD_ra_plot_stepone, file = "fig4c_transplant.eps", width = 160, units = "mm", dpi = 300)
  # ggsave(DS_ra_plot_stepone, file = "fig4d_transplant.eps", width = 160, units = "mm", dpi = 300)

transplant_abundance_table[,c(3,4,5)] <- round(transplant_abundance_table[c(3,4,5)], 2)
transplant_abundance_table
  #   sample_name   time A1z-shallow C1z-shallow C1z-deep transplant
  # 1         DD31   0W        0.00        0.00     1.00         DD
  # 2         DD31   4W        0.00        0.00     1.00         DD
  # 3         DD32   0W        0.00        0.00     1.00         DD
  # 4         DD32   4W        0.00        0.00     1.00         DD
  # 5         DD33   0W        0.00        0.00     1.00         DD
  # 6         DD33   4W        0.00        0.00     1.00         DD
  # 7         DD34   0W        0.00        0.00     1.00         DD
  # 8         DD34   4W        0.00        0.00     1.00         DD
  # 9         DD35   0W        0.00        0.00     1.00         DD
  # 10        DD35   4W        0.00        0.00     1.00         DD
  # 11        DD36   0W        0.00        0.00     1.00         DD
  # 12        DD36   4W        0.00        0.00     1.00         DD
  # 13        DD37   0W        0.00        0.00     1.00         DD
  # 14        DD37   4W        0.00        0.00     1.00         DD
  # 15        DD38   0W        0.00        0.01     0.99         DD
  # 16        DD38   4W        0.00        0.00     1.00         DD
  # 17        DD39   0W        0.00        0.00     1.00         DD
  # 18        DD39   4W        0.03        0.01     0.95         DD
  # 19        DS41   0W        0.00        0.00     1.00         DS
  # 20        DS41   1W        0.00        0.00     1.00         DS
  # 21        DS42   0W        0.00        0.00     1.00         DS
  # 22        DS42   1W        0.00        0.00     1.00         DS
  # 23        DS43   0W        0.00        0.00     1.00         DS
  # 24        DS43   1W        0.00        0.00     1.00         DS
  # 25        DS44   0W        0.00        0.00     1.00         DS
  # 26        DS44   1W        0.00        0.00     1.00         DS
  # 27        DS45   0W        0.00        0.00     1.00         DS
  # 28        DS45   1W        0.00        0.02     0.98         DS
  # 29        DS46   0W        0.00        0.00     1.00         DS
  # 30        DS46   1W        0.00        0.08     0.92         DS
  # 31        DS47   0W        0.00        0.00     1.00         DS
  # 32        DS47   1W        0.00        0.08     0.92         DS
  # 33        DS48   0W        0.00        0.00     1.00         DS
  # 34        DS48   1W        0.01        0.19     0.80         DS
  # 35        DS49   0W        0.00        0.00     1.00         DS
  # 36        DS49   1W        0.32        0.26     0.41         DS
  # 37        SD21   0W        0.00        0.90     0.10         SD
  # 38        SD21   4W        0.90        0.10     0.00         SD
  # 39        SD22   0W        0.99        0.00     0.00         SD
  # 40        SD22   4W        1.00        0.00     0.00         SD
  # 41        SD23   0W        1.00        0.00     0.00         SD
  # 42        SD23   4W        1.00        0.00     0.00         SD
  # 43        SD24   0W        1.00        0.00     0.00         SD
  # 44        SD24   4W        1.00        0.00     0.00         SD
  # 45        SD25   0W        1.00        0.00     0.00         SD
  # 46        SD25   4W        1.00        0.00     0.00         SD
  # 47        SD26   0W        0.99        0.01     0.00         SD
  # 48        SD26   4W        0.99        0.01     0.00         SD
  # 49        SD27   0W        0.98        0.01     0.01         SD
  # 50        SD27   4W        0.94        0.04     0.02         SD
  # 51        SD28   0W        0.95        0.05     0.00         SD
  # 52        SD28   4W        0.95        0.05     0.00         SD
  # 53        SD29   0W        0.00        1.00     0.00         SD
  # 54        SD29   4W        0.00        1.00     0.00         SD
  # 55        SS11   0W        0.00        1.00     0.00         SS
  # 56        SS11   4W        0.98        0.02     0.00         SS
  # 57        SS12   0W        0.00        1.00     0.00         SS
  # 58        SS12   4W        0.99        0.01     0.00         SS
  # 59        SS13   0W        0.83        0.10     0.07         SS
  # 60        SS13   4W        0.86        0.14     0.00         SS
  # 61        SS14   0W        0.83        0.17     0.00         SS
  # 62        SS14   4W        0.98        0.02     0.00         SS
  # 63        SS15   0W        0.95        0.03     0.01         SS
  # 64        SS15   4W        1.00        0.00     0.00         SS
  # 65        SS16   0W        0.99        0.01     0.00         SS
  # 66        SS16   4W        1.00        0.00     0.00         SS
  # 67        SS17   0W        1.00        0.00     0.00         SS
  # 68        SS17   4W        0.99        0.01     0.00         SS
  # 69        SS18   0W        1.00        0.00     0.00         SS
  # 70        SS18   4W        0.96        0.04     0.00         SS
  # 71        SS19   0W        0.92        0.05     0.03         SS
  # 72        SS19   4W        0.00        1.00     0.00         SS

transplant_abundance_table_melt <- melt(transplant_abundance_table) %>%
  transmute(transplant = transplant, time = time, species = variable,
            proportional_abundance = value)

plot_transplant <-
  ggplot(transplant_abundance_table_melt,
         aes(x = "", y = proportional_abundance, fill = species)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("#fc8d59", "#91bfdb", "#4575b4")) + 
  # first pink (A), then light blue (S), then dark blue (D)
  coord_polar("y", start = 0) +
  theme_void() +
  facet_wrap(~ transplant + time, nrow = 4)
plot_transplant 
  # Print plot (not used in publication)

# Plot qpcr data (experiment B: intra-colony)----
var_colony <- c("colony_one","colony_two","colony_three","colony_four")
  # Initiate list of colonies for the loop

colony_abundance_table <- data.frame()
  # Create data frame to which loop output can be added

for( i in var_colony) {
  
  A_data <- data_stepone_expb %>% 
    filter(assay=="A") %>% 
    filter(colony==i) %>% 
    ungroup() %>%
    select(colony, orientation, assay, sample_name, quant_mean) %>% 
    rename(A = quant_mean)
  # Isolate A (shallow Symbiodinium sp., or A1z) assay data
  
  S_data <- data_stepone_expb %>% 
    filter(assay=="S") %>% 
    filter(colony==i) %>%
    ungroup() %>%
    select(quant_mean) %>% 
    rename(S = quant_mean)
  # Isolate S (shallow Cladocopium sp., or C1z_intertidal) assay data
  
  D_data <- data_stepone_expb %>% 
    filter(assay=="D") %>% 
    filter(colony==i) %>%
    ungroup() %>%
    select(quant_mean) %>% 
    rename(D = quant_mean)
  # Isolate D (deep Cladocopium sp., or C1z_subtidal) assay data
  
  count_ASD <- cbind(A_data, S_data, D_data)
  count_ASD <- count_ASD %>% 
    mutate(total_count = A + S + D)
  # Merge assay data to calculate abundance across all three assays per sample
  
  percent_ASD <- count_ASD %>% 
    select(A, S, D) %>% 
    rename(prop_A = A, prop_S = S,prop_D = D)
  percent_ASD <- prop.table(as.matrix(percent_ASD), margin = 1)
  percent_ASD <- bind_cols(count_ASD,as.data.frame(percent_ASD)) %>%
    select(colony, orientation, assay, sample_name, prop_A, prop_S, prop_D)
  percent_ASD <- melt(percent_ASD, variable.name = "Species",
                      value.name = "Proportional_Abundance")
  percent_ASD <- percent_ASD %>% filter(Proportional_Abundance != "NaN")
  # Convert to proportional abundance by dividing assay count by total count
  
  levels(percent_ASD$Species)[levels(percent_ASD$Species)=="prop_A"] <- "A1z-shallow"
  levels(percent_ASD$Species)[levels(percent_ASD$Species)=="prop_S"] <- "C1z-shallow"
  levels(percent_ASD$Species)[levels(percent_ASD$Species)=="prop_D"] <- "C1z-deep"
  # Re-level to rename assays by species
  
  table_species_prop <- percent_ASD %>%
    select(colony, orientation, Species, Proportional_Abundance) %>%
    dcast(orientation ~ Species) %>%
    cbind(colony = i)
  
  colony_abundance_table <- rbind(colony_abundance_table, table_species_prop)
  
  assign(paste(i, "ra_plot_stepone", sep = "_"), 
    ggplot(percent_ASD,
      aes(x = orientation,y = Proportional_Abundance, group = Species, fill = Species)) +
      geom_blank() +
        #geom_hline(yintercept = 0, linetype = "solid") +
      geom_bar(stat = "identity", color = "black", width = 1) +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = c("#fc8d59", "#91bfdb", "#4575b4")) + 
        # first pink (A), then light blue (S), then dark blue (D)
      labs(title = paste(i), x = "Orientation", y = "Proportional Abundance") +
      theme_classic() +
      theme(legend.position = "bottom",
             legend.title = element_blank(),
             panel.spacing = unit (0.5, "lines"),
             axis.text.x = element_text(angle = 90, hjust = 0),
             axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
             axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'))
   )
}

multiplot(colony_one_ra_plot_stepone, colony_two_ra_plot_stepone, 
          colony_three_ra_plot_stepone, colony_four_ra_plot_stepone, cols=2)
  # Print plots (not used in publication)

colony_abundance_table[,c(2,3,4)] <- round(colony_abundance_table[c(2,3,4)], 2)
colony_abundance_table
  # orientation   A1z-shallow  C1z-shallow   C1z-deep           colony
  #       north           0.93         0.07          0       colony_one
  #        east           0.99         0.01          0       colony_one
  #       south           0.16         0.84          0       colony_one
  #        west           0.14         0.86          0       colony_one
  #      center           0.98         0.02          0       colony_one
  #       north           0.00         1.00          0       colony_two
  #        east           0.00         1.00          0       colony_two
  #       south           0.00         1.00          0       colony_two
  #        west           0.02         0.98          0       colony_two
  #      center           0.02         0.98          0       colony_two
  #       north           0.99         0.01          0     colony_three
  #        east           0.34         0.66          0     colony_three
  #       south           0.96         0.04          0     colony_three
  #        west           0.99         0.01          0     colony_three
  #      center           0.13         0.87          0     colony_three
  #       north           0.92         0.08          0      colony_four
  #        east           0.00         1.00          0      colony_four
  #       south           0.00         1.00          0      colony_four
  #        west           0.00         1.00          0      colony_four
  #      center           0.84         0.16          0      colony_four

colony_abundance_table_melt <- melt(colony_abundance_table) %>%
  transmute(colony = colony, orientation = orientation, species = variable,
            proportional_abundance = value)

plot_colony <-
  ggplot(colony_abundance_table_melt,
    aes(x = "", y = proportional_abundance, fill = species)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = c("#fc8d59", "#91bfdb", "#4575b4")) + 
    # first pink (A), then light blue (S), then dark blue (D)
    coord_polar("y", start = 0) +
    theme_void() +
    facet_wrap(~ colony + orientation, nrow = 4)
plot_colony 
  #Print plot
  # ggsave(plot_colony, file = "fig6_colony.eps", width = 180, units = "mm", dpi = 300)