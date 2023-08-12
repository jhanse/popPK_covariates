# Load libraries
library(tidyverse)
library(readr)
library(ggplot2)
library(lubridate)
library(epitools)
library(gridExtra)
library(rstudioapi)

# Import dataset
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
cis <- read_csv('cis.csv')
outcomes <- read_csv('outcomes.csv')

# Group datasets for covariate analysis
out_merged <- outcomes %>% 
  left_join(cis, by = 'study_id')

#############################################
## Basic study and outcome characteristics ##
#############################################

# Descriptive: Number of participants, median and IQR
sum(cis$number_of_patients)
summary(cis$number_of_patients)

# Descriptive: Covariate categories frequencies
out_merged %>% 
  group_by(outcome_category) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  mutate(perc = n/nrow(out_merged)*100) %>% 
  print(n = 20)

# Descriptive: Patient characteristics
pt_char <- data.frame(chars = c(cis$patient_characteristics1, cis$patient_characteristics2, cis$patient_characteristics3, cis$patient_characteristics4)) %>% na.omit() %>% 
  group_by(chars) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  mutate(percent = n/nrow(cis)*100) %>% 
  print(n = 40)
pt_char

# Descriptive: Study country economic categories based on study design
cis %>% 
  group_by(study_design, country_economic_category) %>% 
  count() %>% 
  mutate(freq = round(n/nrow(cis)*100, digits = 1))

# Descriptive: Types of blood samples used based on study design
cis %>% 
  select(study_design, samples_used) %>% 
  group_by(study_design, samples_used) %>%  
  count() %>% 
  mutate(freq = round(n/nrow(cis)*100, digits = 1)) %>% 
  arrange(desc(n))

# Descriptive: Non-blood samples used based on study design 
cis %>% 
  select(study_design, other_samples) %>% 
  group_by(study_design, other_samples) %>%  
  count() %>% 
  mutate(freq = round(n/nrow(cis)*100, digits = 1)) %>% 
  arrange(desc(n)) %>% 
  filter(other_samples != 'no')

# Descriptive: Model testing methods
cis %>%
  select(model_evaluation_process) %>% 
  mutate(model_evaluation_process = strsplit(model_evaluation_process, ", ")) %>% 
  unnest(model_evaluation_process) %>% 
  group_by(model_evaluation_process) %>% 
  count() %>% 
  na.omit() %>% 
  mutate(freq = round(n/nrow(cis)*100, digits = 1)) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# Descriptive: Antimicrobials used
cis %>% 
  select(betalactam_studied) %>% 
  separate_rows(betalactam_studied, sep = ', ') %>% 
  group_by(betalactam_studied) %>% 
  mutate(n = n(), abx = str_to_title(betalactam_studied)) %>%
  ungroup() %>% 
  select(abx, n) %>% 
  distinct() %>% 
  mutate(freq = round(n/nrow(cis)*100, digits = 1)) %>% 
  arrange(desc(n)) %>% 
  as.data.frame()

# Descriptive: Extracorporeal techniques used
cis %>% 
  select(extracorporeal) %>% 
  separate_rows(extracorporeal, sep = ', ') %>% 
  group_by(extracorporeal) %>% 
  mutate(n = n()) %>%
  select(extracorporeal, n) %>% 
  distinct() %>% 
  mutate(freq = round(n/nrow(cis)*100, digits = 1)) %>% 
  arrange(desc(n)) %>% 
  as.data.frame()

#######################
## Subgroup analyses ##
#######################

# Income country
income <- out_merged %>% 
  filter(!is.na(significant)) %>% 
  select(country_economic_category1, significant)

# Antimicrobial category
abx <- out_merged %>% 
  filter(!is.na(significant)) %>% 
  select(abx_cat, significant) 

# Renal function
renal <- out_merged %>% 
  filter(!is.na(significant)) %>%
  filter(renal_function != 'unclear') %>% 
  select(renal_function, significant) 

# Odds ratios for above
or_income <- table(income$significant, income$country_economic_category1) %>% oddsratio()
or_abx <- table(abx$significant, abx$abx_cat) %>% oddsratio()
or_renal <- table(renal$significant, renal$renal_function) %>% oddsratio()

or_dataset <- data.frame(
  index = c(1,2,3),
  study_label = c('Renal function', 'Income level', 'Antimicrobial category'),
  OR = c(or_renal$measure[2], or_income$measure[2], or_abx$measure[2]),
  lci = c(or_renal$measure[2,2], or_income$measure[2,2], or_abx$measure[2,2]),
  uci = c(or_renal$measure[2,3], or_income$measure[2,3], or_abx$measure[2,3])) %>% 
  mutate(ci = paste0(round(lci, 3), ", ", round(uci,3)))
or_dataset

# Renal function
covars_renal <- out_merged %>% select(outcome_cleaner, significant) %>% 
  group_by(outcome_cleaner, significant) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  filter(n > 9)
covars_renal <- covars_renal$outcome_cleaner

for (x in covars_renal) {
  outcomes_renal_function <- out_merged %>% 
    select(study_id, significant, outcome_cleaner, renal_function) %>% 
    filter(renal_function != 'augmented' & renal_function != 'unclear') %>% 
    filter(outcome_cleaner == x)
  print(x)
  table(outcomes_renal_function$renal_function, outcomes_renal_function$significant) %>% print()
  table(outcomes_renal_function$significant, outcomes_renal_function$renal_function) %>% summary() %>% print()
}

# Economic area based on renal function
covars_eco <- c("weight", "serum albumin")
for (x in covars_eco) {
  outcomes_renal_function <- out_merged %>% 
    select(study_id, significant, outcome_cleaner, country_economic_category1) %>% 
    filter(outcome_cleaner == x)
  print(x)
  table(outcomes_renal_function$country_economic_category1, outcomes_renal_function$significant) %>% print()
  table(outcomes_renal_function$country_economic_category1, outcomes_renal_function$significant) %>% oddsratio(correction = TRUE) %>% print()
  table(outcomes_renal_function$significant, outcomes_renal_function$country_economic_category1) %>% summary() %>% print()
}

# Antimicrobial category based on renal function
covars_abx <- c("weight", "serum albumin", "serum creatinine")
for (x in covars_abx) {
  outcomes_renal_function <- out_merged %>% 
    select(study_id, significant, outcome_cleaner, abx_cat) %>% 
    filter(outcome_cleaner == x)
  print(x)
  table(outcomes_renal_function$abx_cat, outcomes_renal_function$significant) %>% print()
  table(outcomes_renal_function$abx_cat, outcomes_renal_function$significant) %>% oddsratio() %>% print()
  table(outcomes_renal_function$abx_cat, outcomes_renal_function$significant) %>% summary() %>% print()
}

###########
## Plots ##
###########

# Study frequency as a function of time
hist(cis$year_date, breaks = 39)

### Sunburst plots

# Plot 1
d1 <- outcomes %>% 
  select(outcome_cleaner, outcome_category) %>% 
  group_by(outcome_cleaner) %>% 
  count(num = n()) %>% 
  ungroup()

d2 <- outcomes %>%
  right_join(d1, by = 'outcome_cleaner', unmatched = 'drop') %>% 
  filter(n > 2)

inner <- d2 %>%
  group_by(outcome_category) %>%
  count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  rename(category_n = n)

middle <- d2 %>%
  group_by(outcome_category, outcome_cleaner) %>%
  count() %>%
  ungroup() %>%
  rename(outcome_n = n) %>%
  left_join(inner, by = "outcome_category") %>%
  group_by(outcome_category) %>%
  arrange(desc(category_n), desc(outcome_n)) %>%
  ungroup()

outer <- d2 %>%
  group_by(outcome_category, outcome_cleaner, significant) %>%
  count() %>%
  ungroup() %>%
  rename(significant_n = n) %>%
  left_join(middle, by = c("outcome_category", "outcome_cleaner")) %>%
  group_by(outcome_category, outcome_cleaner) %>%
  arrange(desc(category_n), desc(outcome_n), outcome_cleaner, desc(significant)) %>%
  ungroup()

inner <- inner %>%
  mutate(prev_cat_n = lag(category_n),
         prev_cat_n = replace_na(prev_cat_n, 0),
         sum_cat_n = cumsum(category_n),
         prev_cat_sum = lag(sum_cat_n),
         prev_cat_sum = replace_na(prev_cat_sum, 0)) %>%
  rename(ymin = "prev_cat_sum",
         ymax = "sum_cat_n") %>%
  rowwise %>%
  mutate(mid_y = mean(c(ymin, ymax)))

middle <- middle %>%
  mutate(prev_out_n = lag(outcome_n),
         prev_out_n = replace_na(prev_out_n, 0),
         sum_out_n = cumsum(outcome_n),
         prev_out_sum = lag(sum_out_n),
         prev_out_sum = replace_na(prev_out_sum, 0)) %>%
  rename(ymin = "prev_out_sum",
         ymax = "sum_out_n") %>%
  rowwise %>%
  mutate(mid_y = mean(c(ymin, ymax)))

outer <- outer %>%
  mutate(prev_sig_n = lag(significant_n),
         prev_sig_n = replace_na(prev_sig_n, 0),
         sum_sig_n = cumsum(significant_n),
         prev_sig_sum = lag(sum_sig_n),
         prev_sig_sum = replace_na(prev_sig_sum, 0)) %>%
  rename(ymin = "prev_sig_sum",
         ymax = "sum_sig_n") %>%
  rowwise %>%
  mutate(mid_y = mean(c(ymin, ymax)))

outer <- outer %>%
  mutate(alpha_val = if_else(significant == "no" | is.na(significant), 0.6, 1))

lines <- tibble(inner_xmin = 1, inner_xmax = 5.0, y = inner$ymax)

circular_angle <- function(x) {
  ifelse(
    x/max(x)*360 > 180, 
    -90-(x/max(x)*360), 
    (-90 - x/max(x)*360)+180
  )
} 

p1 <- ggplot() +
  geom_rect(data = inner, mapping = aes(fill = outcome_category, colour = outcome_category, ymax = ymax, ymin = ymin, xmax = 2.8, xmin = 1)) +
  scale_fill_manual(name = "Covariate category", values =c('#204E79', '#EC7D31', '#8496B0', '#61A5C2', '#89C2D9', '#558699', '#2C7DA0', '#468FAF', '#61A5C2', '#A9D6E5')) +
  scale_colour_manual(name = "Covariate category", values =c('#204E79', '#EC7D31', '#8496B0', '#61A5C2', '#89C2D9', '#558699', '#2C7DA0', '#468FAF', '#61A5C2', '#A9D6E5')) +
  geom_rect(data = middle, mapping = aes(fill = outcome_category, ymax = ymax, ymin = ymin, xmax = 3, xmin = 4.4), colour = "white", linewidth = 0.2) +
  geom_rect(data = outer, mapping = aes(fill = outcome_category, ymax = ymax, ymin = ymin, xmax = 4.6, xmin = 5, alpha = alpha_val), colour = "white", linewidth = 0.2) +
  scale_alpha(guide = "none", range = c(0.5, 1)) +
  geom_segment(data = lines, mapping = aes(x = inner_xmin, xend = inner_xmax, y = y, yend = y), colour = "black") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 5.0)) +
  coord_polar(theta = "y") +
  guides(alpha = "none") +
  geom_text(data = inner, aes(x = 2, y = category_n, label = paste0(category_n, "\n", round(category_n/1077*100), "%")), position = position_stack(vjust = 0.5), colour = 'white') +
  geom_text(size = 2, check_overlap = F, data = middle, aes(angle = circular_angle(mid_y), x = 3.7, y = mid_y, label = outcome_cleaner), colour = 'white')

# Plot 2
sunburst <- out_merged %>% 
  filter(significant == 'yes') %>% 
  select(outcome_cleaner, outcome_category, parameter_clean) %>% 
  drop_na() %>% 
  mutate(para = strsplit(parameter_clean, ", ")) %>% 
  unnest(para) %>% 
  filter(para == "V" | para == "CL") %>% 
  select(outcome_cleaner, outcome_category, para)

d3 <- sunburst %>% 
  select(outcome_cleaner, outcome_category) %>% 
  group_by(outcome_cleaner) %>% 
  count(num = n()) %>% 
  ungroup()

d4 <- sunburst %>% 
  right_join(d3, by = "outcome_cleaner", unmatched = 'drop') %>% 
  filter(n > 0)

inner2 <- d4 %>%
  group_by(outcome_category) %>%
  count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  rename(category_n = n)

middle2 <- d4 %>%
  group_by(outcome_category, outcome_cleaner) %>%
  count() %>%
  ungroup() %>%
  rename(outcome_n = n) %>%
  left_join(inner2, by = "outcome_category") %>%
  group_by(outcome_category) %>%
  arrange(desc(category_n), desc(outcome_n)) %>%
  ungroup()

outer2 <- d4 %>%
  group_by(outcome_category, outcome_cleaner, para) %>% 
  count() %>% 
  ungroup() %>% 
  rename(para_n = n)  %>% 
  left_join(middle2, by = c("outcome_category", "outcome_cleaner")) %>% 
  group_by(outcome_category, outcome_cleaner) %>% 
  arrange(desc(category_n), desc(outcome_n), outcome_cleaner, desc(para)) %>%
  ungroup()

inner2 <- inner2 %>%
  mutate(prev_cat_n = lag(category_n),
         prev_cat_n = replace_na(prev_cat_n, 0),
         sum_cat_n = cumsum(category_n),
         prev_cat_sum = lag(sum_cat_n),
         prev_cat_sum = replace_na(prev_cat_sum, 0)) %>%
  rename(ymin = "prev_cat_sum",
         ymax = "sum_cat_n") %>%
  rowwise %>%
  mutate(mid_y = mean(c(ymin, ymax)))

middle2 <- middle2 %>%
  mutate(prev_out_n = lag(outcome_n),
         prev_out_n = replace_na(prev_out_n, 0),
         sum_out_n = cumsum(outcome_n),
         prev_out_sum = lag(sum_out_n),
         prev_out_sum = replace_na(prev_out_sum, 0)) %>%
  rename(ymin = "prev_out_sum",
         ymax = "sum_out_n") %>%
  rowwise %>%
  mutate(mid_y = mean(c(ymin, ymax)))

outer2 <- outer2 %>%
  mutate(prev_para_n = lag(para_n),
         prev_para_n = replace_na(prev_para_n, 0),
         sum_para_n = cumsum(para_n),
         prev_para_sum = lag(sum_para_n),
         prev_para_sum = replace_na(prev_para_sum, 0)) %>%
  rename(ymin = "prev_para_sum",
         ymax = "sum_para_n") %>%
  rowwise %>%
  mutate(mid_y = mean(c(ymin, ymax)))

outer2 <- outer2 %>%
  mutate(alpha_val = if_else(para == "CL", 1, 0.6))

lines <- tibble(inner_xmin = 1, inner_xmax = 5.0, y = inner2$ymax)

p2 <- ggplot() +
  geom_rect(data = inner2, mapping = aes(fill = outcome_category, colour = outcome_category, ymax = ymax, ymin = ymin, xmax = 2.8, xmin = 1)) +
  scale_fill_manual(name = "Category", values =c('#204E79', '#EC7D31', '#8496B0', '#61A5C2', '#89C2D9', '#558699', '#2C7DA0', '#468FAF', '#61A5C2', '#A9D6E5')) +
  scale_colour_manual(name = "Category", values =c('#204E79', '#EC7D31', '#8496B0', '#61A5C2', '#89C2D9', '#558699', '#2C7DA0', '#468FAF', '#61A5C2', '#A9D6E5')) +
  geom_rect(data = middle2, mapping = aes(fill = outcome_category, ymax = ymax, ymin = ymin, xmax = 3, xmin = 4.4), colour = "white", linewidth = 0.2) +
  geom_rect(data = outer2, mapping = aes(fill = outcome_category, ymax = ymax, ymin = ymin, xmax = 4.6, xmin = 5, alpha = alpha_val), colour = "white", linewidth = 0.2) +
  scale_alpha(guide = "none", range = c(0.5, 1)) +
  geom_segment(data = lines, mapping = aes(x = inner_xmin, xend = inner_xmax, y = y, yend = y), colour = "black") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 5.0)) +
  coord_polar(theta = "y") +
  guides(alpha = "none") +
  geom_text(data = inner2, aes(x = 2, y = category_n, label = paste0(category_n, "\n", round(category_n/243*100), "%")), position = position_stack(vjust = 0.5), colour = 'white') +
  geom_text(size = 2, 
            check_overlap = F, 
            data = middle2, 
            aes(angle = circular_angle(mid_y), x = 3.7, y = mid_y, label = outcome_cleaner), 
            colour = 'white')

# Merged plots
grid.arrange(p1, p2, ncol=2)

### Heatmaps

# Heatmap of countries
countries <- cis %>% 
  select(country) %>% 
  group_by(country) %>% 
  mutate(n = n()) %>% 
  distinct() %>% 
  arrange(desc(n)) %>% 
  as.data.frame()

world_map <- map_data("world")
world_map <- subset(world_map, region!="Antarctica")

ggplot(countries) +
  geom_map(dat = world_map, map = world_map, aes(map_id = region), fill = "white", color = "#7f7f7f", linewidth = 0.25) +
  geom_map(map = world_map, aes(map_id = country, fill = n), linewidth = 0.25) +
  scale_fill_gradient(low = "#DAE3F3", high = "#204E79", name = "Studies") +
  expand_limits(x = world_map$long, y = world_map$lat) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Heatmap of covariate categories vs antimicrobials
hm_data1 <- out_merged %>%
  select(betalactam_studied, outcome_category) %>% 
  mutate(betalactam_studied = strsplit(betalactam_studied, ", ")) %>% 
  unnest(betalactam_studied) %>% 
  group_by(betalactam_studied, outcome_category) %>% 
  count() %>% 
  arrange(desc(n))

ggplot(data = hm_data1, mapping = aes(x = outcome_category, y = betalactam_studied, fill = n, label = n)) +
  geom_tile() +
  geom_text(aes(label=n), alpha = 0.4, size=3) +
  theme_minimal() +
  xlab(label = "Covariate category") + ylab(label = "Antimicrobial agent") + labs(fill = "Covariates (n)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient(low = '#DAE3F3', high = '#204E79')

# Covars in more than one study
hm_data2 <- out_merged %>%
  select(betalactam_studied, outcome_cleaner) %>% 
  mutate(betalactam_studied = strsplit(betalactam_studied, ", ")) %>% 
  unnest(betalactam_studied) %>% 
  group_by(betalactam_studied, outcome_cleaner) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  filter(n > 1)

ggplot(data = hm_data2, mapping = aes(x = reorder(outcome_cleaner, -n), y = betalactam_studied, fill = n)) +
  geom_tile() +
  theme_minimal() +
  xlab(label = "Covariate") + ylab(label = "Antimicrobial agent") + labs(fill = "Number of studies") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "bottom") +
  scale_fill_gradient(low = '#DAE3F3', high = '#204E79')

# Covars in one study only per antimicrobial
hm_data3 <- out_merged %>%
  select(betalactam_studied, outcome_cleaner) %>% 
  mutate(betalactam_studied = strsplit(betalactam_studied, ", ")) %>% 
  unnest(betalactam_studied) %>% 
  group_by(betalactam_studied, outcome_cleaner) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  filter(n < 2)

ggplot(data = hm_data3, mapping = aes(x = reorder(outcome_cleaner, -n), y = betalactam_studied, fill = n)) +
  geom_tile() +
  theme_minimal() +
  xlab(label = "Covariate") + ylab(label = "Antimicrobial agent") + labs(fill = "Number of studies") +
  theme(legend.position = "none") +  
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  scale_fill_gradient(low = '#DAE3F3', high = '#204E79')

### Forest plot of ORs
or_plot <- ggplot(or_dataset, aes(y = index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:3, labels = or_dataset$study_label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

table_base <- ggplot(or_dataset, aes(y=study_label)) +
  ylab(NULL) + xlab("  ") + 
  theme(plot.title = element_text(hjust = 0.5, size=12), 
        axis.text.x = element_text(color="white", hjust = -3, size = 25), ## This is used to help with alignment
        axis.line = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank())

tab1 <- table_base + 
  labs(title = "space") +
  geom_text(aes(y = rev(index), x = 1, label = sprintf("%0.1f", round(OR, digits = 1))), size = 4) + ## decimal places
  ggtitle("OR")

tab2 <- table_base +
  geom_text(aes(y = rev(index), x = 1, label = ci), size = 4) + 
  ggtitle("95% CI")

lay <-  matrix(c(1,1,1,1,1,1,1,1,1,1,2,3,3), nrow = 1)
grid.arrange(or_plot, tab1, tab2, layout_matrix = lay)
