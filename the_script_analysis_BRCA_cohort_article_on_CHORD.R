######################################################################################
######################################################################################

setwd("./ARTICLE_CHORD_SUPPLEMENTARY_FILES_TO_DISPLAY_make_HEATMAPS_BOXPLOTS_SCATTER")

suppressWarnings(library("tidyverse"))
suppressWarnings(library("ggplot2"))
suppressWarnings(library("reshape2"))

######################################################################################
######################################################################################

df = read.table("the_table_CHORD_article.BRCA.txt", header=T, sep="\t", stringsAsFactors=FALSE)

head(df, 2)
colnames(df)
table(df$hr_status)
table(df$hrd_type)

######################################################################################
######################################################################################

df.hr_status = df %>% select("hr_status", "p_BRCA1", "p_BRCA2", "p_none", "p_hrd", 
                              "indel_load", "indel_rep_load", "sv_load" )
head(df.hr_status,2)
dim(df.hr_status)
table(df.hr_status$hr_status)
df.hr_status <- df.hr_status[complete.cases(df.hr_status$hr_status), ]
df.hr_status$hr_status = factor(df.hr_status$hr_status)

df.hrd_type = df %>% select("hrd_type", "p_BRCA1", "p_BRCA2", "p_none", "p_hrd", 
                              "indel_load", "indel_rep_load", "sv_load" )
head(df.hrd_type,2)
dim(df.hrd_type)
table(df.hrd_type$hrd_type)
df.hrd_type <- df.hrd_type[complete.cases(df.hrd_type$hrd_type), ]
df.hrd_type$hrd_type = factor(df.hrd_type$hrd_type)

df.cancer_type = df %>% select("cancer_type", "p_BRCA1", "p_BRCA2", "p_none", "p_hrd", 
                              "indel_load", "indel_rep_load", "sv_load" )
head(df.cancer_type,2)
dim(df.cancer_type)
table(df.cancer_type$cancer_type)
df.cancer_type <- df.cancer_type[complete.cases(df.cancer_type$cancer_type), ]
df.cancer_type$cancer_type = factor(df.cancer_type$cancer_type) 

######################################################################################
######################################################################################
# Working with HR.STATUS

colnames(df.hr_status)
df.hrd_type <- df.hrd_type[complete.cases(df.hrd_type$hrd_type), ]

p.hr_status <- df.hr_status %>%
  gather(key = "variable", value = "value", -hr_status) %>%
  ggplot(aes(x = hr_status, y = value, fill = hr_status)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "HR STATUS", x = "hr_status", y = "Value") +
  theme_minimal()

# Set global options for plot size
options(repr.plot.width=10, repr.plot.height=10)
p.hr_status

summary(df.hr_status)

######################################################################################
######################################################################################
# display HR_STATUS function of p_BRCA1
######################################################################################
######################################################################################

summary_stats <- tapply(df.hr_status$p_BRCA1, df.hr_status$hr_status, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_BRCA1 <- ggplot(df.hr_status, aes(x = hr_status, y = p_BRCA1, fill = hr_status, color = hr_status)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "p_BRCA1 by hr_status", x = "hr_status", y = "p_BRCA1") +
  scale_color_manual(name = "hr_status", 
                     values = c("HR_proficient" = "blue", 
                                "HR_deficient" = "red",
                                "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  scale_fill_manual(name = "hr_status", values = c("HR_proficient" = "lightblue", 
                                                   "HR_deficient" = "lightcoral", 
                                                   "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  )

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(p_BRCA1)

######################################################################################
######################################################################################
# display HR_STATUS function of p_BRCA2
######################################################################################
######################################################################################
# display HR_STATUS function of p_BRCA2

summary_stats <- tapply(df.hr_status$p_BRCA2, df.hr_status$hr_status, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_BRCA2 <- ggplot(df.hr_status, aes(x = hr_status, y = p_BRCA2, fill = hr_status, color = hr_status)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "p_BRCA2 by hr_status", x = "hr_status", y = "p_BRCA2") +
  scale_color_manual(name = "hr_status", 
                     values = c("HR_proficient" = "blue", 
                                "HR_deficient" = "red",
                                "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  scale_fill_manual(name = "hr_status", values = c("HR_proficient" = "lightblue", 
                                                   "HR_deficient" = "lightcoral", 
                                                   "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  )

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(p_BRCA2)


######################################################################################
######################################################################################
# display HR_STATUS function of p_hrd
######################################################################################
######################################################################################
# display HR_STATUS function of p_hrd
######################################################################################
######################################################################################

summary_stats <- tapply(df.hr_status$p_hrd, df.hr_status$hr_status, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))

p_hrd <- ggplot(df.hr_status, aes(x = hr_status, y = p_hrd, fill = hr_status, color = hr_status)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "p_hrd by hr_status", x = "hr_status", y = "p_hrd") +
  scale_color_manual(name = "hr_status", 
                     values = c("HR_proficient" = "blue", 
                                "HR_deficient" = "red",
                                "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  scale_fill_manual(name = "hr_status", values = c("HR_proficient" = "lightblue", 
                                                   "HR_deficient" = "lightcoral", 
                                                   "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  )

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(p_hrd)


######################################################################################
######################################################################################
# Correlation Coefficients
######################################################################################
######################################################################################

cor(df.hr_status$p_BRCA1, df.hr_status$p_BRCA2)
cor(df.hr_status$p_BRCA1, df.hr_status$p_hrd)
cor(df.hr_status$p_BRCA2, df.hr_status$p_hrd)

# display HR_STATUS function of indel_load

summary_stats <- tapply(df.hr_status$indel_load, df.hr_status$hr_status, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_indel_load <- ggplot(df.hr_status, aes(x = hr_status, y = indel_load, fill = hr_status, color = hr_status)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "indel_load by hr_status", x = "hr_status", y = "indel_load") +
  scale_color_manual(name = "hr_status", 
                     values = c("HR_proficient" = "blue", 
                                "HR_deficient" = "red",
                                "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  scale_fill_manual(name = "hr_status", values = c("HR_proficient" = "lightblue", 
                                                   "HR_deficient" = "lightcoral", 
                                                   "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  ) +
  ylim(c(0,6000))

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(p_indel_load)

######################################################################################
######################################################################################
# display HR_STATUS function of indel_rep_load
######################################################################################
######################################################################################

summary_stats <- tapply(df.hr_status$indel_rep_load, df.hr_status$hr_status, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_indel_rep_load <- ggplot(df.hr_status, aes(x = hr_status, y = indel_rep_load, fill = hr_status, color = hr_status)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "indel_rep_load by hr_status", x = "hr_status", y = "indel_rep_load") +
  scale_color_manual(name = "hr_status", 
                     values = c("HR_proficient" = "blue", 
                                "HR_deficient" = "red",
                                "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  scale_fill_manual(name = "hr_status", values = c("HR_proficient" = "lightblue", 
                                                   "HR_deficient" = "lightcoral", 
                                                   "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  ) + 
  ylim(c(0,6000))


# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(p_indel_rep_load)


######################################################################################
######################################################################################
# display HR_STATUS function of sv_load
######################################################################################
######################################################################################


summary_stats <- tapply(df.hr_status$sv_load, df.hr_status$hr_status, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_sv_load <- ggplot(df.hr_status, aes(x = hr_status, y = sv_load, fill = hr_status, color = hr_status)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "sv_load by hr_status", x = "hr_status", y = "sv_load") +
  scale_color_manual(name = "hr_status", 
                     values = c("HR_proficient" = "blue", 
                                "HR_deficient" = "red",
                                "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  scale_fill_manual(name = "hr_status", values = c("HR_proficient" = "lightblue", 
                                                   "HR_deficient" = "lightcoral", 
                                                   "cannot_be_determined" = "brown")
                     ) +  # Adjust colors as needed
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  )

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(p_sv_load)



######################################################################################
######################################################################################
# we are doing the same analysis function of HRD_TYPE
######################################################################################
######################################################################################
# Working with HRD.TYPE

colnames(df.hrd_type)
df.hrd_type <- df.hrd_type[complete.cases(df.hrd_type$hrd_type), ]

p.hrd_type <- df.hrd_type %>%
  gather(key = "variable", value = "value", -hrd_type) %>%
  ggplot(aes(x = hrd_type, y = value, fill = hrd_type)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "HR STATUS", x = "hrd_type", y = "Value") +
  theme_minimal()


summary(df.hrd_type)


######################################################################################
######################################################################################
# display hrd_type function of p_BRCA1
######################################################################################
######################################################################################


summary_stats <- tapply(df.hrd_type$p_BRCA1, df.hrd_type$hrd_type, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_BRCA1 <- ggplot(df.hrd_type, aes(x = hrd_type, y = p_BRCA1, fill = hrd_type, color = hrd_type)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "p_BRCA1 by hrd_type", x = "hrd_type", y = "p_BRCA1") +
  scale_color_manual(name = "hrd_type", 
                     values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                               "none" = "black")
                     ) +  
  scale_fill_manual(name = "hrd_type", 
                    values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                                "none" = "black")
                    ) +         
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  )

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(p_BRCA1)


######################################################################################
######################################################################################
# display hrd_type function of p_BRCA2
######################################################################################
######################################################################################

summary_stats <- tapply(df.hrd_type$p_BRCA2, df.hrd_type$hrd_type, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_BRCA2 <- ggplot(df.hrd_type, aes(x = hrd_type, y = p_BRCA2, fill = hrd_type, color = hrd_type)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "p_BRCA2 by hrd_type", x = "hrd_type", y = "p_BRCA2") +
  scale_color_manual(name = "hrd_type", 
                     values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                               "none" = "black")
                     ) +  
  scale_fill_manual(name = "hrd_type", 
                    values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                                "none" = "black")
                    ) +         
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  )

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(p_BRCA2)

######################################################################################
######################################################################################
# display hrd_type function of p_hrd
######################################################################################
######################################################################################


summary_stats <- tapply(df.hrd_type$p_hrd, df.hrd_type$hrd_type, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_hrd <- ggplot(df.hrd_type, aes(x = hrd_type, y = p_hrd, fill = hrd_type, color = hrd_type)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "p_hrd by hrd_type", x = "hrd_type", y = "p_hrd") +
  scale_color_manual(name = "hrd_type", 
                     values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                               "none" = "black")
                     ) +  
  scale_fill_manual(name = "hrd_type", 
                    values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                                "none" = "black")
                    ) +         
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  )

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(p_hrd)

######################################################################################
######################################################################################
# display hrd_type function of indel_load
######################################################################################
######################################################################################

summary_stats <- tapply(df.hrd_type$indel_load, df.hrd_type$hrd_type, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


indel_load <- ggplot(df.hrd_type, aes(x = hrd_type, y = indel_load, fill = hrd_type, color = hrd_type)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "indel_load by hrd_type", x = "hrd_type", y = "indel_load") +
  scale_color_manual(name = "hrd_type", 
                     values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                               "none" = "black")
                     ) +  
  scale_fill_manual(name = "hrd_type", 
                    values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                                "none" = "black")
                    ) +         
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  ) + 
  ylim(c(0,6000))

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(indel_load)

######################################################################################
######################################################################################
# display hrd_type function of indel_rep_load
######################################################################################
######################################################################################

summary_stats <- tapply(df.hrd_type$indel_rep_load, df.hrd_type$hrd_type, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


indel_rep_load <- ggplot(df.hrd_type, aes(x = hrd_type, y = indel_rep_load, fill = hrd_type, color = hrd_type)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "indel_rep_load by hrd_type", x = "hrd_type", y = "indel_rep_load") +
  scale_color_manual(name = "hrd_type", 
                     values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                               "none" = "black")
                     ) +  
  scale_fill_manual(name = "hrd_type", 
                    values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                                "none" = "black")
                    ) +         
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  )  +
  ylim(c(0,4000))

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(indel_rep_load)


######################################################################################
######################################################################################
# display hrd_type function of sv_load
######################################################################################
######################################################################################

summary_stats <- tapply(df.hrd_type$sv_load, df.hrd_type$hrd_type, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


sv_load <- ggplot(df.hrd_type, aes(x = hrd_type, y = sv_load, fill = hrd_type, color = hrd_type)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "sv_load by hrd_type", x = "hrd_type", y = "sv_load") +
  scale_color_manual(name = "hrd_type", 
                     values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                               "none" = "black")
                     ) +  
  scale_fill_manual(name = "hrd_type", 
                    values = c("BRCA1_type" = "orange", 
                                "BRCA2_type" = "red",
                                "cannot_be_determined" = "brown", 
                                "none" = "black")
                    ) +         
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  )  +
  ylim(c(0,2000))

# Print the ggplot figure
options(repr.plot.width=10, repr.plot.height=8)
print(sv_load)

cor(df.hrd_type$p_BRCA1, df.hrd_type$p_BRCA2)
cor(df.hrd_type$p_BRCA1, df.hrd_type$p_hrd)
cor(df.hrd_type$p_BRCA2, df.hrd_type$p_hrd)

colnames(df.hrd_type)


######################################################################################
######################################################################################
# Working with CANCER_TYPE : df.cancer_type
######################################################################################
######################################################################################


# Working with CANCER.TYPE

colnames(df.cancer_type)
df.cancer_type <- df.cancer_type[complete.cases(df.cancer_type$cancer_type), ]

p.cancer_type <- df.cancer_type %>%
  gather(key = "variable", value = "value", -cancer_type) %>%
  ggplot(aes(x = cancer_type, y = value, fill = cancer_type)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "HR STATUS", x = "cancer_type", y = "Value") +
  theme_minimal()


table(df.cancer_type$cancer_type)
summary(df.cancer_type)

df.cancer_type$cancer_type = factor(df.cancer_type$cancer_type)

# the number of CANCER TYPES
length(unique(df.cancer_type$cancer_type))

# the CANCER TYPES

unique(df.cancer_type$cancer_type)

# to generate random colors for display

NUMBER_CANCER_TYPES = length(unique(df.cancer_type$cancer_type))
NUMBER_CANCER_TYPES

random_colors <- sample(colors(), NUMBER_CANCER_TYPES)
random_colors

# display cancer_type function of p_BRCA1
# options(repr.plot.width=20, repr.plot.height=40)

summary_stats <- tapply(df.cancer_type$p_BRCA1, df.cancer_type$cancer_type, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_BRCA1 <- ggplot(df.cancer_type, aes(x = cancer_type, y = p_BRCA1, fill = cancer_type, color = cancer_type)) +
  geom_boxplot(notch = FALSE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "p_BRCA1 by cancer_type", x = "cancer_type", y = "p_BRCA1") +
  scale_color_manual(name = "cancer_type", 
                                        values = setNames(random_colors, 
                                        unique(df.cancer_type$cancer_type))
                     ) + 
  scale_fill_manual(name = "cancer_type", 
                    values = setNames(random_colors, 
                                      unique(df.cancer_type$cancer_type))  
                    ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  ) + 
  coord_flip() 

# Print the ggplot figure
options(repr.plot.width=18, repr.plot.height=10)
print(p_BRCA1)


######################################################################################
######################################################################################
# display cancer_type function of p_BRCA2
# options(repr.plot.width=20, repr.plot.height=40)
######################################################################################
######################################################################################

summary_stats <- tapply(df.cancer_type$p_BRCA2, df.cancer_type$cancer_type, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_BRCA2 <- ggplot(df.cancer_type, aes(x = cancer_type, y = p_BRCA2, fill = cancer_type, color = cancer_type)) +
  geom_boxplot(notch = FALSE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "p_BRCA2 by cancer_type", x = "cancer_type", y = "p_BRCA2") +
  scale_color_manual(name = "cancer_type", 
                                        values = setNames(random_colors, 
                                        unique(df.cancer_type$cancer_type))
                     ) + 
  scale_fill_manual(name = "cancer_type", 
                    values = setNames(random_colors, 
                                      unique(df.cancer_type$cancer_type))  
                    ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  ) + 
  coord_flip() 

# Print the ggplot figure
options(repr.plot.width=18, repr.plot.height=10)
print(p_BRCA2)


# display cancer_type function of p_hrd
# options(repr.plot.width=20, repr.plot.height=40)

summary_stats <- tapply(df.cancer_type$p_hrd, df.cancer_type$cancer_type, function(x) {
  c(
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
})

# Print the summary statistics
print(data.frame(do.call(rbind, summary_stats)))


p_hrd <- ggplot(df.cancer_type, aes(x = cancer_type, y = p_hrd, fill = cancer_type, color = cancer_type)) +
  geom_boxplot(notch = FALSE, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.7, size = 3) +
  labs(title = "p_hrd by cancer_type", x = "cancer_type", y = "p_hrd") +
  scale_color_manual(name = "cancer_type", 
                                        values = setNames(random_colors, 
                                        unique(df.cancer_type$cancer_type))
                     ) + 
  scale_fill_manual(name = "cancer_type", 
                    values = setNames(random_colors, 
                                      unique(df.cancer_type$cancer_type))  
                    ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Adjust plot title size and style 
    legend.text = element_text(size = 14)  # Adjust the text of the legend  
  ) + 
  coord_flip() 

# Print the ggplot figure
options(repr.plot.width=18, repr.plot.height=10)
print(p_hrd)


######################################################################################
######################################################################################

# To summarizesome of the findings on HR_STATUS

# pBRCA1

#                     Median Min   Max Q1.25. Q3.75.
# cannot_be_determined  0.000   0 0.458   0.00  0.002
# HR_deficient          0.104   0 0.998   0.02  0.762
# HR_proficient         0.000   0 0.454   0.00  0.000


# pBRCA2
#                     Median Min   Max Q1.25. Q3.75.
# cannot_be_determined  0.002   0 0.546  0.000 0.0135
# HR_deficient          0.664   0 0.998  0.152 0.9280
# HR_proficient         0.000   0 0.486  0.000 0.0040

# probab_HRD
#
#                     Median   Min   Max Q1.25. Q3.75.
# cannot_be_determined  0.002 0.000 0.632  0.000  0.020
# HR_deficient          0.954 0.504 1.000  0.868  0.978
# HR_proficient         0.000 0.000 0.498  0.000  0.006

######################################################################################
######################################################################################