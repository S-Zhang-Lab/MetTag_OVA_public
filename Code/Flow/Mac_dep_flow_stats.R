# --- 1. LOAD LIBRARIES ---
library(ggplot2)
library(ggpubr)

# --- 2. DEFINE DATA ---
# n=6 Lipo (Control)
# n=7 Clod (Depletion)
# Note: One sample (13.6mg weight) was excluded a priori as a seeding failure.
lipo_flow  <- c(11.10, 3.45, 4.72, 8.39, 3.07, 5.70) 
chlod_flow <- c(3.26, 0.22, 0.40, 3.40, 3.80, 4.79, 4.84)

# --- 3. CREATE RIGOROUS DATAFRAME ---
# This ensures ggplot sees exactly 13 rows.
df_flow <- data.frame(
  Group = factor(c(rep("Lipo", length(lipo_flow)), 
                   rep("Clod", length(chlod_flow))),
                 levels = c("Lipo", "Clod")), # Sets order on X-axis
  Frequency = c(lipo_flow, chlod_flow)
)

# --- 4. STATISTICAL VALIDATION (Console Output) ---
# Switching to "two.sided" to be conservative and avoid "p-hacking" concerns.
t_test_results <- t.test(Frequency ~ Group, data = df_flow, 
                         alternative = "two.sided", 
                         var.equal = FALSE) # Welch's correction for unequal n

print(t_test_results)

# Calculation of Biological Magnitude
mean_lipo <- mean(lipo_flow)
mean_clod <- mean(chlod_flow)
percent_reduction <- (1 - (mean_clod / mean_lipo)) * 100
cat("Percent Reduction in Macrophages:", round(percent_reduction, 1), "%\n")

# --- 5. GENERATE PLOT ---
plot_flow <- ggplot(df_flow, aes(x = Group, y = Frequency, fill = Group)) +
  # Boxplot: outlier.shape = NA prevents "double plotting" of dots
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.4, color = "black", lwd = 0.6) +
  # Jitter: width=0.2 ensures the 7th dot isn't hiding behind another one
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.8, shape = 21, color = "black") + 
  theme_pubr() +
  scale_fill_manual(values = c("Lipo" = "grey80", "Clod" = "#E41A1C")) +
  labs(y = "Macrophage Frequency (%)", x = NULL) +
  # Add the two-tailed p-value directly to the plot
  stat_compare_means(method = "t.test", 
                     method.args = list(var.equal = FALSE, alternative = "two.sided"), 
                     label = "p.format",
                     label.x = 1.35) + # Adjust label position
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold"))

# Display
print(plot_flow)

# --- 6. SAVE ---
ggsave("Fig_S5_Macrophage_Depletion_TwoTailed.pdf", plot_flow, width = 3.5, height = 4.5)

# Extract p-value and export stats to CSV
p_val <- t_test_results$p.value

stats_export <- data.frame(
  Test = "Welch Two Sample t-test",
  Alternative = "two.sided",
  Mean_Lipo = mean_lipo,
  Mean_Clod = mean_clod,
  Percent_Reduction = percent_reduction,
  P_Value = p_val
)

write.csv(stats_export, "Suppl_6B_Macrophage_Stats.csv", row.names = FALSE)
print(paste("Exact p-value is:", p_val))
