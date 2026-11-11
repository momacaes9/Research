# AGRO-BHARAT — R version of the statistical analysis
# Two-way mixed ANOVA (between: group, within: week) + linear mixed model.
# Run from the repo root:  Rscript analysis/analysis.R
#
# Packages: install.packages(c("afex", "emmeans", "lme4", "lmerTest", "ggplot2"))

suppressPackageStartupMessages({
  library(afex)      # aov_ez: mixed ANOVA with Greenhouse-Geisser correction
  library(emmeans)   # post-hoc contrasts
  library(lmerTest)  # lmer with p-values
  library(ggplot2)
})

df <- read.csv("data/plant_height_long.csv")
df$group    <- factor(df$group, levels = c("control", "neem"))
df$week_f   <- factor(df$week)
df$plant_id <- factor(df$plant_id)

# 1. Descriptive statistics ----------------------------------------------------
desc <- aggregate(height_cm ~ week + group, df,
                  function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
print(do.call(data.frame, desc))

# 2. Mixed ANOVA ----------------------------------------------------------------
fit <- aov_ez(id = "plant_id", dv = "height_cm", data = df,
              between = "group", within = "week_f",
              anova_table = list(es = "pes", correction = "GG"))
print(fit)
# Expected (uncorrected): group F(1,18)=6.38 p=.021; week F(9,162)=259.09;
# group x week F(9,162)=3.38 p=.001 (GG-corrected p = .04)
summary(fit)  # includes Mauchly's test

# 3. Post-hoc: Neem vs Control at each week (Holm) ------------------------------
emm <- emmeans(fit, ~ group | week_f)
# by = NULL pools all 10 weekly contrasts so Holm corrects across weeks
print(test(pairs(emm), by = NULL, adjust = "holm"))

# 4. Linear mixed model: growth rate per group ----------------------------------
df$week_c <- df$week - 1
lmm <- lmer(height_cm ~ week_c * group + (week_c | plant_id), data = df, REML = TRUE)
print(summary(lmm))
# week_c            ~0.390 cm/week (control)
# week_c:groupneem  ~0.084 cm/week extra for Neem
# (Satterthwaite df = 18: p = .051; Python statsmodels Wald z-test: p = .036)

# 5. Figure ---------------------------------------------------------------------
p <- ggplot(df, aes(week, height_cm, colour = group)) +
  geom_line(aes(group = plant_id), alpha = 0.2) +
  stat_summary(fun.data = mean_cl_normal, geom = "ribbon", alpha = 0.2,
               aes(fill = group), colour = NA) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  scale_colour_manual(values = c(control = "#C98A1B", neem = "#2B7A9E")) +
  scale_fill_manual(values = c(control = "#C98A1B", neem = "#2B7A9E")) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Week", y = "Stem height (cm)",
       title = "Plant growth curves by treatment (mean ± 95% CI)") +
  theme_minimal()
ggsave("results/figures/growth_curves_R.png", p, width = 7, height = 4.5, dpi = 200)
