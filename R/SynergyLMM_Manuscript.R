# Code to reproduce the results shown in 'SynergyLMM: A comprehensive statistical framework and interactive web-tool for designing and analyzing in vivo drug combination experiments'

# Install the last version of 'SynergyLMM'

install.packages("SynergyLMM") # version 1.1.1

# Install any other required package not installed:

install.packages("tidyverse")
install.packages("broom")
install.packages("nlme")
install.packages("cowplot")
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("maobinchen/invivoSyn")
install.packages("openxlsx")
install.packages("ggbreak")


# Load Libraries ----

library(SynergyLMM)
library(tidyverse)
library(broom)
library(nlme)
library(cowplot)
library(invivoSyn)
library(openxlsx)
library(ggbreak)

# Figures dimensions (A4)
height <- 11.69
width <- 8.27

# Figure 2 ----

## Fig. 2a ----

U87MG <- read.xlsx("data/Fig2/U87MG.xlsx")

U87MG$Day <- sqrt(U87MG$Day-12)

# Fit model
U87MG_lmm <- lmmModel(data = U87MG, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                      trt_control = "Control", drug_a = "Docetaxel", drug_b = "GNE-317", combination = "Docetaxel+GNE-317",
                      weights = varPower(form = ~Time))

# Generate plots
trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

df <- U87MG_lmm$dt1

a_left <- df %>% ggplot(aes(x = Time, y = RTV, group = SampleID)) + geom_point(aes(colour = Treatment), size = 3) +
  geom_line(aes(colour = Treatment)) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top") +
  labs(title = "U-87-MG FM") +
  ylab("Relative Luminescence Units") + xlab("Days after treatment initiation") +
  scale_x_continuous(breaks = unique(U87MG_lmm$dt1$Time), labels = unique(U87MG_lmm$dt1$Time)^2)
a_left

a_right <- plot_lmmModel(U87MG_lmm, trt_control = "Control", drug_a = "Docetaxel", drug_b = "GNE-317", combination = "Docetaxel+GNE-317") +
  theme(legend.position = "top") +
  labs(title = "") +
  scale_x_continuous(breaks = unique(U87MG_lmm$dt1$Time), labels = unique(U87MG_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") +
  xlab("Days after treatment initiation") + ylab("log (RLU)")
a_right

## Fig. 2b ----

# Calculate synergy
(bliss <- lmmSynergy(U87MG_lmm, robust = F, type = "CR2", min_time = 0, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

(hsa <- lmmSynergy(U87MG_lmm, robust = F, type = "CR2", min_time = 0, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2

# Generate plots

bliss_CI <- bliss$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
bliss_CI$Model <- "SynergyLMM Bliss"
hsa_CI <- hsa$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
hsa_CI$Model <- "SynergyLMM HSA"

Narayan <- data.frame(Model = "Narayan et al", Estimate = c(0.8, 0.56), lwr = NA, upr = NA, Time = c(9, 16), pval = NA)

CI_df <- rbind(Narayan, bliss_CI, hsa_CI)
CI_df$Time <- as.factor(CI_df$Time)

b <- CI_df %>% ggplot(aes(x = Time, y = Estimate)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) + cowplot::theme_cowplot() +
  geom_point(aes(fill  = pval), size = 5, shape = 23, color = "gray60") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", na.value = "firebrick2") +
  ylab("Combination Index") + xlab("Days after treatment initiation") + 
  scale_y_continuous(breaks = c(0 , 0.5 , 0.8, 1, 1.5, 2)) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 0.8, lty = 3) + 
  facet_wrap(~Model) + 
  theme(strip.background = element_rect(fill = "cyan4"), 
        strip.text = element_text(color = "white", face = "bold", size = 16)) +
  labs(title = "U-87-MG FM Docetaxel - GNE-317")
b


## Fig. 2c ----

BV173Gluc <- read.csv("data/Fig2/BV_173.csv")

BV173Gluc$Day <- sqrt(BV173Gluc$Day-7)

# Fit Model
BV173Gluc_lmm <- lmmModel(data = BV173Gluc, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                          trt_control = "Control", drug_a = "Imatinib", drug_b = "Dasatinib", combination = "Imatinib_Dasatinib",
                          weights = nlme::varIdent(form = ~ 1|Treatment))


trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

df <- BV173Gluc_lmm$dt1

c_left <- df %>% ggplot(aes(x = Time, y = RTV, group = SampleID)) + geom_point(aes(colour = Treatment), size = 3) +
  geom_line(aes(colour = Treatment)) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top") +
  labs(title = "BV-173-Gluc") +
  ylab("Relative Luminescence Units (RLU)") + xlab("Days after treatment initiation") +
  scale_x_continuous(breaks = unique(BV173Gluc_lmm$dt1$Time), labels = unique(BV173Gluc_lmm$dt1$Time)^2)
c_left

c_right <- plot_lmmModel(BV173Gluc_lmm, trt_control = "Control", drug_a = "Imatinib", drug_b = "Dasatinib", combination = "Imatinib_Dasatinib") +
  theme(legend.position = "top") +
  labs(title = "")  +
  scale_x_continuous(breaks = unique(BV173Gluc_lmm$dt1$Time), labels = unique(BV173Gluc_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") +
  xlab("Days after treatment initiation") + ylab("log (RLU)")
c_right


## Fig. 2d ----

# Calculate Synergy
(bliss <- lmmSynergy(BV173Gluc_lmm, robust = T, type = "CR2", min_time = 3, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

(hsa <- lmmSynergy(BV173Gluc_lmm, robust = T, type = "CR2", min_time = 3, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2

# Generate Plot

bliss_CI <- bliss$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
bliss_CI$Model <- "SynergyLMM Bliss"
hsa_CI <- hsa$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
hsa_CI$Model <- "SynergyLMM HSA"

Narayan <- data.frame(Model = "Narayan et al", Estimate = c(0.38, 0.24), lwr = NA, upr = NA, Time = c(21-7, 28-7), pval = NA)

CI_df <- rbind(Narayan, bliss_CI, hsa_CI)
CI_df$Time <- as.factor(CI_df$Time)

d <- CI_df %>% ggplot(aes(x = Time, y = Estimate)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) + cowplot::theme_cowplot() +
  geom_point(aes(fill  = pval), size = 5, shape = 23, color = "gray60") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", na.value = "firebrick2") +
  ylab("Combination Index") + xlab("Days after treatment initiation") +
  #scale_y_continuous(breaks = c(0 , 0.8 , 1, seq(2, 20, 2))) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 0.8, lty = 3) + 
  facet_wrap(~Model) + 
  theme(strip.background = element_rect(fill = "cyan4"), 
        strip.text = element_text(color = "white", face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "BV-173-Gluc Imatinib - Dasatinib") + ggbreak::scale_y_break(breaks = c(7.5,10), scales = 0.25)
d

## Fig. 2e ----

CHL1FM <- read.xlsx("data/Fig2/CHL1_FM.xlsx")

CHL1FM$Day <- sqrt(CHL1FM$Day-7)

# Fit model
CHL1FM_lmm <- lmmModel(data = CHL1FM, grwth_model = "gompertz", sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                       trt_control = "Control", drug_a = "Gemcitabine", drug_b = "CGP_082996", combination = "Gemcitabine_CGP",
                       weights = nlme::varIdent(form = ~1|SampleID),
                       control = lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000),
                       start_values = "selfStart")

# Generate plots
trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

df <- CHL1FM_lmm$dt1

e_left <- df %>% ggplot(aes(x = Time, y = RTV, group = SampleID)) + geom_point(aes(colour = Treatment), size = 3) +
  geom_line(aes(colour = Treatment)) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "CHL1-FM") +
  ylab("Relative Luminescence Units (RLU)") + xlab("Days after treatment initiation") +
  scale_x_continuous(breaks = unique(CHL1FM_lmm$dt1$Time), labels = unique(CHL1FM_lmm$dt1$Time)^2)
e_left

e_right <- plot_lmmModel(CHL1FM_lmm, trt_control = "Control", drug_a = "Gemcitabine", drug_b = "CGP_082996", combination = "Gemcitabine_CGP") +
  theme(legend.position = "top") +
  labs(title = "") +
  scale_x_continuous(breaks = unique(CHL1FM_lmm$dt1$Time), labels = unique(CHL1FM_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") +
  xlab("Days after treatment initiation") + ylab("log (RLU)")
e_right

## Fig. 2f ----

# Calculate Synergy
(bliss <- lmmSynergy(CHL1FM_lmm, robust = T, type = "CR2", min_time = 0, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

(hsa <- lmmSynergy(CHL1FM_lmm, robust = T, type = "CR2", min_time = 0, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2

# Generate plots
bliss_CI <- bliss$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
bliss_CI$Model <- "SynergyLMM Bliss"
hsa_CI <- hsa$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
hsa_CI$Model <- "SynergyLMM HSA"

Narayan <- data.frame(Model = "Narayan et al", Estimate = c(0.68, 0.62), lwr = NA, upr = NA, Time = c(21-7, 28-7), pval = NA)

CI_df <- rbind(Narayan, bliss_CI, hsa_CI)
CI_df$Time <- as.factor(CI_df$Time)

f <- CI_df %>% ggplot(aes(x = Time, y = Estimate)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) + cowplot::theme_cowplot() +
  geom_point(aes(fill  = pval), size = 5, shape = 23, color = "gray60") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", na.value = "firebrick2") +
  ylab("Combination Index") + xlab("Days after treatment initiation") + 
  ggbreak::scale_y_break(breaks = c(50,75), scales = "free") +
  #scale_y_continuous(breaks = c(0 , 0.8 , 1, seq(2.5, 12.5, 2.5))) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 0.8, lty = 3) + 
  facet_wrap(~Model) + 
  theme(strip.background = element_rect(fill = "cyan4"), 
        strip.text = element_text(color = "white", face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "CHL-1 FM CGP-082996 - Gemcitabine")
f

## Fig. 2g ----

MDAMD231 <- read.xlsx("data/Fig2/MDA_MD231_Fig5C.xlsx")

MDAMD231$Day <- sqrt(MDAMD231$Day-16) 

# Fit model

MDAMD231_lmm <- lmmModel(data = MDAMD231, grwth_model = "gompertz", sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                         trt_control = "Control", drug_a = "AZ628", drug_b = "Gemcitabine", combination = "AZ628_Gemcitabine",
                         weights = varIdent(form = ~ 1|Treatment))


trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

df <- MDAMD231_lmm$dt1

g_left <- df %>% ggplot(aes(x = Time, y = RTV, group = SampleID)) + geom_point(aes(colour = Treatment), size = 3) +
  geom_line(aes(colour = Treatment)) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MDA-MB-231 FM")+
  ylab("Relative Luminescence Units (RLU)") + xlab("Days after treatment initiation") +
  scale_x_continuous(breaks = unique(MDAMD231_lmm$dt1$Time), labels = unique(MDAMD231_lmm$dt1$Time)^2)
g_left

g_right <- plot_lmmModel(MDAMD231_lmm, trt_control = "Control", drug_a = "AZ628", drug_b = "Gemcitabine", combination = "AZ628_Gemcitabine") +
  theme(legend.position = "top") +
  labs(title = "") +
  scale_x_continuous(breaks = unique(MDAMD231_lmm$dt1$Time), labels = unique(MDAMD231_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") +  xlab("Days after treatment initiation") + ylab("log (RLU)")
g_right


## Fig. 2h ----

# Bliss synergy calculation
(bliss <- lmmSynergy(MDAMD231_lmm, robust = T, type = "CR2", min_time = 0, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

# HSA synergy calculation
(hsa <- lmmSynergy(MDAMD231_lmm, robust = T, type = "CR2", min_time = 0, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2


bliss_CI <- bliss$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
bliss_CI$Model <- "SynergyLMM Bliss"
hsa_CI <- hsa$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
hsa_CI$Model <- "SynergyLMM HSA"

# Generate Plot

Narayan <- data.frame(Model = "Narayan et al", Estimate = c(0.08, 0.11), lwr = NA, upr = NA, Time = c(15, 21), pval = NA)

CI_df <- rbind(Narayan, bliss_CI, hsa_CI)
CI_df$Time <- as.factor(CI_df$Time)

h <- CI_df %>% ggplot(aes(x = Time, y = Estimate)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) + cowplot::theme_cowplot() +
  geom_point(aes(fill  = pval), size = 5, shape = 23, color = "gray60") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", na.value = "firebrick2") +
  ylab("Combination Index") + xlab("Days after treatment initiation") +
  scale_y_continuous(breaks = c(0 , 0.2 , 0.4, 0.6, 0.8, 1)) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 0.8, lty = 3) + 
  facet_wrap(~Model) + 
  theme(strip.background = element_rect(fill = "cyan4"), 
        strip.text = element_text(color = "white", face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MDA-MB-231 FM AZD628 - Gemcitabine")
h

ab <- plot_grid(plot_grid(a_left, a_right), b, rel_widths = c(1,1), labels = "auto", label_size = 18)
cd <- plot_grid(plot_grid(c_left, c_right), d, rel_widths = c(1,1), labels = c("c", "d"), label_size = 18)
ef <- plot_grid(plot_grid(e_left, e_right), f, rel_widths = c(1,1), labels = c("e", "f"), label_size = 18)
gh <- plot_grid(plot_grid(g_left, g_right), h, rel_widths = c(1,1), labels = c("g", "h"), label_size = 18)

Fig2 <- plot_grid(ab, cd, ef, gh, ncol = 1)
Fig2

#ggsave("figures/Figure2.pdf", Fig2, height = 2.25*height, width = 3.9*width)

# Supplementary Figure 1 ----

## Supplementary Fig. 1a ----
U87MG_lmm_ranef <- ranefDiagnostics(U87MG_lmm)
U87MG_lmm_resid <- residDiagnostics(U87MG_lmm)

a <- cowplot::plot_grid(plotlist = c(U87MG_lmm_ranef$Plots[1:4],U87MG_lmm_resid$Plots[1:5]))
a


## Supplementary Fig. 1b ----

BV173Gluc_lmm_ranef <- ranefDiagnostics(BV173Gluc_lmm)
BV173Gluc_lmm_resid <- residDiagnostics(BV173Gluc_lmm)

b <- cowplot::plot_grid(plotlist = c(BV173Gluc_lmm_ranef$Plots[1:4],BV173Gluc_lmm_resid$Plots[1:5]))
b


## Supplementary Fig. 1c ----

CHL1FM_lmm_ranef <- ranefDiagnostics(CHL1FM_lmm)
CHL1FM_lmm_resid <- residDiagnostics(CHL1FM_lmm)

c <- cowplot::plot_grid(plotlist = c(CHL1FM_lmm_ranef$Plots[1:5],CHL1FM_lmm_resid$Plots[1:5]), ncol = 5, nrow = 2)
c

## Supplementary Fig. 1d ----

MDAMD231_lmm_ranef <- ranefDiagnostics(MDAMD231_lmm)
MDAMD231_lmm_resid <- residDiagnostics(MDAMD231_lmm)

d <- cowplot::plot_grid(plotlist = c(MDAMD231_lmm_ranef$Plots[1:5],MDAMD231_lmm_resid$Plots[1:5]), ncol = 5, nrow = 2)
d


SuppFig1 <- cowplot::plot_grid(a,b,c,d, labels = "auto")
SuppFig1

#ggsave(filename = "figures/SuppFig1.pdf", plot = SuppFig1, width = 6*width, height = 4.5*height, limitsize = FALSE)

# Supplementary Figure 2 ----

## Supplementary Fig. 2a ----

ObsvsPred(U87MG_lmm, 3,4)
a <- plot_SynergyLMM(U87MG_lmm, plot_type = "ObsvsPred", nrow = 3, ncol = 4)
a

## Supplementary Fig. 2b ----

ObsvsPred(BV173Gluc_lmm, 4,7)
b <- plot_SynergyLMM(BV173Gluc_lmm, plot_type = "ObsvsPred", nrow = 4, ncol = 7)
b

## Supplementary Fig. 2c ----

ObsvsPred(CHL1FM_lmm, 5,5)
c <- plot_SynergyLMM(CHL1FM_lmm, plot_type = "ObsvsPred", nrow = 5, ncol = 5)
c

## Supplementary Fig. 2d ----

ObsvsPred(MDAMD231_lmm, 6,5)
d <- plot_SynergyLMM(MDAMD231_lmm, plot_type = "ObsvsPred", 6,5)
d


## Supplementary Fig. 2e ----

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")


U87MG_0 <- U87MG_lmm$dt1 |> filter(Time == 0)

U87MG_TV0 <- U87MG_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5) +
  labs(title = "U-87-MG-FM") +
  ylab("Initial Relative Luminescence Units (RLU)")

BV173Gluc_0 <- BV173Gluc_lmm$dt1 |> filter(Time == 0)

BV173Gluc_TV0 <- BV173Gluc_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5) +
  labs(title = "BV-173-Gluc")+
  ylab("Initial Relative Luminescence Units (RLU)")


CHL1FM_0 <- CHL1FM_lmm$dt1 |> filter(Time == 0)

my_comparisons <- list(c("Control", "Gemcitabine"), c("Control", "CGP_082996"), c("Gemcitabine", "CGP_082996"), c("CGP_082996", "Gemcitabine_CGP"))

stat.test <- ggpubr::compare_means(TV ~ Treatment, data = CHL1FM_0, method = "wilcox.test", p.adjust.method = "BH")

stat.test$y.position <- c(4.5e7, 0,0, 5e7, 5.5e7, 6e7)

CHL1FM_TV0 <- CHL1FM_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5, label.y = 6e7) +
  ggpubr::stat_pvalue_manual(data = stat.test, hide.ns = T, label = "p={p.adj}") +
  labs(title = "CHL1-FM") +
  ylab("Initial Relative Luminescence Units (RLU)")

MDAMD231_0 <- MDAMD231_lmm$dt1 |> filter(Time == 0)

MDAMD231_TV0 <- MDAMD231_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = FALSE) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5) +
  labs(title = "MDA-MB-231 FM")+
  ylab("Initial Relative Luminescence Units (RLU)")

e <- cowplot::plot_grid(U87MG_TV0, BV173Gluc_TV0, CHL1FM_TV0, MDAMD231_TV0, nrow = 1, labels = "e")
e

top <- plot_grid(a, b, c, d, ncol = 2, labels = "auto")

SuppFig2 <- plot_grid(top, e, ncol = 1, rel_heights = c(2, 1))
SuppFig2

#ggsave(filename = "figures/SuppFig2.pdf", plot = SuppFig2, width = 3*width, height = 2.25*height, limitsize = FALSE)


# Supplementary Figure 3 ----

## Supplementary Fig. 3a ----

a <- CookDistance(MDAMD231_lmm, type = "fitted")
thra <- 3*mean(a[a !=0])

a_df <- data.frame(Subject = names(a), CookD = a)

a_plot <- a_df %>% ggplot(aes(x = Subject, y = CookD)) +
  geom_segment(aes(x = Subject, y = 0, yend = CookD), color = "slateblue") +
  annotate(geom = "point", x = a_df$Subject[a_df$CookD > thra], y = a_df$CookD[a_df$CookD > thra], color = "slateblue", size = 3) +
  annotate(geom = "text", x = a_df$Subject[a_df$CookD > thra], y = 1.05*a_df$CookD[a_df$CookD > thra],
           label = a_df$Subject[a_df$CookD > thra]) +
  geom_hline(aes(yintercept = thra, linetype = "")) +
  labs(title = "Cook's Distances based on Fitted Values Change") +
  ylab("Cook's Distance") +
  scale_linetype_manual(name = paste("Cook's Distance Threshold: ", round(thra,3)), values = 3) +
  theme_cowplot() + theme(legend.position = "top") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
a_plot
   

## Supplementary Fig. 3b ----

b <- CookDistance(MDAMD231_lmm, type = "fixef")
thr <- 3*mean(b[b !=0])

b_df <- data.frame(Subject = names(b), CookD = b)

b_plot <- b_df %>% ggplot(aes(x = Subject, y = CookD)) +
  geom_segment(aes(x = Subject, y = 0, yend = CookD), color = "#a83260") +
  annotate(geom = "point", x = b_df$Subject[b_df$CookD > thr], y = b_df$CookD[b_df$CookD > thr], color = "#a83260", size = 3) +
  annotate(geom = "text", x = b_df$Subject[b_df$CookD > thr], y = 1.05*b_df$CookD[b_df$CookD > thr],
           label = b_df$Subject[b_df$CookD > thr]) +
  geom_hline(aes(yintercept = thr, linetype = "")) +
  labs(title = "Cook's Distances based on Fixed Effects Values Change") +
  ylab("Cook's Distance") +
  scale_linetype_manual(name = paste("Cook's Distance Threshold: ", round(thr,3)), values = 3) +
  theme_cowplot() + theme(legend.position = "top") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
b_plot

## Supplementary Fig. 3c ----

df <- MDAMD231_lmm$dt1

df$Outliers <- NA

df$Outliers[df$SampleID == "7-4_Control"] <- "7-4_Control"
df$Outliers[!df$SampleID %in% c("7-4_Control")] <- ""

c_plot <- df %>% ggplot(aes(x = Time, y = logRTV, group = SampleID)) + geom_point(aes(color = Outliers, fill = Outliers, shape = Outliers), size = 2) +
  geom_line(aes(color = Outliers), alpha = 0.5) + scale_color_manual(values = c("gray","#d12856")) +
  scale_fill_manual(values = c("gray","#d12856")) +
  scale_shape_manual(values = c(21:24)) +
  facet_wrap(~Treatment) + theme_cowplot() +
  labs(title = "MDA-MB-231 FM") + ylab("log (Relative Luminiscense Units)") +
  theme(strip.text = element_text(size = 18), legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = unique(df$Time), labels = unique(df$Time)^2) +
  xlab("Days since treatment initiation")
c_plot

ab <- plot_grid(a_plot, b_plot, labels = "auto", ncol = 1)
SuppFig3 <- plot_grid(ab, c_plot, rel_widths = c(1,2), labels = c("", "c"))
SuppFig3

#ggsave("figures/SuppFig3.pdf", SuppFig3, height = 1.25*height, width = 2.5*width)


# Figure 3 ----

## Fig. 3b ----

MAS98.06 <- read.csv("data/Fig3/MAS98_06.csv")
colnames(MAS98.06) <- c("Cell", "SampleID","Treatment", "TV","Time", "Figure")
MAS98.06 <- getRTV(MAS98.06, time_start = 0)

MAS98.06$Treatment <- factor(MAS98.06$Treatment, levels = c("Ctrl", "Her", "EGF", "Fulv", "Trast", "Her_Fulv", "EGF_Fulv","Fulv_Trast","Her_Trast"))

hline <- data.frame(yintercept = 0)

b <- MAS98.06 %>% ggplot(aes(x = Time, y = logRTV, color = Treatment)) +
  geom_line(aes(group = SampleID), alpha = 0.33) + geom_point(aes(group = SampleID, shape = Treatment, fill = Treatment)) +
  ylab("Log (RTV)") + 
  xlab("Time since start of treatment") + 
  scale_x_continuous(breaks = c(0,7,14,21,28)) + 
  cowplot::theme_cowplot() + facet_wrap(~Treatment) +
  labs(title = "MAS98.06 ER-dependent PDX") +
  geom_hline(data = hline, aes(yintercept = yintercept), linetype = "dashed") +
  scale_color_manual(values = c("#3c3c3b", "#00a49c", "#fc8d62","#d50c52","#8c510a","#601580","#e78ac3", "#01665e","#386cb0")) +
  scale_fill_manual(values = c("#3c3c3b", "#00a49c", "#fc8d62","#d50c52","#8c510a","#601580","#e78ac3", "#01665e","#386cb0")) +
  scale_shape_manual(values = c(21:25,21:24)) +
  theme(strip.text = element_text(size = 18))
b

## Fig. 3c ----

MAS98.06_Her_Fulv <- read.xlsx("data/Fig3/MAS98.06_Her_Fulv.xlsx")
MAS98.06_Her_Fulv$Mouse <- paste(MAS98.06_Her_Fulv$Treatment, " (",MAS98.06_Her_Fulv$Mouse,")", sep = "")
MAS98.06_Her_Fulv <- MAS98.06_Her_Fulv[order(MAS98.06_Her_Fulv$Mouse),]

MAS98.06_Her_Fulv_lmm <- lmmModel(data = MAS98.06_Her_Fulv, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                  trt_control = "Ctrl", drug_a = "Fulv", drug_b = "Her", combination = "Her_Fulv", time_end = 28,
                                  weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

(Her_Fulv <- lmmSynergy(MAS98.06_Her_Fulv_lmm, robust = T, type = "CR2", min_time = 7, method = "HSA"))
Her_Fulv <- Her_Fulv$Synergy
Her_Fulv$Drug <- "Her_Fulv"

MAS98.06_EGF_Fulv <- read.xlsx("data/Fig3/MAS98.06_EGF_Fulv.xlsx")
MAS98.06_EGF_Fulv$Mouse <- paste(MAS98.06_EGF_Fulv$Treatment, " (",MAS98.06_EGF_Fulv$Mouse,")", sep = "")
MAS98.06_EGF_Fulv <- MAS98.06_EGF_Fulv[order(MAS98.06_EGF_Fulv$Mouse),]


MAS98.06_EGF_Fulv_lmm <- lmmModel(data = MAS98.06_EGF_Fulv, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                  trt_control = "Ctrl", drug_a = "Fulv", drug_b = "EGF", combination = "EGF_Fulv", time_end = 28,
                                  weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

(EGF_Fulv <- lmmSynergy(MAS98.06_EGF_Fulv_lmm, robust = T, type = "CR2", min_time = 7, method = "HSA"))
EGF_Fulv <- EGF_Fulv$Synergy

EGF_Fulv$Drug <- "EGF_Fulv"

MAS98.06_Fulv_Trast <- read.xlsx("data/Fig3/MAS98.06_Fulv_Trast.xlsx")
MAS98.06_Fulv_Trast$Mouse <- paste(MAS98.06_Fulv_Trast$Treatment, " (",MAS98.06_Fulv_Trast$Mouse,")", sep = "")
MAS98.06_Fulv_Trast <- MAS98.06_Fulv_Trast[order(MAS98.06_Fulv_Trast$Mouse),]

MAS98.06_Fulv_Trast_lmm <- lmmModel(data = MAS98.06_Fulv_Trast, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                    trt_control = "Ctrl", drug_a = "Fulv", drug_b = "Trast", combination = "Fulv_Trast", time_end = 28,
                                    weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

(Fulv_Trast <- lmmSynergy(MAS98.06_Fulv_Trast_lmm, robust = T, type = "CR2", min_time = 8, method = "HSA"))
Fulv_Trast <- Fulv_Trast$Synergy
Fulv_Trast$Drug <- "Fulv_Trast"

MAS98.06_Her_Trast <- read.xlsx("data/Fig3/MAS98.06_Her_Trast.xlsx")
MAS98.06_Her_Trast$Sample <- paste(MAS98.06_Her_Trast$Treatment, " (",MAS98.06_Her_Trast$Sample,")", sep = "")
MAS98.06_Her_Trast <- MAS98.06_Her_Trast[order(MAS98.06_Her_Trast$Sample),]

MAS98.06_Her_Trast_lmm <- lmmModel(data = MAS98.06_Her_Trast, sample_id = "Sample", time = "Days", treatment = "Treatment", tumor_vol = "Total_vol",
                                   trt_control = "Ctrl", drug_a = "Her", drug_b = "Trast", combination = "Her_Trast", time_start = 0, time_end = 28,
                                   weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

(Her_Trast <- lmmSynergy(MAS98.06_Her_Trast_lmm, robust = T, type = "CR2", min_time = 7, method = "HSA"))
Her_Trast <- Her_Trast$Synergy
Her_Trast$Drug <- "Her_Trast"

HSA_Syn <- rbind(Her_Fulv, EGF_Fulv, Fulv_Trast, Her_Trast)

HSA_Syn <- HSA_Syn[,-1]

HSA_Syn <- HSA_Syn %>% filter(Time == 28)

HSA_Syn$Drug <- factor(HSA_Syn$Drug, levels = unique(HSA_Syn$Drug))

HSA_CI <- HSA_Syn %>% filter(Metric == "CI")
HSA_SS <- HSA_Syn %>% filter(Metric == "SS")


CI <- HSA_CI %>% ggplot(aes(x = Drug, y = Estimate, group = Drug)) +
  geom_segment(aes(x= Drug, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = pval), size = 5, color = "gray65", shape = 23) +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec") +
  geom_hline(yintercept = 1, lty = "dashed") +
  xlab("Drugs") +
  ylab("Combination Index") + 
  theme_cowplot() +
  labs(title = "Combination Index") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


SS <- HSA_SS %>% ggplot(aes(x = Drug, y = Estimate, group = Drug)) +
  geom_segment(aes(x= Drug, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = pval), size = 5, color = "gray65", shape = 23) +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec") +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Drugs") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Synergy Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


c <- cowplot::plot_grid(CI, SS)

ab <- plot_grid(NULL, b, rel_widths = c(1, 1.3), labels = "auto")
Fig3 <- plot_grid(ab, c, ncol = 1, labels = c("", "c"), rel_heights = c(1,1.2))

#ggsave("figures/Figure3.pdf", Fig3, height = 1.25*height, width = 2*width)

# Supplementary Figure 4 ----

## Supplementary Fig. 4a ----

MAS98.06_TV0 <- MAS98.06 |> filter(Time == 0)

trt_col <- c("#3c3c3b", "#00a49c", "#fc8d62","#d50c52","#8c510a","#601580","#e78ac3", "#01665e","#386cb0")

a <- MAS98.06_TV0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5) +
  labs(title = "MAS98.06 ER-dependent PDX Initial Tumor Volume") +
  ylab("Initial Tumor Volume (mm3)")
a

## Supplementary Fig. 4b ----
ObsvsPred(MAS98.06_Her_Fulv_lmm, 5, 5)
b <- plot_SynergyLMM(MAS98.06_Her_Fulv_lmm, plot_type = "ObsvsPred", nrow = 5, ncol = 5)
b

## Supplementary Fig. 4c ----
ObsvsPred(MAS98.06_EGF_Fulv_lmm, 4, 5)
c <- plot_SynergyLMM(MAS98.06_EGF_Fulv_lmm, plot_type = "ObsvsPred", nrow = 4, ncol = 5)
c

## Supplementary Fig. 4d ----
ObsvsPred(MAS98.06_Fulv_Trast_lmm, 4, 5)
d <- plot_SynergyLMM(MAS98.06_Fulv_Trast_lmm, plot_type = "ObsvsPred", nrow = 4, ncol = 5)
d

## Supplementary Fig. 4e ----
ObsvsPred(MAS98.06_Her_Trast_lmm, 4, 6)
e <- plot_SynergyLMM(MAS98.06_Her_Trast_lmm, plot_type = "ObsvsPred", nrow = 4, ncol = 5)
e

bottom <- plot_grid(b,c,d,e, ncol = 2, labels = list("b","c","d","e"))

SuppFig4 <- plot_grid(a, bottom, ncol = 1, rel_heights = c(1,2.5), labels = c("a", ""))
SuppFig4

#ggsave("figures/SuppFig4.pdf", SuppFig4, height = 2*height, width = 2.5*width)


# Supplementary Figure 5 ----

## Supplementary Fig. 5a ----

MAS98.06_Her_Fulv_lmm_ranef <- ranefDiagnostics(MAS98.06_Her_Fulv_lmm)
MAS98.06_Her_Fulv_lmm_resid <- residDiagnostics(MAS98.06_Her_Fulv_lmm)

a <- cowplot::plot_grid(plotlist = c(MAS98.06_Her_Fulv_lmm_ranef$Plots[1:4],MAS98.06_Her_Fulv_lmm_resid$Plots[1:5]))
a

## Supplementary Fig. 5b ----

MAS98.06_EGF_Fulv_lmm_ranef <- ranefDiagnostics(MAS98.06_EGF_Fulv_lmm)
MAS98.06_EGF_Fulv_lmm_resid <- residDiagnostics(MAS98.06_EGF_Fulv_lmm)

b <- cowplot::plot_grid(plotlist = c(MAS98.06_EGF_Fulv_lmm_ranef$Plots[1:4],MAS98.06_EGF_Fulv_lmm_resid$Plots[1:5]))
b

## Supplementary Fig. 5c ----

MAS98.06_Fulv_Trast_lmm_ranef <- ranefDiagnostics(MAS98.06_Fulv_Trast_lmm)
MAS98.06_Fulv_Trast_lmm_resid <- residDiagnostics(MAS98.06_Fulv_Trast_lmm)

c <- cowplot::plot_grid(plotlist = c(MAS98.06_Fulv_Trast_lmm_ranef$Plots[1:4],MAS98.06_Fulv_Trast_lmm_resid$Plots[1:5]))
c

## Supplementary Fig. 5d ----

MAS98.06_Her_Trast_lmm_ranef <- ranefDiagnostics(MAS98.06_Her_Trast_lmm)
MAS98.06_Her_Trast_lmm_resid <- residDiagnostics(MAS98.06_Her_Trast_lmm)

d <- cowplot::plot_grid(plotlist = c(MAS98.06_Her_Trast_lmm_ranef$Plots[1:4],MAS98.06_Her_Trast_lmm_resid$Plots[1:5]))
d

SuppFig5 <- cowplot::plot_grid(a,b,c,d, labels = "auto")
SuppFig5

#ggsave("figures/SuppFig5.pdf", SuppFig5, height = 2.25*height, width = 3.25*width)

# Figure 4 ----

## Fig. 4a ----

MAS98.06_Her_Fulv_Trast <- read.xlsx("data/Fig4/MAS98.06_Her_Fulv_Trast.xlsx")
MAS98.06_Her_Fulv_Trast$Sample <- paste(MAS98.06_Her_Fulv_Trast$Treatment, " (",MAS98.06_Her_Fulv_Trast$Sample,")", sep = "")
MAS98.06_Her_Fulv_Trast <- MAS98.06_Her_Fulv_Trast[order(MAS98.06_Her_Fulv_Trast$Sample),]

# Fit Model
MAS98.06_Her_Fulv_Trast_lmm <- lmmModel(data = MAS98.06_Her_Fulv_Trast, sample_id = "Sample", time = "Days", treatment = "Treatment", tumor_vol = "Total_vol",
                                        trt_control = "Ctrl", drug_a = "Her", drug_b = "Trast", drug_c = "Fulv", combination = "Her_Fulv_Trast", time_start = 0, time_end = 28,
                                        weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

trt_col <- c("#3c3c3b", "#00a49c", "#fc8d62","#d50c52","#8c510a","#601580","#e78ac3", "#01665e","#386cb0")

a <- plot_lmmModel(MAS98.06_Her_Fulv_Trast_lmm,   trt_control = "Ctrl", drug_a = "Her", drug_b = "Trast", drug_c = "Fulv", combination = "Her_Fulv_Trast") +
  theme(legend.position = "top") +
  labs(title = "MAS98.06 Heregulin - Fulvestrant - Trastuzumab") + scale_x_continuous(breaks = c(0,7,14,21,28)) +
  ylab("log (Relative Tumor Volume)") + xlab("Days after Treatment Initiation") +
  scale_color_manual(values = c(trt_col[c(1,2,4,5)], "red3")) +
  theme(strip.text = element_text(size = 18))
a

## Fig. 4b ----

(bliss <- lmmSynergy(MAS98.06_Her_Fulv_Trast_lmm, robust = T, type = "CR2", min_time = 8, method = "Bliss"))

(hsa <- lmmSynergy(MAS98.06_Her_Fulv_Trast_lmm, robust = T, type = "CR2", min_time = 8, method = "HSA"))

set.seed(123)
(ra <- lmmSynergy(MAS98.06_Her_Fulv_Trast_lmm, robust = T, type = "CR2", min_time = 8, method = "RA", ra_nsim = 1000))

syn <- rbind(bliss$Synergy, hsa$Synergy, ra$Synergy)

syn <- syn %>% filter(Metric == "SS" & Time %in% c(8,12, 16, 20, 24, 28))

b <- syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) +
  geom_point(aes(fill = pval), size = 5, color = "gray65", shape = 23) +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec") +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after Treatment Initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "MAS98.06 Heregulin - Fulvestrant - Trastuzumab", subtitle = "Synergy Score") + facet_wrap(~Model) +
  scale_x_continuous(breaks = c(8,12, 16, 20, 24, 28)) +
  theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold", size = 18))
b

## Fig. 4c ----

triple <- read.csv("data/Fig4/U87MG_Triple_long.csv")

triple$Day <- sqrt(triple$Day)

# Fit model
triple_lmm <- lmmModel(data = triple, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                       trt_control = "None", drug_a = "AZD2014", drug_b = "Tagrisso", drug_c = "DOC",combination = "Triple",
                       weights = nlme::varComb(nlme::varIdent(form = ~1|SampleID), nlme::varIdent(form = ~1|Time)),
                       control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, opt = "optim"))

c <- plot_lmmModel(triple_lmm, trt_control = "None", drug_a = "AZD2014", drug_b = "Tagrisso", drug_c = "DOC",combination = "Triple") +
  theme(legend.position = "top") +
  labs(title = "U87-MG FM Osimertinib-AZD2014-Docetaxel") + xlab("Days after Treatment Initiation") +
  scale_x_continuous(breaks = unique(triple_lmm$dt1$Time), labels = unique(triple_lmm$dt1$Time)^2) + ylab("log (Relative Luminiscence Units)") +
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1))
c

## Fig. 4d ----

(bliss <- lmmSynergy(triple_lmm, robust = T, type = "CR2", min_time = 3, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

(hsa <- lmmSynergy(triple_lmm, robust = T, type = "CR2", min_time = 3, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2

set.seed(123)
(ra <- lmmSynergy(triple_lmm, robust = T, type = "CR2", min_time = 3, method = "RA", ra_nsim = 1000))
ra$Synergy$Time <- ra$Synergy$Time^2

syn <- rbind(bliss$Synergy, hsa$Synergy, ra$Synergy)

syn <- syn %>% filter(Metric == "SS")

d <- syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) +
  geom_point(aes(fill = pval), size = 5, color = "gray65", shape = 23) +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec") +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after Treatment Initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "U87-MG FM Osimertinib - AZD2014 - Docetaxel", subtitle = "Synergy Score") + facet_wrap(~Model) +
  scale_x_continuous(breaks = unique(syn$Time)) +
  theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold", size = 18))
d

ab <- plot_grid(a, b, nrow = 1, rel_widths = c(1.25), labels = "auto")
cd <- plot_grid(c, d, nrow = 1, rel_widths = c(1.25), labels = c("c", "d"))

Fig4 <- plot_grid(ab, cd, ncol = 1)
Fig4
#ggsave("figures/Figure4.pdf", Fig4, height = 1.25*height, width = 2.5*width)

# Supplementary Figure 6 ----

## Supplementary Fig. 6a ----

MAS98.06_Her_Fulv_Trast_lmm_ranef <- ranefDiagnostics(MAS98.06_Her_Fulv_Trast_lmm)
MAS98.06_Her_Fulv_Trast_lmm_resid <- residDiagnostics(MAS98.06_Her_Fulv_Trast_lmm)

a <- cowplot::plot_grid(plotlist = c(MAS98.06_Her_Fulv_Trast_lmm_ranef$Plots[1:4],MAS98.06_Her_Fulv_Trast_lmm_resid$Plots[1:5]))
a

## Supplementary Fig. 6b ----

ObsvsPred(MAS98.06_Her_Fulv_Trast_lmm, 5, 5)
b <- plot_SynergyLMM(MAS98.06_Her_Fulv_Trast_lmm, plot_type = "ObsvsPred", 5, 5)
b


## Supplementary Fig. 6c ----


MAS98.06_Her_Fulv_Trast_0 <- MAS98.06_Her_Fulv_Trast_lmm$dt1 |> filter(Time == 0)

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#fc8d62", "#601580")

c <- MAS98.06_Her_Fulv_Trast_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5, label.y = 650) +
  labs(title = "MAS98.06 Heregulin - Fulvestrant - Trastuzumab") +
  ylab("Initial Tumor Volume (mm3)")

bc <- plot_grid(b,c, nrow = 1, rel_widths = c(1.75,1), labels = c("b", "c"))

SuppFig6 <- cowplot::plot_grid(a,bc, ncol = 1, rel_heights = c(2,1), labels = c("a",""))
SuppFig6

#ggsave("figures/SuppFig6.pdf", SuppFig6, height = 2.5*height, width = 2.5*width)

# Supplementary Figure 7 ----

## Supplementary Fig. 7a ----

triple_lmm_ranef <- ranefDiagnostics(triple_lmm)
triple_lmm_resid <- residDiagnostics(triple_lmm)

a <- cowplot::plot_grid(plotlist = c(triple_lmm_ranef$Plots[1:4],triple_lmm_resid$Plots[1:5]))
a

## Supplementary Fig. 7b ----

ObsvsPred(triple_lmm, 6, 6)
b <- plot_SynergyLMM(triple_lmm, plot_type = "ObsvsPred", 6, 6)
b


## Supplementary Fig. 7c ----

triple_0 <- triple_lmm$dt1 |> filter(Time == 0)

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#fc8d62", "#601580")

c <- triple_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5, label.y = 6.2e8) +
  labs(title = "U-87-MG AZD2014 - Osimertinib - Docetaxel") +
  ylab("Initial Relative Luminiscence Units") +
  ggbreak::scale_y_break(breaks = c(2e8, 6e8), scales = 0.2)

bc <- plot_grid(b,c, nrow = 1, rel_widths = c(1.75,1), labels = c("b", "c"))

SuppFig7 <- cowplot::plot_grid(a,bc, ncol = 1, rel_heights = c(2,1), labels = c("a",""))
SuppFig7

#ggsave("figures/SuppFig7.pdf", SuppFig7, height = 2.5*height, width = 2.5*width)

# Supplementary Fig. 8 ----

## Supplementary Fig. 8b ----
(hsa <- lmmSynergy(MAS98.06_Her_Fulv_lmm, robust = T, type = "CR2", min_time = 7, method = "HSA", padj = "none"))
hsa$Synergy

a <- plot_lmmSynergy(hsa)$SS +
  labs(title = paste("HSA", "Synergy Score (SynergyLMM)", sep = " "), subtitle = "Fulvestrant vs Heregulin + Fulvestrant") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", na.value = "white", limits = c(0,1))
a

df <- MAS98.06_Her_Fulv_lmm$dt1

# Filter for the two treatments of interest
filtered_data <- df %>%
  filter(Treatment %in% c("Her_Fulv", "Fulv") & Time > 0)

# Group by Time and perform the t-test
ttest_results <- filtered_data %>%
  group_by(Time) %>%
  filter(n_distinct(Treatment) == 2 & n_distinct(SampleID) > 4) %>%  # Ensure both groups are present
  summarise(broom::tidy(t.test(logRTV ~ Treatment)), .groups = "drop")

# View the results
print(ttest_results)

ttest_results$padj <- p.adjust(ttest_results$p.value, method = "none")
print(ttest_results)

b <- ttest_results %>% ggplot(aes(x = .data$Time, y = .data$estimate)) +
  geom_segment(aes(x= .data$Time, y = conf.low, yend = conf.high), color = "gray60", lwd = 1,
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) + cowplot::theme_cowplot() +
  geom_point(aes(fill  = .data$padj), size = 5, shape = 23, color = "gray60") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", na.value = "white", limits = c(0,1)) +
  ylab("T-test Estimate") + xlab("Time since start of treatment") +
  scale_x_continuous(breaks = unique(ttest_results$Time)) + 
  geom_hline(yintercept = 0, lty = "dashed") + #facet_wrap(~Metric) + theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold"))
  labs(title = paste("T-test ", sep = " "), subtitle = "Fulvestrant vs Heregulin + Fulvestrant") + 
  annotate(geom = "text", x = (min(ttest_results$Time)-(ttest_results$Time[2]-ttest_results$Time[1])), 
           y = 0.33, angle = 90, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = (min(ttest_results$Time)-(ttest_results$Time[2]-ttest_results$Time[1])), 
           y = -0.33, angle = 90, hjust = 1, label = "Antagonism", fontface = "bold", color = "#c21d2f")
b

SuppFig8 <- cowplot::plot_grid(a, b, labels = "auto")
SuppFig8

#ggsave("figures/SuppFig8.pdf", SuppFig8, height = 0.66 * height, width = 2*width)

# Figure 5 ----

## Fig. 5a ----

CR1197 <- read.csv("data/Fig5/CR1197_sel_tv.csv")
CR1197$Mice.ID <- paste(CR1197$Treatment, " (",CR1197$Mice.ID,")", sep = "")
CR1197 <- CR1197[order(CR1197$Mice.ID),]

CR1197$DaysPostT0 <- sqrt(CR1197$DaysPostT0)

# Fit model
CR1197_lmm <- lmmModel(data = CR1197, sample_id = "Mice.ID", time = "DaysPostT0", treatment = "Treatment", tumor_vol = "Tumor.Volume",
                       trt_control = "Vehicle", drug_a = "Cetuximab", drug_b = "Palbociclib", combination = "Palbociclib+Cetuximab",
                       weights = nlme::varComb(nlme::varIdent(form = ~1|SampleID), nlme::varPower(form = ~Time)),
                       control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000))


# Generate plot
a <- plot_lmmModel(CR1197_lmm, trt_control = "Vehicle", drug_a = "Cetuximab", drug_b = "Palbociclib", combination = "Palbociclib+Cetuximab") +
  theme(legend.position = "top") +
  labs(title = "CR1197 Cetuximab - Palbociclib") +  ylab("log (Relative Tumor Volume)") +
  scale_x_continuous(breaks = unique(CR1197_lmm$dt1$Time), labels = unique(CR1197_lmm$dt1$Time)^2) + 
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Days after treatment initiation")
a  

## Fig. 5b ----

# Calculate Synergy

## SynergyLMM
(bliss <- lmmSynergy(CR1197_lmm, robust = T, type = "CR2", min_time = 2, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

## invivoSyn
CR1197 <- read.csv("data/Fig5/CR1197_sel_tv.csv")

colnames(CR1197) <- c("Cell.line", "Mouse", "Treatment", "TV", "Day")

# Assign groups

CR1197$Group <- NA

CR1197$Group[CR1197$Treatment == "Vehicle"] <- "Group 1"
CR1197$Group[CR1197$Treatment == "Cetuximab"] <- "Group 2"
CR1197$Group[CR1197$Treatment == "Palbociclib"] <- "Group 3"
CR1197$Group[CR1197$Treatment == "Palbociclib+Cetuximab"] <- "Group 4"

CR1197$Group <- as.factor(CR1197$Group)

# We need an extra column with the initial tumor volume. Let's add it.

TV0 <- CR1197 %>% filter(Day == 0) %>% select(Mouse, TV)

CR1197 <- left_join(CR1197, TV0, by = "Mouse", suffix = c("", "0"))

CR1197 <- na.omit(CR1197)

# Calculate Synergy

AUC_lst <- get_mAUCr(CR1197, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 10000)

# Bliss Synergy with Time

CR1197_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(CR1197$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  CR1197_invivoSyn_bliss <- rbind(CR1197_invivoSyn_bliss, tSyn)
  }

#write.csv(CR1197_invivoSyn_bliss, "data/Fig5/CR1197_invivoSyn_bliss_byT.csv")

# The simulations can take a while to run. Load the results running the following line of code:
CR1197_invivoSyn_bliss <- read.csv("data/Fig5/CR1197_invivoSyn_bliss_byT.csv")

CR1197_bliss <- bliss$Synergy %>% filter(Metric == "SS")

CR1197_bliss <- CR1197_bliss %>% select(Estimate, lwr, upr, pval, Time)

CR1197_bliss$Model <- "SynergyLMM"

CR1197_invivoSyn_bliss <- CR1197_invivoSyn_bliss %>% filter(Metric == "Synergy_score")

CR1197_invivoSyn_bliss <- CR1197_invivoSyn_bliss %>% select(Value, lb, ub, p.val,Time)

CR1197_invivoSyn_bliss <- CR1197_invivoSyn_bliss %>% mutate(Value = Value/100, lb = lb/100, ub = ub/100)

CR1197_invivoSyn_bliss$Model <- "invivoSyn"

colnames(CR1197_invivoSyn_bliss) <- colnames(CR1197_bliss)

CR1197_invivoSyn_bliss <- CR1197_invivoSyn_bliss %>% filter(Time > 0)

CR1197_syn <- rbind(CR1197_invivoSyn_bliss, CR1197_bliss)

CR1197_syn$Time <- as.factor(CR1197_syn$Time)

b <- CR1197_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.75)) +
  geom_point(aes(fill = pval, shape = Model), size = 5, color = "gray65", stroke = 1, position = position_dodge(width = 0.75)) +
  #scale_fill_manual(values = c("firebrick", "darkcyan")) +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec") +
  scale_shape_manual(name = "", values = c(22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Bliss Synergy", subtitle = "CR1197 Cetuximab - Palbociclib") +
  annotate(geom = "text", x = 1, 
           y = 4, angle = 0, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 1, 
           y = -1, angle = 0, hjust = 0, label = "Antagonism", fontface = "bold", color = "#c21d2f")
b
  

## Fig. 5c ----

GABA <- read.csv("data/Fig5/GABA_antiPD1_long.csv")

GABA$Day <- sqrt(GABA$Day-7)

# Fit model
GABA_lmm <- lmmModel(data = GABA, sample_id = "MouseID", time = "Day", treatment = "Treatment", tumor_vol = "TumVol",
                     trt_control = "Ctrl", drug_a = "GABA", drug_b = "Anti-PD1", combination = "GABA+anti-PD1",
                     weights = nlme::varComb(nlme::varPower(form = ~Time), nlme::varIdent(form = ~ 1|Treatment)), 
                     control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000)) 

# Generate plots
c <- plot_lmmModel(GABA_lmm, trt_control = "Ctrl", drug_a = "GABA", drug_b = "Anti-PD1", combination = "GABA+anti-PD1") +
  theme(legend.position = "top") +
  labs(title = "4T1 GABA + Anti-PD1") +  ylab("log (Relative Tumor Volume)") +
  scale_x_continuous(breaks = unique(GABA_lmm$dt1$Time), labels = unique(GABA_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18)) +
  xlab("Days after treatment initiation")
c

## Fig. 5d ----

# Synergy calculation

## SynergyLMM

(hsa <- lmmSynergy(GABA_lmm, robust = T, type = "CR2", min_time = 2.5, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2

## invivoSyn

GABA <- read.csv("data/Fig5/GABA_antiPD1_long.csv")

# Assign groups
colnames(GABA) <- c("Mouse", "Treatment", "Day", "TV", "Cell")

GABA$Day <- GABA$Day -7

GABA$Group <- NA

GABA$Group[GABA$Treatment == "Ctrl"] <- "Group 1"
GABA$Group[GABA$Treatment == "GABA"] <- "Group 2"
GABA$Group[GABA$Treatment == "Anti-PD1"] <- "Group 3"
GABA$Group[GABA$Treatment == "GABA+anti-PD1"] <- "Group 4"

GABA$Group <- as.factor(GABA$Group)

# We need an extra column with the initial tumor volume. Let's add it.

TV0 <- GABA %>% filter(Day == 0) %>% select(Mouse, TV)


GABA <- left_join(GABA, TV0, by = "Mouse", suffix = c("", "0"))

GABA <- na.omit(GABA)

# Calculate Synergy

AUC_lst <- get_mAUCr(GABA, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 10000)

# HSA

# Synergy with Time

GABA_invivoSyn_hsa <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(GABA$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "HSA", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  GABA_invivoSyn_hsa <- rbind(GABA_invivoSyn_hsa, tSyn)
  }

#write.csv(GABA_invivoSyn_hsa, file = "data/Fig5/GABA_invivoSyn_hsa_byT.csv")

# The simulations can take a while to run. Load the results running the following line of code:

GABA_invivoSyn_hsa <- read.csv("data/Fig5/GABA_invivoSyn_hsa_byT.csv")

## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Fig5/4T1_GABA_aPD1_CombPDX.Rdata")
GABA_combpdx_hsa <- hsa_effect


# Generate plots

GABA_hsa <- hsa$Synergy %>% filter(Metric == "SS")

GABA_hsa <- GABA_hsa %>% select(Estimate, lwr, upr, pval, Time)

GABA_hsa$Model <- "SynergyLMM"

GABA_combpdx_hsa <- GABA_combpdx_hsa %>% select(CI, L95, U95, p.value, Days)

GABA_combpdx_hsa$Days <- as.character(GABA_combpdx_hsa$Days %>% str_replace(pattern = "D", replacement = ""))
GABA_combpdx_hsa$Model <- "CombPDX"

colnames(GABA_combpdx_hsa) <- colnames(GABA_hsa)

GABA_invivoSyn_hsa <- GABA_invivoSyn_hsa %>% filter(Metric == "Synergy_score")

GABA_invivoSyn_hsa <- GABA_invivoSyn_hsa %>% select(Value, lb, ub, p.val, Time)

GABA_invivoSyn_hsa$Model <- "invivoSyn"

colnames(GABA_invivoSyn_hsa) <- colnames(GABA_hsa)

GABA_invivoSyn_hsa <- GABA_invivoSyn_hsa %>% filter(Time > 0)
GABA_invivoSyn_hsa[,1:3] <- GABA_invivoSyn_hsa[,1:3]/100

GABA_syn <- rbind(GABA_combpdx_hsa[-1,], GABA_invivoSyn_hsa, GABA_hsa)

GABA_syn$Time <- factor(GABA_syn$Time, levels = c("5", "7", "10", "12", "Global"))

d <- GABA_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "HSA Synergy", subtitle = "4T1 GABA - Anti PD-1") +
  annotate(geom = "text", x = 0.5, 
           y = 0.33, angle = 0, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -2, angle = 0, hjust = 0, label = "Antagonism", fontface = "bold", color = "#c21d2f") +
  geom_text(aes(x = Time, y = 0.1+upr ,colour = Model, label = paste0("p=",round(pval, 3))),
            position = position_dodge(0.5), angle = 45, hjust = 0) +
  scale_color_manual(values = c("slateblue1", "firebrick", "darkcyan"))
d  

## Fig. 5e ----

SW837 <- read.xlsx("data/Fig5/SW837_long.xlsx")

# Fit model
SW837_lmm <- lmmModel(data = SW837, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                      trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                      weights = nlme::varIdent(form = ~1|Treatment))

# Generate plot
e <- plot_lmmModel(SW837_lmm, trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib") +
  theme(legend.position = "top") +
  labs(title = "SW837 Rabusertib - Irinotecan") + ylab("log (Relative Tumor Volume)") +
  theme(strip.text = element_text(size = 18)) +
  xlab("Days after treatment initiation")
e

## Fig. 5f ----

# Calculate Synergy

## SynergyLMM

(bliss <- lmmSynergy(SW837_lmm, robust = T, type = "CR2", min_time = 14, method = "Bliss"))

## invivoSyn

SW837 <- read.xlsx("data/Fig5/SW837_long.xlsx")

SW837$Day <- as.numeric(SW837$Day)
SW837$Day <- SW837$Day + 1

# Assign groups

SW837$Group <- NA

SW837$Group[SW837$Treatment == "Ctrl"] <- "Group 1"
SW837$Group[SW837$Treatment == "Rabusertib"] <- "Group 2"
SW837$Group[SW837$Treatment == "Irinotecan"] <- "Group 3"
SW837$Group[SW837$Treatment == "Irinotecan_Rabusertib"] <- "Group 4"

SW837$Group <- as.factor(SW837$Group)

# We need an extra column with the initial tumor volume. Let's add it.

TV0 <- SW837 %>% filter(Day == 0) %>% select(Mouse, TV)

SW837 <- left_join(SW837, TV0, by = "Mouse", suffix = c("", "0"))

SW837 <- na.omit(SW837)

# Calculate Synergy

AUC_lst <- get_mAUCr(SW837, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 10000)

# Bliss Synergy with Time

SW837_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(SW837$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SW837_invivoSyn_bliss <- rbind(SW837_invivoSyn_bliss, tSyn)
  }

#write.csv(SW837_invivoSyn_bliss, file = "data/Fig5/SW837_invivoSyn_bliss_byT.csv")

# The simulations can take a while to run. Load the results running the following line of code:
SW837_invivoSyn_bliss <- read.csv("data/Fig5/SW837_invivoSyn_bliss_byT.csv")

## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:

load("data/Fig5/SW837_CombPDX.Rdata")
SW837_combpdx_bi <- bi_effect

# Generate plots

SW837_bliss <- bliss$Synergy %>% filter(Metric == "SS")

SW837_bliss <- SW837_bliss %>% select(Estimate, lwr, upr, pval, Time)

SW837_bliss$Model <- "SynergyLMM"

SW837_combpdx_bi <- SW837_combpdx_bi %>% select(CI, L95, U95, p.value, Days)

SW837_combpdx_bi$Days <- as.character(SW837_combpdx_bi$Days %>% str_replace(pattern = "D", replacement = ""))
SW837_combpdx_bi$Model <- "CombPDX"

colnames(SW837_combpdx_bi) <- colnames(SW837_bliss)

SW837_invivoSyn_bliss <- SW837_invivoSyn_bliss %>% filter(Metric == "Synergy_score")

SW837_invivoSyn_bliss <- SW837_invivoSyn_bliss %>% select(Value, lb, ub, p.val, Time)

SW837_invivoSyn_bliss$Model <- "invivoSyn"

colnames(SW837_invivoSyn_bliss) <- colnames(SW837_bliss)

SW837_invivoSyn_bliss <- SW837_invivoSyn_bliss %>% filter(Time > 0)
SW837_invivoSyn_bliss[,1:3] <- SW837_invivoSyn_bliss[,1:3]/100

SW837_syn <- rbind(SW837_combpdx_bi[-1,], SW837_invivoSyn_bliss, SW837_bliss)

SW837_syn$Time <- factor(SW837_syn$Time, levels = c("7", "14", "21", "25", "Global"))

f <- SW837_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Bliss Synergy", subtitle = "SW837 Rabusertib - Irinotecan") +
  annotate(geom = "text", x = 0.5, 
           y = 0.75, angle = 0, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -0.5, angle = 0, hjust = 0, label = "Antagonism", fontface = "bold", color = "#c21d2f") +
  geom_text(aes(x = Time, y = 0.1+upr ,colour = Model, label = paste0("p=",round(pval, 3))),
            position = position_dodge(0.5), angle = 45, hjust = 0) +
  scale_color_manual(values = c("slateblue1", "firebrick", "darkcyan"))
f

ab <- plot_grid(a,b, nrow = 1, labels = "auto", rel_widths = c(1,1.4))
cd <- plot_grid(c,d, nrow = 1, labels = c("c","d"), rel_widths = c(1,1.4))
ef <- plot_grid(e,f, nrow = 1, labels = c("e","f"), rel_widths = c(1,1.4))

Fig5 <- plot_grid(ab, cd, ef, ncol = 1)
Fig5

#ggsave("figures/Figure5.pdf", Fig5, width = 2*width, height = 2*height)


# Supplementary Figure 9 ----

## Supplementary Fig. 9a ----

LS1034 <- read.xlsx("data/Supp_Fig9/LS_1034_long.xlsx")

LS1034_lmm <- lmmModel(data = LS1034, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                       trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                       weights = nlme::varIdent(form = ~1|SampleID), 
                       control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000)) 

a <- plot_lmmModel(LS1034_lmm, trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib") +
  theme(legend.position = "top") +
  labs(title = "LS-1034 Rabusertib - Irinotecan") + ylab("log (Relative Tumor Volume)") +
  theme(strip.text = element_text(size = 18)) +
  xlab("Days after treatment initiation")
a

## Supplementary Fig. 9b ----

# Synergy Calculation

## SynergyLMM

(bliss <- lmmSynergy(LS1034_lmm, robust = T, type = "CR2", min_time = 10, method = "Bliss"))

## invivoSyn

LS1034 <- read.xlsx("data/Supp_Fig9/LS_1034_long.xlsx")

LS1034$Day <- as.numeric(LS1034$Day)
LS1034$Day <- LS1034$Day + 1

# Assign groups

LS1034$Group <- NA

LS1034$Group[LS1034$Treatment == "Ctrl"] <- "Group 1"
LS1034$Group[LS1034$Treatment == "Rabusertib"] <- "Group 2"
LS1034$Group[LS1034$Treatment == "Irinotecan"] <- "Group 3"
LS1034$Group[LS1034$Treatment == "Irinotecan_Rabusertib"] <- "Group 4"

LS1034$Group <- as.factor(LS1034$Group)

# We need an extra column with the initial tumor volume. Let's add it.

TV0 <- LS1034 %>% filter(Day == 0) %>% select(Mouse, TV)


LS1034 <- left_join(LS1034, TV0, by = "Mouse", suffix = c("", "0"))

LS1034 <- na.omit(LS1034)

# Calculate Synergy

AUC_lst <- get_mAUCr(LS1034, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 10000)


# Bliss Synergy with Time

LS1034_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(LS1034$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  LS1034_invivoSyn_bliss <- rbind(LS1034_invivoSyn_bliss, tSyn)
}

#write.csv(LS1034_invivoSyn_bliss, file = "data/Supp_Fig9/LS1034_invivoSyn_bliss_byT.csv")

# The simulations can take a while to run. Load the results running the following line of code:

LS1034_invivoSyn_bliss <- read.csv("data/Supp_Fig9/LS1034_invivoSyn_bliss_byT.csv")


## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:

load("data/Supp_Fig9/LS1034_CombPDX.Rdata")
LS1034_combpdx_bi <- bi_effect

# Generate plots

LS1034_bliss <- bliss$Synergy %>% filter(Metric == "SS")

LS1034_bliss <- LS1034_bliss %>% select(Estimate, lwr, upr, pval, Time)

LS1034_bliss$Model <- "SynergyLMM"

LS1034_combpdx_bi <- LS1034_combpdx_bi %>% select(CI, L95, U95, p.value, Days)

LS1034_combpdx_bi$Days <- as.character(LS1034_combpdx_bi$Days %>% str_replace(pattern = "D", replacement = ""))
LS1034_combpdx_bi$Model <- "CombPDX"

colnames(LS1034_combpdx_bi) <- colnames(LS1034_bliss)

LS1034_invivoSyn_bliss <- LS1034_invivoSyn_bliss %>% filter(Metric == "Synergy_score")

LS1034_invivoSyn_bliss <- LS1034_invivoSyn_bliss %>% select(Value, lb, ub, p.val, Time)

LS1034_invivoSyn_bliss$Model <- "invivoSyn"

colnames(LS1034_invivoSyn_bliss) <- colnames(LS1034_bliss)

LS1034_invivoSyn_bliss <- LS1034_invivoSyn_bliss %>% filter(Time > 0)
LS1034_invivoSyn_bliss[,1:3] <- LS1034_invivoSyn_bliss[,1:3]/100

LS1034_syn <- rbind(LS1034_combpdx_bi[-1,], LS1034_invivoSyn_bliss, LS1034_bliss)

LS1034_syn$Time <- factor(LS1034_syn$Time, levels = c("3", "10", "17", "21", "28", "33", "Global"))


b <- LS1034_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Bliss Synergy", subtitle = "LS1034 Rabusertib - Irinotecan")  +
  annotate(geom = "text", x = 0.5, 
           y = 1, angle = 0, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -01, angle = 0, hjust = 0, label = "Antagonism", fontface = "bold", color = "#c21d2f") +
  geom_text(aes(x = Time, y = 0.1+upr ,colour = Model, label = paste0("p=",round(pval, 3))),
            position = position_dodge(0.5), angle = 45, hjust = 0) +
  scale_color_manual(values = c("slateblue1", "firebrick", "darkcyan"))
b

## Supplementary Fig. 9c ----

SNU81 <- read.xlsx("data/Supp_Fig9/SNU81_long.xlsx")

# Fit model
SNU81_lmm <- lmmModel(data = SNU81, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                      trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                      weights = nlme::varIdent(form = ~1|Time),
                      control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000, opt = "optim"))

# Generate plot
c <- plot_lmmModel(SNU81_lmm, trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib") +
  theme(legend.position = "top") +
  labs(title = "SNU-81 Rabusertib - Irinotecan") + ylab("log (Relative Tumor Volume)") +
  theme(strip.text = element_text(size = 18)) +
  xlab("Days after treatment initiation")
c

## Supplementary Fig. 9d ----

# Synergy Calculation

## SynergyLMM

(bliss <- lmmSynergy(SNU81_lmm, robust = T, type = "CR2", min_time = 0, method = "Bliss"))

## invivoSyn

SNU81 <- read.xlsx("data/Supp_Fig9/SNU81_long.xlsx")

SNU81$Day <- as.numeric(SNU81$Day)
SNU81$Day <- SNU81$Day + 1

# Assign groups

SNU81$Group <- NA

SNU81$Group[SNU81$Treatment == "Ctrl"] <- "Group 1"
SNU81$Group[SNU81$Treatment == "Rabusertib"] <- "Group 2"
SNU81$Group[SNU81$Treatment == "Irinotecan"] <- "Group 3"
SNU81$Group[SNU81$Treatment == "Irinotecan_Rabusertib"] <- "Group 4"

SNU81$Group <- as.factor(SNU81$Group)

# We need an extra column with the initial tumor volume. Let's add it.

TV0 <- SNU81 %>% filter(Day == 0) %>% select(Mouse, TV)

SNU81 <- left_join(SNU81, TV0, by = "Mouse", suffix = c("", "0"))

SNU81 <- na.omit(SNU81)

# Calculate Synergy

AUC_lst <- get_mAUCr(SNU81, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 10000)

# Bliss Synergy with Time

SNU81_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(SNU81$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SNU81_invivoSyn_bliss <- rbind(SNU81_invivoSyn_bliss, tSyn)
}

#write.csv(SNU81_invivoSyn_bliss, file = "data/Supp_Fig9/SNU81_invivoSyn_bliss_byT.csv")

# The simulations can take a while to run. Load the results running the following line of code:

SNU81_invivoSyn_bliss <- read.csv("data/Supp_Fig9/SNU81_invivoSyn_bliss_byT.csv")

## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:

load("data/Supp_Fig9/SNU81_CombPDX.Rdata")
SNU81_combpdx_bi <- bi_effect

# Generate plots

SNU81_bliss <- bliss$Synergy %>% filter(Metric == "SS")

SNU81_bliss <- SNU81_bliss %>% select(Estimate, lwr, upr, pval, Time)

SNU81_bliss$Model <- "SynergyLMM"

SNU81_combpdx_bi <- SNU81_combpdx_bi %>% select(CI, L95, U95, p.value, Days)

SNU81_combpdx_bi$Days <- as.character(SNU81_combpdx_bi$Days %>% str_replace(pattern = "D", replacement = ""))
SNU81_combpdx_bi$Model <- "CombPDX"

colnames(SNU81_combpdx_bi) <- colnames(SNU81_bliss)


SNU81_invivoSyn_bliss <- SNU81_invivoSyn_bliss %>% filter(Metric == "Synergy_score")

SNU81_invivoSyn_bliss <- SNU81_invivoSyn_bliss %>% select(Value, lb, ub, p.val, Time)

SNU81_invivoSyn_bliss$Model <- "invivoSyn"

colnames(SNU81_invivoSyn_bliss) <- colnames(SNU81_bliss)

SNU81_invivoSyn_bliss <- SNU81_invivoSyn_bliss %>% filter(Time > 0)
SNU81_invivoSyn_bliss[,1:3] <- SNU81_invivoSyn_bliss[,1:3]/100

SNU81_syn <- rbind(SNU81_combpdx_bi[-1,], SNU81_invivoSyn_bliss, SNU81_bliss)

SNU81_syn$Time <- factor(SNU81_syn$Time, levels = c("8", "15", "19", "30", "36", "Global"))

d <-  SNU81_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Bliss Synergy", subtitle = "SNU-81 Rabusertib - Irinotecan") +
  annotate(geom = "text", x = 0.5, 
           y = 1, angle = 90, hjust = 0.5, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -1, angle = 90, hjust = 0.5, label = "Antagonism", fontface = "bold", color = "#c21d2f") +
  geom_text(aes(x = Time, y = 0.1+upr ,colour = Model, label = paste0("p=",round(pval, 3))),
            position = position_dodge(0.5), angle = 45, hjust = 0) +
  scale_color_manual(values = c("slateblue1", "firebrick", "darkcyan"))

ab <- plot_grid(a,b, nrow = 1, labels = "auto", rel_widths = c(1,1.4))
cd <- plot_grid(c,d, nrow = 1, labels = c("c","d"), rel_widths = c(1,1.4))

SuppFig9 <- plot_grid(ab, cd, ncol = 1)
SuppFig9
#ggsave("figures/SuppFig9.pdf", SuppFig9, width = 2*width, height = 1.33*height)

# Supplementary Figure 10 ----

## Supplementary Fig. 10a ----

CR1197_0 <- CR1197_lmm$dt1 |> filter(Time == 0)

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")


a <- CR1197_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5) +
  labs(title = "CR1197 Cetuximab - Palbociclib") +
  ylab("Initial Tumor Volume (mm3)")

## Supplementary Fig. 10b ----

CR1197_lmm_ranef <- ranefDiagnostics(CR1197_lmm)
CR1197_lmm_resid <- residDiagnostics(CR1197_lmm)

ab <- cowplot::plot_grid(a, plotlist = c(CR1197_lmm_ranef$Plots[1:4],CR1197_lmm_resid$Plots[1:5]), nrow = 2,
                         labels = c("a", "b"))
ab

## Supplementary Fig. 10c ----

GABA_0 <- GABA_lmm$dt1 |> filter(Time == 0)

range(df$TV)

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

stat.test <- ggpubr::compare_means(TV ~ Treatment, data = GABA_0, method = "wilcox.test", p.adjust.method = "BH")

stat.test$y.position <- c(0, 0,0, 33, 35, 0)

c <- GABA_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5, label.y = 40) +
  ggpubr::stat_pvalue_manual(data = stat.test, hide.ns = T, label = "p={p.adj}") +
  labs(title = "4T1 GABA - Anti-PD1") +
  ylab("Initial Tumor Volume (mm3)")

## Supplementary Fig. 10d ----

GABA_lmm_ranef <- ranefDiagnostics(GABA_lmm)
GABA_lmm_resid <- residDiagnostics(GABA_lmm)

cd <- cowplot::plot_grid(c, plotlist = c(GABA_lmm_ranef$Plots[1:4],GABA_lmm_resid$Plots[1:5]), nrow = 2,
                         labels = c("c", "d"))
cd

## Supplementary Fig. 10e ----

SW837_0 <- SW837_lmm$dt1 |> filter(Time == 0)

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

e <- SW837_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5) +
  labs(title = "SW837")+
  ylab("Initial Tumor Volume (mm3)")

## Supplementary Fig. 10f ----

SW837_lmm_ranef <- ranefDiagnostics(SW837_lmm)
SW837_lmm_resid <- residDiagnostics(SW837_lmm)

ef <- cowplot::plot_grid(e, plotlist = c(SW837_lmm_ranef$Plots[1:4],SW837_lmm_resid$Plots[1:5]), nrow = 2,
                         labels = c("e", "f"))
ef

## Supplementary Fig. 10g ----

LS1034_0 <- LS1034_lmm$dt1 |> filter(Time == 0)

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

g <- LS1034_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5) +
  labs(title = "LS-1034")+
  ylab("Initial Tumor Volume (mm3)")

## Supplementary Fig. 10h ----

LS1034_lmm_ranef <- ranefDiagnostics(LS1034_lmm)
LS1034_lmm_resid <- residDiagnostics(LS1034_lmm)

gh <- cowplot::plot_grid(g, plotlist = c(LS1034_lmm_ranef$Plots[1:4],LS1034_lmm_resid$Plots[1:5]), nrow = 2,
                         labels = c("g", "h"))
gh

## Supplementary Fig. 10i ----

SNU81_0 <- SNU81_lmm$dt1 |> filter(Time == 0)

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

i <- SNU81_0 |> ggplot(aes(x = Treatment, y = TV)) + geom_boxplot(aes(color = Treatment, fill = Treatment), alpha = 0.33, outliers = F) +
  geom_point(aes(color = Treatment), position = position_jitter(widt = 0.2), size = 3) +
  scale_fill_manual(values = trt_col) + scale_color_manual(values = trt_col) + theme_cowplot() +
  theme(legend.position = "top") + ggpubr::stat_compare_means(method = "kruskal.test", size = 5) +
  labs(title = "SNU-81")+
  ylab("Initial Tumor Volume (mm3)")

## Supplementary Fig. 10j ----

SNU81_lmm_ranef <- ranefDiagnostics(SNU81_lmm)
SNU81_lmm_resid <- residDiagnostics(SNU81_lmm)

ij <- cowplot::plot_grid(i, plotlist = c(SNU81_lmm_ranef$Plots[1:4],SNU81_lmm_resid$Plots[1:5]), nrow = 2,
                         labels = c("i", "j"))
ij

SuppFig10 <- cowplot::plot_grid(ab,cd,ef,gh,ij, nrow = 5)
SuppFig10

#ggsave("figures/SuppFig10.pdf", SuppFig10, height = 4.2*height, width = 4.2*width)

# Supplementary Figure 11 ----

## Supplementary Fig. 11a ----

ObsvsPred(CR1197_lmm, 6, 7)
a <- plot_SynergyLMM(CR1197_lmm, plot_type = "ObsvsPred", 6,7)

## Supplementary Fig. 11b ----

ObsvsPred(GABA_lmm, 6, 7)
b <- plot_SynergyLMM(GABA_lmm, plot_type = "ObsvsPred", 6,7)

## Supplementary Fig. 11c ----

ObsvsPred(SW837_lmm, 5, 5)
c <- plot_SynergyLMM(SW837_lmm, plot_type = "ObsvsPred", 5, 5)

## Supplementary Fig. 11d ----

ObsvsPred(LS1034_lmm, 6, 7)
d <- plot_SynergyLMM(LS1034_lmm, plot_type = "ObsvsPred", 6, 7)

## Supplementary Fig. 11e ----

ObsvsPred(SNU81_lmm, 5, 6)
e <- plot_SynergyLMM(SNU81_lmm, plot_type = "ObsvsPred", 5, 6)

SuppFig11 <- cowplot::plot_grid(a,b,c,d,e, labels = "auto", nrow = 3)
SuppFig11

#ggsave("figures/SuppFig11.pdf", SuppFig11, height = 1.5*height, width = 1.5*width)

# Supplementary Figure 12 ----

## Supplementary Fig. 12a ----

# Create simulated tumor growths

simulatedSyn <- data.frame(Estimate = numeric(0), lwr = numeric(0), upr = numeric(0), pval = numeric(0), Time = numeric(0), Simulation = numeric(0))

for (i in 1:1000) {
  set.seed(i)
  grwth_data_1 <- simulateTumorGrowth(npg = 5, timepoints = seq(0,30,3), 
                                      initial_volume = 200, grwrControl = 0.08,
                                      grwrA = 0.07, grwrB = 0.065, 
                                      grwrComb = 0.03, sd = 0.15)
  
  
  ## Build Models
  
  
  lmm1 <- lmmModel(data = grwth_data_1, sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                   trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination", show_plot = F)
  
  
  ## Synergy
  
  lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T, show_plot = F)
  lmm1_Bliss <- lmm1_Bliss$Synergy
  lmm1_Bliss$Simulation <- i
  
  simulatedSyn <- rbind(simulatedSyn, lmm1_Bliss)
  
}

#write.csv(simulatedSyn, file = "data/Supp_Fig12/SimulatedSyn1000.csv")

# The simulations can take some time to run. The results can be loaded as:

simulatedSyn <- read.csv("data/Supp_Fig12/SimulatedSyn1000.csv", row.names = 1)

# Generate plot of summary results

lwrSS <- simulatedSyn %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(lwr))
uprSS <- simulatedSyn %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(upr))
estimateSS <- simulatedSyn %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(Estimate))
pvalSS <- simulatedSyn %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(pval))

simSS <- cbind(lwrSS, uprSS[,-1], estimateSS[,-1], pvalSS[,-1])

simSS$Metric <- "SS"

lwrCI <- simulatedSyn %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(lwr))
uprCI <- simulatedSyn %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(upr))
estimateCI <- simulatedSyn %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(Estimate))
pvalCI <- simulatedSyn %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(pval))

simCI <- cbind(lwrCI, uprCI[,-1], estimateCI[,-1], pvalCI[,-1])
simCI$Metric <- "CI"

simSimulated <- rbind(simSS, simCI)


a <- simSS %>% ggplot(aes(x = Time, y = `median(Estimate)`)) +
  geom_segment(aes(x= Time, y = `median(lwr)`, yend = `median(upr)`), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = `median(pval)`), shape = 23, size = 5, color = "gray65") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", limits = c(0,0.5)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  scale_x_continuous(breaks = unique(simSS$Time)) +
  labs(title = "Exponential Growth Simulated Data - Synergistic Effect")
a

## Supplementary Fig. 12b ----

simulatedAdd <- data.frame(Estimate = numeric(0), lwr = numeric(0), upr = numeric(0), pval = numeric(0), Time = numeric(0), Simulation = numeric(0))

for (i in 1:1000) {
  print(i)
  set.seed(i)
  grwth_data_1 <- simulateTumorGrowth(npg = 5, timepoints = seq(0,30,3), 
                                      initial_volume = 200, grwrControl = 0.08,
                                      grwrA = 0.07, grwrB = 0.065, 
                                      grwrComb = 0.055, sd = 0.15)
  
  
  ## Build Models
  
  
  lmm1 <- lmmModel(data = grwth_data_1, sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                   trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination", show_plot = F)
  
  
  ## Synergy
  
  lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T, show_plot = F)
  lmm1_Bliss <- lmm1_Bliss$Synergy
  lmm1_Bliss$Simulation <- i
  
  simulatedAdd <- rbind(simulatedAdd, lmm1_Bliss)
  
}

#write.csv(simulatedAdd, file = "data/Supp_Fig12/SimulatedAdd1000.csv")

# The simulations can take some time to run. The results can be loaded as:

simulatedAdd <- read.csv("data/Supp_Fig12/SimulatedAdd1000.csv", row.names = 1)

lwrSS <- simulatedAdd %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(lwr))
uprSS <- simulatedAdd %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(upr))
estimateSS <- simulatedAdd %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(Estimate))
pvalSS <- simulatedAdd %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(pval))

simSS <- cbind(lwrSS, uprSS[,-1], estimateSS[,-1], pvalSS[,-1])
simSS$Metric <- "SS"


lwrCI <- simulatedAdd %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(lwr))
uprCI <- simulatedAdd %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(upr))
estimateCI <- simulatedAdd %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(Estimate))
pvalCI <- simulatedAdd %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(pval))

simCI <- cbind(lwrCI, uprCI[,-1], estimateCI[,-1], pvalCI[,-1])
simCI$Metric <- "CI"

addSimulated <- rbind(simSS, simCI)

b <- simSS %>% ggplot(aes(x = Time, y = `median(Estimate)`)) +
  geom_segment(aes(x= Time, y = `median(lwr)`, yend = `median(upr)`), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = `median(pval)`), shape = 23, size = 5, color = "gray65") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", limits = c(0, 0.5)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  scale_x_continuous(breaks = unique(simSS$Time)) +
  labs(title = "Exponential Growth Simulated Data - Additive Effect")
b

## Supplementary Fig. 12c ----

# Generate data
set.seed(123)
grwth_data_1 <- simulateTumorGrowth(npg = 5, timepoints = seq(0,30,3), 
                                    initial_volume = 200, grwrControl = 0.08,
                                    grwrA = 0.07, grwrB = 0.065, 
                                    grwrComb = 0.03, sd = 0.15)

grwth_data_1$Cell <- "None"

# Fit model
lmm1 <- lmmModel(data = grwth_data_1, sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                 trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination")

# Generate plot
c <- plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Exponential Growth Simulated Data - Synergistic Effect") +
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Days after treatment initiation")
c

## Supplementary Fig. 12d ----

# Calculate Synergy

## SynergyLMM

lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T)

## invivoSyn

colnames(grwth_data_1) <- c("Mouse", "Day", "Treatment", "TV", "Cell")

# Assign groups

grwth_data_1$Group <- NA

unique(grwth_data_1$Treatment)

grwth_data_1$Group[grwth_data_1$Treatment == "Control"] <- "Group 1"
grwth_data_1$Group[grwth_data_1$Treatment == "DrugA"] <- "Group 2"
grwth_data_1$Group[grwth_data_1$Treatment == "DrugB"] <- "Group 3"
grwth_data_1$Group[grwth_data_1$Treatment == "Combination"] <- "Group 4"

grwth_data_1$Group <- as.factor(grwth_data_1$Group)

# We need an extra column with the initial tumor volume. Let's add it.

TV0 <- grwth_data_1 %>% filter(Day == 0) %>% select(Mouse, TV)


grwth_data_1 <- left_join(grwth_data_1, TV0, by = "Mouse", suffix = c("", "0"))

grwth_data_1 <- na.omit(grwth_data_1)

# Calculate Synergy

AUC_lst <- get_mAUCr(grwth_data_1, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 10000)

# Bliss Synergy with Time

SimData1_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(grwth_data_1$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss",ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SimData1_invivoSyn_bliss <- rbind(SimData1_invivoSyn_bliss, tSyn)
}

#write.csv(SimData1_invivoSyn_bliss, file = "data/Supp_Fig12/SimData1_invivoSyn_bliss_byT.csv")

# The simulations can take a while to run. Load the results running the following line of code:

SimData1_invivoSyn_bliss <- read.csv("data/Supp_Fig12/SimData1_invivoSyn_bliss_byT.csv", row.names = 1)


## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Supp_Fig12/SimData1_CombPDX.Rdata")
SimData1_combpdx_bi <- bi_effect

# Generate plots

SimData1_bliss <- lmm1_Bliss$Synergy %>% filter(Metric == "SS")

SimData1_bliss <- SimData1_bliss %>% select(Estimate, lwr, upr, pval, Time)

SimData1_bliss$Model <- "SynergyLMM"

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% filter(Metric == "Synergy_score")

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% select(Value, lb, ub, p.val, Time)

SimData1_invivoSyn_bliss$Model <- "invivoSyn"

colnames(SimData1_invivoSyn_bliss) <- colnames(SimData1_bliss)

SimData1_combpdx_bi <- SimData1_combpdx_bi %>% select(CI, L95, U95, p.value, Days)

SimData1_combpdx_bi$Days <- as.character(SimData1_combpdx_bi$Days %>% str_replace(pattern = "D", replacement = ""))
SimData1_combpdx_bi$Model <- "CombPDX"
colnames(SimData1_combpdx_bi) <- colnames(SimData1_bliss)

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% filter(Time > 0)
SimData1_invivoSyn_bliss[,1:3] <- SimData1_invivoSyn_bliss[,1:3]/100

SimData1_syn <- rbind(SimData1_combpdx_bi[-1,], SimData1_invivoSyn_bliss, SimData1_bliss)


SimData1_syn$Time <- factor(SimData1_syn$Time, levels = c("3", "6", "9", "12","15", "18", "21", "24", "27", "30", "Global"))

d <- SimData1_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Bliss Synergy") +
  annotate(geom = "text", x = 0.5, 
           y = 5, angle = 0, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -1.5, angle = 0, hjust = 0, label = "Antagonism", fontface = "bold", color = "#c21d2f") +
  geom_text(aes(x = Time, y = 0.1+upr ,colour = Model, label = paste0("p=",round(pval, 3))),
            position = position_dodge(0.5), angle = 45, hjust = 0) +
  scale_color_manual(values = c("slateblue1", "firebrick", "darkcyan"))
d

## Supplementary Fig. 12e ----

# Generate data
set.seed(123)
grwth_data_1 <- simulateTumorGrowth(npg = 5, timepoints = seq(0,30,3), 
                                    initial_volume = 200, grwrControl = 0.08,
                                    grwrA = 0.07, grwrB = 0.065, 
                                    grwrComb = 0.03, sd = 0.15)


outlier1 <- data.frame(subject = "outlier", Time = seq(0,30,3),
                       Treatment = "Combination")
outlier1$TumorVolume <- 200*exp(0.08*seq(0,30,3))
set.seed(123)
outlier1$TumorVolume <- outlier1$TumorVolume * rnorm(nrow(outlier1), mean = 1, sd = 0.15)

grwth_data_1 <- rbind(grwth_data_1, outlier1)

grwth_data_1$Cell <- "None"

# Fit model
lmm1 <- lmmModel(data = grwth_data_1, sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                 trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination")

# Generate plot
e <- plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Exponential Growth Simulated Data - Synergistic Effect") +
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Days after treatment initiation")
e

## Supplementary Fig. 12f ----

# Calculate Synergy

## SynergyLMM

lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T)

## invivoSyn

colnames(grwth_data_1) <- c("Mouse", "Day", "Treatment", "TV", "Cell")

# Assign groups

grwth_data_1$Group <- NA

grwth_data_1$Group[grwth_data_1$Treatment == "Control"] <- "Group 1"
grwth_data_1$Group[grwth_data_1$Treatment == "DrugA"] <- "Group 2"
grwth_data_1$Group[grwth_data_1$Treatment == "DrugB"] <- "Group 3"
grwth_data_1$Group[grwth_data_1$Treatment == "Combination"] <- "Group 4"

grwth_data_1$Group <- as.factor(grwth_data_1$Group)

# We need an extra column with the initial tumor volume. Let's add it.

TV0 <- grwth_data_1 %>% filter(Day == 0) %>% select(Mouse, TV)

grwth_data_1 <- left_join(grwth_data_1, TV0, by = "Mouse", suffix = c("", "0"))

grwth_data_1 <- na.omit(grwth_data_1)

# Calculate Synergy

AUC_lst <- get_mAUCr(grwth_data_1, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 10000)

# Bliss Synergy with Time

SimData1_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(grwth_data_1$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss",ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SimData1_invivoSyn_bliss <- rbind(SimData1_invivoSyn_bliss, tSyn)
}

#write.csv(SimData1_invivoSyn_bliss, file = "data/Supp_Fig12/SimData1_invivoSyn_bliss_byT_outlier.csv")

# The simulations can take a while to run. Load the results running the following line of code:

SimData1_invivoSyn_bliss <- read.csv("data/Supp_Fig12/SimData1_invivoSyn_bliss_byT_outlier.csv", row.names = 1)


## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Supp_Fig12/SimData1_CombPDX_outlier.Rdata")
SimData1_combpdx_bi <- bi_effect

# Generate plots

SimData1_bliss <- lmm1_Bliss$Synergy %>% filter(Metric == "SS")

SimData1_bliss <- SimData1_bliss %>% select(Estimate, lwr, upr, pval, Time)

SimData1_bliss$Model <- "SynergyLMM"

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% filter(Metric == "Synergy_score")

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% select(Value, lb, ub, p.val, Time)

SimData1_invivoSyn_bliss$Model <- "invivoSyn"

colnames(SimData1_invivoSyn_bliss) <- colnames(SimData1_bliss)

SimData1_combpdx_bi <- SimData1_combpdx_bi %>% select(CI, L95, U95, p.value, Days)

SimData1_combpdx_bi$Days <- as.character(SimData1_combpdx_bi$Days %>% str_replace(pattern = "D", replacement = ""))
SimData1_combpdx_bi$Model <- "CombPDX"
colnames(SimData1_combpdx_bi) <- colnames(SimData1_bliss)

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% filter(Time > 0)
SimData1_invivoSyn_bliss[,1:3] <- SimData1_invivoSyn_bliss[,1:3]/100

SimData1_syn <- rbind(SimData1_combpdx_bi[-1,], SimData1_invivoSyn_bliss, SimData1_bliss)


SimData1_syn$Time <- factor(SimData1_syn$Time, levels = c("3", "6", "9", "12","15", "18", "21", "24", "27", "30", "Global"))

f <- SimData1_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Bliss Synergy") +
  annotate(geom = "text", x = 0.5, 
           y = 5, angle = 0, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -1.5, angle = 0, hjust = 0, label = "Antagonism", fontface = "bold", color = "#c21d2f") +
  geom_text(aes(x = Time, y = 0.1+upr ,colour = Model, label = paste0("p=",round(pval, 3))),
            position = position_dodge(0.5), angle = 45, hjust = 0) +
  scale_color_manual(values = c("slateblue1", "firebrick", "darkcyan"))
f

## Supplementary Fig. 12g ----

# Generate data
set.seed(123)
grwth_data_1 <- simulateTumorGrowth(npg = 5, timepoints = seq(0,30,3), 
                                    initial_volume = 200, grwrControl = 0.08,
                                    grwrA = 0.07, grwrB = 0.065, 
                                    grwrComb = 0.055, sd = 0.15)


grwth_data_1$Cell <- "None"

# Fit model
lmm1 <- lmmModel(data = grwth_data_1, sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                 trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination")

# Generate plot
g <- plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Exponential Growth Simulated Data - Additive Effect") +
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Days after treatment initiation")
g

## Supplementary Fig. 12h ----

# Calculate Synergy

## SynergyLMM

lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T)

## invivoSyn

colnames(grwth_data_1) <- c("Mouse", "Day", "Treatment", "TV", "Cell")

# Assign groups

grwth_data_1$Group <- NA

grwth_data_1$Group[grwth_data_1$Treatment == "Control"] <- "Group 1"
grwth_data_1$Group[grwth_data_1$Treatment == "DrugA"] <- "Group 2"
grwth_data_1$Group[grwth_data_1$Treatment == "DrugB"] <- "Group 3"
grwth_data_1$Group[grwth_data_1$Treatment == "Combination"] <- "Group 4"

grwth_data_1$Group <- as.factor(grwth_data_1$Group)

# We need an extra column with the initial tumor volume. Let's add it.

TV0 <- grwth_data_1 %>% filter(Day == 0) %>% select(Mouse, TV)


grwth_data_1 <- left_join(grwth_data_1, TV0, by = "Mouse", suffix = c("", "0"))

grwth_data_1 <- na.omit(grwth_data_1)

# Calculate Synergy

AUC_lst <- get_mAUCr(grwth_data_1, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 10000)

# Bliss Synergy with Time

SimData1_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(grwth_data_1$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss",ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SimData1_invivoSyn_bliss <- rbind(SimData1_invivoSyn_bliss, tSyn)
}

#write.csv(SimData1_invivoSyn_bliss, file = "data/Supp_Fig12/SimData1_invivoSyn_bliss_byTAdditive.csv")

# The simulations can take a while to run. Load the results running the following line of code:

SimData1_invivoSyn_bliss <- read.csv("data/Supp_Fig12/SimData1_invivoSyn_bliss_byTAdditive.csv", row.names = 1)


## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Supp_Fig12/SimData1_CombPDX_Additive.Rdata")
SimData1_combpdx_bi <- bi_effect

# Generate plots

SimData1_bliss <- lmm1_Bliss$Synergy %>% filter(Metric == "SS")

SimData1_bliss <- SimData1_bliss %>% select(Estimate, lwr, upr, pval, Time)

SimData1_bliss$Model <- "SynergyLMM"

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% filter(Metric == "Synergy_score")

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% select(Value, lb, ub, p.val, Time)

SimData1_invivoSyn_bliss$Model <- "invivoSyn"

colnames(SimData1_invivoSyn_bliss) <- colnames(SimData1_bliss)

SimData1_combpdx_bi <- SimData1_combpdx_bi %>% select(CI, L95, U95, p.value, Days)

SimData1_combpdx_bi$Days <- as.character(SimData1_combpdx_bi$Days %>% str_replace(pattern = "D", replacement = ""))
SimData1_combpdx_bi$Model <- "CombPDX"
colnames(SimData1_combpdx_bi) <- colnames(SimData1_bliss)

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% filter(Time > 0)
SimData1_invivoSyn_bliss[,1:3] <- SimData1_invivoSyn_bliss[,1:3]/100

SimData1_syn <- rbind(SimData1_combpdx_bi[-1,], SimData1_invivoSyn_bliss, SimData1_bliss)


SimData1_syn$Time <- factor(SimData1_syn$Time, levels = c("3", "6", "9", "12","15", "18", "21", "24", "27", "30", "Global"))

h <- SimData1_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Bliss Synergy") +
  annotate(geom = "text", x = 0.5, 
           y = 3, angle = 0, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -2, angle = 0, hjust = 0, label = "Antagonism", fontface = "bold", color = "#c21d2f") +
  geom_text(aes(x = Time, y = 0.1+upr ,colour = Model, label = paste0("p=",round(pval, 3))),
            position = position_dodge(0.5), angle = 45, hjust = 0) +
  scale_color_manual(values = c("slateblue1", "firebrick", "darkcyan"))
h

## Supplementary Fig. 12i ----

# Generate data
set.seed(123)
grwth_data_1 <- simulateTumorGrowth(npg = 5, timepoints = seq(0,30,3), 
                                    initial_volume = 200, grwrControl = 0.08,
                                    grwrA = 0.07, grwrB = 0.065, 
                                    grwrComb = 0.055, sd = 0.15)


outlier1 <- data.frame(subject = "outlier", Time = seq(0,30,3),
                       Treatment = "Combination")
outlier1$TumorVolume <- 200*exp(-0.01*seq(0,30,3))
set.seed(123)
outlier1$TumorVolume <- outlier1$TumorVolume * rnorm(nrow(outlier1), mean = 1, sd = 0.15)

grwth_data_1 <- rbind(grwth_data_1, outlier1)

grwth_data_1$Cell <- "None"

# Fit model
lmm1 <- lmmModel(data = grwth_data_1, sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                 trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination")

# Generate plot
i <- plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Simulated Data - Additive Effect") +
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Days after treatment initiation")
i

## Supplementary Fig. 12j ----

# Calculate Synergy

## SynergyLMM

lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T)

## invivoSyn

colnames(grwth_data_1) <- c("Mouse", "Day", "Treatment", "TV", "Cell")

# Assign groups

grwth_data_1$Group <- NA

grwth_data_1$Group[grwth_data_1$Treatment == "Control"] <- "Group 1"
grwth_data_1$Group[grwth_data_1$Treatment == "DrugA"] <- "Group 2"
grwth_data_1$Group[grwth_data_1$Treatment == "DrugB"] <- "Group 3"
grwth_data_1$Group[grwth_data_1$Treatment == "Combination"] <- "Group 4"

grwth_data_1$Group <- as.factor(grwth_data_1$Group)

# We need an extra column with the initial tumor volume. Let's add it.

TV0 <- grwth_data_1 %>% filter(Day == 0) %>% select(Mouse, TV)


grwth_data_1 <- left_join(grwth_data_1, TV0, by = "Mouse", suffix = c("", "0"))

grwth_data_1 <- na.omit(grwth_data_1)

# Calculate Synergy

AUC_lst <- get_mAUCr(grwth_data_1, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 10000)

# Bliss Synergy with Time

SimData1_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(grwth_data_1$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss",ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SimData1_invivoSyn_bliss <- rbind(SimData1_invivoSyn_bliss, tSyn)
}

#write.csv(SimData1_invivoSyn_bliss, file = "data/Supp_Fig12/SimData1_invivoSyn_bliss_byTAdditive_outlier.csv")

# The simulations can take a while to run. Load the results running the following line of code:

SimData1_invivoSyn_bliss <- read.csv("data/Supp_Fig12/SimData1_invivoSyn_bliss_byTAdditive_outlier.csv", row.names = 1)

## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Supp_Fig12/SimData1_CombPDX_Additive_outlier.Rdata")
SimData1_combpdx_bi <- bi_effect

# Generate plots

SimData1_bliss <- lmm1_Bliss$Synergy %>% filter(Metric == "SS")

SimData1_bliss <- SimData1_bliss %>% select(Estimate, lwr, upr, pval, Time)

SimData1_bliss$Model <- "SynergyLMM"

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% filter(Metric == "Synergy_score")

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% select(Value, lb, ub, p.val, Time)

SimData1_invivoSyn_bliss$Model <- "invivoSyn"

colnames(SimData1_invivoSyn_bliss) <- colnames(SimData1_bliss)

SimData1_combpdx_bi <- SimData1_combpdx_bi %>% select(CI, L95, U95, p.value, Days)

SimData1_combpdx_bi$Days <- as.character(SimData1_combpdx_bi$Days %>% str_replace(pattern = "D", replacement = ""))
SimData1_combpdx_bi$Model <- "CombPDX"
colnames(SimData1_combpdx_bi) <- colnames(SimData1_bliss)

SimData1_invivoSyn_bliss <- SimData1_invivoSyn_bliss %>% filter(Time > 0)
SimData1_invivoSyn_bliss[,1:3] <- SimData1_invivoSyn_bliss[,1:3]/100

SimData1_syn <- rbind(SimData1_combpdx_bi[-1,], SimData1_invivoSyn_bliss, SimData1_bliss)


SimData1_syn$Time <- factor(SimData1_syn$Time, levels = c("3", "6", "9", "12","15", "18", "21", "24", "27", "30", "Global"))

j <- SimData1_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Bliss Synergy") +
  annotate(geom = "text", x = 0.5, 
           y = 3.5, angle = 0, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -1.75, angle = 0, hjust = 0, label = "Antagonism", fontface = "bold", color = "#c21d2f") +
  geom_text(aes(x = Time, y = 0.1+upr ,colour = Model, label = paste0("p=", round(pval, 3))),
            position = position_dodge(0.5), angle = 45, hjust = 0) +
  scale_color_manual(values = c("slateblue1", "firebrick", "darkcyan"))
j

ab <- plot_grid(a, b, nrow = 1, labels = "auto")
ab
cd <- plot_grid(c,d, nrow = 1, rel_widths = c(1,1.66), labels = c("c", "d"))
cd
ef <- plot_grid(e,f, nrow = 1, rel_widths = c(1,1.66), labels = c("e", "f"))
ef
gh <- plot_grid(g,h, nrow = 1, rel_widths = c(1,1.66), labels = c("g", "h"))
gh
ij <- plot_grid(i,j, nrow = 1, rel_widths = c(1,1.66), labels = c("i", "j"))
ij

SuppFig12 <- plot_grid(ab,cd,ef,gh,ij, ncol = 1)
SuppFig12

#ggsave("figures/SuppFig12.pdf", SuppFig12, height = 2.66*height, width = 2.66*width)

# Supplementary Fig. 13 ----

# Define logistic model
log_fun <- function(V0, K, r, t){
  K/(1+((K-V0)/V0)*exp(-r*t))
}


LogTumorGrowth <- function(npg = 5,
                           timepoints = c(0, 3, 5, 10),
                           initial_volume = 100,
                           K = 2000,
                           r_Ctrl = 0.08,
                           r_A = 0.07,
                           r_B = 0.06,
                           r_Comb = 0.03,
                           sd = 0.1) {
  
  subject <- 1:(4 * npg) # Subjects' ids
  Treatment <- gl(4, npg, labels = c("Control", "DrugA", "DrugB", "Combination")) # Treatment for each subject
  dts <- data.frame(subject, Treatment) # Subject-level data
  
  dtL <- list(Time = timepoints, subject = subject)
  dtLong <- expand.grid(dtL) # Long format
  mrgDt <- merge(dtLong, dts, sort = FALSE) # Merged
  
  mrgDt$TumorVolume <- NULL
  
  # Simulate exponential growth for each subject
  mrgDt$TumorVolume[mrgDt$Treatment == "Control"] <- log_fun(V0 = initial_volume, K = K, r = r_Ctrl, t = timepoints)
  mrgDt$TumorVolume[mrgDt$Treatment == "DrugA"] <- log_fun(V0 = initial_volume, K = K, r = r_A, t = timepoints)
  mrgDt$TumorVolume[mrgDt$Treatment == "DrugB"] <- log_fun(V0 = initial_volume, K = K, r = r_B, t = timepoints)
  mrgDt$TumorVolume[mrgDt$Treatment == "Combination"] <- log_fun(V0 = initial_volume, K = K, r = r_Comb, t = timepoints)
  
  # Add random noise to simulate variability
  mrgDt$TumorVolume <- mrgDt$TumorVolume * rnorm(nrow(mrgDt), mean = 1, sd = sd)
  mrgDt <- as.data.frame(mrgDt)
  return(mrgDt)
}


## Supplementary Fig. 13a ----

# Example plot
set.seed(123)
grwth_data_1 <- LogTumorGrowth(npg = 5,
                               timepoints = seq(0,30,3),
                               initial_volume = 200,
                               K = 1000,
                               r_Ctrl = 0.25,
                               r_A = 0.15,
                               r_B = 0.1,
                               r_Comb = 0.03,
                               sd = 0.1)

set.seed(123)
grwth_data_1_line <- LogTumorGrowth(npg = 1,
                                    timepoints = seq(0,30,1),
                                    initial_volume = 200,
                                    K = 1000,
                                    r_Ctrl = 0.25,
                                    r_A = 0.15,
                                    r_B = 0.1,
                                    r_Comb = 0.03,
                                    sd = 0)

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")


a <- grwth_data_1 %>% ggplot(aes(x = Time, y = TumorVolume, group = subject)) + geom_point(aes(colour = Treatment), size = 2, alpha = 0.33) +
  geom_line(aes(colour = Treatment), alpha = 0.33) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top") +
  geom_line(data = grwth_data_1_line, aes(colour = Treatment), lwd = 1.5) +
  scale_x_continuous(breaks = unique(grwth_data_1$Time)) +
  xlab("Days after treatment initiation") + ylab("Tumor Volume (mm^3)") +
  labs(title = "Logistic Growth Simulated Data - Synergistic Effect")
a

## Supplementary Fig. 13b ---- 

# Exponential model
lmm1 <- lmmModel(data = grwth_data_1, grwth_model = "exp", sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                 trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination", show_plot = T)

exp_p <- plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Exponential Model") + guides(color = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gompertz model

gomp1 <- lmmModel(data = grwth_data_1, grwth_model = "gompertz", start_values = "selfStart",
                 sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                 trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination",
                 control = nlme::nlmeControl(maxIter = 1000), show_plot = T)
gomp_p <- plot_lmmModel(gomp1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Gompertz Model") + guides(color = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

b <- cowplot::plot_grid(exp_p, gomp_p)
b

## Supplementary Fig. 13c ----

lmm1_resid <- residDiagnostics(lmm1)

c <- cowplot::plot_grid(plotlist = c(lmm1_resid$Plots[c(1,4,5)]), nrow = 1)
c

## Supplementary Fig. 13d ----

gomp1_resid <- residDiagnostics(gomp1)

d <- cowplot::plot_grid(plotlist = c(gomp1_resid$Plots[c(1,4,5)]), nrow = 1)
d

## Supplementary Fig. 13e ----

# Create simulated tumor growths and use exponential growth model

simulatedSyn <- data.frame(Estimate = numeric(0), lwr = numeric(0), upr = numeric(0), pval = numeric(0), Time = numeric(0), Simulation = numeric(0))

for (i in 1:1000) {
  print(i)
  set.seed(i)
  grwth_data_1 <- LogTumorGrowth(npg = 5,
                                 timepoints = seq(0,30,3),
                                 initial_volume = 200,
                                 K = 1000,
                                 r_Ctrl = 0.25,
                                 r_A = 0.15,
                                 r_B = 0.1,
                                 r_Comb = 0.03,
                                 sd = 0.1)
  
  
  ## Build Models
  
  
  lmm1 <- lmmModel(data = grwth_data_1, grwth_model = "exp", sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                   trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination", show_plot = F)
  
  
  ## Synergy
  
  lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T, show_plot = F)
  lmm1_Bliss <- lmm1_Bliss$Synergy
  lmm1_Bliss$Simulation <- i
  
  simulatedSyn <- rbind(simulatedSyn, lmm1_Bliss)
  
}

#write_csv(simulatedSyn, "data/Supp_Fig13/SimulatedSyn1000.csv")

simulatedSyn <- read.csv("data/Supp_Fig13/SimulatedSyn1000.csv")


# Generate plot of summary results

lwrSS <- simulatedSyn %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(lwr))
uprSS <- simulatedSyn %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(upr))
estimateSS <- simulatedSyn %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(Estimate))
pvalSS <- simulatedSyn %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(pval))

simSS <- cbind(lwrSS, uprSS[,-1], estimateSS[,-1], pvalSS[,-1])

simSS$Metric <- "SS"

lwrCI <- simulatedSyn %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(lwr))
uprCI <- simulatedSyn %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(upr))
estimateCI <- simulatedSyn %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(Estimate))
pvalCI <- simulatedSyn %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(pval))

simCI <- cbind(lwrCI, uprCI[,-1], estimateCI[,-1], pvalCI[,-1])
simCI$Metric <- "CI"

simSimulated <- rbind(simSS, simCI)


e <- simSS %>% ggplot(aes(x = Time, y = `median(Estimate)`)) +
  geom_segment(aes(x= Time, y = `median(lwr)`, yend = `median(upr)`), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = (`median(pval)`)), shape = 23, size = 5, color = "gray65") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", limits = c(0,0.5)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  scale_x_continuous(breaks = unique(simSS$Time)) +
  labs(title = "Bliss Sinergy for Logistic Growth Simulated Data", subtitle = "Exponential Model") +
  geom_text(aes(x = Time, y = `median(upr)`+0.25, label = paste0("p=",round(`median(pval)`, 3)), angle = 45,
                hjust = 0)) +
  annotate(geom = "text", x = 0.5, 
           y = 0.5, angle =90, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -0.5, angle = 90, hjust = 1, label = "Antagonism", fontface = "bold", color = "#c21d2f")
e
format(simSS$`median(pval)`, scientific = F)

## Supplementary Fig. 13f ----

# Create simulated tumor growths and use Gompertz growth

simulatedGomp <- data.frame(Estimate = numeric(0), lwr = numeric(0), upr = numeric(0), pval = numeric(0), Time = numeric(0), Simulation = numeric(0))

error_iteration <- c()

i <- 1

while (nrow(simulatedGomp) < 1000*20) {
  print(i)
  set.seed(i)
  grwth_data_1 <- LogTumorGrowth(
    npg = 5,
    timepoints = seq(0, 30, 3),
    initial_volume = 200,
    K = 1000,
    r_Ctrl = 0.25,
    r_A = 0.15,
    r_B = 0.1,
    r_Comb = 0.03,
    sd = 0.1
  )
  
  
  ## Build Models
  
  result <- tryCatch({
    lmm1 <- lmmModel(
      data = grwth_data_1,
      grwth_model = "gompertz",
      start_values = c(0.1, 0.01),
      sample_id = "subject",
      time = "Time",
      treatment = "Treatment",
      tumor_vol = "TumorVolume",
      trt_control = "Control",
      drug_a = "DrugA",
      drug_b = "DrugB",
      combination = "Combination",
      control = nlme::nlmeControl(maxIter = 1000),
      show_plot = F
    )
    
    
    ## Synergy
    
    lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", show_plot = F)
    lmm1_Bliss <- lmm1_Bliss$Synergy
    lmm1_Bliss$Simulation <- i
    
    simulatedGomp <- rbind(simulatedGomp, lmm1_Bliss)
    
    NULL # No error, return NULL
  }, error = function(e) {
    message(paste("Error in iteration", i, ":", conditionMessage(e)))
    error_iteration <- c(error_iteration, i)
    NULL # Continue loop
  })
  i <- i + 1 # update seed
}

# 1009 iterations
# Failure rate

9/1009*100

#write_csv(simulatedGomp, "data/Supp_Fig13/SimulatedSynGomp1000.csv")

SimulatedSynGomp <- read.csv("data/Supp_Fig13/SimulatedSynGomp1000.csv")


lwrSS <- SimulatedSynGomp %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(lwr))
uprSS <- SimulatedSynGomp %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(upr))
estimateSS <- SimulatedSynGomp %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(Estimate))
pvalSS <- SimulatedSynGomp %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(pval))

simSS <- cbind(lwrSS, uprSS[,-1], estimateSS[,-1], pvalSS[,-1])

simSS$Metric <- "SS"

lwrCI <- SimulatedSynGomp %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(lwr))
uprCI <- SimulatedSynGomp %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(upr))
estimateCI <- SimulatedSynGomp %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(Estimate))
pvalCI <- SimulatedSynGomp %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(pval))

simCI <- cbind(lwrCI, uprCI[,-1], estimateCI[,-1], pvalCI[,-1])
simCI$Metric <- "CI"

simSimulated <- rbind(simSS, simCI)


f <- simSS %>% ggplot(aes(x = Time, y = `median(Estimate)`)) +
  geom_segment(aes(x= Time, y = `median(lwr)`, yend = `median(upr)`), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = (`median(pval)`)), shape = 23, size = 5, color = "gray65") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", limits = c(0, 0.5)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  scale_x_continuous(breaks = unique(simSS$Time)) +
  labs(title = "Bliss Sinergy for Logistic Growth Simulated Data", subtitle = "Gompertz Model") +
  geom_text(aes(x = Time, y = `median(upr)`+0.25, label = paste0("p=",round(`median(pval)`, 3)), angle = 45,
                hjust = 0)) +
  annotate(geom = "text", x = 0.5, 
           y = 0.5, angle =90, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -0.5, angle = 90, hjust = 1, label = "Antagonism", fontface = "bold", color = "#c21d2f")
f

ab <- cowplot::plot_grid(a, b, nrow = 1, rel_widths = c(1,2), labels = "auto")
ef <- cowplot::plot_grid(e,f, nrow = 1, labels = c("e", "f"))

SuppFig13 <- cowplot::plot_grid(plotlist = list(ab, c, d, ef), ncol = 1, rel_heights = c(1.2,1,1,1.2),
                                labels = c("", "c", "d", ""))
SuppFig13
#ggsave("figures/SuppFig13.pdf", SuppFig13, height = 2*height, width = 2.33*width)

# Supplementary Fig. 14 ----

## Supplementary Fig. 14a ----

set.seed(123)
grwth_data_1 <- LogTumorGrowth(npg = 5,
                               timepoints = seq(0,30,3),
                               initial_volume = 200,
                               K = 1000,
                               r_Ctrl = 0.25,
                               r_A = 0.15,
                               r_B = 0.1,
                               r_Comb = 0.06,
                               sd = 0.1)

set.seed(123)
grwth_data_1_line <- LogTumorGrowth(npg = 1,
                                    timepoints = seq(0,30,1),
                                    initial_volume = 200,
                                    K = 1000,
                                    r_Ctrl = 0.25,
                                    r_A = 0.15,
                                    r_B = 0.1,
                                    r_Comb = 0.06,
                                    sd = 0)

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

a <- grwth_data_1 %>% ggplot(aes(x = Time, y = TumorVolume, group = subject)) + geom_point(aes(colour = Treatment), size = 2, alpha = 0.33) +
  geom_line(aes(colour = Treatment), alpha = 0.33) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top") +
  geom_line(data = grwth_data_1_line, aes(colour = Treatment), lwd = 1.5) +
  scale_x_continuous(breaks = unique(grwth_data_1$Time)) +
  xlab("Days after treatment initiation") + ylab("Tumor Volume (mm^3)") +
  labs(title = "Logistic Growth Simulated Data - Additive Effect")
a

## Supplementary Fig. 14b ----

# Exponential model
lmm1 <- lmmModel(data = grwth_data_1, grwth_model = "exp", sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                 trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination", show_plot = T)

exp_p <- plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Exponential Model") + guides(color = "none")

# Gompertz model

gomp1 <- lmmModel(data = grwth_data_1, grwth_model = "gompertz", start_values = "selfStart",
                  sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                  trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination",
                  control = nlme::nlmeControl(maxIter = 1000), show_plot = T)
gomp_p <- plot_lmmModel(gomp1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Gompertz Model") + guides(color = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

b <- cowplot::plot_grid(exp_p, gomp_p)
b

## Supplementary Fig. 14c ----

lmm1_resid <- residDiagnostics(lmm1)

c <- cowplot::plot_grid(plotlist = c(lmm1_resid$Plots[c(1,4,5)]), nrow = 1)
c

## Supplementary Fig. 14d ----

gomp1_resid <- residDiagnostics(gomp1)

d <- cowplot::plot_grid(plotlist = c(gomp1_resid$Plots[c(1,4,5)]), nrow = 1)
d

## Supplementary Fig. 14e ----

simulatedAdd <- data.frame(Estimate = numeric(0), lwr = numeric(0), upr = numeric(0), pval = numeric(0), Time = numeric(0), Simulation = numeric(0))

for (i in 1:1000) {
  print(i)
  set.seed(i)
  grwth_data_1 <- LogTumorGrowth(npg = 5,
                                 timepoints = seq(0,30,3),
                                 initial_volume = 200,
                                 K = 1000,
                                 r_Ctrl = 0.25,
                                 r_A = 0.15,
                                 r_B = 0.1,
                                 r_Comb = 0.06,
                                 sd = 0.1)
  
  
  ## Build Models
  
  
  lmm1 <- lmmModel(data = grwth_data_1, grwth_model = "exp", sample_id = "subject", time = "Time", treatment = "Treatment", tumor_vol = "TumorVolume",
                   trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination", show_plot = F)
  
  
  ## Synergy
  
  lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T, show_plot = F)
  lmm1_Bliss <- lmm1_Bliss$Synergy
  lmm1_Bliss$Simulation <- i
  
  simulatedAdd <- rbind(simulatedAdd, lmm1_Bliss)
  
}

#write_csv(simulatedAdd, "data/Supp_Fig14/SimulatedAdd1000.csv")

simulatedAdd <- read.csv("data/Supp_Fig14/SimulatedAdd1000.csv")


# Generate plot of summary results

lwrSS <- simulatedAdd %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(lwr))
uprSS <- simulatedAdd %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(upr))
estimateSS <- simulatedAdd %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(Estimate))
pvalSS <- simulatedAdd %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(pval))

simSS <- cbind(lwrSS, uprSS[,-1], estimateSS[,-1], pvalSS[,-1])

simSS$Metric <- "SS"

lwrCI <- simulatedAdd %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(lwr))
uprCI <- simulatedAdd %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(upr))
estimateCI <- simulatedAdd %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(Estimate))
pvalCI <- simulatedAdd %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(pval))

simCI <- cbind(lwrCI, uprCI[,-1], estimateCI[,-1], pvalCI[,-1])
simCI$Metric <- "CI"

simSimulated <- rbind(simSS, simCI)


e <- simSS %>% ggplot(aes(x = Time, y = `median(Estimate)`)) +
  geom_segment(aes(x= Time, y = `median(lwr)`, yend = `median(upr)`), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = (`median(pval)`)), shape = 23, size = 5, color = "gray65") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", limits = c(0, 0.5)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  scale_x_continuous(breaks = unique(simSS$Time)) +
  labs(title = "Bliss Sinergy for Logistic Growth Simulated Data", subtitle = "Exponential Model") +
  geom_text(aes(x = Time, y = `median(upr)`+0.25, label = paste0("p=",round(`median(pval)`, 3)), angle = 45,
                hjust = 0)) +
  annotate(geom = "text", x = 0.5, 
           y = 0.5, angle =90, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -0.5, angle = 90, hjust = 1, label = "Antagonism", fontface = "bold", color = "#c21d2f")
e

## Supplementary Fig. 14f ----

# Create simulated tumor growths and use Gompertz growth

simulatedAdd_Gomp <- data.frame(Estimate = numeric(0), lwr = numeric(0), upr = numeric(0), pval = numeric(0), Time = numeric(0), Simulation = numeric(0))

error_iteration <- c()

i <- 1

while (nrow(simulatedAdd_Gomp) < 1000*20) {
  print(i)
  set.seed(i)
  grwth_data_1 <- LogTumorGrowth(
    npg = 5,
    timepoints = seq(0, 30, 3),
    initial_volume = 200,
    K = 1000,
    r_Ctrl = 0.25,
    r_A = 0.15,
    r_B = 0.1,
    r_Comb = 0.06,
    sd = 0.1
  )
  
  
  ## Build Models
  
  result <- tryCatch({
    lmm1 <- lmmModel(
      data = grwth_data_1,
      grwth_model = "gompertz",
      start_values = c(0.1, 0.01),
      sample_id = "subject",
      time = "Time",
      treatment = "Treatment",
      tumor_vol = "TumorVolume",
      trt_control = "Control",
      drug_a = "DrugA",
      drug_b = "DrugB",
      combination = "Combination",
      control = nlme::nlmeControl(maxIter = 1000),
      show_plot = F
    )
    
    
    ## Synergy
    
    lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", show_plot = F)
    lmm1_Bliss <- lmm1_Bliss$Synergy
    lmm1_Bliss$Simulation <- i
    
    simulatedAdd_Gomp <- rbind(simulatedAdd_Gomp, lmm1_Bliss)
    
    NULL # No error, return NULL
  }, error = function(e) {
    message(paste("Error in iteration", i, ":", conditionMessage(e)))
    error_iteration <- c(error_iteration, i)
    NULL # Continue loop
  })
  i <- i + 1 # update seed
}

# Failure rate 

8/1008*100

#write_csv(simulatedAdd_Gomp, "data/Supp_Fig14/SimulatedAdd_Gomp.csv")

SimulatedAdd_Gomp <- read.csv("data/Supp_Fig14/SimulatedAdd_Gomp.csv")

# Generate plot of summary results

lwrSS <- SimulatedAdd_Gomp %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(lwr))
uprSS <- SimulatedAdd_Gomp %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(upr))
estimateSS <- SimulatedAdd_Gomp %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(Estimate))
pvalSS <- SimulatedAdd_Gomp %>% filter(Metric == "SS") %>% group_by(Time) %>% summarise(median(pval))

simSS <- cbind(lwrSS, uprSS[,-1], estimateSS[,-1], pvalSS[,-1])

simSS$Metric <- "SS"

lwrCI <- SimulatedAdd_Gomp %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(lwr))
uprCI <- SimulatedAdd_Gomp %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(upr))
estimateCI <- SimulatedAdd_Gomp %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(Estimate))
pvalCI <- SimulatedAdd_Gomp %>% filter(Metric == "CI") %>% group_by(Time) %>% summarise(median(pval))

simCI <- cbind(lwrCI, uprCI[,-1], estimateCI[,-1], pvalCI[,-1])
simCI$Metric <- "CI"

simSimulated <- rbind(simSS, simCI)


f <- simSS %>% ggplot(aes(x = Time, y = `median(Estimate)`)) +
  geom_segment(aes(x= Time, y = `median(lwr)`, yend = `median(upr)`), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = (`median(pval)`)), shape = 23, size = 5, color = "gray65") +
  scale_fill_gradient(name = "p-value", high = "darkorchid4", low = "#d0f5ec", limits = c(0, 0.5)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  scale_x_continuous(breaks = unique(simSS$Time)) +
  labs(title = "Bliss Sinergy for Logistic Growth Simulated Data", subtitle = "Gompertz Model") +
  geom_text(aes(x = Time, y = `median(upr)`+0.25, label = paste0("p=",round(`median(pval)`, 3)), angle = 45,
                hjust = 0)) +
  annotate(geom = "text", x = 0.5, 
           y = 0.5, angle =90, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
  annotate(geom = "text", x = 0.5, 
           y = -0.5, angle = 90, hjust = 1, label = "Antagonism", fontface = "bold", color = "#c21d2f")
f

ab <- cowplot::plot_grid(a, b, nrow = 1, rel_widths = c(1,2), labels = "auto")
ef <- cowplot::plot_grid(e,f, nrow = 1, labels = c("e", "f"))

SuppFig14 <- cowplot::plot_grid(plotlist = list(ab, c, d, ef), ncol = 1, rel_heights = c(1.2,1,1,1.2),
                                labels = c("", "c", "d", ""))
SuppFig14
#ggsave("figures/SuppFig14.pdf", SuppFig14, height = 2*height, width = 2.33*width)

# Figure 6 ----

### Post-Hoc Power ----
# Note: the Post-Hoc power calculation can take some time to run

# BV-173-Gluc

set.seed(123)
PostHocPwr(BV173Gluc_lmm, method = "Bliss")
# [1] 0.812

set.seed(123)
PostHocPwr(BV173Gluc_lmm, method = "HSA")
# [1] 0.96

BV173Gluc <- lmmModel_estimates(BV173Gluc_lmm)

BV173Gluc$PwrBliss <- 0.812
BV173Gluc$PwrHSA <- 0.96
#write.csv(BV173Gluc, "data/Fig6/BV173Gluc_model_estimates.csv")

BV173Gluc <- read.csv("data/Fig6/BV173Gluc_model_estimates.csv")
BV173Gluc <- BV173Gluc[,c(2,4,6,8,10:13)]
BV173Gluc$Cell <- "BV-173-Gluc"
colnames(BV173Gluc) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# U87-MG FM

set.seed(123)
PostHocPwr(U87MG_lmm, method = "Bliss")
# [1] 0.271

set.seed(123)
PostHocPwr(U87MG_lmm, method = "HSA")
# [1] 0.993

U87MG <- lmmModel_estimates(U87MG_lmm)

U87MG$PwrBliss <- 0.271
U87MG$PwrHSA <- 0.993
#write.csv(U87MG, "data/Fig6/U87MG_model_estimates.csv")


U87MG <- read.csv("data/Fig6/U87MG_model_estimates.csv")
U87MG <- U87MG[,c(2,4,6,8,10:13)]
U87MG$Cell <- "U-87-MG-FM"
colnames(U87MG) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# CR1197

set.seed(123)
PostHocPwr(CR1197_lmm, method = "Bliss")
# [1] 0.658

set.seed(123)
PostHocPwr(CR1197_lmm, method = "HSA")
# [1] 0.998

CR1197 <- lmmModel_estimates(CR1197_lmm)

CR1197$PwrBliss <- 0.658
CR1197$PwrHSA <- 0.998
#write.csv(CR1197, "data/Fig6/CR1197_model_estimates.csv")


CR1197 <- read.csv("data/Fig6/CR1197_model_estimates.csv")
CR1197 <- CR1197[,c(2,4,6,8,10:13)]
CR1197$Cell <- "CR1197"
colnames(CR1197) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# LS-1034
set.seed(123)
PostHocPwr(LS1034_lmm, method = "Bliss")
# [1] 0.896
set.seed(123)
PostHocPwr(LS1034_lmm, method = "HSA")
# [1] 0.994

LS1034 <- lmmModel_estimates(LS1034_lmm)

LS1034$PwrBliss <- 0.896
LS1034$PwrHSA <- 0.994
#write.csv(LS1034, "data/Fig6/LS1034_model_estimates.csv")


LS1034 <- read.csv("data/Fig6/LS1034_model_estimates.csv")
LS1034 <- LS1034[,c(2,4,6,8,10:13)]
LS1034$Cell <- "LS-1034"
colnames(LS1034) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# SW837

set.seed(123)
PostHocPwr(SW837_lmm, method = "Bliss")
# [1] 0.458

set.seed(123)
PostHocPwr(SW837_lmm, method = "HSA")
# [1] 0.74

SW837 <- lmmModel_estimates(SW837_lmm)

SW837$PwrBliss <- 0.458
SW837$PwrHSA <- 0.740
#write.csv(SW837, "data/Fig6/SW837_model_estimates.csv")


SW837 <- read.csv("data/Fig6/SW837_model_estimates.csv")
SW837 <- SW837[,c(2,4,6,8,10:13)]
SW837$Cell <- "SW837"
colnames(SW837) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# SNU81

set.seed(123)
PostHocPwr(SNU81_lmm, method = "Bliss")
# [1] 0.443

set.seed(123)
PostHocPwr(SNU81_lmm, method = "HSA")
# [1] 0.597

SNU81 <- lmmModel_estimates(SNU81_lmm)

SNU81$PwrBliss <- 0.443
SNU81$PwrHSA <- 0.597
#write.csv(SNU81, "data/Fig6/SNU81_model_estimates.csv")


SNU81 <- read.csv("data/Fig6/SNU81_model_estimates.csv")
SNU81 <- SNU81[,c(2,4,6,8,10:13)]

SNU81$Cell <- "SNU-81"
colnames(SNU81) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# 4T1

set.seed(123)
PostHocPwr(GABA_lmm, method = "Bliss")
# [1] 0.385
set.seed(123)
PostHocPwr(GABA_lmm, method = "HSA")
# [1] 0.667

GABA <- lmmModel_estimates(GABA_lmm)

GABA$PwrBliss <- 0.385
GABA$PwrHSA <- 0.667
#write.csv(GABA, "data/Fig6/GABA_model_estimates.csv")


GABA <- read.csv("data/Fig6/GABA_model_estimates.csv")
GABA <- GABA[,c(2,4,6,8,10:13)]

GABA$Cell <- "4T1"
colnames(GABA) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# MAS98.06 Heregulin-Fulvestrant

set.seed(123)
PostHocPwr(MAS98.06_Her_Fulv_lmm, method = "Bliss")
# [1] 0.151

set.seed(123)
PostHocPwr(MAS98.06_Her_Fulv_lmm, method = "HSA")
# [1] 0.655

Her_Fulv <- lmmModel_estimates(MAS98.06_Her_Fulv_lmm)

Her_Fulv$PwrBliss <- 0.151
Her_Fulv$PwrHSA <- 0.655
#write.csv(Her_Fulv, "data/Fig6/MAS98.06_Her_Fulv_model_estimates.csv")


Her_Fulv <- read.csv("data/Fig6/MAS98.06_Her_Fulv_model_estimates.csv")
Her_Fulv <- Her_Fulv[,c(2,4,6,8,10:13)]

Her_Fulv$Cell <- "MAS98.06 (Her+Fulv)"
colnames(Her_Fulv) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# MAS98.06 EGF-Fulvestrant

set.seed(123)
PostHocPwr(MAS98.06_EGF_Fulv_lmm, method = "Bliss")
# [1] 0.194

set.seed(123)
PostHocPwr(MAS98.06_EGF_Fulv_lmm, method = "HSA")
# [1] 0.489

EGF_Fulv <- lmmModel_estimates(MAS98.06_EGF_Fulv_lmm)

EGF_Fulv$PwrBliss <- 0.194
EGF_Fulv$PwrHSA <- 0.489
#write.csv(EGF_Fulv, "data/Fig6/MAS98.06_EGF_Fulv_model_estimates.csv")


EGF_Fulv <- read.csv("data/Fig6/MAS98.06_EGF_Fulv_model_estimates.csv")
EGF_Fulv <- EGF_Fulv[,c(2,4,6,8,10:13)]

EGF_Fulv$Cell <- "MAS98.06 (EGF+Fulv)"
colnames(EGF_Fulv) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# MAS98.06 Trastuzumab - Fulvestrant

set.seed(123)
PostHocPwr(MAS98.06_Fulv_Trast_lmm, method = "Bliss")
# [1] 0.157

set.seed(123)
PostHocPwr(MAS98.06_Fulv_Trast_lmm, method = "HSA")
# [1] 0.09


Trast_Fulv <- lmmModel_estimates(MAS98.06_Fulv_Trast_lmm)

Trast_Fulv$PwrBliss <- 0.157
Trast_Fulv$PwrHSA <- 0.09
#write.csv(Trast_Fulv, "data/Fig6/MAS98.06_Trast_Fulv_model_estimates.csv")


Trast_Fulv <- read.csv("data/Fig6/MAS98.06_Trast_Fulv_model_estimates.csv")
Trast_Fulv <- Trast_Fulv[,c(2,4,6,8,10:13)]

Trast_Fulv$Cell <- "MAS98.06 (Trast+Fulv)"
colnames(Trast_Fulv) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# MAS98.06 Heregulin - Trastuzumab
set.seed(123)
PostHocPwr(MAS98.06_Her_Trast_lmm, method = "Bliss")
# [1] 0.145

set.seed(123)
PostHocPwr(MAS98.06_Her_Trast_lmm, method = "HSA")
# [1] 0.079

Her_Trast <- lmmModel_estimates(MAS98.06_Her_Trast_lmm)

Her_Trast$PwrBliss <- 0.145
Her_Trast$PwrHSA <- 0.079
#write.csv(Her_Trast, "data/Fig6/MAS98.06_Her_Trast_model_estimates.csv")


Her_Trast <- read.csv("data/Fig6/MAS98.06_Her_Trast_model_estimates.csv")
Her_Trast <- Her_Trast[,c(2,4,6,8,10:13)]

Her_Trast$Cell <- "MAS98.06 (Her+Trast)"
colnames(Her_Trast) <- c("Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# MAS98.06 Heregulin - Trastuzumab - Fulvestrant
set.seed(123)
PostHocPwr(MAS98.06_Her_Fulv_Trast_lmm, method = "Bliss")
# [1] 0.084

set.seed(123)
PostHocPwr(MAS98.06_Her_Fulv_Trast_lmm, method = "HSA")
# [1] 0.21

Her_Trast_Fulv <- lmmModel_estimates(MAS98.06_Her_Fulv_Trast_lmm)

Her_Trast_Fulv$PwrBliss <- 0.084
Her_Trast_Fulv$PwrHSA <- 0.21
#write.csv(Her_Trast_Fulv, "data/Fig6/MAS98.06_Her_Trast_Fulv_model_estimates.csv")


Her_Trast_Fulv <- read.csv("data/Fig6/MAS98.06_Her_Trast_Fulv_model_estimates.csv")
Her_Trast_Fulv <- Her_Trast_Fulv[,c(2,4,6,8,10, 12:15)]

Her_Trast_Fulv$Cell <- "MAS98.06 (Her+Fulv+Trast)"
colnames(Her_Trast_Fulv) <- c("Ctrl", "DrugA", "DrugB", "DrugC", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# U87-MG FM AZD - OSi - Doc

set.seed(123)
PostHocPwr(triple_lmm, method = "Bliss")
# 0.217

set.seed(123)
PostHocPwr(triple_lmm, method = "HSA")
#  0.107

U87MG_Triple <- lmmModel_estimates(triple_lmm)

U87MG_Triple$PwrBliss <- 0.217
U87MG_Triple$PwrHSA <- 0.107
#write.csv(U87MG_Triple, "data/Fig6/U87MG_Triple_model_estimates.csv")


U87MG_Triple <- read.csv("data/Fig6/U87MG_Triple_model_estimates.csv")
U87MG_Triple <- U87MG_Triple[,c(2,4,6,8,10, 12:15)]

U87MG_Triple$Cell <- "U-87-MG-FM (AZD+Osi+Doc)"
colnames(U87MG_Triple) <-c("Ctrl", "DrugA", "DrugB", "DrugC", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

unit_norm_l1 <- function(x) {
  x / sum(abs(x))  # Divide by the L1 norm
  }


Pwr_df <- rbind(BV173Gluc, CR1197, EGF_Fulv, GABA, Her_Fulv, Her_Trast, LS1034, SNU81, SW837, Trast_Fulv,U87MG)

rownames(Pwr_df) <- Pwr_df$Cell
Pwr_bliss_hsa <- Pwr_df[,c(7:9)]
Pwr_df <- Pwr_df[,-c(7:9)]

Pwr_df_norm <- t(apply(Pwr_df, 1, unit_norm_l1))

Pwr_df_norm <- Pwr_df_norm[,c(4:6)]

Pwr_df_triple <- rbind(Her_Trast_Fulv, U87MG_Triple)
rownames(Pwr_df_triple) <- Pwr_df_triple$Cell
Pwr_bliss_hsa <- rbind(Pwr_bliss_hsa, Pwr_df_triple[,c(8:10)])
Pwr_df_triple <- Pwr_df_triple[,-c(8:10)]

Pwr_df_triple_norm <- t(apply(Pwr_df_triple, 1, unit_norm_l1))

Pwr_df_triple_norm <- Pwr_df_triple_norm[,5:7]


Pwr_df_norm <- rbind(Pwr_df_norm, Pwr_df_triple_norm)
Pwr_df_norm <- as.data.frame(Pwr_df_norm)

Pwr_df_norm$PwrBliss <- Pwr_bliss_hsa$PwrBliss
Pwr_df_norm$PwrHSA <- Pwr_bliss_hsa$PwrHSA
Pwr_df_norm$Cell <- Pwr_bliss_hsa$Cell


## Fig. 6a ----

a <- Pwr_df_norm %>% ggplot(aes(x = sd_ranef, y = sd_resid, group = Cell)) +
  geom_point(aes(size = PwrBliss, colour = Combination)) + 
  ggrepel::geom_text_repel(aes(label = Cell),size = 3, box.padding = 0.5,segment.color = "gray25", max.overlaps = Inf) + theme_cowplot() +
  scale_color_gradient2(name = "Normalized\nCombination Growth Rate",low = "purple3",mid = "#f7f7f7", high = "#e66101") +
  scale_size_area(name = "Statistical Power", max_size = 10, limits = c(0,1)) +
  xlab("Normalized SD of Random Effects") + ylab("Normalized SD of Residuals") +
  labs(title = "Power for Bliss Synergy") #+ ylim(0, 0.6)
a

## Fig. 6b ----

df <- SW837_lmm$dt1

plot_df <- SW837_lmm$data %>%
  dplyr::mutate(fixed_effect = predict(SW837_lmm, level = 0))  # level = 0 = fixed only

# Add manually data at time zero
initial_points <- SW837_lmm$data %>%
  dplyr::distinct(.data$SampleID, .data$Treatment) %>%
  dplyr::mutate(Time = 0, logRTV = 0)
initial_points <- initial_points %>%
  dplyr::mutate(fixed_effect = 0)

plot_df <- dplyr::bind_rows(plot_df, initial_points)

b <- df %>% ggplot(aes(x = Time, y = logRTV)) + geom_point(aes(colour = Treatment, group = SampleID), shape = 19, alpha = 0.33) +
  geom_line(aes(colour = Treatment, group = SampleID), alpha = 0.25) +
  scale_color_manual(values = c("#3c3c3b", "#d50c52", "#00a49c","#601580")) +
  geom_line(data = plot_df, aes(y = fixed_effect, group = SampleID, color = Treatment), 
            lwd = 1.5, alpha = 1)  + theme_cowplot() + theme(legend.position = "top") +
  labs(title = "SW837 Rabusertib-Irinotecan") +
  scale_x_continuous(breaks = unique(df$Time)) +
  xlab("Days after treatment initiation") + ylab("log (Relative Luminescence Units)")

b <- plot_grid(b, NULL, ncol = 1, rel_heights = c(1,0.25))
b

## Fig. 6c and 6d ----

length(unique(SW837_lmm$dt1$SampleID))/4
est <- lmmModel_estimates(SW837_lmm)

APrioriPwr(
  npg = 6,
  time = unique(SW837_lmm$dt1$Time),
  grwrControl = est$Ctrl,
  grwrA = est$Rabusertib,
  grwrB = est$Irinotecan,
  grwrComb = est$Combination,
  sd_ranef = est$sd_ranef,
  sgma = est$sd_resid,
  sd_eval = seq(0.1 * est$sd_ranef, 2 * est$sd_ranef, 0.05 * est$sd_ranef),
  sgma_eval = seq(0.1 * est$sd_resid, 2 * est$sd_resid, 0.05 * est$sd_resid),
  grwrComb_eval = seq(-5 * est$Combination, 16 * est$Combination, 0.1 * est$Combination)
)
#ggsave("figures/Figure6cd.pdf", height = 0.8 * height, width = 3*width)

## Fig. 6e ----

PwrSampleSize(
  npg = 4:20,
  time = unique(SW837_lmm$dt1$Time),
  grwrControl = est$Ctrl,
  grwrA = est$Rabusertib,
  grwrB = est$Irinotecan,
  grwrComb = est$Combination,
  sd_ranef = est$sd_ranef,
  sgma = est$sd_resid
)
#ggsave("figures/Figure6e.pdf", height = 0.8 * height, width = 2*width)

## Fig. 6f ----

df <- LS1034_lmm$dt1

plot_df <- LS1034_lmm$data %>%
  dplyr::mutate(fixed_effect = predict(LS1034_lmm, level = 0))  # level = 0 = fixed only

# Add manually data at time zero
initial_points <- LS1034_lmm$data %>%
  dplyr::distinct(.data$SampleID, .data$Treatment) %>%
  dplyr::mutate(Time = 0, logRTV = 0)
initial_points <- initial_points %>%
  dplyr::mutate(fixed_effect = 0)

plot_df <- dplyr::bind_rows(plot_df, initial_points)

f <- df %>% ggplot(aes(x = Time, y = logRTV)) + geom_point(aes(colour = Treatment, group = SampleID), shape = 19, alpha = 0.33) +
  geom_line(aes(colour = Treatment, group = SampleID), alpha = 0.25) +
  scale_color_manual(values = c("#3c3c3b", "#d50c52", "#00a49c","#601580")) +
  geom_line(data = plot_df, aes(y = fixed_effect, group = SampleID, color = Treatment), 
               lwd = 1.5, alpha = 1)  + 
  theme_cowplot() + theme(legend.position = "top") +
  labs(title = "LS-1031 Rabusertib - Irinotecan") + xlab("Days after treatment initiation") + ylab("log (Relative Tumor Volume)") +
  scale_x_continuous(breaks = unique(df$Time))
f

## Fig. 6g ----

PwrTime(
  npg = 10,
  time = list(
    seq(0, 9, 3),
    seq(0, 12, 3),
    seq(0, 15, 3),
    seq(0, 18, 3),
    seq(0, 21, 3),
    seq(0, 24, 3),
    seq(0, 27, 3),
    seq(0, 30, 3),
    seq(0, 33, 3)
  ),
  grwrControl = 0.097,
  grwrA = 0.103,
  grwrB = -0.006,
  grwrComb = -0.037,
  sd_ranef = 0.028,
  sgma = 0.402
)
#ggsave("figures/Figure6g.pdf", height = 0.8 * height, width = 2*width)

## Fig. 6h ----


PwrTime(
  npg = 10,
  time = list(
    seq(0, 18, 1),
    seq(0, 18, 2),
    seq(0, 18, 3),
    seq(0, 18, 6),
    seq(0, 18, 9),
    seq(0, 18, 18)
  ),
  grwrControl = 0.097,
  grwrA = 0.103,
  grwrB = -0.006,
  grwrComb = -0.037,
  sd_ranef = 0.028,
  sgma = 0.402,
  type = "freq"
)
#ggsave("figures/Figure6h.pdf", height = 0.8 * height, width = 2*width)

ab <- plot_grid(a, b, nrow = 1, labels = "auto")
ab

cde <- plot_grid(NULL, NULL, NULL, nrow = 1, labels = c("c", "d", "e"))

fgh <- plot_grid(f, NULL, NULL, nrow = 1, labels = c("f", "g", "h"), rel_widths = c(1.33,1,1))

Fig6 <- plot_grid(ab,cde, fgh, ncol = 1, rel_heights = c(1.2,1,1))
Fig6
#ggsave("figures/Figure6.pdf", Fig6, height = 1.66*height, width = 2*width)

# Supplementary Table 1 ----

# The code below generates the exponential and Gompertz models for all datasets,
# together with the model diagnostics and metrics used for the comparison of
# models shown in Supplementary Table 1.

### Fig. 2a ----

U87MG <- read.xlsx("data/Fig2/U87MG.xlsx")
U87MG$Day <- sqrt(U87MG$Day-12)

U87MG_lmm <- lmmModel(data = U87MG, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                      trt_control = "Control", drug_a = "Docetaxel", drug_b = "GNE-317", combination = "Docetaxel+GNE-317",
                      weights = varPower(form = ~Time))
U87MG_lmm_gomp <- lmmModel(data = U87MG, grwth_model = "gompertz", sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                      trt_control = "Control", drug_a = "Docetaxel", drug_b = "GNE-317", combination = "Docetaxel+GNE-317",
                      weights = varPower(form = ~Time))

ObsvsPred(U87MG_lmm)

### Fig. 2c ----

BV173Gluc <- read.csv("data/Fig2/BV_173.csv")
BV173Gluc$Day <- sqrt(BV173Gluc$Day-7)

BV173Gluc_lmm <- lmmModel(data = BV173Gluc, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                          trt_control = "Control", drug_a = "Imatinib", drug_b = "Dasatinib", combination = "Imatinib_Dasatinib",
                          weights = nlme::varIdent(form = ~ 1|Treatment))

BV173Gluc_lmm_gomp <- lmmModel(data = BV173Gluc, grwth_model = "gompertz", sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                          trt_control = "Control", drug_a = "Imatinib", drug_b = "Dasatinib", combination = "Imatinib_Dasatinib",
                          weights = nlme::varIdent(form = ~ 1|Treatment),
                          control = nlmeControl(maxIter = 5000, msMaxIter = 1000),
                          start_values = "selfStart")

ranefDiagnostics(BV173Gluc_lmm)
ranefDiagnostics(BV173Gluc_lmm_gomp)

residDiagnostics(BV173Gluc_lmm)
residDiagnostics(BV173Gluc_lmm_gomp)

ObsvsPred(BV173Gluc_lmm)
ObsvsPred(BV173Gluc_lmm_gomp)

### Fig. 2e ----

CHL1FM <- read.xlsx("data/Fig2/CHL1_FM.xlsx")
CHL1FM$Day <- sqrt(CHL1FM$Day-7)

CHL1FM_lmm <- lmmModel(data = CHL1FM, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                       trt_control = "Control", drug_a = "Gemcitabine", drug_b = "CGP_082996", combination = "Gemcitabine_CGP",
                       weights = nlme::varIdent(form = ~1|SampleID),
                       control = lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000))

CHL1FM_lmm_gomp <- lmmModel(data = CHL1FM, grwth_model = "gompertz", sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                            trt_control = "Control", drug_a = "Gemcitabine", drug_b = "CGP_082996", combination = "Gemcitabine_CGP",
                            weights = nlme::varIdent(form = ~1|SampleID),
                            control = lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000),
                            start_values = "selfStart")

ranefDiagnostics(CHL1FM_lmm)
ranefDiagnostics(CHL1FM_lmm_gomp)

residDiagnostics(CHL1FM_lmm)
residDiagnostics(CHL1FM_lmm_gomp)

ObsvsPred(CHL1FM_lmm)
ObsvsPred(CHL1FM_lmm_gomp)

### Fig. 2g ----

MDAMD231 <- read.xlsx("data/Fig2/MDA_MD231_Fig5C.xlsx")
MDAMD231$Day <- sqrt(MDAMD231$Day-16) 

MDAMD231_lmm <- lmmModel(data = MDAMD231, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                         trt_control = "Control", drug_a = "AZ628", drug_b = "Gemcitabine", combination = "AZ628_Gemcitabine",
                         weights = nlme::varIdent(form = ~ 1|Treatment))

MDAMD231_lmm_gomp <- lmmModel(data = MDAMD231, grwth_model = "gompertz", sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                              trt_control = "Control", drug_a = "AZ628", drug_b = "Gemcitabine", combination = "AZ628_Gemcitabine",
                              weights = nlme::varIdent(form = ~ 1|Treatment))
ranefDiagnostics(MDAMD231_lmm)
residDiagnostics(MDAMD231_lmm_gomp)


ObsvsPred(MDAMD231_lmm_gomp)

### Fig. 3 ----

MAS98.06_Her_Fulv <- read.xlsx("data/Fig3/MAS98.06_Her_Fulv.xlsx")
MAS98.06_Her_Fulv$Mouse <- paste(MAS98.06_Her_Fulv$Treatment, " (",MAS98.06_Her_Fulv$Mouse,")", sep = "")
MAS98.06_Her_Fulv <- MAS98.06_Her_Fulv[order(MAS98.06_Her_Fulv$Mouse),]

MAS98.06_Her_Fulv_lmm <- lmmModel(data = MAS98.06_Her_Fulv, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                  trt_control = "Ctrl", drug_a = "Fulv", drug_b = "Her", combination = "Her_Fulv", time_end = 28,
                                  weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 
MAS98.06_Her_Fulv_gomp <- lmmModel(data = MAS98.06_Her_Fulv, grwth_model = "gompertz", start_values = "selfStart",
                                   sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                   trt_control = "Ctrl", drug_a = "Fulv", drug_b = "Her", combination = "Her_Fulv", time_end = 28,
                                   weights = nlme::varIdent(form = ~1|SampleID),
                                   control = nlme::nlmeControl(maxIter = 1000), min_observations = 3) 
ObsvsPred(MAS98.06_Her_Fulv_lmm)


MAS98.06_EGF_Fulv <- read.xlsx("data/Fig3/MAS98.06_EGF_Fulv.xlsx")
MAS98.06_EGF_Fulv$Mouse <- paste(MAS98.06_EGF_Fulv$Treatment, " (",MAS98.06_EGF_Fulv$Mouse,")", sep = "")
MAS98.06_EGF_Fulv <- MAS98.06_EGF_Fulv[order(MAS98.06_EGF_Fulv$Mouse),]

MAS98.06_EGF_Fulv_lmm <- lmmModel(data = MAS98.06_EGF_Fulv, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                  trt_control = "Ctrl", drug_a = "Fulv", drug_b = "EGF", combination = "EGF_Fulv", time_end = 28,
                                  weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

MAS98.06_EGF_Fulv_gomp <- lmmModel(data = MAS98.06_EGF_Fulv,grwth_model = "gompertz", start_values = c(0.00001, 0.001),
                                   sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                   trt_control = "Ctrl", drug_a = "Fulv", drug_b = "EGF", combination = "EGF_Fulv", time_end = 28,
                                   weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 
ranefDiagnostics(MAS98.06_EGF_Fulv_lmm)
ranefDiagnostics(MAS98.06_EGF_Fulv_gomp)
residDiagnostics(MAS98.06_EGF_Fulv_lmm)
residDiagnostics(MAS98.06_EGF_Fulv_gomp)

ObsvsPred(MAS98.06_EGF_Fulv_lmm)
ObsvsPred(MAS98.06_EGF_Fulv_gomp)


MAS98.06_Fulv_Trast <- read.xlsx("data/Fig3/MAS98.06_Fulv_Trast.xlsx")
MAS98.06_Fulv_Trast$Mouse <- paste(MAS98.06_Fulv_Trast$Treatment, " (",MAS98.06_Fulv_Trast$Mouse,")", sep = "")
MAS98.06_Fulv_Trast <- MAS98.06_Fulv_Trast[order(MAS98.06_Fulv_Trast$Mouse),]

MAS98.06_Fulv_Trast_lmm <- lmmModel(data = MAS98.06_Fulv_Trast, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                    trt_control = "Ctrl", drug_a = "Fulv", drug_b = "Trast", combination = "Fulv_Trast", time_end = 28,
                                    weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 
MAS98.06_Fulv_Trast_gomp <- lmmModel(data = MAS98.06_Fulv_Trast, grwth_model = "gompertz", start_values = "selfStart",
                                     sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                     trt_control = "Ctrl", drug_a = "Fulv", drug_b = "Trast", combination = "Fulv_Trast", time_end = 28,
                                     weights = nlme::varIdent(form = ~1|SampleID), 
                                     min_observations = 3) 
ObsvsPred(MAS98.06_Fulv_Trast_lmm)


MAS98.06_Her_Trast <- read.xlsx("data/Fig3/MAS98.06_Her_Trast.xlsx")
MAS98.06_Her_Trast$Sample <- paste(MAS98.06_Her_Trast$Treatment, " (",MAS98.06_Her_Trast$Sample,")", sep = "")
MAS98.06_Her_Trast <- MAS98.06_Her_Trast[order(MAS98.06_Her_Trast$Sample),]

MAS98.06_Her_Trast_lmm <- lmmModel(data = MAS98.06_Her_Trast, sample_id = "Sample", time = "Days", treatment = "Treatment", tumor_vol = "Total_vol",
                                   trt_control = "Ctrl", drug_a = "Her", drug_b = "Trast", combination = "Her_Trast", time_start = 0, time_end = 28,
                                   weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

MAS98.06_Her_Trast_gomp <- lmmModel(data = MAS98.06_Her_Trast, grwth_model = "gompertz", start_values = c(0.00001, 0.001),
                                    sample_id = "Sample", time = "Days", treatment = "Treatment", tumor_vol = "Total_vol",
                                    trt_control = "Ctrl", drug_a = "Her", drug_b = "Trast", combination = "Her_Trast", time_start = 0, time_end = 28,
                                    weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3,
                                    control = nlme::nlmeControl(maxIter = 1000)) 

ranefDiagnostics(MAS98.06_Her_Trast_lmm)
ranefDiagnostics(MAS98.06_Her_Trast_gomp)
residDiagnostics(MAS98.06_Her_Trast_lmm)
residDiagnostics(MAS98.06_Her_Trast_gomp)

ObsvsPred(MAS98.06_Her_Trast_lmm)
ObsvsPred(MAS98.06_Her_Trast_gomp)

### Fig. 4a ----
MAS98.06_Her_Fulv_Trast <- read.xlsx("data/Fig4/MAS98.06_Her_Fulv_Trast.xlsx")
MAS98.06_Her_Fulv_Trast$Sample <- paste(MAS98.06_Her_Fulv_Trast$Treatment, " (",MAS98.06_Her_Fulv_Trast$Sample,")", sep = "")
MAS98.06_Her_Fulv_Trast <- MAS98.06_Her_Fulv_Trast[order(MAS98.06_Her_Fulv_Trast$Sample),]

MAS98.06_Her_Fulv_Trast_lmm <- lmmModel(data = MAS98.06_Her_Fulv_Trast, sample_id = "Sample", time = "Days", treatment = "Treatment", tumor_vol = "Total_vol",
                                        trt_control = "Ctrl", drug_a = "Her", drug_b = "Trast", drug_c = "Fulv", combination = "Her_Fulv_Trast", time_start = 0, time_end = 28,
                                        weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

MAS98.06_Her_Fulv_Trast_lmm <- lmmModel(data = MAS98.06_Her_Fulv_Trast, grwth_model = "gompertz", start_values = "selfStart",
                                        sample_id = "Sample", time = "Days", treatment = "Treatment", tumor_vol = "Total_vol",
                                        trt_control = "Ctrl", drug_a = "Her", drug_b = "Trast", drug_c = "Fulv", combination = "Her_Fulv_Trast", time_start = 0, time_end = 28,
                                        weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 
ObsvsPred(MAS98.06_Her_Fulv_Trast_lmm)

### Fig. 4c ----

triple <- read.csv("data/Fig4/U87MG_Triple_long.csv")

triple$Day <- sqrt(triple$Day)

# Fit model
triple_lmm <- lmmModel(data = triple, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                       trt_control = "None", drug_a = "AZD2014", drug_b = "Tagrisso", drug_c = "DOC",combination = "Triple",
                       weights = nlme::varComb(nlme::varIdent(form = ~1|SampleID), varIdent(form = ~1|Time)))

triple_gomp <- lmmModel(data = triple, grwth_model = "gompertz", start_values = c(0.0001, 0.001),
                        sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                        trt_control = "None", drug_a = "AZD2014", drug_b = "Tagrisso", drug_c = "DOC",combination = "Triple",
                        weights = nlme::varComb(nlme::varIdent(form = ~1|SampleID), varIdent(form = ~1|Time)))
ObsvsPred(triple_lmm)

### Fig. 5a ----

CR1197 <- read.csv("data/Fig5/CR1197_sel_tv.csv")
CR1197$Mice.ID <- paste(CR1197$Treatment, " (",CR1197$Mice.ID,")", sep = "")
CR1197 <- CR1197[order(CR1197$Mice.ID),]

CR1197$DaysPostT0 <- sqrt(CR1197$DaysPostT0)

# Fit model
CR1197_lmm <- lmmModel(data = CR1197, sample_id = "Mice.ID", time = "DaysPostT0", treatment = "Treatment", tumor_vol = "Tumor.Volume",
                       trt_control = "Vehicle", drug_a = "Cetuximab", drug_b = "Palbociclib", combination = "Palbociclib+Cetuximab",
                       weights = nlme::varComb(nlme::varIdent(form = ~1|SampleID), nlme::varPower(form = ~Time)),
                       control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000))


CR1197 <- read.csv("data/Fig5/CR1197_sel_tv.csv")
CR1197$Mice.ID <- paste(CR1197$Treatment, " (",CR1197$Mice.ID,")", sep = "")
CR1197_or <- CR1197[order(CR1197$Mice.ID),]

CR1197_gomp <- lmmModel(data = CR1197_or, grwth_model = "gompertz",
                        sample_id = "Mice.ID", time = "DaysPostT0", treatment = "Treatment", tumor_vol = "Tumor.Volume",
                        trt_control = "Vehicle", drug_a = "Cetuximab", drug_b = "Palbociclib", combination = "Palbociclib+Cetuximab",
                        weights = nlme::varComb(nlme::varIdent(form = ~1|SampleID), nlme::varPower(form = ~Time)),
                        control = nlme::nlmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000))

ranefDiagnostics(CR1197_lmm)
ranefDiagnostics(CR1197_gomp)
residDiagnostics(CR1197_lmm)
residDiagnostics(CR1197_gomp)

ObsvsPred(CR1197_lmm)
ObsvsPred(CR1197_gomp)

### Fig. 5c ----

GABA <- read.csv("data/Fig5/GABA_antiPD1_long.csv")

GABA$Day <- sqrt(GABA$Day-7)

GABA_lmm <- lmmModel(data = GABA, sample_id = "MouseID", time = "Day", treatment = "Treatment", tumor_vol = "TumVol",
                     trt_control = "Ctrl", drug_a = "GABA", drug_b = "Anti-PD1", combination = "GABA+anti-PD1",
                     weights = nlme::varComb(nlme::varPower(form = ~Time), nlme::varIdent(form = ~ 1|Treatment)), 
                     control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000)) 

GABA_or <- read.csv("data/Fig5/GABA_antiPD1_long.csv")
GABA_gomp <- lmmModel(data = GABA_or, grwth_model = "gompertz", start_values = "selfStart",
                      sample_id = "MouseID", time = "Day", treatment = "Treatment", tumor_vol = "TumVol",
                      trt_control = "Ctrl", drug_a = "GABA", drug_b = "Anti-PD1", combination = "GABA+anti-PD1",
                      weights = nlme::varComb(nlme::varPower(form = ~Time), nlme::varIdent(form = ~ 1|Treatment)), 
                      control = nlme::nlmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000)) 


ranefDiagnostics(GABA_lmm)
ranefDiagnostics(GABA_gomp)
residDiagnostics(GABA_lmm)
residDiagnostics(GABA_gomp)

ObsvsPred(GABA_lmm)
ObsvsPred(GABA_gomp)

### Fig. 5e ----

SW837 <- read.xlsx("data/Fig5/SW837_long.xlsx")

# Fit model
SW837_lmm <- lmmModel(data = SW837, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                      trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                      weights = nlme::varIdent(form = ~1|Treatment))

SW837_gomp <- lmmModel(data = SW837, grwth_model = "gompertz", start_values = "selfStart",
                       sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                       trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                       weights = nlme::varIdent(form = ~1|Treatment),
                       control = nlme::nlmeControl(maxIter = 10000))

ObsvsPred(SW837_lmm)

### Supplementary Fig. 9a ----
LS1034 <- read.xlsx("data/Supp_Fig9/LS_1034_long.xlsx")

LS1034_lmm <- lmmModel(data = LS1034, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                       trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                       weights = nlme::varIdent(form = ~1|SampleID), 
                       control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000)) 

LS1034_gomp <- lmmModel(data = LS1034, grwth_model = "gompertz", start_values = c(0.0001, 0.0001),
                        sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                        trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                        weights = nlme::varIdent(form = ~1|SampleID), 
                        control = nlme::nlmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000)) 

ranefDiagnostics(LS1034_lmm)
ranefDiagnostics(LS1034_gomp)
residDiagnostics(LS1034_lmm)
residDiagnostics(LS1034_gomp)

ObsvsPred(LS1034_lmm)
ObsvsPred(LS1034_gomp)

### Supplementary Fig. 9c ----

SNU81 <- read.xlsx("data/Supp_Fig9/SNU81_long.xlsx")

SNU81_lmm <- lmmModel(data = SNU81, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                      trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                      weights = nlme::varIdent(form = ~1|Time),
                      control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000, opt = "optim"))
SNU81_gomp <- lmmModel(data = SNU81, grwth_model = "gompertz", start_values = "selfStart",
                       sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                       trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                       weights = nlme::varIdent(form = ~1|Time),
                       control = nlme::nlmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000))

ObsvsPred(SNU81_lmm)


# Version information about R, the OS and attached or loaded packages. ----

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
