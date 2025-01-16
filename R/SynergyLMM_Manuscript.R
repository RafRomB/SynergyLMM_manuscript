# Code to reproduce the results shown in 'SynergyLMM: A comprehensive statistical framework and interactive web-tool for designing and analyzing in vivo drug combination experiments'

# Libraries

library(SynergyLMM)
library(tidyverse)
library(nlme)
library(cowplot)
library(invivoSyn)
library(openxlsx)
library(ggbreak)

# Figure 2 ----


MDAMD231 <- read.xlsx("data/Fig2/MDA_MD231_Fig5C.xlsx")

MDAMD231$Day <- sqrt(MDAMD231$Day-16) 

# Fit model

MDAMD231_lmm <- lmmModel(data = MDAMD231, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                         trt_control = "Control", drug_a = "AZ628", drug_b = "Gemcitabine", combination = "AZ628_Gemcitabine",
                         weights = varIdent(form = ~ 1|Treatment))

## Fig. 2a ----

trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

df <- MDAMD231_lmm$dt1

df %>% ggplot(aes(x = Time, y = RTV, group = SampleID)) + geom_point(aes(colour = Treatment), size = 3) +
  geom_line(aes(colour = Treatment)) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top") +
  labs(title = "MDA-MB-231 FM")+
  ylab("Relative Luminescence Units (RLU)") + xlab("Days after treatment initiation") +
  scale_x_continuous(breaks = unique(MDAMD231_lmm$dt1$Time), labels = unique(MDAMD231_lmm$dt1$Time)^2)


plot_lmmModel(MDAMD231_lmm, trt_control = "Control", drug_a = "AZ628", drug_b = "Gemcitabine", combination = "AZ628_Gemcitabine") +
  theme(legend.position = "top") +
  labs(title = "MDA-MB-231 FM") +
  scale_x_continuous(breaks = unique(MDAMD231_lmm$dt1$Time), labels = unique(MDAMD231_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18)) +
  xlab("Days after treatment initiation") + ylab("log (RLU)")

## Fig. 2b ----

# Bliss synergy calculation
(bliss <- lmmSynergy(MDAMD231_lmm, robust = T, type = "CR2", min_time = 3, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

# HSA synergy calculation
(hsa <- lmmSynergy(MDAMD231_lmm, robust = T, type = "CR2", min_time = 3, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2


bliss_CI <- bliss$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
bliss_CI$Model <- "SynergyLMM Bliss"
hsa_CI <- hsa$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
hsa_CI$Model <- "SynergyLMM HSA"

# Generate Plot

Narayan <- data.frame(Model = "Narayan et al", Estimate = c(0.08, 0.11), lwr = NA, upr = NA, Time = c(15, 21), pval = NA)

CI_df <- rbind(Narayan, bliss_CI, hsa_CI)
CI_df$Time <- as.factor(CI_df$Time)

CI_df %>% ggplot(aes(x = Time, y = Estimate)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) + cowplot::theme_cowplot() +
  geom_point(aes(fill  = -log10(pval)), size = 5, shape = 23, color = "gray60") +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3, na.value = "firebrick2") +
  ylab("Combination Index") + xlab("Days after treatment initiation") +
  scale_y_continuous(breaks = c(0 , 0.2 , 0.4, 0.6, 0.8, 1)) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 0.8, lty = 3) + 
  facet_wrap(~Model) + 
  theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold", size = 16)) +
  labs(title = "MDA-MB-231 FM\nAZD628 - Gemcitabine")

## Fig. 2c ----


BV173Gluc <- read.delim("data/Fig2/BV_173.csv", sep = "\t")

BV173Gluc$Day <- sqrt(BV173Gluc$Day-7)

# Fit Model
BV173Gluc_lmm <- lmmModel(data = BV173Gluc, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                          trt_control = "Control", drug_a = "Imatinib", drug_b = "Desatinib", combination = "Imatinib_Desatinib")


trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

df <- BV173Gluc_lmm$dt1

df %>% ggplot(aes(x = Time, y = RTV, group = SampleID)) + geom_point(aes(colour = Treatment), size = 3) +
  geom_line(aes(colour = Treatment)) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top") +
  labs(title = "BV-173-Gluc") +
  ylab("Relative Luminescence Units (RLU)") + xlab("Days after treatment initiation") +
  scale_x_continuous(breaks = unique(BV173Gluc_lmm$dt1$Time), labels = unique(BV173Gluc_lmm$dt1$Time)^2)


plot_lmmModel(BV173Gluc_lmm, trt_control = "Control", drug_a = "Imatinib", drug_b = "Desatinib", combination = "Imatinib_Desatinib") +
  theme(legend.position = "top") +
  labs(title = "BV-173-Gluc")  +
  scale_x_continuous(breaks = unique(BV173Gluc_lmm$dt1$Time), labels = unique(BV173Gluc_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18)) +
  xlab("Days after treatment initiation") + ylab("log (RLU)")

## Fig. 2d ----

# Calculate Synergy
(bliss <- lmmSynergy(BV173Gluc_lmm, robust = T, type = "CR2", min_time = 0, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

(hsa <- lmmSynergy(BV173Gluc_lmm, robust = T, type = "CR2", min_time = 0, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2

# Generate Plot

bliss_CI <- bliss$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
bliss_CI$Model <- "SynergyLMM Bliss"
hsa_CI <- hsa$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
hsa_CI$Model <- "SynergyLMM HSA"

Narayan <- data.frame(Model = "Narayan et al", Estimate = c(0.38, 0.24), lwr = NA, upr = NA, Time = c(21-7, 28-7), pval = NA)

CI_df <- rbind(Narayan, bliss_CI, hsa_CI)
CI_df$Time <- as.factor(CI_df$Time)

CI_df %>% ggplot(aes(x = Time, y = Estimate)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) + cowplot::theme_cowplot() +
  geom_point(aes(fill  = -log10(pval)), size = 5, shape = 23, color = "gray60") +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3, na.value = "firebrick2") +
  ylab("Combination Index") + xlab("Days after treatment initiation") +
  #scale_y_continuous(breaks = c(0 , 0.8 , 1, seq(2, 20, 2))) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 0.8, lty = 3) + 
  facet_wrap(~Model) + 
  theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold", size = 16)) +
  labs(title = "BV-173-Gluc\nImatinib - Desatinib") + ggbreak::scale_y_break(breaks = c(7.5,10), scales = 0.25)


## Fig. 2e ----

CHL1FM <- read.xlsx("data/Fig2/CHL1_FM.xlsx")

CHL1FM$Day <- sqrt(CHL1FM$Day-7)

# Fit model
CHL1FM_lmm <- lmmModel(data = CHL1FM, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                       trt_control = "Control", drug_a = "Gemcitabine", drug_b = "CGP_082996", combination = "Gemcitabine_CGP",
                       weights = nlme::varIdent(form = ~1|SampleID),
                       control = lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000))

# Generate plots
trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

df <- CHL1FM_lmm$dt1

df %>% ggplot(aes(x = Time, y = RTV, group = SampleID)) + geom_point(aes(colour = Treatment), size = 3) +
  geom_line(aes(colour = Treatment)) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top") +
  labs(title = "CHL1-FM") +
  ylab("Relative Luminescence Units (RLU)") + xlab("Days after treatment initiation") +
  scale_x_continuous(breaks = unique(CHL1FM_lmm$dt1$Time), labels = unique(CHL1FM_lmm$dt1$Time)^2)


plot_lmmModel(CHL1FM_lmm, trt_control = "Control", drug_a = "Gemcitabine", drug_b = "CGP_082996", combination = "Gemcitabine_CGP") +
  theme(legend.position = "top") +
  labs(title = "CHL1-FM") +
  scale_x_continuous(breaks = unique(CHL1FM_lmm$dt1$Time), labels = unique(CHL1FM_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18)) +
  xlab("Days after treatment initiation") + ylab("log (RLU)")

## Fig. 2f ----

# Calculate Synergy
(bliss <- lmmSynergy(CHL1FM_lmm, robust = T, type = "CR2", min_time = 3, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

(hsa <- lmmSynergy(CHL1FM_lmm, robust = T, type = "CR2", min_time = 3, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2

# Generate plots
bliss_CI <- bliss$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
bliss_CI$Model <- "SynergyLMM Bliss"
hsa_CI <- hsa$Synergy %>% filter(Metric == "CI") %>% select(Model, Estimate, lwr, upr, Time, pval)
hsa_CI$Model <- "SynergyLMM HSA"

Narayan <- data.frame(Model = "Narayan et al", Estimate = c(0.68, 0.62), lwr = NA, upr = NA, Time = c(21-7, 28-7), pval = NA)

CI_df <- rbind(Narayan, bliss_CI, hsa_CI)
CI_df$Time <- as.factor(CI_df$Time)

CI_df %>% ggplot(aes(x = Time, y = Estimate)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) + cowplot::theme_cowplot() +
  geom_point(aes(fill  = -log10(pval)), size = 5, shape = 23, color = "gray60") +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3, na.value = "firebrick2") +
  ylab("Combination Index") + xlab("Days after treatment initiation") +
  scale_y_continuous(breaks = c(0 , 0.8 , 1, seq(2.5, 12.5, 2.5))) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 0.8, lty = 3) + 
  facet_wrap(~Model) + 
  theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold", size = 16)) +
  labs(title = "CHL-1 FM\nCGP-082996 - Gemcitabine")

## Fig. 2g ----

U87MG <- read.xlsx("data/Fig2/U87MG.xlsx")

U87MG$Day <- sqrt(U87MG$Day-12)

# Fit model
U87MG_lmm <- lmmModel(data = U87MG, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                      trt_control = "Control", drug_a = "Docetaxel", drug_b = "GNE-317", combination = "Docetaxel+GNE-317",
                      weights = varPower(form = ~Time))

# Generate plots
trt_col <- c("#3c3c3b", "#d50c52", "#00a49c","#601580")

df <- U87MG_lmm$dt1

df %>% ggplot(aes(x = Time, y = RTV, group = SampleID)) + geom_point(aes(colour = Treatment), size = 3) +
  geom_line(aes(colour = Treatment)) + theme_cowplot() +
  scale_color_manual(values = trt_col) +
  theme(legend.position = "top") +
  labs(title = "U-87-MG-FM") +
  ylab("Relative Luminescence Units (RLU)") + xlab("Days after treatment initiation") +
  scale_x_continuous(breaks = unique(U87MG_lmm$dt1$Time), labels = unique(U87MG_lmm$dt1$Time)^2)

plot_lmmModel(U87MG_lmm, trt_control = "Control", drug_a = "Docetaxel", drug_b = "GNE-317", combination = "Docetaxel+GNE-317") +
  theme(legend.position = "top") +
  labs(title = "U-87-MG-FM") +
  scale_x_continuous(breaks = unique(U87MG_lmm$dt1$Time), labels = unique(U87MG_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18)) +
  xlab("Days after treatment initiation") + ylab("log (RLU)")

## Fig. 2h ----

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

CI_df %>% ggplot(aes(x = Time, y = Estimate)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) + cowplot::theme_cowplot() +
  geom_point(aes(fill  = -log10(pval)), size = 5, shape = 23, color = "gray60") +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3, na.value = "firebrick2") +
  ylab("Combination Index") + xlab("Days after treatment initiation") + 
  scale_y_continuous(breaks = c(0 , 0.5 , 0.8, 1, 1.5, 2)) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 0.8, lty = 3) + 
  facet_wrap(~Model) + 
  theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold", size = 16)) +
  labs(title = "U-87-MG-FM\nDocetaxel - GNE-317")


# Supplementary Figure 1 ----

## Supplementary Fig. 1a ----

ranefDiagnostics(MDAMD231_lmm)$Plots[1]
residDiagnostics(MDAMD231_lmm)$Plots[1]
residDiagnostics(MDAMD231_lmm)$Plots[4]

## Supplementary Fig. 1b ----

ranefDiagnostics(BV173Gluc_lmm)$Plots[1]
residDiagnostics(BV173Gluc_lmm)$Plots[1]
residDiagnostics(BV173Gluc_lmm)$Plots[4]

## Supplementary Fig. 1c ----

ranefDiagnostics(CHL1FM_lmm)$Plots[1]
residDiagnostics(CHL1FM_lmm)$Plots[1]
residDiagnostics(CHL1FM_lmm)$Plots[4]

## Supplementary Fig. 1d ----

ranefDiagnostics(U87MG_lmm)$Plots[1]
residDiagnostics(U87MG_lmm)$Plots[1]
residDiagnostics(U87MG_lmm)$Plots[4]

# Supplementary Figure 2 ----

## Supplementary Fig. 2a ----

ObsvsPred(MDAMD231_lmm, 6,5)

## Supplementary Fig. 2b ----

ObsvsPred(BV173Gluc_lmm, 4,7)

## Supplementary Fig. 2c ----

ObsvsPred(CHL1FM_lmm, 5,5)

## Supplementary Fig. 2d ----

ObsvsPred(U87MG_lmm, 3,4)


# Supplementary Figure 3 ----

## Supplementary Fig. 3a ----

logLikSubjectDisplacements(MDAMD231_lmm, var_name = "Treatment")

## Supplementary Fig. 3b ----

CookDistance(MDAMD231_lmm)

## Supplementary Fig. 3c ----

df <- MDAMD231_lmm$dt1

df$Outliers <- NA

df$Outliers[df$SampleID == "3-1_AZ628"] <- "3-1_AZ628"
df$Outliers[df$SampleID == "7-4_Control"] <- "7-4_Control"
df$Outliers[df$SampleID == "9-3_AZ628"] <- "9-3_AZ628"
df$Outliers[!df$SampleID %in% c("3-1_AZ628","7-4_Control","9-3_AZ628")] <- ""

df %>% ggplot(aes(x = Time, y = logRTV, group = SampleID)) + geom_point(aes(color = Outliers, fill = Outliers, shape = Outliers), size = 2) +
  geom_line(aes(color = Outliers), alpha = 0.5) + scale_color_manual(values = c("gray","#d12856","gray15","#f20747")) +
  scale_fill_manual(values = c("gray","#d12856","gray15","#f20747")) +
  scale_shape_manual(values = c(21:24)) +
  facet_wrap(~Treatment) + theme_cowplot() +
  labs(title = "MDA-MB-231 FM") + ylab("log (Relative Luminiscense Units)") +
  theme(strip.text = element_text(size = 18), legend.position = "top") +
  scale_x_continuous(breaks = unique(df$Time), labels = unique(df$Time)^2) +
  xlab("Days since treatment initiation")


# Figure 3 ----

## Fig. 3a ----

CR1197 <- read.csv("data/Fig3/CR1197_sel_tv.csv")
CR1197$Mice.ID <- paste(CR1197$Treatment, " (",CR1197$Mice.ID,")", sep = "")
CR1197 <- CR1197[order(CR1197$Mice.ID),]

CR1197$DaysPostT0 <- sqrt(CR1197$DaysPostT0)

# Fit model
CR1197_lmm <- lmmModel(data = CR1197, sample_id = "Mice.ID", time = "DaysPostT0", treatment = "Treatment", tumor_vol = "Tumor.Volume",
                       trt_control = "Vehicle", drug_a = "Cetuximab", drug_b = "Palbociclib", combination = "Palbociclib+Cetuximab",
                       weights = nlme::varComb(nlme::varIdent(form = ~1|SampleID), nlme::varPower(form = ~Time)),
                       control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000))


# Generate plot
plot_lmmModel(CR1197_lmm, trt_control = "Vehicle", drug_a = "Cetuximab", drug_b = "Palbociclib", combination = "Palbociclib+Cetuximab") +
  theme(legend.position = "top") +
  labs(title = "CR1197 Cetuximab - Palbociclib") +  ylab("log (Relative Tumor Volume)") +
  scale_x_continuous(breaks = unique(CR1197_lmm$dt1$Time), labels = unique(CR1197_lmm$dt1$Time)^2) + 
  theme(strip.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1))

## Fig. 3b ----

# Calculate Synergy

## SynergyLMM
(bliss <- lmmSynergy(CR1197_lmm, robust = T, type = "CR2", min_time = 2, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

## invivoSyn
CR1197 <- read.csv("data/Fig3/CR1197_sel_tv.csv")

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

AUC_lst <- get_mAUCr(CR1197, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 1000)

# Bliss Synergy with Time

CR1197_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(CR1197$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  CR1197_invivoSyn_bliss <- rbind(CR1197_invivoSyn_bliss, tSyn)
  }

# The simulations can take a while to run. Load the results running the following line of code:
CR1197_invivoSyn_bliss <- read.csv("data/Fig3/CR1197_invivoSyn_bliss_byT.csv")

CR1197_bliss <- bliss$Synergy %>% filter(Metric == "CI")

CR1197_bliss <- CR1197_bliss %>% select(Estimate, lwr, upr, pval, Time)

CR1197_bliss$Model <- "SynergyLMM"

CR1197_invivoSyn_bliss <- CR1197_invivoSyn_bliss %>% filter(Metric == "CI")

CR1197_invivoSyn_bliss <- CR1197_invivoSyn_bliss %>% select(Value, lb, ub, p.val,Time)

CR1197_invivoSyn_bliss$Model <- "invivoSyn"

colnames(CR1197_invivoSyn_bliss) <- colnames(CR1197_bliss)

CR1197_invivoSyn_bliss <- CR1197_invivoSyn_bliss %>% filter(Time > 0)

CR1197_syn <- rbind(CR1197_invivoSyn_bliss, CR1197_bliss)

CR1197_syn$Time <- as.factor(CR1197_syn$Time)

CR1197_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.75)) +
  geom_point(aes(fill = -log10(pval), shape = Model), size = 5, color = "gray65", stroke = 1, position = position_dodge(width = 0.75)) +
  #scale_fill_manual(values = c("firebrick", "darkcyan")) +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
  scale_shape_manual(values = c(22,23)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Combination Index") + 
  theme_cowplot() +
  labs(title = "CR1197 Cetuximab - Palbociclib Bliss Synergy Comparison", subtitle = "invivoSyn - SynergyLMM")

## Fig. 3c ----

GABA <- read.csv("data/Fig3/GABA_antiPD1_long.csv")

GABA$Day <- sqrt(GABA$Day-7)

# Fit model
GABA_lmm <- lmmModel(data = GABA, sample_id = "MouseID", time = "Day", treatment = "Treatment", tumor_vol = "TumVol",
                     trt_control = "Ctrl", drug_a = "GABA", drug_b = "Anti-PD1", combination = "GABA+anti-PD1",
                     weights = nlme::varComb(nlme::varPower(form = ~Time), nlme::varIdent(form = ~ 1|Treatment)), 
                     control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000)) 

# Generate plots
plot_lmmModel(GABA_lmm, trt_control = "Ctrl", drug_a = "GABA", drug_b = "Anti-PD1", combination = "GABA+anti-PD1") +
  theme(legend.position = "top") +
  labs(title = "4T1 GABA + Anti-PD1") +  ylab("log (Relative Tumor Volume)") +
  scale_x_continuous(breaks = unique(GABA_lmm$dt1$Time), labels = unique(GABA_lmm$dt1$Time)^2) +
  theme(strip.text = element_text(size = 18))

## Fig. 3d ----

# Synergy calculation

## SynergyLMM

(hsa <- lmmSynergy(GABA_lmm, robust = T, type = "CR2", min_time = 2.5, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2

## invivoSyn

GABA <- read.csv("data/Fig3/GABA_antiPD1_long.csv")

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

AUC_lst <- get_mAUCr(GABA, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 1000)

# HSA

# Synergy with Time

GABA_invivoSyn_hsa <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(GABA$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "HSA", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  GABA_invivoSyn_hsa <- rbind(GABA_invivoSyn_hsa, tSyn)
  }

# The simulations can take a while to run. Load the results running the following line of code:

GABA_invivoSyn_hsa <- read.csv("data/Fig3/GABA_invivoSyn_hsa_byT.csv")

## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Fig3/4T1_GABA_aPD1_CombPDX.Rdata")
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

GABA_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "4T1 GABA - Anti PD-1 hsa Synergy Comparison", subtitle = "CombPDX - invivoSyn - SynergyLMM") #+ ggbreak::scale_y_break(breaks = c(10,100), scales = 0.25)

## Fig. 3e ----

SW837 <- read.xlsx("data/Fig3/SW837_long.xlsx")

# Fit model
SW837_lmm <- lmmModel(data = SW837, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                      trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                      weights = nlme::varIdent(form = ~1|Treatment))

# Generate plot
plot_lmmModel(SW837_lmm, trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib") +
  theme(legend.position = "top") +
  labs(title = "SW837") + ylab("log (Relative Tumor Volume)") +
  theme(strip.text = element_text(size = 18))

## Fig. 3f ----

# Calculate Synergy

## SynergyLMM

(bliss <- lmmSynergy(SW837_lmm, robust = T, type = "CR2", min_time = 14, method = "Bliss"))

## invivoSyn

SW837 <- read.xlsx("data/Fig3/SW837_long.xlsx")

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

AUC_lst <- get_mAUCr(SW837, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 1000)

# Bliss Synergy with Time

SW837_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(SW837$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SW837_invivoSyn_bliss <- rbind(SW837_invivoSyn_bliss, tSyn)
  }

# The simulations can take a while to run. Load the results running the following line of code:
SW837_invivoSyn_bliss <- read.csv("data/Fig3/SW837_invivoSyn_bliss_byT.csv")

## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:

load("data/Fig3/SW837_CombPDX.Rdata")
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

SW837_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "SW837 Rabusertib - Irinotecan Bliss Synergy Comparison", subtitle = "CombPDX - invivoSyn - SynergyLMM") #+ ggbreak::scale_y_break(breaks = c(10,100), scales = 0.25)

# Supplementary Figure 4 ----

## Supplementary Fig. 4a ----

LS1034 <- read.xlsx("data/Supp_Fig4/LS_1034_long.xlsx")

LS1034_lmm <- lmmModel(data = LS1034, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                       trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                       weights = nlme::varIdent(form = ~1|SampleID), 
                       control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000)) 

plot_lmmModel(LS1034_lmm, trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib") +
  theme(legend.position = "top") +
  labs(title = "LS-1034 Rabusertib - Irinotecan") + ylab("log (Relative Tumor Volume)") +
  theme(strip.text = element_text(size = 18))

## Supplementary Fig. 4b ----

# Synergy Calculation

## SynergyLMM

(bliss <- lmmSynergy(LS1034_lmm, robust = T, type = "CR2", min_time = 10, method = "Bliss"))

## invivoSyn

LS1034 <- read.xlsx("data/Supp_Fig4/LS_1034_long.xlsx")

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

AUC_lst <- get_mAUCr(LS1034, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 1000)


# Bliss Synergy with Time

LS1034_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(LS1034$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  LS1034_invivoSyn_bliss <- rbind(LS1034_invivoSyn_bliss, tSyn)
}

# The simulations can take a while to run. Load the results running the following line of code:

LS1034_invivoSyn_bliss <- read.csv("data/Supp_Fig4//LS1034_invivoSyn_bliss_byT.csv")


## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:

load("data/Supp_Fig4/LS1034_CombPDX.Rdata")
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


LS1034_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "LS1034 Rabusertib - Irinotecan Bliss Synergy Comparison", subtitle = "CombPDX - invivoSyn - SynergyLMM") #+ ggbreak::scale_y_break(breaks = c(10,100), scales = 0.25)

## Supplementary Fig. 4c ----

SNU81 <- read.xlsx("data/Supp_Fig4/SNU81_long.xlsx")

# Fit model
SNU81_lmm <- lmmModel(data = SNU81, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                      trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib",
                      weights = nlme::varIdent(form = ~1|Time),
                      control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, msMaxEval = 1000, opt = "optim"))

# Generate plot
plot_lmmModel(SNU81_lmm, trt_control = "Ctrl", drug_a = "Rabusertib", drug_b = "Irinotecan", combination = "Irinotecan_Rabusertib") +
  theme(legend.position = "top") +
  labs(title = "SNU-81 Rabusertib - Irinotecan") + ylab("log (Relative Tumor Volume)") +
  theme(strip.text = element_text(size = 18))

## Supplementary Fig. 4d ----

# Synergy Calculation

## SynergyLMM

(bliss <- lmmSynergy(SNU81_lmm, robust = T, type = "CR2", min_time = 0, method = "Bliss"))

## invivoSyn

SNU81 <- read.xlsx("data/Supp_Fig4/SNU81_long.xlsx")

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

AUC_lst <- get_mAUCr(SNU81, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 1000)

# Bliss Synergy with Time

SNU81_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(SNU81$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss", ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SNU81_invivoSyn_bliss <- rbind(SNU81_invivoSyn_bliss, tSyn)
}

# The simulations can take a while to run. Load the results running the following line of code:

SNU81_invivoSyn_bliss <- read.csv("data/Supp_Fig4/SNU81_invivoSyn_bliss_byT.csv")

## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:

load("data/Supp_Fig4/SNU81_CombPDX.Rdata")
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

SNU81_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "SNU81 Rabusertib - Irinotecan Bliss Synergy Comparison", subtitle = "CombPDX - invivoSyn - SynergyLMM") #+ ggbreak::scale_y_break(breaks = c(10,100), scales = 0.25)

# Supplementary Figure 5 ----

## Supplementary Fig. 5a ----

ranefDiagnostics(CR1197_lmm)$Plots[1]
residDiagnostics(CR1197_lmm)$Plots[1]
residDiagnostics(CR1197_lmm)$Plots[4]

## Supplementary Fig. 5b ----

ranefDiagnostics(GABA_lmm)$Plots[1]
residDiagnostics(GABA_lmm)$Plots[1]
residDiagnostics(GABA_lmm)$Plots[4]

## Supplementary Fig. 5c ----

ranefDiagnostics(SW837_lmm)$Plots[1]
residDiagnostics(SW837_lmm)$Plots[1]
residDiagnostics(SW837_lmm)$Plots[4]

## Supplementary Fig. 5d ----

ranefDiagnostics(LS1034_lmm)$Plots[1]
residDiagnostics(LS1034_lmm)$Plots[1]
residDiagnostics(LS1034_lmm)$Plots[4]

## Supplementary Fig. 5e ----

ranefDiagnostics(SNU81_lmm)$Plots[1]
residDiagnostics(SNU81_lmm)$Plots[1]
residDiagnostics(SNU81_lmm)$Plots[4]

# Supplementary Figure 6 ----

## Supplementary Fig. 6a ----

ObsvsPred(CR1197_lmm, 6, 7)

## Supplementary Fig. 6b ----

ObsvsPred(GABA_lmm, 6, 7)

## Supplementary Fig. 6c ----

ObsvsPred(SW837_lmm, 5, 5)

## Supplementary Fig. 6d ----

ObsvsPred(LS1034_lmm, 6, 7)

## Supplementary Fig. 6e ----

ObsvsPred(SNU81_lmm, 5, 6)

# Supplementary Figure 7 ----

## Supplementary Fig. 7a ----

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
                   trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination")
  
  
  ## Synergy
  
  lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T)
  lmm1_Bliss <- lmm1_Bliss$Synergy
  lmm1_Bliss$Simulation <- i
  
  simulatedSyn <- rbind(simulatedSyn, lmm1_Bliss)
  
}

# The simulations can take some time to run. The results can be loaded as:

simulatedSyn <- read.csv("data/Supp_Fig7/SimulatedSyn1000.csv", row.names = 1)

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


simSS %>% ggplot(aes(x = Time, y = `median(Estimate)`)) +
  geom_segment(aes(x= Time, y = `median(lwr)`, yend = `median(upr)`), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = -log10(`median(pval)`)), shape = 23, size = 5, color = "gray65") +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  scale_x_continuous(breaks = unique(simSS$Time)) +
  labs(title = "Bliss Sinergy for Simulated Data")

## Supplementary Fig. 7b ----

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
                   trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination")
  
  
  ## Synergy
  
  lmm1_Bliss <- lmmSynergy(lmm1, method = "Bliss", robust = T)
  lmm1_Bliss <- lmm1_Bliss$Synergy
  lmm1_Bliss$Simulation <- i
  
  simulatedAdd <- rbind(simulatedAdd, lmm1_Bliss)
  
}

# The simulations can take some time to run. The results can be loaded as:

simulatedAdd <- read.csv("data/Supp_Fig7/SimulatedAdd1000.csv", row.names = 1)

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

simSS %>% ggplot(aes(x = Time, y = `median(Estimate)`)) +
  geom_segment(aes(x= Time, y = `median(lwr)`, yend = `median(upr)`), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = -log10(`median(pval)`)), shape = 23, size = 5, color = "gray65") +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  scale_x_continuous(breaks = unique(simSS$Time)) +
  labs(title = "Bliss Sinergy (Additivity) for Simulated Data")

## Supplementary Fig. 3c ----

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
plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Simulated Data - Synergistic Effect") +
  theme(strip.text = element_text(size = 18))

## Supplementary Fig. 7d ----

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

AUC_lst <- get_mAUCr(grwth_data_1, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 1000)

# Bliss Synergy with Time

SimData1_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(grwth_data_1$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss",ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SimData1_invivoSyn_bliss <- rbind(SimData1_invivoSyn_bliss, tSyn)
}

# The simulations can take a while to run. Load the results running the following line of code:

SimData1_invivoSyn_bliss <- read.csv("data/Supp_Fig7/SimData1_invivoSyn_bliss_byT.csv", row.names = 1)


## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Supp_Fig7/SimData1_CombPDX.Rdata")
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

SimData1_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Simulated Data Bliss Synergy Comparison", subtitle = "CombPDX - invivoSyn - SynergyLMM") #+ ggbreak::scale_y_break(breaks = c(10,100), scales = 0.25)


## Supplementary Fig. 7e ----

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
plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Simulated Data - Synergistic Effect") +
  theme(strip.text = element_text(size = 18))

## Supplementary Fig. 7f ----

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

AUC_lst <- get_mAUCr(grwth_data_1, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 1000)

# Bliss Synergy with Time

SimData1_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(grwth_data_1$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss",ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SimData1_invivoSyn_bliss <- rbind(SimData1_invivoSyn_bliss, tSyn)
}

# The simulations can take a while to run. Load the results running the following line of code:

SimData1_invivoSyn_bliss <- read.csv("data/Supp_Fig7/SimData1_invivoSyn_bliss_byT_outlier.csv", row.names = 1)


## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Supp_Fig7/SimData1_CombPDX_outlier.Rdata")
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

SimData1_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Simulated Data Bliss Synergy Comparison - Outlier", subtitle = "CombPDX - invivoSyn - SynergyLMM") #+ ggbreak::scale_y_break(breaks = c(10,100), scales = 0.25)

## Supplementary Fig. 7g ----

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
plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Simulated Data - Additive Effect") +
  theme(strip.text = element_text(size = 18))

## Supplementary Fig. 7h ----

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

AUC_lst <- get_mAUCr(grwth_data_1, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 1000)

# Bliss Synergy with Time

SimData1_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(grwth_data_1$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss",ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SimData1_invivoSyn_bliss <- rbind(SimData1_invivoSyn_bliss, tSyn)
}

# The simulations can take a while to run. Load the results running the following line of code:

SimData1_invivoSyn_bliss <- read.csv("data/Supp_Fig7/SimData1_invivoSyn_bliss_byTAdditive.csv", row.names = 1)


## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Supp_Fig7/SimData1_CombPDX_Additive.Rdata")
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

SimData1_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Simulated Data Bliss Synergy Comparison - Additive Effect", subtitle = "CombPDX - invivoSyn - SynergyLMM") #+ ggbreak::scale_y_break(breaks = c(10,100), scales = 0.25)

## Supplementary Fig. 7i ----

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
plot_lmmModel(lmm1, trt_control = "Control", drug_a = "DrugA", drug_b = "DrugB", combination = "Combination") +
  labs(title = "Simulated Data - Additive Effect") +
  theme(strip.text = element_text(size = 18))

## Supplementary Fig. 7j ----

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

AUC_lst <- get_mAUCr(grwth_data_1, ci = 0.95, ci_type = "bca", ref_group = "Group 1", nrep = 1000)

# Bliss Synergy with Time

SimData1_invivoSyn_bliss <- data.frame(Metric = NA, Value = NA, std.err = NA, lb = NA, ub = NA, p.val = NA, Time = NA)

for (t in unique(grwth_data_1$Day)) {
  set.seed(123)
  tSyn <- AUC_synergy(auc_lst = AUC_lst, boot_n = 10000, method = "Bliss",ci_type = "bca", t = t, display = T, save = F)
  tSyn$Time <- t
  SimData1_invivoSyn_bliss <- rbind(SimData1_invivoSyn_bliss, tSyn)
}

# The simulations can take a while to run. Load the results running the following line of code:

SimData1_invivoSyn_bliss <- read.csv("data/Supp_Fig7/SimData1_invivoSyn_bliss_byTAdditive_outlier.csv", row.names = 1)


## CombPDX

# The analysis were done using the web-app: 
# Load the results by running:
load("data/Supp_Fig7/SimData1_CombPDX_Additive_outlier.Rdata")
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

SimData1_syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both"), position = position_dodge(0.5)) +
  geom_point(aes(fill = Model, shape = Model), size = 5, color = "gray65", position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("slateblue1", "firebrick", "darkcyan")) +
  scale_shape_manual(values = c(21,22,23)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after treatment initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Simulated Data Bliss Synergy Comparison - Additive Effect - Outlier", subtitle = "CombPDX - invivoSyn - SynergyLMM") #+ ggbreak::scale_y_break(breaks = c(10,100), scales = 0.25)


# Figure 4 ----

## Fig. 4b ----

MAS98.06 <- read.csv("data/Fig4/MAS98_06.csv")
colnames(MAS98.06) <- c("Cell", "SampleID","Treatment", "TV","Time", "Figure")
MAS98.06 <- getRTV(MAS98.06, time_start = 0)

MAS98.06$Treatment <- factor(MAS98.06$Treatment, levels = c("Ctrl", "Her", "EGF", "Fulv", "Trast", "Her_Fulv", "EGF_Fulv","Fulv_Trast","Her_Trast"))

hline <- data.frame(yintercept = 0)

MAS98.06 %>% ggplot(aes(x = Time, y = logRTV, color = Treatment)) +
  geom_line(aes(group = SampleID), alpha = 0.33) + geom_point(aes(group = SampleID, shape = Treatment, fill = Treatment)) +
  ylab("Log (RTV)") + 
  xlab("Time since start of treatment") + 
  scale_x_continuous(breaks = c(0,7,14,21,28)) + 
  cowplot::theme_cowplot() + facet_wrap(~Treatment) +
  labs(title = "MAS98.06") +
  geom_hline(data = hline, aes(yintercept = yintercept), linetype = "dashed") +
  scale_color_manual(values = c("#3c3c3b", "#00a49c", "#fc8d62","#d50c52","#8c510a","#601580","#e78ac3", "#01665e","#386cb0")) +
  scale_fill_manual(values = c("#3c3c3b", "#00a49c", "#fc8d62","#d50c52","#8c510a","#601580","#e78ac3", "#01665e","#386cb0")) +
  scale_shape_manual(values = c(21:25,21:24)) +
  theme(strip.text = element_text(size = 18), legend.position = "top")

## Fig. 4c ----

MAS98.06_Her_Fulv <- read.xlsx("data/Fig4/MAS98.06_Her_Fulv.xlsx")
MAS98.06_Her_Fulv$Mouse <- paste(MAS98.06_Her_Fulv$Treatment, " (",MAS98.06_Her_Fulv$Mouse,")", sep = "")
MAS98.06_Her_Fulv <- MAS98.06_Her_Fulv[order(MAS98.06_Her_Fulv$Mouse),]

MAS98.06_Her_Fulv_lmm <- lmmModel(data = MAS98.06_Her_Fulv, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                  trt_control = "Ctrl", drug_a = "Fulv", drug_b = "Her", combination = "Her_Fulv", time_end = 28,
                                  weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

(Her_Fulv <- lmmSynergy(MAS98.06_Her_Fulv_lmm, robust = T, type = "CR2", min_time = 7, method = "HSA"))
Her_Fulv <- Her_Fulv$Synergy
Her_Fulv$Drug <- "Her_Fulv"

MAS98.06_EGF_Fulv <- read.xlsx("data/Fig4/MAS98.06_EGF_Fulv.xlsx")
MAS98.06_EGF_Fulv$Mouse <- paste(MAS98.06_EGF_Fulv$Treatment, " (",MAS98.06_EGF_Fulv$Mouse,")", sep = "")
MAS98.06_EGF_Fulv <- MAS98.06_EGF_Fulv[order(MAS98.06_EGF_Fulv$Mouse),]


MAS98.06_EGF_Fulv_lmm <- lmmModel(data = MAS98.06_EGF_Fulv, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                  trt_control = "Ctrl", drug_a = "Fulv", drug_b = "EGF", combination = "EGF_Fulv", time_end = 28,
                                  weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

(EGF_Fulv <- lmmSynergy(MAS98.06_EGF_Fulv_lmm, robust = T, type = "CR2", min_time = 7, method = "HSA"))
EGF_Fulv <- EGF_Fulv$Synergy

EGF_Fulv$Drug <- "EGF_Fulv"

MAS98.06_Fulv_Trast <- read.xlsx("data/Fig4/MAS98.06_Fulv_Trast.xlsx")
MAS98.06_Fulv_Trast$Mouse <- paste(MAS98.06_Fulv_Trast$Treatment, " (",MAS98.06_Fulv_Trast$Mouse,")", sep = "")
MAS98.06_Fulv_Trast <- MAS98.06_Fulv_Trast[order(MAS98.06_Fulv_Trast$Mouse),]

MAS98.06_Fulv_Trast_lmm <- lmmModel(data = MAS98.06_Fulv_Trast, sample_id = "Mouse", time = "Day", treatment = "Treatment", tumor_vol = "TV",
                                    trt_control = "Ctrl", drug_a = "Fulv", drug_b = "Trast", combination = "Fulv_Trast", time_end = 28,
                                    weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

(Fulv_Trast <- lmmSynergy(MAS98.06_Fulv_Trast_lmm, robust = T, type = "CR2", min_time = 8, method = "HSA"))
Fulv_Trast <- Fulv_Trast$Synergy
Fulv_Trast$Drug <- "Fulv_Trast"

MAS98.06_Her_Trast <- read.xlsx("data/Fig4/MAS98.06_Her_Trast.xlsx")
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
  geom_point(aes(fill = -log10(pval)), size = 5, color = "gray65", shape = 23) +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
  geom_hline(yintercept = 1, lty = "dashed") +
  xlab("Drugs") +
  ylab("Combination Index") + 
  theme_cowplot() +
  labs(title = "Combination Index")


SS <- HSA_SS %>% ggplot(aes(x = Drug, y = Estimate, group = Drug)) +
  geom_segment(aes(x= Drug, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) +
  geom_point(aes(fill = -log10(pval)), size = 5, color = "gray65", shape = 23) +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Drugs") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Synergy Score")

cowplot::plot_grid(CI, SS)

# Supplementary Figure 8 ----

## Supplementary Fig. 8a ----

ranefDiagnostics(MAS98.06_Her_Fulv_lmm)$Plots[1]
residDiagnostics(MAS98.06_Her_Fulv_lmm)$Plots[1]
residDiagnostics(MAS98.06_Her_Fulv_lmm)$Plots[4]

## Supplementary Fig. 8b ----

ranefDiagnostics(MAS98.06_EGF_Fulv_lmm)$Plots[1]
residDiagnostics(MAS98.06_EGF_Fulv_lmm)$Plots[1]
residDiagnostics(MAS98.06_EGF_Fulv_lmm)$Plots[4]

## Supplementary Fig. 8c ----

ranefDiagnostics(MAS98.06_Fulv_Trast_lmm)$Plots[1]
residDiagnostics(MAS98.06_Fulv_Trast_lmm)$Plots[1]
residDiagnostics(MAS98.06_Fulv_Trast_lmm)$Plots[4]

## Supplementary Fig. 8d ----

ranefDiagnostics(MAS98.06_Her_Trast_lmm)$Plots[1]
residDiagnostics(MAS98.06_Her_Trast_lmm)$Plots[1]
residDiagnostics(MAS98.06_Her_Trast_lmm)$Plots[4]


# Supplementary Figure 9 ----

## Supplementary Fig. 9a ----
ObsvsPred(MAS98.06_Her_Fulv_lmm, 5, 5)

## Supplementary Fig. 9b ----
ObsvsPred(MAS98.06_EGF_Fulv_lmm, 4, 5)

## Supplementary Fig. 9c ----
ObsvsPred(MAS98.06_Fulv_Trast_lmm, 4, 5)

## Supplementary Fig. 9d ----
ObsvsPred(MAS98.06_Her_Trast_lmm, 4, 6)

# Figure 5 ----

## Fig. 5a ----

MAS98.06_Her_Fulv_Trast <- read.xlsx("data/Fig5/MAS98.06_Her_Fulv_Trast.xlsx")
MAS98.06_Her_Fulv_Trast$Sample <- paste(MAS98.06_Her_Fulv_Trast$Treatment, " (",MAS98.06_Her_Fulv_Trast$Sample,")", sep = "")
MAS98.06_Her_Fulv_Trast <- MAS98.06_Her_Fulv_Trast[order(MAS98.06_Her_Fulv_Trast$Sample),]

# Fit Model
MAS98.06_Her_Fulv_Trast_lmm <- lmmModel(data = MAS98.06_Her_Fulv_Trast, sample_id = "Sample", time = "Days", treatment = "Treatment", tumor_vol = "Total_vol",
                                        trt_control = "Ctrl", drug_a = "Her", drug_b = "Trast", drug_c = "Fulv", combination = "Her_Fulv_Trast", time_start = 0, time_end = 28,
                                        weights = nlme::varIdent(form = ~1|SampleID), min_observations = 3) 

trt_col <- c("#3c3c3b", "#00a49c", "#fc8d62","#d50c52","#8c510a","#601580","#e78ac3", "#01665e","#386cb0")

plot_lmmModel(MAS98.06_Her_Fulv_Trast_lmm,   trt_control = "Ctrl", drug_a = "Her", drug_b = "Trast", drug_c = "Fulv", combination = "Her_Fulv_Trast") +
  theme(legend.position = "top") +
  labs(title = "MAS98.06 Heregulin - Fulvestrant - Trastuzumab") + scale_x_continuous(breaks = c(0,7,14,21,28)) +
  ylab("log (Relative Tumor Volume)") +
  scale_color_manual(values = c(trt_col[c(1,2,4,5)], "red3")) +
  theme(strip.text = element_text(size = 18))

## Fig. 5b ----

(bliss <- lmmSynergy(MAS98.06_Her_Fulv_Trast_lmm, robust = T, type = "CR2", min_time = 8, method = "Bliss"))

(hsa <- lmmSynergy(MAS98.06_Her_Fulv_Trast_lmm, robust = T, type = "CR2", min_time = 8, method = "HSA"))

set.seed(123)
(ra <- lmmSynergy(MAS98.06_Her_Fulv_Trast_lmm, robust = T, type = "CR2", min_time = 8, method = "RA", ra_nsim = 1000))

syn <- rbind(bliss$Synergy, hsa$Synergy, ra$Synergy)

syn <- syn %>% filter(Metric == "SS" & Time %in% c(8,12, 16, 20, 24, 28))

syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) +
  geom_point(aes(fill = -log10(pval)), size = 5, color = "gray65", shape = 23) +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after Treatment Initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Synergy Score") + facet_wrap(~Model) +
  scale_x_continuous(breaks = c(8,12, 16, 20, 24, 28)) +
  theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold", size = 18))

## Fig. 5c ----

triple <- read.csv("data/Fig5/U87MG_Triple_long.csv")

triple$Day <- sqrt(triple$Day)

# Fit model
triple_lmm <- lmmModel(data = triple, sample_id = "MOUSE", time = "Day", treatment = "Treatment", tumor_vol = "RLU",
                       trt_control = "None", drug_a = "AZD2014", drug_b = "Tagrisso", drug_c = "DOC",combination = "Triple",
                       weights = nlme::varComb(nlme::varIdent(form = ~1|SampleID), varIdent(form = ~1|Time)))

plot_lmmModel(triple_lmm, trt_control = "None", drug_a = "AZD2014", drug_b = "Tagrisso", drug_c = "DOC",combination = "Triple") +
  theme(legend.position = "top") +
  labs(title = "U87-MG-FM Tagrisso-AZD2014-DOC") +
  scale_x_continuous(breaks = unique(triple_lmm$dt1$Time), labels = unique(triple_lmm$dt1$Time)^2) + ylab("log (RLU)") +
  theme(strip.text = element_text(size = 18))

## Fig. 5d ----

(bliss <- lmmSynergy(triple_lmm, robust = T, type = "CR2", min_time = 3, method = "Bliss"))
bliss$Synergy$Time <- bliss$Synergy$Time^2

(hsa <- lmmSynergy(triple_lmm, robust = T, type = "CR2", min_time = 3, method = "HSA"))
hsa$Synergy$Time <- hsa$Synergy$Time^2

set.seed(123)
(ra <- lmmSynergy(triple_lmm, robust = T, type = "CR2", min_time = 3, method = "RA", ra_nsim = 1000))
ra$Synergy$Time <- ra$Synergy$Time^2

syn <- rbind(bliss$Synergy, hsa$Synergy, ra$Synergy)

syn <- syn %>% filter(Metric == "SS")

syn %>% ggplot(aes(x = Time, y = Estimate, group = Model)) +
  geom_segment(aes(x= Time, y = lwr, yend = upr), color = "gray60", lwd = 1, 
               arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both")) +
  geom_point(aes(fill = -log10(pval)), size = 5, color = "gray65", shape = 23) +
  scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab("Days after Treatment Initiation") +
  ylab("Synergy Score") + 
  theme_cowplot() +
  labs(title = "Synergy Score") + facet_wrap(~Model) +
  scale_x_continuous(breaks = unique(syn$Time)) +
  theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold", size = 18))

# Supplementary Figure 10 ----

## Supplementary Fig. 10a ----

ranefDiagnostics(MAS98.06_Her_Fulv_Trast_lmm)$Plots[1]
residDiagnostics(MAS98.06_Her_Fulv_Trast_lmm)$Plots[1]
residDiagnostics(MAS98.06_Her_Fulv_Trast_lmm)$Plots[4]

## Supplementary Fig. 10b ----

ObsvsPred(MAS98.06_Her_Fulv_Trast_lmm, 5, 5)

## Supplementary Fig. 10c ----

ranefDiagnostics(triple_lmm)$Plots[1]
residDiagnostics(triple_lmm)$Plots[1]
residDiagnostics(triple_lmm)$Plots[4]

## Supplementary Fig. 10d ----

ObsvsPred(triple_lmm, 6, 6)

# Figure 6 ----

### Post-Hoc Power ----
# Note: the Post-Hoc power calculation can take some time to run

# MDA-MB-231
set.seed(123)
PostHocPwr(MDAMD231_lmm, method = "Bliss")
# [1] 0.52

set.seed(123)
PostHocPwr(MDAMD231_lmm, method = "HSA")
# [1] 0.916

MDAMD231 <- lmmModel_estimates(MDAMD231_lmm)

MDAMD231$PwrBliss <- 0.52
MDAMD231$pwrHSA <- 0.916

MDAMD231 <- read.csv("data/Fig6/MDAMD231_model_estimates.csv")
MDAMD231$Cell <- "MDA-MB-231"
colnames(MDAMD231) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# BV-173-Gluc

set.seed(123)
PostHocPwr(BV173Gluc_lmm, method = "Bliss")
# [1] 0.781

set.seed(123)
PostHocPwr(BV173Gluc_lmm, method = "HSA")
# [1] 0.894

BV173Gluc <- lmmModel_estimates(BV173Gluc_lmm)

BV173Gluc$PwrBliss <- 0.781
BV173Gluc$PwrHSA <- 0.894

BV173Gluc <- read.csv("data/Fig6/BV173Gluc_model_estimates.csv")
BV173Gluc$Cell <- "BV-173-Gluc"
colnames(BV173Gluc) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

# CHL1 FM

set.seed(123)
PostHocPwr(CHL1FM_lmm, method = "Bliss")
# [1] 0.66

set.seed(123)
PostHocPwr(CHL1FM_lmm, method = "HSA")
# [1] 0.384

CHL1FM <- lmmModel_estimates(CHL1FM_lmm)

CHL1FM$PwrBliss <- 0.66
CHL1FM$PwrHSA <- 0.384

CHL1FM <- read.csv("data/Fig6/CHL1FM_model_estimates.csv")
CHL1FM$Cell <- "CHL1-FM"
colnames(CHL1FM) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

U87MG <- read.csv("data/Fig6/U87MG_model_estimates.csv")
U87MG$Cell <- "U-87-MG-FM"
colnames(U87MG) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

CR1197 <- read.csv("data/Fig6/CR1197_model_estimates.csv")
CR1197$Cell <- "CR1197"
colnames(CR1197) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

LS1034 <- read.csv("data/Fig6/LS1034_model_estimates.csv")
LS1034$Cell <- "LS-1034"
colnames(LS1034) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

SW837 <- read.csv("data/Fig6/SW837_model_estimates.csv")
SW837$Cell <- "SW837"
colnames(SW837) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

SNU81 <- read.csv("data/Fig6/SNU81_model_estimates.csv")
SNU81$Cell <- "SNU-81"
colnames(SNU81) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

GABA <- read.csv("data/Fig6/GABA_model_estimates.csv")
GABA$Cell <- "4T1"
colnames(GABA) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

Her_Fulv <- read.csv("data/Fig6/MAS98.06_Her_Fulv_model_estimates.csv")
Her_Fulv$Cell <- "MAS98.06 (Her+Fulv)"
colnames(Her_Fulv) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

EGF_Fulv <- read.csv("data/Fig6/MAS98.06_EGF_Fulv_model_estimates.csv")
EGF_Fulv$Cell <- "MAS98.06 (EGF+Fulv)"
colnames(EGF_Fulv) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

Trast_Fulv <- read.csv("data/Fig6/MAS98.06_Fulv_Trast_model_estimates.csv")
Trast_Fulv$Cell <- "MAS98.06 (Trast+Fulv)"
colnames(Trast_Fulv) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

Her_Trast <- read.csv("data/Fig6/MAS98.06_Her_Trast_model_estimates.csv")
Her_Trast$Cell <- "MAS98.06 (Her+Trast)"
colnames(Her_Trast) <- c("X","Ctrl", "DrugA", "DrugB", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

Her_Trast_Fulv <- read.csv("data/Fig6/MAS98.06_Her_Fulv_Trast_model_estimates.csv")
Her_Trast_Fulv$Cell <- "MAS98.06 (Her+Fulv+Trast)"
colnames(Her_Trast_Fulv) <- c("X","Ctrl", "DrugA", "DrugB", "DrugC", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

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

U87MG_Triple <- read.csv("data/Fig6/U87_Triple_model_estimates.csv")
U87MG_Triple$Cell <- "U-87-MG-FM (AZD+Osi+Doc)"
colnames(U87MG_Triple) <-c("X","Ctrl", "DrugA", "DrugB", "DrugC", "Combination", "sd_ranef", "sd_resid", "PwrBliss", "PwrHSA", "Cell")

unit_norm_l1 <- function(x) {
  x / sum(abs(x))  # Divide by the L1 norm
  }


Pwr_df <- rbind(BV173Gluc, CHL1FM, CR1197, EGF_Fulv, GABA, Her_Fulv, Her_Trast, LS1034,MDAMD231, SNU81, SW837, Trast_Fulv,U87MG)

rownames(Pwr_df) <- Pwr_df$Cell
Pwr_bliss_hsa <- Pwr_df[,c(8:10)]
Pwr_df <- Pwr_df[,-c(1,8:10)]

Pwr_df_norm <- t(apply(Pwr_df, 1, unit_norm_l1))

Pwr_df_norm <- Pwr_df_norm[,c(4:6)]

Pwr_df_triple <- rbind(Her_Trast_Fulv, U87MG_Triple)
rownames(Pwr_df_triple) <- Pwr_df_triple$Cell
Pwr_bliss_hsa <- rbind(Pwr_bliss_hsa, Pwr_df_triple[,c(9:11)])
Pwr_df_triple <- Pwr_df_triple[,-c(1,9:11)]

Pwr_df_triple_norm <- t(apply(Pwr_df_triple, 1, unit_norm_l1))

Pwr_df_triple_norm <- Pwr_df_triple_norm[,5:7]


Pwr_df_norm <- rbind(Pwr_df_norm, Pwr_df_triple_norm)
Pwr_df_norm <- as.data.frame(Pwr_df_norm)

Pwr_df_norm$PwrBliss <- Pwr_bliss_hsa$PwrBliss
Pwr_df_norm$PwrHSA <- Pwr_bliss_hsa$PwrHSA
Pwr_df_norm$Cell <- Pwr_bliss_hsa$Cell


## Fig. 6a ----

Pwr_df_norm %>% ggplot(aes(x = sd_ranef, y = sd_resid, group = Cell)) +
  geom_point(aes(size = PwrBliss, colour = Combination)) + 
  ggrepel::geom_text_repel(aes(label = Cell),size = 3, box.padding = 0.5,segment.color = "gray25") + theme_cowplot() +
  scale_color_gradient2(name = "Normalized\nCombination Growth Rate",low = "purple3",mid = "#f7f7f7", high = "#e66101") +
  scale_size_area(name = "Statistical Power", max_size = 10, limits = c(0,1)) +
  xlab("Normalized SD of Random Effects") + ylab("Normalized SD of Residuals") +
  labs(title = "Power for Bliss Synergy")

## Fig. 6b ----

df <- MDAMD231_lmm$dt1

segment_data <- data.frame(x = rep(0,4), 
                           xend = df %>% dplyr::group_by(Treatment) %>% dplyr::summarise(Max = max(Time)) %>% dplyr::select(Max),
                           y = rep(0, 4), 
                           yend = nlme::fixef(MDAMD231_lmm))

segment_data$yend <- segment_data$Max*segment_data$yend
colnames(segment_data) <- c("x", "xend", "y", "yend")

segment_data$Treatment <- factor(x = c("Control", "AZ628", "Gemcitabine", "AZ628_Gemcitabine"))

df %>% ggplot(aes(x = Time, y = logRTV)) + geom_point(aes(colour = Treatment, group = SampleID), shape = 19, alpha = 0.33) +
  geom_line(aes(colour = Treatment, group = SampleID), alpha = 0.25) +
  scale_color_manual(values = c("#3c3c3b", "#d50c52", "#00a49c","#601580")) +
  geom_segment(data = segment_data, 
               aes(x = x, xend = xend, y = y, yend = yend), 
               lwd = 1.5, alpha = 1, colour = c("#3c3c3b", "#d50c52", "#00a49c","#601580"))  + theme_cowplot() +
  labs(title = "MDA-MB-231-FM AZD628-Gemcitabine") +
  scale_x_continuous(breaks = unique(df$Time), labels = unique(df$Time)^2) +
  xlab("Days after treatment initiation") + ylab("log (Relative Luminescence Units)")


## Fig. 6c and 6d ----

APrioriPwr(npg = 7, time = unique(MDAMD231_lmm$dt1$Time), grwrControl = 1.54, grwrA = 1.48,
           grwrB = 0.79, grwrComb = 0.46, sd_ranef = 0.25, sgma = 0.7, sd_eval = seq(0.1*0.25, 2*0.25, 0.05*0.25),
           sgma_eval = seq(0.1*0.7, 2*0.7, 0.05*0.7), grwrComb_eval = seq(0.01*0.46, 3*0.46, 0.05*0.46))


## Fig. 6e ----

PwrSampleSize(npg = 4:21, time = unique(MDAMD231_lmm$dt1$Time), grwrControl = 1.54, grwrA = 1.48,
              grwrB = 0.79, grwrComb = 0.46, sd_ranef = 0.25, sgma = 0.7)


## Fig. 6f ----

df <- LS1034_lmm$dt1

segment_data <- data.frame(x = rep(0,4), 
                           xend = df %>% dplyr::group_by(Treatment) %>% dplyr::summarise(Max = max(Time)) %>% dplyr::select(Max),
                           y = rep(0, 4), 
                           yend = nlme::fixef(LS1034_lmm))

segment_data$yend <- segment_data$Max*segment_data$yend
colnames(segment_data) <- c("x", "xend", "y", "yend")

segment_data$Treatment <- factor(x = c("Ctrl", "Rabusertib", "Irinotecan", "Irinotecan_Rabusertib"))

df %>% ggplot(aes(x = Time, y = logRTV)) + geom_point(aes(colour = Treatment, group = SampleID), shape = 19, alpha = 0.33) +
  geom_line(aes(colour = Treatment, group = SampleID), alpha = 0.25) +
  scale_color_manual(values = c("#3c3c3b", "#d50c52", "#00a49c","#601580")) +
  geom_segment(data = segment_data, 
               aes(x = x, xend = xend, y = y, yend = yend), 
               lwd = 1.5, alpha = 1, colour = c("#3c3c3b", "#d50c52", "#00a49c","#601580"))  + theme_cowplot() +
  labs(title = "LS-1031 Rabusertib - Irinotecan") + xlab("Days since start of treatment") + ylab("log (Relative Tumor Volume)") +
  scale_x_continuous(breaks = unique(df$Time))

## Fig. 6g ----

PwrTime(npg = 10, time = list(seq(0,9,3), seq(0,12,3), seq(0,15,3),seq(0,18,3), seq(0,21,3), seq(0,24,3), seq(0,27,3), seq(0,30,3),seq(0,33,3)), 
        grwrControl = 0.097, grwrA = 0.103, grwrB = -0.006, grwrComb = -0.037, sd_ranef = 0.028, sgma = 0.402)

## Fig. 6h ----

PwrTime(npg = 10, time = list(seq(0,18,1), seq(0,18,2), seq(0,18,3), seq(0,18,6), seq(0,18,9), seq(0,18,18)), 
        grwrControl = 0.097, grwrA = 0.103, grwrB = -0.006, grwrComb = -0.037, sd_ranef = 0.028, sgma = 0.402,
        type = "freq")

# Version information about R, the OS and attached or loaded packages. ----

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
