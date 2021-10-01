# Tutorial code chunk for BiocompR's ggvolcano.test function.
# Author: Verena Bitto

library(BiocompR)
library(data.table)

dt.test <- data.table(
  "gene" = c("Sox2", "DDIT4", "CA9", "ANKRD37", "NDRG1", "ANGPTL4", "ENO2",
             "EGLN3", "INSIG2", "MYCBP", "BNIP3", "RPA2", "HK2", "MPRS17",
             "KDM3A", "ADM", "MKi67", "RFC4", "BRCA1", "PCNA", "PDK1", "KCTD11",
             "CDKN3", "ATR", "BNIP3L", "P4HA1"),
  "logFC" = c(1.2692153707, -2.2290638339, 2.4717095237, 2.8707390695,
              2.2085696781, 0.3848021, -0.5088270, -0.2909698, -0.2061311,
              0.6963465, 0.3052798, 0.2936876, 0.7998389, -0.2318962, 0.2822826,
              -0.3196105, -0.3173267, 0.3468562, -0.3316569, 0.5381214,
              0.2010866, -0.6567128, 0.1185155, -0.3362876, -0.3014294,
              -0.3184143),
  "P-value" = c( 1.476961e-19, 3.234006e-19, 8.833711e-26, 1.730797e-24,
                 2.465876e-21, 0.003812423, 0.004239498, 0.004392182,
                 0.005220615, 0.005357136, 0.005703521, 0.005756199,
                 0.007229292, 0.011714971, 0.015301721, 0.016189800,
                 0.017896844, 0.017976486, 0.019790513, 0.019868352,
                 0.020064957, 0.020361514, 0.021650503, 0.026751092,
                 0.033636934, 0.038799781),
  "gene group" = rep(c("A", "B", "C", "D"), times = c(6, 6, 6, 8)),
  "ontology weight" = sample.int(n = 26))

# basic
ggvolcano.test(data = dt.test[,c(1:5)], p.cutoff = 0.00001, l2fc.cutoff = 1.5,
               label.cutoff = 2, y.col.sign = TRUE)

# # basic with label
# ggvolcano.test(data[,c(1:3)], p.cutoff = 0.00001, label = c("CA9", "Sox2"))
#
#
# # coloring by groups
# data$groups <- c("base")
# data$groups[2:10] <- c("special")
# ggvolcano.test(data, p.cutoff = 0.0001)
#
# # coloring by multiple groups
# colnames(data)[1] <- "gene"
# cutoff <- data$P.Value<=0.0001 & abs(data$logFC) >= 1
# data$groups <- ifelse(data$P.Value<=0.0001, ifelse(abs(data$logFC) >= 2, "logFC >= |2|", "logFC < |2|"), "not significant")
#
# ggvolcano.test(data, p.cutoff = 0.0001)
# ggvolcano.test(data, p.cutoff = 0.0001, label = c("logFC >= |2|"), group.label = TRUE)
# ggvolcano.test(data, p.cutoff = 0.0001, label = c("Sox2"), group.label = 0)
#
#
#
# ## negative checks
# data_fail <- data.table(c("Sox2"), c(1.2692153707), c(1.2692153707))
# ggvolcano.test(data_fail)
# data_fail <- data.table(c("Sox2"),
#                         c("Sox2"),
#                         c(0.05))
# ggvolcano.test(data_fail)
