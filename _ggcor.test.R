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

dt.test[, FC := 2^logFC]
nologdt <- dt.test[, c("gene", "FC", "P-value", "gene group", "ontology weight"), ]

# basic
ggvolcano.test(data = nologdt, p.cutoff = 0.00001, l2fc.cutoff = 1.5,
               label.cutoff = 2, y.col.sign = TRUE, l2.transform = TRUE)

# Manually tell which point should be labeled with force.label.
ggvolcano.test(data = dt.test, p.cutoff = 0.00001, l2fc.cutoff = 1.5,
               label.cutoff = 2, x.col.sign = FALSE, l2.transform = FALSE)
# Create a function that supports Free Y-axis values
ggvolcano.free(data = dt.test, label.cutoff = 1, p.cutoff = 0.00001,
               x.cutoff = 1, title.x.cutoff = "Delta methylation cut-off",
               x.col.sign = FALSE, force.label = c(
                 "Sox2", "DDIT4", "CA9", "ANKRD37", "NDRG1", "ANGPTL4", "ENO2",
                 "EGLN3", "INSIG2", "MYCBP", "BNIP3", "RPA2", "HK2", "MPRS17",
                 "KDM3A", "ADM", "MKi67", "RFC4", "BRCA1", "PCNA", "PDK1",
                 "KCTD11", "CDKN3", "ATR", "BNIP3L", "P4HA1"))

#Test panel correlation
dt.test2 <- data.table(
  "gene" = c("Sox2", "DDIT4", "CA9", "ANKRD37", "NDRG1", "ANGPTL4", "ENO2",
             "EGLN3", "INSIG2", "MYCBP", "BNIP3", "RPA2", "HK2", "MPRS17",
             "KDM3A", "ADM", "MKi67", "RFC4", "BRCA1", "PCNA", "PDK1", "KCTD11",
             "CDKN3", "ATR", "BNIP3L", "P4HA1"),
  "correlation" = rnorm(n = 26, mean = 0, sd = 0.3),
  "P-value" = c( 1.476961e-19, 3.234006e-19, 8.833711e-26, 1.730797e-24,
                 2.465876e-21, 0.003812423, 0.004239498, 0.004392182,
                 0.005220615, 0.005357136, 0.005703521, 0.005756199,
                 0.007229292, 0.011714971, 0.015301721, 0.016189800,
                 0.017896844, 0.017976486, 0.019790513, 0.019868352,
                 0.020064957, 0.020361514, 0.021650503, 0.026751092,
                 0.033636934, 0.038799781),
  "gene group" = rep(c("A", "B", "C", "D"), times = c(6, 6, 6, 8)),
  "ontology weight" = sample.int(n = 26))

ggpanel.corr(data = dt.test2, label.cutoff = 0.2, jitter.height = 0.5, )

