# Tutorial code chunk for BiocompR's ggcoverage function.
# Author: Niklas Beumer

plot <- ggcoverage(cell_numbers, round.unit = 1, decreasing.order = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(breaks = c("Subset", "Remaining"), values = c("red", "green"), name = "Status", labels = c("Removed", "Retained")) +
  xlab("Sample") +
  ylab("Cell number") +
  theme_classic() +
  theme(axis.text.y = element_text(colour= "black", size = 14), 
        axis.text.x = element_text(colour= "black", angle = 45, hjust = 1, size = 14),
        axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
png("cell_numbers_by_sample.png", width = 500, height = 300)
print(plot)
dev.off()