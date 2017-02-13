library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

name=args[1]

d <- read.table( file('stdin'), header=TRUE, comment.char="%", check.names=FALSE)
print('')
print(head(d))
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

p1 < -ggplot(d, aes(ref_allele_ratio), check.names=False)+
  geom_histogram(bins=41) +
  theme_bw(base_size=11) +
  labs(title=name, x="ref allele ratio")
p1

ggsave(paste(gsub("/","_", name),"_ref_allele_ratios.pdf", sep=''), p1, width=5, height=5)

