import sys
import matplotlib
matplotlib.use('pdf')
import pandas as pd
#import plotnine as pn
from plotnine import *

df = pd.read_table(sys.stdin)

print df[1:5]
print df.shape

p0 = ggplot(df, aes("ref_allele_ratio")) + \
	geom_histogram(bins=41) + \
	xlab("reference allele ratio")+ \
	ggtitle(sys.argv[1].replace('_',' '))

ggsave(p0, sys.argv[1]+"_ref_allele_ratios.pdf", width=5, height=5)

