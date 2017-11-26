import sys
import matplotlib
matplotlib.use('pdf')
import pandas as pd
#import plotnine as pn
from plotnine import *

df = pd.read_table(sys.stdin)

#print df[1:5]
#print df.shape

if len(sys.argv) > 2:
	colmn = sys.argv[2]
else: colmn = "ref_allele_ratio"

p0 = ggplot(df, aes(colmn)) + \
	geom_histogram(bins=41) + \
	xlab("reference allele ratio")+ \
	ggtitle(sys.argv[1].replace('_',' '))

ggsave(p0, sys.argv[1]+"_"+colmn+"s.pdf", width=5, height=5)

