import sys
import matplotlib
matplotlib.use('pdf')
import pandas as pd
#import plotnine as pn
from plotnine import *

df = pd.read_table(sys.stdin)

#print df[1:5]
#print df.shape

if len(sys.argv) > 3:
	colmn = sys.argv[3]
	label = ' '.join(sys.argv[3].split('_'))
else: 
	colmn = "ref_allele_ratio"
	label = "reference allele ratio"

p0 = ggplot(df, aes(colmn)) + \
	geom_histogram(bins=41) + \
	xlab(label) + \
	ggtitle(sys.argv[1].replace('_',' '))

ggsave(p0, sys.argv[1] + "_" + colmn + "s." + sys.argv[2] + ".pdf", width=5, height=5)

