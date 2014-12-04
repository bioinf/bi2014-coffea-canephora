setwd('/media/drozdovapb/big/Studies/bioinformaticsinstitute/bioinfo-lab-14/bi2014-coffea-canephora/Gene_clusters/')

library(ape)
library(phangorn)
?tree
plot


###
install.packages('ade4')
library(ade4)

tr <- read.tree(text="(>Arabidopsis_thaliana:5,(>Coffea_canephora:1,>Solanum_lycopersicon:2):3);")
plot(tr)

#CAFE results

## IDs of nodes:(>AT<0>,(>Cc<2>,>So<4>)<3>)<1>
# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', 
#and 'Branch-specific P-values' = (node ID, node ID): (0,3) (2,4) 
# Output format for 'Branch cutting P-values' and 'Likelihood Ratio Test': (0, 1, 2, 3, 4)
#Average Expansion:    (-0.449018,-0.311996)	(-0.341409,-0.150791)
#Expansion :	(4343,407)	(1245,3467)
#Remain :	(15160,32766)	(28726,33023)
#Decrease:	(28504,14834)	(18036,11517)

?text
png('tree.png', width=480, height=320)
plot(tr)
text(x=1, y=1.5, labels="(4343,407)\n(15160,32766)\n(28504,14834)", cex=.9)
text(x=4, y=2.5, labels="(1245,3467)\n(28726,33023)\n(18036,11517)", cex=.9)
text(x=0.5, y=2.8, labels="Expansion\nRemain\nDecrease", cex=.9)
dev.off()

plot(tr)
text(x=1, y=1.5, labels="(4343,407)\n(15160,32766)\n(28504,14834)", cex=.9)
text(x=4, y=2.5, labels="(1245,3467)\n(28726,33023)\n(18036,11517)", cex=.9)
text(x=0.5, y=2.8, labels="Expansion\nRemain\nDecrease", cex=.9)

png('tree.png', width=480, height=320)
plot(tr)
text(x=4.1, y=2.9, labels="4343 / 15160 / 28504", cex=.8)
text(x=4, y=1.85, labels="3467 / 33023 / 11517", cex=.8)
text(x=1.8, y=2.6, labels="407 / 32766 / 14834", cex=.8)
text(x=1.8, y=1.1, labels="1245 / 28726 / 18036", cex=.8)
dev.off()
