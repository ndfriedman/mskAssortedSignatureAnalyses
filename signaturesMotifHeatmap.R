#written by Noah Friedman

install.packages("reshape");
library(reshape)
library(scales)

signatures <- readRDS("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/pensona/dmp_sigs/diff_test/signatures.rds")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

meltedSigs = melt(signatures)
#meltedSigs[c("value")] <- lapply(meltedSigs[c("value")], function(x) log(x))
#meltedSigs[c("value")] <- lapply(meltedSigs[c("value")], 
#                                 function(x) ifelse(x > -4.564348, x, -4.564348))

write.table(meltedSigs, file='/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/sigMotifsMelted.tsv', quote=FALSE, sep='\t', col.names = NA)

sigsDf <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/sigMotifsMelted_adj.tsv', sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

plt <- ggplot(sigsDf, 
       aes(reorder(X2, trinucOrderingCol), reorder(X1, sigOrderingCol), fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours=c("white","#ffffe0","purple", "black"),
  #scale_fill_gradientn(colours=c("white","#f2f2f2","black"),
                       values  = rescale(c(0, .03333, .1, maxVal)))+
  #scale_fill_gradient(low = "blue", high = "red")+
  theme(axis.text.x = element_text(angle = 90, size=3, vjust=.2),
        axis.text.y = element_text(size=5))+
  ggtitle('Trinucleotide Percentages Across Signatures')

plt <- plt + labs(x = "Trinucleotide Change") 
plt <- plt + labs(y = "Signature")

ggsave('~/Desktop/noahDogTe.pdf', plt)










sigsDfOther <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/sigMotifsMelted_adjSignif.tsv', sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

sigsDfOther$value

plt <- ggplot(sigsDfOther, 
              aes(reorder(X2, trinucOrderingCol), reorder(X1, sigOrderingCol), fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours=c("white","#ffffe0","purple", "black"),
                       #scale_fill_gradientn(colours=c("white","#f2f2f2","black"),
                       values  = rescale(c(0, .03333, .1, maxVal)))+
  #scale_fill_gradient(low = "blue", high = "red")+
  theme(axis.text.x = element_text(angle = 90, size=3, vjust=.2),
        axis.text.y = element_text(size=5))+
  ggtitle('Trinucleotide Percentages: between 0.005 and 0.05')

plt <- plt + labs(x = "Trinucleotide Change") 
plt <- plt + labs(y = "Signature")

ggsave('~/Desktop/otherVersion.pdf', plt)

max(sigsDfOther$value)

sigsDfOther <- read.table('/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/sigMotifsMelted_adjSignif.tsv', sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
plt <- ggplot(sigsDfOther, 
              aes(reorder(X2, trinucOrderingCol), reorder(X1, sigOrderingCol), fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours=c("white","#DCDCDC", "black"),
                       values  = rescale(c(0, .01, .5)))+
  #scale_fill_gradient(low = "blue", high = "red")+
  theme(axis.text.x = element_text(angle = 90, size=3, vjust=.2),
        axis.text.y = element_text(size=5))+
  ggtitle('Trinucleotide Percentages Across Signatures in KRAS Hotspots')

plt<- plt + theme(axis.title.x = element_text(vjust=-2))
plt <- plt + labs(x = "Trinucleotide Change") 
plt <- plt + labs(y = "Signature")

ggsave('~/Desktop/noahDogTe.pdf', plt)



