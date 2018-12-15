#homology

postflank <- read.delim("C:/Users/dmull/bme232/tRNA/homology-cp-postflank.txt")

preflank <- read.delim("C:/Users/dmull/bme232/tRNA/homology-cp-preflank.txt")

sequence <- read.delim("C:/Users/dmull/bme232/tRNA/homology-cp-sequence.txt")


sequence$position = sequence$position+1
postflank$position = postflank$position+2

allseq = rbind(preflank,sequence,postflank)

ggplot(allseq)+
  geom_point(aes(position,identity,color=chromosome),alpha="0.5")#+
  geom_smooth(aes(position,identity,group=chromosome),
              method="loess",span=0.02,se=FALSE,
              color="black")+
  labs(title="Chlorplast tRNA Sequence Identity to Consensus - Gene and Flanks",
       x="Distance Along: Preflank-Sequence-Postflank",
       y="Proportion of Sequences that Match Consensus")+
  geom_vline(xintercept = 1,alpha=0.3)+geom_vline(xintercept=2,alpha=0.3)#+
    theme(legend.position = "none")
