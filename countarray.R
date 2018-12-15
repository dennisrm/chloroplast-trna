#count plotter ? 
library(ggplot2)

countarray <- read.delim("C:/Users/dmull/bme232/tRNA/countarray-GAA-Phe-postflank.txt", header=FALSE)
names(countarray)=c("a","c","g","t","gap")

ggplot(countarray)+
  geom_point(aes(row.names(countarray),a),col="red")+
  geom_point(aes(row.names(countarray),c),col="blue")+
  geom_point(aes(row.names(countarray),g),col="purple")+
  geom_point(aes(row.names(countarray),t),col="darkblue")+
  geom_point(aes(row.names(countarray),gap))+
  ylim(1,120)
#########################


temp.identity <- read.delim("C:/Users/dmull/bme232/tRNA/temp-identity.txt")

temp.identity = temp.identity[temp.identity$seqid>0.9 | (temp.identity$seqid>0.7 & (temp.identity$preid >0.9 | temp.identity$postid>0.9)),]

temp.identity = temp.identity[
  (temp.identity$seqid>0.90| temp.identity$preid >0.90 | temp.identity$postid>0.90) &
    (temp.identity$seqid>0.6 & temp.identity$preid >0.6 & temp.identity$postid>0.6),]
temp.identity$rank = rank(temp.identity$seqid,ties.method = "first")

subset = temp.identity[temp.identity$chromosome!='chloroplast' & temp.identity$chromosome != 'bacteria',]
mean(subset$postid/subset$seqid)

mean(temp.identity$preid)

sum(temp.identity$chromosome != "chloroplast" & temp.identity$chromosome != "bacteria")

######
ggplot(temp.identity)+
  geom_point(aes(postid,preid,color=genus))+
  scale_color_manual(values=c("gray15","gray20","purple","gray25","gray30","gray35","gray40",
                              "gray45","gray50","gray55","gray60","gray65","gray70","blue","gray75"))+
  labs(title="Chloroplast Flanking Region Identity",x="Post-flank Identity",y="Pre-flank Identity")
  
#####

######3
ggplot(temp.identity)+
  geom_point(aes(rank,preid,color=chromosome),shape=8)+
  geom_point(aes(rank,seqid,color=chromosome))+
  geom_point(aes(rank,postid,color=chromosome),shape=1)+
  #scale_color_manual(values=c("red","green3","blue","darkblue","gray","black"))+
  scale_color_manual(values=c("red","green3","gray20","gray25","gray30","gray35","gray40",
                              "gray45","gray50","gray55","gray60","gray65","gray70","gray75"))+
  labs(title="ALL tRNA Gene Identity: tRNA Gene, Pre- and Post-Flanks",y="Identity",x="Ordered by tRNA Gene Sequence Identity")
  theme(legend.position = "none")
  

####################
ggplot(temp.identity)+geom_point(aes(rank,identity,color=chromosome,shape=genus))+
  labs(title="CAA tRNA Gene - No Flanking Sequences")+
  scale_color_manual(values=c("red","green3","blue","darkblue","gray","black"))+
  theme(legend.position = "top")+
  theme(legend.title=element_blank())
  
##############  
temp.identity <- read.delim("C:/Users/dmull/bme232/tRNA/temp-identity.txt")
  
ggplot(temp.identity)+geom_histogram(aes(preid,fill=chromosome),color="black",binwidth = 0.01)+
  scale_fill_manual(values=c("red","green3","gray75","gray70","gray65","gray60","gray55",
                              "gray50","gray45","gray40","gray35","gray30","gray25","gray20"))+
  labs(x="Sequence Identity",y="Number of Sequences",
       title="Number of Sequences vs Identity Value - Preflank region - By Chromosome")+
  theme(legend.position="right",legend.title = element_blank())

ggplot(temp.identity)+geom_histogram(aes(seqid,fill=amino),color="black",binwidth = 0.01)+
  labs(x="Sequence Identity",y="Number of Sequences",
       title="Number of Sequences vs Identity Level - Gene region - By Amino Acid")+
  theme(legend.position="right",legend.title = element_blank())

temp.identity <- read.delim("C:/Users/dmull/bme232/tRNA/temp-identity.txt")
temp.identity$spechro=paste(temp.identity$species,temp.identity$chromosome)
temp.identity$spechro = sub('[0-9]',"",temp.identity$spechro)
temp.identity$spechro = sub('[0-9]',"",temp.identity$spechro)

ggplot(temp.identity)+geom_histogram(aes(preid,fill=spechro),color="black",binwidth = 0.01)+
  scale_fill_manual(values=c("green3","gray60","orange","red2","green2","gray40","red3",
                             "green4","purple2","darkgreen","purple4","orange"))+
  labs(x="Sequence Identity",y="Number of Sequences",
       title="Number of Sequences vs Identity Value - Postflank region - By Species")+
  theme(legend.position="right",legend.title = element_blank())

 