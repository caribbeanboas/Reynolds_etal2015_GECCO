data<-read.csv("Book1.csv",header=T)
names<-c("1","3","6","9","12")

par(mfrow=c(2,1))
boxplot(data$Ho1,data$Ho3,data$Ho6,data$Ho9,data$Ho12,xlab="Number of Loci", names=names, ylab="Observed Heterozygosity")
boxplot(data$allele1,data$allele3,data$allele6,data$allele9,data$allele12,names=names,xlab="Number of Loci", ylab="Number of Alleles")
