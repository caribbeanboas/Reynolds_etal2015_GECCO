setwd("C:/Users/Graham/Dropbox/Manuscripts/VI Boas/MANUSCRIPT/New MS/Simulations/MetaSim")
getwd()


data<-read.table("data_boxplots.txt",header=T)
names<-c("Ne1000","Ne500","Ne250","Ne1000","Ne500","Ne250")


boxplot(data$Ho1000_72,data$Ho500_72,data$Ho250_72,data$He1000_72,data$He500_72,data$He250_72, data$Ho250, data$He250,ylim=c(0,1), names=names, xlab="Effective Population Size")
points(1,.35, col="black", pch=18, cex=1.5)
points(2,.35, col="black", pch=18, cex=1.5)
points(3,.35, col="black", pch=18, cex=1.5)
points(4,.40, col="black", pch=18, cex=1.5)
points(5,.40, col="black", pch=18, cex=1.5)
points(6,.40, col="black", pch=18, cex=1.5)
abline(v=3.5, lty=3)
legend("topright", "Cayo Diablo", pch=18, cex=1.2)


names3<-c("Ne=72","Ne=8","Ne=1000")
boxplot(data$Na1000_72,data$Na1000_8,data$Na250, names=names3, ylab="# Alleles", xlab="Effective Population Size")
points(1,2.78, col="black", pch=18, cex=1.5)
points(2,2.78, col="black", pch=18, cex=1.5)
points(3,2.78, col="black", pch=18, cex=1.5)
abline(v=1.5, lty=3)
abline(v=2.5, lty=3)
legend("topleft", "St. Thomas", pch=18, cex=1.2)



names<-c("Ho","He","Ho","He","Ho","He")
boxplot(data$Ho1000_72,data$He1000_72, data$Ho1000_8,data$He1000_8,data$Ho250,data$He250,ylim=c(0,1), names=names)
points(1,.35, col="black", pch=18, cex=1.5)
points(2,.40, col="black", pch=18, cex=1.5)
points(3,.35, col="black", pch=18, cex=1.5)
points(4,.40, col="black", pch=18, cex=1.5)
points(5,.35, col="black", pch=18, cex=1.5)
points(6,.40, col="black", pch=18, cex=1.5)
abline(v=2.5, lty=3)
abline(v=4.5, lty=3)
legend("topright", "St. Thomas", pch=18, cex=1.2)
