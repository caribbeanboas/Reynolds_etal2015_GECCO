setwd("C:/Users/Graham/Dropbox/Manuscripts/VI Boas/MANUSCRIPT/New MS/Simulations/Easypop")
getwd()

data1<-read.table("100_1600_001.equ",header=T)
data2<-read.table("100_1600_002.equ",header=T)
data3<-read.table("100_1600_003.equ",header=T)
data4<-read.table("100_1600_004.equ",header=T)
data5<-read.table("100_1600_005.equ",header=T)
data6<-read.table("100_1600_006.equ",header=T)
data7<-read.table("100_1600_007.equ",header=T)
data8<-read.table("100_1600_008.equ",header=T)
data9<-read.table("100_1600_009.equ",header=T)
data10<-read.table("100_1600_010.equ",header=T)

data11<-read.table("50_1600_001.equ",header=T)
data12<-read.table("50_1600_002.equ",header=T)
data13<-read.table("50_1600_003.equ",header=T)
data14<-read.table("50_1600_004.equ",header=T)
data15<-read.table("50_1600_005.equ",header=T)
data16<-read.table("50_1600_006.equ",header=T)
data17<-read.table("50_1600_007.equ",header=T)
data18<-read.table("50_1600_008.equ",header=T)
data19<-read.table("50_1600_009.equ",header=T)
data20<-read.table("50_1600_010.equ",header=T)

data21<-read.table("25_1600_001.equ",header=T)
data22<-read.table("25_1600_002.equ",header=T)
data23<-read.table("25_1600_003.equ",header=T)
data24<-read.table("25_1600_004.equ",header=T)
data25<-read.table("25_1600_005.equ",header=T)
data26<-read.table("25_1600_006.equ",header=T)
data27<-read.table("25_1600_007.equ",header=T)
data28<-read.table("25_1600_008.equ",header=T)
data29<-read.table("25_1600_009.equ",header=T)
data30<-read.table("25_1600_010.equ",header=T)

plot(data21$Ho,type="l",ylab="Ho", xlab="# Generations", ylim=c(0,1))
lines(data21$Ho, col="gray80")
lines(data22$Ho, col="gray80")
lines(data23$Ho, col="gray80")
lines(data24$Ho, col="gray80")
lines(data25$Ho, col="gray80")
lines(data26$Ho, col="gray80")
lines(data27$Ho, col="gray80")
lines(data28$Ho, col="gray80")
lines(data29$Ho, col="gray80")
lines(data30$Ho, col="gray80")

lines(data11$Ho, col="gray50")
lines(data12$Ho, col="gray50")
lines(data13$Ho, col="gray50")
lines(data14$Ho, col="gray50")
lines(data15$Ho, col="gray50")
lines(data16$Ho, col="gray50")
lines(data17$Ho, col="gray50")
lines(data18$Ho, col="gray50")
lines(data19$Ho, col="gray50")
lines(data20$Ho, col="gray50")

lines(data1$Ho)
lines(data2$Ho)
lines(data3$Ho)
lines(data4$Ho)
lines(data5$Ho)
lines(data6$Ho)
lines(data7$Ho)
lines(data8$Ho)
lines(data9$Ho)
lines(data10$Ho)

legend("topright", c("Ne=25","Ne=50","Ne=100"), fill=c("gray80","gray50","black"))
abline(h=0.52, lty=2)




plot(data21$allel,type="l", ylab="# Alleles", xlab="# Generations", ylim=c(1,20))
lines(data21$allel, col="gray80")
lines(data22$allel, col="gray80")
lines(data23$allel, col="gray80")
lines(data24$allel, col="gray80")
lines(data25$allel, col="gray80")
lines(data26$allel, col="gray80")
lines(data27$allel, col="gray80")
lines(data28$allel, col="gray80")
lines(data29$allel, col="gray80")
lines(data30$allel, col="gray80")

lines(data1$allel)
lines(data2$allel)
lines(data3$allel)
lines(data4$allel)
lines(data5$allel)
lines(data6$allel)
lines(data7$allel)
lines(data8$allel)
lines(data9$allel)
lines(data10$allel)

lines(data11$allel, col="gray50")
lines(data12$allel, col="gray50")
lines(data13$allel, col="gray50")
lines(data14$allel, col="gray50")
lines(data15$allel, col="gray50")
lines(data16$allel, col="gray50")
lines(data17$allel, col="gray50")
lines(data18$allel, col="gray50")
lines(data19$allel, col="gray50")
lines(data20$allel, col="gray50")

legend("topright", c("Ne=25","Ne=50","Ne=100"), fill=c("gray80","gray50","black"))
abline(h=2.33, lty=2)



##RESAMPLE
data<-read.table("genalex.txt",header=F)

for (i in 1:100){
resamp=data[sample(nrow(data),size=3,replace=TRUE),]
write.table(resamp,file="resamp.txt",append=TRUE)
}


#Mutation=0.01
data1<-read.table("high_mut_001.equ",header=T)
data2<-read.table("high_mut_002.equ",header=T)
data3<-read.table("high_mut_003.equ",header=T)
data4<-read.table("high_mut_004.equ",header=T)
data5<-read.table("high_mut_005.equ",header=T)
data6<-read.table("high_mut_006.equ",header=T)
data7<-read.table("high_mut_007.equ",header=T)
data8<-read.table("high_mut_008.equ",header=T)
data9<-read.table("high_mut_009.equ",header=T)
data10<-read.table("high_mut_010.equ",header=T)

data11<-read.table("50_mut_001.equ",header=T)
data12<-read.table("50_mut_002.equ",header=T)
data13<-read.table("50_mut_003.equ",header=T)
data14<-read.table("50_mut_004.equ",header=T)
data15<-read.table("50_mut_005.equ",header=T)
data16<-read.table("50_mut_006.equ",header=T)
data17<-read.table("50_mut_007.equ",header=T)
data18<-read.table("50_mut_008.equ",header=T)
data19<-read.table("50_mut_009.equ",header=T)
data20<-read.table("50_mut_010.equ",header=T)

data21<-read.table("25_mut_001.equ",header=T)
data22<-read.table("25_mut_002.equ",header=T)
data23<-read.table("25_mut_003.equ",header=T)
data24<-read.table("25_mut_004.equ",header=T)
data25<-read.table("25_mut_005.equ",header=T)
data26<-read.table("25_mut_006.equ",header=T)
data27<-read.table("25_mut_007.equ",header=T)
data28<-read.table("25_mut_008.equ",header=T)
data29<-read.table("25_mut_009.equ",header=T)
data30<-read.table("25_mut_010.equ",header=T)

plot(data21$Ho,type="l",ylab="Ho", xlab="# Generations", ylim=c(0,1))
lines(data21$Ho, col="gray80")
lines(data22$Ho, col="gray80")
lines(data23$Ho, col="gray80")
lines(data24$Ho, col="gray80")
lines(data25$Ho, col="gray80")
lines(data26$Ho, col="gray80")
lines(data27$Ho, col="gray80")
lines(data28$Ho, col="gray80")
lines(data29$Ho, col="gray80")
lines(data30$Ho, col="gray80")

lines(data11$Ho, col="gray50")
lines(data12$Ho, col="gray50")
lines(data13$Ho, col="gray50")
lines(data14$Ho, col="gray50")
lines(data15$Ho, col="gray50")
lines(data16$Ho, col="gray50")
lines(data17$Ho, col="gray50")
lines(data18$Ho, col="gray50")
lines(data19$Ho, col="gray50")
lines(data20$Ho, col="gray50")

lines(data1$Ho)
lines(data2$Ho)
lines(data3$Ho)
lines(data4$Ho)
lines(data5$Ho)
lines(data6$Ho)
lines(data7$Ho)
lines(data8$Ho)
lines(data9$Ho)
lines(data10$Ho)

legend("topright", c("Ne=25","Ne=50","Ne=100"), fill=c("gray80","gray50","black"))
abline(h=0.52,  lty=2)




plot(data21$allel,type="l", ylab="# Alleles", xlab="# Generations", ylim=c(1,20))
lines(data21$allel, col="gray80")
lines(data22$allel, col="gray80")
lines(data23$allel, col="gray80")
lines(data24$allel, col="gray80")
lines(data25$allel, col="gray80")
lines(data26$allel, col="gray80")
lines(data27$allel, col="gray80")
lines(data28$allel, col="gray80")
lines(data29$allel, col="gray80")
lines(data30$allel, col="gray80")

lines(data1$allel)
lines(data2$allel)
lines(data3$allel)
lines(data4$allel)
lines(data5$allel)
lines(data6$allel)
lines(data7$allel)
lines(data8$allel)
lines(data9$allel)
lines(data10$allel)

lines(data11$allel, col="gray50")
lines(data12$allel, col="gray50")
lines(data13$allel, col="gray50")
lines(data14$allel, col="gray50")
lines(data15$allel, col="gray50")
lines(data16$allel, col="gray50")
lines(data17$allel, col="gray50")
lines(data18$allel, col="gray50")
lines(data19$allel, col="gray50")
lines(data20$allel, col="gray50")

legend("topright", c("Ne=25","Ne=50","Ne=100"), fill=c("gray80","gray50","black"))
abline(h=2.33, lty=2)




data<-read.table("summarystats_resamp_highmu.txt",header=T)
names<-c("Ne100","Ne50","Ne25","Ne100","Ne50","Ne25")


boxplot(data$Ho100,data$Ho50,data$Ho25,data$He100,data$He50,data$He25,ylim=c(0,1), names=names, xlab="Effective Population Size")
points(1,.54, col="black", pch=18, cex=1.5)
points(2,.54, col="black", pch=18, cex=1.5)
points(3,.54, col="black", pch=18, cex=1.5)
points(4,.53, col="black", pch=18, cex=1.5)
points(5,.53, col="black", pch=18, cex=1.5)
points(6,.53, col="black", pch=18, cex=1.5)
abline(v=3.5, lty=3)
legend("topright", "Cayo Diablo", pch=18, cex=1.2)


names3<-c("Ne100","Ne50","Ne25")
boxplot(data$Na100,data$Na50,data$Na25, names=names3, ylab="# Alleles", xlab="Effective Population Size")
points(1,2.33, col="black", pch=18, cex=1.5)
points(2,2.33, col="black", pch=18, cex=1.5)
points(3,2.33, col="black", pch=18, cex=1.5)
legend("topright", c("Cayo Diablo", "Simulated Mean"), col=c("black","gray"),pch=18, cex=1.2)







