setwd("D:\\Desktop\\temporary")
library(bio3d)
library(nortest)
pdbf <- read.pdb('4eiycut.pdb')
pdbm <- read.pdb('4ug2cut.pdb')
pdbt <- read.pdb("5g53cut.pdb")
pdbf <- trim.pdb(pdbf, atom.select(pdbf, chain="A"))
pdbm <- trim.pdb(pdbm, atom.select(pdbm, chain="A"))
pdbt <- trim.pdb(pdbt, atom.select(pdbt, chain="A"))
seqf <- pdbseq(pdbf);locf <- names(seqf)
ca.indsf <- atom.select(pdbf,"calpha")
bff<-pdbf$atom[ca.indsf$atom,"b"]
seqm <- pdbseq(pdbm);locm <- names(seqm)
ca.indsm <- atom.select(pdbm,"calpha")
bfm<-pdbm$atom[ca.indsm$atom,"b"]
seqt <- pdbseq(pdbt);loct <- names(seqt)
ca.indst <- atom.select(pdbt,"calpha")
bft<-pdbt$atom[ca.indst$atom,"b"]

#Q-Q图
library(ggpubr)
ggqqplot(bff,title = 'Q-Q plot')
ggqqplot(bfm,title = 'Q-Q plot')
ggqqplot(bft,title = 'Q-Q plot')
###Lilliefor test和S-W检验
options(digits = 3)
a1<-scale(bff)[,1]
a2<-scale(bfm)[,1]
a3<-scale(bft)[,1]
aa<-list(a1,a2,a3)
for (i in 1:3) {
  print(lillie.test(aa[[i]])) 
}
shapiro.test(bff)
shapiro.test(bfm)
shapiro.test(bft)

#斯皮尔曼秩相关系数
fo<-function(){
  pp<-c()
  for (i in 7:50) {
    n6h<-nma.pdb(pdbf, ff="anm", mass=F, cutoff=i)
    fluct0 <- fluct.nma(n6h,mode.inds=seq(7,length(n6h[["L"]])))
    p<-cor(bff,fluct0,method="spearman")
    p<-c(i,p)
    pp<-rbind(pp,p)
  }
  print(pp)
}
ppf<-fo()
mo<-function(){
  pp<-c()
  for (i in 7:50) {
    n6h<-nma.pdb(pdbm, ff="anm", mass=F, cutoff=i)
    fluct0 <- fluct.nma(n6h,mode.inds=seq(7,length(n6h[["L"]])))
    p<-cor(bfm,fluct0,method="spearman")
    p<-c(i,p)
    pp<-rbind(pp,p)
  }
  print(pp)
}
ppm<-mo()
to<-function(){
  pp<-c()
  for (i in 7:50) {
    n6h<-nma.pdb(pdbt, ff="anm",mass = F,cutoff=i)
    fluct0 <- fluct.nma(n6h,mode.inds=seq(7,length(n6h[["L"]])))
    p<-cor(bft,fluct0,method="spearman")
    p<-c(i,p)
    pp<-rbind(pp,p)
  }
  print(pp)
}
ppt<-to()
#截断半径与相关性分布图
pnum <- data.frame(dis=ppf[,1],ppf=ppf[,2],ppm=ppm[,2],ppt=ppt[,2])
tiff("趋势1.tif",width = 5.5,height = 4,units="cm",res=800,compression = "lzw")
par(ps=7.5,family="serif",mai=c(0.25,0.25,0.02,0.02),tck=0.015,mgp=c(0.55,0,0),lwd=0.5)
plot(pnum$dis,pnum$ppf, type = "l", pch = 15, lty = 1,yaxt="n" ,xaxt="n",
     ylim = c(0.4,0.85),col = "black", xlab = "", ylab = "")/
  mtext('截断半径(Å)',side = 1,line = 0.24,cex = 1.12)/
  mtext('斯皮尔曼相关系数',side = 2.2,line = 0.5,cex=1.12)/
  axis(1, at=c(7,10,20,30,40,50), 
       labels=c(7,10,20,30,40,50), las=0,lwd=0.5,padj=-1)/
  axis(2,at=c(0.4,0.5,0.6,0.7,0.8), 
       labels=c(0.4,0.5,0.6,0.7,0.8),las=0,lwd=0.5)/
  lines(pnum$dis,pnum$ppm, type = "l", pch = 16, lty = 2, col = "red")/
  lines(pnum$dis,pnum$ppt, type = "l", pch = 18, lty = 4, col = "purple")/
  legend("topleft", inset =-0.03, c("4EIY","4UG2","5G53"), lty = c(1,2,3),x.intersp = 0, 
         y.intersp =c(0,0.25,0.30) ,col = c("black","red","purple"),bty="n")
dev.off()

####计算B-factor
anmf<-nma.pdb(pdbf, ff="anm", mass=F,cutoff=30)
fluct0f <- fluct.nma(anmf,mode.inds=seq(7,length(anmf[["L"]])))#所有模式下的波动
fluctf<-(8/3)*(pi^2)*5.09*(fluct0f)
cor(bff,fluctf,method="spearman")
oned<-data.frame(V1=1:length(bff),V2=bff,V3=scale(fluctf),V4=fluctf)
#write.table(oned,"4eiy_b.txt")
anmm<-nma.pdb(pdbm, ff="anm", mass=F,cutoff=26)
fluct0m <- fluct.nma(anmm,mode.inds=seq(7,length(anmm[["L"]])))#所有模式下的波动
fluctm<-(8/3)*(pi^2)*7.25*(fluct0m)
cor(bfm,fluctm,method="spearman")
twod<-data.frame(V1=1:length(bfm),V2=bfm,V3=scale(fluctm),V4=fluctm)
#write.table(twod,"4ug2_b.txt")
anmt<-nma.pdb(pdbt, ff="anm", mass=F,cutoff=28)
fluct0t <- fluct.nma(anmt,mode.inds=seq(7,length(anmt[["L"]])))#所有模式下的波动
fluctt<-(8/3)*(pi^2)*11.44*(fluct0t)
cor(bft,fluctt,method="spearman")
threed<-data.frame(V1=1:length(bft),V2=bft,V3=scale(fluctt),V4=fluctt)
#write.table(threed,"5g53_b.txt")

#做图B-factor
###4EIY温度因子图
data1<-read.table("4eiy_b.txt",header = F)
v4<-append(data1$V4, c(rep(NA,10)), after = 208);v2<-append(data1$V2, c(rep(NA,10)), after = 208)
tiff("Bfactor1.tif",width = 5.5,height = 4,units="cm",res=800,compression = "lzw")
par(ps=7.5,family="serif",mai=c(0.25,0.25,0.02,0.02),tck=0.015,mgp=c(0.55,0,0),lwd=0.5)
plot(v4,type = "l",xaxt="n",yaxt="n",xlab = "",ylab="",ylim = c(10,80))
mtext('残基位置',side = 1,line = 0.24,cex = 1.12)
mtext('温度因子',side = 2,line = 0.48,cex=1.12)
axis(1, at=c(1,35,70,105,140,175,210,245,280,307), 
     labels=c(1,35,NA,105,140,175,210,245,280,307),las=0,lwd=0.5,padj=-1)
axis(2,at=c(10,20,30,40,50,60,70,80), 
     labels=c(NA,20,30,40,50,60,70,80),las=0,lwd=0.5)
lines(v2, type="l", col="red",lty=2)
legend("topleft", inset=-.03, c("实验值","计算值"),
       lty=c(2, 1), col=c("red", "black"),bty="n",x.intersp = 0.01, 
       y.intersp =c(0.01,0.31))
dev.off()

###4UG2温度因子图
data2<-read.table("4ug2_b.txt",header = F)
v4<-append(data2$V4, c(rep(NA,12)), after = 203);v2<-append(data2$V2, c(rep(NA,12)), after = 203)
tiff("Bfactor2.tif",width = 5.5,height = 4,units="cm",res=800,compression = "lzw")
par(ps=7.5,family="serif",mai=c(0.25,0.25,0.02,0.02),tck=0.015,mgp=c(0.55,0,0),lwd=0.5)
plot(v4,type = "l",xaxt="n",yaxt="n",xlab = "",ylab="",ylim = c(20,147))
mtext('残基位置',side = 1,line = 0.24,cex = 1.12)
mtext('温度因子',side = 2,line = 0.48,cex=1.12)
axis(1, at=c(1,30,65,100,135,170,205,240,275,304), 
     labels=c(6,35,NA,105,140,175,210,245,280,309),las=0,lwd=0.5,padj=-1)
axis(2,at=c(20,40,60,80,100,120,140), 
     labels=c(20,NA,60,NA,100,NA,140),las=0,lwd=0.5)
lines(v2, type="l", col="red",lty=2)
legend("topleft", inset=-.03, c("实验值","计算值"),
       lty=c(2, 1), col=c("red", "black"),bty="n",x.intersp = 0.01, 
       y.intersp =c(0.01,0.31))
dev.off()

###5G53温度因子图
data3<-read.table("5g53_b.txt",header = F)
v4<-append(data3$V4, c(rep(NA,12)), after = 141); v2<-append(data3$V2, c(rep(NA,12)), after = 141)
v4<-append(v4, c(rep(NA,12)), after = 206);v2<-append(v2, c(rep(NA,12)), after = 206)
tiff("Bfactor3.tif",width = 5.5,height = 4,units="cm",res=800,compression = "lzw")
par(ps=7.5,family="serif",mai=c(0.25,0.25,0.02,0.02),tck=0.015,mgp=c(0.55,0,0),lwd=0.5)
plot(v4,type = "l",xaxt="n",yaxt="n",xlab = "",ylab="",ylim = c(25,255))
mtext('残基位置',side = 1,line = 0.24,cex = 1.12)
mtext('温度因子',side = 2,line = 0.48,cex=1.12)
axis(1, at=c(1,30,65,100,135,170,205,240,275,307), 
     labels=c(6,35,70,105,140,175,210,245,280,312),las=0,lwd=0.5,padj=-1)
axis(2,at=c(50,100,150,200,250), 
     labels=c(50,100,150,200,250),las=0,lwd=0.5)
lines(v2, type="l", col="red",lty=2)
legend("topleft", inset=-.03, c("实验值","计算值"),
       lty=c(2, 1), col=c("red", "black"),bty="n",x.intersp = 0.01, 
       y.intersp =c(0.01,0.31))
dev.off()


###变形能、柔性度
pdbf <- read.pdb('4eiycut.pdb')#4eiycut.pdb分别替换为4ug2cut.pdb和5g53cut.pdb
#min-max标准化
normalize <- function(x){
  return ((x - min(x)) / (max(x) - min(x)))
}
aa<-c()
for (i in 10:35) {
  anmf<-nma.pdb(pdbf, ff="anm", mass=T, temp=NULL,cutoff=i)
  defe <- deformation.nma(anmf)
  defsums <- rowSums(defe$ei[,1:20])
  defsums<-normalize(defsums)
  fluct<-1/2*normalize(fluct.nma(anmf, mode.inds=sample(7:length(anmf[["frequencies"]]),50)))+1/2*defsums
  aa=rbind(aa,fluct)
  aa<-as.data.frame(aa)
}
fluct<-colSums(aa)
fluct<-apply(aa,2,mean)
write.pdb(pdb=NULL, xyz=anmf$xyz, file="R-fluct柔性度.pdb", b=fluct)

###*****快慢运动模式图*******
pdbf <- read.pdb('4eiycut.pdb')
pdbm <- read.pdb('4ug2cut.pdb')
pdbt <- read.pdb("5g53cut.pdb")
anmf<-nma.pdb(pdbf, ff="anm", mass=F, temp=NULL,cutoff=30)
anmm<-nma.pdb(pdbm, ff="anm", mass=F, temp=NULL,cutoff=26)
anmt<-nma.pdb(pdbt, ff="anm", mass=F, temp=NULL,cutoff=28)
slowf <- fluct.nma(anmf, mode.inds=seq(7,8))#前2种慢模式下的波动
fastf <- fluct.nma(anmf, mode.inds=seq((length(anmf[["L"]])-1),length(anmf[["L"]])))
slowm <- fluct.nma(anmm, mode.inds=seq(7,8))
fastm <- fluct.nma(anmm, mode.inds=seq((length(anmm[["L"]])-1),length(anmm[["L"]])))
slowt <- fluct.nma(anmt, mode.inds=seq(7,8))
fastt <- fluct.nma(anmt, mode.inds=seq((length(anmt[["L"]])-1),length(anmt[["L"]])))
slowf<-append(slowf, c(rep(NA,10)), after = 208);slowf<-append(slowf, c(rep(NA,5)), after = 307)
fastf<-append(fastf, c(rep(NA,10)), after = 208);fastf<-append(fastf, c(rep(NA,5)), after = 307)
slowm<-append(slowm, c(rep(NA,5)), after = 0);slowm<-append(slowm, c(rep(NA,12)), after = 208);slowm<-append(slowm, c(rep(NA,3)), after = 309)
fastm<-append(fastm, c(rep(NA,5)), after = 0);fastm<-append(fastm, c(rep(NA,12)), after = 208);fastm<-append(fastm, c(rep(NA,3)), after = 309)
slowt<-append(slowt, c(rep(NA,5)), after = 0);slowt<-append(slowt, c(rep(NA,12)), after = 146);slowt<-append(slowt, c(rep(NA,12)), after = 211)
fastt<-append(fastt, c(rep(NA,5)), after = 0);fastt<-append(fastt, c(rep(NA,12)), after = 146);fastt<-append(fastt, c(rep(NA,12)), after = 211)
#所有结合位点
naall<-c(24,27,48,49,51,52,91,95,242,246,280,281,284,285,288)#钠离子位点
sitef<-c(66,81,84,85,88,153,167,168,169,170,174,177,181,246,
         249,250,252,253,256,264,265,267,270,271,274,277)#配体结合位点
sitem<-c(9,13,59,63,66,67,68,81,84,85,86,88,89,90,92,167,168,
         169,170,174,177,181,185,186,246,249,250,252,253,256,
         264,265,266,267,270,271,274,277,278)
sitet<-c(9,13,59,63,66,81,84,85,86,88,89,92,168,169,174,177,181,185,
         186,246,249,250,252,253,256,264,270,273,274,277,278)
bing<-sort(union(union(sitef,sitem),sitet))
allr<-sort(union(naall,bing));length(allr)

#****************快运动图****************
tiff("fig1.tif",width = 5.5,height = 4,units="cm",res=800,compression = "lzw")
par(ps=7.5,family="serif",mai=c(0.25,0.25,0.02,0.02),tck=0.015,mgp=c(0.55,0,0),lwd=0.5)
plot(fastf,type = "l",ylim = c(0.000,0.00055),xaxt="n",yaxt="n",xlab = "",ylab="")/
  mtext('残基位置',side = 1,line = 0.24,cex = 1.12)/
  mtext('均方涨落',side = 2,line = 0.48,cex=1.12)/
  axis(2,at=c(0,0.0001,0.0002,0.0003,0.0004,0.0005), 
       labels=c(0,NA,0.0002,NA,0.0004,0.0005),las=0,lwd=0.5)/
  lines(fastm, type="l", col='blue',lty=2)/
  lines(fastt, type="l", col=rainbow(50)[7],lty=3)/
  axis(1, at=c(1,50,100,150,200,250,300,312), 
       labels=c(1,50,100,150,200,250,NA,312),las=0,lwd=0.5,padj=-1)/
  legend("topright", inset=-0.01, c("4EIY","4UG2","5G53"),x.intersp = 0.01,y.intersp =c(0.01,0.26,0.31),
         lty=c(1,2,3), col=c( "black","blue",rainbow(50)[7]),bty="n")
points(c(24,48,52,91,281,284,285,63,84,88,185,92,168,181,186,278),
       fastm[c(24,48,52,91,281,284,285,63,84,88,185,92,168,181,186,278)],col="black",cex=0.3)/
  points(c(27,242,246),fastt[c(27,242,246)],col="black",cex=0.3)
points(c(55,128),fastm[c(55,128)],col="red",cex=0.3)
points(c(17,20,132,189,197,239),fastt[c(17,20,132,189,197,239)],col="red",cex=0.3)
dev.off()

###****************慢运动图****************
tiff("fig2.tif",width = 5.5,height = 4,units="cm",res=800,compression = "lzw")
par(ps=7.5,family="serif",mai=c(0.25,0.25,0.055,0.02),tck=0.015,mgp=c(0.55,0,0),lwd=0.5)
plot(slowf,type = "l",ylim = c(0.000,0.02),xaxt="n",yaxt="n",xlab = "",ylab="",lwd=0.35)/
  mtext('残基位置',side = 1,line = 0.24,cex = 1.12)/
  mtext('均方涨落',side = 2,line = 0.48,cex=1.12)/
  axis(2,at=c(0,0.005,0.01,0.015,0.02), 
       labels=c(0,0.005,0.01,0.015,0.02),las=0,lwd=0.5)/
  lines(slowm, type="l", col='blue',lty=2)/
  lines(slowt, type="l", col=rainbow(50)[7],lty=3)/
  axis(1, at=c(1,50,100,150,200,250,300,312), 
       labels=c(1,50,100,150,200,250,NA,312), las=0,lwd=0.5,padj=-1)/
  legend("topleft", inset=0, c("4EIY","4UG2","5G53"),x.intersp = 0.01,y.intersp =c(0.01,0.26,0.31),
         lty=c(1,2,3), col=c( "black","blue",rainbow(50)[7]),bty="n") 
points(naall,slowt[naall],col="black",cex=0.3)
points(bing,slowt[bing],col="black",cex=0.3)
points(c(20,55,128,132,189,197,239),slowt[c(20,55,128,132,189,197,239)],col="red",cex=0.3)
points(c(17),slowf[c(17)],col="red",cex=0.3)
points(c(246),slowt[c(246)],col="black",cex=0.3)
dev.off()

#**************相关性热图****************
#基于GPCRdb数据库同源建模结构
pdbf <- read.pdb('4EIYDB.pdb')
pdbm <- read.pdb('4UG2DB.pdb')
pdbt <- read.pdb("5G53DB.pdb")
anmf<-nma.pdb(pdbf, ff="anm", mass=F, temp=NULL,cutoff=30)
anmm<-nma.pdb(pdbm, ff="anm", mass=F, temp=NULL,cutoff=26)
anmt<-nma.pdb(pdbt, ff="anm", mass=F, temp=NULL,cutoff=28)
cm <- dccm(anmf)
#write.csv(cm,"关联性分析.csv")
colnames(cm)=c(seq(1,314,1))
rownames(cm)=c(seq(1,314,1))
plot(cm,xlab = "Residue",ylab="Residue",axes=F,
     contour=F, col.regions=bwr.colors(2000), at=seq(-1,1,0.01) )
# View the correlations in the structure
pymol.dccm(cm, pdb, launch=TRUE)

#将GPCRdb数据库的相互作用接触残基对拆分
pos<-read.csv("接触对残基.csv",header=T)
v<-unlist(strsplit(pos[1,],split ='-'))
a<-c()
for (i in 1:nrow(pos)) {
  v<-unlist(strsplit(pos[i,],split ='-'))
  a<-rbind(a,v)
}
a
#write.csv(a,"caifen的表.csv")
#**************绘制相互作用图****************
san<-read.csv("pos散点数据.csv",header=T)
pos1<-san$pos1
pos2<-san$pos2
pos3<-c(pos1,pos2)
pos4<-c(pos2,pos1)
tiff("散点1.tif",width = 3.68,height = 3.7,units="cm",res=900,compression = "lzw")
par(ps=7.5,family="serif",mai=c(0.25,0.13,0.02,0.07),tck=0.015,mgp=c(0.55,0,0),lwd=0.5)
#plot(slowf,type = "l",ylim = c(0.000,0.02),xaxt="n",yaxt="n",xlab = "",ylab="均方涨落",lwd=0.35)/
plot(pos4~pos3,xaxt="n",yaxt="n",xlab = "",ylab="",cex=0.09,pch=19,col="red")/
  mtext('残基位置',side = 1,line = 0.23,cex = 1.12)/
  axis(2,at=c(1,50,100,150,200,250,300,314), 
       labels=c(1,NA,100,150,200,250,300,NA),las=0,lwd=0.5)/
  axis(1, at=c(1,50,100,150,200,250,300,314), 
       labels=c(1,50,100,150,200,250,NA,314), las=0,lwd=0.5,padj=-1)
dev.off()




