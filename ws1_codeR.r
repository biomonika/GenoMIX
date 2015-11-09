#! /bin/Rscript

#Makes PNAS submission Figure 1 (age correlation)

#rm(list=ls())
args = commandArgs(trailingOnly = TRUE)
data=read.delim(args[1],header=T)
#data=read.delim("~/Desktop/genomix/Fig2_main_data.txt",header=T)

colnames(data)=c("sample","class","hqc","hqb","hqt","qtc","qtb","qtt","age","con")

# sample: M-identifier of sample. i.e. M172
# class: mother (M) or child (C)
# hqc: high quality cheek (number of sites MAF>=1%)
# hqb: high quality blood (number of sites MAF>=1%)
# hqt: high quality total (number of sites MAF>=1%)
# qtc: quartet cheek (number of sites MAF>=0.2%)
# qtb: quartet blood (number of sites MAF>=0.2%)
# qtt: quartet total (number of sites MAF>=0.2%)
# age: age of individual in days
# con: conception (age of mother at fertilization of child in days)

#############################################################
#              THE PLOT IS AGE (yr) vs. HQT                 #
#############################################################

# Define colors

black="black"
m_col="royalblue1"
c_col="tomato1"
mother=rgb(matrix(col2rgb(m_col),1,3),alpha=120,maxColorValue=255)
child=rgb(matrix(col2rgb(c_col),1,3),alpha=120,maxColorValue=255)
borders=c(m_col,c_col)

# Plot to PDF

pdf("Fig2_main.pdf",height=3.42,width=3.42)

# Margins

par(oma=c(0,0,0,0))
par(mar=c(0.5,2.1,0,0))

# Blank plot
plot(
  1:10,1:10,
  xlim=c(15,60),ylim=c(-1.5,5),
  type="n",frame=F,
  axes=F,xlab="",ylab="",main="")

# Plot mothers data, calculate glm for the mothers, and plot a fitted line

mothers=data[data$class=="M",]
points(mothers$age/365,mothers$hqt,pch=23,col=borders[1],lwd=1,cex=1,bg=mother)
glm_mothers=glm(hqt~age,data=mothers,family="poisson")     
fit_mothers=data.frame(age=glm_mothers$data$age/365,f=glm_mothers$fitted.values)
fit_mothers=fit_mothers[order(fit_mothers$age),]
lines(fit_mothers$age,fit_mothers$f,col=m_col,lwd=2.5)
#pval_mothers=round(summary(glm_mothers)$coefficients[2,4],2)
pval_mothers=summary(glm_mothers)$coefficients[2,4]

# Plot children data on top of current window
par(new=T)
children=data[data$class=="C",]
points(children$con/365,children$hqt,pch=21,col=borders[2],bg=child,lwd=1,cex=1)
glm_children=glm(hqt~con,data=children,family="poisson")
fit_children=data.frame(con=glm_children$data$con/365,f=glm_children$fitted.values)
fit_children=fit_children[order(fit_children$con),]
lines(fit_children$con,fit_children$f,col=c_col,lwd=2.5)
#pval_children=round(summary(glm_children)$coefficients[2,4],3)
pval_children=summary(glm_children)$coefficients[2,4]

# Axis format

axis(side=1,at=seq(15,60,by=5),lab=NA,lwd=2,pos=-0.25,cex.axis=0.75)
axis(side=2,at=0:5,lwd=2,pos=14.5,las=2,cex.axis=0.75)
mtext("number of point heteroplasmies",2,1.3,at=2.5,cex=0.75,font=2)
mtext("maternal age (years)",1,-0.5,cex=0.75,font=2)
for (i in seq(15,60,by=5)){text(i,-0.75,lab=i,cex=0.75)}

# Legend
M=sprintf("mother p=%.3f",pval_mothers)
C=sprintf("child     p=%.3f",pval_children)
legend(15,4.6,legend=c(M,C),col=borders,
       pt.bg=c(mother,child),pch=c(23,21),pt.cex=1.2,pt.lwd=1.2,
       cex=0.6,title="Poisson model",box.col="gray",box.lwd=1.5,bg=NA)

# Age ranges
x1=min(fit_mothers$age)
x2=max(fit_mothers$age)
x3=min(fit_children$con)
x4=max(fit_children$con)

lines(c(x1,x2),c(-1,-1),lwd=1.5,col=m_col)
lines(c(x3,x4),c(-1.1,-1.1),lwd=1.5,col=c_col)
mtext("fertilization",1,-1.25,at=27,cex=0.75,col=c_col,font=2)
mtext("collection",1,-1.25,at=45,cex=0.75,col=m_col,font=2)

dev.off()
