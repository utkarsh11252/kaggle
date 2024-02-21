---
  title: "Prac"
output:
  pdf_document:
  latex_engine: xelatex # or pdflatex
---


#install.packages("gplots")
#install.packages("factoextra")
#install.packages("MVN")

#library(gplots)
#library(factoextra)
#library(MVN)

source(
  "http://www1.maths.leeds.ac.uk/~john/3772/clusterings.r")


nutrients = read.csv("http://www1.maths.leeds.ac.uk/~john/3772/nutrients.csv")
nutrients = nutrients[,2:6]

head(nutrients)
attach(nutrients)

#testing independence
chisq.test(nutrients)

summary(nutrients)

#inspecting univarate normality

# par(mfrow=c(2,3))
# qqnorm(Calcium, main = "Q-Q Plot for Calcium")
# qqline(Calcium, col = "red")
# hist(Calcium,breaks=30)
# boxplot(Calcium, main = "Box plot for Calcium")
# 
# qqnorm(Iron, main = "Q-Q Plot for Iron")
# qqline(Iron, col = "red")
# hist(Iron,breaks=30)
# boxplot(Iron, main = "Box plot for Iron")
# 
# qqnorm(Protein, main = "Q-Q Plot for Protein")
# qqline(Protein, col = "red")
# hist(Protein,breaks=30)
# boxplot(Protein, main = "Box plot for Protein")
# 
# qqnorm(VitaminA, main = "Q-Q Plot for Vitamin A")
# qqline(VitaminA, col = "red")
# hist(VitaminA,breaks=30)
# boxplot(VitaminA, main = "Box plot for Vitamin A")
# 
# qqnorm(VitaminC, main = "Q-Q Plot for Vitamin C")
# qqline(VitaminC, col = "red")
# hist(VitaminC,breaks=30)
# boxplot(VitaminC, main = "Box plot for Vitamin C")

#inspecting bivariate normality
#mvn(data=nutrients, mvnTest="mardia")

# pairs(nutrients, main="Pairwise scatter plot of all variables")
# 
# cor(nutrients)
# heatmap.2(cor(nutrients), 
#           col = colorRampPalette(c("red", "white", "green"))(100), 
#           key = TRUE,
#           srtRow=34,
#           srtCol = 34,
#           density.info = "density",  
#           trace = "none",
#           dendrogram = "none",
#           key.title="",
#           margins = c(5, 10),
#           main="Heatmap for the data"
#           )


#hotellings t2 test

n = nrow(nutrients)
rda = c(Cal=1000, Iro=15, Pro=60, VitA=800, VitC=75)
xbar = colMeans(nutrients)

diff = xbar - rda
S = round(cov(nutrients),3)

tsq = n*t(diff)%*%solve(S)%*%diff
tsq

p=5
m=n-1

fstat = tsq*(m-p+1)/(m*p)
fstat
pf(fstat, df1=p, df2=m-p+1, lower.tail=FALSE)

xbar
#calculating sci

inv.tsq = qf(0.95, df1=p, df2=m-p+1) * (p*m)/(m-p+1)
inv.tsq

sci1 = sqrt(S[1,1] * inv.tsq/n)
c(xbar[1]-sci1,xbar[1]+sci1)
rda[1]

sci2 = sqrt(S[2,2] * inv.tsq/n)
c(xbar[2]-sci2,xbar[2]+sci2)
rda[2]

sci3 = sqrt(S[3,3] * inv.tsq/n)
c(xbar[3]-sci3,xbar[3]+sci3)
rda[3]

sci4 = sqrt(S[4,4] * inv.tsq/n)
c(xbar[4]-sci4,xbar[4]+sci4)
rda[4]

sci5 = sqrt(S[5,5] * inv.tsq/n)
c(xbar[5]-sci5,xbar[5]+sci5)
rda[5]

#var high for cal, vit C
#vit A and vit C not significant diff
#asssumptions: multivariate normality, independence


#testing equality of means

st.nut = data.frame(mapply("/",nutrients,rda))
#st.nut
attach(st.nut)

A = matrix(c(1,-1,0,0,0,0,1,-1,0,0,0,0,1,-1,0,0,0,0,1,-1),nrow=4,byrow=T)
A
dif.nut = data.frame(as.matrix(st.nut)%*%t(A))    
colnames(dif.nut) = c("Calcium-Iron","Iron-Protein","Protein-VitaminA","VitaminA-VitaminC")
#dif.nut

#bivariate plot for normality
# par(mfrow = c(2,2))
# hist(dif.nut[,1], main="Calcium-Iron",xlab="")
# hist(dif.nut[,2], main="Iron-Protein",xlab="")
# hist(dif.nut[,3], main="Protein-Vitamin A",xlab="")
# hist(dif.nut[,4], main="Protein-Vitamin C",xlab="")


xbar.s = colMeans(dif.nut)
mu0 = rep(0,4)
dif.s = xbar.s-mu0
S.s = round(cov(dif.nut),3)

tsq.s = n*t(dif.s)%*%solve(S.s)%*%dif.s
tsq.s

p.s = 4
fstat.s = tsq.s*(m-p.s+1)/(m*p.s)
fstat.s
pf(fstat.s, df1=p.s, df2=m-p.s+1, lower.tail=FALSE)

#sci 

inv.tsq.s = qf(0.95, df1=p.s, df2=m-p.s+1) * (p.s*m)/(m-p.s+1)
inv.tsq.s

sci1.s = sqrt(S.s[1,1] * inv.tsq.s/n)
c(xbar.s[1]-sci1.s,xbar.s[1]+sci1.s)

sci2.s = sqrt(S.s[2,2] * inv.tsq.s/n)
c(xbar.s[2]-sci2.s,xbar.s[2]+sci2.s)

sci3.s = sqrt(S.s[3,3] * inv.tsq.s/n)
c(xbar.s[3]-sci3.s,xbar.s[3]+sci3.s)

sci4.s = sqrt(S.s[4,4] * inv.tsq.s/n)
c(xbar.s[4]-sci4.s,xbar.s[4]+sci4.s)


#kmeans

sc.data = scale(nutrients)
#sc.data
x.cl = clustering(sc.data,"kmeans",k=2)

# fviz_cluster(kmeans(sc.data, centers = 2),data=sc.data)
# 
# plot(sc.data, col = kmeans(sc.data, centers = 2)$cluster, pch = 16, main = "K-Means Clusters")
# points(kmeans(sc.data, centers = 2)$means, col = 1:2, pch = 8, cex = 4)


cal.m = x.cl$means[,1]*sd(nutrients$Calcium)+mean(nutrients$Calcium)
ir.m = x.cl$means[,2]*sd(nutrients$Iron)+mean(nutrients$Iron)
pr.m = x.cl$means[,3]*sd(nutrients$Protein)+mean(nutrients$Protein)
va.m = x.cl$means[,4]*sd(nutrients$VitaminA)+mean(nutrients$VitaminA)
vc.m = x.cl$means[,5]*sd(nutrients$VitaminC)+mean(nutrients$VitaminC)

# barplot(matrix(c(cal.m,ir.m,pr.m,va.m,vc.m),nrow=2,byrow=F), beside=T, col = c("lightblue", "lightgreen"),names.arg=c("Calcium","Iron","Protein","Vitamin A","Vitamin C"),main="Mean values for the clusters")
# legend("topright", legend = c("c1", "c2"), fill = c("lightblue", "lightgreen"), title = "Clusters")




#knitr::stitch_rhtml('prac.r')        
knitr::stitch('prac.r')
browseURL('prac.pdf')

#tinytex::install_tinytex()
#tinytex:::install_prebuilt()
