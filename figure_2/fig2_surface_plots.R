## Loading libraries
library(mgcv)
library(fields)

## load data
dat <- read.csv("data/fig2_gfw_urea_food_intake.csv")

k1<-9
k2<-6

model1 <- gam(urea ~ s(prot, k=k1, bs="tp") + s(carb, k=k1, bs="tp") + s(prot, carb, k=k2, bs="tp"),
            gamma = 1.0, method = "REML", select = TRUE, data = dat)

model2 <- gam(urea ~ s(prot, k=k1, bs="tp") + s(fat, k=k1, bs="tp") + s(prot, fat, k=k2, bs="tp"),
              gamma = 1.0, method = "REML", select = TRUE, data = dat)

surf.te <- Tps(cbind(dat$prot, dat$carb), dat$urea, lambda=0.01)
summary(surf.te)
surf.te.out=predictSurface(surf.te)

image(surf.te.out, col=tim.colors(128), lwd=5, las=1, font.lab=2, cex.lab=1.3, mgp=c(2.7,0.5,0), 
      font.axis=1, lab=c(4,5,6), xlim=c(0.0,0.095), ylim=c(0,1.5), xlab=expression("Protein intake (kcal/g)"), 
      ylab=expression("Carbohydrate intake (kcal/g)"), asp=0)
title("Urea")
contour(surf.te.out, lwd=2, labcex=1, add=T)
points(dat$prot,dat$carb, cex=1.0, pch=1, col="black")

surf.te <- Tps(cbind(dat$prot, dat$fat), dat$urea, lambda=0.01)
summary(surf.te)
surf.te.out=predictSurface(surf.te)

image(surf.te.out, col=tim.colors(128), lwd=5, las=1, font.lab=2, cex.lab=1.3, mgp=c(2.7,0.5,0), 
      font.axis=1, lab=c(4,5,6), xlim=c(0.0,0.095), ylim=c(0,0.06), xlab=expression("Protein intake (kcal/g)"), 
      ylab=expression("Fat intake (kcal/g)"), asp=0)
title("Urea")
contour(surf.te.out, lwd=2, labcex=1, add=T)
points(dat$prot,dat$intake.fat, cex=1.0, pch=1, col="black")
