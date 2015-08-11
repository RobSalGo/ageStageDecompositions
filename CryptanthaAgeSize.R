#Building an age by size mega-matrix IPM
#Code created by Rob Salguero-Gomez (r.salguero@uq.edu.au)
#Creation date: Jan 27 2015
#Last modified: Jul 31 2015


Hal: get a set of matrices for the J Ecol paper.
Combine individuals > 50 rosettes.
Reconstruct age only and size only matrices based on vital rate functions.
Anonymous reproduction.



#Clear all previous content
rm(list=ls(all=TRUE))
      
library(IPMpack)
library(plotrix)
library(scales)

###Functions

picGrow2 <- function (dataf, growObj, mainTitle = "Growth", ...) 
{
    colfunc <- colorRampPalette(c("Goldenrod 1", "Dark green"))
    colors=colfunc(length(unique(dataf$covariate)))
    colfunc2 <- colorRampPalette(c("Goldenrod 1", "Dark green"))
    colors2=colfunc2(length(unique(dataf$covariate)))
    predVar <- attr(growObj@fit$terms, "predvars")[[2]]
    plot(jitter(dataf$size,0.1), jitter(dataf$sizeNext,0.1), pch = 19, xlab = "Size at t", 
        ylab = "Size at t+1", col="white",main = mainTitle)
    abline(a = 0, b = 1)
    
    dataf$covariate <- as.factor(dataf$covariate)
    levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
    ud <- unique(dataf$covariate)
    ud <- ud[!is.na(ud)]
    for (k in 1:length(ud)) {
        tp <- dataf$covariate == ud[k]
        points(jitter(dataf$size[tp],0.1), jitter(dataf$sizeNext[tp],0.1), pch = 19, 
            col = alpha(colors[k],0.5))
                #col = paste("gray",round(105-6.6*k,0)))
    }
    
        ud <- as.factor(ud)
        
    sizes <- dataf$size[!is.na(dataf$size)]
    sizes <- sizes[order(sizes)]
    for (k in 1:length(ud)) {
        newd <- data.frame(size = sizes, size2 = sizes^2, size3 = sizes^3, 
            covariate = as.factor(rep(ud[k], length(sizes))))
        if (length(grep("expsize", names(growObj@fit$coefficients))) == 
            1) 
            newd$expsize = exp(sizes)
        if (length(grep("logsize", names(growObj@fit$coefficients))) == 
            1) 
            newd$logsize = log(sizes)
        if (length(grep("logsize2", names(growObj@fit$coefficients))) == 
            1) 
            newd$logsize = (log(sizes))^2
        if (length(grep("decline", tolower(as.character(class(growObj))))) > 
            0 | length(grep("trunc", tolower(as.character(class(growObj))))) > 
            0) {
            pred.size <- .predictMuX(growObj, newd, covPred = k)
        }
        else {
            pred.size <- predict(growObj@fit, newd, type = "response")
        }
        if (length(grep("incr", tolower(as.character(class(growObj))))) == 
            0) {
            points(sizes, pred.size, type = "l",
                col = colors2[k],lwd=3)
                #col = paste("gray",round(105-6.6*k,0)),lwd=3)
        }
        else {
            if (length(grep("logincr", tolower(as.character(class(growObj))))) > 
                0) {
                points(sizes, sizes + exp(pred.size), type = "l", 
                  col = colors2[k],lwd=3)
                #col = paste("gray",round(105-6.6*k,0)),lwd=3)
            }
            else {
                lines(sizes, pred.size, col = colors2[k])
            }
        }
    }

legend("bottomright",legend=levels(ud),col=colors,pch=16,cex=.7,bg="white")
}





picSurv2 <- function (dataf, survObj, ncuts = 20, makeTitle = "Survival", 
    ...) 
{
    os <- order(dataf$size)
    os.surv <- (dataf$surv)[os]
    os.size <- (dataf$size)[os]
    psz <- tapply(os.size, as.numeric(cut(os.size, ncuts)), mean, 
        na.rm = TRUE)
    ps <- tapply(os.surv, as.numeric(cut(os.size, ncuts)), mean, 
        na.rm = TRUE)
    if (length(grep("covariate", names(survObj@fit$model))) == 
        0) {
        plot(as.numeric(psz), as.numeric(ps), pch = 19, xlab = "Size at t", 
            ylab = "Survival to t+1", main = makeTitle, ...)
        points(dataf$size[order(dataf$size)], surv(dataf$size[order(dataf$size)], 
            data.frame(covariate = 1), survObj), type = "l", 
            col = 2)
    }
    else {
        plot(as.numeric(psz), as.numeric(ps), type = "n", pch = 19, 
            xlab = "Size at t", ylab = "Survival to t+1", main = "Survival", 
            ...)
        dataf$covariate <- as.factor(dataf$covariate)
        levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
        os.cov <- (dataf$covariate)[os]
        sizes <- dataf$size[!is.na(dataf$size)]
        sizes <- sizes[order(sizes)]
        ud <- unique(dataf$covariate)
        ud <- ud[!is.na(ud)]
        for (k in 1:length(ud)) {
            tp <- os.cov == ud[k]
            psz <- tapply(os.size[tp], as.numeric(cut(os.size[tp], 
                ncuts)), mean, na.rm = TRUE)
            ps <- tapply(os.surv[tp], as.numeric(cut(os.size[tp], 
                ncuts)), mean, na.rm = TRUE)
            points(as.numeric(psz), as.numeric(ps), pch = 19, 
                col = paste("gray",round(100-6.6*k,0)))
            newd <- data.frame(size = sizes, size2 = sizes^2, 
                size3 = sizes^3, covariate = rep(as.factor(ud[k]), 
                  length(sizes)))
            if (length(grep("expsize", survObj@fit$formula)) == 
                1) 
                newd$expsize = exp(sizes)
            if (length(grep("logsize", survObj@fit$formula)) == 
                1) 
                newd$logsize = log(sizes)
            if (length(grep("logsize2", survObj@fit$formula)) == 
                1) 
                newd$logsize = (log(sizes))^2
            pred.surv <- predict(survObj@fit, newd, type = "response")
            points(newd$size, pred.surv, type = "l", col = paste("gray",round(100-6.6*k,0)))
        }
    }
}


#Read data
setwd("~/Dropbox/ageStageDecompositions")

#Read longterm data of Cryptantha flava. Make sure you understand waht the metadata are (also in dropbox)
d <- read.table("Population_dynamics_data.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)

#Some reformatting needs to happen. Also, note that in 2002 no data was collected.

d$sizeNext=NA
for (i in 1998:2012){
  d$sizeNext[which(d$Year==i-1)]=d$Size[which(d$Year==i)]
}

d$surv=NA
d$surv[which(!is.na(d$Size) & !is.na(d$sizeNext))]=1
d$surv[which(!is.na(d$Size) & is.na(d$sizeNext))]=0

d$fec0=NA
d$fec0[which(d$Fert>0)]=1
d$fec0[which(d$Fert==0)]=0

d$Fert[(d$Fert==0)]=NA
colnames(d)

d=d[,c("ID","Treatment","Year","Age","Size","sizeNext","surv","fec0","Fert")]
colnames(d)=c("ID","treatment","year","covariate","size","sizeNext","surv","fec0","fec1")

#logging sizes
d$size=log(d$size+1)
d$sizeNext=log(d$sizeNext+1)

#Ignoring individuals for which no good data exists
d=d[-which(d$covariate==999),]
d=d[-which(d$year%in%c(2001:2003)),]
d=d[-which(d$year==2012),]

go = makeGrowthObj(dataf = d,Formula=sizeNext~size+size2+size3+covariate)
picGrow2(d,go)
#Rob: check individuals of small size growing too quickly
#Add legend

so = makeSurvObj(d, Formula = surv~size+size2+size3+covariate)
picSurv2(d,so,ylim=c(0,1))


#Building age transition matrix
  nAgeClasses=max(d$covariate,na.rm=TRUE)
  ageMat=new("envMatrix", nEnvClass = nAgeClasses)
  ageMat@.Data <- matrix(0,nAgeClasses,nAgeClasses)
  ageMat@.Data[cbind(2:nAgeClasses,1:(nAgeClasses-1))] <- 1
  #ageMat@.Data[n.age.classes,n.age.classes] <- 1
  ageMat

Pmatrix=makeCompoundPmatrix(nBigMatrix = 100, minSize = 1,maxSize = 100,envMatrix = ageMat,growObj = go,survObj = so,correction = "changingVar")
dim(Pmatrix)

image(as(Pmatrix[,],'sparseMatrix'),
	xlab = "Continuous stage (e.g. size) at t",
		ylab = "Continuous stage (e.g. size) at t+1", axes = FALSE)

#Fecundity
fo=makeFecObj(d, Formula = c(fec0~size+covariate,fec1~size+size2+size3+covariate),Family = c("binomial","gaussian"),Transform = c("none","none"))

nAgeClasses=max(d$covariate,na.rm=TRUE)
 ageMat1=new("envMatrix", nEnvClass = nAgeClasses)
 ageMat1@.Data <- matrix(0,nAgeClasses,nAgeClasses)
ageMat1@.Data[1,2:nAgeClasses] <- 1
Fmatrix <- makeCompoundFmatrix(nBigMatrix = 100, minSize = 1,maxSize = 100,envMatrix = ageMat1,fecObj = fo,correction = "constant")

#Plot it all together
IPM <- Pmatrix+Fmatrix

image(as(log(IPM),'sparseMatrix'),
	xlab = "Size within age at t",
		ylab = "Size within age at t+1", axes = FALSE)

library(popbio)
elast=elasticity(IPM)
image(as(log(elast),'sparseMatrix'),
  xlab = "Size within age at t",
    ylab = "Size within age at t+1", axes = FALSE)
