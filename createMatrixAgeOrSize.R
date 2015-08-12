#Building an age by size mega-matrix IPM
#Code created by Rob Salguero-Gomez (r.salguero@uq.edu.au)
#Creation date: Jan 27 2015
#Last modified: Jul 31 2015

#Clear all previous content
rm(list=ls(all=TRUE))


#Function to create a population matrix model


makeMatrix <- function(d=d,discretize=c(1,2,3,4), variable="Size", seedlingSizeThreshold=NA){

	if (variable == "Age" & !is.na(seedlingSizeThreshold)) {print("seedlingSizeThreshold reset to NA for age-based models")}
	
	d1=d			#d[which(!is.na(d[,variable])),]
	d1$catNow=d1$catNext=NA

#Number of individuals transitioning
	#Evaluate only those that are transitioning
		indexTrans=which(!is.na(d[,variable]) & !is.na(d[,paste(variable,"Next",sep="")]))

		for (i in rev(discretize)) {
			for (j in indexTrans) {
				if (d1[,variable][j]<i) {d1$catNow[j]=i}
				if (d1$Surv[j]==0) {d1$catNext[j]=NA} else {
				if (d1[,paste(variable,"Next",sep="")][j]<i) {d1$catNext[j]=i}}
				print(paste(variable, "class", i, "of individual", j, "in dataset"))
			}
		}

	#d1$catNow[which(d1$catNow==max(c(d$Size,d$SizeNext),na.rm=T)+1)]=50
		d1$catNow[which(d1$catNow==max(c(d[,variable],d[,paste(variable,"Next",sep="")]),na.rm=T)+1)]=discretize[length(discretize)-1]

	#d1$catNext[which(d1$catNext==max(c(d$Size,d$SizeNext),na.rm=T)+1)]=50
		d1$catNext[which(d1$catNext==max(c(d[,variable],d[,paste(variable,"Next",sep="")]),na.rm=T)+1)]=discretize[length(discretize)-1]

	#Making sure that both catNow and catNext have all range of values for discretize
		#missingCatNow=which(!discretize %in% unique(d1$catNow))
		#missingCatNext=which(!discretize %in% unique(d1$catNext))

	#Frequency of transitions
		matU = as.matrix(table(d1$catNext,d1$catNow))

	#Making sure that all values of discretize are present in matU, even if not recorded transitions appear
	#for (i in discretize){
	#	if (!i %in% colnames(matU)) {
	#		matU=cbind(matU,rep(0,dim(matU)[1]))
	#		colnames(matU)[which(colnames(matU)=="")]=i
	#	}
	#	if (!i %in% rownames(matU)) {
	#		matU=rbind(matU,rep(0,dim(matU)[2]))
	#		rownames(matU)[which(rownames(matU)=="")]=i
	#	}
	#}

	#Reorganizing the transitions in order of discretize
	#matU=matU[order(as.numeric(rownames(matU))),order(as.numeric(colnames(matU)))]
	
	#Dividing each stage by its specific survival value
	surv = as.data.frame(colSums(table(d1$Surv,d1$catNow)))
	if (variable=="Age") {for (i in as.numeric(rownames(surv))){matU[,i]=matU[,i]/surv[which(rownames(surv)==as.character(i)),]}}
	if (variable=="Size") {matU=matU/surv[,1]}

#Reproduction through anonymous motherhood
	if (variable == "Size") {
		recruits=d$SizeNext[which(is.na(d$Size) & !is.na(d$SizeNext) & d$SizeNext<seedlingSizeThreshold)]
	}
	if (variable == "Age") {
		recruits=d$AgeNext[which(d$AgeNext==1)]
	}

	#Total number of recruits to attribute reproductive output from
	recruitNumber=length(recruits)
	
	#Total number of reproductive structures per stage
	rep=rep1=table(d1$Fert,d1$catNow)
	for (i in 1:dim(rep)[2]) {
		rep1[,i]=rep[,i]*as.numeric(rownames(rep))
	}
	
	#Total number of reproductive structures per capita in each stage
	rep2=as.matrix(colSums(rep1[2:dim(rep1)[1],]))
	#Proportion of reproductive effort per capita per stage
	rep3=rep2/sum(rep2)
	#Total of seedlings being produced per capita in each stage
	rep4=recruitNumber*rep3
	#Creating 0 F matrix
	matF=matrix(0,length(discretize)-1,length(discretize)-1)
	colnames(matF)=colnames(matF)=discretize[-length(discretize)]

	matF[1,]=rep4
	#for (i in as.numeric(rownames(rep4))){
	#	matF[1,i]=as.numeric(rep4[which(as.numeric(rownames(rep4))==i)])
	#}

	
	#Making A matrix
		matA = matU + matF
		out=list(matA,matU,matF)

	Re(eigen(matA)$values[1])



#Read data
setwd("~/Dropbox/ageStageDecompositions/")

data <- read.table("Population_dynamics_data.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
data$SizeNext=NA
data$AgeNext=NA

for (i in 1997:(2012-1)){
	data$SizeNext[which(data$Year==i)]=data$Size[which(data$Year==i+1)]
	data$AgeNext[which(data$Year==i)]=data$Age[which(data$Year==i+1)]
}

data$Age[which(d$Age==999)]=NA
data$AgeNext[which(d$AgeNext==999)]=NA

data$Surv=NA
data$Surv[which(is.na(data$Size))]=NA
data$Surv[which(!is.na(data$Size) & is.na(data$SizeNext))]=0
data$Surv[which(!is.na(data$Size) & !is.na(data$SizeNext))]=1


data$ID[which(data$Age>data$AgeNext)]

plot(data$Age,data$AgeNext)
plot(data$Size,data$SizeNext)


#Choosing to discretize the data by size or age:
	#By size:
		discretize=c(2,5,10,15,20,25,50,max(c(data$Size,data$SizeNext),na.rm=T)+1)
		seedlingSizeThreshold=5
	#By age
		discretize=1:max(c(data$Age,data$AgeNext),na.rm=T)
		seedlingSizeThreshold=NA

	#Take subsets by treament level:
		d=data
		d=data[which(data$Treatment=="C"),]		#For Control
		d=data[which(data$Treatment=="D1"),]	#For Drought in 1998
		d=data[which(data$Treatment=="D2"),]	#For Drought in 1999







