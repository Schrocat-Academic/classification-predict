#
# klaR/classifier.R
#
# Created by Michael Walker, 22/05/2014
# Emails: walkerm1@student.unimelb.edu.au, m.walker@aip.org.au
#
# Functions to interface with R naive Bayes function NaiveBayes()
# from package klaR. This is an upgrade of that from e1071.
#

library(klaR)
library(ROCR)
#library(pROC)

source("work/code/classifier/regression/regression.R")
#source("work/code/classifier/rank_normal.R")    # No longer using rank_normal transform
#source("work/code/classifier/classfunctions.R")    # Called by regression.R

call_klaR_assessment <- function(filenames = c("work/data/CAS_Immunology_db.xlsx"),observables=c(),response="crwz5",plot_title="",runLR=TRUE,runNB=TRUE,runIRN=FALSE,runARC=FALSE,runs=100,testfraction=3,overfitcheck=FALSE,conditions="",auconly=FALSE,grouping="",cutoffs=c(0.35,0.5,0.75),clt=FALSE){
	
	# Prepare data, observables, and priors
	if (length(observables)==0) observables = readobserves()	
	print(observables)
	filedata = readframe(filenames[[1]])
	datafile = preparedata(filedata,observables,response,conditions=conditions)
	priors = findpriors(datafile[response])
	
	klaR_assessment(datafile,observables,runLR=runLR,runNB=runNB,runIRN=runIRN,runARC=runARC,runs=runs,testfraction=testfraction,priors=priors,response=response,plot_title=plot_title,overfitcheck=overfitcheck,conditions=conditions,auconly=auconly,cutoffs=cutoffs,clt=clt)
}

klaR_assessment <- function(datafile,observables=c(),response = "crwz5",plot_title=paste0(c(response,"from",observables),collapse=" "),runLR=TRUE,runNB=TRUE,runIRN=FALSE,runARC=TRUE,runs = 100,testfraction = 3,priors = 0,overfitcheck=FALSE,conditions="",auconly=FALSE,datavars=TRUE,grouping=c("subjid","SUBJ_ID"),cutoffs=c(0.35,0.5,0.75),clt=FALSE){
	# Runs cross-validation on datafile based on observables
	# clt refers to whether output should be suitable for studying the Central Limit Theorem
	#in relation to the AUC distribution from cross-validation.
	
	initdatarows = rownames(datafile)
	#cat("DEBUG: rownames",initdatarows,'\n')
	#print("responseframe")
	responseframe = datafile[response][[1]]
	if (length(observables)==0){
		if (datavars) observables = names(datafile)[names(datafile) != response]
		else observables = readobserves()	
	}
	#print(observables)
	grouping = grouping[grouping %in% names(datafile)]
	if (length(grouping)>0) groupingsframe = datafile[grouping]
	else grouping = c("")
	
	datafile = datafile[c(observables,response)]
	# Generate asinh transformed data
	asinhdata = asinh(datafile[observables])
	asinhdata = cbind(asinhdata,datafile[response])
	
	#print(responseframe)
	names(responseframe) = rownames(datafile)
	if (length(priors) == 1) priors = findpriors(datafile[response])

	posteriorslistNB = list()		# List of posteriors from each run
	posteriorslistLR = list()		# List of posteriors from each run
	#posteriorslistNBparam = list()		# List of posteriors from each run
	posteriorslistNBnorm = list()		# List of posteriors from each run
	asinhposteriorslistNB = list()		# List of posteriors from each run
	asinhposteriorslistLR = list()		# List of posteriors from each run
	responselist = list()		# List of correct values of response
	
	datarows = initdatarows
	for (run in 1:runs){
		# Do we have correlated samples?
		posrows = positive_rows(responseframe,1)
		if (grouping != "") {
			datarows = initdatarows[declusterize(groupingsframe)] #Choose single sample from each group
			posrows = posrows[posrows %in% datarows]
		}
		negrows = datarows[!datarows %in% posrows]
		testrows = sample_choice(posrows,fraction = testfraction)
		testrows = c(testrows,sample_choice(negrows,fraction = testfraction))
		trainrows = datarows[!datarows %in% testrows]
		if (overfitcheck) testrows = trainrows	# Checking for overfitting
	
		#train the models
		trainingdata = datafile[trainrows,]
		#print(trainingdata)
		if (runNB) modelNB = trainingNB(trainingdata,response)
		if (runLR) modelLR = call_glm(observables,trainingdata,response)
		# Inverse rank-normalisaed
		if (runIRN) modelNBnorm = rank_training(trainingdata,response)
		#modelNBparam = paramtraining(observables,trainingdata,response=response)
		
		# train with transformed data
		if (runARC){
			trainingasinh = asinhdata[trainrows,]
			#print(asinhtraining)
			if (runNB) asinhNB = trainingNB(trainingasinh,response)
			if (runLR) asinhLR = call_glm(observables,trainingasinh,response)
		}
		#print("logistic regression")
		#print(modelLR)
		#print("naive Bayes")
		#print(modelNB)
	
		#prepare test set
		testdata = datafile[testrows,]
		#testdata = replace_na(testdata)
		responsevec = responseframe[testrows] # testdata[response]
		responselist[[length(responselist)+1]] = responsevec

		# find NB posteriors
		if (runNB){
			posteriorsNB = multiple_testing(modelNB,testdata)
			asthmaposteriors = posteriorsNB[,2]
			posteriorslistNB[[length(posteriorslistNB)+1]] = asthmaposteriors
		}
	
	if (runLR){
		posteriorsLR = testing(modelLR,testdata)
		posteriorslistLR[[length(posteriorslistLR)+1]] = posteriorsLR
	}
	
	if (runIRN){
		posteriorsNBnorm = rank_testing(modelNBnorm,testdata)
		posteriorslistNBnorm[[run]] = posteriorsNBnorm
	}
	
	if (runARC){
		#prepare asinh-transformed test set
		testasinh = asinhdata[testrows,]
		#testasinh = replace_na(testasinh)

		# find NB posteriors for asinh-transformed data
		if (runNB){
			asinhposteriorsNB = multiple_testing(asinhNB,testasinh)
			asinhasthmaposteriors = asinhposteriorsNB[,2]
			asinhposteriorslistNB[[run]] = asinhasthmaposteriors
		}
	
		if (runLR){
			asinhposteriorsLR = testing(asinhLR,testasinh)
			asinhposteriorslistLR[[run]] = asinhposteriorsLR
		}
	}
	
	#TODO: Reenable param testing
	#posteriorsNBparam = multiple_testingparams(testdata,priors,response=response,observables)
	#posteriorslistNBparam[[length(posteriorslistNBparam)+1]] = posteriorsNBparam
	
	}
	posteriorslist = list()
	if (runLR){
			posteriorslist[["LR-raw"]] = posteriorslistLR
			if (runARC) posteriorslist[["LR-arc"]] = asinhposteriorslistLR
	}
	if (runNB){
			posteriorslist[["NB-raw"]] = posteriorslistNB
			if (runARC) posteriorslist[["NB-arc"]] = asinhposteriorslistNB
	}
	if (runIRN) posteriorslist[["NB-IRN"]] = posteriorslistNBnorm
#posteriorslist = list("LR-raw" = posteriorslistLR, "NB-raw" = posteriorslistNB, "NB-IRN" = posteriorslistNBnorm, "LR-arc" = asinhposteriorslistLR, "NB-arc" = asinhposteriorslistNB)

	assesslist = call_ROCR(posteriorslist,responselist=responselist,predictors=observables,runs=runs,plotassess=!auconly,plot_title=plot_title,cutoffs=cutoffs,clt=clt)

	assesslist = fpovertp(assesslist)
	#if (!auconly) print(assesslist$auc)
	assesslist
}

#refine_assessment <- function(assesslist,assess1="fpr",assess2="tpr",op)

call_delta_auc <- function(filenames = c("work/data/CAS_Immunology_db.xlsx"),observables = c(),priors = 0,runs = 100,testfraction = 3,response = "crwz5",conditions="",auconly=TRUE){
	# Calls klaR_assessment() twice, via delta_auc(), 
	# w and wo overfitcheck, finds delta AUC
	
	if (length(observables) == 0) observables = readobserves()	
	datafile = readframe(filenames[[1]])	# Read in data
	datafile = preparedata(datafile,observables,response,conditions=conditions)
	delta_auc(datafile,observables = observables,priors = priors,runs = runs,testfraction = testfraction,response = response,conditions=conditions,auconly=auconly)
}

delta_auc <- function(datafile,observables = names(datafile),priors = 0,runs = 100,testfraction = 3,response = "crwz5",conditions="",auconly=TRUE){
	# Calls klaR_assessment() twice, w and wo overfitcheck, finds delta AUC
	#print(observables)
	#responseframe = datafile[response][[1]]	#TODO: Implement responseframe
	#names(responseframe) = rownames(datafile)	#TODO: Implement this
	
	if (length(priors) == 1) priors = findpriors(datafile[response])

	#Find AUC without overfit check
	auc1 = klaR_assessment(datafile,observables=observables,runs = 100,testfraction = 3,response = response,overfitcheck=FALSE,conditions="",auconly=auconly)

	#Find AUC with overfit check
	auc2 = klaR_assessment(datafile,observables=observables,runs = 100,testfraction = 3,response = response,overfitcheck=TRUE,conditions="",auconly=auconly)	

	deltaAUC = c(auc2[["auc-lr"]] - auc1[["auc-lr"]],auc2[["auc-nb"]] - auc1[["auc-nb"]]) 
	#print(deltaAUC)
	#print(auc1)
	#print(auc2)
	deltavec = c(deltaAUC,auc1[["auc-lr"]],auc1[["auc-nb"]],auc2[["auc-lr"]],auc2[["auc-nb"]])
	names(deltavec) = c("LRdeltaAUC","NBdeltaAUC","aucLR1","aucNB1","aucLR2","aucNB2")
	deltavec
}

multisingleAUC <- function(filedata,observables=names(filedata),response="crwz5",runs = 100,testfraction = 3,priors=0,overfitcheck=FALSE,conditions="",auconly=FALSE,datavars=TRUE,plot_title=""){
	# Cycles through every observable and finds the AUC for that observable alone.
	# Returns list/vec of AUCs for each observable
	
	predictors = observables[!observables %in% response]
	predictorlist = list()
	
	for (predict in predictors){
		#cat("DEBUG: predict",predict,'\n')
		cleandata = remove_na(filedata[unique(c(observables,response))])
		auclist = klaR_assessment(filedata,observables=predict,runs=runs,testfraction=testfraction,priors=priors,response=response,overfitcheck=overfitcheck,conditions=conditions,auconly=auconly,datavars=datavars,plot_title=predict,runNB=FALSE,runARC=FALSE)
		predictorlist[[predict]] = auclist[['auc']][['AUC-LR-raw']][['aucmean']]
	}
	predictorlist
}

ascdescAUC <- function(nb_predict=100,filenames = c("work/data/CAS_Immunology_db.xlsx"),response="crwz5"){
	# Lists AUC with increasing numbers of predictors in both ascending
	# and descending order of Shapiro-Wilks index
	
	Wauc = read.table("work/code/distributions/presentations/Wauc.txt")
	predictors = names(Wauc)[1:nb_predict]
	
	# Prepare filedata for klaR_assessment
	filedata = call_preparedata(filenames=filenames,observables=predictors,response=response)
	
	
	# List in order of descending W
	descending = as.data.frame(array(data=0,dim=c(1,2)))
	names(descending) = c('auc-lr','auc-nb')
	for (nb in 1:nb_predict){
		print(predictors[1:nb])
		newdesc = klaR_assessment(filedata,observables=predictors[1:nb],auconly=TRUE)
		descending = rbind(descending,newdesc)
	}
	descending = descending[2:nrow(descending),]
	rownames(descending) = predictors

	# List in order of ascending W
	ascending = as.data.frame(array(data=0,dim=c(1,2)))
	names(ascending) = c('auc-lr','auc-nb')

	for (nb in nb_predict:1){
		print(predictors[nb_predict:nb])
		newasc = klaR_assessment(filedata,observables=predictors[nb_predict:nb],auconly=TRUE)
		ascending = rbind(ascending,newasc)
	}
	ascending = ascending[2:nrow(ascending),]
	#names(ascending) = c('auc-lr','auc-nb')
	rownames(ascending) = rev(predictors)
	
	list('desc' = descending,'asc' = ascending)
}

call_pROC <- function(classvec,posteriorsNB,posteriorsLR){
	#Perform analysis using pROC tools
	
	rocNB = roc(classvec,posteriorsNB,plot=FALSE,ci=TRUE)
	#ciNB = ci.thresholds(rocNB)
	#print(ciNB)
	
	# find LR posteriors
	rocLR = roc(classvec,posteriorsLR,plot=FALSE,ci=TRUE)
	#ciLR = ci.thresholds(rocLR)
	#print(ciLR)
	roctest = roc.test(rocNB,rocLR,method="delong")
	#Zvalues = c(Zvalues,roctest$statistic)
	#pvalues = c(pvalues,roctest$p.value)

	#plot(rocNB,col='blue')
	#plot(ciNB,col='blue')
	#plot(rocLR,col='green',add=TRUE)
	#plot(ciLR,col='green')
	#print(c("Mean Z-value",mean(Zvalues)))
	#print(c("Mean p-value",mean(pvalues)))
	roctest
}

	
interface_training <- function(filenames = c("work/data/CAS_Immunology_db.xlsx"),class = c("crwz5"),samples = NA){	# Train classifier
	#class = c(class)
	trainingdata = preparedata(filenames,class)	# Read in data
	if (!is.na(samples)) trainingdata = trainingdata[samples,]

	model = trainingNB(trainingdata,class)
	model
}

assess <- function(predictions,classvec){	# Prediction accuracy
	predictnames = names(predictions)
	posipred = names(predictions[predictions == 1])	# predicted positive
	negapred = names(predictions[predictions == 0])	# predicted negative
	positives = rownames(classvec)[classvec == 1]		# total positives
	negatives = rownames(classvec)[classvec == 0]		# total negatives
	posipred = posipred[!is.na(posipred)]
	negapred = negapred[!is.na(negapred)]
	
	truepos = positives[positives %in% posipred]		# true positives
	falseneg = positives[positives %in% negapred] 		# false negatives
	falsepos = negatives[negatives %in% posipred] 		# false positives
	trueneg = negatives[negatives %in% negapred]		# true negatives
	
	nb_truepos = length(truepos)
	nb_trueneg = length(trueneg)
	nb_falsepos = length(falsepos)
	nb_falseneg = length(falseneg)
	
	sensitivity = nb_truepos/(length(positives))
	specificity = nb_trueneg/(length(negatives))
	if (nb_falsepos == 0) pvpve = 1		# avoid division by zero
	else pvpve = nb_truepos/(nb_truepos + nb_falsepos)	# pred. value +ve
	if (nb_falseneg == 0) pvnve = 1
	else pvnve = nb_trueneg/(nb_trueneg + nb_falseneg)	# pred. value -ve
	sspp = c(sensitivity,specificity,pvpve,pvnve)
	sspp
}

interface_testing <- function(model,record,filenames = c("work/data/CAS_Immunology_db.xlsx")){	
	#make prediction for single record
	testdata = readframe(filenames[[1]])
	sample = testdata[record,]
	posteriors = testingNB(model,sample)
	#print(posteriors)
}

multiple_testing <- function(model,testdata){	
	#Predict for testing set, testdata
	
	firstloop = TRUE
	for (record in 1:nrow(testdata)){
		posterior = testingNB(model,testdata[record,])
		if (firstloop){
			posteriors = posterior
			firstloop = FALSE
			#print(posteriors)
		}
		else{
			posteriors = rbind(posteriors,posterior)
		}
		
	}
	rownames(posteriors) = rownames(testdata)
	posteriors
}

trainingNB <- function(trainingdata,response = "crwz5"){	
	#Generate model from NaiveBayes
	
	response = c(response)	# Ensure response is a vector
	grouping = as.factor(trainingdata[[response]])	#TODO: insert response after debug
	observables = names(trainingdata)
	observables = observables[!observables == response]
	trainingdata = trainingdata[observables]
	model <- NaiveBayes(trainingdata,grouping,usekernel = FALSE, fL=0)	#, na.action = na.omit)
	model
}

testingNB <- function(model,sample){	
	#Predict sample given model
	
	predictions = predict(model,sample)[2][[1]]	#TODO: Check if handles vector samples
	predictions
}
