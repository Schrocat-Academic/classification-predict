# functions2.R
#
# Created by Michael Walker, 4/6/15
# Emails: walkerm1@student.unimelb.edu.au, m.walker@aip.org.au
#
# Contains functions of general use which are not called by call_ARTIVA()
#

padrows <- function(thinframe,fatframe,variables = names(thinframe),ref_var = c("SUBJ_ID","subjid")){
	# Pads out thin frame so it has all the same rows as fatframe.
	# The variables used to match are in ref_var
	
	reffat = ref_var[ref_var %in% names(fatframe)]
	refthin = ref_var[ref_var %in% names(thinframe)]
	cat("DEBUG:",refthin,'\n',variables,'\n')
	variables = variables[!variables %in% refthin]

	outputframe = fatframe[reffat]
	for (var in variables){
		varvec = thinframe[[var]]
		names(varvec) = thinframe[[refthin]]
		varcol = varvec[fatframe[[reffat]]]
		cat("DEBUG: varcol",varcol,'\n')
		outputframe = cbind(outputframe,varcol)
	}
	print(outputframe)
	outputframe = outputframe[,2:ncol(outputframe)]
	names(outputframe) = variables
	outputframe
}

cv_average <- function(y_val,x_val,cutoffs=c(0.35,0.5,0.75)){
	# Averages over vectors y_val, parameterised by x_val,
	#evaluated at cutoffs
	# y_val, x_val are lists of vectors
	
	outputvec = c()
	for (cutoff in cutoffs){
		cutoffvec = c()
		for (run in 1:length(y_val)){
			x_valrun = as.vector(unlist(x_val[run]))
			y_valrun = y_val[[run]]
			y_valrun[is.na(y_valrun)] = 0
			# Elements of x_valrun in descending order
			lwrind = length(x_valrun[x_valrun>=cutoff])
			upprind = lwrind+1
			if (lwrind==length(x_valrun)) yweight = y_valrun[length(x_valrun)]
			else if (lwrind < 2) yweight = y_valrun[1]
			else {
				upprweight = x_valrun[lwrind] - cutoff
				lwrweight = cutoff - x_valrun[upprind]
				yweight = (y_valrun[upprind]*upprweight + y_valrun[lwrind]*lwrweight)/(lwrweight + upprweight)
			}
			cutoffvec = c(cutoffvec,yweight)
		}
		outputvec = c(outputvec,mean(cutoffvec))
	}
	names(outputvec) = cutoffs
	outputvec
}

remove_na_OLD <- function(rawdata){	
	# Remove missing data and NA from rawdata
	
	nadata = is.na(rawdata)
	notRows = rowSums(nadata)
	cleandata = rawdata[notRows==0,]
	cleandata
}

remove_na <- function(rawdata){
	# Remove records with missing data
	
	# Select rows with missing data to remove
	removerows = apply(rawdata,1,anyNA) 
	# Keep other rows (records)
	keeprows = !removerows
		
	# Remove records with missing data
	rawdata[keeprows,]
}
	
replace_na <- function(rawdata){ 
	# Replace NA with variable mean
	
	no_nadata = remove_na(rawdata)		# Identify missing data
	cleandata = rawdata
	#print(rawdata)
	colMeanData = colMeans(rawdata,na.rm=TRUE)
	#print(colMeanData)
	for (var in names(colMeanData)) {
		narows = is.na(rawdata[var])
		isinteger = is_integer(no_nadata[,var])
		replaceval = colMeanData[var]
		if (isinteger) replaceval = round(replaceval)
		cleandata[narows==TRUE,var] = replaceval
	}
	#print(cleandata)
	cleandata
}
