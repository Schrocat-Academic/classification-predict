# functions.R
#
# Created by Michael Walker, mid-2013
# Emails: walkerm1@student.unimelb.edu.au, m.walker@aip.org.au
#
# Contains functions of general use which are called, perhaps indirectly, by call_ARTIVA()
#

library(gdata)

#source("work/code/functions2.R")

readframe <- function(filenames=c("work/data/CAS_Immunology_db.xlsx")){  
	# Reads excel file and removes excluded entries
	
	namelength = nchar(filenames[1])
	if (strsplit(filenames[1],"[.]")[[1]][2]=="txt") {
		print("read.table")
		xlsdata <- read.table(filenames[1],sep=" ")
	}
	else {
		print("read.xlsx")
		#print(filenames)
		xlsdata <- read.xls(filenames[1])
	}
	if (length(filenames) > 1){
		filenames <- filenames[2:length(filenames)] 
		for (filename in  filenames){
			if (strsplit(filename,"[.]")[[1]][2]=="txt") 	new.table = as.data.frame(read.table(filename))
			else new.table = read.xls(filename)
			xlsdata <- appendframe(xlsdata,new.table)	# bind all specified data files
		}
	}
	nullstring = "#NULL!"
	if ("EXCL_AGE" %in% names(xlsdata)) {
		excl_records = as.numeric(rownames(subset(xlsdata,EXCL_AGE != nullstring)))+1
		#print(c("Removing",length(excl_records),"excluded entries",excl_records),quote = FALSE)
		as.data.frame(subset(xlsdata,EXCL_AGE == nullstring)) 
	}else {
	#print("No excluded data entries to remove",quote = FALSE)
	as.data.frame(xlsdata)
	}
}

appendframe <- function(olddata,new.table){
	# Appends new.table to frame olddata, keeping only those samples 
	# that are common to both
	newtableIDs = as.vector(unlist(new.table$SUBJ_ID))
	subjids = as.vector(unlist(olddata$SUBJ_ID))
	subjidsINnewtableIDs = subjids %in% newtableIDs
	if (FALSE %in% subjidsINnewtableIDs & 'SUBJ_ID' %in% names(new.table)){
		subjids = subjids[subjidsINnewtableIDs]
		olddata = subset(olddata,SUBJ_ID %in% subjids,select=names(olddata))
		new.table = subset(new.table,SUBJ_ID %in% subjids,select=names(new.table))
	}
	# TODO: Ensure that SUBJ_IDs are in the same order for both tables.
	olddata <- cbind(olddata,new.table)	# bind all specified data files
	olddata
}

is_integer <- function(datacol){
	#Test a column or vector of data to see if it is exclusively integer valued.
	testcol = datacol == round(datacol)
	!FALSE %in% testcol
}

readtest <- function(){  #Reads testnodes.txt and prepares ARTIVA simulation
	parentfolder = ""
	targetfolder = ""
	bothfolder = ""
	conditionfolder = ""
	suboutfolder = ""
	checkCPprior = FALSE
	checkEDGESprior = FALSE
	datafiles = c()
	testfile <- file('/Users/walkerm1/work/code/testnodes.txt','r') 
    fileline <- readLines(testfile,1)
	while (fileline[[1]][1] != 'conditions'){    # read lines until conditions
		fileline <- strsplit(readLines(testfile,1),'\t')
		#print(fileline)
		if (fileline[[1]][1] == 'datafile') datafiles = c(datafiles,fileline[[1]][2])
		else if (fileline[[1]][1] == 'outfolder') {
			outfolder = fileline[[1]][2] 
			if (substring(outfolder,nchar(outfolder)) != "/") outfolder = paste0(outfolder,"/")
		}
		else if (length(fileline[[1]]) > 1 & fileline[[1]][1] == 'suboutfolder') {
			suboutfolder = fileline[[1]][2]
			if (substring(suboutfolder,nchar(suboutfolder)) != "/") suboutfolder = paste0(suboutfolder,"/") # end with '/'
			cat("DEBUG: suboutfolder",suboutfolder,'\n')
		}
		else if (fileline[[1]][1] == 'timesteps') timesteps = as.numeric(fileline[[1]][2])
		else if (fileline[[1]][1] == 'firsttime') firsttime = as.numeric(fileline[[1]][2])
		else if (fileline[[1]][1] == 'parents') parents = fileline[[1]][2:length(fileline[[1]])]
		else if (fileline[[1]][1] == 'targets') targets = fileline[[1]][2:length(fileline[[1]])]
		else if (fileline[[1]][1] == 'alphaCP') {
			alphacp = as.numeric(fileline[[1]][2])
			checkCPprior = TRUE
			}
		else if (fileline[[1]][1] == 'betaCP') betacp = as.numeric(fileline[[1]][2])
		else if (fileline[[1]][1] == 'alphaEdges'){
			alphaedges = as.numeric(fileline[[1]][2])
			checkEDGESprior = TRUE
			}
		else if (fileline[[1]][1] == 'betaEdges') betaedges = as.numeric(fileline[[1]][2])
	}
	parentstring = "parents"
	targetstring = "targets"
	simnodes = sort(unique(c(parents,targets)))
	simstring = "simnodes"
	for (node in simnodes) {
		simstring <- paste(simstring,node,sep='\t')
		#print(c("node",node))
		if (!node %in% targets) {
			parentstring <- paste(parentstring,node,sep='\t') 
			parentfolder = paste(parentfolder,node,sep='-') # if not in targets
		}
		else if (!node %in% parents) {
			targetstring = paste(targetstring,node,sep='\t')
			targetfolder = paste(targetfolder,node,sep='-') # if not in parents
		}
		else {
			targetstring <- paste(targetstring,node,sep='\t')
			parentstring <- paste(parentstring,node,sep='\t') 
			bothfolder = paste(bothfolder,node,sep='-') # both target and parent
		}
	}
	#print(c("bothfolder",bothfolder))
	if (substr(bothfolder,1,1)=="." | substr(bothfolder,1,1)=="-") bothfolder = substr(bothfolder,2,nchar(bothfolder))
	if (substr(parentfolder,1,1)=="." | substr(parentfolder,1,1)=="-") parentfolder = substr(parentfolder,2,nchar(parentfolder))
	if (substr(targetfolder,1,1)=="." | substr(targetfolder,1,1)=="-") targetfolder = substr(targetfolder,2,nchar(targetfolder))
	#print(c("bothfolder",bothfolder))
	simfolder = paste("B",bothfolder,"P",parentfolder,"T",targetfolder,sep="")
	#print("Reading conditions")
	fileline <- readLines(testfile)	# Read conditions on remaining lines
	close(testfile)
	#print(c("fileline",fileline))
	nullstring = "#NULL!"
	conditions = "TRUE"
	conditionlist = sort(fileline)
	conditionfolder = paste0(conditionlist,collapse="-")
	simfolder = paste(simfolder,"C",conditionfolder,sep="") # Contribution to folder name
	#checkCPprior = (TRUE %in% grepl("priors",strplit(suboutfolder,"/"))) 
	#checkEDGESprior = (TRUE %in% grepl("edges",strplit(suboutfolder,"/"))) 
	checkprior = (checkCPprior | checkEDGESprior) 
	if (checkprior) simfolder = paste0(simfolder,'.prior')
	if (checkCPprior) simfolder = paste0(simfolder,'.alphaCP',alphacp,'.betaCP',betacp)
	if (checkEDGESprior) simfolder = paste0(simfolder,'.alphaEDGES',alphaedges,'.betaEDGES',betaedges)
	outfolder = paste0(outfolder,suboutfolder)
    if(! outfolder %in% system("ls" ,intern=TRUE)){	# Add folder if missing
   		syscmnd = paste("mkdir",outfolder)
	   	print(c("Creating folder",outfolder))
		system(syscmnd)
	}
	outfolder = paste0(outfolder,simfolder)
	print(c("outfolder",outfolder))
	#if (length(conditionlist)) {
	for (condition in conditionlist) conditions <- paste(conditions,condition,sep = ' & ')
	conditions = paste(conditions,'TRUE',sep=" & ")
	simoutdir = paste0("/Users/walkerm1/work/data/simulations/AlbertChibb/",suboutfolder)
	simfilename = paste0(simoutdir,simfolder,'.txt')
	#print(c("simfilename",simfilename))
    if(! simoutdir %in% system("ls" ,intern=TRUE)){	# Add folder if missing
    	syscmnd = paste("mkdir",simoutdir)
		system(syscmnd)
    }
	figoutdir = paste0("/Users/walkerm1/work/data/figures/AlbertChibb/",suboutfolder)
	figfilename = paste0(figoutdir,simfolder,'.txt')
	#print(c("simfilename",simfilename))
    if(! figoutdir %in% system("ls" ,intern=TRUE)){	# Add folder if missing
    	syscmnd = paste("mkdir",figoutdir)
		system(syscmnd)
    }

	simfile = file(simfilename,'w')
	writeLines(paste(outfolder,simfolder,sep='\t'),simfile)
	writeLines(paste("outfolder",outfolder,sep='\t'),simfile)
	writeLines(paste("suboutfolder",suboutfolder,sep='\t'),simfile)
	writeLines(paste("simfolder",simfolder,sep='\t'),simfile)
	writeLines(paste("datafiles",datafiles,sep='\t'),simfile)
	writeLines(paste("firsttime",firsttime,sep='\t'),simfile)
	writeLines(paste("timesteps",timesteps,sep='\t'),simfile)
	writeLines(simstring,simfile)
	writeLines(parentstring,simfile)
	writeLines(targetstring,simfile)
	if (checkCPprior) {
		writeLines(paste("alphaCP",alphacp,sep='\t'),simfile)
		writeLines(paste("betaCP",betacp,sep='\t'),simfile)
	}
	if (checkEDGESprior) {
		writeLines(paste("alphaEdges",alphaedges,sep='\t'),simfile)
		writeLines(paste("betaEdges",betaedges,sep='\t'),simfile)
	}
	writeLines(paste("conditions",conditions,sep='\t'),simfile)
	close(simfile)
	simname = paste0(suboutfolder,simfolder)
	simname
}

readsim <- function(simname){	#,simfolder = "/Users/walkerm1/work/data/simulations/AlbertChibb"
	simfile = paste0(simname,".txt")
	testfile <- file(simname,'r') 
    fileline <- readLines(testfile,1)
	while (fileline[[1]][1] != ''){    # read lines until EOF
		fileline <- strsplit(readLines(simfile,1),'\t')
		if (fileline[[1]][1] == 'datafile') datafile = fileline[[1]][2] 
		else if (fileline[[1]][1] == 'outfolder') {
			outfolder = fileline[[1]][2] 
			if (substring(outfolder,nchar(outfolder)) != "/") outfolder = paste0(outfolder,"/")
		}
		else if (fileline[[1]][1] == 'simnodes') outfolder = fileline[[1]][2] 
		else if (fileline[[1]][1] == 'parents') parents = fileline[[1]][2] 
		else if (fileline[[1]][1] == 'targets') targets = fileline[[1]][2] 
		else if (fileline[[1]][1] == 'conditions') conditions = fileline[[1]][2] 
	}
	close(testfile)
    c(datafile,outfolder,simnodes,parents,targets,conditions)
}

call_findvars <- function(filename = c("work/data/CAS_Immunology_db.xlsx"),years=c()){
	# Generates list of variables names given filename
	datafile = readframe(filename[[1]])
	findvars(datafile,years)
}

findvars <- function(filedata,years=c()){
	# Generates list of variables from names of filedata
	returnvars = names(filedata)
	siftvars(returnvars,years)
}

siftvars <- function(returnvars,years=c()){
	# Extract desired variables from allnames. Called by findvars
	
	excluded = c("ID","study","code","EXCL","dob","coll")
	for (excl in excluded){
		returnvars = grep(excl,returnvars,invert=TRUE,value=TRUE,ignore.case=TRUE)
	}
	if (length(years) != 0){
		yearvar = list()
		# Handle birth year (year 0) data
		if (0 %in% years){ 
			yearvar[["0"]] = returnvars[!substring(returnvars,nchar(returnvars),nchar(returnvars)) %in% 1:5]
			years = years[years != 0]
		}
		for (year in years)	yearvar[[as.character(year)]] = returnvars[substring(returnvars,nchar(returnvars),nchar(returnvars)) == year]
		returnvars = unlist(yearvar)
	}
	returnvars
}

timesteps <- function(adddynaroots,xlsnames,steps=c('cd','6m','1','2','3','4','5')){	#6m
	# Returns timesteps available for given factors
	
	newsteps = steps
	dynamic = c()
	for (step in steps){	# Construct vector of timepoint name roots
		cat("DEBUG: newsteps",newsteps,'\n')
		dyna <- paste(adddynaroots,step,sep="")
		for (dynaloop in dyna){
			if (!dynaloop %in% xlsnames){
				newsteps = newsteps[newsteps != step]
				cat("DEBUG: MISSING dynaloop",dynaloop,'\n')
				# Guard against variables with no cd entry but with 1,2x6m entries
				#if (step == 'cd') newsteps = newsteps[newsteps != '6m']
				break
			}
			else dynamic <- c(dynamic,dynaloop)
		}
	}
	list("dynamic" = dynamic, "timesteps" = newsteps)
}

addfnctn <- function(olddata,fieldnames,fnctn){	
	#Apply fnctn to factors in fieldnames
	
	newnodes = readnewnodes(filestring)
	newnames = names(newnodes)
	fnctnnodes = olddata
	fnctnnames = c()
	olddata = cbind(olddata,newnodes)
	for (fieldname in fieldnames){
		fnctnname = paste0(fnctn,fieldname)
		if (!TRUE %in% grepl(fnctnname,names(olddata))){	#avoid recalculation
			grepnames = grep("fieldname",olddata)
			for (grepname in grepnames){
				fngrep = paste0(fnctn,grepname)
				#if (grepname == fieldname) newdata = switch(fnctn,'lg' = log10(olddata[grepname]),'log10' = log10(olddata[grepname]))
				if (nchar(grepname) <= nchar(fieldname)+1 | fieldname == strsplit(grepname," ")[[1]][1]){
					newdata = switch(fnctn,'lg' = log10(olddata[grepname]),'log10' = log10(olddata[grepname]))
					fnctnnodes = cbind(fnctnnodes,newdata)
					fnctnnames = c(fnctnnames,fngrep)
				}else newdata = switch(fnctn,'lg' = log10(addbooldata[grepname]),'log10' = log10(addbooldata[grepname]))
			}
		}
	} 
	newdata
}

minus <- function(xlsdata,node1,node2){	
	# Finds difference between fields
	field1 <- xlsdata[node1]  #reads included data field
	field2 <- xlsdata[node2]  #reads excluded data field
	field1 - field2
	}

plus <- function(xlsdata,node1,node2){	
	# Finds sum of two fields
	field1 <- xlsdata[node1]  #reads included data field
	field2 <- xlsdata[node2]  #reads excluded data field
	field1 + field2
	}

without <- function(xlsdata,node1,node2){	
	# node1 with node2 == 0
	field1 <- xlsdata[node1]  #reads included data field
	field2 <- xlsdata[node2]  #reads excluded data field
	field1 * !field2
	}

with <- function(xlsdata,node1,node2){	
	# node1 with node2 != 0
	field1 <- xlsdata[node1]  #reads included data field
	field2 <- xlsdata[node2]  #reads excluded data field
	field1 * !!field2
	}

exclude <- function(xlsdata,nodeincl,nodeexcl){ 
	# Generates subset of nodeincl which excludes nodeexcl.
	# Intended for infection data with different, overlapping infection types.
	
	fieldexcl <- xlsdata[nodeexcl]  #reads excluded data field
	nodeunion = findunion(nodeincl,nodeexcl)
	fieldunion <- xlsdata[nodeunion]  #reads union data field
	fieldunion - fieldexcl
	#dataincl
}
	
intersect <- function(xlsdata,node1,node2){ 
	#Generates field giving members of fields node1, node2.
	# Intended for infection data with different, overlapping infection types.
	field1 <- xlsdata[node1]  #reads included data field
	field2 <- xlsdata[node2]  #reads excluded data field
	nodeunion = findunion(node1,node2)
	fieldunion <- xlsdata[nodeunion]  #reads union data field
    field1 + field2 - fieldunion  #intersection field
    #dataintersect
}

findunion <- function(name1,name2){	
	#Finds the union of fields name1, name2.
	# Intended for infection data with different, overlapping infection types.
	prefix1 = substring(name1,1,1)	# First prefix for comparison
	suffix = substring(name1,2,nchar(name1))    # Suffix of name1 is union suffix
	prefix2 = substring(name2,1,1)   # Second prefix for comparison
	prefixes = c(prefix1,prefix2)
	unionprefix = ""
	if ('f' %in% prefixes & 'w' %in% prefixes) unionprefix = 's'
	if (unionprefix == "") name1
	else paste0(unionprefix,suffix)
	}
    
isnewboundary <- function(simnode){ # Determines if new field is boundary
	filestring = "work/data/newnodes.txt"
	readfile = file.access(filestring,4)  # is file readable
	if (readfile) newnodes <- read.table(filestring)
	if (readfile & TRUE %in% grepl(simnode,names(newnodes))) isboundary = TRUE
	else {
		nodeparse = strsplit(simnode,"_")
		isboundary = (TRUE %in% grepl(nodeparse[1],names(xlsdata)) & TRUE %in% grepl(nodeparse[3],names(xlsdata))) # only two boundaries make a boundary
	}
	isboundary
}

makebooldata <- function(testframe,simnode){  
	#Adds datafields based on exclusion/intersection and _plus_, _minus_ to testframe
	
	testnames = names(testframe)
	
	# Extract nodes to be combined, nodeparse[1], nodeparse[3]
	nodeparse = strsplit(simnode,"_")[[1]]
	#TODO: find empty strings in nodeparse indicating trailing '-' in nodes
	if (substr(simnode,nchar(simnode),nchar(simnode))=='_') nodeparse[3] = paste0(nodeparse[3],'_')	# in case nodeparse[3] ends in '_'
	operation = nodeparse[2]
	while (substr(nodeparse[3],nchar(nodeparse[3]),nchar(nodeparse[3]))==" "){
		nodeparse[3] = substr(nodeparse[3],1,nchar(nodeparse[3])-1)	# in case nodeparse[3] ends in ' '
	}
	nodes = c(nodeparse[1],c(nodeparse[3]))
	nodelist1 = c()
	nodelist2 = c()
	datanames = c()
	
	# Combine two boundary nodes
	if (!FALSE %in% (nodes %in% testnames)){
		nodelist1 = c(nodes[1])
		nodelist2 = c(nodes[2])
		datanames = c(simnode)	 # list field names
		seedfields = nodes[1]
	} 
	# Combine two dynamic nodes
	else if (!TRUE %in% (nodes %in% testnames)){
		grepout = grep(nodes[1],testnames,value=TRUE)
		for (nodes1 in grepout){	# For datafields of node[1]
			if (nchar(nodes1) == nchar(nodes[1])+1 | nodes[1] == strsplit(nodes1,'6m')[[1]][1] | nodes[1] == strsplit(nodes1," ")[[1]][1]){  #avoid unwanted suffices
				
				# Identify time-period suffix
				# TODO: Handles only years and '6m'. Handle '24m'
				suffix = substr(nodes1,nchar(nodes1),nchar(nodes1))
				if (suffix == 'm') suffix = '6m'
				
				#Corresponding field of node[2]
				nodes2 = paste0(nodes[2],suffix)
				if (!nodes2 %in% testnames) nodes2 = paste(nodes[2],suffix,sep=" ")
				if (nodes2 %in% testnames){	# if both values at timeslice
					nodelist1 = c(nodelist1,nodes1)	# list fields to add
					nodelist2 = c(nodelist2,nodes2)
					datanames = c(datanames,paste(simnode,suffix,sep=" ")) # list field names
					#print(c("DEBUG: datanames",datanames))
				}
			}
		}
		seedfields = nodelist1	# seed datafields
	} 
	#nodes[1] only in testnames
	else if (!nodes[2] %in% testnames){	
		grepout = grep(nodes[2],testnames,value=TRUE)
		for (nodes2 in grepout){	# For datafields of node[1]
			if (nchar(nodes2) == nchar(nodes[2])+1 | nodes[2] == strsplit(nodes2,'6m')[[1]][1] | nodes[2] == strsplit(nodes2," ")[[1]][1]){  #avoid unwanted suffices
				
				# Identify time-period suffix
				# TODO: Handles only years and '6m'. Handle '24m'
				suffix = substr(nodes2,nchar(nodes2),nchar(nodes2))
				if (suffix == 'm') suffix = '6m'
				
				#Corresponding field of node[2]
				nodelist1 = c(nodelist1,nodes[1])	# list fields to add
				nodelist2 = c(nodelist2,nodes2)
				datanames = c(datanames,paste(simnode,suffix,sep=" ")) # list field names
			}
		}		
		seedfields = nodes[1]	# seed datafields
	} 
	#nodes[2] only in testnames
	else{		# !nodes[1] %in% testnames	
		grepout = grep(nodes[1],testnames,value=TRUE)
		for (nodes1 in grepout){	# For datafields of node[1]
			if (nchar(nodes1) == nchar(nodes[1])+1 | nodes[1] == strsplit(nodes1,'6m')[[1]][1] | nodes[1] == strsplit(nodes1," ")[[1]][1]){  #avoid unwanted suffices
				
				# Identify time-period suffix
				# TODO: Handles only years and '6m'. Handle '24m'
				suffix = substr(nodes1,nchar(nodes1),nchar(nodes1))
				if (suffix == 'm') suffix = '6m'
				
				#Corresponding field of node[1]
				nodelist1 = c(nodelist1,nodes1)	# list fields to add
				nodelist2 = c(nodelist2,nodes[2])
				datanames = c(datanames,paste(simnode,suffix,sep=" ")) # list field names
			}
		}		
		seedfields = nodes[2]    # seed datafields
	}
	datafields = testframe[seedfields]
	fieldnames <- c(names(datafields),datanames)
	for (i in 1:length(nodelist1)){	
		datafield = switch(operation,'and' = intersect(testframe,nodelist1[i],nodelist2[i]),'not' = exclude(testframe,nodelist1[i],nodelist2[i]),'without' = without(testframe,nodelist1[i],nodelist2[i]),'with' = with(testframe,nodelist1[i],nodelist2[i]),'minus' = minus(testframe,nodelist1[i],nodelist2[i]),'plus' = plus(testframe,nodelist1[i],nodelist2[i]))
		datafields <- cbind(datafields,datafield)
	}
	names(datafields) = fieldnames  # Assign field names to frame datafields
	datafields <- datafields[datanames]
	#print(fieldnames)
	#print(length(rownames(datafields)))
	datafields
}

appendnewdata <- function(oldframe,fieldnames,filestring = "work/data/newnodes.txt"){
	# Appends after generating, if necessary, new variables onto oldframe.
	# Newly generated variables are stored in newnodes.txt
	
	fileexists = file.exists(filestring)
	if (fileexists)	newnodes <- readnewnodes(filestring)
	else{
		newnodes = makebooldata(oldframe,fieldnames[1])	# Create newnodes
		if (length(fieldnames)>1) fieldnames <- fieldnames[2:length(fieldnames)]	# don't double count
		#print(names(newnodes))
	}
	if (length(fieldnames)>0){	# in case one fieldname and no newnodes
		newnames = names(newnodes)
		for (name in fieldnames){
			if (name %in% newnames){
				oldframe <-cbind(oldframe,newnodes[name])
				next	#foundname = TRUE
			}
			else{ 
				foundname = FALSE
				grepout = grep(name,newnames,value=TRUE)
				if (TRUE %in% grepl(name,substr(newnames,1,nchar(name)))){
					for (grepname in grepout) if (nchar(grepname) <= nchar(name)+2 | name == strsplit(grepname," ")[[1]][1]){
							oldframe <-cbind(oldframe,newnodes[grepname]) #newnodes dynamic data   
							foundname = TRUE
					}
				}
			}
			if (!foundname){
				boolframe <- cbind(oldframe,newnodes)
				newnode = makebooldata(boolframe,name)
				#print(c("newnode from makebooldata",newnode))
				newnodes <- cbind(newnodes,newnode)
				oldframe <- appendframe(oldframe,newnode)
			}
		}
	}
	print("Writing table newnodes")
	write.table(newnodes,filestring)		# save newnodes
	oldframe
}

renamenewnode <- function(oldname,newname){	
	lenold = nchar(oldname)
	filestring = "work/data/newnodes.txt"
	newnodes = readnewnodes(filestring)
	newnames = names(newnodes)
	namepos = grep(oldname,newnames)
	if (length(namepos)==0) print(c(oldname," not found in newnames"),quote=FALSE)
	else for (ind in namepos) {
			old = newnames[[ind]]
			if (nchar(old) == lenold +1|oldname == strsplit(old," ")[[1]][1]){	#Avoid unwanted suffices
				new = paste0(newname,substr(old,nchar(old),nchar(old)))
				newnames[[ind]] = new
			} else if (oldname == old) newnames[[ind]] = newname
			names(newnodes) = newnames
			write.table(newnodes,filestring)
		}
	#print(names(newnodes))
}

readnewnames <- function(filestring="work/data/newnodes.txt") names(readnewnodes(filestring))

readnewnodes <- function(filestring="work/data/newnodes.txt"){	
	# Reads in text table with no dots in headers
	newnodes = read.table(filestring)
	newnames = names(newnodes)
	hasdots = grep("\\.",newnames)	# Which headers did R replace ' ' with '.'
	for (ind in hasdots){
		dotsplit = strsplit(newnames[ind],"\\.")[[1]]	# Replace '.'
		newnames[ind] = paste(dotsplit[1],dotsplit[2])	# with ' '
	}
	names(newnodes) = newnames
	newnodes
}

removenewnodes <- function(fieldnames,filestring="work/data/newnodes.txt"){	
	# Removes unwanted fields from newnodes
	
	newnodes = readnewnodes(filestring)
	newnames = names(newnodes)
	for (fieldname in fieldnames){
		keptnames = c()	# Fields to retain				
		for (name in newnames){
			keepname = TRUE
			if (TRUE %in% grepl(fieldname,name)) if (fieldname == name | nchar(name) == nchar(fieldname)+1 | fieldname == strsplit(name," ")[[1]][1]) keepname = FALSE
			if (keepname) keptnames = c(keptnames,name)
		}
		newnames = keptnames
	}
	newnodes = newnodes[newnames]
	write.table(newnodes,filestring)
}

find_binaries <- function(filedata){
	# Identifies which variables in data frame filedata have binary 0,1 data
	# Returns logical vector
	binaries = c()
	for (observe in names(filedata)){
		binaryvals = (filedata[observe] == 0) | (filedata[observe] == 1) | (is.na(filedata[observe]))
		notbinary = FALSE %in% binaryvals
		binaries = c(binaries,!notbinary)
	}
	names(binaries) = names(filedata)
	binaries
}

find_integers <- function(filedata){
	# Identifies which variables in data frame filedata have integer data
	# Returns logical vector
	integers = c()
	for (observe in names(filedata)){
		integervals = (filedata[observe] == round(filedata[observe])) | (is.na(filedata[observe]))
		notinteger = FALSE %in% integervals
		integers = c(integers,!notinteger)
	}
	names(integers) = names(filedata)
	integers
}

declusterize <- function(initdatarows){
	# Generates random vector of datarows so that each has a unique value and all
	#original values are represented.
	
	selections = c()
	initrownames = rownames(initdatarows)
	initdatarows = as.character(unlist(initdatarows[[1]]))
	uniquedata = unique(initdatarows)
	for (entry in uniquedata){
		indices = grep(entry,initdatarows,fixed=TRUE)
		if (length(indices)==1) select = indices[1]
		else select = sample(x=indices,size=1)
		selections = c(selections,select)	# initrownames[]
	}
	selections
}
