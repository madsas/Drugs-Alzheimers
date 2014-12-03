#read the data in#################### 
if (!exists('data')) {
	data <- read.csv("altmann06122013.csv", as.is=T)
	data.df <- data.frame(data)
}
#Get drug columns
dc <- grep("DRUG",colnames(data))
#Get NACCID column
n <- grep("NACCID",colnames(data))

#SWITCH####################
#switch <- 0
switch <- 1

if (!switch) {
	#Take all drugs, regardless of visits and make histogram####################
	all.drugs <- data[,dc]
	all.drugs.l <- as.list(all.drugs)
	all.drugs.l.flat <- unlist(all.drugs.l)
	#remove empty entries 
	all.drugs.l.flat.clean <- all.drugs.l.flat[as.logical(lapply(x,function(x) x!=""))]
	#make frequency table
	all.drugs.table <- sort(table(all.drugs.l.flat.clean))
	#print most and least common
	tail(all.drugs.table)
	head(all.drugs.table)
	#output to csv
	write.table(all.drugs.table,file="sorted_all_drugs_table.csv",sep=",",col.names=NA,qmethod="double")
}

if (switch) {
	#Now repeat but count subjects on drugs insteado of simply drug occurence####################
	id.drugs <- data[,c(n,dc)]
	#get list of ids and make empty data frame
	ids<-unique(id.drugs$NACCID)
	unique.drugs <- rep(list(c()),length(ids))
	names(unique.drugs) <- ids
	#parameters
	d.cols <- ncol(id.drugs)
	#collect drugs from each row
	for (i in 1:nrow(id.drugs)) {
		x <- id.drugs[i,]
		idx <- x$NACCI
		unique.drugs[[idx]] <- c(unique.drugs[[idx]],unlist(c(unname(x[2:d.cols]))))
		#get rid of repeats and blanks
		unique.drugs[[idx]] <- unique(unique.drugs[[idx]][unique.drugs[[idx]] != ""])
	}
	#flatten list and make table
	unique.drugs.table <- sort(table(unlist(unname(unique.drugs))))
	#print most and least common
	tail(unique.drugs.table)
	head(unique.drugs.table)
	#output to csv
	write.table(unique.drugs.table,file="sorted_unique_drugs_table.csv",sep=",",col.names=NA,qmethod="double")
}

