#read the data in 
if (!exists('data')) {
	data <- read.csv("altmann06122013.csv", as.is=T)
	data.df <- data.frame(data)
}

#Get drug columns
dc <- grep("DRUG",colnames(data))

#Take all drugs, regardless of visits and make histogram
all.drugs <- data[,dc]
all.drugs.l <- as.list(all.drugs)
#all.drugs.l.flat <- as.list(unlist(all.drugs.l))
all.drugs.l.flat <- unlist(all.drugs.l)
#remove empty entries (takes a lot of time for reason)
all.drugs.l.flat.clean <- all.drugs.l.flat[as.logical(lapply(x,function(x) x!=""))]
#make frequency table
all.drugs.table <- table(all.drugs.l.flat.clean)
#print most and least common
tail(sort(all.drugs.table))
head(sort(all.drugs.table))

if (FALSE) {
#make one row per subject
##get row number of first instance of each id
ids <- unique(data$NACCID)
rows <- c()
for (i in ids) {
	rows <- c(rows,which(data$NACCID==i)[1])
}

##make dataframe with just drugs and first instance of id


#get drug column numbers and put in new df
dc <- grep("DRUG",colnames(data))
}
