######Add num_drugs
ppp <- t(sapply(ids, function(y) {
	message(y)
    idx <- data$NACCID == y
    ttt <- data[idx,]
    first_row <- ttt[1,]
    drugs <- first_row[c(26:65)]
    num_drugs <- 0
    for (i in 1:40) {
    drug <- drugs[i]
       if (drug != "") {
          num_drugs <- num_drugs + 1
       }
    }
    return (c(y, num_drugs))
}))

colnames(ppp) <- c("NACCID", "NUMDRUGS")

#NUMDRUGS <- ppp[,2]
#NUMDRUGS <- as.numeric(NUMDRUGS)
NUMDRUGS <- c()
for (row_name in rownames(mydata.use)) {
    row <- grep(row_name, rownames(ppp))
    if (length(row) != 0) {
        n_d <- as.numeric(ppp[row,2])
        NUMDRUGS <- c(NUMDRUGS, n_d)
    }
}
