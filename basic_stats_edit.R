### scripts for producing basic stats###

fname <- "altmann06122013.csv"

mydata <- read.csv(fname)


## number of different subjects ##
nsub <- length(table(mydata$NACCID))

## median num. of visits ##
med.vis <- median(table(mydata$NACCID))
barplot(table(table(mydata$NACCID)))

## APOE info ##
#generate map NACCID 2 APOE
naccid2first <- sapply(names(table(mydata$NACCID)), function(x){
  idx <- mydata$NACCID == x
  return(which(idx)[1])
})



naccid2apoe <- mydata[naccid2first,"NACCAPOE"]

table(naccid2apoe)
## 9 - unknown	- 2786
## 1 - e3/e3	- 4863
## 2 - e3/e4	- 2400
## 3 - e3/e2	-  976
## 4 - e4/e4	-  358
## 5 - e4/e2	-  233
## 6 - e2/e2	-   38
#-----------------------
#		11,654

convertAPOE <- function(x){
  if (x == 1)
    return(33)
  if (x == 2)
    return(34)
  if (x == 3)
    return(23)
  if (x == 4)
    return(44)
  if (x == 5)
    return(24)
  if (x == 6)
    return(22)

  return(NA)
}

myapoe <- sapply(mydata$NACCAPOE,convertAPOE)

##convert visit datas to Date datatype
visit.dates <- apply(mydata,1,function(x){
  as.Date(paste(x["visityr"],x["visitmo"], x["visitday"], sep="-"))
})

##convert birth days to Date datatype
birth.dates <- apply(mydata,1,function(x){
  as.Date(paste(x["BIRTHYR"],x["BIRTHMO"], 15, sep="-"))  
})

##age @ visit
age.at.visit <- round((visit.dates - birth.dates) / 365.25,1)

##Disease Status at enrolment

naccid2normcog   <- mydata[naccid2first,"NORMCOG"]
naccid2demented  <- mydata[naccid2first,"DEMENTED"]
naccid2mciamem   <- mydata[naccid2first,"MCIAMEM"]
naccid2mciaplus  <- mydata[naccid2first,"MCIAPLUS"]
naccid2mcinon1   <- mydata[naccid2first,"MCINON1"]
naccid2mcinon2   <- mydata[naccid2first,"MCINON2"]
naccid2impnomci  <- mydata[naccid2first,"IMPNOMCI"]

demrating <- apply(mydata[,c("NORMCOG","IMPNOMCI","MCINON1","MCINON2","MCIAMEM","MCIAPLUS","DEMENTED","PROBAD","PROBADIF","POSSAD","POSSADIF")], 1, function(x){
  if(x[1] == 1)
    return(1)
  if(x[2] == 1)
    return(2)
  if(x[3] == 1 | x[4] == 1)
    return(3)
  if(x[5] == 1 | x[6] == 1)
    return(4)
  #PRIMARY AD DEMENTIA
  if(x[7] == 1 & ( (x[8] == 1 & x[9] == 1) | x[10] == 1 & x[11] == 1) )
    return(5)
  #PRIMARY AD DEMENTIA
  if(x[7] == 1 &  (x[8] == 1 | x[10] == 1 ))
    return(6)
  #NON (PRIMARY) AD DEMENTIA
  if (x[7] == 1)
    return(7)
  return(NA)
})

#HC at enrolment: 7197
#AD at enrolment: 0
##--amnestic MCI
#MCIA at enrolment: 1523
#MCIA+ at enrolment: 1391

#MCI non-a: 468
#MCI non-a2: 250

#impaired no MCI: 825

#sum: 11,654

##summary code
#1 - HC
#2 - impaired no MCI
#3 - non-amnestic MCI
#4 - amnestic MCI
#5 - primary AD dementia
#6 - contributing AD
#7 - non-AD dementia


extractConversionData <- function(x){

  ids <- unique(x$NACCID)

  ppp <- t(sapply(ids, function(y){
    idx <- x$NACCID == y
    
    #first visit (i.e., the one with minimal age)
    ttt <- x[idx,]
    maxage <- max(ttt$age.at.visit)
    first.v <- which.min(ttt$age.at.visit)
    age <- ttt[first.v, "age.at.visit"]
    edu <- ttt[first.v, "EDUC"]
    gender <- ttt[first.v, "SEX"]
    mmse <- ttt[first.v, "MMSE"]
    apoe   <- ttt[first.v, "myapoe"]
    irating <- ttt[first.v, "demrating"]
    nvis   <- nrow(ttt)
    #event <- 0
    #event.time <- age
    orating <- irating
    race <- ttt[first.v, "RACE"]
    hisp <- ttt[first.v, "HISPANIC"] 

    ## non amnestic MCI
    event.naMCI <- 0
    event.naMCI.time <- max(ttt[,"age.at.visit"])
    idx2 <- ttt[,"demrating"] > irating & ttt[,"demrating"] == 3
    if (sum(idx2) > 0){
      event.naMCI <- 1
      event.id <- which.min(ttt[idx2,"age.at.visit"])
      event.naMCI.time <- ttt[idx2,][event.id,"age.at.visit"]
      orating    <- ttt[idx2,][event.id,"demrating"]
    }

    ## amnestic MCI
    event.aMCI <- 0
    event.aMCI.time <- max(ttt[,"age.at.visit"])
    idx2 <- ttt[,"demrating"] > irating & ttt[,"demrating"] == 4
    if (sum(idx2) > 0){
      event.aMCI <- 1
      event.id <- which.min(ttt[idx2,"age.at.visit"])
      event.aMCI.time <- ttt[idx2,][event.id,"age.at.visit"]
      orating    <- ttt[idx2,][event.id,"demrating"]
    }

    ## primary AD - dementia
    event.pAD <- 0
    event.pAD.time <- max(ttt[,"age.at.visit"])
    idx2 <- ttt[,"demrating"] > irating & ttt[,"demrating"] == 5
    if (sum(idx2) > 0){
      event.pAD <- 1
      event.id <- which.min(ttt[idx2,"age.at.visit"])
      event.pAD.time <- ttt[idx2,][event.id,"age.at.visit"]
      orating    <- ttt[idx2,][event.id,"demrating"]
    }

    ## contributing AD - dementia
    event.cAD <- 0
    event.cAD.time <- max(ttt[,"age.at.visit"])
    idx2 <- ttt[,"demrating"] > irating & ttt[,"demrating"] == 6
    if (sum(idx2) > 0){
      event.cAD <- 1
      event.id <- which.min(ttt[idx2,"age.at.visit"])
      event.cAD.time <- ttt[idx2,][event.id,"age.at.visit"]
      #orating    <- ttt[idx2,][event.id,"demrating"]
    }

    ## primary AD - dementia
    event.nonADdem <- 0
    event.nonADdem.time <- max(ttt[,"age.at.visit"])
    idx2 <- ttt[,"demrating"] > irating & ttt[,"demrating"] == 7
    if (sum(idx2) > 0){
      event.nonADdem <- 1
      event.id <- which.min(ttt[idx2,"age.at.visit"])
      event.nonADdem.time <- ttt[idx2,][event.id,"age.at.visit"]
      #orating    <- ttt[idx2,][event.id,"demrating"]
    }

    event      <- max(event.naMCI, event.aMCI, event.pAD)
    event.time <- min(event.naMCI.time, event.aMCI.time, event.pAD.time)

    ##Counts the number of number of drugs for each patient
    first_row <- ttt[1,]
    drugs <- first_row[c(26:65)] #Indices of DRUG1 through DRUG40
    
    num_drugs <- 0

    for (i in 1:40) {
    drug <- drugs[i]
        if (drug != "") {
            num_drugs <- num_drugs + 1
        }
    }

    return(c(y, nvis, age, maxage, gender,race, hisp, apoe, edu, mmse, irating, event, event.time, event.naMCI, event.naMCI.time, event.aMCI, event.aMCI.time, event.pAD, event.pAD.time, event.cAD, event.cAD.time, event.nonADdem, event.nonADdem.time, orating, num_drugs))   
  }))
  colnames(ppp) <- c("ID","NVIS", "AGE","MAXAGE","SEX","RACE","HISP","APOE","EDU","MMSE","demA","event","event.time","conv.naMCI","naMCI.time","conv.aMCI","aMCI.time","conv.pAD","pAD.time","conv.cAD", "cAD.time","conv.nonADdem","nonADdem.time", "demB", "NUMDRUGS")
  return(ppp)
}


#tmpdata <- subset(cbind(mydata, age.at.visit, demrating, myapoe), subset=(NACCAPOE==1 | NACCAPOE==2))
tmpdata <- cbind(mydata, age.at.visit, demrating, myapoe)


all.conversion <- extractConversionData(tmpdata)
all.conversion <- data.frame(all.conversion)

ofname="xxx.RData"
save(all.conversion, file=ofname)
