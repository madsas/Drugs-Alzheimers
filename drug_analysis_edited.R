require(survival)

#read.data

data <- read.csv("altmann06122013.csv", as.is=T)
data.df <- data.frame(data)

#prepared data
xtab.fname <- "xxx.RData"
#xtab.fname <- "/mnt/mapricot/musk2/home/altmann/data/NACC/RData/drug_conversion.RData"

#Read in each of the drugs#
drug_list <- read.csv("sorted_unique_drugs_table.csv", header = TRUE, sep = ",")
drug_list <- drug_list$Drug
drug_list <- as.character(drug_list)
drug_list <- t(drug_list)
drug_list <- as.vector(drug_list)
drug_list <- drug_list[1489:1584]

for (cur_drug in drug_list) {

drugA <- "LANSOPRAZOLE"
drugB <- "DOXAZOSIN"


drugA <- cur_drug
#out <- capture.output(drugB)
cat(drugA, file = "drug_coxph_summ.txt", sep="\n", append=TRUE)

drugs <- c(drugA)

#drug columns
dc <- grep("DRUG",colnames(data))
#zzz <- apply(data[,dc],2,function(x){
#  idx <- grep(drug, x)
#  return(data[idx,"NACCID"])
#})
#drug.users <- unique(unlist(zzz))

naccapoe2apoe <- function(x){
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

resolveDX <- function(bbb){

  norm <- bbb["NORMCOG"]
  if (norm == 1)
    return(1)

  mci <- sum(bbb[c("MCIAMEM","MCIAPLUS","MCINON1","MCINON2")] ==1)
  if(mci > 0)
    return(2)

  ad <- sum(bbb[c("PROBAD","POSSAD")] == 1)
  #if we want to distinguish between primary and contributing check variables:
  #"PROBADIF","POSSADIF"
  if(ad > 0)
    return(3)

  #other impairment non-MCI/AD
  return(0)

}

getInfo <- function(dd, id){
  #print(id)
  #baseline info
  myinfo <- c()

  tmp  <- subset(dd, subset=NACCID==id)
  tmp2 <- subset(tmp, subset=vnumber==1)

  #SEX
  myinfo <- c(myinfo, tmp2[1,"SEX"])
  #MMSE
  myinfo <- c(myinfo, tmp2[1,"MMSE"])
  #EDUC
  myinfo <- c(myinfo, tmp2[1,"EDUC"])
  #APOE
  ncapoe <- tmp2[1,"NACCAPOE"]
  myinfo <- c(myinfo, naccapoe2apoe(ncapoe))
  #entry data (Y-M-D)
  edate <- as.Date(paste(tmp2[1,c("visityr","visitmo","visitday")], collapse="-"))  
  myinfo <- c(myinfo, edate)
  #birth day
  xxx <- paste(paste(tmp2[1,c("BIRTHYR","BIRTHMO")], collapse="-"),15,sep="-")
  bdate <- as.Date(xxx)
  myinfo <- c(myinfo, bdate)

  #Baseline Dx
  base.dx <- resolveDX(tmp2)
  myinfo <- c(myinfo, base.dx)

  #number of visits
  nvis <- max(tmp[,"vnumber"])
  #time covered
  tmp3 <- subset(tmp, subset=vnumber==nvis)
  ldate <- as.Date(paste(tmp3[1,c("visityr","visitmo","visitday")], collapse="-"))
  myinfo <- c(myinfo, nvis, ldate)
  
  all.dx <- apply(tmp, 1, resolveDX)
  ##first switch to MCI (if more than baseline)
  idx.mci <- all.dx == 2 & base.dx < 2
  vis.mci <- NA
  time.mci <- ldate
  conv.mci <- sum(idx.mci) > 0
  if (conv.mci){
    vis.mci <- min(tmp[idx.mci,"vnumber"])
    tmpx <- subset(tmp,subset=vnumber==vis.mci)
    time.mci <- as.Date(paste(tmpx[1,c("visityr","visitmo","visitday")], collapse="-"))
  }
  myinfo <- c(myinfo, conv.mci, vis.mci, time.mci)

  ##firest switch to AD (if more than baseline)
  idx.ad <- all.dx == 3 & base.dx < 3
  vis.ad <- NA
  time.ad <- ldate
  conv.ad <- sum(idx.ad) > 0
  if (conv.ad){
    vis.ad <- min(tmp[idx.ad,"vnumber"])
    tmpx <- subset(tmp,subset=vnumber==vis.ad)
    time.ad <- as.Date(paste(tmpx[1,c("visityr","visitmo","visitday")], collapse="-"))
  }
  myinfo <- c(myinfo, conv.ad, vis.ad, time.ad)

  return(myinfo)

}

RXlength <- function(dat, id, dc, dn){
  #message(id)

  rxlen <- 0

  tmp <- subset(dat, subset=NACCID==id)
  use.idx <- apply(tmp[,dc],1,function(x){ length(grep(dn, x)) > 0})
  vmax <- max(tmp[,"vnumber"])
  for(z in which(use.idx)){
    vn <- tmp[z,"vnumber"]

    vnm1 <- max(subset(tmp, subset=vnumber<vn)[,"vnumber"])
    vnp1 <- min(subset(tmp, subset=vnumber>vn)[,"vnumber"])

    #middle visit
    if (vn > 1 && vn < vmax){
      #dist (vn n+1 - vn-1)/2
      k1 <- subset(tmp, subset=vnumber==vnm1)
      vn1 <- as.Date(paste(k1[1,c("visityr","visitmo","visitday")], collapse="-"))
      k2 <- subset(tmp, subset=vnumber==vnp1)
      vn2 <- as.Date(paste(k2[1,c("visityr","visitmo","visitday")], collapse="-"))
      rxlen <- rxlen + (vn2 - vn1)/2
    }

    #first visit
    if (vn == 1){
      #dist vn 2 - vn 1
      k1 <- subset(tmp, subset=vnumber==1)
      vn1 <- as.Date(paste(k1[1,c("visityr","visitmo","visitday")], collapse="-"))
      k2 <- subset(tmp, subset=vnumber==vnp1)
      vn2 <- as.Date(paste(k2[1,c("visityr","visitmo","visitday")], collapse="-"))
      rxlen <- rxlen + (vn2 - vn1)
    }

    #last visit
    if (vn == vmax){
      #dist vn N - vn N-1
      k1 <- subset(tmp, subset=vnumber==vnm1)
      vn1 <- as.Date(paste(k1[1,c("visityr","visitmo","visitday")], collapse="-"))
      k2 <- subset(tmp, subset=vnumber==vmax)
      vn2 <- as.Date(paste(k2[1,c("visityr","visitmo","visitday")], collapse="-"))
      rxlen <- rxlen + (vn2 - vn1)
    }
  }
  return(rxlen)

}

cnames <- c("SEX","MMSE","EDU","APOE","ENTRY","BIRTH","DX","NVIS","LAST","CONVMCI","VISMCI","TIMEMCI","CONVAD","VISAD","TIMEAD")

ids <- sort(unique(data[,"NACCID"]))

if (!file.exists(xtab.fname)){
  message("preparing table... (takes a while)")
  xtab <- t(sapply(ids,function(x){ getInfo(data,x)}))
  colnames(xtab) <- cnames
  save(xtab, file=xtab.fname)
} else {
  message("loading table from: ", xtab.fname)
  load(xtab.fname)
  #xtab <- readRDS("xxx.RData") #Added temporarily for specific file
}


message("computing drug users and usage duration")
#create Drug indicators
drugcol <- c()
drugdur <- c()
for(drug in drugs){
  message(drug)
  zzz <- apply(data[,dc],2,function(x){
    idx <- grep(drug, x)
    return(data[idx,"NACCID"])
  })
  drug.users <- unique(unlist(zzz))

  tmp <- rep(0, nrow(xtab))
  names(tmp) <- rownames(xtab)
  tmp[drug.users] <- 1
  drugcol <- cbind(drugcol, tmp)

  tmp <- rep(0, nrow(xtab))
  names(tmp) <- rownames(xtab)
  tmp2 <- sapply(drug.users, function(x){ RXlength(data,x,dc,drug)})
  tmp[drug.users] <- tmp2
  tmp <- tmp/365.25
  drugdur <- cbind(drugdur,tmp)  
}
colnames(drugcol) <- drugs
drugs_dur <- c("DOXAZOSIN")
colnames(drugdur) <- paste(drugs_dur,"dur",sep="_")

mydata <- data.frame(cbind(xtab, drugcol, drugdur))

##select subset to work with
#HC @ baseline
mydata.use <- subset(mydata, subset=DX==1 & !is.na(APOE))

#MCI @ baseline
mydata.use <- subset(mydata, subset=DX==2 & !is.na(APOE))

#HC or MCI @ baseline
mydata.use <- subset(mydata, subset=(DX==1 | DX==2) & !is.na(APOE))

AGE <- (mydata.use$ENTRY-mydata.use$BIRTH)/365.25

### analysis
APOE4      <- (mydata.use$APOE > 33 | mydata.use$APOE == 24) * 1
APOE2      <- (mydata.use$APOE < 33 ) * 1

##MMSE binning
#for normals OK
MMSE<-mydata.use$MMSE
mmse.bin <- MMSE
mmse.bin[MMSE <= 24] <- 4
mmse.bin[MMSE <= 26 & MMSE > 24] <- 3
mmse.bin[MMSE <= 28 & MMSE > 26] <- 2
mmse.bin[MMSE  > 28] <- 1

EDUX <- mydata.use$EDU
EDUX[EDUX >=99] <- NA

agestart   <- (mydata.use$ENTRY - mydata.use$BIRTH)/365.25
ageconvert <-  apply(mydata.use,1,function(x){
  bb <- min(x["TIMEMCI"],x["TIMEAD"])
  (bb - x["BIRTH"])/365.25
})
convert    <- apply(mydata.use[,c("CONVMCI","CONVAD")],1,max)

mysurv <- Surv(agestart, ageconvert, convert)

######Add num_drugs
ppp <- t(sapply(ids, function(y) {
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

#all subjects 
coxph_val <- coxph(mysurv ~ DOXAZOSIN_dur + NUMDRUGS + MMSE + EDUX + APOE2 + SEX * APOE4, data=mydata.use, subset = AGE >= 55)
###coxph(mysurv ~ LANSOPRAZOLE_dur + MMSE + EDUX + APOE2 + SEX * APOE4, data=mydata.use, subset= AGE >= 55)

#Write coxph values to file
out <- capture.output(coxph_val)
cat(out, file = "drug_coxph_summ.txt", sep = "\n", append=TRUE)

#limit to males
###coxph(mysurv ~ DOXAZOSIN_dur + MMSE + EDUX + APOE2 + APOE4, data=mydata.use, subset=SEX==1)
###coxph(mysurv ~ LANSOPRAZOLE_dur + MMSE + EDUX + APOE2 + APOE4, data=mydata.use, subset=SEX==1)


#limit to APOE4 carriers
###coxph(mysurv ~ DOXAZOSIN_dur + SEX + EDUX + MMSE + APOE2, data=mydata.use, subset=APOE4==1)
###coxph(mysurv ~ LANSOPRAZOLE_dur + SEX + EDUX + MMSE+ APOE2, data=mydata.use, subset=APOE4==1)

#in addition limit to males
###coxph(mysurv ~ DOXAZOSIN_dur + EDUX + MMSE + APOE2, data=mydata.use, subset=APOE4==1 & SEX==1)
###coxph(mysurv ~ LANSOPRAZOLE_dur + EDUX + MMSE +APOE2, data=mydata.use, subset=APOE4==1 & SEX==1)


##survival curve##
bh2scurve <- function(x){

  res <- basehaz(x,F)[,c(2,1)]
  res[,2] <- exp(-res[,2])
  return(res)
}


AGE <- (mydata.use$ENTRY-mydata.use$BIRTH)/365.25
#DOXAZOSIN
#mycox.1 <- coxph(mysurv ~ MMSE + EDUX + APOE2 + SEX + APOE4, data=mydata.use, subset=DOXAZOSIN==1 & AGE>=55)
#mycox.2 <- coxph(mysurv ~ MMSE + EDUX + APOE2 + SEX + APOE4, data=mydata.use, subset=DOXAZOSIN==0 & AGE>=55)
#
#plot(bh2scurve(mycox.1),type="s", lwd=2, lty=1)
#points(bh2scurve(mycox.2),type="s", lwd=2, lty=1, col=2)

###par(mfrow=c(1,2))
###plot(survfit(mysurv ~ DOXAZOSIN, data=mydata.use, subset=AGE>=55), col=c(1,2), xlim=c(50,110), main="DOXAZOSIN")
###plot(survfit(mysurv ~ LANSOPRAZOLE, data=mydata.use, subset=AGE>=55), col=c(1,2), xlim=c(50,110), main="LANSOPRAZOLE")


#durguration in users
###coxph(mysurv ~ DOXAZOSIN_dur + MMSE + EDUX + APOE2 + APOE4 + SEX, data=mydata.use, subset=DOXAZOSIN==1 & AGE >= 55)
###coxph(mysurv ~ LANSOPRAZOLE_dur + MMSE + EDUX + APOE2 + APOE4 + SEX, data=mydata.use, subset=LANSOPRAZOLE==1 & AGE >= 55)

}
