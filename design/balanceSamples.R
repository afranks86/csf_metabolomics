rm(list=ls())
library(ggplot2)

ADTracking <- read.csv("../data/ADTracking.csv")
PDTracking <- read.csv("../data/PDTracking.csv")

N <- nrow(ADTracking) + nrow(PDTracking)

maxBatch <- 30

nbatch <- ceiling(N / maxBatch)

table(ADTracking$Gender)
table(PDTracking$Sex)

## Add AD and PD data
Type <- c(as.character(ADTracking$PrimaryDx), as.character(PDTracking$Disease.Diagnosis))

PDTracking$Sex <- ifelse(PDTracking$Sex=="Female", "F", "M")
Gender <- c(as.character(ADTracking$Gender), as.character(PDTracking$Sex))

## Add ids
PDTracking$Id <- paste(paste(PDTracking$SubjectID,
                             PDTracking$Aliquot.Number, sep="#"),
                       PDTracking$PaNUCID, sep="; ")
ADTracking$Id <- ADTracking$Lab.Label.LPDate
Id <- c(as.character(ADTracking$Id), as.character(PDTracking$Id))
                       
Age <- c(ADTracking$LPAge, PDTracking$AgeAtDraw)
APOE <- c(as.character(ADTracking$APOEGen), rep("NA", nrow(PDTracking)))


NDTracking <- data.frame(list(Type=Type, Gender=Gender, Age=Age, APOE=APOE, Id=Id))

## APOE <- as.character(NDTracking$APOE)
## APOE[APOE == "NA"]  <- "Unknown"
## NDTracking$APOE <-  as.factor(APOE)

## 
## Based on Fite Selection Model (Morris, 1979)

set.seed(888)

assigned <- c()
remaining <- 1:N
groupAssignment <- rep(0, N)
count <- 1

######################################33
## Compute design matrix and rescale so all variables have unit variance
######################################33

formula <- "~  Age + Gender + Type + APOE + Type:Gender + Type:Age - 1"
mat <- model.matrix(as.formula(formula), dat=NDTracking)
tokeep <- which(colSums(mat) > 0)
mat <- mat[, tokeep]
mat <- apply(mat, 2, scale)

group.means <- matrix(0, nrow=nbatch, ncol=ncol(mat))

###########################################
## Core algorithm  
###########################################

for(round in 1:floor(N / nbatch)){

    ## for g in random group order
    for(g in sample(1:nbatch)){
        
        ## Find the observation that is _furthest_ from the current mean
        ## And select it for current group g
        MSE <- apply(sweep(mat[remaining, , drop=FALSE], 2, group.means[g, ], "-")^2, 1, sum)
        choice <- remaining[which.max(MSE)]

        ## add the choice, update group means
        groupAssignment[choice] <- g
        assigned <- c(assigned, choice)
        group.means[g, ] <- group.means[g, ]*(round-1)/(round)+mat[choice,] / round
        count <- count+1
        remaining <- setdiff(1:N, assigned)
    }

}

####################################################
## For remainder, find which minimize between group var
## Or more simply randomize the remainder
####n################################################

## function to calculate between group variance
calcB <- function(X, asgn.vec, nbatch){

    xbar <- apply(X, 2, mean)
    
    groupMeans <- sapply(1:nbatch,
                         function(g) apply(X[which(asgn.vec==g), , drop=FALSE], 2, mean))
    B <- (nrow(X)/nbatch)*(groupMeans-xbar) %*% t(groupMeans-xbar)
    return(B)
}

gasgn <- c()
for(r in remaining) {

    betweenGroupSS <- 
        sapply(setdiff(1:nbatch, gasgn), function(g) {
            tmpAsgn <- groupAssignment
            tmpAsgn[r] <- g
            B <- calcB(mat, tmpAsgn, nbatch)
            sum(diag(B %*% t(B)))
        })
    gchoice <- setdiff(1:nbatch, gasgn)[which.min(betweenGroupSS)]
    groupAssignment[r] <- gchoice
    assigned <- c(assigned, r)
    gasgn <- c(gasgn, gchoice)
}

## Uncomment below to evaluate for a true randomization
## groupAssignment <- c(sapply(1:28, function(x) sample(1:7)), sample(1:7, size=3))

NDTracking$Batch <- groupAssignment

## Check counts per batch
table(NDTracking$Batch)

## Check age balance
pdf("../figs/design/age_balance.pdf")
t(sapply(1:nbatch, function(g) quantile(NDTracking$Age[NDTracking$Batch==g])))
ggplot(data=NDTracking, aes(y = Age, x=as.factor(Batch))) +
    geom_boxplot(fill="grey") +
    theme_bw(base_size=20) + xlab("Batch")
dev.off()

## Check gender balance
pdf("../figs/design/gender_balance.pdf")
table(NDTracking$Gender, NDTracking$Batch)
ggplot(data=NDTracking, aes(Batch)) + geom_bar() + facet_wrap(~ Gender) +
    theme_bw(base_size=20) + ylab("Count")
dev.off()

## Check balance of disease types
pdf("figs/design/type_balance.pdf")
table(NDTracking$Type, NDTracking$Batch)
ggplot(data=NDTracking, aes(Batch)) + geom_bar() + facet_wrap(~ Type)
dev.off()

## Check balance of APOE status
pdf("figs/design/apoe_balance.pdf")
table(NDTracking$APOE, NDTracking$Batch)
ggplot(data=NDTracking, aes(Batch)) + geom_bar() + facet_wrap(~ APOE)
dev.off()

## Interaction between gender and type
pdf("figs/design/gender_type_balance.pdf")
table(paste(as.character(NDTracking$Type), as.character(NDTracking$Gender), sep="-"), NDTracking$Batch)
ggplot(data=NDTracking, aes(Batch)) + geom_bar() + facet_wrap(~ Gender + Type)
dev.off()

## APOE and Gender
pdf("figs/design/apoe_gender_balance.pdf")
table(paste(as.character(NDTracking$APOE), as.character(NDTracking$Gender), sep="-"), NDTracking$Batch)
dev.off()

## Age and Type
pdf("figs/design/age_type_balance.pdf")
ggplot(data=NDTracking[-163, ], aes(y = Age, x=as.factor(Batch), fill=Type)) + geom_boxplot()
dev.off()

## write.csv(NDTracking[, c("Id", "Batch")], file="SampleBatches.csv", row.names=FALSE)

if(FALSE) {
    write.csv(NDTracking, file="../data/NDTracking.csv", row.names=FALSE)
}
