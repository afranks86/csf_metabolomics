rm(list=ls())
NDTracking <- read.csv("../data/NDBatchTracking.csv")
SampleOrder <- read.csv("../data/SampleOrder.csv")

merged <- merge(NDTracking, SampleOrder, by=c("Id"))

merged <- merged[, c("Id", "Type", "Gender", "Age", "APOE", "Batch.y", "Index")]
colnames(merged)[6] <- "Batch"

write.csv(merged, file="../data/NDTracking.csv", row.names=FALSE)
