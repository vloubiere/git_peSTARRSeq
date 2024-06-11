test <- as.data.table(readRDS("Rdata/vl_library_twist015_112022.rds"))
test <- test[sublib=="A"]

test[, c("class", "mut", "ID"):= tstrsplit(ID, "_", keep= c(1,2,4))]
test[, start:= substr(oligo_full_sequence, 1, 140)]
test[, dist:= stringdist(start, start[mut=="WT"]), .(class, ID)]
table(test[mut!="WT", dist])

table(test$dist)
