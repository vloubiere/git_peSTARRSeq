if(!exists("test"))
{
  test <- CJ(L= letters[1:4],
             R= letters[1:4])
  set.seed(1)
  umis <- c(sapply(seq(500), function(x) paste0(sapply(seq(10), function(y) sample(c("A", "T", "C", "G"), 1)), collapse= "")))
  set.seed(1)
  umis <- rep(umis, sample(3, 500, replace = T))
  test <- test[, .(UMI= umis), .(L, R)]
}
.c <- copy(test)
# Count UMIs and order
.c <- data.table(L= "a", R= "b", UMI= c("AAAAAAAAAA", "AAAAAAAAAT", "AAAAAAAATT", "AAAAAAATTT", "AAAAAATTTT", "AAAAATTTTT"))
.c <- .c[, .(umi_N= .N), .(L, R, UMI)]
length(unique(.c$UMI))
.c <- .c[!agrepl("GGGGGGGGG", UMI)]
setorderv(.c, "umi_N", order = -1)
# Identify UMIs that might require collapsing
t1 <- Sys.time()
.c[, done:= T]
# done <- rep(T, nrow(.c))
for(i in 0:9) 
  .c[(done), done:= ifelse(.N>1, F, T), .(L, R, sub(paste0("(.{", i, "})."), "\\1", UMI))]
# Advanced UMI collapsing (>1 diff)
while(any(!.c$done))
{
  .c[!(done), c("UMI", "done") := {
    idx <- agrepl(UMI[1], UMI)
    .(ifelse(idx, UMI[1], UMI), idx)
  }, .(L, R)]
}
t2 <- Sys.time()
t2-t1
length(unique(.c$UMI))
