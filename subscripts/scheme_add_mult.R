

tab <- data.table(tCP= 1,
                  ctlL= c(T,T,F,F),
                  ctlR= c(T,F,F,T),
                  indL= c(1,4,1,4),
                  indR= c(1,1,4,4))
tab[, add:= indL+indR-tCP]
tab[, mult:= indL*indR/tCP]

par(mfrow=c(2,1))
vl_plot_table(tab)
vl_plot_table(log2(tab))