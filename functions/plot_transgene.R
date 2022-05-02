transgene <- function(x,
                      y,
                      CP= "dev",
                      rotate= F,
                      cex= 0.6)
{
  width <- strwidth(" dev ", cex= cex)*6
  CP_height <- strheight("dev", cex= cex)*1.5
  if(rotate)
  {
    width <- diff(grconvertY(c(0, diff(grconvertX(c(0, width), "user", "inches"))), "inches", "user"))
    CP_height <- diff(grconvertX(c(0, diff(grconvertY(c(0, CP_height), "user", "inches"))), "inches", "user"))
  }
  RE.w <- width/6
  if(rotate)
  {
    # TG
    segments(x0= x,
             y0= y,
             x1= x,
             y1= y+width,
             xpd= T,
             lend= 2)
    # TSS
    lines(x-CP_height*c(0,1,1),
          y+RE.w+c(0,0,RE.w/2),
          xpd= T,
          lend= 2)
    polygon(x-CP_height*c(1.2, 0.8, 1),
            y+RE.w*1.5+c(0,0,RE.w*0.25),
            col= "black",
            border= NA,
            xpd= T,
            lend= 2)
    # CP
    rect(x-CP_height/2,
         y+RE.w/1.8,
         x+CP_height/2,
         y+RE.w*1.5,
         col= switch(CP, 
                     "dev"= "#74C27A",
                     "hk"= "tomato",
                     "white"),
         xpd= T,
         lend= 2)
    text(x,
         y+RE.w,
         CP,
         offset= 0,
         srt= 90,
         cex= cex,
         xpd= T)
    # Enhancers
    segments(x,
             y+RE.w*c(2.5,4.5),
             x,
             y+RE.w*c(2.5,4.5)+RE.w,
             lwd= 4,
             col= "darkgrey",
             xpd= T,
             lend= 2)
  }else
  {
    # TG
    segments(x0= x,
             y0= y,
             x1= x+width,
             y1= y,
             xpd= T,
             lend= 2)
    # TSS
    lines(x+RE.w+c(0,0,RE.w/2),
          y+CP_height*c(0,1,1),
          xpd= T,
          lend= 2)
    polygon(x+RE.w*1.5+c(0,0,RE.w*0.25),
            y+CP_height*c(1.2, 0.8, 1),
            col= "black",
            border= NA,
            xpd= T,
            lend= 2)
    # CP
    rect(x+RE.w/2,
         y-CP_height/2,
         x+RE.w*1.5,
         y+CP_height/2,
         col= switch(CP, 
                     "dev"= "#74C27A",
                     "hk"= "tomato",
                     "white"),
         xpd= T,
         lend= 2)
    text(x+RE.w,
         y,
         CP,
         offset= 0,
         cex= cex,
         xpd= T)
    # Enhancers
    segments(x+RE.w*c(2.5,4.5),
             y,
             x+RE.w*c(2.5,4.5)+RE.w,
             y,
             lwd= 4,
             col= "darkgrey",
             xpd= T,
             lend= 2)
  }
}
