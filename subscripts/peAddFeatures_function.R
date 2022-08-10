peAddFeatures_list <- function(feat_path= "Rdata/final_300bp_enhancer_features.txt")
  names(fread(feat_path, nrows = 1))

peAddFeatures <- function(dat,
                          vars,
                          vars.regexpr= F,
                          aggregate_LR_FUN,
                          feat_path= "Rdata/final_300bp_enhancer_features.txt")
{
  if(vars.regexpr)
    sel <- grep(vars, names(fread(feat_path, nrows = 1))) else
      sel <- match(vars, names(fread(feat_path, nrows = 1)))
  add <- fread(feat_path, 
               sel= c(1, sel))
  cols <- names(add)[-1]
  if(!missing(aggregate_LR_FUN))
  {
    if(is.function(aggregate_LR_FUN))
    {
      L <- add[dat, cols, on= "ID==L", with= F][, idx:= .I]
      R <- add[dat, cols, on= "ID==L", with= F][, idx:= .I]
      .m <- rbind(L, R)
      .m <- .m[, lapply(.SD, aggregate_LR_FUN), idx]
      dat[, (cols):= .m[, -1]]
    }else
      stop("aggregate_LR_FUN should be a function")
  }else
  {
    dat[, paste0(cols, "_L"):= add[dat, cols, on= "ID==L", with= F]]
    dat[, paste0(cols, "_R"):= add[dat, cols, on= "ID==R", with= F]]
  }
}
