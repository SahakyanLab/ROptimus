ex.r.fun <- function(K){
  K.new <- K
  # Indices of two randomly chosen js to swap
  new.ind <- sample(x=1:length(K.new[,"j"]), size=2, replace=FALSE)
  # Swap the js
  K.new[new.ind,"j"] <- K.new[rev(new.ind),"j"]
  return(K.new)
}
usethis::use_data(ex.r.fun, overwrite = TRUE)
