rescale <- function(x) {
  rng <-range(x)
  (x - rng[1]) / (rng[2] - rng[1])
}

rescale2 <- function(x, na.rm=TRUE, plot=TRUE) {
  rng <-range(x, na.rm=na.rm)
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  if(plot) {
    plot(answer, typ="b")
  }
  return(answer)
}
rescale2( c(1,2,NA,3,10), plot=FALSE )


rescale <- function(x, na.rm=TRUE, plot=FALSE) {
  if(na.rm) {
    rng <-range(x, na.rm=na.rm)
  } else {
    rng <-range(x)
  }
  print("Hello")
  
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  
  print("is it me you are looking for?")
  
  if(plot) { 
    plot(answer, typ="b", lwd=4) 
  }
  print("I can see it in ...")
}

