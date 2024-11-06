## Min-max scaling of x using the mean of top and bottom n values as the range
scaleDistro <- function(x, tail.size=20){
  top <- mean(sort(x, decreasing=T)[1:tail.size])
  bot <- mean(sort(x, decreasing=F)[1:tail.size])
  
  res <- ( x - bot ) / (top-bot)
  return(res)
}