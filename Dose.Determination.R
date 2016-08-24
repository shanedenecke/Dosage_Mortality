
gm=function(min,max){
  a=c(0,(max-min))
  b=c()
  for (i in 1:4){
    b[i]=(a/(2*i))+min
  }
  b=rev(b)
  b=c(min,b,max)
  return(b)
}


gm2=function(min,max){
   a=c(0,(max-min))
   b=a[2]/8
   c=c(min,min+b,min+(b*2),min+(b*4),max)
  return(c)
   }
   
dd <- function(min,max,Lc50){
    b <- c(min,mean(c(min,Lc50)),Lc50,mean(c(max,Lc50)),max)
    return(b)
}