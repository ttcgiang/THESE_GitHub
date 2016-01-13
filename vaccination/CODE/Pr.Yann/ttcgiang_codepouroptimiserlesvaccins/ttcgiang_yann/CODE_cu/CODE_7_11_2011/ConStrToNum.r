ConStrToNum<-function(x)
{
x = sapply(strsplit(x, split = "/", fixed="TRUE"), as.numeric)
if (length(x)==2)
{
x = apply(x,2, function(x) x[1] / x[2])
}
else
x=x
}
