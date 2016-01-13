stringFixe <- function(x,n)
{
l= nchar(x)
k=n-l+1
tp=""
if(k<0) 
	x = substr(x,0,n)
else if(k>0)
{
for(i in 0 : (k-1)) tp <- paste(tp,"");
x = paste(x,tp)
x=substr(x,0,n)
}
else x = x;
}
