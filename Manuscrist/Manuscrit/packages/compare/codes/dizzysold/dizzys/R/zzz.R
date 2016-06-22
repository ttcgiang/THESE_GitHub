is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
callpackg <- function(pckgname="deSolve"){
	if(!is.installed(pckgname)) {
		invisible(install.packages(pckgname, repos="http://cran.r-project.org"))
		if(require(pckgname,character.only = TRUE)){
       			 print(paste("",pckgname," installed and loaded"))
   		} else {
    		    stop(paste("could not install ",pckgname))
  		}

	}
	else
		invisible(require(pckgname,character.only = TRUE))
}

.onLoad <- function(libname, pkgname){
     # do whatever needs to be done when the package is loaded
     # some people use it to bombard users with 
     # messages using 
	packageStartupMessage( "My package is so cool" )
	packageStartupMessage( "so I will print these lines each time you load it")
	callpackg("deSolve")
	callpackg("rgl")
	callpackg("survival")
	callpackg("KMsurv")
	
}

