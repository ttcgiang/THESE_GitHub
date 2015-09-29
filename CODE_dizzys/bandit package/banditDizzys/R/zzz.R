####
# The functions below help us install automatically the R dependant packages
# that the package 'dizzys' needs
#
# checking a R package installed or not yet
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
# installing the dependent packages automatically
# the dependant packages are called from "http://cran.r-project.org"
callpackg <- function(pckgname="deSolve"){
	if(!is.installed(pckgname)) {
		invisible(install.packages(pckgname, repos="http://cran.r-project.org"))
		#checking packages
		if(require(pckgname,character.only = TRUE)){
       			 print(paste("",pckgname," installed and loaded"))
   		} else {#error
    		    stop(paste("could not install ",pckgname))
  		}

	}
	else
		invisible(require(pckgname,character.only = TRUE))
}

# names of the dependant packages, we need
.onLoad <- function(libname, pkgname){
     # do whatever needs to be done when the package is loaded
     # some people use it to bombard users with 
     # messages using 
	packageStartupMessage( "My package is so cool" )
	packageStartupMessage( "so I will print these lines each time you load it")
	# package 'deSolve', we call the fucntion 'ODE'
	callpackg("dizzysNEWYANN")
}
##############################################################################

