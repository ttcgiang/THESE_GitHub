#all: EvolVacc InterVacc BatchVacc
all: EvolVacc BatchVacc
################################################
## This Makefile allows to create the Targets  #
## EvolVacc and InterVacc. EvolVacc supports   #
## finding automatically the optimum parms     #
## InterVacc allows to create the rpogram      #
## that interacts with the user.               #
##                                             #
## Version 1.0 of 5/12/2011 by Giang           #
################################################
	
#InterVacc: EvoStrat.o
#	g++ $(INCLUDE) $(LIB) -o $@ $<
	
BatchVacc: BatchVacc.o
	g++ $(INCLUDE) $(LIB) -o $@ $<
BatchVacc.o:BatchVacc.cpp
	g++ -c $(INCLUDE) -o $@ BatchVacc.cpp



EvolVacc: EvolVacc.o
	g++ $(INCLUDE) $(LIB) -o $@ $<

EvolVacc.o:EvolVacc.cpp
	g++ -c $(INCLUDE) -o $@ EvolVacc.cpp

# more general option
.cpp:
	g++ $(INCLUDE) $(LIB) -o $@ $<
	g++ -lgraph

clean:
	rm -f *.o *~ core*
	ls -l

touch:
	touch  BatchVacc.cpp  EvolVacc.cpp  

