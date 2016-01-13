all: EvoStrat
################################################
## This Makefile allows to create the Targets  #
## EvolVacc and InterVacc. EvolVacc supports   #
## finding automatically the optimum parms     #
## InterVacc allows to create the rpogram      #
## that interacts with the user.               #
##                                             #
## Version 1.0 of 5/12/2011 by Giang           #
################################################
	
EvoStrat:EvoStrat.o
	g++ $(INCLUDE) $(LIB) -o $@ $<
EvoStrat.o:EvoStrat.cpp
	g++ -c $(INCLUDE) -o $@ EvoStrat.cpp

# more general option
.cpp:
	g++ $(INCLUDE) $(LIB) -o $@ $<
	g++ -lgraph

clean:
	rm -f *.o *~ core*
	ls -l

touch:
	touch  EvoStrat.cpp  

