all: VaccinationSEIR_EvoStrat_version1
# one option
	
EvoStrat:EvoStrat.o
	g++ $(INCLUDE) $(LIB) -o $@ $<
EvoStrat.o:VaccinationSEIR_EvoStrat_version1.cpp
	g++ -c $(INCLUDE) -o $@ VaccinationSEIR_EvoStrat_version1.cpp

# more general option
.cpp:
	g++ $(INCLUDE) $(LIB) -o $@ $<
	g++ -lgraph

clean:
	rm -f *.o *~ core*
	ls -l

touch:
	touch  VaccinationSEIR_EvoStrat_version1.cpp  

