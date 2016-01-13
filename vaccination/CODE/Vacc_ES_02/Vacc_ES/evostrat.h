/*
 *  evostrat.h
 *  EvoStrat
 *
 *  Created by yann chevaleyre on 24/11/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <vector>
#include <string>

using namespace std;

typedef vector<double> double_vector;

double normal_distribution();
void test_normal_distribution();

void test_evo_strat();
double_vector evo_strat(double_vector v0,double (*eval)(double_vector*),int ntrials,double sigma,bool adaptive,double_vector* fitness_sequence);

