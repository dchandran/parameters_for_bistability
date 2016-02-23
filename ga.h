#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mtrand.h"

#ifndef GA_MAIN_LOOP
#define GA_MAIN_LOOP

#define NO_FPRINTF_OUTPUT 1

/* Population of individuals -- an individual is represented by a user struct (void*) */
typedef void ** Population;

/***********************************************************************
    The following functions MUST be defined somewhere
***********************************************************************/

/*
 * Free an individual from memory
 * @param: a single individual
*/
extern void deleteIndividual(void *);
/*
 * Make a copy of an individual and return the memory pointer
 * @param: target individual
*/
extern void * clone(void *);
/*
 * Compute fitness of an individual. Fitness must be positive if default selection function is used
 * @param: target individual
 * @ret: fitness (double) of the individual (MUST BE POSITIVE if default selecte)
*/
typedef double(*GAFitnessFnc)(void *);

/************************************************************************************************************
    At least one of the following two functions SHOULD (not must) be defined, otherwise the GA will not evolve
****************************************************************************************************************/
/*
 * combine two individuals to generate a new individual
 * @param: parent individual 1
 * @param: parent individual 2
 * @ret: pointer to an individual (can be the same as one of the parents)
*/
typedef void* (*GACrossoverFnc)(void *, void *);
/*
 * Change an individual randomly to generate a new individual
 * @param: parent individual
 * @ret: pointer to an individual (can be the same as one of the parents)
*/
typedef void* (*GAMutateFnc)(void *);

/************************************************************************************************************
  The following two functions are entirely optional. They may or may not affect the GA performance 
****************************************************************************************************************/
/*
 * Selection function. GA provides a default selection function that linearly converts fitness values to probabilities
 * @param: Population of individuals
 * @param: array of fitness values for the individuals
 * @param: total fitness (sum of all fitness values)
 * @param: number of individuals in the population
 * @ret: index (in population vector) of the individual to select
*/
typedef int(*GASelectionFnc)(Population , double * , double , int );
/*
 * Callback function. This function is called during each iteration of the GA.
 * @param: iteration
 * @param: Population of individuals
 * @param: number of individuals in the population
 * @ret: 0 = continue GA, 1 = stop GA. This can be used to stop the GA before it reaches max iterations
*/
typedef int(*GACallbackFnc)(int,Population,int);

/************************************************************************************************************
  The central functions of the genetic algorithm
****************************************************************************************************************/

/*
 * Selects an individual at random, with probability of selection ~ fitness
 * @param: array of individuals
 * @param: array of corresponding fitness values
 * @param: sum of all fitness values
 * @param: number of individual
*/
int GAselect(Population , double * , double , int );
/*
 * Get next population from current population
 * @param: array of individuals
 * @param: number of individual in population currently
 * @param: number of individual in the new population (returned array)
 * @param: fitness function pointer
 * @param: crossover function pointer
 * @param: mutation function pointer
 * @param: selection function pointer
 * @param: 0 = delete old population, 1 = keep old population (warning: user must delete it later)
 * @ret: new array of individual (size = 3rd parameter)
*/
Population GAnextGen(Population,int,int,GAFitnessFnc,GACrossoverFnc,GAMutateFnc,GASelectionFnc,short);

/*
 * The main GA loop
 * @param: array of individuals
 * @param: number of individual in the initial population
 * @param: number of individual to be kept in the successive populations
 * @param: total number of generations
 * @param: fitness function pointer
 * @param: crossover function pointer
 * @param: mutation function pointer
 * @param: selection function pointer
 * @param: 0 = delete old population, 1 = keep old population (warning: user must delete it later)
 * @ret: final array of individuals (sorted)
*/
Population GArun(Population,int,int,int,GAFitnessFnc,GACrossoverFnc,GAMutateFnc, GACallbackFnc);

/*
 * sort (Quicksort) a population by its fitness
 * @param: population to sort
 * @param: fitness function
 * @param: size of population
 * @ret: void
*/
void GAsort(Population population, GAFitnessFnc, int);

#endif

