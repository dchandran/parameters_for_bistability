/****************************************************************************
 **
 ** Copyright (C) 2008 Deepak Chandran
 ** Contact: Deepak Chandran (dchandran1@gmail.com)
 **
 ****************************************************************************/
 
#include "ga.h"
#include <stdio.h>
/*
 * Selects an individual at random, with probability of selection ~ fitness
 * @param: array of individuals
 * @param: array of corresponding fitness values
 * @param: sum of all fitness values
 * @param: number of individual
 * @ret: index of selected individual in the population
*/
int GAselect(Population population, double * fitnessValues, double sumOfFitness, int popSz)
{
   int i;
   double randNum = mtrand() * sumOfFitness, 
          total = 0;
   for (i=0; i < popSz-1; ++i)
       if (total < randNum && randNum < (total+fitnessValues[i]))
          return (i);
       else
          total += fitnessValues[i];
   return (i);
}

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
Population GAnextGen(Population currentPopulation, int oldPopSz, int newPopSz,
                     GAFitnessFnc fitness, GACrossoverFnc crossover, GAMutateFnc mutate,
                     GASelectionFnc select,
                     short keepOldPopulation)
{
   //allocate memory for next generation
   Population nextPopulation = malloc( newPopSz * sizeof(void*) );
   if (nextPopulation == NULL) 
   {
      return (0);
   }
   int i,k;
   //make array of fitness values
   double * fitnessArray = malloc ( oldPopSz * sizeof(double) );
   double totalFitness = 0;
   int best = 0;  //save best's index

   for (i = 0; i < oldPopSz; ++i)
   {
      fitnessArray[i] = fitness(currentPopulation[i]);
      if (fitnessArray[i] < 0) fitnessArray[i] = 0;   //negative fitness not allowed

      totalFitness += fitnessArray[i];
      if (fitnessArray[i] > fitnessArray[best]) 
         best = i;
   }

   //keep the best
   nextPopulation[0] = clone(currentPopulation[best]);

   //select the fit individuals
   void * x1 = NULL, * x2 = NULL;
   for (i = 1; i < newPopSz; ++i)
   {
      k = select(currentPopulation,fitnessArray,totalFitness,oldPopSz);

      x1 = currentPopulation[k];
      if (crossover != NULL) 
      {
         double temp = fitnessArray[k];
         fitnessArray[k] = 0;   //this is to prevent self-self crossover
         int k2 = select(currentPopulation,fitnessArray,totalFitness,oldPopSz);
         fitnessArray[k] = temp;
         x2 = currentPopulation[k2];
         x1 = crossover(x1,x2);
      }
      else
      {
         x1 = clone(x1);
      }

      if (mutate != NULL) 
      {
         x1 = mutate(x1);
      }
      nextPopulation[i] = x1; //add to the new population
   }
   /*free the memory from the old population*/
   if (keepOldPopulation == 0)
   {
 
     for (i = 0; i < oldPopSz; ++i)
        if (currentPopulation[i] != NULL)
           deleteIndividual(currentPopulation[i]);
     free(currentPopulation);
   }
   free(fitnessArray);
   return (nextPopulation);
}

/*
 * The main GA loop
 * @param: array of individuals
 * @param: number of individual initially
 * @param: number of individual in successive populations
 * @param: total number of generations
 * @param: fitness function pointer
 * @param: crossover function pointer
 * @param: mutation function pointer
 * @param: selection function pointer
 * @param: 0 = delete old population, 1 = keep old population (warning: user must delete it later)
 * @ret: final array of individuals (sorted)
*/
Population GArun(Population initialPopulation, int initPopSz, int popSz, int numGenerations,
                 GAFitnessFnc fitness, GACrossoverFnc crossover, GAMutateFnc mutate, 
                 GACallbackFnc callback)
{
   FILE * errfile = freopen("GArun_errors.log", "w", stderr);

   initMTrand(); /*initialize seeds for MT random number generator*/
   int i = 0, stop = 0;
   Population population = initialPopulation;

   while (stop == 0)
   { 
      if (i == 0)
         population = GAnextGen(population, initPopSz, popSz, fitness, crossover, mutate, &GAselect, 0);
      else
         population = GAnextGen(population, popSz, popSz, fitness, crossover, mutate, &GAselect, 0);

      if (callback != NULL)
         stop = callback(i,population,popSz);

     ++i;
     if (i >= numGenerations) stop = 1;
   }
   GAsort(population,fitness,popSz);

   fclose(errfile);
   return (population);
}

/***********************************************************************
    *  Quicksort code from Sedgewick 7.1, 7.2.
***********************************************************************/

    // is x < y ?
    int less(double x, double y) {
        return (x > y);
    }

    // exchange a[i] and a[j]
    void exch(Population population, double* a, int i, int j) {
        double swap = a[i];
        a[i] = a[j];
        a[j] = swap;

        void * temp = population[i];
        population[i] = population[j];
        population[j] = temp;
    }

    // partition a[left] to a[right], assumes left < right
    int partition(Population population, double* a, int left, int right) {
        int i = left - 1;
        int j = right;
        while (1) {
            while (less(a[++i], a[right]))      // find item on left to swap
                ;                               // a[right] acts as sentinel
            while (less(a[right], a[--j]))      // find item on right to swap
                if (j == left) break;           // don't go out-of-bounds
            if (i >= j) break;                  // check if pointers cross
            exch(population, a, i, j);         // swap two elements into place
        }
        exch(population, a, i, right);                      // swap with partition element
        return i;
    }

    // quicksort helper a[left] to a[right]
    void quicksort(Population population, double* a, int left, int right) {
        if (right <= left) return;
        int i = partition(population, a, left, right);
        quicksort(population, a, left, i-1);
        quicksort(population, a, i+1, right);
    }

    //quicksort
    void GAsort(Population population, GAFitnessFnc fitness, int populationSz) 
    {
        double * a = malloc ( populationSz * sizeof(double) );
        int i;
        for (i=0; i < populationSz; ++i)
        {
            a[i] = fitness(population[i]);
        }
        if (a != NULL)
        {
           quicksort(population, a, 0, populationSz - 1);
           free(a);
        }
    }

