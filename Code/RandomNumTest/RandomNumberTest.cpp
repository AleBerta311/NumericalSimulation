//We want to test a linear random generator by estimating the average of a large sample of randomly generated numbers
//It should be 1/2 and then we can also estimate the variance with the block average method

#include "random.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

double error(double av, double av2)
    {
        return sqrt(av2 - av*av);
    }

int main()
{
    Random myGen;

    //We need to extract some primes and seed(s) from the files listed for the "advanced" random generator to work

    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open())
    {
        Primes >> p1 >> p2;
    }
    else cerr << "Unable to open file Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open())
    {
        while(! input.eof() )
        {
            input >> property;
            if(property == "RANDOMSEED")
            {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            }
        }

        myGen.SetRandom(seed, p1, p2);
        input.close();
    }
    else cerr << "Cannot open seed.in" << endl;

    ofstream myresults;

    myresults.open("myresults.txt");

    //Ok now we can use the generator

    int M = 100000;         // number of "throws" i.e randomly generated nums
    int N = 100;            // number of blocks
    int L = M/N ;           // number of throws in each block

    // We could generate M random numbers already and put them in a vector but there's no need to waste space
    // We'll generate one at a time, and work with cumulative sums and averages    

    //Now we do the actual stuff and we'll write all the results for the averages and errors in a file to be plotted in jupyter notebook

    double ave = 0;
    double av2 = 0;
    double sum_prog = 0;
    double su2_prog = 0;
    double err_prog = 0;

    double var_ave = 0;
    double var_2_ave = 0;
    double var_sum_prog = 0;
    double var_2_sum_prog = 0;
    double var_error_prog = 0;

    for (int i = 0; i < N; i++)
    {
        double sum = 0;
        double var_progr = 0;

        for (int j = 0; j < L; j++)
        {
            double r = myGen.Rannyu();  //random number between 0 ( included ) and 1 (excluded)
            sum+= r;
            var_progr = var_progr + (r-0.5)*(r-0.5);
        }
        ave = sum/L;
        av2 = ave*ave;

        sum_prog+= ave;
        su2_prog+= av2;


        var_ave = var_progr/L;
        var_2_ave = var_ave*var_ave;

        var_sum_prog+= var_ave;
        var_2_sum_prog+= var_2_ave;


        if (i != 0) 
        {
            err_prog = error(sum_prog/(i+1), su2_prog/(i+1))/sqrt(i);
            var_error_prog = error (var_sum_prog/(i+1), var_2_sum_prog/(i+1))/sqrt(i);
        }
            

        myresults << ave << " ";
        myresults << err_prog << " ";
        myresults << sum_prog/(i+1) << " ";
        myresults << var_ave << " ";
        myresults << var_sum_prog/(i+1) << " ";
        myresults << var_error_prog;
        myresults << endl;
        
    } 

    myresults.close();


    return 0;
}