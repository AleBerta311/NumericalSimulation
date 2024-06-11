#include "random.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

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

    ofstream my_chi2;
    my_chi2.open("chi2.txt");

    int M = 1000000;         // number of "throws" i.e randomly generated nums
    double N = 100.00;       // number of blocks
    int L = M/N ;           // number of throws in each block

    // We could generate M random numbers already and put them in a vector but there's no need to waste space
    // We'll generate one at a time, and work with cumulative sums and averages    

  
    double r;
    double chi_2_progr = 0;
    int count = 0;

    for(int k = 0; k < 10000; k++)
    {
        chi_2_progr = 0;
        for (int i = 0; i < N; i++)
        {
            count = 0;
            double min = i/N;
            double max = (i+1)/N;

            for (int j = 0; j < L; j++)
            {
                r = myGen.Rannyu();
                
                if(r > min and r < max)
                {
                    count+= 1;
                }
                
            }

            chi_2_progr+= ((count - L/N)* (count -L/N));
        
        }

        chi_2_progr/=(L/N);
        my_chi2 << chi_2_progr << endl;
    }
    
    return 0;
}