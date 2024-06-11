//A program to compute an integral using montecarlo techniques
//the function we want to integrate is f(x) = pi/2*cos(pi*x/2) in [0, 1]

#include "random.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

using namespace std;

void errors(const vector<double> &av, const vector<double> &av2, vector<double> &stds)
    {
        double progr;
        double progr2;
        for(int i = 1; i<100; i++)
        {
            progr = 0;
            progr2 = 0;
            for (int j = 0; j<=i; j++)
            {
                progr+= av[j];
                progr2+= av2[j];
                
            }
            stds[i] = (progr2/(i+1))- pow( progr/(i+1) ,2); 
            stds[i] = sqrt(stds[i]/i);
        }
    }

double myFunction(double x)
    {
        return M_PI/2*cos(M_PI*x/2);
    }

double myNewPx(double x)
    {
        double C = 1 - 1/M_E;       //norm cost for p(x)
        return 1/C*exp(-x);
    }

double myNewFunction(double x)
    {
        
        
        //cout << myFunction(x) << endl;
        return myFunction(x)/ myNewPx(x);

    }

int main()
{
    ofstream myoutput;
    myoutput.open("results.txt");

    ofstream mynewoutput;

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

    int n_throws = 200;       // throws per evaluation
    int n_blocks = 100;
    int n_per_block = 100;      // evaluations per block
    double interval = 1;        // (b-a)
    double f_ave;
    double integral;
    double int2;
    double x;

    vector<double> integrals; 
    vector<double> ints2;

    for (int j = 0; j< n_blocks; j++)
    {
        integral = 0;
        int2 = 0;
        for (int k = 0; k < n_per_block; k++)
        {
            f_ave = 0;
            for (int i = 0; i<n_throws; i++)
            {
                x = myGen.Rannyu();

                f_ave+= myFunction(x); 
            }
        
            f_ave/=(n_throws);

            integral+= f_ave/(interval);
            

        }
        
        integral/= n_per_block;
        int2+= integral*integral;

        integrals.push_back(integral);
        ints2.push_back(int2);
        
    }

    vector<double> std_devs(100);
    
    errors(integrals, ints2, std_devs);
    double cumul_int = 0;

    for (int i = 0; i<100; i++)
    {
        //cout << integrals[i] << endl;
        cumul_int = 0;
        for (int j= 0; j<=i; j++)
            cumul_int+= integrals[j];

        if (i != 0)
        {
            cumul_int/=(i+1);
            myoutput << cumul_int << " " << std_devs[i] << endl;
        }
            

        
        
    }
     
    myoutput.close();

    // same stuff but with importance sampling and p(x) = e^-x

    mynewoutput.open("newresults.txt");

    vector<double> new_integrals; 
    vector<double> new_ints2;
    double y;

    for (int j = 0; j< n_blocks; j++)
    {
        integral = 0;
        int2 = 0;
        for (int k = 0; k < n_per_block; k++)
        {
            f_ave = 0;
            for (int i = 0; i<n_throws; i++)
            {
                
                
                x = myGen.Rannyu(0, 1- 1/M_E);
                y = -log(1- x);   
        
                f_ave+= myNewFunction(y); 
            }
        
            f_ave/=(n_throws);

            integral+= f_ave/(interval);
            
        }
        
        integral/= n_per_block;
        int2+= integral*integral;

        new_integrals.push_back(integral);
        new_ints2.push_back(int2);
        
    }

    vector<double> new_std_devs(100);
    
    errors(new_integrals, new_ints2, new_std_devs);
    cumul_int = 0;

    for (int i = 0; i<100; i++)
    {
        cumul_int = 0;
        for (int j= 0; j<=i; j++)
            cumul_int+= new_integrals[j];

        if (i != 0)
        {
            cumul_int/=(i+1);
            mynewoutput << cumul_int << " " << new_std_devs[i] << endl;
        }
            

    }
     
    mynewoutput.close();
    
    return 0;
}
