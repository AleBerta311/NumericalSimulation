// We want to generate average values ofrom different distributions and compare them

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

    ofstream myout;
    myout.open("distroresults_unif.txt");

    //First we work with the uniform distro

    double average = 0;

    for (int i = 0; i< 10000; i++)
    {   
        average = 0;
       // N = 1

        average+= myGen.Rannyu();
        myout << average << " ";

        average = 0;

        // N = 2

        average+= myGen.Rannyu();
        average+= myGen.Rannyu();

        average/=(2);

        myout << average << " ";

        // N = 10

        average = 0;
        for (int i = 0; i <10 ; i++)
        {
            average+= myGen.Rannyu();
        }

        average/=10;
        myout << average << " ";

        // N = 100

        average = 0;

        for (int i = 0; i<100; i++)
        {
            average+= myGen.Rannyu();
        }
        average/=100;
        myout << average << endl;
    
    }

    myout.close();

    //Now we work with the expo and its inverse is f^-1 (x) = ln(lambda) -ln(x)/lambda
    double x;

    myout.open("distroresults_expo.txt");
    for (int i = 0; i< 10000; i++)
    {
        average = 0;
        // N = 1

        x = myGen.Rannyu();
        average+= -log(x);
        myout << average << " ";

        average = 0;

        // N = 2

        x = myGen.Rannyu();
        average+= -log(x);
        x = myGen.Rannyu();
        average+= -log(x);

        average/=(2);

        myout << average << " ";

        // N = 10

        average = 0;
        for (int i = 0; i <10 ; i++)
        {
            x= myGen.Rannyu();
            average+= -log(x);
            
        }
        average/=10;
        myout << average << " ";

        // N = 100

        average = 0;

        for (int i = 0; i<100; i++)
        {
            x= myGen.Rannyu();
            average+= -log(x);
            
        }
        average/=100;
        myout << average << endl;
    
    }

    myout.close();

    //Now we work with the Lorentz. and its inverse is f^-1 (x) = sqrt (1/pi*x -1) and  x must be in [0, 1/pi]

    myout.open("distroresults_lorentz.txt");

    for (int i = 0; i< 10000; i++)
    {
        // N = 1
        average = 0;

        x = myGen.Rannyu(0, 1/M_PI);
        average+= sqrt(1/(M_PI*x)-1);
        myout << average << " ";

        average = 0;

        // N = 2

        x = myGen.Rannyu(0, 1/M_PI);
        average+= sqrt(1/(M_PI*x)-1);
        x = myGen.Rannyu(0, 1/M_PI);
        average+= sqrt(1/(M_PI*x)-1);

        average/=(2);

        myout << average << " ";

        // N = 10

        average = 0;
        for (int i = 0; i <10 ; i++)
        {
            x = myGen.Rannyu(0, 1/M_PI);
            average+= sqrt(1/(M_PI*x)-1);
            
        }
        average/=10;
        myout << average << " ";

        // N = 100

        average = 0;

        for (int i = 0; i<100; i++)
        {
            x = myGen.Rannyu(0, 1/M_PI);
            average+= sqrt(1/(M_PI*x)-1);
            
        }
        average/=100;
        myout << average << endl;
    
    }

    myout.close();
    


    return 0;
}
