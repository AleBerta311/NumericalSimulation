//A code to simulate the Buffon experiment to evaluate pi

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
 
    double L = 1.5 ;     //length of the needle in arbitrary units 
    double d = 2.0 ;     //distance between the lines in a.u.
    double table_length = 100.0 ;

    int n_throws = 10000;      //  # of times we throw the needle per evaluation of pi
    int n_blocks = 100;
    int n_perblock = 100; 
    double theta;
    double cos_theta;           // cos of the angle the needle makes w/ the perpend. to the table's length
    double lower_corner;        // lower corner of the needle
    double upper_corner_x;
    int N_hits = 0;             // #of times the needle hits a line
    double n_lines = double(table_length/d);
    double sum_prog = 0;
    double su2_prog = 0;
    double pi_block_av_2 = 0;
    double pi_av_2 = 0;
    double err = 0;


    ofstream myPi;
    myPi.open("mypi's.txt");

    for (int j = 0; j < n_blocks; j++)
    {
        pi_block_av_2 = 0;
        pi_av_2 = 0;
        for (int k = 0; k < n_perblock; k++)
        {
            N_hits =0;
            for (int i = 0; i<n_throws; i++)
            {
                lower_corner = myGen.Rannyu(0, table_length);       // Generate x0, consider y0 = 0 
                //theta = myGen.Rannyu(-M_PI_2, M_PI_2);      // NOOOO bcs you cant use pi to make pi :(
                //cos_theta = cos(theta);
                //upper_corner_x = lower_corner + L*cos_theta;
                
                // Method to generate angle without using everyone's favourite number
                bool in_circle = false;
                double x, y;
                do
                {
                   x = myGen.Rannyu(0,1);
                   y = myGen.Rannyu(-1, 1);         // These represent a random point that should be on the unit half-circle (RHS)
                   
                   if (x*x + y*y <= 1)
                    in_circle = true;
                    
                } while (!in_circle);
                
                theta = atan2(y, x);

                cos_theta = cos(theta);
                upper_corner_x = lower_corner + L*cos_theta;
            
                for (int j = 1; j <= n_lines; j++)
                {
                    if (lower_corner <= j*d && upper_corner_x >= j*d)
                        N_hits++;
                }
               
            }
            double pi_est = 2.0*L*n_throws/(N_hits*d);
            pi_block_av_2+= pi_est;
          
        }
            pi_block_av_2/= n_perblock;     //this is the pi estimate in the block
            pi_av_2+= pi_block_av_2*pi_block_av_2;
            
            
            sum_prog+= pi_block_av_2;
            su2_prog+= pi_av_2;

            

        myPi << pi_block_av_2 << " ";
        myPi << sum_prog/(j+1) << " ";
        if (j != 0 )
            err = error(sum_prog/(j+1), su2_prog/(j+1))/sqrt(j);

        myPi << err << endl;
    }

    return 0;

}