//A program to simulate a 3D random walk with a uniform random continuous angle chosen at each step

#include "random.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

using namespace std;

void errors_with_matrix(const vector<vector<double>> mat, vector<double> &stds)
    {
        double progr;
        double progr2;
        for(int i = 1; i<100; i++)      //for each step
        {
            progr = 0;
            progr2 = 0;
            for (int j = 0; j< 100; j++)  //for each block num
            {
                progr+=sqrt(mat[j][i]);
                progr2+=mat[j][i];
            }

            progr/=100;
            progr2/=100;
            

            stds[i] = (progr2- pow(progr,2)); 
            stds[i] = sqrt(stds[i]/99);            //error for each step after 100 blocks
        }
    }



struct position
{
    double x , y , z ;
    
    void print_pos()
    {
        cout <<  " " << x << " " << y << " " << z << endl;
    }

    double norm2()
    {
        return double(x*x + y*y + z*z);
    }

};

void do_step(Random&, vector<position>&, int);

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

    //double spacing = 1.0;
    int n_steps = 100;          // steps in each random walk
    int n_blocks = 100;         // number of blocks
    int n_perblock = 100;       // random walks per block

    
    vector <double> stds(n_steps);

    vector<vector<double>> my_matrix(n_blocks);


    for (int i = 0; i< n_blocks; i++)
    {
        vector<double> average_dist(n_blocks);

        for(int j = 0; j<n_perblock; j++)
        {
            vector<position> positions(n_steps);    //reset initial position to origin
            for (int q = 0; q<n_steps; q++)
            {
                positions[q].x = 0;
                positions[q].y = 0;
                positions[q].z = 0;
            }

            for (int k = 0; k<n_steps; k++)         // do the random walk
            {
                do_step(myGen, positions, k);

                average_dist[k]+= positions[k].norm2(); //average dist in a block after k steps
                //if(i == 0 && j== 0)
                    //cout << positions[k].x << positions[k].y << positions[k].z << endl;
            }
                  
        }
        for (int k = 0; k< 100; k++)
        {
            average_dist[k]/=100; 
            
        }
                

        my_matrix[i] = average_dist;


    }

    vector<double> cumul_average_per_step(n_blocks);

    for(int j = 0; j<n_steps; j++)
    {
        for (int i = 0; i<n_blocks; i++)
        {
            cumul_average_per_step[j]+= my_matrix[i][j];
        }

        cumul_average_per_step[j]/=n_blocks;        
        cumul_average_per_step[j] = sqrt(cumul_average_per_step[j]);    //average distances after j steps after 100 blocks
        
    }


    errors_with_matrix(my_matrix, stds);

    ofstream myout;
    myout.open("results_100_blocks.txt");

    for (int i = 0; i< n_steps; i++)
        myout << cumul_average_per_step[i] << " " << stds[i] << endl;

    myout.close();

    return 0;
}


void do_step(Random & myGen, vector<position> &positions, int i) 
{
    
    double rho = 1.0;         // radius of the sphere we can move on for each step
    double theta, phi;

    theta = myGen.Rannyu(0, M_PI);
    phi = myGen.Rannyu(0, 2*M_PI);

    if (i != 0)
        positions[i] = positions[i-1];
    
    positions[i].x+= rho*sin(theta)*cos(phi);
    positions[i].y+= rho*sin(theta)*sin(phi);
    positions[i].z+= rho*cos(theta);
    
}