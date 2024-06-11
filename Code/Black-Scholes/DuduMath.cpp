

#include "DuduMath.h"

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


void errors(const vector<double> &av, const vector<double> &av2, vector<double> &stds, int leng)
    {
        double progr;
        double progr2;
        for(int i = 1; i<leng; i++)
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

double error(double av, double av2)
    {
        return sqrt(av2 - av*av);
    }

void Initialize_Rannyu(Random& myGen)
{
    
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
}
