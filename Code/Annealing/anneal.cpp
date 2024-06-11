#include "random.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <iomanip>

using namespace std;

void InitializeRannyu(Random& );

double PsiTrial(double, double, double);            // Wave function with parameters to be optimized
double Psi2Mod(double, double, double);            // Square modulus of wave function

double Metropolis(double, double, Random&, double, double);   //This returns the n+1 x position after metropolis accept/reject using a uniform transition step
double Hamilt(double, double, double);
double Potential(double);
double Kinetic(double, double, double);
double error (vector<double>&, vector<double>&, int);
double global_av(vector<double>&, int);
double energy_calc(vector<double>&, vector<double>&, vector<double>&, double, double, int, int, int, Random&, bool);
double Boltz(double, double, double);


int main()
{
    Random rng;

    int n_eval = 100000;
    int block_length = 2500;

    InitializeRannyu(rng);
    double mu = 1.1, sigma = 0.5, old_mu, old_sigma;       // Initial guess + variables for SA
    double step = 1.8;    // step for metropolis sampling of psi
    double anneal_step_mu = 0.5;       // step for updating mu and sigma in SA
    double anneal_step_sigma = 0.3;    

    int n_temp = 40;    // num of temps to try for annealing
    double T , beta = 0.1, cooling_rate = 0.91, T_0 = 1.4;     
    vector<int> steps_per_temp;     //how many steps to do at each T

    vector<double> energies;            // These are for progr.av eval for a single mu, sigma
    vector <double> en2;
    vector<double> errors;

    vector <double> final_energies;        // These are for the final evals after all the Psi blocks
    vector <double> final_errors;

    double old_energy, new_energy = 1, old_error;       // stuff needed for SA
    bool accept = true;
    
    ofstream coutf;
    coutf.open("<H>_estim.dat");
    coutf << "      T:         <H>:          err:          mu:           sigma:        acc?"  << endl;
    coutf.close();
    
    vector <double> average_en_per_temp;

    bool to_print = false;
    
    for (int k = 0; k<n_temp; k++)
    {

        anneal_step_mu*=0.99;
        anneal_step_sigma*=0.99;
        double accum = 0; 
        T = pow(cooling_rate,k)*T_0;
        steps_per_temp.push_back(floor(T*100));

        final_energies.clear();
        final_errors.clear();
        for (int n = 0; n < steps_per_temp[k]; n++)
        {
            ////////////  This is the energy sampling for a given (mu, sigma) config /////////////////
            if(accept)
            {
                energies.clear();
                en2.clear();
                errors.clear();
                old_energy = energy_calc(energies, en2, errors, mu, sigma, n_eval, block_length, step, rng, to_print);
                old_error = errors[errors.size()-1];
            }
            
            energies.clear();
            en2.clear();
            errors.clear();
            
            
            /////// Here we try to update mu and sigma for a given temp /////
            old_mu = mu;
            old_sigma = sigma;

            mu = rng.Rannyu(mu-anneal_step_mu, mu+anneal_step_mu);
            sigma = rng.Rannyu(sigma-anneal_step_sigma, sigma+anneal_step_sigma);

            new_energy = energy_calc(energies, en2, errors, mu, sigma, n_eval, block_length, step, rng, to_print);

            if(new_energy < old_energy)
            {
                final_energies.push_back(new_energy);
                final_errors.push_back(errors[errors.size()-1]);
                accept = true;
            }
            else if (new_energy >= old_energy)
            {
                double alpha = min(1.0, Boltz(old_energy, new_energy, T));
                double r = rng.Rannyu(0,1);

                if (r < alpha)
                {
                    final_energies.push_back(new_energy);
                    final_errors.push_back(errors[errors.size()-1]);
                    accept = true;
                    // mu and sigma are left as the new ones
                }

                else
                {
                    mu = old_mu;
                    sigma = old_sigma;
                    final_energies.push_back(old_energy);
                    final_errors.push_back(old_error);
                    accept = false;

                }
            }  
            
             
            
            coutf.open("<H>_estim.dat", ios::app);
            coutf << std::fixed << std::setprecision(6);
            coutf << setw(12) << T
                  << setw(12) << final_energies[n]
                  << setw(12) << final_errors[n]
                  << setw(12) << mu
                  << setw(12) << sigma 
                  << setw(12) << accept << '\n';
            coutf.close();

            accum+= final_energies[n];
            
        }
            
        // Then we move on to a new temperature 
        average_en_per_temp.push_back(accum/steps_per_temp[k]);       
        
    }
    
    coutf.open("Averages_per_temp.dat");
    for (int i = 0; i< n_temp; i++)
    {
        coutf << std::fixed << std::setprecision(4);
        coutf << setw(12) << pow(cooling_rate,i)*T_0
              << setw(12) << average_en_per_temp[i] << endl;
    }
    coutf.close();

    to_print = true;
    energies.clear();
    en2.clear();
    errors.clear();
    double store = energy_calc(energies, en2, errors, mu, sigma, n_eval, block_length, step, rng, to_print);

    return 0;
}

void InitializeRannyu(Random& rng)
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

        rng.SetRandom(seed, p1, p2);
        input.close();
    }
    else cerr << "Cannot open seed.in" << endl;

}

double PsiTrial(double x, double mu, double sigma)
{
    return exp(-pow(x-mu, 2)/(2*pow(sigma,2))) + exp(-pow(x+mu, 2)/(2*pow(sigma,2)));
}

double Psi2Mod(double x, double mu, double sigma)
{
    return pow( fabs(PsiTrial(x, mu, sigma)) , 2);
}

double Metropolis(double x_curr, double step, Random& rng, double mu, double sigma)
{
    double x_new = rng.Rannyu(x_curr-step, x_curr+step);

    bool acceptance = false;

    double alpha = min(1.0, Psi2Mod(x_new, mu, sigma)/Psi2Mod(x_curr, mu, sigma) );

    double accept = rng.Rannyu(0,1);

    if (accept < alpha)
        acceptance = true;

    if (acceptance)
        return x_new;
    else
        return x_curr;
}

double Hamilt (double x, double mu, double sigma)
{
    return Kinetic (x, mu, sigma) +Potential(x) ; 
}

double Kinetic(double x, double mu, double sigma)
{
    double dd2x =  1.0/pow(sigma,4) * (pow(x-mu,2) * exp(-pow(x-mu,2)/(2.0*sigma*sigma)) + pow(x+mu,2) * exp(-pow(x+mu,2)/(2.0*sigma*sigma))) -1.0/pow(sigma,2) * PsiTrial(x,mu,sigma);
    return -0.5 * dd2x / PsiTrial(x,mu,sigma);
}

double Potential(double x)
{
    return pow(x, 4) -2.5 * pow(x, 2);
    
}

double error (vector<double>& acc, vector<double>& acc2, int blk)
{
    double sum_prog = 0, sum2= 0;
    for (int i = 0; i< blk; i++)
    {
        sum_prog+= acc[i];
        sum2+= acc2[i];
    }

    sum_prog/=blk;
    sum2/=blk;

    return sqrt(fabs((sum2 - pow(sum_prog,2)))/blk);
}

double global_av(vector<double>& acc, int blk)
{
    double sum = 0;
    for (int i = 0; i< blk; i++)
        sum+= acc[i];
    
    sum/=blk;
    
    return sum;
}

double energy_calc(vector<double>& energies, vector<double>& en2, vector<double>& errors, double mu, double sigma, int n_eval, int block_length, int step, Random& rng, bool to_print)
{
            ofstream coutf;
            int n_launch, n_accept;
            int n_blocks = n_eval/block_length;
            double x = mu, x_old;
            double energy_est;

            for(int i = 0; i< n_blocks; i++)
            {
                n_launch = 0, n_accept = 0;
                energy_est = 0;
                
                for (int j = 0; j< block_length; j++)
                {
                    energy_est+=Hamilt(x, mu, sigma);
                    x_old = x; 
                    x = Metropolis(x, step, rng, mu, sigma);
                    n_launch++;
                    if (x_old != x)
                        n_accept++;
                }
            ////////////////////////// This is where the block results are processed and printed ///////////////////

                energies.push_back(energy_est/block_length);        //current block estim.
                en2.push_back(pow(energies[i],2));
                errors.push_back(error(energies,en2,i+1));
                
                if(to_print)
                {
                    coutf.open("<H>_calc.dat", ios::app);
                    coutf << setw(12) << i+1            //block num
                        << setw(12) << energies[i]      // block av.
                        << setw(12) << global_av(energies, i+1)     //rolling av
                        << setw(12) << errors[i] << endl;       //error on rolling av
                    coutf.close();
                }
                
                // cout << double(n_accept)/double(n_launch) << endl;      //gotta check that acceptance
            }        
            /////////////////    End of the config sampling     /////////
            return global_av(energies,energies.size());
}

double Boltz(double old, double nov, double temp)
{
    double beta = 1/temp;

    return exp(-beta * (nov-old));
}
