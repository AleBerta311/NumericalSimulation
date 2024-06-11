//A program to simulate the Black-Scholes equations

#include "DuduMath.h"

using namespace std;

double EvaluateFinalPrice(double, double , double, double, Random&);


int main()
{
    ofstream myoutput;
    myoutput.open("call_prices.txt");

    Random myGen;

    Initialize_Rannyu(myGen);

    //first we try sampling the final asset price S(T) directly and then estimate the price after sampling S(T) a lot of times

    
    double initial_asset_price = 100.;
    double delivery_time = 1.0;     // T in equations
    double strike_price = 100.;       // K in equations
    double no_risk_interest_rate = 0.1;      // r = 10%
    double volatility = 0.25;                //sigma
    int n_eval = 1000000;
    int n_in_block = 5000;

    vector<double> final_prices(n_eval/n_in_block);         // this will contain the various sampled final call prices
    vector<double> fin_2(n_eval/n_in_block);               // same vals squared for error calculation
    vector<double> final_put_prices(n_eval/n_in_block);      // this will contain the various sampled final put prices
    vector<double> fin_put_2(n_eval/n_in_block);  
    
    double  final_asset_price;
    double exp_call_price_cumul_in_block;
    double exp_put_price_cumul_in_block;

    for (int j = 0; j< n_eval/n_in_block; j++)
    {
        exp_call_price_cumul_in_block = 0;
        exp_put_price_cumul_in_block = 0;
        for (int i = 0; i< n_in_block; i++)
        {
            final_asset_price = EvaluateFinalPrice(initial_asset_price, delivery_time, no_risk_interest_rate, volatility, myGen);
            exp_call_price_cumul_in_block+= max(0.0, final_asset_price-strike_price) * exp(-delivery_time*no_risk_interest_rate);
            exp_put_price_cumul_in_block+= max(0.0, strike_price - final_asset_price) * exp(-delivery_time*no_risk_interest_rate);
        
        }

        exp_call_price_cumul_in_block/=(n_in_block);
        exp_put_price_cumul_in_block/=(n_in_block);

        final_prices[j] = exp_call_price_cumul_in_block;
        fin_2[j] = pow(final_prices[j],2);   
        final_put_prices[j] = exp_put_price_cumul_in_block;
        fin_put_2[j] = pow(final_put_prices[j],2);   
    }
    
    
    vector <double> stds_call(n_eval/n_in_block);
    vector <double> stds_put(n_eval/n_in_block);

    errors(final_prices, fin_2, stds_call, n_eval/n_in_block);
    errors(final_put_prices, fin_put_2, stds_put, n_eval/n_in_block);

    double cumul_average;

    double cumul_average_put;

    for (int i = 0; i < n_eval/n_in_block; i++)
    {
        cumul_average = 0;
        myoutput << final_prices[i] << " " << stds_call[i] << " ";

        for(int k = 0; k < i; k++)
        {
            cumul_average+= final_prices[k];
        }
        if (i!= 0)
            cumul_average/=i;
        myoutput << cumul_average << " " << endl;
    }
    myoutput.close();

    myoutput.open("put_prices.txt");

    for (int i = 0; i < n_eval/n_in_block; i++)
    {
        cumul_average_put = 0;
        myoutput << final_put_prices[i] << " " << stds_put[i] << " ";

        for(int k = 0; k < i; k++)
        {
            cumul_average_put+= final_put_prices[k];
        }
        if (i!= 0)
            cumul_average_put/=i;
        myoutput << cumul_average_put << endl;
    }
    
    myoutput.close();

    return 0;
}

    

double EvaluateFinalPrice(double S_init, double final_T, double mu, double sigma, Random& myGen)
{
    double W = myGen.Gauss(0, final_T);       //W(t) of the GBM

    double final_price = S_init*exp((mu-0.5*sigma*sigma)*final_T + sigma* W*sqrt(final_T));

    return final_price;
}
    
    
