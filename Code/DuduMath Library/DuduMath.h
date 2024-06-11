//      This library contains a bunch of various functions related to mathematical or statistical operations and methods
//      It was developed for personal use mainly so it's not "professional"
//      Made by Alessandro Bertarelli
//---------------------------------------------------------------------------------------------------------------------------------------


//This function fills the stds vector with the uncertainties after all the blocks are calculated using the block method from a matrix mat that MUST be a vector of vectors mat[j][i]
// where mat[j] is a vector containing average values in ONE BLOCK for each step i and mat[j][i] is the average of the i-th calculation in the j-th block
void errors_with_matrix(const vector<vector<double>> mat, vector<double> &stds);    


// This function fills the stds vector with uncertainties calculated from the av and av2 vectors that are assumed to contain the average
// in each block and its square in av2.  i.e. : av[i] is the average from the i-th block
void errors(const vector<double> &av, const vector<double> &av2, vector<double> &stds);


// This function returns the error estimation as av2 - av*av
double error(double av, double av2);


//  This function simply does all the stuff needed to setup the linear random generator to work correctly
void Initialize_Rannyu(Random& myGen);