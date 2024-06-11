#include "Path.h"     //already includes city.h and random.h
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "Population.h"

using namespace std;

void Initialize_Rannyu(Random&); //necessary
const vector<City> Create_city_list();   //creates the necessary amount of cities with unique positions and identifiers
Path Generate_Path(vector<City>, Random&, int);


int main()        // 3 arguments : program_name, #_of_cities, figure_type
        {
    
            Random rng;
            Initialize_Rannyu(rng);     //self-explanatory

            double lin = 1;     //character. length of the problem
            bool convergence = false;

            // Gonna have to generate the cities' list with their positions and assign index
            vector<City> city_list = Create_city_list();
            
            //Now we gotta generate the first generation of paths

            Population my_pop;

            int n_paths = 1000; // number of  individuals of the population

            for (int i = 0; i<n_paths; i++)
            {
                Path new_path = Generate_Path(city_list, rng, i);
                if (new_path.quick_check())        //honestly should be redundant but its quick enough
                {
                    my_pop.Add_individual(new_path);
                    
                }     
                   
                else 
                    i--;
            }

           // Now we have a Gen 0 of (probably bad) paths

           // We want to order this first Gen by fitness

            my_pop.Order_by_Fitness();

            string filename = "Gen_Length.dat";
            string best_averages = "Best_Averages.dat";
            string best_result = "Best_path.dat";

            my_pop.Write_ordered_costs(filename);
            
            // ------------------------------------------Genetic Algo for optimiz ------------------------------------------------

            bool unchanged_generations = false;
            // Define a threshold for considering lengths as unchanged (e.g., tolerance for change)
            double length_tolerance = 0.1; 
            
            ofstream coutf;
            coutf.open(best_averages);
            coutf << " #    Gen: " << "     Top 50% av." << endl;
            coutf.close();
            int number_of_gens = 500;       //numbers of generations to make (the original one is 0th)
            int i = 0;
            do
            {
                
                cout << "\n\tStarting reproduction of " << i << " -th gen" << endl << endl << endl;
                // First we gotta select #individ/2 unions to be made, each producing 2 offspring, save the offspring
                Population* next_gen = new Population();
                for (int j = 0; j <n_paths/2; j++)
                {
                    parent_pair parents = my_pop.Select_Parents(rng);
                    children kids = my_pop.Make_Kids(parents, i+1, j);
                    next_gen->Add_individual(kids.first);
                    next_gen->Add_individual(kids.second);
                }
                
                //Replace the old pop with the new pop

                my_pop.Create_new_gen(*next_gen);
                delete next_gen;

                //apply random mutations to the new pop: voilà new generation!
                
                my_pop.Apply_Mutations(rng);//ok now


                my_pop.Order_by_Fitness();

                // Check if lengths remain unchanged
                double length_difference = abs(my_pop.Get_path(0).Get_Length() - my_pop.Get_path(n_paths - 1).Get_Length());
                if (length_difference < length_tolerance) 
                {
                    cout << length_difference;
                    unchanged_generations = true;
                    if (unchanged_generations) 
                    {
                        cout << "CONVERGENCE REACHED. Stopping algorithm." << endl;
                        my_pop.Write_ordered_costs(filename);
                        my_pop.Write_best_half_average(best_averages, i+1);
                        
                        cout << "Gen " << i+1 << " done" << endl;
                        convergence = true;
                        break; // Exit the loop to stop the algorithm
                    }
                } 
                else 
                    unchanged_generations = false; // Reset if lengths change
                
                            

                my_pop.Write_ordered_costs(filename);
                my_pop.Write_best_half_average(best_averages, i+1);
                cout << "   ##########   Gen " << i+1 << " done     #########" << endl;

              
                if (i == number_of_gens )
                {
                    cout << "MAX NUM OF GENS REACHED. Stopping algorithm." << endl;
                    convergence = true;
                    break;
                }
                i++;
            } while (!convergence);

            my_pop.Write_best_individual(best_result);
            return 0;
        }        

void Initialize_Rannyu(Random& rng)
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

const vector<City> Create_city_list() 
{
    vector<City> city_list;
    ifstream in("cap_prov_ita.dat");
    if (!in.is_open()) {
        cerr << "Error opening file.\n";
        return city_list;
    }
    
    double lon, lat;
    int i = 0;
    while (in >> lon >> lat) {
        City new_city;
        new_city.Generate_City(lon, lat, i);
        city_list.push_back(new_city);
        ++i;
    }
    
    in.close();
    return city_list;
}


Path Generate_Path(vector<City> city_list, Random& rng, int id)
{
    Path new_path;
    int final_size = city_list.size();

    // Always put the first city of the list as the starting point
    new_path.add_City(city_list[0]);
    city_list.erase(city_list.begin());

    // Debugging output to log the coordinates of the first city
    //cout << "First city coordinates: (" << new_path.Get_city(0).Get_x() << ", " << new_path.Get_city(0).Get_y() << ")" << endl;


    for (int i = 0; i < final_size - 1; i++) {
        // Generate a random index within the current size of city_list
        int rand_index = static_cast<int>(rng.Rannyu(0, city_list.size()));
        // Add the city at the random index to the path
        new_path.add_City(city_list[rand_index]);
        // Erase the city from city_list
        city_list.erase(city_list.begin() + rand_index);
         
        // Debugging output to log the coordinates of the newly added city
        //cout << "New city coordinates: (" << new_path.Get_city(i + 1).Get_x() << ", " << new_path.Get_city(i + 1).Get_y() << ")" << endl;
    }

    new_path.Set_ID(id);

    return new_path;
}

