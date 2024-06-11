#include "Path.h"     //already includes city.h and random.h
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "Population.h"
#include "mpi.h"

using namespace std;

void Initialize_Rannyu(Random&, int); //necessary
void Initialize_Rannyu(Random&); //necessary
const vector<City> Create_city_list(Random&, int, string, double);   //creates the necessary amount of cities with unique positions and identifiers
Path Generate_Path(vector<City>, Random&, int);
Path Deserialize_index_list(vector<int>&, const vector<City>&, Random&);
Path Olympics(int, int, Population&, const vector<City>&, Random&);       // gets top individ from processes, compares and sends the best back to all


int main(int argc, char* argv[])      
{
    if (argc< 1)
    {
        cerr << "Wrong or insufficient number of arguments passed to .exe" << endl;
        return -404;

    }
        
    else            // the actual program
    {
            int size, rank;
            MPI_Init(&argc, &argv);         //Stuff needed for parallel processing
            MPI_Comm_size(MPI_COMM_WORLD, &size);   
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);       //homeboy needs to get its own rank 

            Random rng;
            Initialize_Rannyu(rng);     //first they all need the same list generated
            
            int n_cities = 34;
            string figure = "Square";

            double lin = 1;     //character. length of the problem
            bool convergence = false;

            // Gonna have to generate the cities' list with their positions and assign index
            vector<City> city_list = Create_city_list(rng, n_cities, figure, lin);

            Initialize_Rannyu(rng, rank);  // now we need each MPI process to re-initialize with different seed
            
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

            string filename = "Gen_Length" + to_string(rank) + ".dat";
            string best_averages = "Best_Averages" + to_string(rank) + ".dat";
            string best_result = "Best_path" + to_string(rank) + ".dat";

            my_pop.Write_ordered_costs(filename);
            
            // ------------------------------------------Genetic Algo for optimiz ------------------------------------------------

            bool unchanged_generations = false;
            // Define a threshold for considering lengths as unchanged (e.g., tolerance for change)
            double length_tolerance = 0.1; 
            
            ofstream coutf;
            coutf.open(best_averages);
            coutf << " #    Gen: " << "     Top 50% av." << endl;
            coutf.close();
            int number_of_gens = 150;       //numbers of generations to make (the original one is 0th)
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

                Path new_champion = Olympics(size, rank, my_pop, city_list, rng);

                my_pop.Remove_last_individual();
                my_pop.Add_individual(new_champion);

                my_pop.Order_by_Fitness();

                //////////////////////////////////////////////////////////////////////////////////////////////

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
                
                ////////////////////////////////////////////////////////////////////////////////////////////////            

                my_pop.Write_ordered_costs(filename);
                my_pop.Write_best_half_average(best_averages, i+1);
                cout << "   ##########   Gen " << i+1 << " done     #########" << endl;

                ////////////////////////////////////////////////////////////////////////////////////////////////
                if(figure == "Circle" && my_pop.Get_best_individual().Get_Length() < 0.995*2*M_PI*lin)
                {
                    cout << "CONVERGENCE REACHED. Stopping algorithm." << endl;
                    convergence = true;
                    break;
                }
                if (figure == "Square" && i == number_of_gens )
                {
                    cout << "MAX NUM OF GENS REACHED. Stopping algorithm." << endl;
                    convergence = true;
                    break;
                }
                i++;

                /////////////////////////////////////////////////////////////////////////////////////////////////////
            } while (!convergence);

            vector <int> send_array = my_pop.Get_best_individual().Get_sequence_ind();
            vector<int> recv_array;
            if (rank == 0)
            {
                
                recv_array.resize(size*n_cities);
        
            }

            //Gather best-individual arrays to root process
            MPI_Gather(send_array.data(), n_cities, MPI_INT,
                recv_array.data(), n_cities, MPI_INT, 0, MPI_COMM_WORLD);
            
            if (rank == 0)
            {
                for (int i = 0; i < recv_array.size(); i++)
                {
                    if (i == 0 || i%n_cities == 0)
                        cout << "\n\tBest from continent:     ";
                    cout << recv_array[i] << "   ";
                }
            }
            

            my_pop.Write_best_individual(best_result);

            MPI_Finalize();
            return 0;
    }        
}

void Initialize_Rannyu(Random& rng, int rank)
{
    //We need to extract some primes and seed(s) from the files listed for the "advanced" random generator to work

    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");

    std::string dummy;
    for (int i = 0; i < rank; ++i) 
    {
        std::getline(Primes, dummy); // Just read and discard the line
    }

    if(Primes.is_open())
        Primes >> p1 >> p2;

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

    cout << "Process " << rank << " initialized with primes " << p1 << " and " << p2 << endl; 

}

const vector<City> Create_city_list(Random& rng, int n_cities, string figure, double lin) 
{
    vector<City> city_list;
    for (int i = 0; i<n_cities; i++)
            {
                City new_city;
               
                do
                {
                    new_city.Generate_City(rng, figure, lin);
                    new_city.SetIndex(i);
                    new_city.SetVisited(false);

                } while ( count(city_list.begin(), city_list.end(), new_city) > 0);  // all the cities have different coords and index

               city_list.push_back(new_city);
                
            }
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

Path Deserialize_index_list(vector<int>& index_list, const vector<City>& city_list, Random& rng)
{
    Path equivalent_path;
    
    for (int i = 0; i< index_list.size(); i++)
    {
        auto it = std::find_if(city_list.begin(), city_list.end(), [i](const City& city) 
        {
            return city.Get_index() == i; 
        });
        if (it != city_list.end())
        {
            int loc = distance(city_list.begin(), it);
            equivalent_path.add_City(city_list[loc]);
        }
        else
        {
            cout << "\n\tCan't do index " << i << endl;
            throw out_of_range("Cannot deserialize list");
        }
    }

    equivalent_path.Set_ID(300000+ rng.Rannyu(0, 10000));
    return equivalent_path;
    
}

Path Olympics(int size, int rank, Population& my_pop, const vector<City>& city_list, Random& rng)
{
    // all the processes except root send their #1 individual to root to compare them

    vector<int> send_ids;
    vector<int> recv_ids;

    int n_cities = my_pop.Get_best_individual().Get_size();

    vector <int> best_sequence(n_cities);

    if (rank == 0)
    {
        send_ids = my_pop.Get_best_individual().Get_sequence_ind();
        recv_ids.resize(size*n_cities);
    
        MPI_Gather(send_ids.data(), n_cities, MPI_INT,
                recv_ids.data(), n_cities, MPI_INT, 0, MPI_COMM_WORLD);

        
            Population champions;

            for (int i = 0; i < recv_ids.size(); i++)
            {
                if (i % n_cities == 0)
                {
                    vector<int> serialized_list(recv_ids.begin() + i, recv_ids.begin() + i + n_cities);
                    champions.Add_individual(Deserialize_index_list(serialized_list, city_list, rng));
                }       
            }

            best_sequence = champions.Get_best_individual().Get_sequence_ind();        

    }
    else
    {
        // Non-root processes send their data to the root process
        MPI_Gather(my_pop.Get_best_individual().Get_sequence_ind().data(), n_cities, MPI_INT,
                   NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }

    
    
    
    MPI_Bcast(best_sequence.data(), n_cities, MPI_INT, 0, MPI_COMM_WORLD);
    
    //Construct best path on all processes so the function can return it

    Path best_path = Deserialize_index_list(best_sequence, city_list, rng);

    return best_path;

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