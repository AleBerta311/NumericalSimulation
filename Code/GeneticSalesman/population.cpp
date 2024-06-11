#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <algorithm>
#include "Population.h"


using namespace std;

Population :: Population(){m_pop.clear(); m_pop_id.clear(); m_pop_count = 0; m_IsSorted = false;}
Population :: ~Population(){}

void Population::Add_individual(const Path& path)
{
    m_pop.push_back(path);
    m_pop_id.push_back(path.Get_ID());
    m_pop_count++;
    chromosome_size = path.Get_size();
}      
int Population::Get_path_position(const Path& path)
{
    auto it = find(m_pop.begin(), m_pop.end(), path);

    if (it != m_pop.end())
        return distance(m_pop.begin(), it);
    else
        return -10;
}
vector<Path> Population::Get_population() const
{
    return m_pop;
}  
vector<int> Population::Get_population_ids () const
{
    return m_pop_id;
} 

int Population::Get_Pop_Count()
{
    return m_pop_count;
}

Path Population::Get_path(int pos)const
{
    if(pos <0)
        throw out_of_range("Invalid individual position in m_pop");
    else
        return m_pop[pos];
}

void Population::Order_by_Fitness()
{
    cout << "\n\n\t Evaluating new pop's fitness..." << endl << endl;;
    sort(m_pop.begin(), m_pop.end(), &Population::compare);
    for (int i = 0; i< m_pop.size(); i++)
    {
        m_pop_id[i] = m_pop[i].Get_ID();
        //cout << "\n\tOrdered ID: " << m_pop[i].Get_ID() << endl;
    }

    m_IsSorted = true;
        
}

void Population::Write_ordered_costs(string filename)
{
    ofstream coutf;
    coutf.open(filename);
    //coutf << "  Length:   " << endl;
    for (int i = 0 ; i< m_pop_count; i++)
    {
        coutf << setw(12) << m_pop[i].Get_Length() << endl;
    }
    coutf.close();
}

parent_pair Population::Select_Parents(Random& rng)
{
    double rand = rng.Rannyu(0,1);
    int parent_1 = static_cast<int>(m_pop_count*pow(m_selection_curve_steep, rand) );
    int parent_2;
    do
    {   
        rand = rng.Rannyu(0,1);
        parent_2 = static_cast<int>(m_pop_count*pow(m_selection_curve_steep, rand) );
    } while (parent_2 == parent_1 || abs(m_pop[parent_1].Get_Length() == m_pop[parent_2].Get_Length()) );
    double length_1 =m_pop[parent_1].Get_Length();
    double length_2 = m_pop[parent_2].Get_Length();
    /* Debugging output to log the selected parent indices
    cout << "Selected parent indices: " << parent_1 << ", " << parent_2  
         << "lenghts: " << length_1 << " , " << length_2 << endl;
    */
    parent_pair parents;
    parents.first = parent_1;
    parents.second = parent_2;

    return parents;
    
}

children Population::Make_Kids(const parent_pair parents, int gen, int couple)
{
    Path father = m_pop[parents.first];
    Path mother = m_pop[parents.second];

    m_crossover_threshold = static_cast<int>(father.Get_size()*0.5);

    children kiddos;

    //-------------------------------------
    for (int i = 0; i< m_crossover_threshold; i++)
    {
        
        kiddos.first.add_City(father.Get_city(i));
        kiddos.second.add_City(mother.Get_city(i));
    }
    
    //This is the not so easy part

    vector<int> missing_from_first;
    vector<int> missing_from_second;

    /*std::cout << "Debug: Checking missing cities for kiddo 1" << std::endl;
    cout << "Debug: Father size: " << father.Get_size() << endl;
    cout << "Debug: Mother size: " << mother.Get_size() << endl;
    */
    for (int i = 0; i< father.Get_size(); i++)
    {
        //std::cout << "Debug: Checking city at index " << i << std::endl;
        if(!kiddos.first.IsPresent(i))
        {
            missing_from_first.push_back(i);
            //std::cout << "Debug: City at index " << i << " not present in kiddo 1" << std::endl;
        }
        

        if(!kiddos.second.IsPresent(i))
            missing_from_second.push_back(i);
    } 

    /* Debugging output to log the positions of missing cities
    cout << "Missing cities from first offspring: ";
    for (int i : missing_from_first) {
        cout << i << " ";
    }
    cout << endl;

    cout << "Missing cities from second offspring: ";
    for (int i : missing_from_second) {
        cout << i << " ";
    }
    cout << endl;
    */

    sort(missing_from_first.begin(), missing_from_first.end(), [&mother](int a, int b)
    {
        return mother.Get_city_position(a) < mother.Get_city_position(b);
    });

   
    sort(missing_from_second.begin(), missing_from_second.end(), [&father](int a, int b)
    {
        return father.Get_city_position(a) < father.Get_city_position(b);
    });

    
    /* Debugging output to log the sorted positions of missing cities
    cout << "Sorted positions of missing cities from first offspring: ";
    for (int i : missing_from_first) {
        cout << i << " ";
    }
    cout << endl;

    cout << "Sorted positions of missing cities from second offspring: ";
    for (int i : missing_from_second) {
        cout << i << " ";
    }
    cout << endl;
    */

    

    for (int j = 0; j< missing_from_first.size(); j++)
        kiddos.first.add_City(mother.Get_city_from_ind(missing_from_first[j]));
    for (int j = 0; j< missing_from_second.size(); j++)
        kiddos.second.add_City(father.Get_city_from_ind(missing_from_second[j]));

    kiddos.first.Set_ID(m_pop_count + (gen+1)*(couple*2));
    kiddos.second.Set_ID(m_pop_count + (gen+1)*(couple*2)+1);

    //cout << "\nDebug: kids' size: " << kiddos.first.Get_size()
        // << " and " << kiddos.second.Get_size() << endl << "--------------------" << endl;
    
    
    return kiddos;

}

void Population::Apply_Mutations(Random& rng)
{

    cout << "\n\n\tApplying mutations.." << endl << endl;
    for (int i = 0; i < m_pop_count; i++)       // for each member decide whether to apply various mutations
    {
        double rand1, rand2, rand3;
        rand1 = rng.Rannyu(0,1);
        rand2 = rng.Rannyu(0,1);
        rand3 = rng.Rannyu(0,1);

        if(rand1 < m_pair_permutation_chance)
        {
            int pos1 = static_cast<int>(rng.Rannyu(1, chromosome_size));
            int pos2;
            do
            {
                pos2 = static_cast<int>(rng.Rannyu(1, chromosome_size));
            }while (pos2 == pos1);

            //cout << "Pair Permutation: pos1 = " << pos1 << ", pos2 = " << pos2 << endl;
            m_pop[i].Pair_permutation(pos1,pos2);
            
        }
        if(rand2 < m_group_invert_chance)
        {
            int pos1 = static_cast<int>(rng.Rannyu(1, chromosome_size-1));
            int pos2 = static_cast<int>(rng.Rannyu(pos1+1, chromosome_size));


            //cout << "Inversion: pos1 = " << pos1 << ", pos2 = " << pos2 << endl;
            m_pop[i].Inversion(pos1, pos2);

        }
        if (rand3 < m_shift_contig_chance)
        {
            int pos1 = static_cast<int>(rng.Rannyu(1, chromosome_size/2));
            
             
            int num = static_cast<int>(rng.Rannyu(2, 6));
            
            int pos2 = static_cast<int>(rng.Rannyu(pos1+num+1, chromosome_size-num));
            
            
            //cout << "chromosize:" << chromosome_size << endl;
            //cout << "Contiguous Permutation: pos1 = " << pos1 << ", num = " << num << ", pos2 = " << pos2 << endl;
            m_pop[i].Contiguous_permutation(pos1, num, pos2);
        }
    }
    
}

void Population::Create_new_gen(const Population& new_pop)
{
    cout << "\n\n\tSubstituting old pop..." << endl << endl;;
    
    // Print sizes before any operation
    //std::cout << "Before clear: m_pop size: " << m_pop.size() << ", m_pop_id size: " << m_pop_id.size() << std::endl;
    //std::cout << "New population size: " << new_pop.Get_population().size() << ", New population IDs size: " << new_pop.Get_population_ids().size() << std::endl;

    // Clear the current population
    m_pop.clear();
    m_pop_id.clear();
    
    // Get population data from new_pop
    vector<Path> new_population = new_pop.Get_population();
    vector<int> new_population_ids = new_pop.Get_population_ids();
    
    // Copy population data manually
    for (const auto& path : new_population) {
        m_pop.push_back(path);
    }
    
    for (int id : new_population_ids) {
        m_pop_id.push_back(id);
    }

    // Verify data integrity by comparing IDs
    for (size_t i = 0; i < m_pop.size(); ++i) 
    {
        if (m_pop[i].Get_ID() != m_pop_id[i]) 
        {
            cerr << "Error: ID mismatch for path at index " << i << endl;
        }
    }
}

void Population::Write_best_half_average(string filename, int gen)
{
    if (m_IsSorted)
    {
        ofstream coutf;
        coutf.open(filename, ios::app);

        double accum = 0;
        for (int i = m_pop_count/2 ; i < m_pop_count; i++)
        {
            accum+= m_pop[i].Get_Length();
        }

        accum/= double(m_pop_count/2);

        coutf << setw(12) << gen 
              << setw(12) << accum << endl;

        coutf.close();
    }
    else
    {
        Order_by_Fitness();
        Write_best_half_average(filename, gen);
    }
        
}

void Population::Write_best_individual(string filename)
{
    if (m_IsSorted)
    {
        Path best_one = m_pop[m_pop_count-1];
        ofstream coutf;
        coutf.open(filename);
        coutf << "#    Position:        Index:           x_coord:          y_coord:" << endl;
        for (int i = 0; i<chromosome_size; i++)
        {
            City current_city = best_one.Get_city(i);

            coutf << setw(12) << i 
                  << setw(12) << current_city.Get_index()
                  << setw(12) << current_city.Get_x()
                  << setw(12) << current_city.Get_y() << endl;
        }

        coutf.close();
    }
    else
    {
        Order_by_Fitness();
        Write_best_individual(filename);    
    }
}

Path Population::Get_best_individual()
{
    if (m_IsSorted)
        return m_pop[m_pop_count-1];
    else
    {
        Order_by_Fitness();
        return Get_best_individual();
    }
}