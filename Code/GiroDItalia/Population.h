#ifndef _POP_
#define _POP_

#include "Path.h"

struct parent_pair
{
    int first, second;
};

struct children
{
    Path first, second;
};

class Population 
{
    public:
        Population();
        ~Population();


        void Add_individual(const Path&);      //Adds an individual path to the first available position of m_pop
        int Get_path_position(const Path&);    // get the position within m_pop of an individual (path) in the population 
        Path Get_path(int)const;                    // Get the path at the i-th position in m_pop
        vector<Path> Get_population() const;  // Get the population of paths as paths
        vector<int> Get_population_ids()const;     // Get the pop listed as their path-id
        int Get_Pop_Count();

        static bool compare (const Path& a, const Path& b)
        {
            return a.Get_Length() > b.Get_Length();
        };
        
        void Order_by_Fitness();        // orders the paths in the m_pop and m_pop_id based on cost function 

        void Write_ordered_costs(string);       // write a file containing all lenghts of the pop members in descending order

        parent_pair Select_Parents(Random&);

        children Make_Kids(const parent_pair, int, int);

        void Apply_Mutations(Random&);

        void Create_new_gen(const Population&);

        void Write_best_half_average(string, int);       // write the average length of the best (shortest) half of the pop

        void Write_best_individual(string);         // log the coords list as it appears in the fittest path
        Path Get_best_individual();
       


    private:
        int m_pop_count;
        vector<Path> m_pop;         // contains list of the population members at a given time
        vector<int> m_pop_id;       // same list but only with the individual ID
        bool m_IsSorted;            // Has it been sorted by fitness already?

        // These are for the actual genetic part
        double m_selection_curve_steep = 0.74; // in the function f = N*r^p the value of r (p is Unif in 0-1 and N is m_pop_count) Cant be too high or problem in parent selection arise
        int m_crossover_threshold;      // Where to perform the "cromosome cut" during reproduction
        double m_pair_permutation_chance = 0.25;    //this works ok
        double m_shift_contig_chance = 0.25;     // seems to work ok
        double m_group_invert_chance = 0.2;        // this works ok
        int chromosome_size;
};

#endif