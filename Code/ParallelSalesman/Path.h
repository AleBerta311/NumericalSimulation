#ifndef __Path__
#define __Path__

#include "city.h"
#include <vector>


using namespace std;

class Path
{

    public:

        Path();
        ~Path();
        void Set_size(int);
        int Get_size()const;
        void Set_ID(int);
        int Get_ID() const;

        City Get_city(int);     // get a city based on its position in the sequence
        City Get_city_from_ind(int); // get city based on its index
        int Get_city_position(const City); // get the position within the sequence of a city
        int Get_city_position(int);   // same as previous but pass only its identifier index to search the sequence
        vector<City> Get_sequence() const; // Get the sequence of cities as cities
        vector<int> Get_sequence_ind() const;     // Get the sequence of cities listed as their identification index
        void Set_Index_in_pos(int, int);        //set the i-th index to be equal to n
        void add_City(const City);        //adds the city to the first available position in m_sequence 

        double Get_Length() const;           // Get the COST FUNCTION aka the sum of distances between cities as ordered
        
        bool IsPresent(double, double);       // are the given coords (aka city) already present in the path?
        bool IsPresent(int);                // is the given city-index present in the path (useful when making children, not before)

        void write_Path_to_file(const char*);

        bool quick_check();         //easiest check to see that all cities remain present in the path (not guarantee but quick)

        const bool operator==(const Path& other) const
        {
            return (m_sequence == other.Get_sequence());
        };

        // ------------------These are the possible mutations----------------------

        void Pair_permutation(int, int);        //swap two cities based on a selection of position in the sequence
        void Inversion(int, int);            // take a sub-sequence from pos a to pos b and invert it
        void Contiguous_permutation(int, int, int); //swap a number m of contiguous cities from an initial pos with the same num from a final pos
        void Repair_Indices();          // THIS MUST BE EXECUTED AFTER APPLYING MUTATIONS

    private:

        int m_size;     //number of cities in the path
        vector<City> m_sequence;        // This sequence represent the order of travel for the salesman
        vector<int> m_sequence_ind;      // the same as previous but only contains the identification indices
        int m_ID;       // to distinguish paths within the population

};


#endif 