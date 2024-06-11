#ifndef __City__
#define __City__

#include <string>
#include "random.h"

using namespace std;


class City
{
    public:
        City();
        ~City();
        void Generate_City(Random&, string, double);        // Will take an rng in the main to generate a city's random coords on a figure (square/circ)
        void Generate_City(double, double, int);        // ONLY IMPLEMENTED WITH LATITUDE/LONGITUDE
        void SetIndex(int);                    // Will be used to set the index within a path
        void SetVisited(bool);

        double Get_x() const;
        double Get_y()const;
        double Get_lat() const;
        double Get_lon()const;
        int Get_index()const;
        bool Is_Visited();

        const bool operator==(const City& other) const
        {
            return (m_x == other.Get_x() && m_y == other.Get_y());
        };

        double GetDistance(City)const;       // distance from another city

        double GetDistanceSphere(City) const;


    private:
        double m_x, m_y;    // coordinates of the city
        double m_lon, m_lat; // coordinates but on a sphere because you know..
        int m_index;        // identifier of a city (e.g. Los Angeles, Paris etc)
        bool m_visited;     // has it already been visited by the salesman? y/n
        
};

#endif