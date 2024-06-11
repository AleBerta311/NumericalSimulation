#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include "random.h"
#include "city.h"


using namespace std;

City :: City(){}
City :: ~City(){}

void City::Generate_City(Random& rng, string figure_type, double lin)       // Generate the coordinates based on the type of "enclosure"
{
    if(figure_type == "Circle")     // Generate a city ON a circumference
    {
        // lin here is the radius of the circle
        double theta = rng.Rannyu(0, 2*M_PI);
        m_x = lin*cos(theta);
        m_y = lin*sin(theta);
    }
    else if (figure_type == "Square")    // Generate a city INSIDE the square
    {
        // lin here is the side of the square

        m_x = rng.Rannyu(-lin, lin);
        m_y = rng.Rannyu(-lin, lin);

    }
    else
        cerr << "Invalid input" << endl;

    // Validate coordinates
    if (std::isnan(m_x) || std::isnan(m_y) || std::isinf(m_x) || std::isinf(m_y)) 
    {
        std::cerr << "Invalid coordinates generated: (" << m_x << ", " << m_y << ")" << std::endl;
        throw std::runtime_error("Invalid coordinates generated");
    }
}

double City::Get_x() const
{
    return m_x;
}
double City::Get_y() const
{
    return m_y;
}

int City::Get_index() const
{
    return m_index;
}

double City::GetDistance(City other) const
{
    double dx = m_x - other.Get_x();
    double dy = m_y - other.Get_y();
    double distance = sqrt(dx * dx + dy * dy);

    // Check for invalid distance values
    if (std::isnan(distance) || std::isinf(distance)) {
        std::cerr << "Invalid distance calculation between (" << m_x << ", " << m_y 
                  << ") and (" << other.Get_x() << ", " << other.Get_y() << "): " << distance 
                  << " between city indices " << m_index  << " and " << other.Get_index() <<endl;
                   throw std::runtime_error("Invalid distance encountered");
    }

    return distance;
}

bool City::Is_Visited()
{
    return m_visited;
}

void City::SetIndex(int ind)
{
    if(ind >= 0)
        m_index = ind;
    else 
        cerr << "Invalid index" << ind << endl;
}

void City::SetVisited(bool vis)
{
    m_visited = vis;
}
