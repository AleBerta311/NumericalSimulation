#include "Path.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <math.h>

Path::Path(){m_sequence.clear(); m_sequence_ind.clear(); m_size = 0;}
Path::~Path(){}


void Path::Set_size(int size)
{
    if (size >0)
        m_size = size;
    else
        cerr << "Invalid size" << endl;
}

int Path::Get_size() const
{
    return m_size;
}

City Path::Get_city(int position)
{
    if (position >=0)
        return m_sequence[position];
    else
    {
        throw out_of_range("Invalid position");
    }
}

int Path::Get_city_position(const City city) 
{
    auto it = find(m_sequence.begin(), m_sequence.end(), city);

    if (it != m_sequence.end())
        return distance(m_sequence.begin(), it);
    else
        return -1;
}

int Path::Get_city_position(int index)
{
    auto it = find(m_sequence_ind.begin(), m_sequence_ind.end(), index);

    if (it != m_sequence_ind.end())
        return distance(m_sequence_ind.begin(), it);
    else
        return -2;
}

City Path::Get_city_from_ind(int index)
{
    auto it = find(m_sequence_ind.begin(), m_sequence_ind.end(), index);
    if (it != m_sequence_ind.end())
    {
        int loc = distance(m_sequence_ind.begin(), it);
        return m_sequence[loc];
    }
    else
    {
        cout << "\n\tCan't do index " << index << endl;
        throw out_of_range("Cannot get city from invalid index");
    }
       
        
}

vector<City> Path:: Get_sequence() const
{
    return m_sequence;
}

vector<int> Path:: Get_sequence_ind() const
{
    return m_sequence_ind;
}

void Path::Set_Index_in_pos(int pos, int ind)
{
    if (pos >= 0 && pos <m_size && ind >=0 && ind<m_size)
    {
        m_sequence_ind[pos] = ind;
    }
    else
        throw out_of_range("Cannot assign index to desired position");
}

double Path:: Get_Length() const
{
    double cost = 0;
    for (int i = 0; i < m_size - 1; i++) {
        double distance = m_sequence[i].GetDistance(m_sequence[i + 1]);
        if (std::isnan(distance) || std::isinf(distance)) {
            std::cerr << "Invalid distance between city " << i << " and city " << (i + 1) << ": " << distance << std::endl;
            throw std::runtime_error("Invalid distance encountered");
        }
        cost += distance;
    }

    double wrap_around_distance = m_sequence[m_size - 1].GetDistance(m_sequence[0]);
    if (std::isnan(wrap_around_distance) || std::isinf(wrap_around_distance)) {
        std::cerr << "Invalid wrap-around distance between city " << (m_size - 1) << " and city 0: " << wrap_around_distance << std::endl;
        throw std::runtime_error("Invalid wrap-around distance encountered");
    }
    cost += wrap_around_distance;

    if (std::isnan(cost) || std::isinf(cost)) {
        std::cerr << "Invalid total cost calculated: " << cost << std::endl;
        throw std::runtime_error("Invalid total cost encountered");
    }

    return cost;
}

bool Path::IsPresent(double x, double y)
{
    bool present = false;
    for (int i = 0; i< m_size; i++)
    {
        if (m_sequence[i].Get_x() == x && m_sequence[i].Get_y() == y)
        {
            present = true;
            break;
        }
        else present = false;
    }

    return present;
}
bool Path::IsPresent(int index)
{
    bool present = false;
    for (int i = 0; i< m_size; i++)
    {
        if (m_sequence_ind[i] == index)
        {
            present = true;
            break;
        }
        else present = false;
    }

    return present;
}

void Path::add_City(const City new_city)
{
    double new_x= new_city.Get_x();
    double new_y= new_city.Get_y();
    if (!IsPresent(new_x, new_y) || !IsPresent(new_city.Get_index()))
    {
        m_sequence.push_back(new_city);
        m_sequence_ind.push_back(new_city.Get_index());
        m_size++;
    }
    
}

void Path::write_Path_to_file(const char* filename)
{
    ofstream coutf;
    coutf.open(filename);
    coutf << "# Index:        x_coord:          y_coord: " << endl;
    for (int i = 0 ; i< m_size; i++)
    {
        coutf << setw(12) << m_sequence_ind[i] <<
                 setw(12) << m_sequence[i].Get_x() <<
                 setw(12) << m_sequence[i].Get_y() << endl;
    }
    coutf.close();
}

bool Path::quick_check()
{
    int checksum = 0;
    for (int i = 0; i< m_size; i++)
        checksum+=m_sequence_ind[i];
    int correct_sum = (m_size-1)*floor(m_size/2);

    if(checksum == correct_sum)
        return true;
    else
        return false;

}

void Path::Set_ID(int id)
{
    if(id >= 0)
        m_ID = id;
    else
        throw out_of_range("Invalid ID number");
}

int Path::Get_ID() const
{
    return m_ID;
}

void Path::Pair_permutation(int pos_a, int pos_b)
{
    if (pos_a >= 0 && pos_a < m_sequence.size() && pos_b >= 0 && pos_b < m_sequence.size()) 
    {
        City first_city  = m_sequence[pos_a];
        City second_city  = m_sequence[pos_b];
        m_sequence[pos_a] = second_city;
        m_sequence[pos_b] = first_city;
        m_sequence_ind[pos_a] = second_city.Get_index();
        m_sequence_ind[pos_b] = first_city.Get_index();
    }
    else
        throw out_of_range("Index out of range in Pair Permutation");
    
}

void Path::Inversion(int pos_in, int pos_fin)
{
    if (pos_in >= 0 && pos_fin < m_sequence.size() && pos_in < pos_fin) 
    {
        //cout << "Inverting positions from " << pos_in << " to " << pos_fin << endl;
        while (pos_in < pos_fin) 
        {
            swap(m_sequence[pos_in], m_sequence[pos_fin]);
            swap(m_sequence_ind[pos_in], m_sequence_ind[pos_fin]);
            pos_in++;
            pos_fin--;
        }
    } 
    else 
    {
        std::cerr << "Invalid indices in Inversion: pos_in = " << pos_in << ", pos_fin = " << pos_fin << std::endl;
        throw out_of_range("Index out of range in Inversion");
    }
        
}

void Path::Contiguous_permutation(int pos_a, int num, int pos_b)
{
    // Check if positions and lengths are valid
    if (pos_a >= 0 && pos_b >= 0 &&
        pos_a + num < m_sequence.size() && 
        pos_b + num <= m_sequence.size() &&
        pos_a + num < pos_b) 
    {
        
        // Temporary storage for the elements to be swapped
        vector<City> to_exchange(m_sequence.begin() + pos_a, m_sequence.begin() + pos_a + num);
        vector<City> to_exchange_back(m_sequence.begin() + pos_b, m_sequence.begin() + pos_b + num);

        // Perform the swap
        for (int i = 0; i < num; i++) 
        {
            m_sequence[pos_b + i] = to_exchange[i];
            m_sequence[pos_a + i] = to_exchange_back[i];
        }

        this->Repair_Indices();
    } 
    else 
        throw out_of_range("Invalid indices or blocks too close in Contiguous_permutation");
    
    
}
void Path::Repair_Indices()
{
     // Clear the existing indices
    m_sequence_ind.clear();

    // Reconstruct the indices based on the current sequence
    for (const auto& city : m_sequence) {
        m_sequence_ind.push_back(city.Get_index());
    }       
    
}