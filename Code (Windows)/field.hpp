#ifndef FIELD_HPP
#define FIELD_HPP

#include <typeinfo>
#include <cstdlib>
#include <cassert>

#include "vec2.hpp"

#include <string>
using std::string;

#include <cstddef>
using std::size_t;

#include <complex>
using std::complex;
using std::real;
using std::abs;
using std::arg;

#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;

#include <fstream>
using std::ifstream;
using std::ofstream;

template <class Type> class field
{
private:
    int xpos();
    int ypos();
private:
    template<class T>
    inline int xpos(T i) { return i % N; }        //Extract position on x-axis from index
    template<class T>
    inline int ypos(T i) { return i / N; }        //Extract position on y-axis from index
    template<class T>
    inline size_t index(T i, T j) { return ((i + N) % N) + ((j + N) % N)*N; }       //Periodic boundary condition for cell update
    /*{
     if(i >= N or i < 0 or j >= N or j < 0)
     return (N2);
     else
     return i + j*N;                                 //Alternative boundary condition, hard open boundary condition
     }*/

public:
    size_t N;
    size_t N2;
    double dx;
    double dy;
    vector<Type> data;
public:
    field() : N(0), N2(0), dx(0), dy(0), data() {}
    field(double L, size_t N_) : N(N_), N2(N_*N_), dx(L/N_), dy(L/N_), data(N_*N_+1){}		//Create a discrete field of (N^2) cells of size (L/N_) with one additional dummy cell for boundary conditions.
    
    Type& operator[](size_t i)
    {
        return data[i];
    }
    
    const Type& operator[](size_t i) const
    {
        return data[i];
    }

    field<Type>& operator=(const field<Type> &B) 	//assign values from another field
    {
        for (size_t i = 0; i < B.data.size(); ++i)
        {
            data[i] = B.data[i];
        }
        return *this;
    }
    
	void report_data(const string filename)     //output data
	{
        string newname = filename + ".dat";
        ofstream fout;
        fout.precision(6);
        fout.open(newname.c_str());
        for (size_t i = 0; i < N2; ++i)
        {
            fout << xpos(i) << "\t" << ypos(i) << "\t" << data[i] << endl;
        }
        fout.close();
	}

	void map_data(const string filename)		//output data as a density map
	{
		string newname = filename + "_map.dat";
		ofstream fout;
		fout.precision(6);
		fout.open(newname.c_str());
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t i = 0; i < N; ++i)
			{
				fout << data[i + j*N] << "\t";
			}
			fout << endl;
		}
        fout.close();
	}
    
    void import_data(const string filename)		//import data from a file
    {
        string newname = filename + ".dat";
        ifstream fin(newname.c_str());
        if (!fin.is_open())
        {
            cerr << "ERROR: Can't open the file" << endl;
            exit(1);
        }
        size_t i = 0;
        while (true)
        {
            size_t x;
            size_t y;
            Type tamp;
            fin >> x >> y >> tamp;
            data[index(x,y)] = tamp;
            ++i;
            if (fin.eof())
            {
                break;
            }
        }
        fin.close();
    }
    
    void import_map(const string filename)		//import data from a density map file
    {
        string newname = filename + "_map.dat";
        ifstream fin(newname.c_str());
        if (!fin.is_open())
        {
            cerr << "ERROR: Can't open the ifle" << endl;
            exit(1);
        }
        for (size_t j = 0; j < N; ++j)
        {
            for (size_t i = 0; i < N; ++i)
            {
                fin >> data[i + j*N];
            }
        }
        fin.close();
    }
};
#endif
