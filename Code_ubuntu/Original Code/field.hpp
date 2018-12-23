#ifndef FIELD_HPP
#define FIELD_HPP

#include <cassert>

#include "vec2.hpp"

#include <complex>
using std::complex;
using std::arg;
using std::conj;
using std::abs;

#include <vector>
using std::vector;

template <class Type> class field
{
public:
    size_t N;
    double dx;
    double dy;
    vector<Type> data;
public:
    field() : N(0), dx(0), dy(0), data() {}
    field(double L, size_t N_) : N(N_), dx(L/N_), dy(L/N_), data(N_*N_+1){}
    
    Type& operator[](size_t i)
    {
        return data[i];
    }
    
    const Type& operator[](size_t i) const
    {
        return data[i];
    }

    field<Type>& operator=(const field<Type> &B)
    {
        for (size_t i = 0; i < B.data.size(); ++i)
        {
            data[i] = B.data[i];
        }
        return *this;
    }
};
#endif
