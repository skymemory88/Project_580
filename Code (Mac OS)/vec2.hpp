#ifndef VEC2_HPP
#define VEC2_HPP

#include <iostream>
using std::ostream;
using std::istream;

#include <iomanip>
using std::setw;

class vec2
{
public:
	// member data
	double x,y;

	// constructor
	vec2() : x(0), y(0) {}
	vec2(double x_, double y_) : x(x_), y(y_) {}
	
	vec2& operator+=(const vec2 &v)
	{
		x += v.x;
		y += v.y;
		return *this;
	}

	vec2& operator-=(const vec2 &v)
	{
		x -= v.x;
		y -= v.y;
		return *this;
	}
	
	double norm2(void) const
	{
		return x*x + y*y;
	}
	
	friend vec2 operator+(const vec2 &u, const vec2 &v)
	{
		return vec2(u.x+v.x,u.y+v.y);
	}

	friend vec2 operator-(const vec2 &u, const vec2 &v)
	{
		return vec2(u.x-v.x,u.y-v.y);
	}
	
	friend vec2 operator*(const vec2 &v, double scale)
	{
		return vec2(scale*v.x,scale*v.y);
	}
	
	friend vec2 operator*(double scale, const vec2 &v)
	{
		return vec2(scale*v.x,scale*v.y);
	}
    
    friend vec2 operator/(const vec2 &v, double scale)
	{
		return vec2(v.x/scale,v.y/scale);
	}

	friend double dot(const vec2 &u, const vec2 &v)
	{
		return u.x*v.x + u.y*v.y;
	}
	
	// a friend function does not belong to the class but is affiliated with it
	friend ostream& operator<<(ostream& os, const vec2 &v)
	{
		os.precision(8);
		os << setw(20) << v.x
		   << setw(20) << v.y;
		return os;
	}
	friend istream& operator>>(istream& is, vec2 &v)
	{
		is >> v.x >> v.y;
		return is;
	}
};

#endif // VEC2_HPP

