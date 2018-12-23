#include <cassert>
#include "field.hpp"
#include "mtrand.hpp"   //Random number generator
#include <ctime>
#include <cstdlib>
#include <omp.h>        //Parallel computing package

#include <cstddef>
using std::size_t;

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

#include <iomanip>
using std::setw;

#include <cmath>
using std::fabs;
using std::tanh;
using std::atan2;
using std::sqrt;
using std::cos;
using std::sin;

mtrand Rand(time(0));
double damp = 1.0;                //choose damping factor for cell update
size_t N = 500;                   //choose discretizing size (see "field.h").
size_t N2 = N*N;
complex<double> ci(0.0,1.0);      //create mathematical imaginary number "i"

template<class T>
inline int xpos(T i) {return i%N;}        //Extract position on x-axis from index

template<class T>
inline int ypos(T i) {return i/N;}        //Extract position on y-axis from index

template<class T>
inline bool ODD(T x) {return (x%2 == 0? false:true);}       //Determine whether or not a number is odd

template<class T>                                    //Determine the minimum number between two values of the same type
inline T min(T x, T y){return (x>y?y:x);}

template<class T>
inline size_t index(T i, T j) { return ((i + N) % N) + ((j + N) % N)*N; }       //Periodic boundary condition for cell update
/*{
   if(i >= N or i < 0 or j >= N or j < 0)
    return (N2);
 else
     return i + j*N;                                 //Alternative boundary condition, hard open boundary condition
}*/

void ermessage()
{
    cerr << "Please input (Sample size(L) | Cycles of evolution | Number of vortices)" << endl;
    exit(1);
}

void ermessage2()
{
    cerr << "invalid size of sample." << endl;
    exit(1);
}

void ermessage3()
{
	cerr << "The number of vortices can't be negative." << endl;
	exit(1);
}

void Decompose_wavefunc(field< complex<double> >& wavefunc, field<double>& ProbDen, field<double>& Phase)
{
    for (size_t i = 0; i < N2; ++i)
    {
        ProbDen[i] = abs(wavefunc[i]);
        Phase[i] = arg(wavefunc[i]);
    }
}

void print_data(field<vec2> &vecP, field<complex<double> >& wavefunc, field<vec2> &supC, field<double> &ProbDen, field<double> &Phase)
{
    vecP.report_data("Vector_potential");
    wavefunc.report_data("Wave_function");
    supC.report_data ("Supercurrent");
    
    Decompose_wavefunc(wavefunc, ProbDen, Phase);
    ProbDen.map_data("Probability");
    Phase.map_data("Phase");    //Data for density plot
}

void Import_Data(field<vec2> &vecP, field<complex<double> >& wavefunc, field<vec2> &supC, field<double> &ProbDen, field<double> &Phase)
{
    wavefunc.import_data("Initial_wave_function");
    vecP.import_data("Initial_Vector_potential");
    supC.import_data("Initial_supercurrent");
    ProbDen.import_map("Initial_Probability");
    Phase.import_map("Initial_Phase");
}

void initia_data(field<vec2> &vecP, field<complex<double> > &wavefunc, field<vec2> &supC, field<double> &ProbDen, field<double> &Phase, const size_t &vnum, const double &xsi)        //sample initializaiton: discretizing, fill initial values
{    
    const double dx = wavefunc.dx;
    const double dy = wavefunc.dy;
   
    for (size_t i = 0; i < N2; ++i)
    {
        wavefunc[i] = complex<double>(1.0,0.0);       //Make the amplitude of wavefunction throughout the system to unity for convenience of multiplication later
    
        //wavefunc[i] = (2*Rand()-1,2*Rand()-1)*(1./sqrt(2.));     //Alternative wave of set up the system with random fluctuation of phase and amplitude
    }
//const size_t vnum = rand()%10;            //Alternative way of generate random number of vortices to be placed in

vector<size_t> vindex(vnum);                   //create an arry of size in accordance with the number of vortices to store their positions
    
    for (size_t i = 0; i < vindex.size(); ++i)
    {
        //vindex[i] = N2/2+N/2- 29 + 48*i;                   //Put the vortex (vortices) in the center, multiple vortices overlap here physically means giant vortex state
        vindex[i] = Rand(N) + Rand(N)*N;				     //Alternative way of putting vortices in with random postion generated
		cout << "Vortex generated at " << "(" << xpos(vindex[i]) << ", " << ypos(vindex[i]) << ")" << "\t" << endl;              //check point
    }

	if (vnum <= 1)
		cout << vnum << " vortex generated." << endl;
	else
		cout << vnum << " vortices generated." << endl;      //check point
    
    for (size_t i = 0; i < vindex.size(); ++i)
    {
        vecP[vindex[i]] = vec2();
        wavefunc[vindex[i]] = complex<double>(0.0,0.0);      //Vanish wavefunction and vector potential at the enter of a vortex to be physical sensible
    }
    
    for (size_t i = 0; i < N2; ++i)
    {
        
        for (size_t j = 0; j < vindex.size(); ++j)
            if (i != vindex[j])
            {
                const double halfL = N*dx/2;
                const double L = 2*halfL;
                double x = (xpos(i)-xpos(vindex[j])) * dx;
                double y = (ypos(i)-ypos(vindex[j])) * dy;        //Compute relative distance between a vortex and the current cell
                
                if (x > halfL)
                    x -= L;
                if (x <= -halfL)
                    x += L;
                
                if (y > halfL)
                    y -= L;
                if (y <= -halfL)
                    y += L;         //Apply truncated rule to only take into account the influence from the most nearby vortex under periodic condition.
                
                double r2 = x*x + y*y;
                double r = sqrt(r2);
                double amp0 = r/274.0/(xsi*xsi + r2);      //determine the amplitude of local vector potential
                
                vec2 Alocal0(-y,x);                        //determine the direction of the local vector potential
                vecP[i] += amp0 * Alocal0;                 //take into account all the vortices within the cut-off range
				

            }
    }

	for (size_t i = 0; i < N2; ++i)
	{
		for (size_t j = 0; j < vindex.size(); ++j)
			if (i != vindex[j])
			{
				const double halfL = N*dx / 2;
				const double L = 2 * halfL;
				double x = (xpos(i) - xpos(vindex[j])) * dx;
				double y = (ypos(i) - ypos(vindex[j])) * dy;

				//double xx = (fabs(xpos(i)-xpos(vindex[j])) - N) * dx;
				//double yy = (fabs(ypos(i)-ypos(vindex[j])) - N) * dy;

				if (x > halfL)
					x -= L;
				if (x <= -halfL)
					x += L;
             
				if (y > halfL)
					y -= L;
				if (y <= -halfL)
					y += L;         //Apply truncated rule to only take into account the influence from the most nearby vortex under periodic condition.

				const double r = sqrt(x*x + y*y);
				double theta0 = atan2(y, x);            //determine the phase angle of the local wavefunction
				double phase = vnum*theta0;            //put in information about the number of total vortices
				complex<double> LocalPh(cos(phase), sin(phase));   //determine the phase of the local wavefunction
				complex<double> psiLocal = tanh(r / xsi)*LocalPh;       //Determine the amplitude of local wavefunction
				wavefunc[i] *= psiLocal;                        //take into account all vortices within the cut-off range
			}
    }
    
    for (size_t i = 0; i < N2; ++i) //calculating supercurrent
    {
        double jx = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i) + 1, ypos(i))] - wavefunc[index(xpos(i) - 1, ypos(i))]) - wavefunc[i] * (conj(wavefunc[index(xpos(i) + 1, ypos(i))]) - conj(wavefunc[index(xpos(i) - 1, ypos(i))]))) / dx - 274.0*conj(wavefunc[i])*wavefunc[i] * vecP[i].x);
        double jy = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i), ypos(i) + 1)] - wavefunc[index(xpos(i), ypos(i) - 1)]) - wavefunc[i] * (conj(wavefunc[index(xpos(i), ypos(i) + 1)]) - conj(wavefunc[index(xpos(i), ypos(i) - 1)]))) / dy - 274.0*conj(wavefunc[i])*wavefunc[i] * vecP[i].y);
        vec2 current(jx, jy);
        supC[i] = current;
    }
    
    vecP.report_data("Initial_Vector_potential");
    wavefunc.report_data("Initial_wave_function");
    supC.report_data("Initial_supercurrent");
    
    Decompose_wavefunc(wavefunc, ProbDen, Phase);
    ProbDen.map_data("Initial_Probability");
    Phase.map_data("Initial_Phase");
}

void linear_update_data(field<vec2>& vecP, field< complex<double> >& wavefunc, field<vec2> &supC, const size_t& rev, const double &xsi)       //updating function using linear GL equations
{
    const double dx = wavefunc.dx;
    const double dy = wavefunc.dy;
    const double alpha = -0.25/xsi/xsi;     //GL expansion coefficient (determins temperature in our case)
    for (size_t time = 0; time < rev; ++time)
    {
    for (size_t i = 0; i < N2; ++i)
    {
        if (ODD(xpos(i)+ypos(i)) == true)   //Chess board updating
        {
            const double ax = vecP[i].x;
            const double ay = vecP[i].y;
            
            double axu = vecP[index(xpos(i),ypos(i)+1)].x;   // up neighbour
            double axd = vecP[index(xpos(i),ypos(i)-1)].x;   // down neighbour
            double axl = vecP[index(xpos(i)-1,ypos(i))].x;   // left neighbour
            double axr = vecP[index(xpos(i)+1,ypos(i))].x;   // right neighbour
            
            double axur = vecP[index(xpos(i)+1,ypos(i)+1)].x; // up right neighbour
            double axul = vecP[index(xpos(i)-1,ypos(i)+1)].x; // up left neighbour
            double axdr = vecP[index(xpos(i)+1,ypos(i)-1)].x; // down right neighbour
            double axdl = vecP[index(xpos(i)-1,ypos(i)-1)].x; // down left neighbour

            double ayu = vecP[index(xpos(i),ypos(i)+1)].y;
            double ayd = vecP[index(xpos(i),ypos(i)-1)].y;
            double ayl = vecP[index(xpos(i)-1,ypos(i))].y;
            double ayr = vecP[index(xpos(i)+1,ypos(i))].y;
            
            double ayur = vecP[index(xpos(i)+1,ypos(i)+1)].y;
            double ayul = vecP[index(xpos(i)-1,ypos(i)+1)].y;
            double aydr = vecP[index(xpos(i)+1,ypos(i)-1)].y;
            double aydl = vecP[index(xpos(i)-1,ypos(i)-1)].y;


            
            complex<double> messR0 = xsi*xsi*((wavefunc[index(xpos(i)+1,ypos(i))] + wavefunc[index(xpos(i)-1,ypos(i))])/dx/dx + (wavefunc[index(xpos(i),ypos(i)+1)] + wavefunc[index(xpos(i),ypos(i)-1)])/dy/dy) + ci*274.0*(ax*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx + ay*(wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy);
            
            complex<double> messL0 = xsi*xsi*(2.0/dx/dx+2.0/dy/dy - 274.0*274.0 * dot(vecP[i], vecP[i]) - 274.0*ci*((axr-axl)/2.0/dx + (ayu-ayd)/2.0/dy)) + 1.0;
            
            wavefunc[i] = (1-damp)*wavefunc[i] + damp*messR0/messL0;      //upadate local wavefunction with damping factor
            
            complex<double> messRx = -0.25*(ayur-aydr-ayul+aydl)/dx/dy + (axu+axd)/dy/dy - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i])*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx - wavefunc[i]*(conj(wavefunc[index(xpos(i)+1,ypos(i))])-conj(wavefunc[index(xpos(i)-1,ypos(i))]))/2.0/dx );
            
            complex<double> messRy = -0.25*(axur-axdr-axul+axdl)/dx/dy + (ayr+ayl)/dx/dx - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i]) * (wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy - wavefunc[i]*(conj(wavefunc[index(xpos(i),ypos(i)+1)])-conj(wavefunc[index(xpos(i),ypos(i)-1)]))/2.0/dy);
            
            complex<double> messLx = 2.0/dy/dy + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            complex<double> messLy = 2.0/dx/dx + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            
            double axx = (1-damp)*ax + damp*real(messRx/messLx);
            double ayy = (1-damp)*ay + damp*real(messRy/messLy);
            
            vecP[i] = vec2(axx,ayy);       //update local vector potential with damping factor

        }
    }
    
    for (size_t i = 0; i < N2; ++i)
    {
        //const int l = xpos(i) + ypos(i);
        if (ODD(xpos(i)+ypos(i)) == false)   //Chess board update
        {
            const double ax = vecP[i].x;
            const double ay = vecP[i].y;
            
            double axu = vecP[index(xpos(i),ypos(i)+1)].x;   // up neighbour
            double axd = vecP[index(xpos(i),ypos(i)-1)].x;   // down neighbour
            double axl = vecP[index(xpos(i)-1,ypos(i))].x;   // left neighbour
            double axr = vecP[index(xpos(i)+1,ypos(i))].x;   // right neighbour
            
            double axur = vecP[index(xpos(i)+1,ypos(i)+1)].x; // up right neighbour
            double axul = vecP[index(xpos(i)-1,ypos(i)+1)].x; // up left neighbour
            double axdr = vecP[index(xpos(i)+1,ypos(i)-1)].x; // down right neighbour
            double axdl = vecP[index(xpos(i)-1,ypos(i)-1)].x; // down left neighbour
            
            double ayu = vecP[index(xpos(i),ypos(i)+1)].y;
            double ayd = vecP[index(xpos(i),ypos(i)-1)].y;
            double ayl = vecP[index(xpos(i)-1,ypos(i))].y;
            double ayr = vecP[index(xpos(i)+1,ypos(i))].y;
            
            double ayur = vecP[index(xpos(i)+1,ypos(i)+1)].y;
            double ayul = vecP[index(xpos(i)-1,ypos(i)+1)].y;
            double aydr = vecP[index(xpos(i)+1,ypos(i)-1)].y;
            double aydl = vecP[index(xpos(i)-1,ypos(i)-1)].y;

            
            complex<double> messR0 = xsi*xsi*((wavefunc[index(xpos(i)+1,ypos(i))] + wavefunc[index(xpos(i)-1,ypos(i))])/dx/dx + (wavefunc[index(xpos(i),ypos(i)+1)] + wavefunc[index(xpos(i),ypos(i)-1)])/dy/dy) + ci*274.0*(ax*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx + ay*(wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy);
            
            complex<double> messL0 = xsi*xsi*(2.0/dx/dx+2.0/dy/dy - 274.0*274.0 * dot(vecP[i], vecP[i]) - 274.0*ci*((axr-axl)/2.0/dx + (ayu-ayd)/2.0/dy)) + 1.0;
            
            wavefunc[i] = (1-damp)*wavefunc[i] + damp*messR0/messL0;
            
            complex<double> messRx = -0.25*(ayur-aydr-ayul+aydl)/dx/dy + (axu+axd)/dy/dy - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i])*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx - wavefunc[i]*(conj(wavefunc[index(xpos(i)+1,ypos(i))])-conj(wavefunc[index(xpos(i)-1,ypos(i))]))/2.0/dx );
            
            complex<double> messRy = -0.25*(axur-axdr-axul+axdl)/dx/dy + (ayr+ayl)/dx/dx - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i]) * (wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy - wavefunc[i]*(conj(wavefunc[index(xpos(i),ypos(i)+1)])-conj(wavefunc[index(xpos(i),ypos(i)-1)]))/2.0/dy);
            
            complex<double> messLx = 2.0/dy/dy + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            complex<double> messLy = 2.0/dx/dx + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            
            double axx = (1-damp)*ax + damp*real(messRx/messLx);
            double ayy = (1-damp)*ay + damp*real(messRy/messLy);
            
            vecP[i] = vec2(axx,ayy);

        }
    }
    }
    
    for (size_t i = 0; i < N2; ++i) //calculating supercurrent
    {
        double jx = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i) + 1, ypos(i))] - wavefunc[index(xpos(i) - 1, ypos(i))]) - wavefunc[i] * (conj(wavefunc[index(xpos(i) + 1, ypos(i))]) - conj(wavefunc[index(xpos(i) - 1, ypos(i))]))) / dx - 274.0*conj(wavefunc[i])*wavefunc[i] * vecP[i].x);
        double jy = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i), ypos(i) + 1)] - wavefunc[index(xpos(i), ypos(i) - 1)]) - wavefunc[i] * (conj(wavefunc[index(xpos(i), ypos(i) + 1)]) - conj(wavefunc[index(xpos(i), ypos(i) - 1)]))) / dy - 274.0*conj(wavefunc[i])*wavefunc[i] * vecP[i].y);
        vec2 current(jx, jy);
        supC[i] = current;
    }
}

void update_data(field<vec2>& vecP, field< complex<double> >& wavefunc, field<vec2> &supC, const size_t& rev, const double &xsi)
{
    const double dx = wavefunc.dx;
    const double dy = wavefunc.dy;
    const double alpha = -0.25/xsi/xsi;

    for (size_t time = 0; time < rev; ++time)
    {
    for (size_t i = 0; i < N2; ++i)
    {
        //const int l = xpos(i) + ypos(i);
        if (ODD(xpos(i)) == true && ODD(ypos(i)) == true)         //a variety of chess board updating
       {
            double ax = vecP[i].x;             //splite the local two dimensional vector potential into two directions
            double ay = vecP[i].y;
            
            double axu = vecP[index(xpos(i),ypos(i)+1)].x;
            double axd = vecP[index(xpos(i),ypos(i)-1)].x;
            double axl = vecP[index(xpos(i)-1,ypos(i))].x;
            double axr = vecP[index(xpos(i)+1,ypos(i))].x;
            
            double axur = vecP[index(xpos(i)+1,ypos(i)+1)].x;
            double axul = vecP[index(xpos(i)-1,ypos(i)+1)].x;
            double axdr = vecP[index(xpos(i)+1,ypos(i)-1)].x;
            double axdl = vecP[index(xpos(i)-1,ypos(i)-1)].x;
            
            
            double ayu = vecP[index(xpos(i),ypos(i)+1)].y;
            double ayd = vecP[index(xpos(i),ypos(i)-1)].y;
            double ayl = vecP[index(xpos(i)-1,ypos(i))].y;
            double ayr = vecP[index(xpos(i)+1,ypos(i))].y;
            
            double ayur = vecP[index(xpos(i)+1,ypos(i)+1)].y;
            double ayul = vecP[index(xpos(i)-1,ypos(i)+1)].y;
            double aydr = vecP[index(xpos(i)+1,ypos(i)-1)].y;
            double aydl = vecP[index(xpos(i)-1,ypos(i)-1)].y;

            
            complex<double> messR0 = xsi*xsi*((wavefunc[index(xpos(i)+1,ypos(i))] + wavefunc[index(xpos(i)-1,ypos(i))])/dx/dx + (wavefunc[index(xpos(i),ypos(i)+1)] + wavefunc[index(xpos(i),ypos(i)-1)])/dy/dy) + ci*274.0*(ax*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx + ay*(wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy);
            
            complex<double> messL0 = xsi*xsi*(2.0/dx/dx+2.0/dy/dy - 274.0*274.0 * dot(vecP[i],vecP[i]) - 274.0*ci*((axr-axl)/2.0/dx + (ayu-ayd)/2.0/dy)) + 1.0 - conj(wavefunc[i])*wavefunc[i];
            
            wavefunc[i] = (1-damp)*wavefunc[i] + damp*messR0/messL0;
            
            complex<double> messRx = -0.25*(ayur-aydr-ayul+aydl)/dx/dy + (axu+axd)/dy/dy - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i])*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx - wavefunc[i]*(conj(wavefunc[index(xpos(i)+1,ypos(i))])-conj(wavefunc[index(xpos(i)-1,ypos(i))]))/2.0/dx );
            
            complex<double> messRy = -0.25*(axur-axdr-axul+axdl)/dx/dy + (ayr+ayl)/dx/dx - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i]) * (wavefunc[index(xpos(i),ypos(i)+1)] - wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy - wavefunc[i]*(conj(wavefunc[index(xpos(i),ypos(i)+1)]) - conj(wavefunc[index(xpos(i),ypos(i)-1)]))/2.0/dy);
            
            complex<double> messLx = 2.0/dy/dy + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            complex<double> messLy = 2.0/dx/dx + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
           
            double axx = (1-damp)*ax + damp*real(messRx/messLx);
            double ayy = (1-damp)*ay + damp*real(messRy/messLy);
            
            vecP[i] = vec2(axx,ayy);       //computed respectively, and here recombine the two components of local vector potential.
            
        }
    }
    for (size_t i = 0; i < N2; ++i)
    {
        //const int l = xpos(i) + ypos(i);
        if (ODD(xpos(i)) == false && ODD(ypos(i)) == false)    //a variety of chess board update
        {
            const double ax = vecP[i].x;
            const double ay = vecP[i].y;
            
            double axu = vecP[index(xpos(i),ypos(i)+1)].x;
            double axd = vecP[index(xpos(i),ypos(i)-1)].x;
            double axl = vecP[index(xpos(i)-1,ypos(i))].x;
            double axr = vecP[index(xpos(i)+1,ypos(i))].x;
            
            double axur = vecP[index(xpos(i)+1,ypos(i)+1)].x;
            double axul = vecP[index(xpos(i)-1,ypos(i)+1)].x;
            double axdr = vecP[index(xpos(i)+1,ypos(i)-1)].x;
            double axdl = vecP[index(xpos(i)-1,ypos(i)-1)].x;
            
            
            double ayu = vecP[index(xpos(i),ypos(i)+1)].y;
            double ayd = vecP[index(xpos(i),ypos(i)-1)].y;
            double ayl = vecP[index(xpos(i)-1,ypos(i))].y;
            double ayr = vecP[index(xpos(i)+1,ypos(i))].y;
            
            double ayur = vecP[index(xpos(i)+1,ypos(i)+1)].y;
            double ayul = vecP[index(xpos(i)-1,ypos(i)+1)].y;
            double aydr = vecP[index(xpos(i)+1,ypos(i)-1)].y;
            double aydl = vecP[index(xpos(i)-1,ypos(i)-1)].y;

            
            complex<double> messR0 = xsi*xsi*((wavefunc[index(xpos(i)+1,ypos(i))] + wavefunc[index(xpos(i)-1,ypos(i))])/dx/dx + (wavefunc[index(xpos(i),ypos(i)+1)] + wavefunc[index(xpos(i),ypos(i)-1)])/dy/dy) + ci*274.0*(ax*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx + ay*(wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy);
            
            complex<double> messL0 = xsi*xsi*(2.0/dx/dx+2.0/dy/dy - 274.0*274.0 * dot(vecP[i], vecP[i]) - 274.0*ci*((axr-axl)/2.0/dx + (ayu-ayd)/2.0/dy)) + 1.0 - conj(wavefunc[i])*wavefunc[i];
           
            wavefunc[i] = (1-damp)*wavefunc[i] + damp*messR0/messL0;
           
            complex<double> messRx = -0.25*(ayur-aydr-ayul+aydl)/dx/dy + (axu+axd)/dy/dy - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i])*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx - wavefunc[i]*(conj(wavefunc[index(xpos(i)+1,ypos(i))])-conj(wavefunc[index(xpos(i)-1,ypos(i))]))/2.0/dx );
            
            complex<double> messRy = -0.25*(axur-axdr-axul+axdl)/dx/dy + (ayr+ayl)/dx/dx - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i]) * (wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy - wavefunc[i]*(conj(wavefunc[index(xpos(i),ypos(i)+1)])-conj(wavefunc[index(xpos(i),ypos(i)-1)]))/2.0/dy);
            
            complex<double> messLx = 2.0/dy/dy + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            complex<double> messLy = 2.0/dx/dx + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            
            double axx = (1-damp)*ax + damp*real(messRx/messLx);
            double ayy = (1-damp)*ay + damp*real(messRy/messLy);
            
            vecP[i] = vec2(axx,ayy);

        }
    }
    for (size_t i = 0; i < N2; ++i)
    {
        //const int l = xpos(i) + ypos(i);
        if (ODD(xpos(i)) == false && ODD(ypos(i)) == true)       //a variety of chess board updating
        {
            double ax = vecP[i].x;
            double ay = vecP[i].y;
            
            double axu = vecP[index(xpos(i),ypos(i)+1)].x;
            double axd = vecP[index(xpos(i),ypos(i)-1)].x;
            double axl = vecP[index(xpos(i)-1,ypos(i))].x;
            double axr = vecP[index(xpos(i)+1,ypos(i))].x;
            
            double axur = vecP[index(xpos(i)+1,ypos(i)+1)].x;
            double axul = vecP[index(xpos(i)-1,ypos(i)+1)].x;
            double axdr = vecP[index(xpos(i)+1,ypos(i)-1)].x;
            double axdl = vecP[index(xpos(i)-1,ypos(i)-1)].x;
            
            
            double ayu = vecP[index(xpos(i),ypos(i)+1)].y;
            double ayd = vecP[index(xpos(i),ypos(i)-1)].y;
            double ayl = vecP[index(xpos(i)-1,ypos(i))].y;
            double ayr = vecP[index(xpos(i)+1,ypos(i))].y;
            
            double ayur = vecP[index(xpos(i)+1,ypos(i)+1)].y;
            double ayul = vecP[index(xpos(i)-1,ypos(i)+1)].y;
            double aydr = vecP[index(xpos(i)+1,ypos(i)-1)].y;
            double aydl = vecP[index(xpos(i)-1,ypos(i)-1)].y;
            
            
            complex<double> messR0 = xsi*xsi*((wavefunc[index(xpos(i)+1,ypos(i))] + wavefunc[index(xpos(i)-1,ypos(i))])/dx/dx + (wavefunc[index(xpos(i),ypos(i)+1)] + wavefunc[index(xpos(i),ypos(i)-1)])/dy/dy) + ci*274.0*(ax*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx + ay*(wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy);
            
            complex<double> messL0 = xsi*xsi*(2.0/dx/dx+2.0/dy/dy - 274.0*274.0 * dot(vecP[i],vecP[i]) - 274.0*ci*((axr-axl)/2.0/dx + (ayu-ayd)/2.0/dy)) + 1.0 - conj(wavefunc[i])*wavefunc[i];
            
            wavefunc[i] = (1-damp)*wavefunc[i] + damp*messR0/messL0;
            
            complex<double> messRx = -0.25*(ayur-aydr-ayul+aydl)/dx/dy + (axu+axd)/dy/dy - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i])*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx - wavefunc[i]*(conj(wavefunc[index(xpos(i)+1,ypos(i))])-conj(wavefunc[index(xpos(i)-1,ypos(i))]))/2.0/dx );
            
            complex<double> messRy = -0.25*(axur-axdr-axul+axdl)/dx/dy + (ayr+ayl)/dx/dx - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i]) * (wavefunc[index(xpos(i),ypos(i)+1)] - wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy - wavefunc[i]*(conj(wavefunc[index(xpos(i),ypos(i)+1)]) - conj(wavefunc[index(xpos(i),ypos(i)-1)]))/2.0/dy);
            
            complex<double> messLx = 2.0/dy/dy + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            complex<double> messLy = 2.0/dx/dx + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            
            double axx = (1-damp)*ax + damp*real(messRx/messLx);
            double ayy = (1-damp)*ay + damp*real(messRy/messLy);
            
            vecP[i] = vec2(axx,ayy);
            
        }
    }
    for (size_t i = 0; i < N2; ++i)
    {
        //const int l = xpos(i) + ypos(i);
        if (ODD(xpos(i)) == true && ODD(ypos(i)) == false)     //a variety of chess board updating
        {
            double ax = vecP[i].x;
            double ay = vecP[i].y;
            
            double axu = vecP[index(xpos(i),ypos(i)+1)].x;
            double axd = vecP[index(xpos(i),ypos(i)-1)].x;
            double axl = vecP[index(xpos(i)-1,ypos(i))].x;
            double axr = vecP[index(xpos(i)+1,ypos(i))].x;
            
            double axur = vecP[index(xpos(i)+1,ypos(i)+1)].x;
            double axul = vecP[index(xpos(i)-1,ypos(i)+1)].x;
            double axdr = vecP[index(xpos(i)+1,ypos(i)-1)].x;
            double axdl = vecP[index(xpos(i)-1,ypos(i)-1)].x;
            
            
            double ayu = vecP[index(xpos(i),ypos(i)+1)].y;
            double ayd = vecP[index(xpos(i),ypos(i)-1)].y;
            double ayl = vecP[index(xpos(i)-1,ypos(i))].y;
            double ayr = vecP[index(xpos(i)+1,ypos(i))].y;
            
            double ayur = vecP[index(xpos(i)+1,ypos(i)+1)].y;
            double ayul = vecP[index(xpos(i)-1,ypos(i)+1)].y;
            double aydr = vecP[index(xpos(i)+1,ypos(i)-1)].y;
            double aydl = vecP[index(xpos(i)-1,ypos(i)-1)].y;
            
            
            complex<double> messR0 = xsi*xsi*((wavefunc[index(xpos(i)+1,ypos(i))] + wavefunc[index(xpos(i)-1,ypos(i))])/dx/dx + (wavefunc[index(xpos(i),ypos(i)+1)] + wavefunc[index(xpos(i),ypos(i)-1)])/dy/dy) + ci*274.0*(ax*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx + ay*(wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy);
            
            complex<double> messL0 = xsi*xsi*(2.0/dx/dx+2.0/dy/dy - 274.0*274.0 * dot(vecP[i],vecP[i]) - 274.0*ci*((axr-axl)/2.0/dx + (ayu-ayd)/2.0/dy)) + 1.0 - conj(wavefunc[i])*wavefunc[i];
            
            wavefunc[i] = (1-damp)*wavefunc[i] + damp*messR0/messL0;
            
            complex<double> messRx = -0.25*(ayur-aydr-ayul+aydl)/dx/dy + (axu+axd)/dy/dy - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i])*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))])/2.0/dx - wavefunc[i]*(conj(wavefunc[index(xpos(i)+1,ypos(i))])-conj(wavefunc[index(xpos(i)-1,ypos(i))]))/2.0/dx );
            
            complex<double> messRy = -0.25*(axur-axdr-axul+axdl)/dx/dy + (ayr+ayl)/dx/dx - ci * 274.0 * M_PI * fabs(alpha) * (conj(wavefunc[i]) * (wavefunc[index(xpos(i),ypos(i)+1)] - wavefunc[index(xpos(i),ypos(i)-1)])/2.0/dy - wavefunc[i]*(conj(wavefunc[index(xpos(i),ypos(i)+1)]) - conj(wavefunc[index(xpos(i),ypos(i)-1)]))/2.0/dy);
            
            complex<double> messLx = 2.0/dy/dy + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            complex<double> messLy = 2.0/dx/dx + 8.0*137.0*137.0*M_PI*fabs(alpha)*conj(wavefunc[i])*wavefunc[i];
            
            double axx = (1-damp)*ax + damp*real(messRx/messLx);
            double ayy = (1-damp)*ay + damp*real(messRy/messLy);
            
            vecP[i] = vec2(axx,ayy);
            
        }
    }
        if (time%200 == 0)
            cout << rev-time << " rounds left" << endl;

    }
    
    for (size_t i = 0; i < N2; ++i) //calculating supercurrent
    {
        double jx = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i) + 1, ypos(i))] - wavefunc[index(xpos(i) - 1, ypos(i))]) - wavefunc[i] * (conj(wavefunc[index(xpos(i) + 1, ypos(i))]) - conj(wavefunc[index(xpos(i) - 1, ypos(i))]))) / dx - 274.0*conj(wavefunc[i])*wavefunc[i] * vecP[i].x);
        double jy = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i), ypos(i) + 1)] - wavefunc[index(xpos(i), ypos(i) - 1)]) - wavefunc[i] * (conj(wavefunc[index(xpos(i), ypos(i) + 1)]) - conj(wavefunc[index(xpos(i), ypos(i) - 1)]))) / dy - 274.0*conj(wavefunc[i])*wavefunc[i] * vecP[i].y);
        vec2 current(jx, jy);
        supC[i] = current;
    }
}


int main(int argc, char *argv[])
{
    if (argc != 4) ermessage();
    
    const double L = atof(argv[1]);     //take the first argument input as the size of sample.
    const double rev = atof(argv[2]);   //take the second argument input as the count of relaxation time.
	const size_t vnum = atof(argv[3]);  //take the third argument input as vortex number.
    
    if (L <= 0) ermessage2();
	if (vnum <= 0) ermessage3();
 
    field<vec2> vecP(L,N);  //create a vector field throughout the sample (see "field.h")
    vecP[N2] = vec2(0.0, 0.0);  //make the additional one cell zero to behave as a dummy margin when usin open boundary condition
    field<complex<double> > wavefunc(L,N);  //create a complex scaler field throughout the sample (see "field.h")
    wavefunc[N2] = complex<double>(0.0,0.0);    //make the additional one cell zero to behave as a dummy margin when usin open boundary condition
    field<double> ProbDen(L,N);
    ProbDen[N2] = 0;
    
    field<double> Phase(L,N);
    Phase[N2] = 0;
    
	field<vec2> supC(L, N); //create a supercurrent field
	supC[N2] = vec2(0.0, 0.0);

	const double xsi = 200.0 * L / N;    //Determine coherent length according to discretization

	cout << "Sample size: " << L << " x " << L << endl;
	cout << "Discrete cell size: " << L / N << endl;
	cout << "coherent length: " << xsi * N / L << endl;

    initia_data(vecP, wavefunc, supC, ProbDen, Phase, vnum, xsi);           //initialize the fields (a report of initial conditions are included in this function)
    
    //Import_Data(vecP, wavefunc, supC, ProbDen, Phase);                //Alternative way of initialization by reading a complete set of data from files.
    /*
    field<vec2> LvecP = vecP;
    field<complex<double> > Lwavefunc = wavefunc;
    field<vec2> LsupC = supC;   //copy the initial condition for linear GL equations
    
    linear_update_data(LvecP, Lwavefunc, LsupC, rev, xsi);       //update the copy using linear GL equations
     */
    update_data(vecP, wavefunc, supC, rev, xsi);                //update the original copy using non-linear equations
    print_data(vecP, wavefunc, supC, ProbDen, Phase);
    return 0;
}
