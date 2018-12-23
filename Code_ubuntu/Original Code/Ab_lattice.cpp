#include <cstddef>
using std::size_t;

#include <cassert>

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

#include <complex>
using std::complex;
using std::real;
using std::abs;
using std::arg;

#include "field.hpp"

#include "mtrand.hpp"
mtrand Rand(2000);

double damp = 0.3;                //choose damping factor for cell update
int N = 299;                   //choose discretizing step (see "field.h").
int N2 = N*N;
complex<double> ci(0.0,1.0);      //create mathematical imaginary number "i"

template<class T>
inline int xpos(T i) {return i%N;}        //Extract position on x-axis from index

template<class T>
inline int ypos(T i) {return i/N;}        //Extract position on y-axis from index

template<class T>
inline bool ODD(T x){return (x%2 == 0? false:true);}       //Determine whether or not a number is odd

template<class T>
inline size_t index(T i, T j)//{return ((i+N)%N)+((j+N)%N)*N;}       //Periodic boundary condition for cell update
{
   if(i >= N or i < 0 or j >= N or j < 0)
    return (N2);
 else
     return i + j*N;                                 //Alternative boundary condition, hard open boundary condition
}

template<class T>                                    //Determine the minimum number between two values of the same type
inline T min(T x, T y){return (x>y?y:x);}



void ermessage()
{
    cerr << "(size(L) | cycles of revolution)" << endl;
    exit(1);
}

void ermessage2()
{
    cerr << "invalid size of sample" << endl;
    exit(1);
}


void read_grid(field<vec2> &vecP, field<complex<double> >& wavefunc)
{
    ifstream fin("filename(vector potential).dat");    //read files containing initialization information of vector potential
    if (!fin.is_open())
    {
        cerr << "ERROR: Can't open the file" << endl;
        exit(1);
    }
    size_t i = 0;
    while (true)
    {
        fin >> vecP[i].x >> vecP[i].y >> vecP[i];
        ++i;
        if (fin.eof())
        {
            break;
        }
    }
    fin.close();
    
    ifstream fin2("filename(wavefunction).dat");    //read files containing initialization information of wavefunction (stored differently from that of vector potential)
    if (!fin2.is_open())
    {
        cerr << "ERROR: Can't open the file" << endl;
        exit(1);
    }
    
    for (size_t i = 0; i < N2; ++i)
    {
        fin2 >> wavefunc[index(xpos(i),ypos(i))];
        if (fin2.eof())
        {
            break;
        }
    }
    fin2.close();
}

void initiallize_grid(field<vec2> &vecP, field<complex<double> > &wavefunc)        //sample initializaiton: discretizing, fill initial values
{    
    const double dx = wavefunc.dx;
    const double dy = wavefunc.dy;
    const double xsi = 100.0 * dx;                    //Determine coherent length according to discretizing
    
    for (int i = 0; i < N2; ++i)
        wavefunc[i] = complex<double>(1.0,0.0);       //Make the amplitude of wavefunction throughout the system to unity for convenience of multiplication later
        //wavefunc[i] = (2*Rand()-1,2*Rand()-1)*(1./sqrt(2.));     //Alternative wave of set up the system with random fluctuation of phase and amplitude
    
    const int vnum = 3;                     //Choose number of vortices to be placed in manually
    //const size_t vnum = rand()%10;            //Alternative way of generate random number of vortices to be placed in
    cout << vnum << endl;                                      //check point
    vector<int> vindex(vnum);                   //create an arry of size in accordance with the number of vortices to store their positions
    
    for (size_t i = 0; i < vindex.size(); ++i)
    {
        //vindex[i] = N2/2+N/2- 29 + 48*i;                   //Put the vortex (vortices) in the center, multiple vortices overlap here physically means giant vortex state
        vindex[i] = Rand(N) + Rand(N)*N;      //Alternative way of putting vortices in with random postion generated
        cout << i << "\t" << vindex[i] << endl;                //check point
    }
    
    for (size_t i = 0; i < vindex.size(); ++i)
    {
        vecP[vindex[i]] = vec2();
        wavefunc[vindex[i]] = complex<double>(0.0,0.0);      //Vanish wavefunction and vector potential at the enter of a vortex to be physical sensible
    }
    
    for (int i = 0; i < N2; ++i)
    {
        
        for (size_t j = 0; j < vindex.size(); ++j)
            if (i != vindex[j])
            {
                const double halfL = N*dx/2;
                const double L = 2*halfL;
                double x = (xpos(i)-xpos(vindex[j])) * dx;
                double y = (ypos(i)-ypos(vindex[j])) * dy;        //Compute relative distance        between a vortex and the current cell
                
                if (x > halfL)
                    x -= L;
                if (x <= -halfL)
                    x += L;
                
                if (y > halfL)
                    y -= L;
                if (y <= -halfL)                            //Apply truncated rule to only take into account the influence from th most nearby vortex
                    y += L;
                
                double r2 = x*x + y*y;
                double r = sqrt(r2);
                double amp0 = r/274.0/(xsi*xsi + r2);      //determine the amplitude of local vector potential
                
                vec2 Alocal0(-y,x);                        //determine the direction of the local vector potential
                vecP[i] += amp0 * Alocal0;                 //take into account all the vortices within the cut-off range
            }
    }
    for (int i = 0; i < N2; ++i)
    {
        for (size_t j = 0; j < vindex.size(); ++j)
            if (i != vindex[j])
            {
                const double halfL = N*dx/2;
                const double L = 2*halfL;
                double x = (xpos(i)-xpos(vindex[j])) * dx;
                double y = (ypos(i)-ypos(vindex[j])) * dy;
                //double xx = (fabs(xpos(i)-xpos(vindex[j])) - N) * dx;
                //double yy = (fabs(ypos(i)-ypos(vindex[j])) - N) * dy;
                
               if (x > halfL)
                    x -= L;
                if (x <= -halfL)
                    x += L;
                
                if (y > halfL)
                    y -= L;
                if (y <= -halfL)
                    y += L;
               
                //double x = min(x0,xx);
                //double y = min(y0,yy);

                const double r = sqrt(x*x + y*y);
                double theta0 = atan2(y,x);            //determine the phase angle of the local wavefunction
                 
                double phase = vnum*theta0;            //put in information about the number of total vortices
                complex<double> LocalPh (cos(phase),sin(phase));   //determine the phase of the local wavefunction
                complex<double> psiLocal = tanh(r/xsi)*LocalPh;       //Determine the amplitude of local wavefunction
                wavefunc[i] *= psiLocal;                        //take into account all vortices within the cut-off range
            }
    }

    
    ofstream fout;                    //generate a report of initial conditions
    fout.precision(6);
    fout.open("Vfield_init_cond.dat");    //vector potential
    for (int i = 0; i < N2; ++i)
    {
        fout << xpos(i) << "\t" << ypos(i) << "\t" << vecP[i] << endl;
    }
    
    fout.close();
    
    fout.open("wavefunc_init_cond.dat");    //amplitude part of the wavefunction (order parameter)
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            fout << abs(wavefunc[i + j*N]) << "\t";
        }
        fout << endl;
    }
    fout.close();
    
    fout.open("phase_init.dat");            // Phase part of the wavefunction (order parameter)
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            fout << arg(wavefunc[i + j*N]) << "\t";
        }
        fout << endl;
    }
    fout.close();

    
    fout.open("supercurrent_init.dat");      //supercurrent generated using info about local wavefunctions and vector potentials
    for (int i = 0; i < N2; ++i)
    {
        double jx = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))]) - wavefunc[i]*(conj(wavefunc[index(xpos(i)+1,ypos(i))])-conj(wavefunc[index(xpos(i)-1,ypos(i))])))/dx - 274.0*conj(wavefunc[i])*wavefunc[i]*vecP[i].x );

        double jy = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)]) - wavefunc[i]*(conj(wavefunc[index(xpos(i),ypos(i)+1)])-conj(wavefunc[index(xpos(i),ypos(i)-1)])))/dy - 274.0*conj(wavefunc[i])*wavefunc[i]*vecP[i].y );
        
        vec2 current(jx,jy);           //Determing the local supercurrent density according to quantum mechanics convention, using data from non-linear GL equations
        
        fout << xpos(i) << "\t" << ypos(i) << "\t" << current << endl;
    }
    fout.close();

}

void linear_update_grid(field<vec2>& vecP, field< complex<double> >& wavefunc, const size_t& rev)       //updating function using linear GL equations
{
    const double dx = wavefunc.dx;
    const double dy = wavefunc.dy;
    const double xsi = 100.0 * dx;          //coherent length
    const double alpha = -0.25/xsi/xsi;     //GL expansion coefficient (determins temperature in our case)
    for (size_t time = 0; time < rev; ++time)
    {
    for (int i = 0; i < N2; ++i)
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
    
    for (int i = 0; i < N2; ++i)
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
}

void update_grid(field<vec2>& vecP, field< complex<double> >& wavefunc, const size_t& rev)
{
    const double dx = wavefunc.dx;
    const double dy = wavefunc.dy;
    const double xsi = 100.0 * dx;
    const double alpha = -0.25/xsi/xsi;

    for (size_t time = 0; time < rev; ++time)
    {
    for (int i = 0; i < N2; ++i)
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
            
            vecP[i] = vec2(axx,ayy);       //compute respectively, here recombine the two components of local vector potential together to make one.
            
        }
    }
    for (int i = 0; i < N2; ++i)
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
    for (int i = 0; i < N2; ++i)
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
    for (int i = 0; i < N2; ++i)
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
}

void report_data(field<vec2>& vecP, field< complex<double> >& wavefunc, field<vec2>& LvecP, field< complex<double> >& Lwavefunc)
{
    const double dx = wavefunc.dx;
    const double dy = wavefunc.dy;
    
    ofstream fout;
    ofstream fout2;
    fout.precision(6);
    fout2.precision(6);
    fout.open("Vfield_results.dat");
    fout2.open("Vfield_Lresults.dat");
    for (int i = 0; i < N2; ++i)
    {
        fout << xpos(i) << "\t" << ypos(i) << "\t" << vecP[i] << endl;
        fout2 << xpos(i) << "\t" << ypos(i) << "\t" << LvecP[i] << endl;
    }
    
    fout.close();
    fout2.close();
    
    fout.open("wavefunc_results.dat");
    fout2.open("wavefunc_Lresults.dat");
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            fout << abs(wavefunc[i + j*N]) << "\t";
            fout2 << abs(Lwavefunc[i + j*N]) << "\t";
        }
        fout << endl;
        fout2 << endl;
    }
    fout.close();
    fout2.close();
    
    fout.open("supercurrent_results.dat");
    for (int i = 0; i < N2; ++i)
    {
        double jx = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i)+1,ypos(i))]-wavefunc[index(xpos(i)-1,ypos(i))]) - wavefunc[i]*(conj(wavefunc[index(xpos(i)+1,ypos(i))])-conj(wavefunc[index(xpos(i)-1,ypos(i))])))/dx - 274.0*conj(wavefunc[i])*wavefunc[i]*vecP[i].x );
        
        double jy = real(-0.25 * ci * (conj(wavefunc[i])*(wavefunc[index(xpos(i),ypos(i)+1)]-wavefunc[index(xpos(i),ypos(i)-1)]) - wavefunc[i]*(conj(wavefunc[index(xpos(i),ypos(i)+1)])-conj(wavefunc[index(xpos(i),ypos(i)-1)])))/dy - 274.0*conj(wavefunc[i])*wavefunc[i]*vecP[i].y );
       
        vec2 current(jx,jy);
        
        fout << xpos(i) << "\t" << ypos(i) << "\t" << current << endl;
    }
    fout.close();
    
    fout.open("phase_results.dat");
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            fout << arg(wavefunc[i + j*N]) << "\t";
        }
        fout << endl;
    }
    fout.close();


}

int main(int argc, char *argv[])
{
    if (argc != 3) ermessage();
    
    const double L = atof(argv[1]);     //take the first argument input as the size of sample
    const double rev = atof(argv[2]);   //take the second argument input as the the count of relaxation time.
    
    if (L <= 0) ermessage2();
    
    field<vec2> vecP(L,N);              //create a vector field throughout the sample (see "field.h")
    vecP[N2] = vec2(0.0,0.0);              //make the additional one cell zero to behave as a dummy margin when usin open boundary condition
    field<complex<double> > wavefunc(L,N);          //create a complex scaler field throughout the sample (see "field.h")
    wavefunc[N2] = complex<double>(0.0,0.0);            //make the additional one cell zero to behave as a dummy margin when usin open boundary condition
    
    field<vec2> LvecP(L,N);                     //create a copy of the vector field to be updated using linear GL equations
    field<complex<double> > Lwavefunc(L,N);     //create a copy of the scalar field to be updated using libear GL equations

    initiallize_grid(vecP, wavefunc);           //initialize the fields (a report of initial conditions are included in this function)
    //read_grid(vecP, wavefunc);                //Alternative way of initialization by reading a complete set of data from a file.

    LvecP = vecP;                         //copy the initial condition
    Lwavefunc = wavefunc;                 //copy the initial condition
    

    linear_update_grid(LvecP, Lwavefunc, rev);       //update the copy using linear GL equations
    update_grid(vecP, wavefunc, rev);                //update the original copy using non-linear equations
      
    report_data(vecP, wavefunc, LvecP, Lwavefunc);  //report the results from both calculations of both fields and their copies
    return 0;
}
