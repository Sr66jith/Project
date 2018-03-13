/* 1d slab with phase change */

/* Properties of ICE (Engineering Tool Box)     Properties from Zivkovic's Paper
    k   = 2.22  W/m.K                           den.Cp  =   2.5x10^6    J/m^3.K
    den = 916.2 kg/m^3                          den.L   =   10^8        J/m^3
    Tm  = 0     C                               delx    =   0.125       m
    Cp  = 2.050 kJ/kg.K                         k       =   2           W/m.K
                                                Timestep=   50 
                                                delt    =   4.32 * 10^4 s 
                                                h       =   16          W/m^2.K */

# include <iostream>
# include <math.h>

# define t      2.5                     // in metres
# define delt   43200.0                 // in seconds
# define delx   0.125                   // in metres
# define h      16.0                    // in W/m^2.K
# define k      2.0                     // in W/m.K
# define tm     0.0                     // degree celsius - melting temperature
# define ti     -2.0                    // degree celsius - initial temperature
# define tf     10.0                    // degree celsius - fluid temperature
# define tstep  5.0

using namespace std;

int node = t/delx + 1;                  // 21
int size = node + 1;                    // 22

/* Dimensionless constants */

float Tm    = (tm - tm) / (tf - tm);    // 0
float Tf    = (tf - tm) / (tf - tm);    // 1
float Ti    = (ti - tm) / (tf - tm);    // -0.2
float alpha = 2.5 * pow(10,6) / 2;      // 1.25 x 10^6
float delT  = alpha * delt / (t * t);   // 8.64 * 10^9
float delX  = delx / t;                 // 0.05
float F     = delT / (delX * delX);     // 2.7648 x 10^6
float Bi    = h * t / k;                // 20

void print(float);
void output(float*);
void initialize(float*, float*, float*, float*, float*, float*);
void tdma(float*, float*, float*, float*, float*);
void transient(float*, float*, float*, float*, float*, float*);
void update(float*, float*,float*);

int main()
{
    float ld[size], md[size], ud[size], rs[node];
    float T[size], To[size];

    initialize(T, To, ld, md, ud, rs);
    transient(T, To, ld, md, ud, rs);
}

void initialize(float T[], float To[], float ld[], float md[], float ud[], float rs[])
{
    float ah, ai, aj;

    ah = -F;
    ai = 1 + F;
    aj = -F;

    for(int i = 1; i <= node; i++) {
        ld[i] = ah;
        md[i] = ai;
        ud[i] = aj;
        To[i] = Ti;
        rs[i] = To[i];
    }
    
    ld[1] = 0;
    ud[1] = ah + aj;
    ld[node] = ah + aj;
    md[node] = ai - (ah * Bi * delX);
    ud[node] = 0;
    rs[node] = Ti - (aj * Bi * delX);
}

void transient(float T[], float To[], float ld[], float md[], float ud[], float rs[])
{
    for(int i = 1; i <= tstep; i++) {
        tdma(T, rs, ld, md, ud);
        update(T, To, rs);
        output(T);
    }
}

void tdma(float T[], float rs[], float ld[], float md[], float ud[])
{
    for(int i = 2; i <= node; i++) {
        ld[i] = ld[i] / md[i-1];
        md[i] = md[i] - (ld[i] * ud[i-1]);
        rs[i] = rs[i] - (ld[i] * rs[i-1]);
    }

    T[node] = rs[node] / md[node];

    for(int i = node -1; i > 0; i--) 
        T[i]    = (rs[i] - (ud[i] * T[i+1])) / md[i];
}
    
void update(float T[], float To[], float rs[])
{
    for(int i = 1; i < size; i++) {
        To[i] = T[i];
        rs[i] = T[i];
    }
}

void print(float a)
{
    cout << a << "\n";
}

void output(float a[])
{
    for(int i = 1; i <= node; i++)
        cout << a[i] << "\n";
    cout << "\n";
}
