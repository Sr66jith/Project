/* 1 D conduction dominated phase change */
# include <iostream>
# define node 5
# define time 1
using namespace std;

// Geometry variable
int g = 0;

float Ste = 0.05, Bi = 2.0, theta_m = 0.2; // Assuming water as the heat transfer fluid
float delx = 1.0/(node-1), delt = 0.5; 

void output(float *);

int main()
{
    // Variable declarations
    float ld[node], md[node], ud[node], rhs[node], temp[node], liq_frac[node];

    // Function declarations
    void input(float *, float *, float *, float *, float *);
    void transient(float *, float *, float *, float *, float *, float *);

    // Function calls
    input(ld, md, ud, rhs, liq_frac);
    transient(ld, md, ud, rhs, temp, liq_frac);
    output(temp);

    return 0;
}

// To assign the initial values for the variables
void input(float ld[], float md[], float ud[], float rhs[], float liq_frac[])
{
    int i;
    float c, a;

    c = -delt / (delx * delx);
    a = 1 - (2 * c);

    for(i = 0; i < node; i++) {
        ld[i] = c;
        md[i] = a;
        ud[i] = c;
        rhs[i] = 1.0;
        liq_frac[i] = 0.0;
    }

    ud[0] = 2 * c;
    md[node-1] = a + (2 * Bi * delx);
    ld[node-1] = 2 * c;
}

// To execute the transient function for determining the temp. and liquid
// fraction
void transient(float ld[], float md[], float ud[], float rhs[], float temp[], float liq_frac[])
{
    void tdma(float *, float *, float *, float *, float *);
    void update(float *, float *);
    void phase(float *, float *);

    for(float i = 0; i <= time; i += delt) {
        tdma(ld, md, ud, rhs, temp);
        update(temp, rhs);
        phase(temp, liq_frac);
    }
}

// To execute the tdma algorithm
void tdma(float ld[], float md[], float ud[], float rhs[], float x[])
{
    int k;

    for(k = 1; k < node; k++) {
        ld[k]  = ld[k] / md[k-1];
        md[k]  = md[k] - (ld[k] * ud[k-1]);
        rhs[k] = rhs[k] - (ld[k] * rhs[k-1]);
    }

    x[node-1] = rhs[node-1] / md[node-1];

    for(k = node-2; k >=0; k--)
        x[k] = (rhs[k] - (ud[k] * x[k+1])) / md[k];
}

// To update the previous temperatures with the new one
void update(float x[], float rhs[])
{
    for(int i = 0; i < node; i++)
        rhs[i] = x[i];
}

// To update the liquid fraction term
void phase(float temp[], float liq_frac[])
{
    int i;
    float c = (Ste * delt) / (delx * delx);

    for(i = 1; i < (node-1); i++) {
        if(temp[i] >= theta_m) { // Needs checking
            liq_frac[i] = liq_frac[i] + 2 * c * theta_m - c * temp[i+1] - c * temp[i-1];
        }


// To output any array of size node
void output(float a[])
{
    for(int i = 0; i < node; i++)
        cout << a[i] << " ";
}
