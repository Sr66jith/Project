/* 1 D unsteady heat conduction */
# include <iostream>
# define node 5
# define time 10
using namespace std;

// Geometry variable
int g = 0;

float Bi = 2.0, delx = 1.0/(node-1), delt = 0.5;

void output(float *);

int main()
{
    // Variable declarations
    float ld[node], md[node], ud[node], rhs[node], temp[node];

    // Function declarations
    void input(float *, float *, float *, float *);
    void transient(float *, float *, float *, float *, float *);

    // Function calls
    input(ld, md, ud, rhs);
    transient(ld, md, ud, rhs, temp);
    output(temp);

    return 0;
}

void input(float ld[], float md[], float ud[], float rhs[])
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
    }

    ud[0] = 2 * c;
    md[node-1] = a + (2 * Bi * delx);
    ld[node-1] = 2 * c;
}

void transient(float ld[], float md[], float ud[], float rhs[], float temp[])
{
    void tdma(float *, float *, float *, float *, float *);
    void update(float *, float *);

    for(float i = 0; i <= time; i += delt) {
        tdma(ld, md, ud, rhs, temp);
        update(temp, rhs);
    }
}

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

void update(float x[], float rhs[])
{
    for(int i = 0; i < node; i++)
        rhs[i] = x[i];
}

void output(float a[])
{
    for(int i = 0; i < node; i++)
        cout << a[i] << " ";
}
