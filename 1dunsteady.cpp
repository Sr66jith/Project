/* 1 D unsteady heat conduction */
# include <iostream>
# define node 5

using namespace std;

// Geometry variable
int g = 0;

float Bi = 2.0, delx = 1.0/node;

int main()
{
    // Function definitions
    void input(float *, float *, float *, float *);
    void tdma(float *, float *, float *, float *, float *);
    void output(float *);

    // Variable declarations
    float ld[node], md[node], ud[node], rhs[node], temp[node];

    // Function calls
    input(ld, md, ud, rhs);
    tdma(ld, md, ud, rhs, temp);
    output(temp);
    return 0;
}

void input(float ld[], float md[], float ud[], float rhs[])
{
    int i;

    for(i = 1; i < node; i++) {
        ld[i] = 1.0;
        md[i] = -2.0;
        ud[i] = 1.0;
        rhs[i] = 0.0;
    }

    rhs[node-1] = -1;
}

void tdma(float ld[], float md[], float rhs[], float x[])
{
    int k;

    for(k = 1; k < node; k++) {
        ld[k] = ld[k] / md[k-1];
        md[k] = md[k] - (ld[k] * ud[k-1]);
    }

