/* For 1 dimensional steady state conduction */

# include <iostream>
# define node 5
using namespace std;

// Geometry variable (0 for slab and 1 for cylinder)
int g = 0;

float Bi = 2.0, delx = 1.0/node;

int main()
{
    // Function definitions
    void input(float *, float *, float *, float *);
    void tdma(float *, float *, float *, float *, float *);
    void output(float *);

    // Variable Declarations
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

    for(i = 0; i < node; i++) {
        ld[i] = 1.0;
        md[i] = -2.0;
        ud[i] = 1.0;
        rhs[i] = 0.0;
    }

    //ld[node-1] = 2.0;
    rhs[node-1] = -1;
    //md[node-1] = -(2 + (Bi * delx));
    //ud[0] = 2.0;
}

void tdma(float ld[], float md[], float ud[], float rhs[], float x[])
{
    int k;

    for(k = 1; k < node; k++) {
        ld[k] = ld[k] / md[k-1];
        md[k] = md[k] - (ld[k] * ud[k-1]);
    }

    for(k = 1; k < node; k++)
        rhs[k] = rhs[k] - (ld[k] * rhs[k-1]);

    x[node-1] = rhs[node-1] / md[node-1];
    for(k = node-2; k >= 0; k--)
        x[k] = (rhs[k] - (ud[k] * x[k+1])) / md[k];
}

void output(float a[])
{
    for(int i = 0; i < node; i++)
        cout << a[i] << " ";
}
