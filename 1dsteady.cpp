/* For 1 dimensional steady state conduction */

# include <iostream>
# define node 5
using namespace std;

// Geometry variable (0 for slab and 1 for cylinder)
int g = 0;

float Bi = 2.0, delx = 1.0/(node-1);

int main()
{
    // Function definitions
    void input(float *, float *, float *, float *);
    void tdma(float *, float *, float *, float *, float *);
    void output(float *);

    // Variable Declarations
    float ld[node+1], md[node+1], ud[node+1], rhs[node+1], temp[node+1];

    // Function calls
    input(ld, md, ud, rhs);
    tdma(ld, md, ud, rhs, temp);
    output(temp);
    return 0;
}

void input(float ld[], float md[], float ud[], float rhs[])
{
    int i;

    for(i = 2; i < node-1; i++) {
        ld[i] = (2*i - 3);
        md[i] = -4*(i-1);
        ud[i] = (2*i - 1);
        rhs[i] = 0.0;
    }

    rhs[node-1] = -2*node + 3;
}

void tdma(float ld[], float md[], float ud[], float rhs[], float x[])
{
    int k;

    x[1] = 0;

    for(k = 3; k < node; k++) {
        ld[k] = ld[k] / md[k-1];
        md[k] = md[k] - (ld[k] * ud[k-1]);
    }

    for(k = 3; k < node; k++)
        rhs[k] = rhs[k] - (ld[k] * rhs[k-1]);

    x[node-1] = rhs[node-1] / md[node-1];
    for(k = node-2; k > 1; k--)
        x[k] = (rhs[k] - (ud[k] * x[k+1])) / md[k];

    x[node] = 1;
}

void output(float a[])
{
    for(int i = 1; i <= node; i++)
        cout << a[i] << " ";
}
