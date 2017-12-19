/* Code to run TDMA algorithm */

# include <iostream>
# define size 4
using namespace std;

int main()
{
    // Variable declaration
    float ld[size], md[size], ud[size], rhs[size], x[size];

    // Function declaration
    void input(float *, float *, float *, float *);
    void tdma(float *, float *, float *, float *, float *);
    void output(float *);

    // Function calls
    input(ld, md, ud, rhs);
    tdma(ld, md, ud, rhs, x);
    output(x);

    return 0;
}

void input(float ld[], float md[], float ud[], float rhs[])
{
    int i;

    //for(i = 1; i < size; i++) {
    //    ld[i] = -1.00;
    //}

    for(i = 0; i < size; i++) {
        md[i]  = 2.04;
        ld[i]  = -1.00;
        ud[i]  = -1.00;
        rhs[i] = 0.8;
    }

    //for(i = 0; i < size-1; i++) {
    //    ud[i] = -1.00;
    //}

    rhs[0] = 40.8;
    //for(i = 1; i < size-1; i++) {
    //   rhs[i] = 0.8;
    //}
    rhs[size-1] = 200.8;
}

void tdma(float ld[], float md[], float ud[], float rhs[], float x[])
{
    int k;
    
    // Forward Decomposition
    for(k = 1; k < size; k++) {
        ld[k] = ld[k] / md[k-1];
        md[k] = md[k] - (ld[k] * ud[k-1]);
    }

    // Forward Substituition
    for(k = 1; k < size; k++) 
        rhs[k] = rhs[k] - (ld[k] * rhs[k-1]);

    // Backward Substituition
    x[size-1] = rhs[size-1] / md[size-1];
    for(k = size-2; k >= 0; k--)
        x[k] = (rhs[k] - (ud[k] * x[k+1])) / md[k];
}

void output(float a[])
{
    for(int i = 0; i < size; i++)
        cout << a[i] << " ";
}
