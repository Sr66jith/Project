/#* 1 D conduction dominated phase change */
# include <iostream>
# define node 3
# define time 1
using namespace std;

// Geometry variable
int g = 0;

float Ste = 0.05, Bi = 2.0, theta_m = 0; 
float delx = 1.0/(node-1), delt = 0.5;

void output(float *);

int main()
{
    // Variable declarations
    float ld[node], md[node], ud[node], rhs[node];
    float theta[node], theta_old[node];
    float liq_frac_old[node], liq_frac[node];

    // Function declarations
    void input(float *, float *, float *, float *, float *, float *);
    void transient(float *, float *, float *, float *, float *, float *, float *, float *);

    // Function calls
    input(ld, md, ud, rhs, liq_frac, liq_frac_old);
    transient(ld, md, ud, rhs, theta, liq_frac, theta_old, liq_frac_old);
    output(theta);

    return 0;
}

// To assign the initial values for the variables
void input(float ld[], float md[], float ud[], float rhs[], float liq_frac[], float liq_frac_old[])
{
    int i;
    float c, a;

    c = -delt / (delx * delx);
    a = 1 - (2 * c);

    for(i = 0; i < node; i++) {
        ld[i] = c;
        md[i] = a;
        ud[i] = c;
        rhs[i] = 6.0;
        liq_frac[i] = 0.0;
        liq_frac_old[i] = 0.0;
    }

    ud[0] = 2 * c;
    md[node-1] = a + (2 * Bi * delx);
    ld[node-1] = 2 * c;
}

// To execute the transient function for determining the temp. and liquid
// fraction
void transient(float ld[], float md[], float ud[], float rhs[], float theta[], float liq_frac[], float theta_old[], float liq_frac_old[])
{
    void tdma(float *, float *, float *, float *, float *);
    void update(float *, float *, float *);
    void phase(float *, float *, float *, float *, float *);

    for(float i = 0; i <= time; i += delt) {
        tdma(ld, md, ud, rhs, theta);
        update(theta, rhs, theta_old);
//        phase(theta, liq_frac, ld, md, ud, theta_old, liq_frac_old);
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
void update(float x[], float rhs[], float theta_old[])
{
    for(int i = 0; i < node; i++) {
        theta_old[i] = rhs[i];
        rhs[i] = x[i];
    }
    output(theta_old);
    cout << "\n";
}

// To update the liquid fraction term
void phase(float theta[], float liq_frac[], float ld[], float md[], float ud[], float theta_old[], float liq_frac_old[])
{
    int i;
    float c = (Ste * delt) / (delx * delx);
    float c1 = -delt / (delx * delx);

    for(i = 1; i < (node-1); i++) {

        if((theta[i] >= theta_m) && (theta_old[i] < theta_m)) { //Checks start of melting
            liq_frac[i] = liq_frac_old[i] - (delt * Ste * (theta[i-1] - 2 * theta_m + theta[i+1])) + Ste * (theta_m - theta_old[i]);
        }
        else {
            liq_frac[i] = liq_frac_old[i] + c * (theta[i+1] - (2 * theta_m) + theta[i-1]);
        }

/*        if(theta[i] <= theta_m) { // Needs checking
            liq_frac[i] = liq_frac[i] - 2 * c * theta_m + c * theta[i+1] + c * theta[i-1];

            if((liq_frac[i] < 1) && (liq_frac[i] > 0)) {
                ld[i] = 0.0;
                md[i] = 1.0;
                ud[i] = 0.0;
                theta[i] = theta_m; // Needs checking
            }
            else { // Did not consider the variation in coefficients when node = 1 and node = n
                ld[i] = c1;
                md[i] = 1 - 2 * c1;
                ud[i] = c1;
            }
        }
*/
    }
}


// To output any array of size node
void output(float a[])
{
    for(int i = 0; i < node; i++)
        cout << a[i] << " ";
}
