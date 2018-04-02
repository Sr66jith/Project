/* For the cylindrical geometry */

# include <iostream>
# define size 11
using namespace std;

void input();
void constant();
void tdma();
void output();
void initialize();

float tf, ro, tm, L, rh, cs, cl, ks, kl;
float to[size+1], tn[size+1], ld[size+1], md[size+1], ud[size+1]; 

int main() {
    constant();


    return 0;
}

void constant() {
    tf = 60.0;
    ro = 0.01;
    tm = 29.9;
    L  = 187000;
    rh = 1710;
    cs = 1400;
    cl = 2200;
    ks = 1.09;
    kl = 0.53;
}
