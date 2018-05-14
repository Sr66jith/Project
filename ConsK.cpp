/* Constant K with phase change */
# include <iostream>
# define NODE 11
# define TIME 1
using namespace std;

float cS, cL, kS, kL, dS, dL, tM, b, h, tI, tF, L;
float aS, aL;
float BiS, BiL, X, T, z;
void constant();
void TDMA();
void matrix();
void output(float);

int main() {
    constant();
    output(aS);
    output(aL);
    output(BiS);
    output(BiL);
}

void constant() {
   cS = 1400.0;
   cL = 2200.0;
   kS = 1.09;
   kL = 0.53;
   dS = 1710.0;
   dL = 1530.0;
   L  = 187000;
   tM = 29.0;
   b  = 0.01;
   h  = 16.0;
   tF = 60.0;

   aS = kS/(dS*cS);
   aL = kL/(dL*cL);
   BiS = h*b/kS;
   BiL = h*b/kL;
}

void output(float a) {
    cout << a << " ";
}
