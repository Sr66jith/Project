# include <iostream>
# define SIZE 6
# define N 5
# define TIME 5
using namespace std;

float BiS, BiL, delX, delT, b, lambda;
float thetaN[NODE],thetaO[NODE];
float ld[SIZE], md[SIZE], ud[SIZE], rs[SIZE];

void constants();
void matrix();

int main() {
    constants();
//    cout << delT << endl << delX << endl << lambda;
}

void constants() {
    BiS = 0.1468;
    BiL = 0.3019;
    delX = 1.0 / (NODE-1);
    delT = 4.553E-3*TIME;
    lambda = delT / (delX*delX);
}

void matrix() {
    ld[0] = 0.0;
    md[0] = (1 + 2*lambda);
    ud[0] = -2*lambda;

    for(int i = 0; i < N; i++) {
        ld[i] = -lambda;
        md[i] = (1 + 2*lambda);
        ud[i] = -lambda;

