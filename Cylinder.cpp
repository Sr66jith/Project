/* Final run of my Project */
# include <iostream>
# include <stdio.h>
# define N 21
using namespace std;

void grid();
void matrix();
void tdma();
void phase();
void melt_start();
void melt_end();
void display();
void input();

double delR, ro, t, Tm, BiS, BiL, alphaS, alphaL, deltS, deltL;

int main() {
    grid();
    input();
}

void grid() {
    delR = (1.0 - 0.7072136)/(N-1);
}

void input() {
    ro = 28.28E-3;
    t = 5.0;
    Tm = 29.0;
    BiS = 0.41512;
    BiL = 0.85374;
    alphaS = 4.553E-7;
    alphaL = 1.5746E-7;
    deltS = (alpha_s*t)/(ro*ro);
    deltL = (alpha_l*t)/(ro*ro);
    lambdaS = deltS/(delR*delR);
    lambdaL = deltL/(delR*delR);
    deltaS = deltS/(2*delR);
    deltaL = deltL/(2*delR);
}

