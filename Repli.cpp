/* To replicate the Zivkovic's Paper */
# include <iostream>
# include <math.h>
# define size 11
# define tstep 3
using namespace std;

float tm, k, h, rho, c, L, delt, delx, alpha, Fo, ti, tf;
float ld[size+1], md[size+1], ud[size+1], rs[size+1], t[size+1], to[size+1];
float fl[size+1], flo[size+1];

void constant();
void output(float*, int);
void initialize();
void matrix();
void update(float*, float*);
void tdma(float*, float*, float*, float*, float*);
int phase();
void tempcheck();
void check();

int main() {
    int i;
    constant();
    initialize();
    for(i = 1; i <= tstep; i++) {
        cout << "\nTemperature : \n";
        output(to, i);
//        cout << "Liquid Fraction : \n";
//        output(flo, i);
        matrix();
        tdma(ld, md, ud, rs, t);
//        if(phase() == 1) {
//            i--;
//            continue;
//       }
//        check();
//        tempcheck();
       update(to, t);
//       if(i%25 == 0)
//        cout << t[1] << "\n";
//        update(flo, fl);
    }

//    cout << "\nTemperature : \n";
//    output(t, i);
    
//    cout << "Liquid Fraction : \n";
//    output(fl, i);

    cout << "\nGrid Size: " << delx;
    return 0;
}

void tempcheck() {
    float ht;
    ht = -k*(t[size]-t[size-1])/delx/(t[size]-tf);
    cout <<"\nHeat transfer coefficient : " << ht;
    cout << "\n";
}


void constant() {
    tm = 29.0;
    ti = 18.0;
    tf = 60.0;
    k  = 1.09;
    rho = 1710;
    c = 1400;
    L = 187000;
    delt = 5.0;
    delx = 0.02 / (size-1);
    alpha = k / (rho*c);
    Fo = alpha * delt / (delx * delx);
    h = 16.0;
}

void initialize() {
    for(int i = 1; i <= size; i++) {
        to[i] = ti;
        t[i]  = ti;
        fl[i] = 0.0;
        flo[i]= 0.0;
    }
}

void matrix() {
    float ai, aj, ah;

    ah = -Fo;
    aj = -Fo;
    ai = 1 + (2 * Fo);

    for(int i = 1; i <= size; i++) {
        ld[i] = ah;
        md[i] = ai;
        ud[i] = aj;
        rs[i] = to[i];
    }

    ud[1] = ah + aj;

    ld[size] = ah + aj;
    md[size] = ai - (2 * h * (delx / k) * aj);
    rs[size] = to[size] - (2 * h * (delx / k) * tf * aj);

   /* for( int i = 1; i <= size; i++) {
        if(t[i] >= tm) {
            ld[i] = 0;
            md[i] = 1;
            ud[i] = 0;
            rs[i] = tm;
        }
    }*/
}

void tdma(float ld[], float md[], float ud[], float rs[], float x[]) {
    int k;

    cout << "\n \n TDMA: \n\n";
    printf("      %+0.2f %+0.2f\t%+0.2f\n", md[1], ud[1], rs[1]);

    for( int i = 2; i < size; i++)
        printf("%+0.2f %+0.2f %+0.2f \t%+0.2f\n", ld[i], md[i], ud[i], rs[i]);

     printf("%+0.2f %+0.2f\t\t%+0.2f\n\n", ld[size], md[size], rs[size]);


    for(k = 2; k <= size; k++) {
        ld[k] = ld[k] / md[k-1];
        md[k] = md[k] - (ld[k] * ud[k-1]);
        rs[k] = rs[k] - (ld[k] * rs[k-1]);
    }

    x[size] = rs[size] / md[size];

    for(k = size-1; k > 0; k--)
        x[k] = (rs[k] - (ud[k] * x[k+1])) / md[k];

}

int phase() {
    int flag = 0;

    if(t[1] >= tm) {
        fl[1] = flo[1] + ((k/(rho*L))*(delt/(delx*delx))*(t[2]-2*tm+t[2]));
        if(to[1] < tm) {
            fl[1] = fl[1] - (c/L)*(tm-to[1]);
            flag = 1;
            to[1] = tm;
        }
   }

    for(int i = 2; i <= size; i++) { 
        if(t[i] >= tm) {
            fl[i] = flo[i] + ((k/(rho*L))*(delt/(delx*delx))*(t[i-1]-(2*tm)+t[i+1]));
            if(to[i] < tm){
                fl[i] = fl[i]-(c/L)*(tm-to[i]);
                flag = 1;
                to[i] = tm;
            }
        }
    }
    
    if(t[size] >= tm) {
        fl[size] = flo[size] + (k/(rho*L))*(delt/(delx*delx))*(2*t[size-1]-2*tm-(2*h*delx/k*(t[size]-tf)));
        if(to[size] < tm){
            fl[size] = fl[size]-(c/L)*(tm-to[size]);
            flag = 1;
            to[size] = tm;
        }
    }

    return flag;
}

void update(float a[], float b[]) {
    for(int i = 1; i <= size; i++)
        a[i] = b[i];
}

void output(float a[], int step) {
    cout << "Time : "<< (step-1) * 5 <<" s\n\n";
    for(int i = 1; i <= size; i++)
        cout << a[i] << "\n";
    cout << "\n";
}

void check() {
    float temp;
    cout << "\n Checking temperature profile: \n";
    for(int i = 2; i <= size-1; i++) {
        temp = alpha * (t[i-1]-2*t[i]+t[i+1])/(delx*delx);
        temp-=(t[i]-to[i])/delt;
        cout << temp <<"\n";
    }
}


