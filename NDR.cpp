/* Implemented Crank Nicolson Method */
# include <iostream>
# define NODE 5
# define TIME 2
using namespace std;

void constant();
void TDMA();
void initialize();
void output(float a[]);
void update();

float tI, t, c, tN[NODE+1], tO[NODE+1];
float u[NODE+1], m[NODE+1], l[NODE+1], r[NODE+1];

int main() {
    constant();
    while(t) {
        initialize();
        output(tO);
        output(r);
        TDMA();
        output(tN);
        update();
        t--;
    }
}

void constant() {
    tI = 0.0;
    t = TIME;
    c = 0.020875;  
    for(int i = 0; i <= NODE; i++) {
        tN[i] = tI;
        tO[i] = tI;
    }
    tO[0] = 100;
    tN[0] = 100;
    tO[NODE] = 50;
    tN[NODE] = 50;
}

void initialize() {
    l[1] = 0.0;
    m[1] = 2*(1+c);
    u[1] = -c;
    r[1] = 2*c*tO[0] + 2*(1-c)*tO[1] + c*tO[2]; 

    for(int i = 2; i < NODE; i++) {
        l[i] = -c;
        m[i] = 2*(1+c);
        u[i] = -c;
        r[i] = c*(tO[i+1]+tO[i-1]) + 2*tO[i]*(1-c);
    }
    l[NODE-1] = -c;
    m[NODE-1] = 2*(1+c);
    r[NODE-1] = 2*c*tO[NODE] + c*tO[NODE-2] + 2*tO[NODE-1]*(1-c);
}

void TDMA() {
    for(int i = 2; i <= NODE; i++) {
        l[i] = l[i] / m[i-1];
        m[i] = m[i] - (l[i]*u[i-1]);
        r[i] = r[i] - (l[i]*r[i-1]);
    }
    tN[NODE-1] = r[NODE-1]/m[NODE-1];

    for(int i = NODE-2; i > 0; i--) 
        tN[i] = (r[i] - u[i]*tN[i+1])/m[i];
}

void output(float a[]) {
    for(int i = 0; i <= NODE; i++)
        cout << a[i] << " ";
    cout << endl;
}

void update() {
    for(int i = 0; i <= NODE; i++) 
        tO[i] = tN[i];
}
