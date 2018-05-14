# include <iostream>
# define N 101
using namespace std;

void grid();
void input();
void conscalc();
void initialize();
void matrixS();
void matrixL();
int meltstart();
void phase();
void update();
void updatephase();
void TDMA();
void output();
void print();

double delX, alphaS, alphaL, BiS, BiL, deltS, deltL, t;
double ldS, ldL, thetaI, thetaF, thetaM, tI, tF, tM, SteS, SteL;
double thetaN[N], thetaO[N], fN[N], fO[N];
double l[N], m[N], u[N], r[N];

int main() {
    grid();
    input();
    initialize();
    for(int i = 0; i < 500; i++) {
        matrixS();
        TDMA();
        update();
//        phase();
//       if(meltstart() == 1) {
//           i--;
//            continue;
//        }
//        output();
        print();
        update();
//        updatephase();
    }
}

void grid() {
    delX = 1.0/(N-1);
}

void input() {
    t = 5.0;
    tF = 60.0;
    tM = 29.0;
    tI = 17.5;
    alphaS = 4.5530E-7;
    alphaL = 1.5746E-7;
    BiS = 0.1468;
    BiL = 0.3019;
    deltS = 4.5530E-3 * t;
    deltL = 1.5746E-3 * t;
    ldS = deltS/delX/delX;
    ldL = deltL/delX/delX;
    thetaF = 1.0;
    thetaM = 0.0;
    thetaI = (tI - tM) / (tF - tM);
    SteS = 0.2321;
    SteL = 0.3647;
}

void initialize() {
    for(int i = 0; i < N; i++) {
        thetaO[i] = thetaI;
        fN[i] = 0.0;
        fO[i] = 0.0;
    }
}

void matrixS() {
    l[0] = 0.0;
    m[0] = 1 + ldS;
    u[0] = -ldS;
    r[0] = (ldS*thetaO[1]) + ((1-ldS)*thetaO[0]);
    if(fN[0]>0 && fN[0]<1){
		m[0]=1;
		u[0]=0;
		r[0]=thetaM;
	}

    for(int i = 1; i < (N-1); i++) {
        l[i] = -ldS;
        m[i] = 2*(1+ldS);
        u[i] = -ldS;
        r[i] = (ldS*thetaO[i-1]) + (2*(1-ldS)*thetaO[i]) + (ldS*thetaO[i+1]);
        if(fN[i]>0 && fN[i]<1){
			l[i]=0;
			m[i]=1;
			u[i]=0;
			r[i]=thetaM;
		}
    }

    l[N-1] = -ldS;
    m[N-1] = ldS + 1 + (ldS*BiS*delX);
    u[N-1] = 0.0;
    r[N-1] = (ldS*thetaO[N-2]) + ((1-ldS-(ldS*BiS*delX))*thetaO[N-1]) + (2*BiS*delX*ldS);
    if(fN[N-1]>0 && fN[N-1]<1){
		l[N-1]=0;
		m[N-1]=1;	
		r[N-1]=thetaM;
	}
}

void matrixL() {
}

void phase() {
    double sum;
    sum = (SteS*ldS)*(2*thetaN[1] + 2*thetaM);
    if(sum > 0)
        fN[0] = fO[0] + sum;    

    for(int i = 1; i < N-1; i++) {
        sum = (SteS*ldS)*(thetaN[i+1] - 2*thetaM + thetaN[i-1]);
        if(sum > 0)
            fN[i] = fO[i] + sum;
   }
    sum = (SteS*ldS)*(thetaN[N-2] - 2*thetaM + thetaN[N-2] + 2*Bi*delX*(1-thetaN[N-1]));
    if(sum > 0)
        fN[N-1] = fO[N-1] + sum;  
}

int meltstart(){
    int flag = 0;
    for(int i = 1; i < N; i++) {
        if(thetaN[i] >= thetaM && thetaO[i] < thetaM) {
            fN[i] = fN[i] - SteS*(thetaM - thetaO[i]);
            flag = 1;
        }
    }
    return flag;
}

void TDMA() {
 //   cout << "\n \n TDMA: \n\n";
 //   printf("         %+0.5lf %+0.5lf\t%+0.5lf\n", m[0], u[0], r[0]);
 //   for( int i = 1; i < N-1; i++)
//        printf("%+0.5lf %+0.5lf %+0.5lf \t%+0.5lf\n", l[i], m[i], u[i], r[i]);
//
//    printf("%+0.5f %+0.5f\t\t%+0.5f\n\n", l[N-1], m[N-1], r[N-1]);

    for(int i = 1; i < N; i++) {
        l[i] = l[i] / m[i-1];
        m[i] = m[i] - (l[i]*u[i-1]);
        r[i] = r[i] - (l[i]*r[i-1]);
    }
    thetaN[N-1] = r[N-1]/m[N-1];

    for(int i = N-2; i >= 0 ; i--)
        thetaN[i] = (r[i] - u[i]*thetaN[i+1])/m[i];
}


void update() {
    for(int i = 0; i < N; i++)
        thetaO[i] = thetaN[i];
}

void output() {
    for(int i = 0; i < N; i++)
        cout << thetaN[i] << " ";
    cout << endl;
    cout << (thetaN[N-1] - thetaN[N-2]) / delX / (1 - thetaN[N-1]);
}

void print() {
    cout << thetaN[0]*(tF-tM)+tM << "\t" << fN[0] << endl;
}

void updatephase() {
    for(int i = 0; i < N; i++)
        fO[i] = fN[i];
}
