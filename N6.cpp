# include <iostream> 
# include <stdio.h>
# define N 21
using namespace std;
void grid();void input();void conscalc();void initialize();void matrix();void tdma(int neqs);void slope();void l_entropy();void g_entropy();
void display();void copy_temp();void update_lf();int meltstart();void liquidfrac();int meltend();
double delX,alpha_s,alpha_l,Bi_s,Bi_l,delt_s,delt_l,t,lambda_s,lambda_l,theta_f,theta_m,theta_i,T_f,T_m,T_i,theta_inf;
int step;
double ste_l,ste_s;int pos;double delta_s,delta_l,gamma,L,b,ko_l,ko_s;
double a(int type,int pos);
double theta[N],theta_old[N],lf[N],lf_old[N],ldia[N],dia[N],udia[N],rhs[N],var[N],dtheta[N],s_l[N],s_g;


int main(){
	step=0; pos=-1;
	grid();input();initialize();
	for(int i=0;i<1440;i++){
		for(int j=1;j<1000;j++){
			slope();
			matrix();	
			tdma(N);
		//	copy_temp();
		}
		liquidfrac();
		if(meltstart()==1){
			i--;
			copy_temp();
			continue;
		}
		pos=meltend();
		update_lf();
		copy_temp();
		l_entropy();
		g_entropy();
		display();
	}
}
void grid(){
	delX=1.0/(double)(N-1);
}
void input(){
	t=5;alpha_s=4.553E-7;alpha_l=1.5746E-7;
	Bi_s=0.1468;Bi_l=0.3019;
	delt_s=(4.5530E-3)*t;delt_l=(1.5746E-3)*t;
	lambda_s=delt_s/delX/delX;lambda_l=delt_l/delX/delX;
	T_f=60;T_m=29;T_i=17.5;ste_s=0.2321;ste_l=0.3647;
	theta_f=1;theta_m=0;theta_i=(T_i-T_m)/(T_f-T_m);
	delta_s=delt_s/delX/2.0;gamma=0.0;delta_l=delt_l/delX/2.0;
	theta_inf=T_m/(T_f-T_m);
	L=187000;b=0.01;ko_l=0.53;ko_s=1.09;
}
void initialize(){
	for(int i=0;i<N;i++){
		theta_old[i]=theta_i;
		lf_old[i]=lf[i]=0;
	}
}
void matrix(){
	dia[0]=1+(2*a(1,0)*lambda_s);
	udia[0]=-2*lambda_s*a(1,0);
	rhs[0]=theta_old[0];
	if((lf[0]>0 && lf[0]<1)||theta[0]>0){
		dia[0]=1;
		udia[0]=0;
		rhs[0]=theta_m;
	}
	if(lf[0]>=1){
		dia[0]=1+(2*a(1,0)*lambda_l);
		udia[0]=-2*lambda_l*a(1,0);
		rhs[0]=theta_old[0];
	}		
	if(pos==0){
	//	system("Pause");
		rhs[0]-=(1/ste_l)*(1-lf_old[0]);
	}
	for(int i=1;i<N-1;i++){
		ldia[i]=(-a(1,i)*lambda_s)+(a(2,i)*delta_s);
		dia[i]=1+(2*lambda_s*a(1,i));
		udia[i]=(-a(1,i)*lambda_s)-(a(2,i)*delta_s);
		rhs[i]=theta_old[i];
		if((lf[i]>0 && lf[i]<1)||theta[i]>0){
			ldia[i]=0;
			dia[i]=1;
			udia[i]=0;
			rhs[i]=theta_m;
		}
		if(lf[i]>=1){
			ldia[i]=(-a(1,i)*lambda_l)+(a(2,i)*delta_l);
			dia[i]=1+(2*lambda_l*a(1,i));
			udia[i]=(-a(1,i)*lambda_l)-(a(2,i)*delta_l);
			rhs[i]=theta_old[i];
		}

	}
	if(pos>0 && pos<N-1){
		//system("Pause");
		rhs[pos]-=(1/ste_l)*(1-lf_old[pos]);
	}
	ldia[N-1]=-2*a(1,N-1)*lambda_s;
	dia[N-1]=1+(2*lambda_s*a(1,N-1))+((a(1,N-1)*lambda_s+a(2,N-1)*delta_s)*2*delX*Bi_s/a(1,N-1));
	rhs[N-1]=theta_old[N-1]+((2*Bi_s*delX*(a(1,N-1)*lambda_s+a(2,N-1)*delta_s))/a(1,N-1));
	if((lf[N-1]>0 && lf[N-1]<1)||theta[N-1]>0){
		ldia[N-1]=0;
		dia[N-1]=1;	
		rhs[N-1]=theta_m;
	}
	if(lf[N-1]>=1){
		ldia[N-1]=-2*a(1,N-1)*lambda_l;
		dia[N-1]=1+(2*lambda_l*a(1,N-1))+((a(1,N-1)*lambda_l+a(2,N-1)*delta_l)*2*delX*Bi_l/a(1,N-1));
		rhs[N-1]=theta_old[N-1]+((2*Bi_l*delX*(a(1,N-1)*lambda_l+a(2,N-1)*delta_l))/a(1,N-1));
		}
	if(pos==N-1){
	//	system("Pause");
		rhs[N-1]-=(1/ste_l)*(1-lf_old[N-1])+((a(1,N-1)*lambda_l+a(2,N-1)*delta_l)*2*delX*Bi_l/a(1,N-1));
	}
	pos=-1;
}	
void tdma(int neqs){
	//	printf("%lf\t%lf\t%lf\n",dia[0],udia[0],rhs[0]);
	//for(int i=1;i<N-1;i++)
	//	printf("%lf\t%lf\t%lf\t%lf\n",ldia[i],dia[i],udia[i],rhs[i]);
	//printf("%lf\t%lf\t%lf\n",ldia[N-1],dia[N-1],rhs[N-1]);
    for(int k=1;k<neqs;k++){
        double m=ldia[k]/dia[k-1];
        dia[k]-=m*udia[k-1];
        rhs[k]-=m*rhs[k-1];
    }
    var[neqs-1]=rhs[neqs-1]/dia[neqs-1];
    for(int k=neqs-2;k>=0;k--)
        var[k]=(rhs[k]-(udia[k]*var[k+1]))/dia[k];
    for(int i=0;i<N;i++)
    	theta[i]=var[i];
}
void liquidfrac(){
	double sum=a(1,0)*ste_s*2*lambda_s*(theta[1]-theta_m);
	if(sum>0){
		lf[0]=lf_old[0]+sum;
		if(lf[0]>1)lf[0]=1;
	}
	for(int i=1;i<N-1;i++){
		sum=((a(1,i)*ste_s*lambda_s)*(theta[i+1]-(2*theta_m)+theta[i-1]))+(a(2,i)*ste_s*delta_s*(theta[i+1]-theta[i-1]));
		if(sum>0)
			lf[i]=lf_old[i]+sum;
		if(lf[i]>1)lf[i]=1;
	}
	sum=((2*a(1,N-1)*ste_s*lambda_s)*(theta[N-2]-(theta_m*(1+(Bi_s*delX/a(1,N-1))))+(delX*Bi_s/a(1,N-1))))+(2*a(2,N-1)*ste_s*delta_s*Bi_s*delX/a(1,N-1)*(1-theta_m));
	if(sum>0)
		lf[N-1]=lf_old[N-1]+sum;
	if(lf[N-1]>1)lf[N-1]=1;
}

void update_lf(){
	for(int i=0;i<N;i++)
		lf_old[i]=lf[i];
}
int meltstart(){
	int flag=0;
	if(theta[0]>0 && theta_old[0]<0){
		lf[0]=lambda_s*ste_s*2*a(1,0)*(theta[1]-theta_m)-(ste_s*(theta_m-theta_old[0]));
		flag=1;
	}
	for(int i=1;i<N-1;i++)
		if(theta[i]>=0 && theta_old[i]<0){
			lf[i]=lambda_s*ste_s*a(1,i)*(theta[i-1]-(2*theta_m)+theta[i+1])-(ste_s*(theta_m-theta_old[i]))+(a(2,i)*ste_s*delta_s*(theta[i+1]-theta[i-1]));
			flag= 1;
		}
	if(theta[N-1]>=0 && theta_old[N-1]<0){
		lf[N-1]=(lambda_s*ste_s*2*a(1,N-1)*(theta[N-1]-((1+(Bi_s*delX/a(1,N-1)))*theta_m)+(Bi_s*delX/a(1,N-1))))-(ste_s*(theta_m-theta[N-1]))+(a(2,N-1)*ste_s*delta_s*Bi_s*delX*(1-theta_m)/a(1,N-1));
		flag=1;
	}
	return flag;		
}
int meltend(){
	int position=-1;
	for(int i=0;i<N;i++)
		if(lf_old[i]<1 && lf[i]>=1)
			position=i;
	return position;
}
	
void display(){
    //double pos;
    //pos=0.0;
	//if(step==80){
	//	for(int i=0;i<N;i++){
	//		printf("%0.3lf %0.3lf\n",pos*10,theta[i]*(T_f-T_m)+T_m);
     //       pos+=delX;
      //  }
	//}

	/*
	for(int i=0;i<N;i++)
		printf("%0.3lf ",lf[i]);
	printf("\n");
	for(int i=0;i<N;i++)
		printf("%0.3lf ",s_l[i]);
	printf("\n\n");
	}*/
	if(step%5==0)
	printf("%lf %0.3lf\n",((step*5)+5)/60.0,(theta[0]*(T_f-T_m)+T_m));
    step++;
	//printf("\t::%0.6lf\n",s_l[N-1]);
	//if(step%100==0){
	//	for(int i=0;i<N;i++)
	//		printf("%0.3lf ",s_l[i]);
	//	printf("\n");
	//}
}	
void copy_temp(){
	for(int i=0;i<N;i++)
		theta_old[i]=theta[i];
}

double a(int type,int pos){
	if(type==1)
		if(lf[pos]>=1)
			return 1+gamma*(theta[pos]-theta_m);
		else
			return 1+gamma*(theta[pos]-theta_i);
	if(type==2)
		return gamma*dtheta[pos];
}
		
void slope(){
    dtheta[0]=((-3*theta_old[0])+(4*theta_old[1])+(-theta_old[2]))/2.0/delX;
    for(int i=1;i<N-1;i++)
        dtheta[i]=(theta_old[i+1]-theta_old[i-1])/2.0/delX;
    dtheta[N-1]=((3*theta_old[N-1])+(-4*theta_old[N-2])+(theta_old[N-3]))/2.0/delX;
}	
	
void l_entropy(){
	//for(int i=0;i<N;i++)
	//	s_l[i]=0;
	for(int i=0;i<N;i++){
	//	if(lf[i]<1){
			s_l[i]=(a(1,i)*dtheta[i]*dtheta[i]/(theta[i]+theta_inf));//*(1-lf[i]);
			//s_l[i]+=lf[i]*L*b*b/ko_l/T_m;
		//}
		//if(lf[i]>=1)
		//	s_l[i]+=(a(1,i)*dtheta[i]*dtheta[i]/(theta[i]+theta_inf));
		
	}
}

void g_entropy(){
	double area=0;
	for(int i=0;i<N-1;i++)
		area+=0.5*(s_l[i]+s_l[i+1])*delX;
	s_g=area;
}
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
