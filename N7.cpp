# include <iostream>
# include <stdio.h>
# define N 11
using namespace std;
void grid();void input();void conscalc();void initialize();void matrix();void tdma(int neqs);void slope();void l_entropy();void g_entropy();
void display(int step);void copy_temp();void update_lf();int meltstart();void liquidfrac();int meltend();

double delX,alpha_s,alpha_l,Bi_s,Bi_l,delt_s,delt_l,t,lambda_s,lambda_l,theta_f,theta_m,theta_i,T_f,T_m,T_i,step,theta_inf;
double ste_l,ste_s;int pos;double delta_s,delta_l,gamma,L,b,ko_l,ko_s,delR,Ro;
double a(int type,int pos);double A(int type,int pos);
double theta[N],theta_old[N],lf[N],lf_old[N],ldia[N],dia[N],udia[N],rhs[N],var[N],dtheta[N],s_l[N],s_g,Ri[N],psi_s[N],psi_l[N];


int main(){
	step=0; pos=-1;
	grid();input();initialize();
	//printf("**%f**",a(2,1));
	for(int i=0;i<2000;i++){
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
		//l_entropy();
		//g_entropy();
		display(step++);
	}
}
void grid(){
	delR=(1.0-0.707213578)/(double)(N-1);
	//cout <<delR <<"!!";
}
void input(){
	t=5;alpha_s=4.553E-7;alpha_l=1.5746E-7;ko_l=0.53;ko_s=1.09;
	Ro=0.01414;
	Bi_s=16.0*Ro/ko_s;
	Bi_l=16.0*Ro/ko_l;
	delt_s=alpha_s*t/Ro/Ro;delt_l=alpha_l*t/Ro/Ro;
	T_f=60;T_m=29;T_i=17.5;ste_s=0.2321;ste_l=0.3647;
	theta_f=1;theta_m=0;theta_i=(T_i-T_m)/(T_f-T_m);
	delta_s=delt_s/delR/2.0;gamma=0;delta_l=delt_l/delR/2.0;
	theta_inf=T_m/(T_f-T_m);lambda_s=delt_s/delR/delR;lambda_l=delt_l/delR/delR;
	L=187000;b=0.00828;
	for(int i=0;i<N;i++)
		Ri[i]=0.707213578+(delR*i);
	for(int i=0;i<N;i++){
		psi_s[i]=delt_s/2.0/Ri[i]/delR;
		psi_l[i]=delt_l/2.0/Ri[i]/delR;
	}
	//cout <<lambda_s<<"))";
}
void initialize(){
	for(int i=0;i<N;i++){
		theta_old[i]=theta_i;
		lf_old[i]=lf[i]=0;
	}
}
void matrix(){
	dia[0]=A(3,0)+(A(2,0)*Bi_s*delR/a(1,0));
	udia[0]=-A(1,0)-A(2,0);
	rhs[0]=theta_old[0]+(A(2,0)*Bi_s*delR/a(1,0));
	//printf("$$%f$$",a(1,0));
	if((lf[0]>0 && lf[0]<1)||theta[0]>0){
		dia[0]=1;
		udia[0]=0;
		rhs[0]=theta_m;
	}
	if(lf[0]>=1){
		dia[0]=A(3,0)+(A(2,0)*Bi_l*delR/a(1,0));
		udia[0]=-A(1,0)-A(2,0);
		rhs[0]=theta_old[0]+(A(2,0)*Bi_l*delR/a(1,0));
	}		
	if(pos==0){
	//	system("Pause");
		rhs[0]-=((1/ste_l)*(1-lf_old[0]))+(A(2,0)*Bi_l*delR/a(1,0));
	}
	for(int i=1;i<N-1;i++){
		ldia[i]=-A(2,i);
		dia[i]=A(3,i);
		udia[i]=-A(1,i);
		rhs[i]=theta_old[i];
		if((lf[i]>0 && lf[i]<1)||theta[i]>0){
			ldia[i]=0;
			dia[i]=1;
			udia[i]=0;
			rhs[i]=theta_m;
		}
		if(lf[i]>=1){
			ldia[i]=-A(2,i);
			dia[i]=A(3,i);
			udia[i]=-A(1,i);
			rhs[i]=theta_old[i];
		}

	}
	if(pos>0 && pos<N-1){
		//system("Pause");
		rhs[pos]-=(1/ste_l)*(1-lf_old[pos]);
	}
	ldia[N-1]=-A(1,N-1)-A(2,N-1);
	dia[N-1]=A(3,N-1);
	rhs[N-1]=theta_old[N-1];
	if((lf[N-1]>0 && lf[N-1]<1)||theta[N-1]>0){
		ldia[N-1]=0;
		dia[N-1]=1;	
		rhs[N-1]=theta_m;
	}
	if(lf[N-1]>=1){
		ldia[N-1]=-A(1,N-1)-A(2,N-1);
		dia[N-1]=A(3,N-1);
		rhs[N-1]=theta_old[N-1];
		}
	if(pos==N-1){
	//	system("Pause");
		rhs[N-1]-=((1/ste_l)*(1-lf_old[N-1]));
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
	double sum1=(A(1,0)+A(2,0))*theta[1];
	double sum2=(A(1,0)+A(2,0)+(A(2,0)*Bi_s*delR/a(1,0)))*theta[0];
	double sum3=A(2,0)*Bi_s*delR/a(1,0);
	double sum=ste_s*(sum1-sum2+sum3);
	if(sum>0 && theta[0]>=0){
		//printf("@@%lf@@%lf@@%lf@@",sum1,sum2,sum3);
		lf[0]=lf_old[0]+sum;
		if(lf[0]>1)lf[0]=1;
	}
	for(int i=1;i<N-1;i++){
		sum1=A(1,i)*theta[i+1];
		sum2=(A(1,i)+A(2,i))*theta[i];
		sum3=A(2,i)*theta[i-1];
		sum=ste_s*(sum1-sum2+sum3);
		if(sum>0 && theta[i]>=0)
			lf[i]=lf_old[i]+sum;
		if(lf[i]>1)lf[i]=1;
	}
	sum1=(A(1,N-1)+A(2,N-1))*theta[N-2];
	sum2=(A(1,N-1)+A(2,N-1))*theta[N-1];
	sum=ste_s*(sum1-sum2);
	if(sum>0 && theta[N-1]>=0)
		lf[N-1]=lf_old[N-1]+sum;
	if(lf[N-1]>1)lf[N-1]=1;
}

void update_lf(){
	for(int i=0;i<N;i++)
		lf_old[i]=lf[i];
}
int meltstart(){
	int flag=0;
	double sum1,sum2,sum3,sum;
	if(theta[0]>0 && theta_old[0]<0){
		sum1=(A(1,0)+A(2,0))*theta[1];
		sum2=(A(3,0)+(A(2,0)*Bi_s*delR/a(1,0)))*theta[0];
		sum3=A(2,0)*Bi_s*delR/a(1,0);
		sum=ste_s*(sum1-sum2+sum3+theta_old[0]);
		lf[0]=sum;
		flag=1;
	}
	for(int i=1;i<N-1;i++)
		if(theta[i]>=0 && theta_old[i]<0){
			sum1=A(1,i)*theta[i+1];
			sum2=A(3,i)*theta[i];
			sum3=A(2,i)*theta[i-1];
			sum=ste_s*(sum1-sum2+sum3+theta_old[i]);
			lf[i]=sum;
			flag= 1;
		}
	if(theta[N-1]>=0 && theta_old[N-1]<0){
		sum1=(A(1,N-1)+A(2,N-1))*theta[N-2];
		sum2=(A(3,N-1))*theta[N-1];
		sum=ste_s*(sum1-sum2+theta_old[N-1]);
		lf[N-1]=sum;
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
	
void display(int step){
/*	if(step%10==0){
		for(int i=0;i<N;i++)
			printf("%+0.3lf ",theta[i]);

	printf("\n");
	for(int i=0;i<N;i++)
		printf("%0.3lf ",lf[i]);
	}
	/*printf("\n");
	for(int i=0;i<N;i++)
		printf("%0.3lf ",s_l[i]);
	printf("\n\n");
	}*/
	if(step%10==0)
	//printf("%0.3lf\n",theta[0]);
	//printf("\t::%0.6lf\n",s_l[N-1]);
	printf("%d\t%lf\n",step,theta[N-1]);
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
	if(type==2){
		//cout << "$$"<<dtheta[pos]<<"$$";
		return gamma*dtheta[pos];
	}
}
double A(int type,int pos){
	if(type==1)
		if(lf[pos]>=1)
			return (a(1,pos)*lambda_l)+(a(2,pos)*delta_l)+(a(1,pos)*psi_l[pos]);
		else{
			//cout <<"##" <<a(2,pos)<<"##";
			return (a(1,pos)*lambda_s)+(a(2,pos)*delta_s)+(a(1,pos)*psi_s[pos]);
		}
	if(type==2)
		if(lf[pos]>=1)
			return (a(1,pos)*lambda_l)-(a(2,pos)*delta_l)-(a(1,pos)*psi_l[pos]);
		else{
			//cout <<"&&"<<a(2,pos)<<"&&";
			return (a(1,pos)*lambda_s)-(a(2,pos)*delta_s)-(a(1,pos)*psi_s[pos]);
		}
	if(type==3)
		return 1+A(1,pos)+A(2,pos);
}
		
void slope(){
    dtheta[0]=((-3*theta_old[0])+(4*theta_old[1])+(-theta_old[2]))/2.0/delR;
    for(int i=1;i<N-1;i++)
        dtheta[i]=(theta_old[i+1]-theta_old[i-1])/2.0/delR;
    dtheta[N-1]=((3*theta_old[N-1])+(-4*theta_old[N-2])+(theta_old[N-3]))/2.0/delR;
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
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
