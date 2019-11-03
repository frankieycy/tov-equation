#include <stdio.h>
#include <math.h>

/*
# Note
Last: 2/11/2019
solve Tolman–Oppenheimer–Volkoff (TOV) equation with RK4
for neutron stars in equilibrium

# Abbrev
EOS: equation of state

# Definitions
G      grav const
K      EOS const
c      light speed
c2     light speed squared
m_sun  sun mass (SI)
gamma  adiabatic const

n      max iter steps
dr     radius increment
r      radius counter

**** Core params here ****
rho0   density at r=0 (center)
m_tot  star total mass
R      star radius
P[]    pressure at r
M[]    mass up to r
RHO[]  density at r
E[]    internal energy at r

# Units
SI units

# Outputs
data.csv   star data for a fixed rho0
mr.csv     M-R data for a range of rho0
*/

#define G      6.67e-11
#define K      1.11e-15
#define c      3.00e+8
#define c2     9.00e+16
#define m_sun  1.99e+30
#define gamma  2.75

#define n      5000000 // suff large st all iter finish
#define dr     1e-2
#define eps    1e-6 // small num

int    r;
double rho0,R,m_tot,
P[n],M[n],RHO[n],E[n];

double density(double p){
	/* EOS: rho=rho(p) */
	return pow(p/K,1/gamma);
}

double energy(double p, double rho){
	/* EOS: e=e(p,rho) */
	return p/(gamma-1)/rho;
}

double dp(double r,double p,double m,double rho,double e){
	/* TOV eq: dp/dr */
	// +eps to avoid division by zero at r=0
	return -dr*G*(rho*(1+e/c2)+p/c2)*
	(m+4*M_PI*pow(r,3)*p/c2)/(r+eps)/(r+eps-2*G*m/c2);
}

double dm(double r,double rho,double e){
	/* TOV eq: dm/dr */
	return dr*4*M_PI*r*r*rho*(1+e/c2);
}

void init(){
	/* initialisation */
	r=0;
	M[0]=0;
	RHO[0]=rho0;
	P[0]=K*pow(RHO[0],gamma);
	E[0]=energy(P[0],RHO[0]);
}

void step(){
	/* update P,M,RHO,E at r+dr with RK4 */
	double dp0,dp1,dp2,dp3,dp4,dm0,dm1,dm2,dm3,dm4,rho,e;

	dp1 = dp(r*dr,P[r],M[r],RHO[r],E[r]);
	dm1 = dm(r*dr,RHO[r],E[r]);

	rho = density(P[r]+dp1/2);
	e = energy(P[r]+dp1/2,rho);
	dp2 = dp((r+0.5)*dr,P[r]+dp1/2,M[r]+dm1/2,rho,e);
	dm2 = dm((r+0.5)*dr,rho,e);

	rho = density(P[r]+dp2/2);
	e = energy(P[r]+dp2/2,rho);
	dp3 = dp((r+0.5)*dr,P[r]+dp2/2,M[r]+dm2/2,rho,e);
	dm3 = dm((r+0.5)*dr,rho,e);

	rho = density(P[r]+dp3);
	e = energy(P[r]+dp3,rho);
	dp4 = dp((r+1)*dr,P[r]+dp3,M[r]+dm3,rho,e);
	dm4 = dm((r+1)*dr,rho,e);

	// final weighted increments
	dp0 = (dp1+2*dp2+2*dp3+dp4)/6;
	dm0 = (dm1+2*dm2+2*dm3+dm4)/6;

	if(P[r]+dp0>0 && M[r]+dm0>0){
		/* further iter still physical */
		// require +ve P[r+1],M[r+1]
		P[r+1]=P[r]+dp0;
		M[r+1]=M[r]+dm0;
		RHO[r+1]=density(P[r+1]);
		E[r+1]=energy(P[r+1],RHO[r+1]);
	}else
		/* next iter fails; radius found */
		R=r*dr;
	r++;
}

void iter(){
	/* iter till radius found */
	while(r<n && r*dr<R) step();
	m_tot=M[r-1];
}

void data(int num){
	/* print data */
	int radius=r,sep=radius/num;
	FILE *f=fopen("data.csv","w");
	fprintf(f,"radius,pressure,mass,density,energy\n");
	for(int r=0; r<radius; r+=sep){
		fprintf(f,"%.8e,%.8e,%.8e,%.8e,%.8e\n",
		r*dr,P[r],M[r],RHO[r],E[r]);
	}
	fclose(f);
}

void info(){
	/* report after iter */
	printf("rho0 : %.4e kg m^-3 | "
	       "R : %7.4lf km | "
	       "m_tot : %6.4lf * m_sun | "
	       "iter : %d\n",
	rho0,R/1e3,m_tot/m_sun,r);
}

void gen_mr(int a, int b){
	FILE *f=fopen("mr.csv","w");
	fprintf(f,"rho0,mass,radius\n");
	for(int rho=a; rho<b; rho++){
		R=2e6;
		rho0=rho*1e17;
		init();
		iter();
		info();
		fprintf(f,"%.6e,%.6e,%.6e\n",rho0,m_tot,R);
	}
	fclose(f);
}

int main(){
/*
	rho0=5e17; // density at r=0
	R=2e6; // radius upper bound (guess)

	init();
	iter();
	data(1000);
	info();
*/
	gen_mr(1,400);
	return 0;
}
