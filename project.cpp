#include<iostream>
#include<iomanip>
double uae[1001][1001];	// for ee
double uan[1001][1001];	//for nne
double uas[1001][1001];	//for sse
double uaw[1001][1001];	//for w
double uap[1001][1001];	//for e
double vae[1001][1001];
double van[1001][1001];
double vas[1001][1001];
double vaw[1001][1001];
double vap[1001][1001];	
double bu[1001][1001];
double bv[1001][1001];		//coefficients for u and v equations
double u[1001][1001], uold[1001][1001];
double v[1001][1001], vold[1001][1001];
double uh[1001][1001];
double vh[1001][1001];	//u hat and v hat
double p[1001][1001];	//pressure
double bp[1001][1001];
double pae[1001][1001];
double pan[1001][1001];
double pas[1001][1001];
double paw[1001][1001];
double pap[1001][1001];	//coefficients for pressure 
double p_corr[1001][1001];
double p_corrap[1001][1001], p_corrae[1001][1001], p_corraw[1001][1001], p_corran[1001][1001], p_corras[1001][1001], p_corrb[1001][1001];	//coefficients for pressure correction equation. 
double rho, vis, g, patm; //density, viscocity, gravity, atmospheric pressure
double delx;	//cell size
int n; //grid size in any directon. Assuming square grid
int i,j;
int count; //counts the number of iterations of SIMPLER
//NOTE : CONVENTION : INCREASING i IS +VE X DIRECTION AND INCREASING J IS +VE Y DIRECTION. IE.,CODE STARTS AT SOUTH-WEST CORNER OF GRID. 
//STAGGERING CONVENTION : THE CORRESPONDING U VELOCITY CELL FOR A GIVEN PRESSURE CELL IS THE ONE SHIFTED TO THE RIGHT
//AND THE CORRESPONDING V VELOCITY CELL IS THE ONE SHIFTED TO THE NORTH.
//WE ASSUME A FLUID OF VISCOSITY 0.001 AND DENSITY OF 1000 kg/m^3
void gaussseidel(double ap[1001][1001],double ae[1001][1001],double aw[1001][1001],double an[1001][1001],double as[1001][1001],double b[1001][1001], double x[1001][1001]);
void SIMPLER();
int main()
{	double res;
	count = 0;
	patm = 100000;
	vis = 0.001;
	rho = 1000;
	n = 31; 	//grid size is 30*30. Can be changed as required. Max grid size is 1000. We 
	g = 9.81;
		for(i=1;i<n;i++)
	{
		for(j=1;j<n;j++)
		{
			u[i][j] = 50;
			v[i][j] = -1;
		}
	}//guessed velocity field
	do
	{
		for(i=1;i<n;i++)
		{
			for(j=1;j<n;j++)
			{
				uold[i][j] = u[i][j];
				vold[i][j] = v[i][j];
			}
		}
		SIMPLER();	//calls the SIMPLER function. The function will keep getting called till the values converge.
		res = 0;
		for(i=1;i<n;i++)
		{
			for(j=1;j<n;j++)
			{
				res = res + (uold[i][j] - u[i][j])*(uold[i][j] - u[i][j]) + (vold[i][j] - v[i][j])*(vold[i][j] - v[i][j]) ;
			}
		}
		
	}while(res>=0.01);
	std::cout<<"The number of iterations taken by SIMPLER is"<<count;
		std::cout<<"\n the residue is \t"<<res<<"\n";
	std::cout<<"The u velocities are\t \n \n \n";
	int k;
	for(j=1;j<n;j++)
	{
		for(i=1;i<n;i++)
		{
			k = n-j;
			std::cout<<std::setprecision(2)<<u[i][k]<<"\t";
		}
		std::cout<<"\n";
	}
	std::cout<<"\n \n \n \n The v velocities are \t \n \n \n ";
	for(j=1;j<n;j++)
	{
		for(i=1;i<n;i++)
		{	k = n - j;
			std::cout<<std::setprecision(2)<<v[i][k]<<"\t";
		}
		std::cout<<"\n";
	}
	std::cout<<"\n \n \n \n The pressures are \t \n \n \n";
	for(j=1;j<n;j++)
	{
		for(i=1;i<n;i++)
		{
			k = n-j;
			std::cout<<std::setprecision(2)<<p[i][k]<<"\t";
		}
		std::cout<<"\n";
	}
return 0;
}
void SIMPLER()	//CONTAINS STEPS 2 TO 7 OF SIMPLE ALGORITHM.
{
	// SIMPLER ALGORITHM FOR BACK STEP FLOW
	double c;
	c = 0.25;
	int i,j;

	for(j=1;j<=(n/4);j++)
	{
		u[0][j] = 0;
		v[0][j] = 0;
	}
	for(j= (n/4) + 1;j<=n;j++)
	{	
		u[0][j] = (c/2)*(j-25)*(100-j)*0.01;
		v[0][j] = 0;
	}//GENERATE VELOCITIES AT LEFT BOUNDARY - VELOCITY BOUNDARY CONDITION
	for(i=0;i<=n;i++)
	{
		u[i][0] = 0;
		v[i][0] = 0;
		u[i][100] = 0;
		v[i][100] = 0;
	}//GENERATE VELOCITIES AT TOP AND BOTTOM BOUNDARY - VELOCITY BOUNDARY CONDITION
	for(j=0;j<=n;j++)
	{
		v[n][j] = 0;
		p[n-1][j] = 0;
	}// dirichlet boundary condition for v and p at right face
		for(i=1;i< n; i++)
	{	
		for(j= 1; j<n;j++)
		{	//initialise coefficients for x momentum equation
			if(u[i+1][j] + u[i][j] >0)	//UPWINDING IN U
				uae[i][j] = vis;
			else
				uae[i][j] = vis - rho*delx*u[i+1][j];	
			if(u[i-1][j] + u[i][j] >0)
				uaw[i][j] = vis - rho*delx*u[i-1][j];
			else
				uaw[i][j] = vis;
			if(v[i][j] + v[i+1][j] > 0)			//UPWINDING IN V
				uan[i][j] = vis ;
			else
				uan[i][j] = vis -  rho*delx*(v[i][j] + v[i+1][j] + v[i][j+1] + v[i+1][j+1])/4 ;
			if(v[i][j-1] + v[i+1][j-1] > 0)
				uas[i][j] = vis - rho*delx*(v[i][j-1] + v[i+1][j-1] + v[i][j-2] + v[i+1][j-2])/4;
			else
				uas[i][j] = vis;
			bu[i][j] = (p[i][j] - p[i+1][j])*delx;
			uap[i][j] = uan[i][j] + uas[i][j] + uae[i][j] + uaw[i][j] + rho*delx*u[i+1][j] - rho*delx*u[i-1][j] + rho*delx*(v[i][j] + v[i+1][j] + v[i][j+1] + v[i+1][j+1])/4 - rho*delx*(v[i][j-1] + v[i+1][j-1] + v[i][j-2] + v[i+1][j-2])/4 ; 
			//initialise constants for y momentum equation
			if(v[i][j+1] + v[i][j] >0)

				van[i][j] = vis;
			else
				van[i][j] = vis - rho*delx*v[i][j+1];	//UPWINDING IN V
			if(v[i][j-1] + v[i][j] >0)
				vas[i][j] = vis - rho*delx*v[i][j-1];
			else
				vas[i][j] = vis;
			if(u[i][j+1] + u[i][j] > 0)
				vae[i][j] = vis;
			else
				vae[i][j] = vis - rho*delx*(u[i][j] + u[i][j+1] + u[i+1][j] + u[i+1][j+1])/4;
			if(u[i-1][j] + u[i-1][j+1] > 0)
				vaw[i][j] = vis - rho*delx*(u[i-2][j] + u[i-2][j+1] + u[i-1][j] + u[i-1][j+1])/4;
			else
				vaw[i][j] = vis;
			bv[i][j] = (p[i][j] - p[i][j+1])*delx - rho*g;
			vap[i][j] = van[i][j] + vas[i][j] + vae[i][j] + vaw[i][j] + rho*delx*v[i][j+1] - rho*delx*v[i][j-1] + rho*delx*(u[i][j] + u[i][j+1] + u[i+1][j] + u[i+1][j+1])/4 - rho*delx*(u[i-2][j] + u[i-2][j+1] + u[i-1][j] + u[i-1][j+1])/4 ; 
		}	
	}	
	for(j=0;j<=n;j++)
	{
		uap[n-1][j] = uap[n-1][j] - uae[n-1][j];
		uae[n][j] = 0; 
	}	//NEUMANN BOUNDARY CONDITION FOR U AT RIGHT SIDE FOR THE VELOCITY EQUATIONS

	for(i=1;i<n;i++)
	{
		for(j=1;j<n;j++)
		{	if(uap[i][j] != 0 && vap[i][j] != 0)
			{
			uh[i][j] = (uae[i][j]*u[i+1][j] + uaw[i][j]*u[i-1][j] + uas[i][j]*u[i][j-1] + uan[i][j]*u[i][j+1])/uap[i][j] ;
			vh[i][j] = (vae[i][j]*v[i+1][j] + vaw[i][j]*vaw[i-1][j] + vas[i][j]*v[i][j-1] + van[i][j]*v[i][j+1] - rho*g)/vap[i][j] ;	//STEP 2 OF SIMPLER
			}
			else
			{
				uh[i][j] = 0;
				vh[i][j] = 0;
			}
		}
	}
	for(i=1;i<n;i++)
	{
		for(j=1;j< n ;j++)
		{
			pae[i][j] = rho*(delx*delx)/uap[i][j];
			paw[i][j] = rho*(delx*delx)/uap[i-1][j];
			pan[i][j] = rho*(delx*delx)/uap[i][j];
			pas[i][j] = rho*(delx*delx)/uap[i][j-1];
			pap[i][j] = pae[i][j] + paw[i][j] + pan[i][j] + pas[i][j];
			bp[i][j] = delx*rho*(uh[i][j]-uh[i+1][j] + vh[i][j] - vh[i][j+1]) ; // The west face is i, the east face is i+1, north face is j+1, south face is j 
		}			//STEP 3 OF SIMPLER
	}
	for(j=0;j<=n;j++)
	{
		p[0][j] = patm;
		p[n][j] = patm;
	}	//DIRICHLET BOUNDARY CONDITIONS FOR PRESSURE AT TOP AND BOTTOM!!
	for(i=0;i<=n;i++)
	{
		pap[i][1] = pap[i][1] - pas[i][1];
		pas[i][1] = 0;
		pap[i][n-1] = pap[i][n-1] - pan[i][n-1];
		pan[i][n-1] = 0;		//NEUMANN BOUNDARY CONDIIONS FOR PRESSURE AT TOP AND BOTTOM!!
	}
	gaussseidel(pap,pae,paw,pan,pas,bp,p);	//STEP 3 OF SIMPLER - SOLVING THE PRESSURE EQUATION
	gaussseidel(uap,uae,uaw,uan,uas,bu,u);
	gaussseidel(vap,vae,vaw,van,vas,bv,v);		// u* and v* calculated. STEP 4 OF SIMPLER
	for(i=1;i< n;i++)
	{
		for(j=1;j<n;j++)
		{
			p_corrb[i][j] = rho*delx*(u[i][j] - u[i-1][j] + v[i][j] - v[i][j-1]);
			p_corrae[i][j] = rho*delx*delx/uap[i][j];
			p_corraw[i][j] = rho*delx*delx/uap[i-1][j];
			p_corran[i][j] = rho*delx*delx/vap[i][j];
			p_corras[i][j] = rho*delx*delx/vap[i][j-1];	//COEFFICIENTS FOR THE PRESSURE CORRECTION EQUATIONS
		}
	}
	gaussseidel(p_corrap,p_corrae,p_corraw,p_corran,p_corras,p_corrb,p_corr);
	u[i][j] = u[i][j] + (delx/uap[i][j])*(p[i][j] - p[i+1][j]);
	v[i][j] = v[i][j] + (delx/vap[i][j])*(p[i][j] - p[i][j+1]);
	count++;
}

void gaussseidel(double ap[1001][1001],double ae[1001][1001],double aw[1001][1001],double an[1001][1001],double as[1001][1001],double b[1001][1001], double x[1001][1001])
{	double xold[1001][1001];
	double err;
	do
	{	err = 0;
		for(i=1;i<n;i++)
		{
			for(j=1;j<n;j++)
			{
				xold[i][j] = x[i][j];
				x[i][j] = (x[i+1][j]*ae[i][j] + x[i-1][j]*aw[i][j] + x[i][j+1]*an[i][j] + x[i][j-1]*as[i][j])/ap[i][j];
				err = err + (xold[i][j] - x[i][j])*(xold[i][j] - x[i][j]);
			}
		}
	} while( err >= 0.001);
}