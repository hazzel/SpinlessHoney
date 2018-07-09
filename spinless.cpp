//----------------------------------------------------------(includes)
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include <string>
#include <complex>
#include <chrono>
#include <boost/filesystem.hpp>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

using namespace std;

//------------------------------------------------(defined constants)
#define pi M_PI
complex<double> I (0.0,1.0);

#define N 20    //rough angles
#define Nb 2    //band index
#define No 2    //orbital index
#define Ns 1    //finer angle for integration within one patch
#define numring 1 // number of patch circles
#define Nq 6    //corners of the honeycomb times number of patch circles 6 x numring
#define Nt 6    //thread number (local: 4, cluster: 12)

//------------------------------------------------(system parameters)

//temperature
double T=1E-10; //T=1E-5;

//interaction parameters
double Vnn;
double Vnnn;
//double Vqah;

//hopping
double t = 1.0;

//chemical potential
double mu;
double mu_cutoff;

//---------------------------------------------------(flow parameters)
double laminit,lam;
double dlaminit=0.001,dlam;
double dlammin=1E-7;

//------------------------------------------------(radial integration)
double Nr=20;//steps in radial integration


//double rpatchint 0.01; //used in radint //obsolete

double rpatch = 0.01;

double bisectionprec=1E-10;
double mf=0.05;                 //temperature for cutoff function
double fnorm;                   //area of BZ and factor pi^2 for scaling of r and dr
double deltaEmin=1E-12; //1E-12;

double rup[N][Ns][Nq][Nb][2];   //limits of radial integration for 2 shells
double rd[N][Ns][Nq][Nb][2];    //limits of radial integration for 2 shells

//--------------------------------------------(integration parameters)
double summax=1E-12;    //largest change in one RG step
double summaxold;
bool stop=0;
bool tp0tag=1;
double vmax;
double vmaxstop=300;double vminstop=-300;

//----------------------------------------------------(gemoetry of BZ)
double rpatchvar[N][Nq];        //radial coordinate of representative patch point n=1,...,N and q=1,...,N
double angle[N][Nq];            //angle of patch point (n,q) with n=1,...,N and q=1,...,Nq
double sangle[N][Ns][Nq];       //finer angle
double rmaxinit[N][Ns];         //maximum value of radial coordinate for given direction
double patchx[N][Nq];           //x-coordinate of patch point (n,q)
double patchy[N][Nq];           //y-coordinate of patch point (n,q)
//double patchradx[N][Nq];
//double patchrady[N][Nq];
double xcorner[Nq];             //x-coordinate of corner q=1,...,Nq
double ycorner[Nq];             //y-coordinate of corner q=1,...,Nq

double deltax[3],deltay[3];			//nearest neighbor distances
double Dx[6],Dy[6];					//next nearest neighbor distances

double G1x; double G1y; double G2x; double G2y; double G3x; double G3y;	//next neighbor distances in reciprocal lattice

//-----------------------------------------------------------(vertex)
complex<double> V[N][Nq][Nb][N][Nq][Nb][N][Nq][Nb][Nb];  //band
complex<double> O[N][Nq][No][N][Nq][No][N][Nq][No][No];   //orbital
complex<double> dV[N][Nq][Nb][N][Nq][Nb][N][Nq][Nb][Nb];  //increment

//----------------------------------------------------(transformation)
complex<double> patchtrafo[N][Nq][No][Nb]; //transformation for patch points

//----------------------------------------------------------(bubbles)
double Lpp[N][Nq][N][Nq][N][Nq][Nb][Nb]; //3 momenta, 2 band indices
double Lph[N][Nq][N][Nq][N][Nq][Nb][Nb]; //3 momenta, 2 band indices

//---------------------------------------------(momentum conservation)
int impuls[N][Nq][N][Nq][N][Nq][2];//3 momenta, 4th one is computed as (n,q)

//-------------------------------------------------------------------
double looppmax; double loophmax;

//---------------------------------------------------------(functions)
#include "funktionen/geometry.cpp"
#include "funktionen/matrix.cpp"
#include "funktionen/energie.cpp"           //defines dispersion and cut-off functions
#include "funktionen/impulserhaltung.cpp"   //momentum conservation
#include "funktionen/initfunktionen.cpp"    //
#include "funktionen/loops.cpp"             //radial and fine-angle integration for bubbles
#include "funktionen/init.cpp"              //BZ geometry, initial condition
#include "funktionen/vstep.cpp"             //computes increment and new vertex function

//-----------------------------------------------------(main function)
int main(int argc, char *argv[])
{
	if (argc < 6)
		return 1;
	char folder[50];
	sprintf(folder,"data");
	double vmin = std::stod(argv[3]);
	double vdelta = std::stod(argv[4]);
	double vcnt = std::stod(argv[5]);
	double vnnn = std::stod(argv[6]);

	//for-loop for parameter scans
	for(int i=0; i<vcnt; ++i)
	{
		Vnn = vmin + vdelta * i;
		if (argc < 7)
			Vnnn = 0;
		else
			Vnnn = vnnn;
		mu = std::stod(argv[1]);
		mu_cutoff = std::stod(argv[2]);
		
		std::string dir_name("../data/mu="+std::to_string(mu)+"_Vmin="+std::to_string(vmin)
			+ "_Vmax="+std::to_string(vmin + vcnt*vdelta));
		boost::filesystem::path dir(dir_name);
		if(!(boost::filesystem::exists(dir)))
			boost::filesystem::create_directory(dir);
		
		std::ofstream aus(dir_name+"/ausgabe"+std::to_string(i)+".dat");
		std::ofstream flow_out(dir_name+"/flow_run_"+std::to_string(i)+".txt");
		flow_out << "RG step, lambda, vmax, summax, dsummax, loophmax, looppmax" << std::endl;

		aus<<"==========================================================="<<endl;
		aus<<"Nk: "<<N<<"; vmax: "<<vmaxstop<<endl;
		aus<<"==========================================================="<<endl;

		int count=0;
		
		std::ofstream reoout(dir_name+"/ReO"+std::to_string(i)+".dat");
		std::ofstream imoout(dir_name+"/ImO"+std::to_string(i)+".dat");
		std::ofstream revout(dir_name+"/ReV"+std::to_string(i)+".dat");
		std::ofstream imvout(dir_name+"/ImV"+std::to_string(i)+".dat");
		std::ofstream vzwout(dir_name+"/vzw"+std::to_string(i)+".dat");
		std::ofstream lout(dir_name+"/l"+std::to_string(i)+".dat");
		std::ofstream vdiagout(dir_name+"/V_diag"+std::to_string(i)+".dat");
		std::ofstream vflowout(dir_name+"/V_flow"+std::to_string(i)+".dat");
		std::ofstream lppout(dir_name+"/Lpp_flow"+std::to_string(i)+".dat");
		std::ofstream lphout(dir_name+"/Lph_flow"+std::to_string(i)+".dat");

		//set up everything:
		geometry();
		init();

		//output to ausgabe file
		aus<<endl;
		aus<<"==========================================================="<<endl;
		aus<<i<<endl;
		aus<<"t: "<<t<<"; mu: "<<mu<<"; V1: "<<Vnn<<"; V2: "<<Vnnn<<endl;
		aus<<"laminit: "<<laminit<<endl;
		aus<<"==========================================================="<<endl;

		//output to terminal
		cout<<endl;
		cout<<"==========================================================="<<endl;
		cout<<i<<endl;
		cout<<"t: "<<t<<"; mu: "<<mu<<"; V1: "<<Vnn<<"; V2: "<<Vnnn<<endl;
		cout<<"laminit: "<<laminit<<endl;
		cout<<"==========================================================="<<endl;

		//this is where the magic happens:
		while (stop==0)
		{
			count = count + 1;

			//output to ausgabe file
			aus<<"-----------------------------------------------------------"<<endl;
			aus<<"RG-step:"<<count<<endl;
			aus<<"lambda: "<<lam<<"; dlambda: "<<dlam<<endl;
			flow_out << count << "," << lam;

			//output to terminal
			cout<<"-----------------------------------------------------------"<<endl;
			cout<<"RG-step:"<<count<<endl;
			cout<<"lambda: "<<lam<<"; dlambda: "<<dlam<<endl;

			//write flow parameter to file with insertion operator
			lout<<lam<<" ";

			//-----------------------------------------------------(RG step)
			std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
			shellvec(); //set radial integration range in energy shell around lambda
			std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
			loopcalc(); //perform radial integration and fine-angle averaging of bubbles
			std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
			dvstep();   //compute increment
			std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
			vstep();    //compute vertex function
			std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
			std::cout << "time shellvec = " << std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0).count() << std::endl;
			std::cout << "time loopcalc = " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << std::endl;
			std::cout << "time dvstep = " << std::chrono::duration_cast<std::chrono::duration<float>>(t3 - t2).count() << std::endl;
			std::cout << "time vstep = " << std::chrono::duration_cast<std::chrono::duration<float>>(t4 - t3).count() << std::endl;
			//-----------------------------------------------------

			lam-=dlam;  //set new flow parameter
			dlam=min(0.05*lam,dlam/(2*summax)); //set new decrement of flow parameter  //min(0.05*lam,dlam/(2*summax));
			if (dlam<dlammin)dlam=dlammin;

			//store vmax in file
			vzwout<<vmax<<" ";

			//output to ausgabe file
			aus<<"summax: "<<summax<<" dsummax:"<<summax-summaxold<<endl;
			aus<<""<<endl;
			aus<<"loophmax "<< loophmax<<endl;
			aus<<"looppmax "<< looppmax<<endl;
			aus<<""<<endl;
			aus<<"Vmax: "    <<vmax<<endl;
			aus<<""<<endl;
			aus<<"-----------------------------------------------------------"<<endl;
			aus<<""<<endl;
			
			flow_out << "," << vmax << "," << summax << "," << summax-summaxold
							<< "," << loophmax << "," << looppmax << std::endl;

			//output to terminal
			cout<<"summax: "<<summax<<" dsummax:"<<summax-summaxold<<endl;
			cout<<""<<endl;
			cout<<"loophmax "<< loophmax<<endl;
			cout<<"looppmax "<< looppmax<<endl;
			cout<<""<<endl;
			cout<<"Vmax: "    <<vmax<<endl;
			cout<<""<<endl;
			
			cout<<"-----------------------------------------------------------"<<endl;
			cout<<""<<endl;
			
			if ((count-1) % 25 == 0)
			{
				trafozuorbital();
				
				for (int k1=0;k1<N;k1++)
				for (int q1=0;q1<Nq;q1++)
				for (int b1=0;b1<Nb;b1++)
					for (int b2=0;b2<Nb;b2++)
						for (int k3=0;k3<N;k3++)
						for (int q3=0;q3<Nq;q3++)
						for (int b3=0;b3<Nb;b3++)
							for (int b4=0;b4<Nb;b4++)
								vdiagout << real(O[k1][q1][b1][k1][(q1+Nq/2)%Nq][b2][k3][q3][b3][b4]) << " ";
				vdiagout << std::endl;
				
				vflowout << lam << " ";
				for (int k1=0;k1<N;k1++)
				for (int q1=0;q1<Nq;q1++)
				for (int b1=0;b1<Nb;b1++)
					for (int k2=0;k2<N;k2++)
					for (int q2=0;q2<Nq;q2++)
					for (int b2=0;b2<Nb;b2++)
						for (int k3=0;k3<N;k3++)
						for (int q3=0;q3<Nq;q3++)
						for (int b3=0;b3<Nb;b3++)
							for (int b4=0;b4<Nb;b4++)
								vflowout<<real(O[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4])<<" ";
				vflowout<<std::endl;
				
				for (int ks=0;ks<N;ks++)
					for (int qs=0;qs<Nq;qs++)
						for (int k1=0;k1<N;k1++)
							for (int q1=0;q1<Nq;q1++)
								for (int k2=0;k2<N;k2++)
									for (int q2=0;q2<Nq;q2++)
										for (int bs=0;bs<Nb;bs++)
											for (int bs2=0;bs2<Nb;bs2++)
												{
													lppout<<real(Lpp[ks][qs][k1][q1][k2][q2][bs][bs2])<<" ";
													lphout<<real(Lph[ks][qs][k1][q1][k2][q2][bs][bs2])<<" ";
												}
				lppout<<std::endl;
				lphout<<std::endl;
			}

			if(vmax>vmaxstop)break;

			if(count==401)break;

			if(lam < 0.0001*T)break;
		}


		std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
		//output to files
		trafozuorbital();

		std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
		for (int k1=0;k1<N;k1++)
			for (int q1=0;q1<Nq;q1++)
				for (int b1=0;b1<Nb;b1++)
					
					for (int k2=0;k2<N;k2++)
						for (int q2=0;q2<Nq;q2++)
							for (int b2=0;b2<Nb;b2++)
								
								for (int k3=0;k3<N;k3++)
									for (int q3=0;q3<Nq;q3++)
										for (int b3=0;b3<Nb;b3++)
											
											for (int b4=0;b4<Nb;b4++)
											{
												//-------------------------------------------------------------------
												//singlet channel:

												//vertex in band basis
												revout<<real(V[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4])<<" ";
												imvout<<imag(V[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4])<<" ";

												//vertex in orbital basis
												reoout<<real(O[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4])<<" ";
												imoout<<imag(O[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4])<<" ";
												//-------------------------------------------------------------------
											}
		std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
		std::cout << "time transform = " << std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0).count() << std::endl;
		std::cout << "time output = " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << std::endl;
	}
}
//-------------------------------------------------------------------
