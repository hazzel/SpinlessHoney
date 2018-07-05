double find_rmax(int q, double ang)
{
	double rmax;
	double phi = ang - pi/3.*q;
	if (phi < pi/3.)
		rmax = xcorner[1]/std::cos(phi);
	else
		rmax = xcorner[1]/std::cos(2.*pi/3.-phi);
	return rmax;
}

void find_fs_point(int q, double ang, double nb, double& x, double& y)
{
	double rmax = find_rmax(q, ang);
	double es=0;
	double r_fs = bisection(rmax, 0, es, ang, q, nb);
	x = xcorner[q%(Nq/numring)] + r_fs*cos(ang);
	y = ycorner[q%(Nq/numring)] + r_fs*sin(ang);
}

void swap_patch_points()
{
	for (int q=0;q<Nq;q++)
		for (int k=0;k<N/2;k++)
		{
			std::swap(angle[k][q], angle[N-1-k][q]);
			std::swap(sangle[k][0][q], sangle[N-1-k][0][q]);
			std::swap(rmaxinit[k][0], rmaxinit[N-1-k][0]);
			std::swap(patchx[k][q], patchx[N-1-k][q]);
			std::swap(patchy[k][q], patchy[N-1-k][q]);
		}
}

void set_patch_points_K()
{
	double nb = (mu < 0 ? 0. : 1.);
	for (int k=0;k<N;k++)
		for (int q=0;q<Nq;q++)
		{
			sangle[0][0][q] = 0.;
			patchx[k][q] = 0.;
			patchy[k][q] = 0.;
		}
	for (int q=0;q<Nq;q++)
	{
		sangle[0][0][q] = pi/3.*q + 2.*pi/3.;
		rmaxinit[0][0] = find_rmax(q, sangle[0][0][q]);
		find_fs_point(q, sangle[0][0][q], nb, patchx[0][q], patchy[0][q]);
	}
	for (int q=0;q<Nq;q++)
	{
		double total_length = 0.;
		int n_res = 100000;
		double ang = sangle[0][0][q];
		double x_0 = patchx[0][q], y_0 = patchy[0][q];
		for (int n = 1; n < n_res; ++n)
		{
			double x_1, y_1;
			ang -= 2.*pi/3. / n_res;
			find_fs_point(q, ang, nb, x_1, y_1);
			total_length += std::sqrt((x_1 - x_0) * (x_1 - x_0) + (y_1 - y_0) * (y_1 - y_0));
			x_0 = x_1;
			y_0 = y_1;
		}
		
		double partial_length = 0.;
		double interval = (std::abs(mu) < 1. ? total_length / (N-1) : total_length / N);
		int k = 1;
		ang = sangle[0][0][q];
		x_0 = patchx[0][q], y_0 = patchy[0][q];
		for (int n = 1; n < n_res; ++n)
		{
			double x_1, y_1;
			ang -= 2.*pi/3. / n_res;
			find_fs_point(q, ang, nb, x_1, y_1);
			partial_length += std::sqrt((x_1 - x_0) * (x_1 - x_0) + (y_1 - y_0) * (y_1 - y_0));
			x_0 = x_1;
			y_0 = y_1;
			if (partial_length >= interval)
			{
				sangle[k][0][q] = ang;
				rmaxinit[k][0] = find_rmax(q, ang);
				patchx[k][q] = x_1;
				patchy[k][q] = y_1;
				partial_length = 0.;
				++k;
			}
		}
		if (std::abs(mu) < 1.)
		{
			//Last point
			sangle[k][0][q] = sangle[0][0][q]-2.*pi/3.;
			rmaxinit[k][0] = find_rmax(q, sangle[k][0][q]);
			find_fs_point(q, sangle[k][0][q], nb, x_0, y_0);
			patchx[k][q] = x_0;
			patchy[k][q] = y_0;
		}
		
		//test distance
		//std::cout << std::sqrt(std::pow(patchx[k][q] - patchx[k-1][q], 2.) + std::pow(patchy[k][q] - patchy[k-1][q], 2.)) << std::endl;
		//std::cout << std::sqrt(std::pow(patchx[1][q] - patchx[0][q], 2.) + std::pow(patchy[1][q] - patchy[0][q], 2.)) << std::endl;
	}
	
	//swap_patch_points();
}

void set_patch_points_G()
{
	//Koordinaten von Gamma gesehen von corner[0]
	double x0=2.0/(3.0*sqrt(3));
	double y0=2.0/3.0;
	double l0=x0;

	//setze winkel und feinere Winkel
	for (int k=0;k<N;k++)
		for (int q=0;q<Nq;q++)
		{
			angle[k][q]= pi+pi/6. + pi/3.0*q + k*pi/3./N + pi/3./(2.*N);
			for (int s=0;s<Ns;s++)
				sangle[k][s][q]=angle[k][q]+((-Ns+1)+2*s)*pi/(3.0*N*Ns);
		}

	//Setze rmax
	for (int k=0;k<N/2;k++)	//half of the smaller angles (to avoid negative r values)
		for (int s=0;s<Ns;s++)	//fine angles
		{
			double x0=2.0/(3.0*sqrt(3)); // 1/2 of the y coordinate of b2 divided by cos(30 degree)
			double y0=2.0/3.0;           // 1/2 of the y coordinate of b2
			double phi = k*pi/3./N + pi/3./(2.*N);
			rmaxinit[k][s]=sqrt(1.5*1.5*x0*x0 + y0*y0/4.) / cos(phi);                //first half of corner
			rmaxinit[N-k-1][Ns-s-1]=sqrt(1.5*1.5*x0*x0 + y0*y0/4.) / cos(phi);		//second half of corner
		}
	
	//new
	//find rpatch for mu different from 0
	for (int k=0;k<N;k++)
		for (int q=0;q<Nq;q++)
		{
			double ang=angle[k][q];
			double es=0;
			double nb = (mu < 0 ? 0. : 1.);
			
			if (std::fmod(ang, pi/3.) < 1E-3)
				rpatchvar[k][q] = rmaxinit[k][(Ns-1)/2];
			else
				rpatchvar[k][q]=bisection(rmaxinit[k][(Ns-1)/2],0,es,ang,q,nb);
		}

	//set patchpoints
	for (int k=0;k<N;k++)
		for (int q=0;q<Nq;q++)
		{
			patchx[k][q]=rpatchvar[k][q]*cos(angle[k][q]);
			patchy[k][q]=rpatchvar[k][q]*sin(angle[k][q]);
		}
}

void set_patch_points_K_const()
{
	//Koordinaten von Gamma gesehen von corner[0]
	double x0=2.0/(3.0*sqrt(3));
	double y0=2.0/3.0;
	double l0=x0;

	//setze winkel und feinere Winkel
	for (int k=0;k<N;k++)
	for (int q=0;q<Nq;q++)
	{
		angle[k][q]= pi/3.0*q + (2*k+1)*pi/3/N;
		for (int s=0;s<Ns;s++)
			sangle[k][s][q]=angle[k][q]+((-Ns+1)+2*s)*pi/(3.0*N*Ns);
	}

	//Setze rmax
	for (int k=0;k<N/2;k++)	//half of the smaller angles (to avoid negative r values)
		for (int s=0;s<Ns;s++)	//fine angles
		{
			rmaxinit[k][s]=xcorner[1]/cos(sangle[k][s][0]);                //first half of corner
			rmaxinit[N-k-1][Ns-s-1]=xcorner[1]/cos(sangle[k][s][0]);		//second half of corner
		}

	//new
	//find rpatch for mu different from 0
	for (int k=0;k<N;k++)
		for (int q=0;q<Nq;q++)
			rpatchvar[k][q]=rpatch * (q / (Nq/numring) + 1);

	//set patchpoints
	for (int k=0;k<N;k++)
		for (int q=0;q<Nq;q++)
		{
			patchx[k][q]=xcorner[q%(Nq/numring)]+rpatchvar[k][q]*cos(angle[k][q]);
			patchy[k][q]=ycorner[q%(Nq/numring)]+rpatchvar[k][q]*sin(angle[k][q]);
		}
		
	swap_patch_points();
}

/*
void set_patch_points()
{
	//-------------------------------------------------------------------(lattice specific)
	//Koordinaten von Gamma gesehen von corner[0]
	double x0=2.0/(3.0*sqrt(3));
	double y0=2.0/3.0;
	double l0=x0;
	//-------------------------------------------------------------------



	//-------------------------------------------------------------------(lattice specific)
	//setze winkel und feinere Winkel
	for (int k=0;k<N;k++){
	for (int q=0;q<Nq;q++){
//		if(abs(mu)>1){
//			double anglemin;
//			double anglemax;
//			double r0;
//			double nb;

//			if(mu > 0){nb = 0;}
//			if(mu < 0){nb = 1;}

//			r0 = bisectionGM(y0,0,0,nb);
//			anglemin = atan(r0/l0);
//			anglemax = pi/3.0 - anglemin;

//			angle[k][q]= pi/3.0*q + anglemin + (2*k+1)*anglemax/N;

//		}

		if(abs(mu)<mu_cutoff){angle[N-1-k][q]= pi/3.0*q + (2*k+1)*pi/3/N;}
		else
			angle[k][q]= pi+pi/6. + pi/3.0*q + k*pi/3./N + pi/3./(2.*N);

		for (int s=0;s<Ns;s++){

			sangle[k][s][q]=angle[k][q]+((-Ns+1)+2*s)*pi/(3.0*N*Ns);

		}


		}
	}
	//-------------------------------------------------------------------



	//-------------------------------------------------------------------
	//Setze rmax
	for (int k=0;k<N/2;k++){	//half of the smaller angles (to avoid negative r values)
	for (int s=0;s<Ns;s++){	//fine angles

		if(abs(mu)<mu_cutoff)
		{
			rmaxinit[k][s]=xcorner[1]/cos(sangle[k][s][0]);                //first half of corner
			rmaxinit[N-k-1][Ns-s-1]=xcorner[1]/cos(sangle[k][s][0]);		//second half of corner
		}
		else
		{
			double x0=2.0/(3.0*sqrt(3)); // 1/2 of the y coordinate of b2 divided by cos(30 degree)
			double y0=2.0/3.0;           // 1/2 of the y coordinate of b2
			double phi = k*pi/3./N + pi/3./(2.*N);
			rmaxinit[k][s]=sqrt(1.5*1.5*x0*x0 + y0*y0/4.) / cos(phi);                //first half of corner
			rmaxinit[N-k-1][Ns-s-1]=sqrt(1.5*1.5*x0*x0 + y0*y0/4.) / cos(phi);		//second half of corner
		}

	}
	}
	//-------------------------------------------------------------------



	//new
	//-------------------------------------------------------------------
	//find rpatch for mu different from 0
		if(abs(mu) > 1E-6){
			for (int k=0;k<N;k++){
					for (int q=0;q<Nq;q++){

						double ang=angle[k][q];
						double es=0;

						if(mu>0)
						{
							if (std::fmod(ang, pi/3.) < 1E-3)
								rpatchvar[k][q] = rmaxinit[k][(Ns-1)/2];
							else
								rpatchvar[k][q]=bisection(rmaxinit[k][(Ns-1)/2],0,es,ang,q,0);
						}
						if(mu<0)
						{
							if (std::fmod(ang, pi/3.) < 1E-3)
								rpatchvar[k][q] = rmaxinit[k][(Ns-1)/2];
							else
								rpatchvar[k][q]=bisection(rmaxinit[k][(Ns-1)/2],0,es,ang,q,1);
						}
					}
			}
		}

		else {
			for (int k=0;k<N;k++){
					for (int q=0;q<Nq;q++){
						rpatchvar[k][q]=rpatch * (q / (Nq/numring) + 1);
					}
			}
		}

	//-------------------------------------------------------------------
	//set patchpoints
	for (int k=0;k<N;k++){
			for (int q=0;q<Nq;q++){

				if(abs(mu)<mu_cutoff)
				{
					patchx[k][q]=xcorner[q%(Nq/numring)]+rpatchvar[k][q]*cos(angle[k][q]);
					patchy[k][q]=ycorner[q%(Nq/numring)]+rpatchvar[k][q]*sin(angle[k][q]);
				}
				else
				{
					patchx[k][q]=rpatchvar[k][q]*cos(angle[k][q]);
					patchy[k][q]=rpatchvar[k][q]*sin(angle[k][q]);
				}

			}
		}
	//-------------------------------------------------------------------
}
*/

//-------------------------------------------------------------------
void init(){

//-------------------------------------------------------------------
//begin cutoff at maximum energy at gamma point
laminit=max(abs(energy(0,0,0)),abs(energy(0,0,1)));

dlam=dlaminit;
lam=laminit;
//-------------------------------------------------------------------

if (std::abs(mu) <= 1E-3)
	set_patch_points_K_const();
else if (std::abs(mu) < mu_cutoff && std::abs(mu) <= 1.)
	set_patch_points_K();
else
	set_patch_points_G();
//set_patch_points();
std::ofstream out("patching.txt");
for (int q=0;q<Nq;q++)
	for (int k=0;k<N;k++)
		//out <<  rmaxinit[k][0]*cos(angle[k][q]) << " " <<  rmaxinit[k][0]*sin(angle[k][q]) << std::endl;
		//out <<  xcorner[q] + rmaxinit[k][0]*cos(sangle[k][0][q]) << " " <<  ycorner[q] + rmaxinit[k][0]*sin(sangle[k][0][q]) << std::endl;
		out << patchx[k][q] << " " <<  patchy[k][q] << " " << rmaxinit[k][0] << " " << sangle[k][0][q] << std::endl;

//-------------------------------------------------------------------
//berechne 4ten Wellenvektor zu 3 gegebenen
impulserhaltung();
//-------------------------------------------------------------------


    //-------------------------------------------------------------------
    //set initial vertex functions in orbital space
    for (int k1=0;k1<N;k1++){
    for (int q1=0;q1<Nq;q1++){

	double kx1=patchx[k1][q1];
	double ky1=patchy[k1][q1];

    for (int o1=0;o1<No;o1++){

    	for (int k2=0;k2<N;k2++){
    	for (int q2=0;q2<Nq;q2++){

        double kx2=patchx[k2][q2];
        double ky2=patchy[k2][q2];

    	for (int o2=0;o2<No;o2++){

		    for (int k3=0;k3<N;k3++){
		    for (int q3=0;q3<Nq;q3++){

		 	double kx3=patchx[k3][q3];
			double ky3=patchy[k3][q3];

            //momenta that enter lattice (f, f2) and bond modulation (bmf, bmf2) factors (honeycomb)
			double qx=(kx3-kx1);
			double qy=(ky3-ky1);

            double qx2=(kx3-kx2);
            double qy2=(ky3-ky2);

                qx=qx*pi;
				qy=qy*pi;

                qx2=qx2*pi;
				qy2=qy2*pi;

				//nn-modulation factors
                complex<double> f=0;
                for (int s=0;s<3;s++){
                    f+= exp(I*(qx*deltax[s]+qy*deltay[s]));
                }


                complex<double> f2=0;
                for (int s=0;s<3;s++){
                    f2+= exp(I*(qx2*deltax[s]+qy2*deltay[s]));
                }


                //nnn-modulation factor
                complex<double> D=0;
                for (int s=0;s<6;s++){
                    D+= exp(I*(qx*Dx[s]+qy*Dy[s]));
                }


                complex<double> D2=0;
                for (int s=0;s<6;s++){
                    D2+= exp(I*(qx2*Dx[s]+qy2*Dy[s]));
                }

                    for (int o3=0;o3<No;o3++){
                    for (int o4=0;o4<No;o4++){

                            //initial condition in orbital representation

                            O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]=0;

                            // Onsite interaction (Hubbard-U) [forbidden for spinless fermions]

                            // nearest-neighbor-interaction Vnn

                            if(o1==0 && o2==1 && o3==0 && o4==1) O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]+= -Vnn*f;
                            if(o1==1 && o2==0 && o3==0 && o4==1) O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]+= +Vnn*f2;
                            if(o1==0 && o2==1 && o3==1 && o4==0) O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]+= +Vnn*conj(f2);
                            if(o1==1 && o2==0 && o3==1 && o4==0) O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]+= -Vnn*conj(f);

                            // next-nearest-neighbor-interaction Vnnn

					        if(o1==0 && o2==0 && o3==0 && o4==0) O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]+= -Vnnn*(D - D2);
                            if(o1==1 && o2==1 && o3==1 && o4==1) O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]+= -Vnnn*(D - D2);



					    }
				    }


			    }
			    }

		    }
	    	}
			}

    }
    }
    }


    trafocalc();
    trafozuband();
}
//-------------------------------------------------------------------
