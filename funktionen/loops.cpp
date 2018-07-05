//-------------------------------------------------------------------
//finds r, for which the dispersion reaches energy en
double bisection(double rmax,double rmin, double en, double phi, int nq, int nb){
   double cphi=cos(phi);
   double sphi=sin(phi);

	double x0, y0;
	if (std::abs(mu) < mu_cutoff && std::abs(mu) <= 1.)
	{
		x0 = xcorner[nq%(Nq/numring)];
		y0 = ycorner[nq%(Nq/numring)];
	}
	else
	{
		x0 = 0.;
		y0 = 0.;
	}
   double x1=x0+rmax*cphi;		//k value shifted away from corner by rmax
   double y1=y0+rmax*sphi;

   double e1=energy(x1,y1,nb)-en;		//band energy at x1 minus en

   double x2=x0;				//k value at corner
   double y2=y0;

   double e2=energy(x2,y2,nb)-en;  		//band energy at corner minus en

    //en is larger than energy at the boundaries
    if (e1<0 && e2<0){
		return -1;
    }

    //en is smaller than energy at the boundaries
    else if (e1>0 && e2>0){
		return -2;
    }

	//works only for increasing energyfunction -> exchange rmax <-> rmin for decreasing e (lower band)
    else {
        if(e2>e1)
			std::swap(rmin, rmax);

        double dr=rmax-rmin;
        double r=rmin;
        double emid=1;

        while(abs(dr)>=bisectionprec && emid!=0){
		dr=dr/2;
		double rmid=r+dr;

   		double xmid=x0+rmid*cphi;
   		double ymid=y0+rmid*sphi;

		emid=energy(xmid,ymid,nb)-en;

		if(emid<=0){r=rmid;}
	    }

	    return r;
    }
}
//-------------------------------------------------------------------




//-------------------------------------------------------------------
//find radial coordinate for given energy along Gamma-M path:
double bisectionGM(double rmax,double rmin, double en, int nb){


    //coordinates gamma-point
    double x1=0;
    double y1=0;

    double e1=energy(x1,y1,nb)-en;

    //coordinates M-point
    double x2=0;
    double y2=-2.0/3.0;

    double e2=energy(x2,y2,nb)-en;

    //en ist grš§er als die energien am rand des intervalls
    if (e1<0 && e2<0){
		 return -1;
    }

    //en ist kleiner als die energien am rand des intervalls
    else if (e1>0 && e2>0){
		return -2;
    }

    else {
        if(e2>e1){
            //fuer diese Dispersion, falls Energien im unteren Band
            double f=rmin;
            rmin=rmax;
            rmax=f;
        }

        double dr=rmax-rmin;
        double r=rmin;
        double emid=1;

        while(abs(dr)>=bisectionprec && emid!=0){
            dr=dr/2;
            double rmid=r+dr;

            double ymid=-2.0/3.0+rmid;

            emid=energy(0,ymid,nb)-en;

            if(emid<=0){r=rmid;}
	    }

	    return r;
    }
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
//provides boundaries for radial integration
void shellvec(){

     for(int k=0;k<N;k++){
     for(int s=0;s<Ns;s++){
     for(int q=0;q<Nq;q++){
     for(int b=0;b<Nb;b++){


         double ang=sangle[k][s][q];
                
         double es=-1.5*lam;
         rd[k][s][q][b][1]=bisection(rmaxinit[k][s],0,es,ang,q,b);			//r where dispersion = -1.5lambda
              
         es=-0.5*lam;
         rup[k][s][q][b][1]=bisection(rmaxinit[k][s],0,es,ang,q,b);			//r where dispersion = -0.5lambda
          
         es=0.5*lam;
         rd[k][s][q][b][0]=bisection(rmaxinit[k][s],0,es,ang,q,b);			//r where dispersion = 0.5lambda
                            
         es=1.5*lam;
         rup[k][s][q][b][0]=bisection(rmaxinit[k][s],0,es,ang,q,b);         //r where dispersion = 1.5lambda
                }
            }
        }
    }
}
//-------------------------------------------------------------------



//-------------------------------------------------------------------
//calculates radial average of the loopkernels
void radint(double& sumph,double& sumpp,double rmin,double dr,double x1,double y1,double x2,double y2,double cphi,double sphi,int qs,int bs,int bs2){
     for (int i=0;i<Nr+1;i++){

         double rs=rmin+dr*i;
			
			double xs, ys;
			if (std::abs(mu) < mu_cutoff && std::abs(mu) <= 1.)
			{
				xs = xcorner[qs%(Nq/numring)]+rs*cphi;
				ys = ycorner[qs%(Nq/numring)]+rs*sphi;
			}
			else
			{
				xs = rs*cphi;
				ys = rs*sphi;
			}

         //ph-part;
         double xi=xs+x1-x2;
         double yi=ys+y1-y2;

         //pp-part
         double xj=-xs+x1+x2;
         double yj=-ys+y1+y2;

         double e=energy(xs,ys,bs);
         double ei=energy(xi,yi,bs2);
         double ej=energy(xj,yj,bs2);

         //derivative of cutoff function
         double dcut=(fermi((abs(e)-lam-dlam/2),mf*(lam+dlam/2))-fermi((abs(e)-lam+dlam/2),mf*(lam-dlam/2)))/dlam;

         //theta cutoff function
         double cutph=1-fermi((abs(ei)-lam),mf*lam);
         double cutpp=1-fermi((abs(ej)-lam),mf*lam);

         double facph;
         double facpp;

         //ph
         if (abs(e-ei)<deltaEmin){
			double en=(e+ei)/2;
         	facph=dfermi(en,T);						//for small energies -> derivative instead of difference
         }
	 	 else{
			facph=(fermi(e,T)-fermi(ei,T))/(e-ei);
         }

	 	//pp
	 	if(abs(e+ej)<deltaEmin){
        	double en=(e-ej)/2;
        	facpp=-dfermi(en,T);
         }
        else{
            facpp=(1-fermi(e,T)-fermi(ej,T))/(e+ej);
         }

	 sumph+=facph*abs(dr)*cutph*dcut*rs;//dr=delta_r/Nr
	 sumpp+=facpp*abs(dr)*cutpp*dcut*rs;
    }
}
//-------------------------------------------------------------------



//-------------------------------------------------------------------
//calculates loopkernel L = GS + SG
void loopcalc()
{
	loophmax=-1000;
	looppmax=-1000;

	#pragma omp parallel for collapse(10)
	for (int k1=0;k1<N;k1++)
	for (int q1=0;q1<Nq;q1++)
	for (int b1=0;b1<Nb;b1++)
		for (int k2=0;k2<N;k2++)
		for (int q2=0;q2<Nq;q2++)
		for (int b2=0;b2<Nb;b2++)
			for (int ks=0;ks<N;ks++)
			for (int qs=0;qs<Nq;qs++)
			for (int bs=0;bs<Nb;bs++)
			for(int bs2=0;bs2<Nb;bs2++)
			{
				double x1=patchx[k1][q1];
				double y1=patchy[k1][q1];
				double x2=patchx[k2][q2];
				double y2=patchy[k2][q2];
				double sumph=0;
				double sumpp=0;

				for (int s=0;s<Ns;s++)
				{
					double cphi=cos(sangle[ks][s][qs]);
					double sphi=sin(sangle[ks][s][qs]);

					for (int shell=0;shell<2;shell++)
					{
						//bounds of radial integration: +- 0.5 lambda
						double rmax = rup[ks][s][qs][bs][shell];
						double rmin = rd[ks][s][qs][bs][shell];

						int tolarge=-1;
						int tosmall=-2;

						//exchange rmax and rmin in the lower band to ensure increasing energy function
						if (std::abs(mu) < mu_cutoff && std::abs(mu) <= 1.)
						{
							if (bs == 0)
							{
								std::swap(rmin, rmax);
								std::swap(tolarge, tosmall);
							}
						}
						else
						{
							if (bs == 1)
							{
								std::swap(rmin, rmax);
								std::swap(tolarge, tosmall);
							}
						}

						if (rmax != tosmall && rmin != tolarge)
						{
							if(rmax==tolarge) rmax=rmaxinit[ks][s];
							if(rmin==tosmall) rmin=0;

							double dr=(rmax-rmin)/Nr;
							radint(sumph,sumpp,rmin,dr,x1,y1,x2,y2,cphi,sphi,qs,bs,bs2);
						}
					}
				}
				Lph[ks][qs][k1][q1][k2][q2][bs][bs2]=sumph/Ns;
				Lpp[ks][qs][k1][q1][k2][q2][bs][bs2]=sumpp/Ns;

				if(abs(Lph[ks][qs][k1][q1][k2][q2][bs][bs2])>loophmax)
					loophmax=abs(Lph[ks][qs][k1][q1][k2][q2][bs][bs2]);
				if(abs(Lpp[ks][qs][k1][q1][k2][q2][bs][bs2])>looppmax)
					looppmax=abs(Lpp[ks][qs][k1][q1][k2][q2][bs][bs2]);
			}
}
//-------------------------------------------------------------------
