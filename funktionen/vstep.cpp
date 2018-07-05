void dvstep()
{
	summaxold=summax;
	summax=-100;

	#pragma omp parallel for collapse(10)
	for (int k1=0; k1<N; k1++)
	for (int q1=0; q1<Nq; q1++)
	for (int b1=0; b1<Nb; b1++)
		for (int k2=0; k2<N; k2++)
		for (int q2=0; q2<Nq; q2++)
		for (int b2=0; b2<Nb; b2++)
			for (int k3=0; k3<N; k3++)
			for (int q3=0; q3<Nq; q3++)
			for (int b3=0; b3<Nb; b3++)
				for (int b4=0; b4<Nb; b4++)
				{
					complex<double> sum=0;

					for (int ks=0;ks<N;ks++)
					for (int qs=0;qs<Nq;qs++)
					for (int bs=0;bs<Nb;bs++)
					for (int bs2=0;bs2<Nb;bs2++)
					{
						//pp-diagramm:
						int kpp=impuls[k1][q1][k2][q2][ks][qs][0];
						int qpp=impuls[k1][q1][k2][q2][ks][qs][1];

						
						sum-= 0.5*V[k1][q1][b1][k2][q2][b2][kpp][qpp][bs2][bs]*V[ks][qs][bs][kpp][qpp][bs2][k3][q3][b3][b4]*Lpp[ks][qs][k1][q1][k2][q2][bs][bs2];
						sum-= 0.5*V[k1][q1][b1][k2][q2][b2][ks][qs][bs][bs2]*V[kpp][qpp][bs2][ks][qs][bs][k3][q3][b3][b4]*Lpp[ks][qs][k1][q1][k2][q2][bs][bs2];
						
						//if (lam > 0.15)
						{
						//phcross-diagramm:
						int kphc1=impuls[ks][qs][k2][q2][k3][q3][0];
						int qphc1=impuls[ks][qs][k2][q2][k3][q3][1];
						int kphc2=impuls[ks][qs][k3][q3][k2][q2][0];
						int qphc2=impuls[ks][qs][k3][q3][k2][q2][1];
						
						sum-= -V[ks][qs][bs][k2][q2][b2][k3][q3][b3][bs2]*V[k1][q1][b1][kphc1][qphc1][bs2][ks][qs][bs][b4]*Lph[ks][qs][k2][q2][k3][q3][bs][bs2];
						sum-= -V[kphc2][qphc2][bs2][k2][q2][b2][k3][q3][b3][bs]*V[k1][q1][b1][ks][qs][bs][kphc2][qphc2][bs2][b4]*Lph[ks][qs][k3][q3][k2][q2][bs][bs2];
						
						//phdirect-diagramm:
						int kphd1=impuls[ks][qs][k3][q3][k1][q1][0];
						int qphd1=impuls[ks][qs][k3][q3][k1][q1][1];
						int kphd2=impuls[ks][qs][k1][q1][k3][q3][0];
						int qphd2=impuls[ks][qs][k1][q1][k3][q3][1];

						sum-= V[ks][qs][bs][k2][q2][b2][kphd1][qphd1][bs2][b4]*V[k1][q1][b1][kphd1][qphd1][bs2][k3][q3][b3][bs]*Lph[ks][qs][k3][q3][k1][q1][bs][bs2];
						sum-= V[kphd2][qphd2][bs2][k2][q2][b2][ks][qs][bs][b4]*V[k1][q1][b1][ks][qs][bs][k3][q3][b3][bs2]*Lph[ks][qs][k1][q1][k3][q3][bs][bs2];
						}
					}

					double propnorm = 1.0/numring;
					sum=sum* 2.0*pi/(3.0*N) * dlam * fnorm * propnorm; /* factor of 1/numring to account for numring */
					dV[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4]=sum;			//dlambda * diagrams
					if(abs(sum)>summax) summax=abs(sum);
				}
}

//calculates RG step
void vstep()
{
	vmax=-1000;

	for (int k1=0; k1<N; k1++)
	for (int q1=0; q1<Nq; q1++)
	for (int b1=0; b1<Nb; b1++)
		for (int k2=0; k2<N; k2++)
		for (int q2=0; q2<Nq; q2++)
		for (int b2=0; b2<Nb; b2++)
			for (int k3=0; k3<N; k3++)
			for (int q3=0; q3<Nq; q3++)
			for (int b3=0; b3<Nb; b3++)
				for (int b4=0; b4<Nb; b4++)
				{
					V[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4]+=dV[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4]; 		//increase V by dV
					if (abs(V[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4])>vmax)
						vmax=abs(V[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4]);
				}
}
