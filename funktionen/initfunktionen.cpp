// INITFUNCTIONS /////////////////////////////////////////////////////////////////////////////

int sign(double i){
	int a=1;
	if(i<0)a=-1;
	return a;
}

// UNITARY TRANSFORMATION MATRIX /////////////////////////////////////////////////////////////

complex<double> trafo2band(int i, int j,double x,double y){

    complex<double> res;
    complex<double> I (0.0,1.0);
    complex<double> dk= exp(I*(0.5*sqrt(3)*pi*x+0.5*pi*y))+exp(I*(-0.5*sqrt(3)*pi*x+0.5*pi*y))+exp(I*(-pi*y)) ;
    complex<double> dkc=conj( exp(I*(0.5*sqrt(3)*pi*x+0.5*pi*y))+exp(I*(-0.5*sqrt(3)*pi*x+0.5*pi*y))+exp(I*(-pi*y)) );
    double dk2=3+2*cos(sqrt(3)*pi*x)+2*cos(0.5*(sqrt(3)*pi*x-3*pi*y))+2*cos(0.5*(sqrt(3)*pi*x+3*pi*y));
    double absdk=sqrt(dk2);
    
    if(i==1){              
        if(j==0) res=-sqrt(0.5)*sqrt(dkc/absdk);
	if(j==1) res=sqrt(0.5)*sqrt(dk/absdk);
    }
      
    if(i==0){              

        if(j==0) res=sqrt(0.5)*sqrt(dkc/absdk);
	if(j==1) res=sqrt(0.5)*sqrt(dk/absdk);
    }

    return res;                         
}

// EVALUATION OF THE UNITARY TRANSFORMATION ON THE PATCHES ////////////////////////////////////

void trafocalc(){
    for (int k=0;k<N;k++){   
    for (int q=0;q<Nq;q++){ 
    for (int o=0;o<No;o++){ 
    for (int b=0;b<Nb;b++){

	double x=patchx[k][q];
	double y=patchy[k][q];
	
	patchtrafo[k][q][b][o]=trafo2band(b,o,x,y);
    }}}}
}


/*
void trafocalc(){		//i - band index, j - original field index

for (int k=0;k<N;k++){   
for (int q=0;q<Nq;q++){ 

	double kx=patchx[k][q];
	double ky=patchy[k][q];

	double hopmat[2*No*No];
	
	matrix(kx,ky,hopmat);   

    //allocate workspace
    gsl_matrix_complex_view m
    = gsl_matrix_complex_view_array (hopmat, 2, 2);
    
    gsl_vector *eval = gsl_vector_alloc (2);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (2, 2);
    
    gsl_eigen_hermv_workspace * w =
    gsl_eigen_hermv_alloc (2);

    //compute eigensystem
    gsl_eigen_hermv (&m.matrix, eval, evec, w);
    
    //sort eigensystem in ascending order
    gsl_eigen_hermv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

	//choose eigenvector
	for (int b=0;b<Nb;b++){ 
    
    gsl_vector_complex_view evec_b = gsl_matrix_complex_column (evec, b);

    //choose entry            
    for (int o=0;o<No;o++){
    
    //transformation of cb to co/co+ to cb+	
    gsl_complex patch_b_o = gsl_vector_complex_get(&evec_b.vector, o); 
    
	//transfer to array
    complex<double> trafo (GSL_REAL(patch_b_o),GSL_IMAG(patch_b_o));
	patchtrafo[k][q][b][o] = trafo; //Johann: yields u_ob
    
    }
	}	


	//free workspace
    gsl_eigen_hermv_free (w);
    gsl_vector_free(eval);
    gsl_matrix_complex_free(evec);
    
    }
	}                         
}
*/

// TRANSFORM VERTEX FROM ORBITAL REPRESENTATION TO BAND REPRESENTATION //////////////////////////

void trafozuband(){
    for (int k1=0;k1<N;k1++){   
    for (int q1=0;q1<Nq;q1++){  
    for (int b1=0;b1<Nb;b1++){  
                                                       
    for (int k2=0;k2<N;k2++){
    for (int q2=0;q2<Nq;q2++){           
    for (int b2=0;b2<Nb;b2++){
                                      
    for (int k3=0;k3<N;k3++){
    for (int q3=0;q3<Nq;q3++){ 
                  
    int k4=impuls[k1][q1][k2][q2][k3][q3][0];
    int q4=impuls[k1][q1][k2][q2][k3][q3][1];
                             
    for (int b3=0;b3<Nb;b3++){                                                                              
    for (int b4=0;b4<Nb;b4++){
                                
         V[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4]=0;
             			                                      
           for(int o1=0;o1<No;o1++){                                      
	   for(int o2=0;o2<No;o2++){                                             
	   for(int o3=0;o3<No;o3++){
	   for(int o4=0;o4<No;o4++){                                
                                
	   V[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4]+=O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]*conj(patchtrafo[k3][q3][b3][o3])*conj(patchtrafo[k4][q4][b4][o4])*patchtrafo[k2][q2][b2][o2]*patchtrafo[k1][q1][b1][o1];
	   }}}}                                                      
    }}
                                                 
    }}
                                 
    }}}
                  
    }}}
}

// TRANSFORM VERTEX FROM BAND REPRESENTATION TO ORBITAL REPRESENTATION //////////////////////////

void trafozuorbital(){
    for (int k1=0;k1<N;k1++){   
    for (int q1=0;q1<Nq;q1++){  
    for (int o1=0;o1<No;o1++){  
                                                       
    for (int k2=0;k2<N;k2++){
    for (int q2=0;q2<Nq;q2++){           
    for (int o2=0;o2<No;o2++){
                                      
    for (int k3=0;k3<N;k3++){
    for (int q3=0;q3<Nq;q3++){ 
                  
    int k4=impuls[k1][q1][k2][q2][k3][q3][0];
    int q4=impuls[k1][q1][k2][q2][k3][q3][1];
                             
    for (int o3=0;o3<No;o3++){                                                                              
    for (int o4=0;o4<No;o4++){
                                
         O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]=0;
             			                                      
           for(int b1=0;b1<Nb;b1++){                                      
	   for(int b2=0;b2<Nb;b2++){                                             
	   for(int b3=0;b3<Nb;b3++){
	   for(int b4=0;b4<Nb;b4++){                                
                                
	   O[k1][q1][o1][k2][q2][o2][k3][q3][o3][o4]+=V[k1][q1][b1][k2][q2][b2][k3][q3][b3][b4]*patchtrafo[k3][q3][b3][o3]*patchtrafo[k4][q4][b4][o4]*conj(patchtrafo[k2][q2][b2][o2])*conj(patchtrafo[k1][q1][b1][o1]);
	   }}}}                                                      
    }}                                                 
    
    }}
                                 
    }}}
                  
    }}}
}
