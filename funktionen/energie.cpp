//Energy from diagonalizing the Hopping matrix


double energy(double x, double y, int nb){
    double e;
    
    double cx=2*cos(sqrt(3)*pi*x);                  // kx = pi*x and ky = pi*y
    double cxminy=2*cos(0.5*sqrt(3)*pi*x-1.5*pi*y);
    double cxplsy=2*cos(0.5*sqrt(3)*pi*x+1.5*pi*y);
    
    double t2=pow(t,2);
    double enarg=4*t2*(3+cx+cxminy+cxplsy);
    
    if(enarg<0){
        enarg=0;
	   }
    
    if(nb==0){
        e=-0.5*sqrt(enarg) -mu;
    }
    
    if(nb==1){
        e=+0.5*sqrt(enarg) -mu;
    }
    
    return e;
}


/*
double energy(double kx, double ky, int nb){
       
    double e;

	double hopmat[2*No*No];
	
	matrix(kx,ky,hopmat);
	
    //allocate workspace
    gsl_matrix_complex_view m
    = gsl_matrix_complex_view_array (hopmat, 2, 2);
    
    gsl_vector *eval = gsl_vector_alloc (2);
    
    gsl_eigen_herm_workspace * w =
    gsl_eigen_herm_alloc (2);

    //compute eigenvalues and store them in vector eval      
    gsl_eigen_herm (&m.matrix, eval, w);
    
    //sort eigenvalues in ascending order
    gsl_sort_vector(eval);
    
    //pick energy for corresponding band
    e = gsl_vector_get (eval, nb)-mu;
	 //std::cout << "nb = " << nb << ", " << "e = " << e << std::endl;
    
    //free workspace
    gsl_eigen_herm_free (w);
    gsl_vector_free(eval);
    
    //return energy
    return e;              
}
*/


//fermi distribution
double fermi (double e,double T){		
	
	if(e/T>40)return 0;					//for e>>T = 0
	else if(e/T<-40)return 1;			//for e<<T = 1
	else{
		double f=1/(exp(e/T)+1); 
		return f;
	}
}

//derivative of fermi distribution
double dfermi(double e,double T){		
	
	double df=-(1/T)*fermi(e,T)*fermi(-e,T);
	return df;
}   	
