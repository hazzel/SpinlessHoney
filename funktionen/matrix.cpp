//--------------------------------------------(define hopping matrix)
void matrix(double kx,double ky,double (&hopmat)[2*No*No]){

    

//--------------------------------------------------(matrix elements)
double d1r;
double d1im;
    
d1r = t*(2*cos(sqrt(3.0)/2.0*pi*kx)*cos(pi*ky/2.0)+cos(-pi*ky)); //NN hopping, real part
d1im =t*(2*cos(sqrt(3.0)/2.0*pi*kx)*sin(pi*ky/2.0)+sin(-pi*ky)); //NN hopping, imaginary part
//-------------------------------------------------------------------

    


//-----------------------------------------------(2x2 complex matrix)
double hoppingmatrix[2*No*No] = {/*Re 00*/ 0.0, /*Im 00*/ 0.0, /*Re 01*/ -d1r, /*Im 01*/ +d1im, /*Re 10*/ -d1r, /*Im 10*/ -d1im, /* Re 11*/ 0.0, /*Im 11*/ 0.0};
//-------------------------------------------------------------------




//---------------------------------------------(pass matrix elements)
int i;    
    
    
for(i=0; i< 2*No*No; i++){
hopmat[i] = hoppingmatrix[i];
}
//-------------------------------------------------------------------
  


}
//-------------------------------------------------------------------
