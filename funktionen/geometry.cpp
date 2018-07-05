//-------------------------------------------------------------------
void geometry(){


    
//-------------------------------------------------------------------
//volume of BZ, divided by pi^2 (scaling factor for r and dr)
    
     double a=8/(3*sqrt(3));
     fnorm=1.0/a;
//-------------------------------------------------------------------

    

//-------------------------------------------------------------------
//distance gamma to corner[0]
    
     double x0=2.0/(3.0*sqrt(3)); // 1/2 of the y coordinate of b2 divided by cos(30 degree)
     double y0=2.0/3.0;           // 1/2 of the y coordinate of b2
     
//coordinates of corners for Gamma = (0,0)
    
     xcorner[0]=-x0;		//down left
     ycorner[0]=-y0;

     xcorner[1]=x0;			//down right
     ycorner[1]=-y0;

     xcorner[2]=2.0*x0;		//right
     ycorner[2]=0.0;

     xcorner[3]=x0;			//up right
     ycorner[3]=y0;

     xcorner[4]=-x0;		//up left
     ycorner[4]=y0;

     xcorner[5]=-2.0*x0;	//left
     ycorner[5]=0.0;
//-------------------------------------------------------------------
    
    

//-------------------------------------------------------------------     
//distances in real lattice - nearest neighbor
    
		deltax[0]=sqrt(3.0)/2.0;			//up right
		deltay[0]=0.5;

		deltax[1]=-sqrt(3.0)/2.0;			//up left
		deltay[1]=0.5;

		deltax[2]=0;						//down
		deltay[2]=-1.0;
//-------------------------------------------------------------------
    
    
    
//-------------------------------------------------------------------
//distances in real lattice - next nearest neighbor
    
		Dx[0]=deltax[1]-deltax[2];			//down to up left
		Dy[0]=deltay[1]-deltay[2];

		Dx[1]=deltax[0]-deltax[2];			//down to up right
		Dy[1]=deltay[0]-deltay[2];

		Dx[2]=deltax[1]-deltax[0];			//up right to up left
		Dy[2]=deltay[1]-deltay[0];

		Dx[3]=deltax[0]-deltax[1];			//up left to up right
		Dy[3]=deltay[0]-deltay[1];

		Dx[4]=-deltax[0]+deltax[2];			//up right to down
		Dy[4]=-deltay[0]+deltay[2];

		Dx[5]=-deltax[1]+deltax[2];			//up left to down
		Dy[5]=-deltay[1]+deltay[2];
//-------------------------------------------------------------------
    
    
    
//-------------------------------------------------------------------
//reciprocal lattice basis vectors
	
		G1x=2.0/sqrt(3);			//down right = b1
		G1y=2.0/3.0;

		G2x=-G1x;					//down left
		G2y=G1y;

		G3x=0;						//up = b2
		G3y=4.0/3.0;
//-------------------------------------------------------------------
    
    
    
}
//-------------------------------------------------------------------
