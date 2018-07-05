//-------------------------------------------------------------------
void impulserhaltung(){


	//Fixing three given momenta

     for (int k1=0;k1<N;k1++){
     for (int q1=0;q1<Nq;q1++){   

         double x1=patchx[k1][q1];
         double y1=patchy[k1][q1]; 

     for (int k2=0;k2<N;k2++){
     for (int q2=0;q2<Nq;q2++){

         double x2=patchx[k2][q2];
         double y2=patchy[k2][q2]; 

     for (int k3=0;k3<N;k3++){
     for (int q3=0;q3<Nq;q3++){

         double x3=patchx[k3][q3];
         double y3=patchy[k3][q3]; 

   //Sum of the free momenta ks

         double xs=x1+x2-x3;    
         double ys=y1+y2-y3;

   //Bringing ks into the first BZ

	double dist=1000;
	double xg=0;
	double yg=0;

	for (int i=-4;i<5;i++){
	for (int j=-4;j<5;j++){
	for (int k=-4;k<5;k++){
	
	double x0=i*G1x+j*G2x+k*G3x;
	double y0=i*G1y+j*G2y+k*G3y;
	
	double dist2=pow(xs-x0,2)+pow(ys-y0,2);
	if(dist2<dist){
            dist=dist2;
            xg=xs-x0;
	    yg=ys-y0;
        }   
	
	}
	}
	}
           
   //Find nearest patchpoint for k4
     
        dist=1000;
	int resk=-1;
	int resq=-1;
                                 
        for (int k4=0;k4<N;k4++){
        for (int q4=0;q4<Nq;q4++){
                                                
             double xx=patchx[k4][q4];
             double yy=patchy[k4][q4];
                                     
             double dist2=pow(xg-xx,2)+pow(yg-yy,2);
             if(dist2<dist){
                 dist=dist2;
                 resk=k4;
                 resq=q4;
            }    
       }
       }

   //Momentum conservaion function
                                                                      
       impuls[k1][q1][k2][q2][k3][q3][0]=resk;
       impuls[k1][q1][k2][q2][k3][q3][1]=resq;                                         
     }
     }
                
     }
     }
           
     }    
     }
}
//-------------------------------------------------------------------
