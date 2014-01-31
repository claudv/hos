//
//  model_2d.c
//  euler
//
//  Created by Claudio Viotti on 6/24/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//


#include "model_2d.h"
#include "hdf5_routines_2d.h"

void rhs_test(fftw_complex* hrhs, fftw_complex* hu){

    int index, index_shift;
    double kx, ky;
    
    for (int i=0; i<Nx/2; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            kx = i*Kx0;
            ky = j*Ky0;
            index = (Ny/2+1)*i + j;
            hrhs[index] = I*kx*hu[index];

            index_shift = index + Nx*(Ny/2+1);
            hrhs[index_shift] = I*ky*hu[index_shift];
            
        }
        
    }


    for (int i=Nx/2; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            kx = -(Nx-i)*Kx0;
            ky = j*Ky0;
            index = (Ny/2+1)*i + j;
            hrhs[index] = I*kx*hu[index];

            index_shift = index + Nx*(Ny/2+1);
            hrhs[index_shift] = I*ky*hu[index_shift];
            
        }
        
    }


}




void rhs_hos_setup(){

    eta_x = malloc(sizeof(double)*Nx*Ny);
    phi_x = malloc(sizeof(double)*Nx*Ny);
    
    temp1 = malloc(sizeof(double)*Nx*Ny);
    temp2 = malloc(sizeof(double)*Nx*Ny);
    temp1_large = malloc(sizeof(double)*Nxl*Nyl);
    temp2_large = malloc(sizeof(double)*Nxl*Nyl);
    htemp2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    htemp1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    htemp1_large = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nxl*(Nyl/2+1));
    htemp2_large = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nxl*(Nyl/2+1));

    hetan = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1)*(NLevs+1));
    hphin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1)*(NLevs+1));
    
    hwn   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1)*(NLevs+1));
    hwM   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    hwM2  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    hw2M  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    hw2M2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));

    
    Coeff = malloc(sizeof(double)*(NLevs+1));
    Coeff[0] = 1.0;
    Coeff[1] = 1.0;
    
    for (int n=2; n<=NLevs; n++) {
    
        Coeff[n] = Coeff[n-1]*n;
    
    }
    
}


/* HOS scheme (West et al. 1987) */
void rhs_hos(fftw_complex* hrhs, fftw_complex* hu, double t){

    int index, index_shift;
    
    /* eta_x*phi_x + eta_y*phi_y */
    Dx(hu, htemp1);
    Dx(&hu[Nx*(Ny/2+1)], htemp2);
    Mult(htemp1, htemp2, hrhs);
    
    Dy(hu, htemp1);
    Dy(&hu[Nx*(Ny/2+1)], htemp2);
    Mult(htemp1, htemp2, htemp1);
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
    
            index = (Ny/2+1)*i + j;
            hrhs[index] = - hrhs[index] - htemp1[index];
    
        }
        
    }
    
  
    /* Vertical velocity */
    Zvel(hu, hwM, hwM2, hw2M, hw2M2, t);
    
    
    /* eta_x*eta_x + eta_y*eta_y */
    Dx(hu, htemp1);
    Mult(htemp1, htemp1, htemp1);  /* eta_x*eta_x --> htemp1 */

    Dy(hu, htemp2);
    Mult(htemp2, htemp2, htemp2);  /* eta_y*eta_y --> htemp2 */
    
    Sum(htemp1, htemp2, htemp1);   /* eta_x*eta_x + eta_y*eta_y --> htemp1  */
    Mult(htemp1, hwM2, htemp2);    /* W^(M-2)*(eta_x*eta_x + eta_y*eta_y) --> htemp2  */
    Sum(hwM, htemp2, htemp2);      /* W^(M-2)*(eta_x*eta_x + eta_y*eta_y) + W^(M) --> htemp2  */
    
    Sum(hrhs, htemp2, hrhs);




    /* RHS for phi part */
    
    Mult(htemp1, hw2M2, htemp2);   /* W2^(M-2)*(eta_x*eta_x + eta_y*eta_y) --> htemp2  */
    Sum(hw2M, htemp2, htemp2);     /* W2^(M-2)*(eta_x*eta_x + eta_y*eta_y) + W2^(M) --> htemp2  */
    
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            index = (Ny/2+1)*i + j;
            index_shift = index + Nx*(Ny/2+1);
            hrhs[index_shift] = -g*hu[index];
            hrhs[index_shift] = hrhs[index_shift] + 0.5*htemp2[index];

        }
        
    }
    
    
    /* -0.5*(phi_x*phi_x + phi_y*phi_y) */
    Dx(&hu[Nx*(Ny/2+1)], htemp1);
    Mult(htemp1, htemp1, htemp1);
  
    Dy(&hu[Nx*(Ny/2+1)], htemp2);
    Mult(htemp2, htemp2, htemp2);
    
    Sum(htemp1, htemp2, htemp1);   /* phi_x*phi_x + phi_y*phi_y --> htemp1  */
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            index = (Ny/2+1)*i + j;
            index_shift = index + Nx*(Ny/2+1);
            hrhs[index_shift] = hrhs[index_shift] - 0.5*htemp1[index];

        }
        
    }

    /* Dealiasing */
    
    int mx, my;
    
    mx = floor( 0.5*Nx/(1 + 0.5*NLevs) );
    my = floor( 0.5*Ny/(1 + 0.5*NLevs) );
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            if ( (i>=mx && i<(Nx-mx+1)) || (j>=my) ) {
            
                index = (Ny/2+1)*i + j;
                index_shift = index + Nx*(Ny/2+1);
                hrhs[index] = 0.0;
                hrhs[index_shift] = 0.0;
                        
            }
            
        }
        
    }



}



/* Compute vertical velocity on the free surface. */
/* HOS scheme (West et al. 1987) */
void Zvel(fftw_complex* hu, fftw_complex* hwM, fftw_complex* hwM2, fftw_complex* hw2M, fftw_complex* hw2M2, double t){

    int index, Nindex, Mindex;

//    double kx, ky, kz;
//
//
//    for (int i=0; i<Nx/2; i++) {
//        
//        for (int j=0; j<Ny/2+1; j++) {
//            
//            kx = i*Kx0;
//            ky = j*Ky0;
//            kz = kx*kx + ky*ky;
//            kz = sqrt(kz);
//    
//            index = (Ny/2+1)*i + j;
//            index_shift = index + Nx*(Ny/2+1);
//            hw[index] = kz*hu[index_shift];
//
//        }
//        
//    }
//
//
//    for (int i=Nx/2; i<Nx; i++) {
//        
//        for (int j=0; j<Ny/2+1; j++) {
//            
//            kx = -(Nx-i)*Kx0;
//            ky = j*Ky0;
//            kz = kx*kx + ky*ky;
//            kz = sqrt(kz);
//    
//            index = (Ny/2+1)*i + j;
//            index_shift = index + Nx*(Ny/2+1);
//            hw[index] = kz*hu[index_shift];
//
//        }
//        
//    }


    for (int n=0; n<=NLevs; n++) {
    
        for (int i=0; i<Nx; i++) {
        
            for (int j=0; j<Ny/2+1; j++) {
    
                index = (Ny/2+1)*i + j;
                hwn[n*(Nx*(Ny/2+1)) + index] = 0.0;
                hphin[n*(Nx*(Ny/2+1)) + index] = 0.0;
             
            }
        
        }
    
    }
    
    
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
    
            index = (Ny/2+1)*i + j;
            hwM[index] = 0.0;
            hwM2[index] = 0.0;
            hw2M[index] = 0.0;
            hw2M2[index] = 0.0;
            hetan[index] = 0.0;
             
        }
        
    }

    hetan[0] = 1.0;




    Dz( &hu[Nx*(Ny/2+1)], hphin);
    
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
    
            index = (Ny/2+1)*i + j;
            hwn[index] = hphin[index];
    
        }
        
    }




    for (int n=1; n<=NLevs; n++) {
    
        Nindex = n*(Nx*(Ny/2+1));
        
        Mult(hu, &hetan[(n-1)*(Nx*(Ny/2+1))], &hetan[Nindex]);

        /* Phi_n */
        for (int m=0; m<=n-1; m++) {
    
            Mindex = m*(Nx*(Ny/2+1));
            
            Mult( &hphin[Mindex], &hetan[Nindex - Mindex], htemp1);
    
            for (int i=0; i<Nx; i++) {
        
                for (int j=0; j<Ny/2+1; j++) {
    
                    index = (Ny/2+1)*i + j;
                    hphin[Nindex + index] = hphin[Nindex + index] - htemp1[index]/Coeff[n-m];
             
                }
        
            }
    
        }
        
        /* W_n */
        for (int m=0; m<=n; m++) {
        
            Mindex = m*(Nx*(Ny/2+1));
        
            Dz( &hphin[Mindex], &hphin[Mindex]);
            Mult( &hphin[Mindex], &hetan[Nindex - Mindex], htemp1);
            
            for (int i=0; i<Nx; i++) {
        
                for (int j=0; j<Ny/2+1; j++) {
    
                    index = (Ny/2+1)*i + j;
                    hwn[Nindex + index] = hwn[Nindex + index] + htemp1[index]/Coeff[n-m];
             
                }
        
            }

    
        }


    
    }

    
//    ifft_2d( &hphin[2*(Nx*(Ny/2+1))], temp1, ifftp);
//    ifft_2d( &hphin[3*(Nx*(Ny/2+1))], temp2, ifftp);
//    
//    savefileid = open_file_2d(savefile_buff);
//    write_extra_2d(savefileid, temp1, temp2, 1);
//    status = close_file_2d(savefileid);
    
    
    //hwM2
    for (int n=0; n<=(NLevs-2); n++) {
    
        for (int i=0; i<Nx; i++) {
        
            for (int j=0; j<Ny/2+1; j++) {
    
                index = (Ny/2+1)*i + j;
                hwM[index] = hwM[index] + hwn[n*(Nx*(Ny/2+1)) + index];
                hwM2[index] = hwM2[index] + hwn[n*(Nx*(Ny/2+1)) + index];
             
            }
        
        }
    
    }
    
    //hwM
    for (int n=(NLevs-1); n<=NLevs; n++) {
        
        for (int i=0; i<Nx; i++) {
            
            for (int j=0; j<Ny/2+1; j++) {
                
                index = (Ny/2+1)*i + j;
                hwM[index] = hwM[index] + hwn[n*(Nx*(Ny/2+1)) + index];
                
            }
            
        }
        
    }

    
    
    
    // W2M2 = W2^(0) + W2^(1) + ... W2^(M-2)
    for (int n=0; n<=(NLevs-2); n++) {
        
        for (int i=0; i<Nx; i++) {
            
            for (int j=0; j<Ny/2+1; j++) {
                
                index = (Ny/2+1)*i + j;
                htemp2[index] = 0;
                
            }
            
        }

        // W2^(n) = W^(0)*W^(n) + W^(1)*W^(n-1) + ... W^(n)*W^(0)
        // W2^(n) stored in temp2
        for (int m=0; m<=n; m++) {
        
            Mult(&hwn[m*(Nx*(Ny/2+1))], &hwn[(n-m)*(Nx*(Ny/2+1))], htemp1);
            
            for (int i=0; i<Nx; i++) {
                
                for (int j=0; j<Ny/2+1; j++) {
                    
                    index = (Ny/2+1)*i + j;
                    htemp2[index] = htemp2[index] + htemp1[index];
                    
                }
                
            }
            
        }
        
        Sum(hw2M2, htemp2, hw2M2);
        Sum(hw2M, htemp2, hw2M);
        
        
    }
    
    
    // Loop over two more levels to finish assembling: W2M = W2^(0) + W2^(1) + ... W2^(M)
    for (int n=(NLevs-1); n<=NLevs; n++) {
                
        for (int i=0; i<Nx; i++) {
            
            for (int j=0; j<Ny/2+1; j++) {
                
                index = (Ny/2+1)*i + j;
                htemp2[index] = 0;
                
            }
            
        }

        // W2^(n) = W^(0)*W^(n) + W^(1)*W^(n-1) + ... W^(n)*W^(0)
        // W2^(n) stored in temp2
        for (int m=0; m<=n; m++) {
            
            Mult(&hwn[m*(Nx*(Ny/2+1))], &hwn[(n-m)*(Nx*(Ny/2+1))], htemp1);
            
            for (int i=0; i<Nx; i++) {
                
                for (int j=0; j<Ny/2+1; j++) {
                    
                    index = (Ny/2+1)*i + j;
                    htemp2[index] = htemp2[index] + htemp1[index];
                    
                }
                
            }
            
        }
        
        Sum(hw2M, htemp2, hw2M);
        
        
    }

    

}


