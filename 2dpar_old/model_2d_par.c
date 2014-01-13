//
//  model_2d.c
//  euler
//
//  Created by Claudio Viotti on 6/24/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//


#include "model_2d_par.h"

void rhs_test(fftw_complex* hrhs, fftw_complex* hu){

    // Linear advection
    Dx(hu, hrhs);
    Dy(&hu[alloc_local], &hrhs[alloc_local]);
    
}




void rhs_hos_setup(){

    int n;
    
    eta_x = fftw_alloc_real(2 * alloc_local);
    phi_x = fftw_alloc_real(2 * alloc_local);
    
    temp1 = fftw_alloc_real(2 * alloc_local);
    temp2 = fftw_alloc_real(2 * alloc_local);
    htemp1 = fftw_alloc_complex(alloc_local);
    htemp2 = fftw_alloc_complex(alloc_local);
    
    hetan = fftw_alloc_complex((NLevs+1)*alloc_local);
    hphin = fftw_alloc_complex((NLevs+1)*alloc_local);
    
    hwn = fftw_alloc_complex((NLevs+1)*alloc_local);
    hwM = fftw_alloc_complex(alloc_local);
    hwM2 = fftw_alloc_complex(alloc_local);
    hw2M = fftw_alloc_complex(alloc_local);
    hw2M2 = fftw_alloc_complex(alloc_local);
    
    Coeff = malloc(sizeof(double)*(NLevs+1));
    Coeff[0] = 1.0;
    Coeff[1] = 1.0;
    
    for (n=2; n<=NLevs; n++) {
    
        Coeff[n] = Coeff[n-1]*n;
    
    }
    
}


/* HOS scheme (West et al. 1987) */
void rhs_hos(fftw_complex* hrhs, fftw_complex* hu, double t){

    ptrdiff_t index, index_shift;
    
    /* eta_x*phi_x + eta_y*phi_y */
    Dx(hu, htemp1);
    Dx(&hu[alloc_local], htemp2);
    Mult(htemp1, htemp2, hrhs);
    
    Dy(hu, htemp1);
    Dy(&hu[alloc_local], htemp2);
    Mult(htemp1, htemp2, htemp1);
    
//    for (i=0; i<local_Nx; i++) {
//        
//        for (j=0; j<Ny/2+1; j++) {
//    
//            index = (Ny/2+1)*i + j;
//            hrhs[index] = - hrhs[index] - htemp1[index];
//    
//        }
//        
//    }
    
    for (j=0; j<local_Nyhpo; j++) {
    
        for (i=0; i<Nx; i++) {
        
            index = Nx*j + i;
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
    
    
//    for (i=0; i<local_Nx; i++) {
//        
//        for (j=0; j<Ny/2+1; j++) {
//            
//            index = (Ny/2+1)*i + j;
//            index_shift = index + alloc_local;
//            hrhs[index_shift] = -g*hu[index];
//            hrhs[index_shift] = hrhs[index_shift] + 0.5*htemp2[index];
//
//        }
//        
//    }
    
    
    for (j=0; j<local_Nyhpo; j++) {
    
        for (i=0; i<Nx; i++) {
        
            index = Nx*j + i;
            index_shift = index + alloc_local;
            hrhs[index_shift] = -g*hu[index];
            hrhs[index_shift] = hrhs[index_shift] + 0.5*htemp2[index];
            
        }
        
    }    
    
    /* -0.5*(phi_x*phi_x + phi_y*phi_y) */
    Dx(&hu[alloc_local], htemp1);
    Mult(htemp1, htemp1, htemp1);
  
    Dy(&hu[alloc_local], htemp2);
    Mult(htemp2, htemp2, htemp2);
    
    Sum(htemp1, htemp2, htemp1);   /* phi_x*phi_x + phi_y*phi_y --> htemp1  */
    
//    for (i=0; i<local_Nx; i++) {
//        
//        for (j=0; j<Ny/2+1; j++) {
//            
//            index = (Ny/2+1)*i + j;
//            index_shift = index + alloc_local;
//            hrhs[index_shift] = hrhs[index_shift] - 0.5*htemp1[index];
//
//        }
//        
//    }
    
    for (j=0; j<local_Nyhpo; j++) {
    
        for (i=0; i<Nx; i++) {
        
            index = Nx*j + i;
            index_shift = index + alloc_local;
            hrhs[index_shift] = hrhs[index_shift] - 0.5*htemp1[index];

        }
        
    }        

    /* Dealiasing */
    Dealias(hrhs);
    Dealias(&hrhs[alloc_local]);
    
//    int mx, my;
//    
//    mx = floor( 0.5*Nx/(1 + 0.5*NLevs) );
//    my = floor( 0.5*Ny/(1 + 0.5*NLevs) );
//    
//    for (i=0; i<local_Nx; i++) {
//        
//        for (j=0; j<Ny/2+1; j++) {
//            
//            if ( ( (local_0_start + i)>=mx && (local_0_start + i)<(Nx-mx)) || (j>=my) ) {
//            
//                index = (Ny/2+1)*i + j;
//                index_shift = index + alloc_local;
//                hrhs[index] = 0.0;
//                hrhs[index_shift] = 0.0;
//                        
//            }
//            
//        }
//        
//    }



}



/* Compute vertical velocity on the free surface. */
/* HOS scheme (West et al. 1987) */
void Zvel(fftw_complex* hu, fftw_complex* hwM, fftw_complex* hwM2, fftw_complex* hw2M, fftw_complex* hw2M2, double t){

    ptrdiff_t index, Nindex, Mindex;
    ptrdiff_t n, m;
    

    for (n=0; n<=NLevs; n++) {
    
//        for (i=0; i<local_Nx; i++) {
//        
//            for (j=0; j<Ny/2+1; j++) {
//    
//                index = (Ny/2+1)*i + j;
//                hwn[n*alloc_local + index] = 0.0;
//                hphin[n*alloc_local + index] = 0.0;
//             
//            }
//
//        }
        
        
//        for (j=0; j<local_Nyhpo; j++) {
//    
//            for (i=0; i<Nx; i++) {
//        
//                index = Nx*j + i;
//                hwn[n*alloc_local + index] = 0.0;
//                hphin[n*alloc_local + index] = 0.0;
//            
//            }
//            
//        
//        }
        
        
        for (i=0; i<(Nx*local_Nyhpo); i++) {
    
                hwn[n*alloc_local + i] = 0.0;
                hphin[n*alloc_local + i] = 0.0;
            
        }
        

    
    }
    
    
    
//    for (i=0; i<local_Nx; i++) {
//        
//        for (j=0; j<Ny/2+1; j++) {
//    
//            index = (Ny/2+1)*i + j;
//            hwM[index] = 0.0;
//            hwM2[index] = 0.0;
//            hw2M[index] = 0.0;
//            hw2M2[index] = 0.0;
//            hetan[index] = 0.0;
//             
//        }
//        
//    }
    
    for (i=0; i<(Nx*local_Nyhpo); i++) {
        
        hwM[i] = 0.0;
        hwM2[i] = 0.0;
        hw2M[i] = 0.0;
        hw2M2[i] = 0.0;
        hetan[i] = 0.0;
        
    }


    if (mpi_rank==0) {
        hetan[0] = 1.0;
    }
    




    Dz( &hu[alloc_local], hphin);
    
    
//    for (i=0; i<local_Nx; i++) {
//        
//        for (j=0; j<Ny/2+1; j++) {
//    
//            index = (Ny/2+1)*i + j;
//            hwn[index] = hphin[index];
//    
//        }
//        
//    }

    for (i=0; i<(Nx*local_Nyhpo); i++) {

            hwn[i] = hphin[i];
        
    }



    for (n=1; n<=NLevs; n++) {
    
        Nindex = n*alloc_local;
        
        Mult(hu, &hetan[(n-1)*alloc_local], &hetan[Nindex]);

        /* Phi_n */
        for (m=0; m<=n-1; m++) {
    
            Mindex = m*alloc_local;
            
            Mult( &hphin[Mindex], &hetan[Nindex - Mindex], htemp1);
    
//            for (i=0; i<local_Nx; i++) {
//        
//                for (j=0; j<Ny/2+1; j++) {
//    
//                    index = (Ny/2+1)*i + j;
//                    hphin[Nindex + index] = hphin[Nindex + index] - htemp1[index]/Coeff[n-m];
//             
//                }
//        
//            }

                
            for (i=0; i<(Nx*local_Nyhpo); i++) {
                    
                hphin[Nindex + i] = hphin[Nindex + i] - htemp1[i]/Coeff[n-m];
                    
            }
            
    
        }
        
        
        /* W_n */
        for (m=0; m<=n; m++) {
        
            Mindex = m*alloc_local;
        
            Dz( &hphin[Mindex], &hphin[Mindex]);
            Mult( &hphin[Mindex], &hetan[Nindex - Mindex], htemp1);
            
            
//            for (i=0; i<local_Nx; i++) {
//        
//                for (j=0; j<Ny/2+1; j++) {
//    
//                    index = (Ny/2+1)*i + j;
//                    hwn[Nindex + index] = hwn[Nindex + index] + htemp1[index]/Coeff[n-m];
//             
//                }
//        
//            }

            
            for (i=0; i<(Nx*local_Nyhpo); i++) {
                    
                hwn[Nindex + i] = hwn[Nindex + i] + htemp1[i]/Coeff[n-m];
                    
            }
                

    
        }

    
    
    }

   
    
    //hwM2
    for (n=0; n<=(NLevs-2); n++) {
    
//        for (i=0; i<local_Nx; i++) {
//        
//            for (j=0; j<Ny/2+1; j++) {
//    
//                index = (Ny/2+1)*i + j;
//                hwM[index] = hwM[index] + hwn[n*alloc_local + index];
//                hwM2[index] = hwM2[index] + hwn[n*alloc_local + index];
//             
//            }
//        
//        }
        
        for (i=0; i<(Nx*local_Nyhpo); i++) {

            hwM[i] = hwM[i] + hwn[n*alloc_local + i];
            hwM2[i] = hwM2[i] + hwn[n*alloc_local + i];
                
        }
    
    }
    
    //hwM
    for (n=(NLevs-1); n<=NLevs; n++) {
        
//        for (i=0; i<local_Nx; i++) {
//            
//            for (j=0; j<Ny/2+1; j++) {
//                
//                index = (Ny/2+1)*i + j;
//                hwM[index] = hwM[index] + hwn[n*alloc_local + index];
//                
//            }
//            
//        }
        
        for (i=0; i<(Nx*local_Nyhpo); i++) {
        
            hwM[i] = hwM[i] + hwn[n*alloc_local + i];
                
        }

    }


       
    // W2M2 = W2^(0) + W2^(1) + ... W2^(M-2)
    for (n=0; n<=(NLevs-3); n++) {
        
//        for (i=0; i<local_Nx; i++) {
//            
//            for (j=0; j<Ny/2+1; j++) {
//                
//                index = (Ny/2+1)*i + j;
//                htemp2[index] = 0;
//                
//            }
//            
//        }
                
        for (i=0; i<(Nx*local_Nyhpo); i++) {
            
            htemp2[i] = 0.0;
            
        }


        // W2^(n) = W^(0)*W^(n) + W^(1)*W^(n-1) + ... W^(n)*W^(0)
        // W2^(n) stored in temp2
        for (m=0; m<=n; m++) {
        
            Mult(&hwn[m*alloc_local], &hwn[(n-m)*alloc_local], htemp1);
            
//            for (i=0; i<local_Nx; i++) {
//                
//                for (j=0; j<Ny/2+1; j++) {
//                    
//                    index = (Ny/2+1)*i + j;
//                    htemp2[index] = htemp2[index] + htemp1[index];
//                    
//                }
//                
//            }

            for (i=0; i<(Nx*local_Nyhpo); i++) {
                
                htemp2[i] = htemp2[i] + htemp1[i];
                
            }
            
            
            
        }
        
        Sum(hw2M2, htemp2, hw2M2);
        Sum(hw2M, htemp2, hw2M);
        
        
    }
    
    
    // Loop over two more levels to finish assembling: W2M = W2^(0) + W2^(1) + ... W2^(M)
    for (n=(NLevs-2); n<=(NLevs-1); n++) {
                
//        for (i=0; i<local_Nx; i++) {
//            
//            for (j=0; j<Ny/2+1; j++) {
//                
//                index = (Ny/2+1)*i + j;
//                htemp2[index] = 0;
//                
//            }
//            
//        }

        for (i=0; i<(Nx*local_Nyhpo); i++) {
                
                htemp2[i] = 0.0;
                
        }


        // W2^(n) = W^(0)*W^(n) + W^(1)*W^(n-1) + ... W^(n)*W^(0)
        // W2^(n) stored in temp2
        for (m=0; m<=n; m++) {
            
            Mult(&hwn[m*alloc_local], &hwn[(n-m)*alloc_local], htemp1);
            
//            for (i=0; i<local_Nx; i++) {
//                
//                for (j=0; j<Ny/2+1; j++) {
//                    
//                    index = (Ny/2+1)*i + j;
//                    htemp2[index] = htemp2[index] + htemp1[index];
//                    
//                }
//                
//            }
            
            for (i=0; i<(Nx*local_Nyhpo); i++) {
                
                htemp2[i] = htemp2[i] + htemp1[i];
                
            }

            
            
        }
        
        Sum(hw2M, htemp2, hw2M);
        
        
    }

    
}


