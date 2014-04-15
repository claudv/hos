//
//  model_2d.c
//  euler
//
//  Created by Claudio Viotti on 6/24/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//


#include "model_2d_par.h"


//static double*  eta_x;
//static double*  phi_x;
static fftw_complex*   hwn;
static fftw_complex*   hwM;
static fftw_complex*   hwM2;
static fftw_complex*   hw2M;
static fftw_complex*   hw2M2;
static double*         Coeff;

static fftw_complex*   hwindpress;

void rhs_test(fftw_complex* hrhs, fftw_complex* hu){

    // Linear advection
    Dx(hu, hrhs);
    Dy(&hu[alloc_local], &hrhs[alloc_local]);
    
}


void rhs_hos_setup(){

    int n;
    
//    eta_x = fftw_alloc_real(2 * alloc_local);
//    phi_x = fftw_alloc_real(2 * alloc_local);
    
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
    
    hwindpress = fftw_alloc_complex(alloc_local);
    
}


/* HOS scheme (West et al. 1987) */
void rhs_hos(fftw_complex* hrhs, fftw_complex* hu, double t){

    ptrdiff_t index, index_shift;
  
    
    /*------------------*/
    /* RHS for phi part */
    
    /* eta_x*phi_x + eta_y*phi_y */
    Dx(hu, htemp1);
    Dx(&hu[alloc_local], htemp2);
    Mult(htemp1, htemp2, hrhs);
    
    Dy(hu, htemp1);
    Dy(&hu[alloc_local], htemp2);
    Mult(htemp1, htemp2, htemp1);
        
    for (i=0; i<local_N; i++) {
        
        hrhs[i] = - hrhs[i] - htemp1[i];
            
    }

    /* Vertical velocity */
    Zvel(hu, hwM, hwM2, hw2M, hw2M2, t);
    
    Dx(hu, htemp1);                /* eta_x --> htemp1 */
    Mult(htemp1, htemp1, htemp1);  /* eta_x*eta_x --> htemp1 */

    Dy(hu, htemp2);                /* eta_y --> htemp2 */
    Mult(htemp2, htemp2, htemp2);  /* eta_y*eta_y --> htemp2 */
    
    Sum(htemp1, htemp2, htemp1);   /* eta_x*eta_x + eta_y*eta_y --> htemp1  */
    Mult(htemp1, hwM2, htemp2);    /* W^(M-2)*(eta_x*eta_x + eta_y*eta_y) --> htemp2  */
    Sum(hwM, htemp2, htemp2);      /* W^(M-2)*(eta_x*eta_x + eta_y*eta_y) + W^(M) --> htemp2  */
    
    Sum(hrhs, htemp2, hrhs);


    /*------------------*/
    /* RHS for phi part */
    Mult(htemp1, hw2M2, htemp2);   /* W2^(M-2)*(eta_x*eta_x + eta_y*eta_y) --> htemp2  */
    Sum(hw2M, htemp2, htemp2);     /* W2^(M-2)*(eta_x*eta_x + eta_y*eta_y) + W2^(M) --> htemp2  */
    
    for (i=0; i<local_N; i++) {
        
        index = i;
        index_shift = index + alloc_local;
        hrhs[index_shift] = -g*hu[index];
        hrhs[index_shift] = hrhs[index_shift] + 0.5*htemp2[index];

    }
        
    /* -0.5*(phi_x*phi_x + phi_y*phi_y) */
    Dx(&hu[alloc_local], htemp1);
    Mult(htemp1, htemp1, htemp1);
  
    Dy(&hu[alloc_local], htemp2);
    Mult(htemp2, htemp2, htemp2);
    
    /* phi_x*phi_x + phi_y*phi_y --> htemp1  */
    Sum(htemp1, htemp2, htemp1);
    
    for (i=0; i<local_N; i++) {
        
        index = i;
        index_shift = index + alloc_local;
        hrhs[index_shift] = hrhs[index_shift] - 0.5*htemp1[index];

    }


    /* Dommermuth initialiation scheme -----------------------*/
    /* Multiply nonlinear part of rhs by ramping coefficient. */
    if (rampflg==1)
    {
    
        ZvelLinear(hu, htemp1);

        for (i=0; i<local_N; i++)
        {
            hrhs[i] = hrhs[i] - htemp1[i];
            hrhs[i] = RampFun(t)*hrhs[i];
            hrhs[i] = hrhs[i] + htemp1[i];
        }
    
        for (i=0; i<local_N; i++)
        {
            index_shift = i + alloc_local;
            hrhs[index_shift] = hrhs[index_shift] + g*hu[i];
            hrhs[index_shift] = RampFun(t)*hrhs[index_shift];
            hrhs[index_shift] = hrhs[index_shift] - g*hu[i];
        }
    
    }

    /* External forcing, directly added to hrhs --------------*/
    if (windflg==1)
        Wind(hu, hu, hrhs);
    
    /* Dealiasing --------------------------------------------*/
    Dealias(hrhs);
    Dealias(&hrhs[alloc_local]);
    
}



/* Compute vertical velocity on the free surface. */
/* HOS scheme (West et al. 1987) */
void Zvel(fftw_complex* hu, fftw_complex* hZvelM, fftw_complex* hZvelM2, fftw_complex* hZvel2M, fftw_complex* hZvel2M2, double t){

    ptrdiff_t Nindex, Mindex;
    ptrdiff_t n, m;
    

    for (n=0; n<=NLevs; n++) {
        
        for (i=0; i<local_N; i++) {
    
                hwn[n*alloc_local + i] = 0.0;
                hphin[n*alloc_local + i] = 0.0;
            
        }
    
    }
    
    for (i=0; i<local_N; i++) {
        
        hZvelM[i] = 0.0;
        hZvelM2[i] = 0.0;
        hZvel2M[i] = 0.0;
        hZvel2M2[i] = 0.0;
        hetan[i] = 0.0;
        
    }

    if (mpi_rank==0) {
        hetan[0] = 1.0;
    }
    
    Dz( &hu[alloc_local], hphin);
    
    for (i=0; i<local_N; i++) {

            hwn[i] = hphin[i];
        
    }

    for (n=1; n<=NLevs; n++) {
    
        Nindex = n*alloc_local;
        
        Mult(hu, &hetan[(n-1)*alloc_local], &hetan[Nindex]);

        /* Phi_n */
        for (m=0; m<=n-1; m++) {
    
            Mindex = m*alloc_local;
            
            Mult( &hphin[Mindex], &hetan[Nindex - Mindex], htemp1);
                
            for (i=0; i<local_N; i++) {
                    
                hphin[Nindex + i] = hphin[Nindex + i] - htemp1[i]/Coeff[n-m];
                    
            }
    
        }
        
        /* W_n */
        for (m=0; m<=n; m++) {
        
            Mindex = m*alloc_local;
        
            Dz( &hphin[Mindex], &hphin[Mindex]);
            Mult( &hphin[Mindex], &hetan[Nindex - Mindex], htemp1);
            
            for (i=0; i<local_N; i++) {
                    
                hwn[Nindex + i] = hwn[Nindex + i] + htemp1[i]/Coeff[n-m];
                    
            }
    
        }

    }

    /* hZvelM2 */
    for (n=0; n<=(NLevs-2); n++) {
            
        for (i=0; i<local_N; i++) {

            hZvelM2[i] = hZvelM2[i] + hwn[n*alloc_local + i];
                
        }
    
    }
    
    /* hZvelM */
    for (i=0; i<local_N; i++) {

        hZvelM[i] = hZvelM2[i];
                
    }
    
    for (n=(NLevs-1); n<=NLevs; n++) {
                
        for (i=0; i<local_N; i++) {
        
            hZvelM[i] = hZvelM[i] + hwn[n*alloc_local + i];
                
        }

    }

    /* hZvel2M2 */
    for (n=0; n<=(NLevs-3); n++) {
        
//        for (i=0; i<local_N; i++) {
//            
//            htemp2[i] = 0.0;
//            
//        }


        /* W2^(n) = W^(0)*W^(n) + W^(1)*W^(n-1) + ... W^(n)*W^(0) */
        for (m=0; m<=n; m++) {
        
            Mult(&hwn[m*alloc_local], &hwn[(n-m)*alloc_local], htemp1);
            
//            for (i=0; i<local_N; i++) {
//                
//                htemp2[i] = htemp2[i] + htemp1[i];
//                
//            }
          
            Sum(hZvel2M2, htemp1, hZvel2M2);
            
        }
        
        //Sum(hZvel2M2, htemp2, hZvel2M2);
        
    }
    
    
    /* hZvel2M */
    for (i=0; i<local_N; i++) {

        hZvel2M[i] = hZvel2M2[i];
                
    }
    
    for (n=(NLevs-2); n<=(NLevs-1); n++) {
        
//        for (i=0; i<local_N; i++) {
//                
//                htemp2[i] = 0.0;
//                
//        }


        /* W2^(n) = W^(0)*W^(n) + W^(1)*W^(n-1) + ... W^(n)*W^(0) */
        for (m=0; m<=n; m++) {
            
            Mult(&hwn[m*alloc_local], &hwn[(n-m)*alloc_local], htemp1);
            
//            for (i=0; i<local_N; i++) {
//                
//                htemp2[i] = htemp2[i] + htemp1[i];
//                
//            }
            Sum(hZvel2M, htemp1, hZvel2M);

        }
        
//        Sum(hZvel2M, htemp2, hZvel2M);
        
    }

    
}

/* Compute vertical velocity on the free surface at linear order. */
void ZvelLinear(const fftw_complex* hu, fftw_complex* hZvelLinear){

    Dz(&hu[alloc_local], hZvelLinear);
 
}

/* Compute the Hamiltonian for the HOS system, as given in West et al. 1987, eq. 8.*/
double Hamiltonian(const fftw_complex* heta, const fftw_complex* heta_t, const fftw_complex* hphi){

    ptrdiff_t index;
    double H;
    
    H = 0.0;
        
    /* j=0 */
    for (i=0; i < Nx; i++)
    {
        H += g*creal(heta[i])*creal(heta[i]);
        H += g*cimag(heta[i])*cimag(heta[i]);

        H += creal(hphi[i])*creal(heta_t[i]);
        H += cimag(hphi[i])*cimag(heta_t[i]);
    }
    
    for (j=1; j<local_Nyhpo; j++)
    {
        for (i=0; i < Nx; i++)
        {
            index = Nx*j + i;

            H += 2*g*creal(heta[index])*creal(heta[index]);
            H += 2*g*cimag(heta[index])*cimag(heta[index]);

            H += 2*creal(hphi[index])*creal(heta_t[index]);
            H += 2*cimag(hphi[index])*cimag(heta_t[index]);
        }
    }
    
    return H;

}

/* Ramping function used for Dommermuth adjusted initialization (D. Dommermuth, Wave Motion, Vol 32, 2000). */
double RampFun(const double t){
    
    return ( 1 - exp(-t*t/Tramp/Tramp) );

}

/* Add wind forcing to rhs */
void Wind(const fftw_complex* hu, const fftw_complex* hut, fftw_complex* hrhs){

    ptrdiff_t i;
    ptrdiff_t index_shift;

    /* implemetation of wind pressure */
    for (i=0; i<local_N; i++)
    {
        hwindpress[i] = 0.4*hu[i];
    }
    
    /* Adding wind pressure to hrhs */
    for (i=0; i<local_N; i++)
    {
        index_shift = i + alloc_local;
        hrhs[index_shift] = hrhs[index_shift] + hwindpress[i];
    }

}



