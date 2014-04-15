//
//  operators.c
//  euler
//
//  Created by Claudio Viotti on 7/10/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "operators_par.h"

/* Spectral x-derivative */
void Dx(fftw_complex* hu, fftw_complex* hu_x){


    double kx;
    ptrdiff_t index;

#if FFT_TRANSPOSE == 0
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny/2+1; j++) {
            
            if ( (local_0_start + i) < Nx/2 ) {
                
                kx = (local_0_start + i)*Kx0;
            
            }
            else{
            
                kx = -(Nx - local_0_start - i)*Kx0;
            
            }

            hu_x[(Ny/2+1)*i + j] = I*kx*hu[(Ny/2+1)*i + j];
            
        }
        
    }
    
#elif FFT_TRANSPOSE == 1
    for (j=0; j<local_Nyhpo; j++) {
    
        for (i=0; i < Nx; i++) {
        
            if ( i <= (Nx/2 + 1) ) {
                
                kx = i*Kx0;
            
            }
            else{
            
                kx = -(Nx - i)*Kx0;
            
            }
            
            index = Nx*j + i;
            hu_x[index] = I*kx*hu[index];
            
        }
        
    }

#endif

}


/* Spectral y-derivative */
void Dy(fftw_complex* hu, fftw_complex* hu_y){

    double ky;
    ptrdiff_t index;
    
#if FFT_TRANSPOSE == 0
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny/2+1; j++) {
            
            ky = j*Ky0;
            hu_y[(Ny/2+1)*i + j] = I*ky*hu[(Ny/2+1)*i + j];

        }
        
    }
    
#elif FFT_TRANSPOSE == 1
    for (j=0; j<local_Nyhpo; j++) {
    
        for (i=0; i < Nx; i++) {
        
            ky = (local_1_start + j)*Ky0;
            index = Nx*j + i;
            hu_y[index] = I*ky*hu[index];
            
        }
        
    }
    
#endif

}


/* Spectral z-derivative for deep water */
void Dz(const fftw_complex* hu, fftw_complex* hu_z){


    ptrdiff_t index;
    double kx, ky, kz;

#if FFT_TRANSPOSE == 0
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny/2+1; j++) {
            
            if ( (local_0_start + i) <= (Nx/2 + 1) ) {
                
                kx = (local_0_start + i)*Kx0;
            
            }
            else{
            
                kx = -(Nx - local_0_start - i)*Kx0;
            
            }
            ky = j*Ky0;
            kz = kx*kx + ky*ky;
            kz = sqrt(kz);
    
            index = (Ny/2+1)*i + j;
            hu_z[index] = kz*hu[index];

        }
        
    }
    
#elif FFT_TRANSPOSE == 1    
    for (j=0; j<local_Nyhpo; j++) {
    
        for (i=0; i<Nx; i++) {
        
            if ( i <= (Nx/2 + 1) ) {
                
                kx = i*Kx0;
            
            }
            else{
            
                kx = -(Nx - i)*Kx0;
            
            }
            ky = (local_1_start + j)*Ky0;
            kz = kx*kx + ky*ky;
            kz = sqrt(kz);

            index = Nx*j + i;
            hu_z[index] = kz*hu[index];
            
        }
        
    }
    
#endif

}



/* Compute product between two arrays. */
/* The operation can be executed fully in-place hu1=hu2=hprod. */
void Mult(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hprod){
    
    
    ptrdiff_t index;
    
    
    ifft_2d(hu1, temp1, ifftp);
    ifft_2d(hu2, temp2, ifftp);


    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
    
            index = (Ny + 2)*i + j;
            temp1[index] = temp1[index]*temp2[index];
    
        }
    
    }
    
    fft_2d(temp1, hprod, fftp);
    
}


void Sum(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hsum){
    
    
    ptrdiff_t index, Nyhal;
    
    Nyhal = Ny/2+1;
#if FFT_TRANSPOSE == 0    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Nyhal; j++) {
    
            index = Nyhal*i + j;
            hsum[index] = hu1[index] + hu2[index];
    
        }
        
    }
#elif FFT_TRANSPOSE == 1
    for (j=0; j<local_Nyhpo; j++) {
    
        for (i=0; i<Nx; i++) {
        
                index = Nx*j + i;
                hsum[index] = hu1[index] + hu2[index];
            
        }
        
    }
#endif
    
}


void Dealias(fftw_complex* hu){

    //int mx, my;
    //mx = floor( 0.5*Nx/(1 + 0.5*NLevs) );
    //my = floor( 0.5*Ny/(1 + 0.5*NLevs) );

#if FFT_TRANSPOSE == 0
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny/2+1; j++) {
            
            if ( ( (local_0_start + i)>=mx && (local_0_start + i)<(Nx-mx+1)) || (j>=my) ) {

                hu[(Ny/2+1)*i + j] = 0.0;
            
            }
            
        }
        
    }


#elif FFT_TRANSPOSE == 1
    for (j=0; j<local_Nyhpo; j++) {
    
        for (i=0; i<Nx; i++) {
        
            if ( ( i>=mx && i<(Nx-mx+1)) || ((local_1_start + j)>=my) ) {

                hu[Nx*j + i] = 0.0;
            
            }
            
        }
        
    }
#endif

}


/*---------------------------------------------------------------------*/
double ij2kx(int i, int j){

    double kx;

#if FFT_TRANSPOSE == 0
    if ( (local_0_start + i) <= (Nx/2 + 1) ) {
                
        kx = (local_0_start + i)*Kx0;
            
    }
    else{
            
        kx = -(Nx - local_0_start - i)*Kx0;
            
    }
#elif FFT_TRANSPOSE == 1
    if ( i <= (Nx/2 + 1) ) {
                
        kx = i*Kx0;
            
    }
    else{
            
        kx = -(Nx - i)*Kx0;
            
    }
#endif

    return kx;

};


/*---------------------------------------------------------------------*/
/* Filter operator as implemented by Xiao et al. JFM 2013.             */
void Filter(fftw_complex* hu, double k_peak, double fb1, double fb2){

double kx, ky, kmod;
ptrdiff_t index;

#if FFT_TRANSPOSE == 0
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny/2+1; j++) {
            
            kx = ij2kx(i,j);
            ky = j*Ky0;
            kmod = kx*kx;
            kmod += ky*ky;
            kmod = sqrt(kmod);
    
            index = (Ny/2+1)*i + j;
            hu[index] *= exp( - pow(kmod/fb1/k_peak, fb2) );
            
        }
        
    }
#elif FFT_TRANSPOSE == 1
    for (j=0; j<local_Nyhpo; j++) {
    
        for (i=0; i<Nx; i++) {
        
            if ( i <= (Nx/2 + 1) ) {
                
                kx = i*Kx0;
            
            }
            else{
            
                kx = -(Nx - i)*Kx0;
            
            }
            ky = (local_1_start + j)*Ky0;
            kmod = kx*kx;
            kmod += ky*ky;
            kmod = sqrt(kmod);

            index = Nx*j + i;
            hu[index] *= exp( - pow(kmod/fb1/k_peak, fb2) );
        }
        
    }
#endif

};

