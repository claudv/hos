//
//  fft_routines.c
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "fft_routines_2d_par.h"


void fft_2d(const double* u, fftw_complex* hu, fftw_plan plan){
    
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            /* Copy real data inside f skipping padded elements 
            (note the Ny+2 striding over the first dimension) */
            f[(fNy+2)*i + j]=u[(fNy+2)*i + j];
            
        }
        
    }
    
    
    fftw_execute(plan);
    
    
//    for (i=0; i<local_Nx; i++) {
//    
//        for (j=0; j<Ny/2+1; j++) {
//            
//            hu[(fNy/2+1)*i + j]=hf[(fNy/2+1)*i + j]/Nx/Ny;
//            
//        }
//        
//    }
    
    
    
    for (i=0; i<local_N; i++) {
        
        hu[i] = hf[i]/Nx/Ny;
        
    }
    
}


void ifft_2d(const fftw_complex* hu, double* u, fftw_plan plan){
    
    
//    for (i=0; i<local_Nx; i++) {
//        
//        for (j=0; j<Ny/2+1; j++) {
//            
//            hf[(Ny/2+1)*i + j]=hu[(Ny/2+1)*i + j];
//            
//        }
//        
//    }

    for (i=0; i<local_N; i++) {
        
        hf[i] = hu[i];
            
    }
    
    
    fftw_execute(plan);
    
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            /* Copy real data inside f skipping padded elements 
            (note the Ny+2 striding over the first dimension) */
            u[(Ny+2)*i + j]=f[(Ny+2)*i + j];
            
        }
        
    }
    
}


