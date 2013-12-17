//
//  time_schemes.c
//  euler
//
//  Created by Claudio Viotti on 3/9/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "time_schemes.h"


void sol_update_RK(fftw_complex* u,double* t,double dt,char* dtflag){
    
    
    for (int s=0; s<Nstages; s++) {
    
    
        for (int i=0; i<N+2; i++) {
            
            htemp[i] = u[i];
            
        }
      
        
        for (int p=0; p<s; p++) {

            for (int i=0; i<N+2; i++) {
            
                htemp[i] = htemp[i] + dt*b[s][p]*fun[p][i];
            
            }
            
        }
     
        //rhs_confmap(fun[s], htemp);
        //rhs_test(fun[s], htemp);
        rhs_hos(fun[s], htemp, t[0]);
    
    }
    
    
    for (int s=0; s<Nstages; s++) {
    
        for (int i=0; i<N+2; i++) {
                
            u[i] = u[i] + dt*c[s]*fun[s][i];
        
        }
        
    }

//       for (int i=0; i<N+2; i++) {
//                
//            u[i] = 1.1*u[i];
//        
//        }

//    rhs(fun[0],u);
//    
//    for (int i=0; i<N/2+1; i++) {
//        
//        u[i] = u[i] + dt*fun[0][i];
//
//    }
    
    t[0] = t[0] + dt;
    
    
}



void Setup_TimeScheme(int scheme_flg){
    
    
    switch (scheme_flg) {
        case 1:
            Nstages = 7;
        break;

        case 2:
            Nstages = 4;
        break;

        default:
        break;

    }
    
    
    fun     = (fftw_complex**) fftw_malloc(sizeof(fftw_complex)*Nstages);
    funElem = (fftw_complex*)  fftw_malloc(sizeof(fftw_complex)*2*(N/2+1)*Nstages);
            
    htemp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(N/2+1));
            
    a = malloc (sizeof(double)*Nstages);
    b = malloc (sizeof(double*)*Nstages);
    bElem = malloc (sizeof(double)*(Nstages*Nstages));
    c =  malloc (sizeof(double)*Nstages);
    cs = malloc (sizeof(double)*Nstages);
            
    for (int i=0; i<Nstages; i++) {
                
        //b[i] = bElem + (Nstages*i)*sizeof(double);
        b[i] = &bElem[Nstages*i];
        //fun[i] = funElem + (N+2)*i*sizeof(fftw_complex);
        fun[i] = &funElem[(N+2)*i];
        
    }

    
    
    switch (scheme_flg) {
 
        case 1:
            
            /* Dorman-Prince-Shampine RK45 pair */
            a[0] = 0.0;
            a[1] = 1.0/5.0;
            a[2] = 3.0/10.0;
            a[3] = 4.0/5.0;
            a[4] = 8.0/9.0;
            a[5] = 1.0;
            a[6] = 1.0;
            
            b[1][0] = 1.0/5.0;
            b[2][0] = 3.0/40.0;
            b[3][0] = 44.0/45.0;
            b[4][0] = 19372.0/6561.0;
            b[5][0] = 9017.0/3168.0;
            b[6][0] = 35.0/384.0;
            
            b[2][1] = 9.0/40.0;
            b[3][1] = -56.0/15.0;
            b[4][1] = -25360.0/2187.0;
            b[5][1] = -355.0/33.0;
            b[6][1] = 0.0;
            
            b[3][2] = 32.0/9.0;
            b[4][2] = 64448.0/6561.0;
            b[5][2] = 46732.0/5247.0;
            b[6][2] = 500.0/1113.0;
            
            b[4][3] = -212.0/729.0;
            b[5][3] = 49.0/176.0;
            b[6][3] = 125.0/192.0;
            
            b[5][4] = -5103.0/18656.0;
            b[6][4] = -2187.0/6784.0;
            
            b[6][5] = 11.0/84.0;
            
            /* 4th order */
            c[0] = 1951.0/21600.0;
            c[1] = 0.0;
            c[2] = 22642.0/50085.0;
            c[3] = 451.0/720.0;
            c[4] = -12231.0/42400.0;
            c[5] = 649.0/6300.0;
            c[6] = 1.0/60.0;
            
            /* 5th order */
            cs[0] = 35.0/384.0;
            cs[1] = 0.0;
            cs[2] = 500.0/1113.0;
            cs[3] = 125.0/192.0;
            cs[4] = -2187.0/6784.0;
            cs[5] = 11.0/84.0;
            cs[6] = 0.0;

    
        break;
        
        
        case 2:
            
            a[0] = 0.0;
            a[1] = 0.5;
            a[2] = 0.5;
            a[3] = 1.0;
            
            b[1][0] = 0.5;
            b[2][0] = 0.0;
            b[3][0] = 0.0;
            
            b[2][1] = 0.5;
            b[3][1] = 0.0;
            
            b[3][2] = 1.0;
            
            c[0] = 1.0/6.0;
            c[1] = 1.0/3.0;
            c[2] = 1.0/3.0;
            c[3] = 1.0/6.0;


        break;


        default:
    
        break;
    
    }
    
    
        
}
