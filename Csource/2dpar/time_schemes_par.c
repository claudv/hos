//
//  time_schemes.c
//  euler
//
//  Created by Claudio Viotti on 3/9/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "time_schemes_par.h"

/* Coeffs for RK routines and other storage room */
static int             Nstages;
static double*         a;
static double*         bElem;
static double*         c;
static double*         cs;
static double**        b;
static fftw_complex    *htemp;
static fftw_complex    *fun;


void sol_update_RK(fftw_complex* u,double* t,double dt,char* dtflag){
    
    ptrdiff_t s, p;
    
    for (s=0; s<Nstages; s++) {
    
        
//        for (i=0; i<local_Nx; i++) {
//    
//            for (j=0; j<Ny/2+1; j++) {
//            
//                htemp[(Ny/2+1)*i + j]=u[(Ny/2+1)*i + j];
//            
//            }
//        
//        }

        for (i=0; i<local_N; i++) {
            
                htemp[i]=u[i];

        }


//        for (i=0; i<local_Nx; i++) {
//    
//            for (j=0; j<Ny/2+1; j++) {
//            
//                htemp[alloc_local + (Ny/2+1)*i + j]=u[alloc_local + (Ny/2+1)*i + j];
//            
//            }
//        
//        }

        for (i=0; i<local_N; i++) {
            
                htemp[alloc_local + i]=u[alloc_local + i];

        }
        
        
        for (p=0; p<s; p++) {
            
//            for (i=0; i<local_Nx; i++) {
//                
//                for (j=0; j<Ny/2+1; j++) {
//                    
//                    htemp[(Ny/2+1)*i + j] = htemp[(Ny/2+1)*i + j] + dt*b[s][p]*fun[2*alloc_local*p + (Ny/2+1)*i + j];
//                    
//                }
//                
//            }
            

            for (i=0; i<local_N; i++) {
            
                htemp[i] = htemp[i] + dt*b[s][p]*fun[2*alloc_local*p + i];

            }
            
            
//            for (i=0; i<local_Nx; i++) {
//                
//                for (j=0; j<Ny/2+1; j++) {
//                    
//                    htemp[alloc_local + (Ny/2+1)*i + j] = htemp[alloc_local +(Ny/2+1)*i + j] + dt*b[s][p]*fun[2*alloc_local*p + alloc_local + (Ny/2+1)*i + j];
//                    
//                }
//                
//            }


            for (i=0; i<local_N; i++) {
            
                htemp[alloc_local + i] = htemp[alloc_local + i] + dt*b[s][p]*fun[2*alloc_local*p + alloc_local + i];

            }
            
        }
     
        //rhs_test(&fun[2*alloc_local*s], htemp);
        rhs_hos(&fun[2*alloc_local*s], htemp, t[0]);
        
        /* Compute runtime diagnostics */
        if (s==0) {
  
            Ham = Hamiltonian(u, fun, &u[alloc_local]);
            MPI_Reduce(&Ham, &Ham_glob, 1, MPI_FLOAT, MPI_SUM, 0, comm);
            
            if (mpi_rank==0) {
  
                runtime_fid = fopen(runtime_data_buff, "a+");
                if ( runtime_fid == NULL ) {
                
                    perror ( "Unable to open runtime datafile" );
                    exit ( EXIT_FAILURE );
                    
                }
                
                fprintf(runtime_fid, "%14.10f\t%.16g\n", *t, Ham);
                fclose(runtime_fid);
            
            }

        }
    
    }
    
    
    for (s=0; s<Nstages; s++) {
        
//        for (i=0; i<local_Nx; i++) {
//                
//            for (j=0; j<Ny/2+1; j++) {
//                    
//                u[(Ny/2+1)*i + j] = u[(Ny/2+1)*i + j] + dt*c[s]*fun[2*alloc_local*s + (Ny/2+1)*i + j];
//                    
//            }
//                
//        }
        
        for (i=0; i<local_N; i++) {
            
                u[i] = u[i] + dt*c[s]*fun[2*alloc_local*s + i];

        }

            
            
//       for (i=0; i<local_Nx; i++) {
//            
//            for (j=0; j<Ny/2+1; j++) {
//                    
//                u[alloc_local + (Ny/2+1)*i + j] = u[alloc_local + (Ny/2+1)*i + j] + dt*c[s]*fun[2*alloc_local*s + alloc_local + (Ny/2+1)*i + j];
//                    
//            }
//                
//        }
        
        
        for (i=0; i<local_N; i++) {
            
                u[alloc_local + i] = u[alloc_local + i] + dt*c[s]*fun[2*alloc_local*s + alloc_local + i];

        }

        
        
    }

    
    t[0] = t[0] + dt;
    
    
}



void Setup_TimeScheme(int scheme_flg){
    
    int n;
    
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
    

    fun = fftw_alloc_complex(2 * Nstages * alloc_local);
    htemp = fftw_alloc_complex(2 * alloc_local);
        
    a = malloc (sizeof(double)*Nstages);
    b = malloc (sizeof(double*)*Nstages);
    bElem = malloc (sizeof(double)*(Nstages*Nstages));
    c =  malloc (sizeof(double)*Nstages);
    cs = malloc (sizeof(double)*Nstages);
            
    for (n=0; n<Nstages; n++) {
        
        b[n] = &bElem[Nstages*n];
        
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
            
            /* Standard RK4 */
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
