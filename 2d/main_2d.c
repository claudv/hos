//
//  main_2d.c
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "euler_2d.h"

#include "fft_routines_2d.h"
#include "hdf5_routines_2d.h"
#include "model_2d.h"
#include "time_schemes.h"

#define DIM 2
#define M 4
#define BUFSIZE 256
#define INIT_DATA_BUFSIZE 128


int main(int argc, char **argv)
{
    
    
    int         scheme_flg, nfld;
    double      *eta, *phi;
    double      *velwM, *velwM2, *velw2M, *velw2M2;
    double      t, dt, t_old;
    char        init_data_buff[INIT_DATA_BUFSIZE], init_pars_buff[INIT_DATA_BUFSIZE], subid_buff[2], dtflag[20], nfld_buff[5];
    
    fftw_complex    *heta, *hphi;
    fftw_complex    *hvelwM, *hvelwM2, *hvelw2M, *hvelw2M2;
    
    
    
    strncpy(init_pars_buff,"\0",INIT_DATA_BUFSIZE);
    //strncpy(init_pars_buff, Sim_root, strlen(Sim_root));
    strncpy(init_pars_buff,"initpars.h5", strlen("initpars.h5"));
    
    /* Read input parameters file */
    get_params(init_pars_buff);
    
    
    printf("Final time = %f\n",T);
    printf("dtsave = %f\n",dtsave);
    
    snprintf(subid_buff, 2,"%d",runsubid);
    
    strncpy(init_data_buff,"\0",INIT_DATA_BUFSIZE);
    strcat(init_data_buff,"initdata.");
    strcat(init_data_buff,subid_buff);
    strcat(init_data_buff,".h5");
   
    nfld = 0;
    strncpy(savefile_buff,"\0",SAVE_FILE_BUFSIZE);
    strcat(savefile_buff,"data");
    snprintf(nfld_buff, 5,"%d",nfld);
    strcat(savefile_buff,nfld_buff);
    strcat(savefile_buff,".");
    strcat(savefile_buff,subid_buff);
    strcat(savefile_buff,".h5");
    
    strncpy(savefile2_buff,"\0",SAVE_FILE_BUFSIZE);
    strcat(savefile2_buff,"data_extra");
    snprintf(nfld_buff, 5,"%d",nfld);
    strcat(savefile2_buff,nfld_buff);
    strcat(savefile2_buff,".");
    strcat(savefile2_buff,subid_buff);
    strcat(savefile2_buff,".h5");
    
    Nxl = 4*Nx/2;
    Nyl = 4*Ny/2;
    
    g=1;
    NLevs=M;
    
    /* Define a total size to share time-stepping routines with 1d case */
    /* N/2+1 = Nx*(Ny/2 + 1) */
    N = 2*Nx*(Ny/2+1) - 2;
    
    eta = (double*) fftw_malloc(sizeof(double)*2*Nx*Ny);
    heta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*Nx*(Ny/2+1));
    
    phi = &eta[Nx*Ny];
    hphi = &heta[Nx*(Ny/2+1)];
    
    //hdummy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*Nx*(Ny/2+1));
    
    //k = malloc (sizeof(double)*Nx*(Ny/2+1));
    
    f =  malloc(sizeof(double)*Nx*Ny);
    ff = malloc(sizeof(double)*Nxl*Nyl);
    hf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    hff = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nxl*(Nyl/2+1));
    
    velwM    =  malloc(sizeof(double)*Nx*Ny);
    velwM2   =  malloc(sizeof(double)*Nx*Ny);
    velw2M   =  malloc(sizeof(double)*Nx*Ny);
    velw2M2  =  malloc(sizeof(double)*Nx*Ny);
    hvelwM   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    hvelwM2  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    hvelw2M  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    hvelw2M2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
    
    /* Setup fft routines */
    //stat = fftw_init_threads();
    //fftw_plan_with_nthreads(4);
    
    fftp  = fftw_plan_dft_r2c_2d(Nx, Ny,  f, hf, FFTW_ESTIMATE);
    ifftp = fftw_plan_dft_c2r_2d(Nx, Ny, hf,  f, FFTW_ESTIMATE);
    
    fftp_large  = fftw_plan_dft_r2c_2d(Nxl, Nyl, ff, hff, FFTW_ESTIMATE);
    ifftp_large = fftw_plan_dft_c2r_2d(Nxl, Nyl, hff, ff, FFTW_ESTIMATE);
    
    
    /* Read initial data */
    printf("Reading initial condition from %s\n",init_data_buff);
    get_ic_2d(init_data_buff, eta, phi);
    

    /* Setup temporal scheme. */
    scheme_flg=2;
    Setup_TimeScheme(scheme_flg);
    rhs_hos_setup();
    
    t = 0;
       
    fft_2d(eta, heta, fftp);
    fft_2d(phi, hphi, fftp);
    
    ifft_2d(heta, eta, ifftp);
    ifft_2d(hphi, phi, ifftp);
  
    
    int mx, my, index;
    
    mx = floor( 0.5*Nx/(1 + 0.5*NLevs) );
    my = floor( 0.5*Ny/(1 + 0.5*NLevs) );
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            if ( (i>=mx && i<(Nx-mx+1)) || (j>=my) ) {
                
                index = (Ny/2+1)*i + j;
                heta[index] = 0.0;
                hphi[index] = 0.0;
            
            }
            
        }
        
    }

    
    Zvel(heta, hvelwM, hvelwM2, hvelw2M, hvelw2M2, t);
    ifft_2d(hvelwM, velwM, ifftp);
    ifft_2d(hvelwM2, velwM2, ifftp);
    
    
    /* Save initial snapshot */
    savefileid = create_file_2d(savefile_buff);
    write_header_2d(savefileid, t);
    write_field_2d(savefileid, eta, phi);
    status = close_file_2d(savefileid);
            
    savefileid2 = create_file_2d(savefile2_buff);
    write_header_2d(savefileid2, t);
    write_extra_2d(savefileid2, velwM, velwM2);
    status = close_file_2d(savefileid2);
    
    printf("Datafile written at t=%f\n",t);
       
    while (t<T-0.00000001) {
        
        t_old = t;
        
        dt = 0.025;
        
        sol_update_RK(heta,&t,dt,dtflag);
        //t = t + dt;
        
        if ( floor( (t*1.000000001)/dtsave) > floor((t_old*1.000000001)/dtsave) ) {
            
            nfld = nfld + 1;
            ifft_2d(heta, eta, ifftp);
            ifft_2d(hphi, phi, ifftp);
            
            
            /* Print field in output file */
            strncpy(savefile_buff,"\0",SAVE_FILE_BUFSIZE);
            strcat(savefile_buff,"data");
            snprintf(nfld_buff, 5,"%d",nfld);
            strcat(savefile_buff,nfld_buff);
            strcat(savefile_buff,".");
            strcat(savefile_buff,subid_buff);
            strcat(savefile_buff,".h5");
    
            strncpy(savefile2_buff,"\0",SAVE_FILE_BUFSIZE);
            strcat(savefile2_buff,"data_extra");
            snprintf(nfld_buff, 5,"%d",nfld);
            strcat(savefile2_buff,nfld_buff);
            strcat(savefile2_buff,".");
            strcat(savefile2_buff,subid_buff);
            strcat(savefile2_buff,".h5");

            savefileid = create_file_2d(savefile_buff);
            write_header_2d(savefileid, t);
            write_field_2d(savefileid, eta, phi);
            status = close_file_2d(savefileid);
            
            savefileid2 = create_file_2d(savefile2_buff);
            write_header_2d(savefileid2, t);
            write_extra_2d(savefileid2, velwM, velwM2);
            status = close_file_2d(savefileid2);

            printf("Datafile written at t=%f\n",t);
            
        }
        
    }
    
    
//
//    ifft_1d(heta, eta, ifftp);
//    ifft_1d(hphi, phi, ifftp);
//    
//    
//    
    
    fftw_destroy_plan(fftp);
    fftw_destroy_plan(ifftp);
    free(eta);
    fftw_free(heta);
    
    fftw_destroy_plan(fftp_large);
    fftw_destroy_plan(ifftp_large);
    
    return 0;
    
}

