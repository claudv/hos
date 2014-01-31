//
//  main_2d.c
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#define DIM 3
#define BUFSIZE 256
#define INIT_DATA_BUFSIZE 128
#define EPSILON 0.00000001

#define THREADS 0
#define N_THREADS 2


#include "euler_2d_par.h"

#include "fft_routines_2d_par.h"
#include "hdf5_routines_2d_par.h"
#include "model_2d_par.h"
#include "time_schemes_par.h"
#include "operators_par.h"



#if THREADS == 1
   int threads_ok;
#endif

int main(int argc, char **argv)
{
    
    
    int         scheme_flg, nfld;
    double      *eta, *phi;
    double      *velwM, *velwM2;
    //double      *velw2M, *velw2M2;
    double      t, dt, t_old;
    char        init_data_buff[INIT_DATA_BUFSIZE], init_pars_buff[INIT_DATA_BUFSIZE], subid_buff[2], nfld_buff[5], dtflag[20];
    
    fftw_complex    *heta, *hphi;
    fftw_complex    *hvelwM, *hvelwM2, *hvelw2M, *hvelw2M2;

    
    comm  = MPI_COMM_WORLD;
    info  = MPI_INFO_NULL;

#if THREADS == 0
    MPI_Init(&argc, &argv);
#elif THREADS == 1
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    threads_ok = provided >= MPI_THREAD_FUNNELED;
#endif
    
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    
    strncpy(init_pars_buff,"\0",INIT_DATA_BUFSIZE);
    //strncpy(init_pars_buff, Sim_root, strlen(Sim_root));
    strncpy(init_pars_buff,"initpars.h5", strlen("initpars.h5"));
    
    /* Setup runtime data file */
    strncpy(runtime_data_buff,"\0",RUNTIME_DATA_BUFSIZE);
    strncpy(runtime_data_buff,"runtime.dat", strlen("runtime.dat"));
    runtime_fid = fopen(runtime_data_buff, "w+");
    fclose(runtime_fid);
    
    
    /* Read input parameters file */
    get_params(init_pars_buff);
    
    snprintf(subid_buff, 2,"%d",runsubid);
    
    strncpy(init_data_buff,"\0",INIT_DATA_BUFSIZE);
    strcat(init_data_buff,"initdata.");
    strcat(init_data_buff,subid_buff);
    strcat(init_data_buff,".h5");
    
    if (mpi_rank==0) {
        printf("Reading initial data from:\t%s\n",init_data_buff);
    }
    MPI_Barrier(comm);
   
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

    
    fNx=Nx;
    fNy=Ny;
    
    
#if THREADS == 0
    fftw_mpi_init();
#elif THREADS == 1
    if (threads_ok) threads_ok = fftw_init_threads();
#endif
    
    
    /* get local data size and allocate. */
    /* NOTE: alloc_local can be greater than local_Nx*(Ny/2+1) */
#if FFT_TRANSPOSE == 0
    alloc_local = fftw_mpi_local_size_2d(fNx, fNy/2+1, MPI_COMM_WORLD, &local_Nx, &local_0_start);
    local_N = local_Nx*(Ny/2 + 1);
#elif FFT_TRANSPOSE == 1
    alloc_local = fftw_mpi_local_size_2d_transposed(fNx, fNy/2+1, MPI_COMM_WORLD, &local_Nx, &local_0_start, &local_Nyhpo, &local_1_start);
    local_N = Nx*local_Nyhpo;
#endif
    
    
    printf("Process %d:\tlocal size = %td,\tlocal x-size = %td,\tlocal y-size (transp-out) = %td.\n",mpi_rank,alloc_local,local_Nx,local_Nyhpo);


    eta = fftw_alloc_real(4 * alloc_local);
    heta = fftw_alloc_complex(2 * alloc_local);

    phi = &eta[2 * alloc_local];
    hphi = &heta[alloc_local];

    f = fftw_alloc_real(2 * alloc_local);
    hf = fftw_alloc_complex(alloc_local);
    
    
    temp1 = fftw_alloc_real(2 * alloc_local);
    velwM    = fftw_alloc_real(2 * alloc_local);
    velwM2   = fftw_alloc_real(2 * alloc_local);
    //velw2M   = fftw_alloc_real(2 * alloc_local);
    //velw2M2  = fftw_alloc_real(2 * alloc_local);
    hvelwM   = fftw_alloc_complex(alloc_local);
    hvelwM2  = fftw_alloc_complex(alloc_local);
    hvelw2M  = fftw_alloc_complex(alloc_local);
    hvelw2M2 = fftw_alloc_complex(alloc_local);
    
    
    /* Setup fft routines */
#if THREADS == 1
    if (threads_ok) fftw_plan_with_nthreads(N_THREADS);
#endif
    
    
#if FFT_TRANSPOSE == 0
    fftp = fftw_mpi_plan_dft_r2c_2d(fNx, fNy, f, hf, MPI_COMM_WORLD, FFTW_MEASURE);
    ifftp = fftw_mpi_plan_dft_c2r_2d(fNx, fNy, hf, f, MPI_COMM_WORLD, FFTW_MEASURE);
#elif FFT_TRANSPOSE == 1
    fftp = fftw_mpi_plan_dft_r2c_2d(fNx, fNy, f, hf, MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);
    ifftp = fftw_mpi_plan_dft_c2r_2d(fNx, fNy, hf, f, MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_IN);
#endif
    

    /* Read initial data */
    get_ic_2d(init_data_buff, eta, phi);

    /* Setup temporal scheme. */
    scheme_flg=2;
    Setup_TimeScheme(scheme_flg);
    rhs_hos_setup();
    
    t = 0;
    
    fft_2d(eta, heta, fftp);
    fft_2d(phi, hphi, fftp);
    
    Dealias(heta);
    Dealias(hphi);

    Zvel(heta, hvelwM, hvelwM2, hvelw2M, hvelw2M2, t);
    
    ifft_2d(hvelwM, velwM, ifftp);
    ifft_2d(hvelwM2, velwM2, ifftp);
    
    
    /* Save initial snapshot */
    savefileid = create_file_2d(savefile_buff);
    write_header_2d(savefileid, t);
    write_field_2d(savefileid, eta, phi);
    status = close_file_2d(savefileid);
    
    if (mpi_rank == 0) {
            printf("Datafile '%s' written at t=%f\n",savefile_buff,t);
    }
    
    
    if (saveflg > 1) {
        savefileid2 = create_file_2d(savefile2_buff);
        write_header_2d(savefileid2, t);
        write_extra_2d(savefileid2, velwM, velwM2);
        status = close_file_2d(savefileid2);
        
        if (mpi_rank == 0) {
            printf("Datafile '%s' written at t=%f\n",savefile2_buff,t);
        }
        
    }

    
    /* Main time loop */
    while (t<T-EPSILON) {
        
        t_old = t;
        
        dt = 0.2;
        
        sol_update_RK(heta,&t,dt,dtflag);
        //t = t + dt;
    
        if ( floor(t*(1+EPSILON)/dtsave) > floor(t_old*(1+EPSILON)/dtsave) ){
            
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
    
            savefileid = create_file_2d(savefile_buff);
            write_header_2d(savefileid, t);
            write_field_2d(savefileid, eta, phi);
            status = close_file_2d(savefileid);
            
            if (mpi_rank==0) {
                printf("Datafile '%s' written at t=%f\n",savefile_buff,t);
            }
            
            if (saveflg > 1) {
            
                strncpy(savefile2_buff,"\0",SAVE_FILE_BUFSIZE);
                strcat(savefile2_buff,"data_extra");
                snprintf(nfld_buff, 5,"%d",nfld);
                strcat(savefile2_buff,nfld_buff);
                strcat(savefile2_buff,".");
                strcat(savefile2_buff,subid_buff);
                strcat(savefile2_buff,".h5");

                savefileid2 = create_file_2d(savefile2_buff);
                write_header_2d(savefileid2, t);
                write_extra_2d(savefileid2, velwM, velwM2);
                status = close_file_2d(savefileid2);
                
                if (mpi_rank==0) {
                    printf("Datafile '%s' written at t=%f\n",savefile2_buff,t);
                }

            }
            
        }
        
    }
    
    
    
    fftw_destroy_plan(fftp);
    fftw_destroy_plan(ifftp);
    free(eta);
    fftw_free(heta);
    
    MPI_Finalize();
    
    return 0;
    
}

