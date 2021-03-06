//
//  main_2d.c
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

/*============================================================
 pending maintenance issues:
  -find a better arrangement for global scratchwork arrays
  -f and hf should not be visible outside fft_routines_2d_par.c
============================================================*/

#define DIM 3
#define BUFSIZE 256
#define INIT_DATA_BUFSIZE 128
#define EPSILON 0.00000001

#define THREADS 0
#define N_THREADS 2


#include "euler_2d_par.h"

#include "fft_routines_2d_par.h"
#include "hdf5_routines_2d_par.h"
#include "global_params.h"
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
    //double      *velwM, *velwM2;
    //double      *velw2M, *velw2M2;
    double      *eta_t;
    double      *press;
    
    double      t, dt, t_old;
    char        init_data_buff[INIT_DATA_BUFSIZE];
    char        init_pars_buff[INIT_DATA_BUFSIZE];
    char        savefile_buff[SAVE_FILE_BUFSIZE];
    char        savefile2_buff[SAVE_FILE_BUFSIZE];
    
    herr_t      status;
    
    char        subid_buff[2], nfld_buff[5], dtflag[20];
    
    fftw_complex    *heta, *hphi;
    //fftw_complex    *hvelwM, *hvelwM2, *hvelw2M, *hvelw2M2;
    fftw_complex    *heta_t;
    fftw_complex    *hpress;

    hid_t       savefileid;
    hid_t       savefileid2;
    
    /*-----------------------------------------------------*/
    /* Initialize MPI                                      */
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
    
    /*-----------------------------------------------------*/
    /* Setup runtime data file                             */
    strncpy(runtime_data_buff,"\0",RUNTIME_DATA_BUFSIZE);
    strncpy(runtime_data_buff,"runtime.dat", strlen("runtime.dat"));
    runtime_fid = fopen(runtime_data_buff, "w+");
    fclose(runtime_fid);
    
    /*----------------------------------------------------*/
    /* Read input parameters                              */
    strncpy(init_pars_buff,"\0",INIT_DATA_BUFSIZE);
    strncpy(init_pars_buff,"initpars.h5", strlen("initpars.h5"));
    get_global_params(init_pars_buff);
    
    if (mpi_rank==0)
    {
        printf("windflg = %d\n",windflg);
        if(windflg==1)
            printf("Uwind_x = %lf, Uwind_y = %lf\n",Uwind_x,Uwind_y);
    }
    
    /*----------------------------------------------------*/
    /* Setup file names                                   */
    strncpy(init_data_buff,"\0",INIT_DATA_BUFSIZE);
    strcat(init_data_buff,"initdata.");
    snprintf(subid_buff, 2,"%d",runsubid);
    strcat(init_data_buff,subid_buff);
    strcat(init_data_buff,".h5");
    
    if (mpi_rank==0)
        printf("Reading initial data from:\t%s\n",init_data_buff);
    
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
    
    

    /*---------------------------------------------------------*/
    /* Get local data size and allocate.                       */
    /* NOTE: alloc_local can be greater than local_Nx*(Ny/2+1) */
#if THREADS == 0
    fftw_mpi_init();
#elif THREADS == 1
    if (threads_ok) threads_ok = fftw_init_threads();
#endif

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

    eta_t = fftw_alloc_real(4 * alloc_local);
    heta_t = fftw_alloc_complex(2 * alloc_local);

    press = fftw_alloc_real(2 * alloc_local);
    hpress = fftw_alloc_complex(alloc_local);

    phi = &eta[2 * alloc_local];
    hphi = &heta[alloc_local];

    f = fftw_alloc_real(2 * alloc_local);
    hf = fftw_alloc_complex(alloc_local);
    
    
    temp1 = fftw_alloc_real(2 * alloc_local);
    //velwM    = fftw_alloc_real(2 * alloc_local);
    //velwM2   = fftw_alloc_real(2 * alloc_local);
    //velw2M   = fftw_alloc_real(2 * alloc_local);
    //velw2M2  = fftw_alloc_real(2 * alloc_local);
    //hvelwM   = fftw_alloc_complex(alloc_local);
    //hvelwM2  = fftw_alloc_complex(alloc_local);
    //hvelw2M  = fftw_alloc_complex(alloc_local);
    //hvelw2M2 = fftw_alloc_complex(alloc_local);
    
    /*------------------------------------------------*/
    /* Setup fftw plans                               */
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

    /*------------------------------------*/
    /* Read initial data                  */
    get_ic_2d(init_data_buff, eta, phi);

    /*------------------------------------*/
    /* Setup temporal scheme.             */
    scheme_flg=2;
    Setup_TimeScheme(scheme_flg);
    rhs_hos_setup();
    
    t = 0;
    
    fft_2d(eta, heta, fftp);
    fft_2d(phi, hphi, fftp);
    
    Dealias(heta);
    Dealias(hphi);
    
    /* Save initial snapshot */
    savefileid = create_file_2d(savefile_buff);
    write_header_2d(savefileid, t);
    write_field_2d(savefileid, eta, phi);
    status = close_file_2d(savefileid);
    
    //savefileid2 = create_file_2d(savefile2_buff);
	//write_header_2d(savefileid, t);
    //write_field_complex_2d(savefileid2, heta, hphi);
    //status = close_file_2d(savefileid2);
    
    if (mpi_rank == 0)
            printf("Datafile '%s' written at t=%f\n",savefile_buff,t);
    
    
    if (saveflg > 1)
    {
        rhs_hos(heta_t, heta, t);
        ifft_2d(heta_t, eta_t, ifftp);
        Wind(heta, heta_t, hpress);
        ifft_2d(hpress, press, ifftp);
        //Zvel(heta, hvelwM, hvelwM2, hvelw2M, hvelw2M2, t);
        //ifft_2d(hvelwM, velwM, ifftp);
        //ifft_2d(hvelwM2, velwM2, ifftp);
    
        savefileid2 = create_file_2d(savefile2_buff);
        write_header_2d(savefileid2, t);
        //write_extra_2d(savefileid2, velwM, velwM2);
        write_extra_2d(savefileid2, eta_t, press);
        status = close_file_2d(savefileid2);
        
        if (mpi_rank == 0)
            printf("Datafile '%s' written at t=%f\n",savefile2_buff,t);
    }
    
    /*--------------------*/
    /* Main time loop     */
    while (t<T-EPSILON)
    {
        t_old = t;
        
        dt = 0.02;
        
        sol_update_RK(heta,&t,dt,dtflag);
        
        Filter(heta, 4.0, 8.0, 30.0);
        Filter(hphi, 4.0, 8.0, 30.0);
     
    
        if ( floor(t*(1+EPSILON)/dtsave) > floor(t_old*(1+EPSILON)/dtsave) )
        {
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
            
            if (mpi_rank==0)
                printf("Datafile '%s' written at t=%f\n",savefile_buff,t);
            
            if (saveflg > 1)
            {
            
                strncpy(savefile2_buff,"\0",SAVE_FILE_BUFSIZE);
                strcat(savefile2_buff,"data_extra");
                snprintf(nfld_buff, 5,"%d",nfld);
                strcat(savefile2_buff,nfld_buff);
                strcat(savefile2_buff,".");
                strcat(savefile2_buff,subid_buff);
                strcat(savefile2_buff,".h5");

                rhs_hos(heta_t, heta, t);
                ifft_2d(heta_t, eta_t, ifftp);
                Wind(heta, heta_t, hpress);
                ifft_2d(hpress, press, ifftp);
        
                savefileid2 = create_file_2d(savefile2_buff);
                write_header_2d(savefileid2, t);
                write_extra_2d(savefileid2, eta_t, press);
                status = close_file_2d(savefileid2);
                
                if (mpi_rank==0)
                    printf("Datafile '%s' written at t=%f\n",savefile2_buff,t);

            }
        }
    }
    
    
    /*------------------------*/
    /* Clean up and finalize  */
    fftw_destroy_plan(fftp);
    fftw_destroy_plan(ifftp);
    free(eta);
    fftw_free(heta);
    
    MPI_Finalize();
    
    return 0;
    
}

