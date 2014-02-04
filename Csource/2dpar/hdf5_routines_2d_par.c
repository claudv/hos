//
//  hdf5_routines_2d.c
//  euler
//
//  Created by Claudio Viotti on 6/24/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "hdf5_routines_2d_par.h"
#include "fft_routines_2d_par.h"


void get_params(char* filename){
    
    hid_t       h5_file, h5_data;
    herr_t      status;
    double      Nd[1];
    
    
    h5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    /* Nx */
    h5_data = H5Dopen(h5_file, "/Nx", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    Nx = (int) Nd[0];
    
    /* Ny */
    h5_data = H5Dopen(h5_file, "/Ny", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    Ny = (int) Nd[0];
    
    /* Lx */
    h5_data = H5Dopen(h5_file, "/Lx", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lx);
    
    /* Ly */
    h5_data = H5Dopen(h5_file, "/Ly", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Ly);
    
    /* g  */
    h5_data = H5Dopen(h5_file, "/g", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &g);
    
    /* runsubid */
    h5_data = H5Dopen(h5_file, "/runsubid", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);

    runsubid = (int) Nd[0];

    /* T */
    h5_data = H5Dopen(h5_file, "/T", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
    
    /* dtsave */
    h5_data = H5Dopen(h5_file, "/dtsave", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dtsave);
    
    /* saveflg */
    h5_data = H5Dopen(h5_file, "/saveflg", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    saveflg = (int) Nd[0];
    
    /* rampflg */
    h5_data = H5Dopen(h5_file, "/rampflg", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    rampflg = (int) Nd[0];
    
    if ( rampflg!=0 && rampflg!=1 ) {
    
        printf("Illegal value assigned to rampflag, setting to default value (0).\n");
        rampflg = 0;
    
    }
    
    /* Tramp */
    h5_data = H5Dopen(h5_file, "/Tramp", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Tramp);
   
   
    status = H5Dclose(h5_data);
    status = H5Fclose(h5_file);
    
    Kx0 = 2*M_PI/Lx;
    Ky0 = 2*M_PI/Ly;
    
}


void get_ic_2d(char* filename, double* u1, double* u2){
    
    hid_t       h5_file, h5_dataset, h5_memspace, h5_filespace;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    
    
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    h5_file  = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_dataset = H5Dopen(h5_file, "/eta0", H5P_DEFAULT);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    /* Read data and write it inside a temporary array, then copy the data inside the padded target array */
    status = H5Dread(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, H5P_DEFAULT, temp1);
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            u1[(Ny+2)*i + j]=temp1[Ny*i + j];
            
        }
        
    }
    
    h5_dataset = H5Dopen(h5_file, "/phi0", H5P_DEFAULT);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    /* Read data and write it inside a temporary array, then copy the data inside the padded target array */
    status = H5Dread(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, H5P_DEFAULT, temp1);
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            u2[(Ny+2)*i + j]=temp1[Ny*i + j];
            
        }
        
    }
    
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    H5Fclose(h5_file);
    
}


hid_t create_file_2d(char* filename){
    
    hid_t       h5_file;
    hid_t       plist_id;
    
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
    
    h5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    
    
    return h5_file;
    
}


hid_t open_file_2d(char* filename){
    
    hid_t       h5_file;
    //hid_t       plist_id;
    
    //plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //H5Pset_fapl_mpio(plist_id, comm, info);
    
    h5_file = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT );
    //H5Pclose(plist_id);
    
    return h5_file;
    
}


herr_t close_file_2d(hid_t file){
    
    herr_t      status;
    
    
    status = H5Fclose(file);
    
    return status;
    
}


void write_header_2d(hid_t file, double time){
    
    hid_t       data, fid;
    herr_t      status;
    hsize_t     fdim [2];
    double      var;
    
    /* For some reason we get a deadlock of only one process executes this part */
    //if (mpi_rank==0) {
  
        fdim[0]=1;
        fdim[1]=1;
        
        fid = H5Screate_simple(2, fdim, NULL);
    
        var = time;
        data = H5Dcreate(file, "time", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&var);
        status = H5Dclose(data);

        var = g;
        data = H5Dcreate(file, "g", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&g);
        status = H5Dclose(data);

        var = (double) Nx;
        data = H5Dcreate(file, "Nx", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&var);
        status = H5Dclose(data);
        
        var = (double) Ny;
        data = H5Dcreate(file, "Ny", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&var);
        status = H5Dclose(data);
        
        data = H5Dcreate(file, "Lx", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lx);
        status = H5Dclose(data);
        
        data = H5Dcreate(file, "Ly", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Ly);
        status = H5Dclose(data);
        
        status = H5Sclose(fid);


    //}
    
}


void write_field_2d(hid_t h5_file, double* eta, double* phi){
    

    hid_t       h5_dataset, h5_memspace, h5_filespace;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "eta", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);

    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = eta[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
   
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "phi", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = phi[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
       
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
}



/* Function that checks whether the global index "index" lies inside (returns 0) or outside (returns 1) the dealiased region. */
ptrdiff_t IsOutside(ptrdiff_t index){ return (index >= mx)*(index <= Nx-mx-1); };


/* Local static defnition of a fftw plan to be used inside write_field_complex_2d() for untransposing fftw arrays. */
static fftw_plan fftp_untrsp;
static ptrdiff_t alloc_local_untrsp;
static ptrdiff_t local_Nx_untrsp;
static ptrdiff_t local_N_untrsp;
static ptrdiff_t local_0_start_untrsp;

enum initstatus{

    false,
    true

};


/* Initialize the fftw planning for transposed output on the first call of write_field_complex_2d(). */
void write_field_complex_2d_init(){

    static enum initstatus initialized = false;
    
    if (initialized==false) {
        
        fftp_untrsp = fftw_mpi_plan_dft_r2c_2d(fNx, fNy, f, hf, MPI_COMM_WORLD, FFTW_MEASURE);
        alloc_local_untrsp = fftw_mpi_local_size_2d(fNx, fNy/2+1, MPI_COMM_WORLD, &local_Nx_untrsp, &local_0_start_untrsp);
        local_N_untrsp = local_Nx_untrsp*(Ny/2 + 1);
        initialized = true;
    }

}


void write_field_complex_2d(hid_t h5_file, const fftw_complex* heta, const fftw_complex* phi){
    

    hid_t       h5_dataset, h5_memspace, h5_filespace;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    ptrdiff_t   lbd, rbd; // local indexes of the left and right boundaries of the portion of the local complex array contained inside the dealiased region.
    
    hsize_t *   count_0_gathered;
    count_0_gathered = malloc(mpi_size*sizeof(hsize_t));
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    
    /* Untranpose data if FFT_TRANSPOSE == 1, by doing ifft and untransposed fft */
#if FFT_TRANSPOSE == 1
    write_field_complex_2d_init();

    ifft_2d(heta, temp1, ifftp);
    
    for (i=0; i<local_Nx_untrsp; i++) {
        
        for (j=0; j<Ny; j++) {
            
            /* Copy real data inside f skipping padded elements. */
            /* (Note the Ny+2 striding over the first dimension) */
            f[(fNy+2)*i + j]=temp1[(fNy+2)*i + j];
            
        }
        
    }
    
    fftw_execute(fftp_untrsp);
    
    for (i=0; i<local_N_untrsp; i++) {
        
        htemp1[i] = hf[i]/Nx/Ny;
        
    }
#elif FFT_TRANSPOSE == 0
   for (i=0; i<local_N; i++) {
        
        htemp1[i] = heta[i];
        
    }
#endif
    

    lbd = ( Nx - mx - local_0_start)*IsOutside(local_0_start);
    rbd = local_Nx - 1 + (mx - local_0_start - local_Nx)*IsOutside(local_0_start+local_Nx-1);
    
    /* Compute count, the size of the local data slab to write. */
    /* If rbd<lbd this chuck has nothing. */
    //count[0] = (rbd - lbd + 1)*(rbd - lbd > 0);
    count[0] = local_Nx*(rbd - lbd > 0);
    count[1] = my;
    
    /* Compute the offset by receving values of count[0] from other processes.               */
    /* This is necessary because the size is not guaranteed to be the same for all processes.*/
    MPI_Barrier(comm);
    MPI_Allgather(count, 1, MPI_UNSIGNED_LONG_LONG, count_0_gathered, 1, MPI_UNSIGNED_LONG_LONG, comm);
    
    offset[0] = 0;
    
    for (int p_index=0; p_index<mpi_rank; p_index++) {
    
        offset[0] += count_0_gathered[p_index];
    
    }
    
    offset[1] = 0;
    
    /* Compute global size of the filespace */
    dimsf[0] = 0;
    
    for (int p_index=0; p_index<mpi_size; p_index++) {
    
        dimsf[0] += count_0_gathered[p_index];
        
    }
    
    dimsf[1] = my;
    
    //printf("process %d, dimsf %lld, count %lld, offset %lld \n",mpi_rank,dimsf[0],count[0],offset[0]);
    
    /* heta */
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "heta_r", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);

    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile. */
    for (i=0; i<count[0]; i++) {
        
        for (j=0; j<my; j++) {
            
            temp1[my*i + j] = creal(htemp1[(Ny/2+1)*i + j]);
        
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
   
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
//    // hphi //
//    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
//    h5_dataset = H5Dcreate(h5_file, "phi", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//    H5Sclose(h5_filespace);
//    
//    h5_memspace = H5Screate_simple(RANK, count, NULL);
//
//    h5_filespace = H5Dget_space(h5_dataset);
//    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
//    
//    plist_id = H5Pcreate(H5P_DATASET_XFER);
//    
//    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
//    for (i=0; i<local_Nx; i++) {
//        
//        for (j=0; j<Ny; j++) {
//            
//            temp1[Ny*i + j] = phi[(Ny+2)*i + j];
//            
//        }
//        
//    }
//    
//    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
//       
//    H5Dclose(h5_dataset);
//    H5Sclose(h5_filespace);
//    H5Sclose(h5_memspace);
//    H5Pclose(plist_id);
    
    
}


void write_extra_2d(hid_t h5_file, double* array1, double* array2){
    

    hid_t       h5_dataset, h5_memspace, h5_filespace;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "Array1", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);

    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = array1[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
   
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "Array2", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = array2[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
       
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
}


