//
//  hdf5_routines_2d.c
//  euler
//
//  Created by Claudio Viotti on 6/24/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "hdf5_routines_2d.h"


void get_params(char* filename){
    
    hid_t       h5_file, h5_data;
    herr_t      status;
    double      Nd[1];
    
    
    h5_file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    
    h5_data = H5Dopen(h5_file, "/Nx", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    Nx = (int) Nd[0];
    
    h5_data = H5Dopen(h5_file, "/Ny", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    Ny = (int) Nd[0];
    
    h5_data = H5Dopen(h5_file, "/Lx", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lx);
    h5_data = H5Dopen(h5_file, "/Ly", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Ly);
    
    h5_data = H5Dopen(h5_file, "/runsubid", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);

    runsubid = (int) Nd[0];

    h5_data = H5Dopen(h5_file, "/T", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
    
    h5_data = H5Dopen(h5_file, "/dtsave", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dtsave);
    
    status = H5Dclose(h5_data);
    status = H5Fclose(h5_file);
    
    Kx0 = 2*M_PI/Lx;
    Ky0 = 2*M_PI/Ly;
    
}


void get_ic_2d(char* filename, double* u1, double* u2){
    
    hid_t       h5_file, h5_data;
    herr_t      status;
    
    
    h5_file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    
    h5_data = H5Dopen(h5_file, "/eta0", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,u1);
    
    h5_data = H5Dopen(h5_file, "/phi0", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,u2);
    
    status = H5Dclose(h5_data);
    status = H5Fclose(h5_file);
    
    
}

hid_t create_file_2d(char* filename){
    
    hid_t       file; //, data, fid1;
   
    
    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    return file;
    
}


hid_t open_file_2d(char* filename){
    
    hid_t       file;
    
    
    file = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT );
    
    return file;
    
}


herr_t close_file_2d(hid_t file){
    
    herr_t      status;
    
    
    status = H5Fclose(file);
    
    return status;
    
}


void write_header_2d(hid_t file, double t){
    
    hid_t       data, fid;
    herr_t      status;
    hsize_t     fdim [2];
    double      var;
    
    fdim[0]=1;
    fdim[1]=1;
    
    fid = H5Screate_simple(2, fdim, NULL);
    
    data = H5Dcreate(file, "time", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&t);
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
    
}


void write_field_2d(hid_t file, double* eta, double* phi){
    
    hid_t       data, fid1, fid2;
    herr_t      status;
    hsize_t     fdim [2], fdim2 [2];
    char        data_name_eta[10], data_name_phi[10];
    
    fdim[0]=Nx;
    fdim[1]=Ny;
    
    fdim2[0]=1;
    fdim2[1]=1;
    
    fid1 = H5Screate_simple(2, fdim, NULL);
    fid2 = H5Screate_simple(2, fdim2, NULL);
    
    sprintf(data_name_eta, "eta");
    sprintf(data_name_phi, "phi");
    
    
    data = H5Dcreate(file, data_name_eta, H5T_IEEE_F64LE, fid1, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,eta);
    status = H5Dclose(data);
    
    data = H5Dcreate(file, data_name_phi, H5T_IEEE_F64LE, fid1, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,phi);
    status = H5Dclose(data);
    
    
}



void write_extra_2d(hid_t file, double* Array1, double* Array2){
    
    hid_t       data, fid1, fid2;
    herr_t      status;
    hsize_t     fdim [2], fdim2 [2];
    char        data_name_array1[12], data_name_array2[12];
    
    fdim[0]=Nx;
    fdim[1]=Ny;
    
    fdim2[0]=1;
    fdim2[1]=1;
    
    fid1 = H5Screate_simple(2, fdim, NULL);
    fid2 = H5Screate_simple(2, fdim2, NULL);
    
    sprintf(data_name_array1, "Array1");
    sprintf(data_name_array2, "Array2");
    
    
    data = H5Dcreate(file, data_name_array1, H5T_IEEE_F64LE, fid1, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,Array1);
    status = H5Dclose(data);
    
    data = H5Dcreate(file, data_name_array2, H5T_IEEE_F64LE, fid1, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,Array2);
    status = H5Dclose(data);
    
    
}


