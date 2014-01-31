//
//  time_schemes.h
//  euler
//
//  Created by Claudio Viotti on 3/9/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#define DIM 2

#if DIM == 1
    #include "euler.h"
#elif DIM == 2
    #include "euler_2d.h"
#endif

#if DIM == 1
    #include "model.h"
#elif DIM == 2
    #include "model_2d.h"
#endif



#ifndef euler_time_schemes_h
    #define euler_time_schemes_h
#endif


void sol_update_RK(fftw_complex* u_old,double* t,double dt,char* dtflag);
void Setup_TimeScheme(int scheme_flg);