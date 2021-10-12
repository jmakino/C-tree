#pragma once
//
// utils.hpp 
//
#include <unistd.h>
#include <time.h>
#include <errno.h>
namespace NbodyLib{
#ifdef UTILBASE    
    double t0;
#else
    extern double t0;
#endif    
    
    
    inline double SetWtimeBase()
    {
	struct timespec ts;
	auto  iret=clock_gettime(CLOCK_MONOTONIC,&ts );
	if (iret){
	    printf("SetWtimeBase Failed with %d %d \n", iret, errno);
	}
	t0= ((double) ts.tv_sec)+ ts.tv_nsec*1e-9;
	return t0;
    }
    inline double GetWtime()
    {
	struct timespec ts;
	auto  iret=clock_gettime(CLOCK_MONOTONIC,&ts );
	if (iret){
	    printf("GetWtime Failed with %d %d \n", iret, errno);
	}
	return ((double) ts.tv_sec)+ ts.tv_nsec*1e-9 - t0;
    }
}
