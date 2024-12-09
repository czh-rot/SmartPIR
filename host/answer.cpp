#include "answer.hpp"
//#include <define_host.h>

unsigned long diff(const struct timeval *newTime, const struct timeval *oldTime)
{
    return (newTime->tv_sec - oldTime->tv_sec) * 1000000 + (newTime->tv_usec - oldTime->tv_usec);
}
