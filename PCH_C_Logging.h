//
//  PCH_C_Logging.h
//  BlueBookFunctions
//
//  Created by Peter Huber on 2018-01-12.
//  Copyright © 2018 Peter Huber. All rights reserved.
//

#ifndef PCH_C_Logging_h
#define PCH_C_Logging_h

#include <assert.h>

#if !defined(__PRETTY_FUNCTION__)

    #if defined(__FUNCTION__)

        #define __PRETTY_FUNCTION__ __FUNCTION__

    #else

        #define __PRETTY_FUNCTION__ __func__

    #endif

#endif

#define DLog(...) do { if (DEBUG) fprintf(stderr, "%s:%d:%s: ", __FILE__, __LINE__, __PRETTY_FUNCTION__, __VA_ARGS__); } while (0)

#define ALog(fmt, ...) do { if (DEBUG) { fprintf(stderr, "%s", fmt) ; assert(0); } else { fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __PRETTY_FUNCTION__, __VA_ARGS__)} } while (0)

void test(char *fmt)
{
    do { if (DEBUG) fprintf(stderr, "%s: %s:%d:%s(): ", fmt, __FILE__, __LINE__, __PRETTY_FUNCTION__, "va_args"); } while (0);
    
    if (DEBUG) { fprintf(stderr, "%s", fmt) ; assert(0); } else { fprintf(stderr, "%s:%d:%s(): ", fmt, __FILE__, __LINE__, __PRETTY_FUNCTION__, __VA_ARGS__)}
}


#endif /* PCH_C_Logging_h */
