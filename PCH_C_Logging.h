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

// This is Mike Ash's method, from https://www.mikeash.com/pyblog/friday-qa-2009-08-21-writing-vararg-macros-and-functions.html
#define DLog(fmt, ...) do { \
if (DEBUG) \
fprintf(stderr, "%s:%d:%s " fmt "\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, ## __VA_ARGS__); \
} while(0)


#define ALog(fmt, ...) do { \
if (DEBUG) \
{ fprintf(stderr, fmt " ", ## __VA_ARGS__); assert(0); } \
else \
{ fprintf(stderr, "%s:%d:%s " fmt "\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, ## __VA_ARGS__); } \
} while (0)



#endif /* PCH_C_Logging_h */
