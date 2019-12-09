//
//  PCH_BasicStructs.hpp
//  BlueBookFunctions
//
//  Created by Peter Huber on 2019-12-08.
//  Copyright © 2019 Peter Huber. All rights reserved.
//

#ifndef PCH_BasicStructs_hpp
#define PCH_BasicStructs_hpp

#include <stdio.h>

typedef struct pch_Point {
    
    double x;
    double y;
    
} PCH_Point;

struct PCH_Size:PCH_Point {
    
    double width()
    {
        return this->x;
    }
    
    double height()
    {
        return this->y;
    }
    
};

typedef struct pch_Rect {
    
    PCH_Point origin;
    PCH_Size size;
    
    double Area()
    {
        return this->size.height() * this->size.width();
    }
    
} PCH_Rect;

#endif /* PCH_BasicStructs_hpp */
