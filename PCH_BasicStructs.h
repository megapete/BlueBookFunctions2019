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
    
    float x;
    float y;
    
    pch_Point() {this->x = 0.0; this->y = 0.0;}
    pch_Point(double x, double y) {this->x = x; this->y = y;}
    pch_Point(const pch_Point &source) {this->x = source.x; this->y = source.y;}
    
} PCH_Point;

typedef struct pch_Size {
    
    float width;
    float height;
    
    pch_Size() {this->width = 0.0; this->height = 0.0;}
    pch_Size(double w, double h) {this->width = w; this->height = h;}
    pch_Size(const pch_Size &source) {this->width = source.width; this->height = source.height;}
    
} PCH_Size;

typedef struct pch_Rect {
    
    PCH_Point origin;
    PCH_Size size;
    
    pch_Rect() {this->origin = PCH_Point(); this->size = PCH_Size();}
    pch_Rect(PCH_Point origin, PCH_Size size) {this->origin = PCH_Point(origin); this->size = PCH_Size(size);}
    pch_Rect(const pch_Rect &source) {this->origin = PCH_Point(source.origin); this->size = PCH_Size(source.size);}
    
    float Area()
    {
        return this->size.height * this->size.width;
    }
    
} PCH_Rect;

#endif /* PCH_BasicStructs_hpp */
