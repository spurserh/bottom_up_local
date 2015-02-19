//
//  PointSearcher.h
//  BottomUp_Local
//
//  Created by Sean R Purser-Haskell on 2/17/15.
//  Copyright (c) 2015 Sean R Purser-Haskell. All rights reserved.
//

#ifndef BottomUp_Local_PointSearcher_h
#define BottomUp_Local_PointSearcher_h

#include <vector>

#include "Vec2f.h"

struct PointSearcher {
    virtual
    void Build(std::vector<Vec2f> const&points)=0;
    
    virtual
    void SearchNearestK(Vec2f const&pt, size_t *output, int k)const=0;
};

#endif
