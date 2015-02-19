//
//  KDTreePointSearcher.h
//  BottomUp_Local
//
//  Created by Sean R Purser-Haskell on 2/17/15.
//  Copyright (c) 2015 Sean R Purser-Haskell. All rights reserved.
//

#ifndef BottomUp_Local_KDTreePointSearcher_h
#define BottomUp_Local_KDTreePointSearcher_h

#include "PointSearcher.h"
#include "kdtree2.hpp"

#include <memory>

struct KDTreePointSearcher : public PointSearcher {
    void Build(std::vector<Vec2f> const&points)override;
    // Point indices output
    void SearchNearestK(Vec2f const&pt, size_t *output, int k)const override;
  private:
    std::unique_ptr<kdtree2> kd_tree_;
};

#endif
