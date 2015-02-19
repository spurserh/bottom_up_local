//
//  ConstantDensityPointSearcher.h
//  BottomUp_Local
//
//  Created by Sean R Purser-Haskell on 2/17/15.
//  Copyright (c) 2015 Sean R Purser-Haskell. All rights reserved.
//

#ifndef __BottomUp_Local__ConstantDensityPointSearcher__
#define __BottomUp_Local__ConstantDensityPointSearcher__

#include "PointSearcher.h"

#include <memory>
#include <map>

struct LocalSystem;

struct ConstantDensityPointSearcher : public PointSearcher {
    ConstantDensityPointSearcher(float cell_width);
    void Build(std::vector<Vec2f> const&points)override;
    // Point indices output
    // TODO: Randomize a bit
    // Does not search outside of AABB of original points
    void SearchNearestK(Vec2f const&pt, size_t *output, int k)const override;
private:
  friend class LocalSystem;
    const float cell_width_;
    int dbg_max_points_considered;
    std::map<Vec2i, Extrema1i> point_indices_to_consider_for_cell_;
    std::vector<Vec2f> ref_points_;
    std::vector<Vec2f> original_points_;
    
    // Temp
    std::map<Vec2f, size_t> index_by_location_;

    Vec2i CellForPoint(Vec2f const&p)const;
};

#endif /* defined(__BottomUp_Local__ConstantDensityPointSearcher__) */
