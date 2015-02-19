//
//  ConstantDensityPointSearcher.cpp
//  BottomUp_Local
//
//  Created by Sean R Purser-Haskell on 2/17/15.
//  Copyright (c) 2015 Sean R Purser-Haskell. All rights reserved.
//

#include "ConstantDensityPointSearcher.h"
#include <set>
#include <map>
#include <queue>

using namespace std;

ConstantDensityPointSearcher::ConstantDensityPointSearcher(float cell_width)
  : cell_width_(cell_width), dbg_max_points_considered(0) {
}

void ConstantDensityPointSearcher::Build(std::vector<Vec2f> const&points) {
    point_indices_to_consider_for_cell_.clear();
    ref_points_.clear();
    original_points_.clear();
    index_by_location_.clear();
    dbg_max_points_considered = 0;
    
    set<Vec2i> cells;
    multimap<Vec2i, Vec2f> points_by_cell;
    for(Vec2f const&p : points) {
        const Vec2i cell = CellForPoint(p);
        cells.insert(cell);
        points_by_cell.insert(multimap<Vec2i, Vec2f>::value_type(cell, p));
    }
    
    for(size_t pi=0;pi<points.size();++pi)
        index_by_location_[points[pi]] = pi;
    
    for(Vec2i const&cell : cells) {
        Extrema1i cell_indices(Vec1i((int)ref_points_.size()), Vec1i(-1));
        for(int row = -1;row <= 1;++row) {
            for(int col = -1;col <= 1;++col) {
                const Vec2i this_cell(cell + Vec2i(col, row));
                for(multimap<Vec2i, Vec2f>::const_iterator it = points_by_cell.lower_bound(this_cell);
                    it != points_by_cell.upper_bound(this_cell);
                    ++it) {
                    ref_points_.push_back(it->second);
                }
            }
        }
        cell_indices.mMax[0] = (int)ref_points_.size();
        point_indices_to_consider_for_cell_.insert(
            map<Vec2i, Extrema1i>::value_type(cell, cell_indices));
        dbg_max_points_considered = std::max(dbg_max_points_considered, cell_indices.mMax[0] - cell_indices.mMin[0]);
    }
    
    fprintf(stderr, "dbg_max_points_considered %i\n", dbg_max_points_considered);
    original_points_ = points;
}

struct CloserToPt {
    CloserToPt(Vec2f const&pt, std::vector<Vec2f> const&points)
      : pt(pt), points_(points) {
    }
    
    bool operator()(size_t a, size_t b)const {
        return (points_[a]-pt).Length() > (points_[b]-pt).Length();
    }
    
    const Vec2f pt;
    std::vector<Vec2f> const&points_;
};

// Point indices output
// TODO: Randomize a bit
void ConstantDensityPointSearcher::SearchNearestK(Vec2f const&pt, size_t *output, int k)const {
    const CloserToPt compare(pt, original_points_);
    priority_queue<size_t, vector<size_t>, CloserToPt> priority_queue(compare);
    map<Vec2i, Extrema1i>::const_iterator it = point_indices_to_consider_for_cell_.find(CellForPoint(pt));
    if(it == point_indices_to_consider_for_cell_.end()) {
        fprintf(stderr, "Warning: ConstantDensityPointSearcher::SearchNearestK() outside AABB");
        return;
    }
    for(int i = it->second.mMin[0];i < it->second.mMax[0];++i) {
        priority_queue.push(index_by_location_.find(ref_points_[i])->second);
    }

    for(size_t i=0;!priority_queue.empty() && i < k;++i) {
        output[i] = priority_queue.top();
        priority_queue.pop();
    }
}

Vec2i ConstantDensityPointSearcher::CellForPoint(Vec2f const&p)const {
    return Vec2i(int(p.x / cell_width_), int(p.y / cell_width_));
}


