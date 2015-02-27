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
    dbg_max_points_considered = 0;
    
    if(points.size() > 0) {
        Extrema2i cell_extents(Vec2i(INT_MAX, INT_MAX), Vec2i(INT_MIN, INT_MIN));
        for(size_t pi = 0;pi < points.size();++pi) {
            Vec2f const&p = points[pi];
            const Vec2i cell = CellForPoint(p);
            cell_extents.DoEnclose(cell);
            cell_extents.DoEnclose(cell + Vec2i(1,1));
        }
        
        vector<vector<int> > points_by_cell;
        points_by_cell.resize(cell_extents.GetSize().width * cell_extents.GetSize().height);
        
        for(int pi = 0;pi < points.size();++pi) {
            Vec2f const&p = points[pi];
            const Vec2i cell = CellForPoint(p);
            const Vec2i array_cell = cell - cell_extents.mMin;
            points_by_cell[array_cell.y * cell_extents.GetSize().width + array_cell.x].push_back(pi);
        }
        
        for(int cell_row = cell_extents.mMin.y;cell_row < cell_extents.mMax.y;++cell_row) {
            for(int cell_col = cell_extents.mMin.x;cell_col < cell_extents.mMax.x;++cell_col) {
                const Vec2i cell(cell_col, cell_row);
                Extrema1i cell_indices(Vec1i((int)ref_points_.size()), Vec1i(-1));
                for(int row = -1;row <= 1;++row) {
                    for(int col = -1;col <= 1;++col) {
                        const Vec2i this_cell(cell + Vec2i(col, row));
                        const Vec2i this_array_cell = this_cell - cell_extents.mMin;
                        if(this_array_cell.x >= 0 && this_array_cell.y >= 0 &&
                           this_array_cell.x < cell_extents.GetSize().width &&
                           this_array_cell.y < cell_extents.GetSize().height) {
                            for(const int o_pt_idx : points_by_cell[this_array_cell.y * cell_extents.GetSize().width + this_array_cell.x]) {
                                ref_points_.push_back(o_pt_idx);
                            }
                        }
                    }
                }
                cell_indices.mMax[0] = (int)ref_points_.size();
                point_indices_to_consider_for_cell_.insert(
                    map<Vec2i, Extrema1i>::value_type(cell, cell_indices));
                dbg_max_points_considered = std::max(dbg_max_points_considered, cell_indices.mMax[0] - cell_indices.mMin[0]);
            }
        }
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
        priority_queue.push(ref_points_[i]);
    }

    for(size_t i=0;!priority_queue.empty() && i < k;++i) {
        output[i] = priority_queue.top();
        priority_queue.pop();
    }
}

Vec2i ConstantDensityPointSearcher::CellForPoint(Vec2f const&p)const {
    return Vec2i(int(::floor(p.x / cell_width_)), int(::floor(p.y / cell_width_)));
}


