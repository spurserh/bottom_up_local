//
//  KDTreePointSearcher.cpp
//  BottomUp_Local
//
//  Created by Sean R Purser-Haskell on 2/17/15.
//  Copyright (c) 2015 Sean R Purser-Haskell. All rights reserved.
//

#include "KDTreePointSearcher.h"

void KDTreePointSearcher::Build(std::vector<Vec2f> const&points_in) {
    multi_array<float, 2> points;
    points.resize(extents[points_in.size()][2]);
    for(size_t i=0;i<points_in.size();++i) {
        points[i][0] = points_in[i].x;
        points[i][1] = points_in[i].y;
    }
    kd_tree_.reset(new kdtree2(points));
}

void KDTreePointSearcher::SearchNearestK(Vec2f const&pt, size_t *output, int k)const {
    kdtree2_result_vector nearest_k;
    vector<float> qv = {pt.x, pt.y};
    kd_tree_->n_nearest(qv, k, nearest_k);
    for(int i=0;i<k;++i)
        output[i] = nearest_k[i].idx;
}