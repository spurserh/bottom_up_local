//
//  LocalSystem.cpp
//  BottomUp_Local
//
//  Created by Sean R Purser-Haskell on 2/15/15.
//  Copyright (c) 2015 Sean R Purser-Haskell. All rights reserved.
//

#include "LocalSystem.h"

using namespace std;

const unsigned ParticleTypeId_Null = 0;

float ParticleType::max_affected_distance(float accel_epsilon)const {
    const float force_epsilon = accel_epsilon / mass;
    float max_affected_dist = 0;
    for(std::tuple<float, float, bool> const& coeff_power : this->coeff_power) {
        const float coefficient = std::get<0>(coeff_power);
        const float power = std::get<1>(coeff_power);
        assert(power >= 1.0f);
//        const float this_max_affect_dist =  ::pow(force_epsilon / ::fabs(coefficient), 1.0f / power);
        const float this_max_affect_dist = ::pow(::fabs(coefficient), 1.0f / power) / ::pow(::fabs(force_epsilon), 1.0f / power);

        fprintf(stderr, "this_max_affect_dist(%f): %f\n", power, this_max_affect_dist);
        max_affected_dist = std::max(max_affected_dist, this_max_affect_dist);
    }
    return max_affected_dist;
}

void ParticleState::AffectBy(ParticleState const&o, float t, float min_d, ParticleTypeById const&types) {
    ParticleType const&type = types.find(typeId)->second;
    ParticleType const&o_type = types.find(o.typeId)->second;
    const float d = (pos - o.pos).Length();
    const float charge_prod = type.charge * o_type.charge;
    float f = 0.0f;
    for(std::tuple<float, float, bool> const& coeff_power : type.coeff_power) {
        const float coefficient = std::get<0>(coeff_power);
        const float power = std::get<1>(coeff_power);
        assert(power >= 1.0f);
        float charge_prod_term = std::get<2>(coeff_power) ? charge_prod : 1.0f;
        f += charge_prod_term * coefficient * (1.0f / ::pow(std::max(d, min_d), power));
    }
    const float accel = f / o_type.mass;
    const Vec2f accel_vec = (o.pos - pos).Normalized() * accel;
    
 //   if (::fabs(accel) >= accel_epsilon)
        vel += accel_vec * t;
}

LocalSystem::LocalSystem(std::vector<ParticleType> const&particle_types,
                         std::unique_ptr<PointSearcher> &point_searcher,
                         int n_interactions,
                         float min_d,
                         float time_epsilon)
  : min_d_(min_d),
    time_epsilon_(time_epsilon),
    n_interactions_(n_interactions),
    point_searcher_(std::move(point_searcher)),
    time_remainder_(0) {
    for(ParticleType const&type : particle_types) {
        types_by_id_.insert(ParticleTypeById::value_type(type.id, type));
    }
#if 0
    max_distance_ = 0;
    for(ParticleType const&type : particle_types) {
        types_by_id_.insert(ParticleTypeById::value_type(type.id, type));
        
        const float this_max_dist = type.max_affected_distance();
        max_distance_ = std::max(max_distance_, this_max_dist);
    }
        fprintf(stderr, "max distance %f\n", max_distance_);
#endif
}

void LocalSystem::AddParticle(ParticleState const&initial) {
    particles_.push_back(initial);
}

void LocalSystem::GetParticles(std::vector<ParticleState> &output)const {
    std::copy(particles_.begin(), particles_.end(), std::back_inserter(output));
}

void LocalSystem::Iterate(float time, std::vector<ParticleAffecter const*> const&extra_affecters) {
    // Closest is always self
    const size_t n_interact = 1 + n_interactions_;
    
    time_remainder_ += time;
    
    while(time_remainder_ > time_epsilon_) {
        vector<Vec2f> points;
        for(ParticleState const&particle : particles_) {
            points.push_back(particle.pos);
        }
        point_searcher_->Build(points);
        
        const float time_slice = std::min(time_remainder_, time_epsilon_);
        time_remainder_ -= time_slice;
        vector<ParticleState> next_particles;
        vector<size_t> nearest_k;
        
        float total_velocity = 0.0f, next_total_velocity = 0.0f;;
        
        // TODO: Parallelize, use GPU?
        for(ParticleState const&a : particles_) {
            ParticleState next_particle(a);
            
            nearest_k.resize(n_interact);
            point_searcher_->SearchNearestK(a.pos, &nearest_k[0], n_interact);
            
            for(int i=0;i<n_interact;++i) {
                ParticleState const&b = particles_[nearest_k[i]];
                // Closest is always self..
                if(&a == &b)
                    continue;
                // Pos is not changed, so this is okay. Probably should be more explicit.
                next_particle.AffectBy(b, time_slice, min_d_, types_by_id_);
            }
            
            for(ParticleAffecter const*affecter : extra_affecters) {
                affecter->Affect(next_particle, time_slice);
            }
            next_particle.pos += next_particle.vel * time_slice;
            next_particles.push_back(next_particle);
            total_velocity += a.vel.Length();
            next_total_velocity += next_particle.vel.Length();
        }
        // Normalize total energy
//        for(ParticleState &state : next_particles)
  //          state.vel *= total_velocity / next_total_velocity;
        particles_ = next_particles;
    }
}
#if 0
Vec2i LocalSystem::CellIndex(Vec2f const&position)const {
    return Vec2i(int(position.x / max_distance_), int(position.y / max_distance_));
}
#endif
