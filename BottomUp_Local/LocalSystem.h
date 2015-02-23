//
//  LocalSystem.h
//  BottomUp_Local
//
//  Created by Sean R Purser-Haskell on 2/15/15.
//  Copyright (c) 2015 Sean R Purser-Haskell. All rights reserved.
//

#ifndef BottomUp_Local_LocalSystem_h
#define BottomUp_Local_LocalSystem_h

#include <vector>
#include <map>

#include "PointSearcher.h"

#include "Vec2f.h"
#include <cassert>

#include <OpenCL/opencl.h>

typedef unsigned ParticleTypeId;
extern const unsigned ParticleTypeId_Null;

struct ParticleType;

typedef std::map<ParticleTypeId, ParticleType> ParticleTypeById;

struct ParticleType
{
    // These are coefficients for a polynomial expression of d, the distance to another particle.
    // First one is linear, second quadratic, and so on
    // The final force between the particles is charge_a * charge_b * polynomial(d)
    const std::vector<std::tuple<float, float, bool> > coeff_power;
    const float charge;
    const float mass;
    const ParticleTypeId id;
    inline ParticleType(ParticleTypeId id, std::vector<std::tuple<float, float, bool> > const&coeff_power, float charge, float mass)
      : id(id), coeff_power(coeff_power), charge(charge), mass(mass) {
        assert(charge == -1.0f || charge == 1.0f);
    }
    float max_affected_distance(float accel_epsilon)const;
};

struct ParticleState {
    ParticleTypeId typeId;
    Vec2f pos, vel;
    
    inline ParticleState(ParticleTypeId const& typeId, Vec2f const&pos, Vec2f const&vel)
      : typeId(typeId), pos(pos), vel(vel) {
    }
    
    void AffectBy(ParticleState const&o, float t, float min_d, ParticleTypeById const&types);
};

struct ParticleAffecter {
    virtual
    void Affect(ParticleState &particle, float time_slice)const=0;
};

struct ConstantDensityPointSearcher;

struct LocalSystem
{
    LocalSystem(std::vector<ParticleType> const&particle_types,
                std::unique_ptr<PointSearcher> &point_searcher,
                int n_interactions,
                float min_d,
                float time_epsilon,
                bool use_opencl,
                cl_device_id device);
    void AddParticle(ParticleState const&initial);
    void GetParticles(std::vector<ParticleState> &output)const;
    
    void Iterate(float time, std::vector<ParticleAffecter const*> const&extra_affecters);
private:
    
    ParticleTypeById types_by_id_;
    std::vector<ParticleState> particles_;
    std::unique_ptr<PointSearcher> point_searcher_;
    const float min_d_, time_epsilon_;
    float time_remainder_;
    const int n_interactions_;

    const bool use_opencl_;
    // OpenCL
    ConstantDensityPointSearcher *searcher_constant_;
    cl_device_id cl_device_;
    cl_context cl_context_;
    cl_command_queue cl_commands_;
    cl_program cl_program_;
    cl_kernel cl_kernel_;
    void MakeIndicesSquare(std::vector<Extrema1i> &ranges,
                           Vec2i &size_out);
    
    Vec2i CellIndex(Vec2f const&position)const;
};

#endif
