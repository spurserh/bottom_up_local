//
//  LocalSystem.cpp
//  BottomUp_Local
//
//  Created by Sean R Purser-Haskell on 2/15/15.
//  Copyright (c) 2015 Sean R Purser-Haskell. All rights reserved.
//

#include "LocalSystem.h"
#include "ConstantDensityPointSearcher.h"
#include <OpenCL/opencl.h>

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
    const float accel = f / type.mass;
    const Vec2f accel_vec = (o.pos - pos).Normalized() * accel;
    
 //   if (::fabs(accel) >= accel_epsilon)
        vel += accel_vec * t;
}

LocalSystem::LocalSystem(std::vector<ParticleType> const&particle_types,
                         std::unique_ptr<PointSearcher> &point_searcher,
                         int n_interactions,
                         float min_d,
                         float time_epsilon,
                         bool use_opencl,
                         cl_device_id device)
  : min_d_(min_d),
    time_epsilon_(time_epsilon),
    n_interactions_(n_interactions),
    point_searcher_(std::move(point_searcher)),
    time_remainder_(0),
    use_opencl_(use_opencl),
    cl_device_(device),
    cl_context_(0),
    cl_commands_(0),
    cl_program_(0),
    cl_kernel_(0) {
    searcher_constant_ = dynamic_cast<ConstantDensityPointSearcher*>(point_searcher_.get());
    assert(!use_opencl || searcher_constant_);
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
        
    if(use_opencl_) {
        int err;
        cl_context_ = clCreateContext(0, 1, &cl_device_, NULL, NULL, &err);
        assert(cl_context_);
        
        // Create a command commands
        //
        cl_commands_ = clCreateCommandQueue(cl_context_, cl_device_, 0, &err);
        assert(cl_commands_);
        
        // TODO
        const char* KernelSource =
            R"(__kernel void affect(
                __global const float2* pos,
                __global const int2* cell_affected_ranges,
                __global const unsigned int* cell_affected_indices,
                __global float2* vel_offset_out,
                const float t,
                const unsigned int cell_ranges_width,
                const unsigned int particle_count) {
                    const float cell_width = %ff;
                    const float min_d = %ff;
                    const unsigned int particle_index = get_global_id(0);
                    if(particle_index < particle_count) {
                        const float2 a_pos = pos[particle_index];
                        const int2 cell_index = (int2)(floor(a_pos.x / cell_width), floor(a_pos.y / cell_width));
                        const int2 cell_range = cell_affected_ranges[cell_index.y * cell_ranges_width + cell_index.x];
                        float2 vel_offset = (float2)(0.0f,0.0f);
                        /*
                        for(int affected_index = cell_range.x;affected_index < cell_range.y;++affected_index) {
                            // TODO: Avoid interacting with self
                            // TODO: Optimization? Avoid using another layer of indices
                          //  const float2 affected_pos = pos[cell_affected_indices[affected_index]];
                          //  vel_offset += normalize(a_pos - affected_pos) * t * 0.1f;
                            vel_offset += (float2)(-0.2f, 0.0f);
                        }
                         */
                        //vel_offset += (float2)(particle_index / 2000.0f, 0.0f);
                        //vel_offset += (float2)(cell_index.x, cell_index.y);
                        //vel_offset += normalize(a_pos) * t * 0.1f;
                        //vel_offset_out[particle_index] = vel_offset;
                        vel_offset_out[particle_index] = (float2)(cell_index.x, cell_index.y);
                    }
                })";
        
        vector<char> kernel_subbed;
        kernel_subbed.resize(strlen(KernelSource) + 1024);
        sprintf(&kernel_subbed[0], KernelSource, searcher_constant_->cell_width_, min_d);
        char const*const kernel_subbed_ptr = &kernel_subbed[0];
        
        fprintf(stderr, "kernel %s\n", kernel_subbed_ptr);
        
        // Create the compute program from the source buffer
        //
        cl_program_ = clCreateProgramWithSource(cl_context_, 1, (const char **)&kernel_subbed_ptr, NULL, &err);
        if (!cl_program_)
        {
            printf("Error: Failed to create compute program!\n");
            exit(1);
        }
        
        // Build the program executable
        //
        err = clBuildProgram(cl_program_, 0, NULL, NULL, NULL, NULL);
        size_t len;
        char buffer[2048];
        clGetProgramBuildInfo(cl_program_, cl_device_, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);

        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to build program executable!\n");
            exit(1);
        }
        
        // Create the compute kernel in the program we wish to run
        //
        cl_kernel_ = clCreateKernel(cl_program_, "affect", &err);
        if (!cl_kernel_ || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel!\n");
            exit(1);
        }
    }
}

void LocalSystem::AddParticle(ParticleState const&initial) {
    particles_.push_back(initial);
}

void LocalSystem::GetParticles(std::vector<ParticleState> &output)const {
    std::copy(particles_.begin(), particles_.end(), std::back_inserter(output));
}

void LocalSystem::MakeIndicesSquare(std::vector<Extrema1i> &ranges,
                                    Vec2i &size_out) {
    Extrema2i cell_extents(Vec2i(INT_MAX, INT_MAX), Vec2i(-INT_MAX, -INT_MAX));
    for(auto it : searcher_constant_->point_indices_to_consider_for_cell_) {
        cell_extents.DoEnclose(it.first);
        cell_extents.DoEnclose(it.first + Vec2i(1,1));
    }
    
    size_out = cell_extents.GetSize();
    ranges.resize(size_out.width * size_out.height, Extrema1i(Vec1i(0), Vec1i(0)));
    for(int row = cell_extents.mMin.y;row < cell_extents.mMax.y;++row) {
        for(int col = cell_extents.mMin.x;col < cell_extents.mMax.x;++col) {
            auto found = searcher_constant_->point_indices_to_consider_for_cell_.find(Vec2i(col, row));
            if(found != searcher_constant_->point_indices_to_consider_for_cell_.end()) {
                const int row_r = row - cell_extents.mMin.y;
                const int col_r = col - cell_extents.mMin.x;
                assert((row_r*size_out.width+col_r) < ranges.size());
                ranges[row_r*size_out.width+col_r] = found->second;
            }
        }
    }
}

void LocalSystem::Iterate(float time, std::vector<ParticleAffecter const*> const&extra_affecters) {
    // Closest is always self
    const size_t n_interact = 1 + n_interactions_;
    
    time_remainder_ += time;
    
    while(time_remainder_ > time_epsilon_) {
        const float time_slice = std::min(time_remainder_, time_epsilon_);
        time_remainder_ -= time_slice;

        vector<Vec2f> points;
        for(ParticleState const&particle : particles_) {
            points.push_back(particle.pos);
        }
        point_searcher_->Build(points);

        vector<ParticleState> next_particles = particles_;
        if(!use_opencl_) {
            vector<size_t> nearest_k;
            for(size_t particle_index = 0;particle_index < particles_.size();++particle_index) {
                ParticleState const&a = particles_[particle_index];
                ParticleState &next_particle = next_particles[particle_index];
                
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
            }
        } else {
            unsigned int count = (unsigned int)particles_.size();
            cl_mem cl_pos, cl_cell_affected_ranges, cl_cell_affected_indices, cl_vel_offset_out;
            cl_pos = clCreateBuffer(cl_context_,  CL_MEM_READ_ONLY,  sizeof(float) * 2 * particles_.size(), NULL, NULL);
            cl_vel_offset_out = clCreateBuffer(cl_context_,  CL_MEM_WRITE_ONLY,  sizeof(float) * 2 * particles_.size(), NULL, NULL);
            vector<Vec2f> pos_data;
            pos_data.resize(particles_.size());
            
            for(size_t i=0;i<particles_.size();++i)
                pos_data[i] = particles_[i].pos;
            
            vector<Extrema1i> ranges;
            Vec2i ranges_size;
            
            MakeIndicesSquare(ranges, ranges_size);
            
            int err = clEnqueueWriteBuffer(cl_commands_, cl_pos, CL_TRUE, 0, sizeof(float) * 2 * particles_.size(), &pos_data[0], 0, NULL, NULL);
            assert(err == CL_SUCCESS);
            cl_cell_affected_indices = clCreateBuffer(cl_context_,  CL_MEM_READ_ONLY,  sizeof(int) * searcher_constant_->ref_points_.size(), NULL, NULL);
            cl_cell_affected_ranges = clCreateBuffer(cl_context_,  CL_MEM_READ_ONLY,  sizeof(int) * 2 * ranges.size(), NULL, NULL);
            err = clEnqueueWriteBuffer(cl_commands_, cl_cell_affected_indices, CL_TRUE, 0, sizeof(int) * searcher_constant_->ref_points_.size(), &searcher_constant_->ref_points_[0], 0, NULL, NULL);
            assert(err == CL_SUCCESS);
            err = clEnqueueWriteBuffer(cl_commands_, cl_cell_affected_ranges, CL_TRUE, 0, sizeof(int) * 2 * ranges.size(), &ranges[0], 0, NULL, NULL);
            assert(err == CL_SUCCESS);
            
            err = 0;
            err  = clSetKernelArg(cl_kernel_, 0, sizeof(cl_mem), &cl_pos);
            err |= clSetKernelArg(cl_kernel_, 1, sizeof(cl_mem), &cl_cell_affected_ranges);
            err |= clSetKernelArg(cl_kernel_, 2, sizeof(cl_mem), &cl_cell_affected_indices);
            err |= clSetKernelArg(cl_kernel_, 3, sizeof(cl_mem), &cl_vel_offset_out);
            err |= clSetKernelArg(cl_kernel_, 4, sizeof(float), &time_slice);
            err |= clSetKernelArg(cl_kernel_, 5, sizeof(unsigned int), &ranges_size.width);
            err |= clSetKernelArg(cl_kernel_, 6, sizeof(unsigned int), &count);
            assert(err == CL_SUCCESS);

            size_t global;                      // global domain size for our calculation
            size_t local;                       // local domain size for our calculation
            
            err = clGetKernelWorkGroupInfo(cl_kernel_, cl_device_, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
            assert(err == CL_SUCCESS);
            
            local = 16;
            
            // Execute the kernel over the entire range of our 1d input data set
            // using the maximum number of work group items for this device
            //
            global = count;
            err = clEnqueueNDRangeKernel(cl_commands_, cl_kernel_, 1, NULL, &global, &local, 0, NULL, NULL);
            assert(err == CL_SUCCESS);
            
            clFinish(cl_commands_);
            
            vector<Vec2f> vel_offset_data;
            vel_offset_data.resize(count);
            
            // Blocking
            err = clEnqueueReadBuffer(cl_commands_, cl_vel_offset_out, CL_TRUE, 0, sizeof(float) * 2 * count, &vel_offset_data[0], 0, NULL, NULL );
            assert(err == CL_SUCCESS);

            clReleaseMemObject(cl_pos);
            clReleaseMemObject(cl_cell_affected_ranges);
            clReleaseMemObject(cl_cell_affected_indices);
            clReleaseMemObject(cl_vel_offset_out);
            
            for(size_t particle_index = 0;particle_index < particles_.size();++particle_index) {
                ParticleState &next_particle = next_particles[particle_index];
                next_particle.vel += vel_offset_data[particle_index];
            }
        }
        for(size_t particle_index = 0;particle_index < particles_.size();++particle_index) {
            ParticleState &next_particle = next_particles[particle_index];
            for(ParticleAffecter const*affecter : extra_affecters) {
                affecter->Affect(next_particle, time_slice);
            }
            next_particle.pos += next_particle.vel * time_slice;
        }

        particles_ = next_particles;
    }
}
#if 0
Vec2i LocalSystem::CellIndex(Vec2f const&position)const {
    return Vec2i(int(position.x / max_distance_), int(position.y / max_distance_));
}
#endif
