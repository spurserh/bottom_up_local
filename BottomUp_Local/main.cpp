

#include <memory>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <GLUT/glut.h>
#include <cassert>
#include <cstdio>
#include <vector>

#include "Vec2f.h"
#include "LocalSystem.h"
#include "KDTreePointSearcher.h"
#include "ConstantDensityPointSearcher.h"

#include <OpenCL/opencl.h>

#include "kdtree2.hpp"

using namespace std;

const float time_epsilon = 0.005f;
const float min_d = 0.001f;
const int n_interactions = 8;

std::unique_ptr<LocalSystem> local_system;
struct BoxAffecter : public ParticleAffecter {
    BoxAffecter(Extrema2f const&bounds)
      : bounds(bounds) {
        
    }
    
    void Affect(ParticleState &particle, float time_slice)const override {
        // pt_a, pt_b
        static const tuple<Vec2f, Vec2f> boundary_lines[] = {
            {bounds.mMin + bounds.GetSize() * Vec2f(0,0), bounds.mMin + bounds.GetSize() * Vec2f(1,0)},
            {bounds.mMin + bounds.GetSize() * Vec2f(1,0), bounds.mMin + bounds.GetSize() * Vec2f(1,1)},
            {bounds.mMin + bounds.GetSize() * Vec2f(1,1), bounds.mMin + bounds.GetSize() * Vec2f(0,1)},
            {bounds.mMin + bounds.GetSize() * Vec2f(0,1), bounds.mMin + bounds.GetSize() * Vec2f(0,0)},
        };
        for(auto const&line : boundary_lines) {
            const Vec2f o = std::get<0>(line);
            const Vec2f d = (std::get<1>(line) - std::get<0>(line)).Normalized();
            const float t = (particle.pos - o).Dot(d);
            const Vec2f line_pt = o + d * t;
            const float dist = (line_pt - particle.pos).Length();
            const Vec2f line_to_pt = (line_pt - particle.pos).Normalized();
            particle.vel -= line_to_pt * time_slice * (1.0f / ::pow(std::max(min_d, dist), 3)) * 0.000001f;
            // Cannot change pos..
            // TODO
//            particle.pos.x = std::max(bounds.mMin.x, std::min(bounds.mMax.x, particle.pos.x));
  //          particle.pos.y = std::max(bounds.mMin.y, std::min(bounds.mMax.y, particle.pos.y));
        }
    }
    
    const Extrema2f bounds;
    //}box_affecter(Extrema2f(Vec2f(-1, -1) * 0.3f, Vec2f(1, 1) * 0.3f));
}box_affecter(Extrema2f(Vec2f(-1, -1) * 0.15f, Vec2f(1, 1) * 0.15f));

// Kind of like wind resistance or something
struct TooFastAffecter : public ParticleAffecter {
    TooFastAffecter(float speed_limit)
    : speed_limit(speed_limit) {
        
    }
    
    void Affect(ParticleState &particle, float time_slice)const override {
        const float speed = particle.vel.Length();
        if(speed > speed_limit)
            particle.vel *= 1.0f / ::pow(speed / speed_limit, 2);
    }
    
    const float speed_limit;
}speed_limiter(0.1f);
std::vector<ParticleAffecter const*> affecters = {&box_affecter, &speed_limiter};

float runi() {
    return float(rand()) / float(RAND_MAX);
}

Vec2f rvel() {
    const float a = 2.0f * M_PI * runi();
    return Vec2f(cos(a), sin(a));
}

vector<Vec2f> RandomParticles(const size_t n, const float min_radius) {
    srand(500);
    const size_t build_every = 100;
    unique_ptr<kdtree2> kd_tree;
    vector<Vec2f> ret;
    while(ret.size() < n) {
        Vec2f new_pt(runi() * 2.0 - 1.0, runi() * 2.0 - 1.0);
        // See if it's too close to any existing points
        if(kd_tree) {
            vector<float> new_pt_v = {new_pt.x, new_pt.y};
            kdtree2_result_vector nearest_v;
            kd_tree->n_nearest(new_pt_v, 1, nearest_v);
            const Vec2f nearest_tree = ret[nearest_v[0].idx];
            if((nearest_tree - new_pt).Length() < min_radius)
                continue;
        }
        if(ret.size() > 0) {
            bool invalid_point = false;
            for(int i = ret.size() % build_every;i > 0;--i) {
                const Vec2f nearest_remainder = ret[ret.size() - i];
                if((nearest_remainder - new_pt).Length() < min_radius) {
                    invalid_point = true;
                    break;
                }
            }
            if(invalid_point)
                continue;
        }
        ret.push_back(new_pt);
        if(ret.size() > 0 && ret.size() % build_every == 0) {
            multi_array<float, 2> points;
            points.resize(extents[ret.size()][2]);
            for(size_t i=0;i<ret.size();++i) {
                points[i][0] = ret[i].x;
                points[i][1] = ret[i].y;
            }
            kd_tree.reset(new kdtree2(points));
        }
    }
    return ret;
}

static void
Init(void)
{
    // No speed limit
    std::vector<std::tuple<float,float, bool> > default_coeffs = {
        // Like particles attract
        {
            // Coefficient
            0.00002f,
            // Power
            2,
            // Charge product?
            true,
        },
        {
            // Coefficient
            -0.0000001f,
            // Power
            3,
            // Charge product?
            false,
        },
    };
#if 0
    // No speed limit
    std::vector<std::tuple<float,float, bool> > default_coeffs = {
        // Like particles attract
        {
            // Coefficient
            0.000004f,
            // Power
            2,
            // Charge product?
            true,
        },
        /*
         // But also repel a bit closer
         {
         // Coefficient
         -0.00001,
         // Power
         2,
         // Charge product?
         false,
         },
         */
        {
            // Coefficient
            -0.00000002f,
            // Power
            3,
            // Charge product?
            false,
        },
    };
#endif
#if 0
    // Speed limit = 0.1
    std::vector<std::tuple<float,float, bool> > default_coeffs = {
        // Like particles attract
        {
            // Coefficient
            0.0001,
            // Power
            2,
            // Charge product?
            true,
        },
        /*
        // But also repel a bit closer
        {
            // Coefficient
            -0.00001,
            // Power
            2,
            // Charge product?
            false,
        },
         */
        {
            // Coefficient
            -0.0000001,
            // Power
            3,
            // Charge product?
            false,
        },
    };
#endif
    vector<ParticleType> particle_types;
    const ParticleType type_1(1,
                               default_coeffs,
                               1,  // charge
                               1   // mass
                               );
    
    const ParticleType type_2(2,
                              default_coeffs,
                              -1,  // charge
                              1   // mass
                              );
    particle_types.push_back(type_1);
    particle_types.push_back(type_2);
    
    cl_device_id device_id = 0;
    fprintf(stderr, "TODO: Try CL_DEVICE_TYPE_GPU, auto work group size");
    int err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to create a device group!\n");
        exit(1);
    }

//    std::unique_ptr<PointSearcher> point_searcher(new KDTreePointSearcher);
    std::unique_ptr<PointSearcher> point_searcher(new ConstantDensityPointSearcher(0.15f / 25.0f));
    local_system.reset(new LocalSystem(particle_types, point_searcher, n_interactions, min_d, time_epsilon, true, device_id));
    
#if 1
//    vector<Vec2f> pos = RandomParticles(2000, min_d);
    vector<Vec2f> pos = RandomParticles(1024, min_d);
    for(Vec2f const&p : pos) {
        // Need to nudge the particles a bit so they don't "fall into each other" as much
        local_system->AddParticle(ParticleState((runi() < 0.5f) ? 1 : 2, p * 0.15f, rvel() * 0.0f));
    }
#endif
    
#if 0
    const float blob_r = 0.08f;
    // Oil
    for(int i=0;i<500;++i) {
        local_system->AddParticle(ParticleState(1, Vec2f(-2 * blob_r, -blob_r) + Vec2f(blob_r * runi(), blob_r * runi()),
                                                    Vec2f(0.1f,0)));
    }
    // Water
    for(int i=0;i<500;++i) {
        local_system->AddParticle(ParticleState(2, Vec2f(2 * blob_r, -blob_r) + Vec2f(blob_r * runi(), blob_r * runi()),
                                                Vec2f(-0.1f,0)));
    }
#endif
}

/* ARGSUSED1 */
static void
Key(unsigned char key, int x, int y)
{
    switch (key) {
        case 27:
            exit(0);
        case ' ':
            Init();
            glutPostRedisplay();
            break;
        case 'p':
            
            break;
    }
}

static void
Draw(void)
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1, 1, 1, -1);
    glMatrixMode(GL_MODELVIEW);
    
    glBlendFunc(GL_ONE, GL_ONE);
    glEnable(GL_BLEND);
    
    vector<ParticleState> particles;
    local_system->GetParticles(particles);
    glPointSize(2);
    glBegin(GL_POINTS);
    for(ParticleState const&particle : particles) {
        if(particle.typeId == 1)
            glColor4f(0,1,0,0.5);
        else
            glColor4f(1,0,0,0.5);
        glVertex2f(particle.pos.x, particle.pos.y);
    }
    glEnd();
    glLineWidth(1);
    glBegin(GL_LINES);
    for(ParticleState const&particle : particles) {
        if(particle.typeId == 1)
            glColor4f(0,1,0,0.5);
        else
            glColor4f(1,0,0,0.5);
        glVertex2f(particle.pos.x, particle.pos.y);
        const Vec2f p_v = particle.pos + particle.vel * 0.03f;
        glVertex2f(p_v.x, p_v.y);
    }
    glEnd();
    
    glutSwapBuffers();
}

static void Idle(void)
{
    static int last_millis = -1;
    int this_millis = glutGet(GLUT_ELAPSED_TIME);
    if(last_millis < 0)
        last_millis = this_millis;
    else {
        int elapsed = this_millis - last_millis;
        last_millis = this_millis;
        local_system->Iterate(float(elapsed) / 1000.0f, affecters);
        glutPostRedisplay();
    }
}

int
main(int argc, char **argv)
{
    GLenum type;
    
    
    glutInit(&argc, argv);
    
    type = GLUT_RGB;
    type |= GLUT_DOUBLE;
    glutInitDisplayMode(type);
    glutInitWindowSize(1024, 1024);
    glutCreateWindow("ABGR extension");
    if (!glutExtensionSupported("GL_EXT_abgr")) {
        printf("Couldn't find abgr extension.\n");
        exit(0);
    }
#if !GL_EXT_abgr
    printf("WARNING: client-side OpenGL has no ABGR extension support!\n");
    printf("         Drawing only RGBA (and not ABGR) images and textures.\n");
#endif
    Init();
    glutKeyboardFunc(Key);
    glutDisplayFunc(Draw);
    glutIdleFunc(Idle);
    glutMainLoop();
    return 0;             /* ANSI C requires main to return int. */
}