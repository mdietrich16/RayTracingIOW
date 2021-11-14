#ifndef CONSTANT_MEDIUM_H
#define CONSTANT_MEDIUM_H

#include "rtweekend.h"

#include "hittable.h"
#include "material.h"
#include "texture.h"

class constant_medium : public hittable {
    public:
        constant_medium(shared_ptr<hittable> b, double d, shared_ptr<texture> a)
            : boundary(b), neg_inv_density(-1/d), phase_function(make_shared<isotropic>(a)) {}
        
        constant_medium(shared_ptr<hittable> b, double d, color c)
            : boundary(b), neg_inv_density(-1/d), phase_function(make_shared<isotropic>(c)) {}
        
        
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
            return boundary->bounding_box(time0, time1, output_box);
        }
        
        shared_ptr<hittable> boundary;
        shared_ptr<material> phase_function;
        double neg_inv_density;
};


bool constant_medium::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    const bool enableDebug = false;
    const bool debugging = enableDebug && random_double() < 0.00001;

    hit_record rec1, rec2;

    if (!boundary->hit(r, -infinity, infinity, rec1))
        return false;
    
    if (!boundary->hit(r, rec1.t+0.001, infinity, rec2))
        return false;

    if(debugging)
        std::cerr << std::endl << "t_min=" << rec1.t << ", t_max=" << rec1.t << std::endl;
    
    if (rec1.t < t_min) rec1.t = t_min;
    if (rec2.t > t_max) rec2.t = t_max;

    if (rec1.t >= rec2.t)
        return false;
    
    if (rec1.t < 0)
        rec1.t = 0;
    
    const auto ray_lengh = r.direction().length();
    const auto distance_inside_boundary = (rec2.t - rec1.t) * ray_lengh;
    const auto hit_distance = neg_inv_density * log(random_double());

    if(hit_distance > distance_inside_boundary)
        return false;
    
    rec.t = rec1.t + hit_distance / ray_lengh;
    rec.p = r.at(rec.t);

    if (debugging) {
        std::cerr << "hit_distance = " << hit_distance << std::endl << "rec.t = " << rec.t << std::endl << "rec.p = " << rec.p << std::endl;
    }
    
    rec.normal = vec3(1, 0, 0);
    rec.front_face = true;
    rec.mat_ptr = phase_function;

    return true;
}

#endif