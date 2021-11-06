#ifndef MOVING_SPHERE_H
#define MOVING_SPHERE_H

#include "rtweekend.h"

#include "hittable.h"

class moving_sphere : public hittable {
    public:
        moving_sphere() {}
        moving_sphere(
            point3 center0,
            point3 center1,
            double time0,
            double time1,
            double r,
            shared_ptr<material> m
        ) : cen0(center0), cen1(center1), tm0(time0), tm1(time1), radius(r), mat_ptr(m) {};
    
        virtual bool hit(
            const ray &r,
            double t_min,
            double t_max,
            hit_record &rec) const override;
        
        point3 center(double time) const;
            
        point3 cen0, cen1;
        double tm0, tm1;
        double radius;
        shared_ptr<material> mat_ptr;
};

point3 moving_sphere::center(double time) const {
    return cen0 + ((time - tm0) / (tm1 - tm0)) * (cen1 - cen0);
}

bool moving_sphere::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    vec3 oc = r.origin() - center(r.time());
    double a = r.direction().length_squared();
    double half_b = dot(oc, r.direction());
    double c = oc.length_squared() - radius * radius;
    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0)
        return false;
    double sqrtd = sqrt(discriminant);
    double root = (-half_b - sqrtd) / a;
    if (root < t_min || root > t_max)
    {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || root > t_max)
        {
            return false;
        }
    }
    rec.mat_ptr = mat_ptr;
    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center(r.time())) / radius;
    rec.set_face_normal(r, outward_normal);
    return true;
}

#endif