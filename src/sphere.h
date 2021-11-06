#ifndef SPHEREH
#define SPHEREH

#include "rtweekend.h"

#include "hittable.h"

class sphere : public hittable
{
public:
    sphere() {}
    sphere(point3 cen, double r, shared_ptr<material> mat) : center(cen), radius(r), mat_ptr(mat){};
    virtual bool hit(const ray &r, double tmin, double tmax, hit_record &rec) const override;
    point3 center;
    double radius;
    shared_ptr<material> mat_ptr;
};

bool sphere::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    vec3 oc = r.origin() - center;
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
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    return true;
}

#endif