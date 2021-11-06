#ifndef CAMERAH
#define CAMERAH

#include "rtweekend.h"

class camera
{
public:
    camera(
        point3 lookfrom,
        point3 lookat,
        vec3 vup,
        double vfov,//vfov is top to bottom in degrees
        double aspect_ratio,
        double aperture,
        double focus_dist,
        double time0,
        double time1)
    {
        double theta = deegrees_to_radians(vfov);
        double h = tan(theta / 2);
        double viewport_height = 2.0 * h;
        double viewport_width = aspect_ratio * viewport_height;
        
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(u, w);

        origin = lookfrom;
        horizontal = focus_dist * viewport_width * u;
        vertical = focus_dist * viewport_height * v;
        upper_left_corner = origin - horizontal/2 - vertical/2 - focus_dist * w;
        lens_radius = aperture/2;
        tm0 = time0;
        tm1 = time1;
    }
    ray get_ray(double s, double t) {
        vec3 rd = lens_radius * random_in_unit_disk();
        vec3 offset = u * rd.x() + v * rd.y();

        return ray(
            origin + offset,
            upper_left_corner + s*horizontal + t*vertical - origin - offset,
            random_double(tm0, tm1)
        );
    }

private:
    point3 origin;
    point3 upper_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 u, v, w;
    double lens_radius;
    double tm0, tm1; // Shutter open/close times
};

#endif