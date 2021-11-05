#ifndef CAMERAH
#define CAMERAH

#include "rtweekend.h"

class camera
{
public:
    camera() : camera(point3(0, 0, -1), point3(0, 0, 0), vec3(0, 1, 0), 90, 16.0 / 9.0) {}
    camera(point3 lookfrom, point3 lookat, vec3 vup, double vfov, double aspect)
    { //vfov is top to bottom in degrees
        point3 u, v, w;
        double theta = vfov * M_PI / 180;
        double half_height = tan(theta / 2);
        double half_width = aspect * half_height;
        origin = lookfrom;
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(u, w);
        //lower_left_corner = point3(-half_width, -half_height, -1.0);
        upper_left_corner = origin - half_width * u - half_height * v - w;
        horizontal = 2 * half_width * u;
        vertical = 2 * half_height * v;
    }
    ray get_ray(double s, double t) { return ray(origin, upper_left_corner + s * horizontal + t * vertical - origin); }

    point3 origin;
    point3 upper_left_corner;
    vec3 horizontal;
    vec3 vertical;
};

#endif