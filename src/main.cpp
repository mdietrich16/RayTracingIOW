#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "sphere.h"
#include "hitable_list.h"
#include "camera.h"
#include "material.h"

vec3 color(const ray &r, hitable *world, int depth)
{
    hit_record rec;
    if (world->hit(r, 0.001, MAXFLOAT, rec))
    {
        ray scattered;
        vec3 attenuation;
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation*color(scattered, world, depth+1);
        }
        else {
            return vec3(0.0, 0.0, 0.0);
        }
    }
    else
    {
        vec3 unit_direction = unit_vector(r.direction());
        double t = 0.5 * (unit_direction.y() + 1.0);
        return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
    }
}

int main()
{
    // Define Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 1920;
    const int image_height = (int)(image_width / aspect_ratio);
    const int total = image_width*image_height;
    const int ns = 20;
    uint8_t image_data[total*3];
    
    // Define scene
    hitable *list[5];
    double R = cos(M_PI/4);
    // Middle matte sphere - blue
    list[0] = new sphere(vec3(0.0, 0.0, -1.0), 0.5, new lambertian(vec3(0.5, 0.2, 0.5)));
    // Big sphere/ground - green
    list[1] = new sphere(vec3(0.0, -100.5, -1.0), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
    // Right metal sphere - slightly bluish
    list[2] = new sphere(vec3(1.0, 0.0, -1.0), 0.5, new metal(vec3(0.6, 0.6, 0.7), 0.0));
    // Outer shell of glass ball
    list[3] = new sphere(vec3(-1.0, 0.0, -1.0), 0.5, new dielectric(1.5));
    // Inner shell og glass ball
    list[4] = new sphere(vec3(-1.0, 0.0, -1.0), -0.45, new dielectric(1.5));
    hitable *world = new hitable_list(list, 5);
    vec3 lookfrom(-2, 2, 1);
    vec3 lookat(0, 0, -1);
    vec3 vup(0, 1, 0);
    camera cam(lookfrom, lookat, vup, 90, aspect_ratio);
    
    // Render
    for (int j = 0; j < image_height; j++)
    {
        for (int i = 0; i < image_width; i++)
        {
            vec3 col(0.0, 0.0, 0.0);
            
            for (int s = 0; s < ns; s++)
            {
                double u = double(i + drand48()) / (image_width);
                double v = double(j + drand48()) / (image_height);
                ray r = cam.get_ray(u, v);
                
                col += color(r, world, 0);
            }
            
            col /= double(ns);
            col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
            int ir = int(255.99 * col.r());
            int ig = int(255.99 * col.g());
            int ib = int(255.99 * col.b());
            image_data[3*(j*image_width + i)] = ir;
            image_data[3*(j*image_width + i) + 1] = ig;
            image_data[3*(j*image_width + i) + 2] = ib;
            std::cout << "Progress: " << j*image_width + i << "/" << total << " - " << 100*double(j*image_width + i)/total << "%\r";
        }
    }
    stbi_write_jpg("pic.jpg", image_width, image_height, 3, image_data, 100);
    return 0;
}