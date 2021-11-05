#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "rtweekend.h"
#include "sphere.h"
#include "hittable_list.h"
#include "camera.h"
#include "color.h"
#include "material.h"

color ray_color(const ray &r, const hittable &world, int depth)
{
    hit_record rec;

    if (depth <= 0)
        return color(0.0, 0.0, 0.0);
    if (world.hit(r, 0.001, infinity, rec))
    {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth - 1);
        return color(0.0, 0.0, 0.0);
    }
    vec3 unit_direction = unit_vector(r.direction());
    double t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}

int main()
{
    // Define Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 960;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int total = image_width * image_height;
    const int samples = 100;
    const int max_depth = 50;
    uint8_t image_data[total * 3];

    // Define scene

    hittable_list world;

    auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    auto material_matte  = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    //auto material_left   = make_shared<metal>(color(0.8, 0.8, 0.8), 0.0);
    auto material_glass  = make_shared<dielectric>(1.5);
    auto material_metal  = make_shared<metal>(color(0.8, 0.6, 0.2), 0.1);

    // Middle matte sphere - blue
    world.add(make_shared<sphere>(point3( 0.0, -100.5, -1.0),   100, material_ground));
    world.add(make_shared<sphere>(point3( 0.0,    0.0, -1.0),   0.5, material_matte));
    world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0),   0.5, material_glass));
    world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0), -0.45, material_glass));
    world.add(make_shared<sphere>(point3( 1.0,    0.0, -1.0),   0.5, material_metal));

    // Camera

    point3 lookfrom(0, 0, 0);
    point3 lookat(0, 0, -1);
    vec3 vup(0, 1, 0);
    camera cam(lookfrom, lookat, vup, 90, aspect_ratio);

    // Render
    
    for (int j = 0; j < image_height; j++)
    {
        for (int i = 0; i < image_width; i++)
        {
            // Pass rays through the scene
            color col(0.0, 0.0, 0.0);
            for (int s = 0; s < samples; s++)
            {
                double u = double(i + random_double()) / (image_width - 1);
                double v = double(j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);

                col += ray_color(r, world, max_depth);
            }

            // Store color data
            write_color(image_data, i, j, image_width, samples, col);

            // Print progress message
            std::cout << "Progress: " << j * image_width + i << "/" << total << " - " << 100 * double(j * image_width + i) / total << "%\r";
        }
    }
    stbi_write_jpg("pic.jpg", image_width, image_height, 3, image_data, 100);
    std::cout << "Done.                            \n";
    return 0;
}