#include <iostream>
#include <thread>
#include <future>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "rtweekend.h"
#include "sphere.h"
#include "moving_sphere.h"
#include "hittable_list.h"
#include "bvh.h"
#include "camera.h"
#include "color.h"
#include "material.h"
#include "options.h"

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

hittable_list random_scene() {
    hittable_list world;

    auto ground_material = make_shared<lambertian>(color(0.3, 0.5, 0.6));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));
    
    float matte_prob = random_double(0.0, 0.8);
    float metal_prob = random_double(0.0, 1.0);
    float glass_prob = random_double(0.0, 0.3);
    float total_prob = matte_prob + metal_prob + glass_prob;
    matte_prob /= total_prob;
    metal_prob /= total_prob;
    metal_prob += matte_prob;

    for(int a = -11; a < 11; a++) {
        for(int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            float radius = random_double(0.1, 0.4);
            point3 center(a + 0.9*random_double(), radius, b + 0.9*random_double());

            if((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;
            
                if(choose_mat < matte_prob) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, radius, sphere_material));
                } else if(choose_mat < metal_prob) {
                    // metal
                    auto albedo = color::random() * color::random();
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, radius, sphere_material));
                } else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    //Make half of the spheres hollow
                    if(random_double(0.0, 1.0) < 0.5)
                        world.add(make_shared<sphere>(center, -0.95*radius, sphere_material));
                    world.add(make_shared<sphere>(center, radius, sphere_material));
                }
            }
        }
    }
    auto material = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));
    
    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return hittable_list(make_shared<bvh_node>(world, 0.0, 1.0));
}

hittable_list threeSphereScene() {
    hittable_list world;
    
    auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    auto material_matte  = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    auto material_glass  = make_shared<dielectric>(1.5);
    auto material_metal  = make_shared<metal>(color(0.8, 0.6, 0.2), 0.0);

    world.add(make_shared<sphere>(point3( 0.0, -100.5, -1.0),   100, material_ground));
    world.add(make_shared<sphere>(point3( 0.0,    0.0, -1.0),   0.5, material_matte));
    world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0),   0.5, material_glass));
    world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0), -0.45, material_glass));
    world.add(make_shared<sphere>(point3( 1.0,    0.0, -1.0),   0.5, material_metal));
    
    return hittable_list(make_shared<bvh_node>(world, 0.0, 1.0));
}

void render(int image_height, int image_width, int samples, camera *cam, hittable_list *world, int max_depth, int total, std::vector<color>* image_data) {
    
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
                ray r = cam->get_ray(u, v);

                col += ray_color(r, *world, max_depth);
            }
            int index = j * image_width + i;

            // Store color data
            image_data->at(index) = col;
            
            // Print progress message
            //std::cout << "Progress: " << j * image_width + i << "/" << total << " - " << 100 * double(j * image_width + i) / total << "%\r";
        }
    }

}

int main(int argc, char *argv[])
{
    // Define Image from parsed arguments
    
    options_rec opts = {};

    if (!parseArguments(argc, argv, opts)) {
        std::cout << "Error while parsing arguments!" << std::endl;
        return 1;
    }
    
    if (!opts.render) {
        std::cout << "Not rendering!" << std::endl;
        return 0;
    }
    
    srand(time(0));

    const auto aspect_ratio = opts.aspect_ratio;
    const int image_width = opts.image_width;
    const int image_height = opts.image_height;
    const int total = image_width * image_height;
    const int samples = opts.samples_per_pixel;
    const int samples_per_thread = static_cast<int>(samples/opts.cores);
    const int max_depth = 50;
    uint8_t image_data[total * 3];

    // Define scene
    
    hittable_list world;
    switch (opts.scene)
    {
    case 0:
        world = random_scene();
        break;
    
    case 1:
        world = threeSphereScene();
        break;
    
    default:
        std::cout << "Specified non-existent scene. Using scene 0." << std::endl;
        world = random_scene();
        break;
    }

    // Camera

    point3 lookfrom(8, 5, 5);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    double dist_to_focus = 10.;
    double aperture = 0.05;
    camera cam(lookfrom, lookat, vup, 30, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);
    
    // Create threads

    std::vector<std::thread> threads;
    std::vector<std::vector<color>> data_arrays = std::vector<std::vector<color>>(opts.cores, std::vector<color>(total));
    for (int i = 0; i < opts.cores; i++) {
        threads.push_back(std::thread(render, image_height, image_width, samples_per_thread, &cam, &world, max_depth, total, &(data_arrays[i])));
    }
    
    std::vector<color> pixel_colors(total);
    for (int i = 0; i < opts.cores; i++) {
        threads[i].join();
        for (int y = 0; y < image_height; y++) {
            for (int x = 0; x < image_width; x++) {
                int index = y * image_width + x;
                pixel_colors[index] += data_arrays[i][index];
            }
        }
    }
    
    write_color(image_data, image_width, image_height, samples, pixel_colors);

    char c[opts.filename.size() + 1];
    stbi_write_jpg(opts.filename.c_str(), image_width, image_height, 3, image_data, 100);
    std::cout << "Done.                                   \n";
    return 0;
}