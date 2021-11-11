#include <iostream>
#include <iomanip>
#include <thread>
#include <future>
#include <atomic>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb_image_write.h"

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

    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));
    
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
                    //Make half of the spheres hollow -- Don't work because of BVH
                    // if(random_double(0.0, 1.0) < 0.5) {
                    //     std::cout << "Made hollow sphere" << std::endl;
                    // }
                    world.add(make_shared<sphere>(center, radius, sphere_material));
                }
            }
        }
    }
    auto material = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material));
    // Hollow spheres don't work because of BVH?
    //world.add(make_shared<sphere>(point3(0, 1, 0), -0.9, material));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));
    
    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return hittable_list(make_shared<bvh_node>(world, 0.0, 1.0));
}

hittable_list two_spheres() {
    hittable_list objects;

    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
    objects.add(make_shared<sphere>(point3(0, -10, 0), 10, make_shared<lambertian>(checker)));
    objects.add(make_shared<sphere>(point3(0,  10, 0), 10, make_shared<lambertian>(checker)));
    
    return objects;
}

hittable_list two_perlin_spheres() {
    hittable_list objects;

    auto pertext = make_shared<noise_texture>(4);
    auto turbtext = make_shared<turb_noise_texture>(4);
    auto marbletext = make_shared<marble_texture>(4);
    objects.add(make_shared<sphere>(point3(0, -1000,  0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(point3(0,     2, -2),    2, make_shared<lambertian>(turbtext)));
    objects.add(make_shared<sphere>(point3(0,     2,  2),    2, make_shared<lambertian>(marbletext)));
    
    return objects;
}

hittable_list three_spheres() {
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

hittable_list earth() {
    auto earth_texture = make_shared<image_texture>("assets/earthmap.jpg");
    auto earth_surface = make_shared<lambertian>(earth_texture);
    auto globe = make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);

    return hittable_list(globe);
}

void render(int seed, int image_height, int image_width, int samples, camera *cam, hittable_list *world, int max_depth, int total, std::atomic<int>* alive_count, std::atomic<int>* finished_pixel_parts,  std::vector<color>* image_data) {
    
    (*alive_count)++;
    srand(seed);
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
            
            (*finished_pixel_parts)++;
        }
    }
    (*alive_count)--;
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

    std::cout << "Rendering a " << opts.image_width << "x" << opts.image_height << " image of scene " << opts.scene << " in file " << opts.filename << std::endl;
    

    std::cout << "Using baseseed " << opts.seed << " for RNG." << std::endl;
    srand(opts.seed);

    const auto aspect_ratio = opts.aspect_ratio;
    const int image_width = opts.image_width;
    const int image_height = opts.image_height;
    const int total = image_width * image_height;
    const int samples = opts.samples_per_pixel;
    const float thread_count = opts.cores;
    const int samples_per_thread = static_cast<int>(ceil(float(samples)/thread_count));
    const int real_samples = samples_per_thread*opts.cores;
    const int max_depth = 50;
    uint8_t* image_data = new uint8_t[total * 3];
  
    std::cout << "Computing " << samples_per_thread << " rays on " << opts.cores << " core(s) for a total of " << real_samples << " rays per pixel or " << real_samples*total << " total rays." << std::endl;

    // Define scene and camera from arguments

    point3 lookfrom;
    point3 lookat;
    double aperture = 0.0;
    double vfov = 40.0;
    
    hittable_list world;
    switch (opts.scene)
    {
    case 0:
        world = random_scene();
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20;
        aperture = 0.1;
        break;
    
    case 1:
        world = three_spheres();
        lookfrom = point3(0, 1, 1);
        lookat = point3(0, 0, 0);
        vfov = 20;
        break;
    
    default:
        std::cout << "Specified non-existent scene. Using scene 2." << std::endl;
    case 2:
        world = two_spheres();
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20;
        break;

    case 3:
        world = two_perlin_spheres();
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20;
        break;

    case 4:
        world = earth();
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20;
        break;
    }

    // Camera

    vec3 vup(0, 1, 0);
    double dist_to_focus = 10.;
    camera cam(lookfrom, lookat, vup, 30, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);
    
    // Time-keeping and stats
    
    std::chrono::time_point start_time = std::chrono::high_resolution_clock::now();
    std::atomic<int> finished_pixel_parts = 0;
    std::atomic<int> alive_count = 0;

    // Create threads
    
    std::vector<std::thread> threads;
    std::vector<std::vector<color>> data_arrays = std::vector<std::vector<color>>(opts.cores, std::vector<color>(total));
    for (int i = 0; i < opts.cores; i++) {
        threads.push_back(std::thread(render, opts.seed + i, image_height, image_width, samples_per_thread, &cam, &world, max_depth, total, &alive_count, &finished_pixel_parts, &(data_arrays[i])));
    }
    
    float progress = 0.f;
    int str_width = static_cast<int>(ceil(log10(total)));
    while(alive_count > 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        // Print progress message
        progress = float(finished_pixel_parts)/thread_count;
        std::cout << "Progress: " << std::fixed << std::setw(str_width) << std::setprecision(0) << progress << "/" << total << " - " << std::fixed << std::setprecision(1) << std::setw(5) << 100.f * progress / total << "%\r";
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
    
    write_color(image_data, image_width, image_height, real_samples, pixel_colors);

    stbi_write_jpg(opts.filename.c_str(), image_width, image_height, 3, image_data, 100);

    std::chrono::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration elapsed = (end_time - start_time);
    const int elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
    std::cout << std::endl << "Done. Took " << elapsed_sec << " seconds." << std::endl;
    std::cout << "That means:" << std::endl;
    std::cout << "\t-" << std::scientific << std::setprecision(2) << float(total*real_samples)/elapsed_sec << " rays/sec" << std::endl;
    std::cout << "\t-" << std::scientific << std::setprecision(2) << (elapsed_sec)/(total) << " sec/pixel" << std::endl;

    delete[] image_data;
    return 0;
}