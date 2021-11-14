#include <iostream>
#include <iomanip>
#include <thread>
#include <future>
#include <atomic>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb_image_write.h"

#include "rtweekend.h"
#include "sphere.h"
#include "aarect.h"
#include "box.h"
#include "moving_sphere.h"
#include "hittable_list.h"
#include "bvh.h"
#include "camera.h"
#include "color.h"
#include "material.h"
#include "constant_medium.h"
#include "options.h"
#include "scenes.h"

color ray_color(const ray &r, const color& background, const hittable &world, int depth)
{
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0.0, 0.0, 0.0);
    
    // If the ray hits nothing, return the background color.
    if(!world.hit(r, 0.001, infinity, rec))
        return background;

    ray scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
}

void render(int seed, int image_height, int image_width, int samples, camera *cam, hittable_list *world, const color* background, int max_depth, int total, std::atomic<int>* alive_count, std::atomic<int>* finished_pixel_parts,  std::vector<color>* image_data) {
    
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
                col += ray_color(r, *background, *world, max_depth);
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
    
    camera cam = camera(point3(0,0,0), point3(0,0,-1), vec3(0, 1, 0), 40., aspect_ratio, .2, 10., 0, 1);
    hittable_list world;
    color background(0, 0, 0);
    choose_scene(opts.scene, aspect_ratio, cam, world, background);
   
    // Time-keeping and stats
    
    std::chrono::time_point start_time = std::chrono::high_resolution_clock::now();
    std::atomic<int> finished_pixel_parts = 0;
    std::atomic<int> alive_count = 0;

    // Create threads
    
    std::vector<std::thread> threads;
    std::vector<std::vector<color>> data_arrays = std::vector<std::vector<color>>(opts.cores, std::vector<color>(total));
    for (int i = 0; i < opts.cores; i++) {
        threads.push_back(std::thread(render, opts.seed + i, image_height, image_width, samples_per_thread, &cam, &world, &background, max_depth, total, &alive_count, &finished_pixel_parts, &(data_arrays[i])));
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

    stbi_write_bmp(opts.filename.c_str(), image_width, image_height, 3, image_data);

    std::chrono::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration elapsed = (end_time - start_time);
    const int elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
    std::cout << std::endl << "Done. Took " << elapsed_sec << " seconds." << std::endl;
    std::cout << "That means:" << std::endl;
    std::cout << "\t-> " << std::scientific << std::setprecision(2) << float(total*real_samples)/elapsed_sec << " rays/sec" << std::endl;
    std::cout << "\t-> " << std::scientific << std::setprecision(2) << float(elapsed_sec)/(total) << " sec/pixel" << std::endl;

    delete[] image_data;
    return 0;
}