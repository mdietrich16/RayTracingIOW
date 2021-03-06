#ifndef SCENES_H
#define SCENES_H

#include "rtweekend.h"
#include "camera.h"
#include "hittable_list.h"
#include "sphere.h"
#include "moving_sphere.h"
#include "box.h"
#include "bvh.h"
#include "constant_medium.h"

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

hittable_list simple_light() {
    hittable_list objects;

    auto pertext = make_shared<noise_texture>(4);
    objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    auto difflight = make_shared<diffuse_light>(color(4, 4, 4));
    objects.add(make_shared<xy_rect>(3, 4, 1, 3, -2, difflight));

    return objects;
}

hittable_list cornell_box() {
    hittable_list objects;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(15, 15, 15));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));
    objects.add(box1);

    shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(130, 0, 65));
    objects.add(box2);

    return objects;
}

hittable_list volume_cornell_box() {
    hittable_list objects;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(7, 7, 7));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<xz_rect>(113, 443, 127, 432, 554, light));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));

    shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(130, 0, 65));

    objects.add(make_shared<constant_medium>(box1, 0.01, color(0, 0, 0)));
    objects.add(make_shared<constant_medium>(box2, 0.01, color(1, 1, 1)));

    return objects;
}

hittable_list final_scene() {
    hittable_list boxes1;
    auto ground = make_shared<lambertian>(color(0.48, 0.83, 0.53));

    const int boxes_per_side = 20;
    for (int i = 0; i < boxes_per_side; i++) {
        for (int j = 0; j < boxes_per_side; j++) {
            auto w = 100.0;
            auto x0 = -1000.0 + i*w;
            auto z0 = -1000.0 + j*w;
            auto y0 = 0.0;
            auto x1 = x0 + w;
            auto z1 = z0 + w;
            auto y1 = random_double(1, 101);
            
            boxes1.add(make_shared<box>(point3(x0, y0, z0), point3(x1, y1, z1), ground));
        }
    }
    
    hittable_list objects;

    objects.add(make_shared<bvh_node>(boxes1, 0, 1));

    auto light = make_shared<diffuse_light>(color(7, 7, 7));
    objects.add(make_shared<xz_rect>(123, 423, 147, 412, 554, light));

    auto center1 = point3(400, 400, 200);
    auto center2 = center1 + vec3(30, 0, 0);
    auto moving_sphere_material = make_shared<lambertian>(color(0.7, 0.3, 0.1));
    objects.add(make_shared<moving_sphere>(center1, center2, 0, 1, 50, moving_sphere_material));

    objects.add(make_shared<sphere>(point3(260, 150, 45), 50, make_shared<dielectric>(1.5)));
    objects.add(make_shared<sphere>(point3(0, 150, 145), 50, make_shared<metal>(color(0.8, 0.8, 0.9), 1.0)));

    auto boundary = make_shared<sphere>(point3(360, 150, 145), 70, make_shared<dielectric>(1.5));
    objects.add(boundary);
    objects.add(make_shared<constant_medium>(boundary, 0.2, color(0.2, 0.4, 0.9)));

    boundary = make_shared<sphere>(point3(0, 0, 0), 5000, make_shared<dielectric>(1.5));
    objects.add(make_shared<constant_medium>(boundary, 0.0001, color(1, 1, 1)));
    
    auto emat = make_shared<lambertian>(make_shared<image_texture>("assets/earthmap.jpg"));
    objects.add(make_shared<sphere>(point3(400, 200 ,400), 100, emat));

    auto pertext = make_shared<noise_texture>(0.1);
    objects.add(make_shared<sphere>(point3(220, 280, 300), 80, make_shared<lambertian>(pertext)));

    hittable_list boxes2;
    auto white = make_shared<lambertian>(color(0.73, 0.73, 0.73));
    int ns = 1000;
    for (int j = 0; j < ns; j++) {
        boxes2.add(make_shared<sphere>(point3::random(0, 165), 10, white));
    }
    
    objects.add(make_shared<translate>(make_shared<rotate_y>(make_shared<bvh_node>(boxes2, 0.0, 1.0), 15), vec3(-100, 270, 395)));

    return objects;
}

bool choose_scene(int scene, double aspect_ratio, camera& cam, hittable_list& world, color& background) {
    point3 lookfrom;
    point3 lookat = point3(0, 0, 0);
    double aperture = 0.0;
    double vfov = 20.0;
    
    switch (scene)
    {
    case 0:
        world = random_scene();
        background = color(0.7, 0.8, 1.0);
        lookfrom = point3(13, 2, 3);
        break;
    
    case 1:
        world = three_spheres();
        background = color(0.7, 0.8, 1.0);
        lookfrom = point3(0, 1, 1);
        break;
    
    default:
        std::cout << "Specified non-existent scene. Using scene 2." << std::endl;
    case 2:
        world = two_spheres();
        background = color(0.7, 0.8, 1.0);
        lookfrom = point3(13, 2, 3);
        break;

    case 3:
        world = two_perlin_spheres();
        background = color(0.7, 0.8, 1.0);
        lookfrom = point3(13, 2, 3);
        break;

    case 4:
        world = earth();
        background = color(0.7, 0.8, 1.0);
        lookfrom = point3(13, 2, 3);
        break;

    case 5:
        world = simple_light();
        lookfrom = point3(26, 3, 6);
        lookat = point3(0, 2, 0);
        break;

    case 6:
        world = cornell_box();
        lookfrom = point3(278, 278, -800);
        lookat = point3(278, 278, 0);
        vfov = 40;
        break;
 
    case 7:
        world = volume_cornell_box();
        lookfrom = point3(278, 278, -800);
        lookat = point3(278, 278, 0);
        vfov = 40;
        break;
 
    case 8:
        world = final_scene();
        lookfrom = point3(478, 278, -600);
        lookat = point3(278, 278, 0);
        vfov = 40;
        break;
    }

    // Camera

    vec3 vup(0, 1, 0);
    double dist_to_focus = (lookat - lookfrom).length();
    cam = camera(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);
 
    return true;
}

#endif