#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdlib.h>
#include <iostream>
#include <string>
#include <cmath>

struct options_rec {
    bool render;
    int scene;
    int cores;
    float aspect_ratio;
    int image_width;
    int image_height;
    int samples_per_pixel;
    std::string filename;
};

static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << " <option(s)> [SCENE]"
              << "Options:\n"
              << "\t-h,--help\t\t\tShow this help message\n"
              << "\t-c,--core NUMBER\t\tSpecify the number of cores for multithreading\n"
              << "\t-n,--samples NUMBER\t\tSpecify the number of rays to cast per pixel\n"
              << "\t-a,--aspect,--ratio NUMBER\tSpecify the aspect ratio for the image\n"
              << "\t-x,--width WIDTH\t\tSpecify the width of the image\n"
              << "\t-y,--height HEIGHT\t\tSpecify the height of the image. If this is set, aspect ratio is silently ignored!\n"
              << "\t-f,--file STRING\t\tSpecify the filename to store the image into\n"
              << "\t-s,--scene SCENE\t\tSpecify the scene by integer ID"
              << std::endl;
}   

bool parseArguments(int argc, char *argv[], options_rec &options) {
    options.render = false;
    options.aspect_ratio = 16.0/9.0;
    options.image_width = 480;
    options.cores = 1;
    options.scene = 0;
    options.samples_per_pixel = 100;
    options.filename = "pic.jpg";
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return true;
        } else if((arg == "-c") || (arg == "--core") || (arg == "--cores")) {
            if(i + 1 < argc) {
                int cores = std::atoi(argv[++i]);
                options.cores = cores < 1 ? 1 : cores;
            } else {
                show_usage(argv[0]);
                return false;
            }
        } else if((arg == "-s") || (arg == "--scene")) {
            if(i + 1 < argc) {
                int scene = atoi(argv[++i]);
                options.scene = scene < 0 ? 0 : scene;
            } else {
                show_usage(argv[0]);
                return false;
            }
        } else if((arg == "-a") || (arg == "--aspect") || (arg == "--ratio")) {
            if(i + 1 < argc) {
                float ar = std::atof(argv[++i]);
                options.aspect_ratio = ar;
            } else {
                show_usage(argv[0]);
                return false;
            }
        } else if((arg == "-n") || (arg == "--samples")) {
            if(i + 1 < argc) {
                int samples = std::atoi(argv[++i]);
                options.samples_per_pixel = samples < 1 ? 1 : samples;
            } else {
                show_usage(argv[0]);
                return false;
            }
        } else if((arg == "-f") || (arg == "--file")) {
            if(i + 1 < argc) {
                options.filename = argv[++i];
            } else {
                show_usage(argv[0]);
                return false;
            }
        } else if((arg == "-x") || (arg == "--width")) {
            if(i + 1 < argc) {
                int width = std::atoi(argv[++i]);
                options.image_width = width < 1 ? 1 : width;
            } else {
                show_usage(argv[0]);
                return false;
            }
        } else if((arg == "-y") || (arg == "--height")) {
            if(i + 1 < argc) {
                int height = std::atoi(argv[++i]);
                options.image_height = height < 1 ? 1 : height;
            } else {
                show_usage(argv[0]);
                return false;
            }
        } else {
            options.scene = fmax(std::stoi(arg), 0);
        }
    }
    
    if(options.image_height) {
        options.aspect_ratio = ((float)options.image_width) / options.image_height;
    } else {
        options.image_height = static_cast<int>(options.image_width / options.aspect_ratio);
    }
    
    std::cout << "Rendering a " << options.image_width << "x" << options.image_height << " image of scene " << options.scene << " on " << options.cores << " core(s) in file " << options.filename << std::endl;
    options.render = true;
    return true;
}



#endif