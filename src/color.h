#ifndef COLOR_H
#define COLOR_H

#include "rtweekend.h"

#include <iostream>

void write_color(uint8_t *data, const int x, const int y, const int stride, const int samples, color pixel_color)
{

    double r = pixel_color.r();
    double g = pixel_color.g();
    double b = pixel_color.b();

    auto scale = 1.0 / samples;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Fill image data array with new pixel
    data[3 * (y * stride + x)    ] = static_cast<uint8_t>(256 * clamp(r, 0.0, 0.999));
    data[3 * (y * stride + x) + 1] = static_cast<uint8_t>(256 * clamp(g, 0.0, 0.999));
    data[3 * (y * stride + x) + 2] = static_cast<uint8_t>(256 * clamp(b, 0.0, 0.999));
}

#endif