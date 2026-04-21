#define _CRT_SECURE_NO_WARNINGS 1
#include <cmath>
#include <random>
#include <vector>
// #include <omp.h>
// Source - https://stackoverflow.com/a/60831839
// Posted by Dani Frățilă
// Retrieved 2026-04-21, License - CC BY-SA 4.0

// #include "/usr/local/opt/libomp/include/omp.h"


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#ifndef M_PI
#define M_PI 3.14159265358979323856
#endif

static std::default_random_engine engine[32];
static std::uniform_real_distribution<double> uniform(0, 1);

double sqr(double x) { return x * x; };

class Vector {
   public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const { return sqrt(norm2()); }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                  a[0] * b[1] - a[1] * b[0]);
}

class Ray {
   public:
    Ray(const Vector& origin, const Vector& unit_direction)
        : O(origin), u(unit_direction) {};
    Vector O, u;
};

class Object {
   public:
    Object(const Vector& albedo, bool mirror = false, bool transparent = false)
        : albedo(albedo), mirror(mirror), transparent(transparent) {};

    virtual bool intersect(const Ray& ray, Vector& P, double& t,
                           Vector& N) const = 0;

    Vector albedo;
    bool mirror, transparent;
};

class Sphere : public Object {
   public:
    Sphere(const Vector& center, double radius, const Vector& albedo,
           bool mirror = false, bool transparent = false)
        : ::Object(albedo, mirror, transparent), C(center), R(radius) {};

    // returns true iif there is an intersection between the ray and the sphere
    // if there is an intersection, also computes the point of intersection P,
    // t>=0 the distance between the ray origin and P (i.e., the parameter along
    // the ray) and the unit normal N
    bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {
        // TODO (lab 1) : compute the intersection (just true/false at the
        // begining of lab 1, then P, t and N as well) perf first, calculate
        // product
        auto dot_prod = dot(ray.u, ray.O - C);
        auto delta = sqr(dot_prod) - ((ray.O - C).norm2() - sqr(R));

        if (delta < 0) return false;

        auto disc = sqrt(delta);

        auto t1 = -dot_prod - disc;
        auto t2 = -dot_prod + disc;

        if (t2 < 0) return false;  // ray looking in other direction

        t = t1 >= 0 ? t1 : t2;  // keep minimal positive t

        P = ray.O + t * ray.u;
        N = P - C;
        N.normalize();

        return true;
    }

    double R;
    Vector C;
};

// I will provide you with an obj mesh loader (labs 3 and 4)
class TriangleMesh : public Object {
   public:
    TriangleMesh(const Vector& albedo, bool mirror = false,
                 bool transparent = false)
        : ::Object(albedo, mirror, transparent) {};

    bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {
        // TODO (labs 3 and 4)
        return false;
    }
};

class Scene {
   public:
    Scene() {};
    void addObject(const Object* obj) { objects.push_back(obj); }

    // returns true iif there is an intersection between the ray and any object
    // in the scene
    // if there is an intersection, also computes the point of the *nearest*
    // intersection P, t>=0 the distance between the ray origin and P (i.e., the
    // parameter along the ray) and the unit normal N. Also returns the index of
    // the object within the std::vector objects in object_id
    bool intersect(const Ray& ray, Vector& P, double& t, Vector& N,
                   int& object_id) const {
        // TODO (lab 1): iterate through the objects and check the intersections
        // with all of them, and keep the closest intersection, i.e., the one if
        // smallest positive value of t

        double min_t = std::numeric_limits<double>::max();
        Vector min_P = P;
        Vector min_N = N;

        double tmp_t = 0;
		bool intersected = false;
        for (size_t i = 0; i < objects.size(); ++i) {
            if (objects[i]->intersect(ray, min_P, tmp_t, min_N)) {
				if (tmp_t < min_t){
					min_t = tmp_t;
					object_id = i;
					N = min_N;
					P = min_P;
					t = min_t;
				}
				intersected = true;
            }
        }
        return intersected;
    }

    // return the radiance (color) along ray
    Vector getColor(const Ray& ray, int recursion_depth) {
        if (recursion_depth >= max_light_bounce) return Vector(0, 0, 0);

        // TODO (lab 1) : if intersect with ray, use the returned information to
        // compute the color ; otherwise black in lab 1, the color only includes
        // direct lighting with shadows

		// return Vector(255, 255, 255);
        Vector P, N;
        double t;
        int object_id;
        if (intersect(ray, P, t, N, object_id)) {
			// return Vector(255, 255, 255);
            if (objects[object_id]->mirror) {
                // return getColor in the reflected direction, with
                // recursion_depth+1 (recursively)
                auto ref_direction = ray.u - 2 * dot(ray.u, N) * N;
                auto new_ray = Ray(P + 1e-10*N, ref_direction);

                return getColor(new_ray, ++recursion_depth);
            }  // else

            if (objects[object_id]->transparent) {  // optional

                // auto ref_direction = ray.u - 2 * dot(ray.u, N) * N;
                // auto new_ray = Ray(P, ref_direction);

                // return getColor(new_ray, ++recursion_depth);

                // return getColor in the refraction direction, with
                // recursion_depth+1 (recursively)
            }  // else

            if (objects[object_id]->transparent) {
                double n1, n2;
                Vector normal = N;

                // Determine if ray is entering or exiting the medium
                if (dot(ray.u, N) > 0) {
                    // Ray is inside the object, exiting — flip normal and IORs
                    n1 = 1.5;  // glass
                    n2 = 1.0;  // air
                    normal = -1.0 * N;
                } else {
                    // Ray is outside, entering
                    n1 = 1.0;  // air
                    n2 = 1.5;  // glass
                }

                double n_ratio = n1 / n2;
                double cos_i = -dot(ray.u, normal);  // cos of incident angle (>0)
                double sin2_t = sqr(n_ratio) * (1.0 - sqr(cos_i));  // sin²(θ_t) via Snell

                if (sin2_t > 1.0) {
                    // Total internal reflection — treat like a mirror
                    Vector ref_direction = ray.u - 2 * dot(ray.u, N) * N;
                    Ray reflected_ray(P + 1e-10 * N, ref_direction);
                    return getColor(reflected_ray, recursion_depth + 1);
                }

                // Refracted ray direction (Snell's law in vector form)
                double cos_t = sqrt(1.0 - sin2_t);
                Vector refracted_dir = n_ratio * ray.u + (n_ratio * cos_i - cos_t) * normal;
                refracted_dir.normalize();

                // Offset origin inward (against the normal) to avoid self-intersection
                Ray refracted_ray(P - 1e-10 * normal, refracted_dir);
                return getColor(refracted_ray, recursion_depth + 1);
            }

        	// return Vector(255, 255, 255);
			

            // test if there is a shadow by sending a new ray
            Vector direction_to_light = light_position - P;
            auto distance_to_light = direction_to_light.norm();
            direction_to_light.normalize();

            auto shadow_ray = Ray(P + 1e-10*N, direction_to_light);

            Vector shadow_P, shadow_N;
            int shadow_id;
            double shadow_t;
            bool in_shadow = false;
            if (intersect(shadow_ray, shadow_P, shadow_t, shadow_N,
                          shadow_id)) {
                if (shadow_t < distance_to_light) in_shadow = true;
            }
            auto direct_light = Vector(0, 0, 0);
            auto indirect_light = Vector(0, 0, 0);
            if (!in_shadow)
                direct_light = (light_intensity /
                    (4 * M_PI * distance_to_light * distance_to_light)) *
                   (objects[object_id]->albedo / M_PI) *
                   std::max(0.0, dot(N, direction_to_light));

            // TODO (lab 2) : add indirect lighting component with a recursive
            // call
            double r1 = uniform(engine[0]);
            double r2 = uniform(engine[0]);

            // double r1 = uniform(engine[omp_get_thread_num()]);
            // double r2 = uniform(engine[omp_get_thread_num()]);

            // cosine-weighted hemisphere sampling
            double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
            double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
            double z = sqrt(r2);

            // now create tangents
            Vector tangent1, tangent2;
            if ((fabs(N[2]) < fabs(N[1])) && (fabs(N[2]) < fabs(N[0]))) {
                tangent1 = Vector(-N[1], N[0], 0);
            } else if ((fabs(N[1]) < fabs(N[0])) && (fabs(N[1]) < fabs(N[2]))) {
                tangent1 = Vector(0, -N[2], N[1]);
            } else {
                tangent1 = Vector(0, N[2], -N[1]);
            }
            tangent1.normalize();
            tangent2 = cross(N, tangent1);

            Vector new_direction = x * tangent1 + y * tangent2 + z * N;
            new_direction.normalize();
            auto new_ray = Ray(P + 1e-10*N, new_direction);
            indirect_light = objects[object_id]->albedo * getColor(new_ray, recursion_depth + 1);
        
            return direct_light + indirect_light;

        }

        return Vector(0, 0, 0);
    }

    std::vector<const Object*> objects;

    Vector camera_center, light_position;
    double fov, gamma, light_intensity;
    int max_light_bounce;
};

void boxMuller ( double stdev , double &x , double &y ) {
    // since MacOS issues with fopenmp, I will just use the first engine for all threads
    // double r1 = uniform ( engine[omp_get_thread_num()] ) ;
    // double r2 = uniform ( engine[omp_get_thread_num()] ) ;
    double r1 = uniform ( engine[0] ) ;
    double r2 = uniform ( engine[0] ) ;
    x = sqrt(-2 * log ( r1 +1e-15) ) *cos ( 2 * M_PI*r2 ) *stdev ;
    y = sqrt(-2 * log ( r1 + 1e-15) ) *sin ( 2 * M_PI*r2 ) *stdev ;
}

int main() {
    int W = 512;
    int H = 512;

    for (int i = 0; i < 32; i++) {
        engine[i].seed(i);
    }

    Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8), false, false);
    Sphere center_sphere_2(Vector(-10, 0, 10), 10., Vector(0.8, 0.8, 0.8), false, true);
    Sphere center_sphere_3(Vector(10, 0, -10), 10., Vector(0.8, 0.8, 0.8), true, false);
    Sphere wall_left(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
    Sphere wall_right(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
    Sphere wall_front(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
    Sphere wall_behind(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
    Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3), false);
    Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7), false);

    Scene scene;
    scene.camera_center = Vector(0, 0, 40);
    scene.light_position = Vector(-10, 20, 40);
    scene.light_intensity = 15E6;
    scene.fov = 90 * M_PI / 180.;
    scene.gamma =
        2.2;  // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
    scene.max_light_bounce = 5;

    scene.addObject(&center_sphere);
    scene.addObject(&center_sphere_2);
    scene.addObject(&center_sphere_3);


    scene.addObject(&wall_left);
    scene.addObject(&wall_right);
    scene.addObject(&wall_front);
    scene.addObject(&wall_behind);
    scene.addObject(&ceiling);
    scene.addObject(&floor);

    int SAMPLES = 128;
    double FOCAL_LENGTH = 50.;
    double APERTURE_SIZE = 1.;
    double ANTIALIASING = 0.5;

    std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color;

            // TODO (lab 1) : correct ray_direction so that it goes through each
            // pixel (j, i)
            // Vector ray_direction(j - W / 2 + 0.5, H / 2 - i - 0.5,
            //                      -W / (2 * tan(scene.fov / 2)));

			// ray_direction.normalize();
            // Ray ray(scene.camera_center, ray_direction);

            // TODO (lab 2) : add Monte Carlo / averaging of random ray
            // contributions here

            for (int s = 0; s < SAMPLES; s++) {
                // TODO (lab 2) : add antialiasing by altering the ray_direction
                // here
                double ax, ay;
                boxMuller(ANTIALIASING, ax, ay);
                Vector ray_direction(j - W / 2 + 0.5 + ax, H / 2 - i - 0.5 + ay,
                                 -W / (2 * tan(scene.fov / 2)));
                ray_direction.normalize();

                // TODO (lab 2) : add depth of field effect by altering the ray
                // origin (and direction) here

                Vector focal_point = scene.camera_center + FOCAL_LENGTH * ray_direction;
                boxMuller(APERTURE_SIZE, ax, ay);
                Vector new_origin = scene.camera_center + APERTURE_SIZE * Vector(ax, ay, 0);
                Vector new_direction = focal_point - new_origin;
                new_direction.normalize();

                Ray ray(new_origin, new_direction);
                color = color + scene.getColor(ray, 0);
            }
            color = color / SAMPLES;



            image[(i * W + j) * 3 + 0] =
                std::min(255., std::max(0., 255. * std::pow(color[0] / 255.,
                                                            1. / scene.gamma)));
            image[(i * W + j) * 3 + 1] =
                std::min(255., std::max(0., 255. * std::pow(color[1] / 255.,
                                                            1. / scene.gamma)));
            image[(i * W + j) * 3 + 2] =
                std::min(255., std::max(0., 255. * std::pow(color[2] / 255.,
                                                            1. / scene.gamma)));
        }
    }
    stbi_write_png("image_lab2.png", W, H, 3, &image[0], 0);

    return 0;
}