#define _CRT_SECURE_NO_WARNINGS 1
#include <cmath>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
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

// Class only used in labs 3 and 4 
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1) {
		vtx[0] = vtxi; vtx[1] = vtxj; vtx[2] = vtxk;
		uv[0] = uvi; uv[1] = uvj; uv[2] = uvk;
		n[0] = ni; n[1] = nj; n[2] = nk;
		this->group = group;
	};
	int vtx[3]; // indices within the vertex coordinates array
	int uv[3];  // indices within the uv coordinates array
	int n[3];   // indices within the normals array
	int group;  // face group
};


// Class only used in labs 3 and 4 
class TriangleMesh : public Object {
public:
	TriangleMesh(const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent) {};

	// first scale and then translate the current object
	void scale_translate(double s, const Vector& t) {
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i] * s + t;
		}
	}

	// read an .obj file
	void readOBJ(const char* obj) {
		std::ifstream f(obj);
		if (!f) return;

		std::map<std::string, int> mtls;
		int curGroup = -1, maxGroup = -1;

		// OBJ indices are 1-based and can be negative (relative), this normalizes them
		auto resolveIdx = [](int i, int size) {
			return i < 0 ? size + i : i - 1;
		};

		auto setFaceVerts = [&](TriangleIndices& t, int i0, int i1, int i2) {
			t.vtx[0] = resolveIdx(i0, vertices.size());
			t.vtx[1] = resolveIdx(i1, vertices.size());
			t.vtx[2] = resolveIdx(i2, vertices.size());
		};
		auto setFaceUVs = [&](TriangleIndices& t, int j0, int j1, int j2) {
			t.uv[0] = resolveIdx(j0, uvs.size());
			t.uv[1] = resolveIdx(j1, uvs.size());
			t.uv[2] = resolveIdx(j2, uvs.size());
		};
		auto setFaceNormals = [&](TriangleIndices& t, int k0, int k1, int k2) {
			t.n[0] = resolveIdx(k0, normals.size());
			t.n[1] = resolveIdx(k1, normals.size());
			t.n[2] = resolveIdx(k2, normals.size());
		};

		std::string line;
		while (std::getline(f, line)) {
			// Trim trailing whitespace
			line.erase(line.find_last_not_of(" \r\t\n") + 1);
			if (line.empty()) continue;

			const char* s = line.c_str();

			if (line.rfind("usemtl ", 0) == 0) {
				std::string matname = line.substr(7);
				auto result = mtls.emplace(matname, maxGroup + 1);
				if (result.second) {
					curGroup = ++maxGroup;
				} else {
					curGroup = result.first->second;
				}
			} else if (line.rfind("vn ", 0) == 0) {
				Vector v;
				sscanf(s, "vn %lf %lf %lf", &v[0], &v[1], &v[2]);
				normals.push_back(v);
			} else if (line.rfind("vt ", 0) == 0) {
				Vector v;
				sscanf(s, "vt %lf %lf", &v[0], &v[1]);
				uvs.push_back(v);
			} else if (line.rfind("v ", 0) == 0) {
				Vector pos, col;
				if (sscanf(s, "v %lf %lf %lf %lf %lf %lf", &pos[0], &pos[1], &pos[2], &col[0], &col[1], &col[2]) == 6) {
					for (int i = 0; i < 3; i++) col[i] = std::min(1.0, std::max(0.0, col[i]));
					vertexcolors.push_back(col);
				} else {
					sscanf(s, "v %lf %lf %lf", &pos[0], &pos[1], &pos[2]);
				}
				vertices.push_back(pos);
			}
			else if (line[0] == 'f') {
				int i[4], j[4], k[4], offset, nn;
				const char* cur = s + 1;
				TriangleIndices t;
				t.group = curGroup;

				// Try each face format: v/vt/vn, v/vt, v//vn, v
				if ((nn = sscanf(cur, "%d/%d/%d %d/%d/%d %d/%d/%d%n", &i[0], &j[0], &k[0], &i[1], &j[1], &k[1], &i[2], &j[2], &k[2], &offset)) == 9) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d/%d %d/%d %d/%d%n", &i[0], &j[0], &i[1], &j[1], &i[2], &j[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]);
				} else if ((nn = sscanf(cur, "%d//%d %d//%d %d//%d%n", &i[0], &k[0], &i[1], &k[1], &i[2], &k[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d %d %d%n", &i[0], &i[1], &i[2], &offset)) == 3) {
					setFaceVerts(t, i[0], i[1], i[2]);
				}
				else continue;

				indices.push_back(t);
				cur += offset;

				// Fan triangulation for polygon faces (4+ vertices)
				while (*cur && *cur != '\n') {
					TriangleIndices t2;
					t2.group = curGroup;
					if ((nn = sscanf(cur, " %d/%d/%d%n", &i[3], &j[3], &k[3], &offset)) == 3) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d/%d%n", &i[3], &j[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]);
					} else if ((nn = sscanf(cur, " %d//%d%n", &i[3], &k[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d%n", &i[3], &offset)) == 1) {
						setFaceVerts(t2, i[0], i[2], i[3]);
					} else { 
						cur++; 
						continue; 
					}

					indices.push_back(t2);
					cur += offset;
					i[2] = i[3]; j[2] = j[3]; k[2] = k[3];
				}
			}
		}
	}
	

	// TODO ray-mesh intersection (labs 3 and 4)
	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {
		
		// lab 3 : for each triangle, compute the ray-triangle intersection with Moller-Trumbore algorithm
		// lab 3 : once done, speed it up by first checking against the mesh bounding box
		// lab 4 : recursively apply the bounding-box test from a BVH datastructure
        bool intersected = false;
        auto t_min = 1e100;

        for (auto& tri : indices) {
            auto v0 = vertices[tri.vtx[0]];
            auto v1 = vertices[tri.vtx[1]];
            auto v2 = vertices[tri.vtx[2]];
            auto edge1 = v1 - v0;
            auto edge2 = v2 - v0;

            auto beta = dot(edge2, cross(ray.u, (v0 - ray.O))) / dot(ray.u, cross(edge1, edge2));

            auto gamma = - dot(edge1, cross(ray.u, (v0 - ray.O))) / dot(ray.u, cross(edge1, edge2));

            auto alpha = 1 - beta - gamma;
            auto tri_t = dot(cross(edge1, edge2), (v0 - ray.O)) / dot(ray.u, cross(edge1, edge2));

            if (tri_t > 0 && alpha >= 0 && beta >= 0 && gamma >= 0) {
                if (tri_t < t_min) {
                    t_min = tri_t;;
                    t = tri_t;
                    intersected = true;
                    P = ray.O + t * ray.u;
                    N = cross(edge1, edge2);
                    N.normalize();
                    if (dot(ray.u, N) > 0) N = -1.0 * N;
                }
            }
        }

		return intersected;
	}


	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
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

    TriangleMesh cat(Vector(0.6, 0.5, 0.4));
    cat.readOBJ("cat.obj");
    cat.scale_translate(0.6, Vector(0, -10, 0));
    // cat.computeBoundingBox();

    

    Scene scene;
    scene.camera_center = Vector(0, 0, 40);
    scene.light_position = Vector(-10, 20, 40);
    scene.light_intensity = 15E6;
    scene.fov = 90 * M_PI / 180.;
    scene.gamma =
        2.2;  // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
    scene.max_light_bounce = 15;

    // scene.addObject(&center_sphere);
    // scene.addObject(&center_sphere_2);
    // scene.addObject(&center_sphere_3);


    scene.addObject(&wall_left);
    scene.addObject(&wall_right);
    scene.addObject(&wall_front);
    scene.addObject(&wall_behind);
    scene.addObject(&ceiling);
    scene.addObject(&floor);

    scene.addObject(&cat);

    int SAMPLES = 8;
    double FOCAL_LENGTH = 40.;
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
    stbi_write_png("image_lab4.png", W, H, 3, &image[0], 0);

    return 0;
}