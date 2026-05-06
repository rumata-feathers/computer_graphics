#define _CRT_SECURE_NO_WARNINGS 1
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
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

// Class only used in labs 3 and 4
class TriangleIndices {
   public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1,
                    int nj = -1, int nk = -1, int uvi = -1, int uvj = -1,
                    int uvk = -1, int group = -1) {
        vtx[0] = vtxi;
        vtx[1] = vtxj;
        vtx[2] = vtxk;
        uv[0] = uvi;
        uv[1] = uvj;
        uv[2] = uvk;
        n[0] = ni;
        n[1] = nj;
        n[2] = nk;
        this->group = group;
    };
    int vtx[3];  // indices within the vertex coordinates array
    int uv[3];   // indices within the uv coordinates array
    int n[3];    // indices within the normals array
    int group;   // face group
};

struct BBox {
    Vector min;
    Vector max;
    bool is_hit(const Ray& ray, double t_min, double t_max) const {
        auto t_min_vec = Vector{(min[0] - ray.O[0]) / ray.u[0],
                                (min[1] - ray.O[1]) / ray.u[1],
                                (min[2] - ray.O[2]) / ray.u[2]};
        auto t_max_vec = Vector{(max[0] - ray.O[0]) / ray.u[0],
                                (max[1] - ray.O[1]) / ray.u[1],
                                (max[2] - ray.O[2]) / ray.u[2]};

        auto t_mins = Vector(std::min(t_min_vec[0], t_max_vec[0]),
                             std::min(t_min_vec[1], t_max_vec[1]),
                             std::min(t_min_vec[2], t_max_vec[2]));
        auto t_maxs = Vector(std::max(t_min_vec[0], t_max_vec[0]),
                             std::max(t_min_vec[1], t_max_vec[1]),
                             std::max(t_min_vec[2], t_max_vec[2]));

        double t_0 = std::max(std::max(t_mins[0], t_mins[1]), t_mins[2]);
        double t_1 = std::min(std::min(t_maxs[0], t_maxs[1]), t_maxs[2]);

        return (t_0 <= t_1) && (t_1 >= 0) && (t_0 < t_max) && (t_1 > t_min);
    }
};

struct BVHNode {
    BBox bbox;
    int left;
    int right;
    int tri_start;
    int tri_end;

    BVHNode() : left(-1), right(-1), tri_start(-1), tri_end(-1) {}
    bool is_leaf() const { return left == -1 && right == -1; }
    bool is_hit(const Ray& ray, double t_min, double t_max) const {
        return bbox.is_hit(ray, t_min, t_max);
    }
};

// Class only used in labs 3 and 4
class TriangleMesh : public Object {
   public:
    TriangleMesh(const Vector& albedo, bool mirror = false,
                 bool transparent = false)
        : ::Object(albedo, mirror, transparent) {};

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

        // OBJ indices are 1-based and can be negative (relative), this
        // normalizes them
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
                if (sscanf(s, "v %lf %lf %lf %lf %lf %lf", &pos[0], &pos[1],
                           &pos[2], &col[0], &col[1], &col[2]) == 6) {
                    for (int i = 0; i < 3; i++)
                        col[i] = std::min(1.0, std::max(0.0, col[i]));
                    vertexcolors.push_back(col);
                } else {
                    sscanf(s, "v %lf %lf %lf", &pos[0], &pos[1], &pos[2]);
                }
                vertices.push_back(pos);
            } else if (line[0] == 'f') {
                int i[4], j[4], k[4], offset, nn;
                const char* cur = s + 1;
                TriangleIndices t;
                t.group = curGroup;

                // Try each face format: v/vt/vn, v/vt, v//vn, v
                if ((nn = sscanf(cur, "%d/%d/%d %d/%d/%d %d/%d/%d%n", &i[0],
                                 &j[0], &k[0], &i[1], &j[1], &k[1], &i[2],
                                 &j[2], &k[2], &offset)) == 9) {
                    setFaceVerts(t, i[0], i[1], i[2]);
                    setFaceUVs(t, j[0], j[1], j[2]);
                    setFaceNormals(t, k[0], k[1], k[2]);
                } else if ((nn = sscanf(cur, "%d/%d %d/%d %d/%d%n", &i[0],
                                        &j[0], &i[1], &j[1], &i[2], &j[2],
                                        &offset)) == 6) {
                    setFaceVerts(t, i[0], i[1], i[2]);
                    setFaceUVs(t, j[0], j[1], j[2]);
                } else if ((nn = sscanf(cur, "%d//%d %d//%d %d//%d%n", &i[0],
                                        &k[0], &i[1], &k[1], &i[2], &k[2],
                                        &offset)) == 6) {
                    setFaceVerts(t, i[0], i[1], i[2]);
                    setFaceNormals(t, k[0], k[1], k[2]);
                } else if ((nn = sscanf(cur, "%d %d %d%n", &i[0], &i[1], &i[2],
                                        &offset)) == 3) {
                    setFaceVerts(t, i[0], i[1], i[2]);
                } else
                    continue;

                indices.push_back(t);
                cur += offset;

                // Fan triangulation for polygon faces (4+ vertices)
                while (*cur && *cur != '\n') {
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if ((nn = sscanf(cur, " %d/%d/%d%n", &i[3], &j[3], &k[3],
                                     &offset)) == 3) {
                        setFaceVerts(t2, i[0], i[2], i[3]);
                        setFaceUVs(t2, j[0], j[2], j[3]);
                        setFaceNormals(t2, k[0], k[2], k[3]);
                    } else if ((nn = sscanf(cur, " %d/%d%n", &i[3], &j[3],
                                            &offset)) == 2) {
                        setFaceVerts(t2, i[0], i[2], i[3]);
                        setFaceUVs(t2, j[0], j[2], j[3]);
                    } else if ((nn = sscanf(cur, " %d//%d%n", &i[3], &k[3],
                                            &offset)) == 2) {
                        setFaceVerts(t2, i[0], i[2], i[3]);
                        setFaceNormals(t2, k[0], k[2], k[3]);
                    } else if ((nn = sscanf(cur, " %d%n", &i[3], &offset)) ==
                               1) {
                        setFaceVerts(t2, i[0], i[2], i[3]);
                    } else {
                        cur++;
                        continue;
                    }

                    indices.push_back(t2);
                    cur += offset;
                    i[2] = i[3];
                    j[2] = j[3];
                    k[2] = k[3];
                }
            }
        }
    }

    BBox compute_bbox_simple(int start_triangle, int end_triangle) {
        Vector bbox_min( 1e10,  1e10,  1e10);
        Vector bbox_max(-1e10, -1e10, -1e10);
        for (int i = start_triangle; i < end_triangle; i++) {
            for (int k = 0; k < 3; k++) {
                const Vector& v = vertices[indices[i].vtx[k]];
                for (int axis = 0; axis < 3; axis++) {
                    bbox_min[axis] = std::min(bbox_min[axis], v[axis]);
                    bbox_max[axis] = std::max(bbox_max[axis], v[axis]);
                }
            }
        }
        return BBox{bbox_min, bbox_max};
    }

    BBox compute_bbox(int start_triangle, int end_triangle) {
        auto v0_min =
            (*std::min_element(
                 indices.begin() + start_triangle,
                 indices.begin() + end_triangle,
                 [this](const TriangleIndices& a, const TriangleIndices& b) {
                     return std::min(std::min(vertices[a.vtx[0]][0],
                                              vertices[a.vtx[1]][0]),
                                     vertices[a.vtx[2]][0]) <
                            std::min(std::min(vertices[b.vtx[0]][0],
                                              vertices[b.vtx[1]][0]),
                                     vertices[b.vtx[2]][0]);
                 }))
                .vtx[0];
        auto v0_max =
            (*std::max_element(
                 indices.begin() + start_triangle,
                 indices.begin() + end_triangle,
                 [this](const TriangleIndices& a, const TriangleIndices& b) {
                     return std::max(std::max(vertices[a.vtx[0]][0],
                                              vertices[a.vtx[1]][0]),
                                     vertices[a.vtx[2]][0]) <
                            std::max(std::max(vertices[b.vtx[0]][0],
                                              vertices[b.vtx[1]][0]),
                                     vertices[b.vtx[2]][0]);
                 }))
                .vtx[0];

        auto v1_min =
            (*std::min_element(
                 indices.begin() + start_triangle,
                 indices.begin() + end_triangle,
                 [this](const TriangleIndices& a, const TriangleIndices& b) {
                     return std::min(std::min(vertices[a.vtx[0]][1],
                                              vertices[a.vtx[1]][1]),
                                     vertices[a.vtx[2]][1]) <
                            std::min(std::min(vertices[b.vtx[0]][1],
                                              vertices[b.vtx[1]][1]),
                                     vertices[b.vtx[2]][1]);
                 }))
                .vtx[1];
        auto v1_max =
            (*std::max_element(
                 indices.begin() + start_triangle,
                 indices.begin() + end_triangle,
                 [this](const TriangleIndices& a, const TriangleIndices& b) {
                     return std::max(std::max(vertices[a.vtx[0]][1],
                                              vertices[a.vtx[1]][1]),
                                     vertices[a.vtx[2]][1]) <
                            std::max(std::max(vertices[b.vtx[0]][1],
                                              vertices[b.vtx[1]][1]),
                                     vertices[b.vtx[2]][1]);
                 }))
                .vtx[1];

        auto v2_min =
            (*std::min_element(
                 indices.begin() + start_triangle,
                 indices.begin() + end_triangle,
                 [this](const TriangleIndices& a, const TriangleIndices& b) {
                     return std::min(std::min(vertices[a.vtx[0]][2],
                                              vertices[a.vtx[1]][2]),
                                     vertices[a.vtx[2]][2]) <
                            std::min(std::min(vertices[b.vtx[0]][2],
                                              vertices[b.vtx[1]][2]),
                                     vertices[b.vtx[2]][2]);
                 }))
                .vtx[2];
        auto v2_max =
            (*std::max_element(
                 indices.begin() + start_triangle,
                 indices.begin() + end_triangle,
                 [this](const TriangleIndices& a, const TriangleIndices& b) {
                     return std::max(std::max(vertices[a.vtx[0]][2],
                                              vertices[a.vtx[1]][2]),
                                     vertices[a.vtx[2]][2]) <
                            std::max(std::max(vertices[b.vtx[0]][2],
                                              vertices[b.vtx[1]][2]),
                                     vertices[b.vtx[2]][2]);
                 }))
                .vtx[2];

        Vector bbox_min(vertices[v0_min][0], vertices[v1_min][1],
                        vertices[v2_min][2]);
        Vector bbox_max(vertices[v0_max][0], vertices[v1_max][1],
                        vertices[v2_max][2]);
        return BBox{bbox_min, bbox_max};
    }

    // // void compute_bbox
    // BBox compute_bbox(const std::vector<Vector>::const_iterator& begin,
    //                   const std::vector<Vector>::const_iterator& end) {
    //         auto v0_min =
    //             (*std::min_element(begin, end,
    //                                [](const Vector& a, const Vector& b) {
    //                                    return a[0] < b[0];
    //                                }))[0];
    //         auto v0_max =
    //             (*std::max_element(begin, end,
    //                                [](const Vector& a, const Vector& b) {
    //                                    return a[0] < b[0];
    //                                }))[0];

    //         auto v1_min =
    //             (*std::min_element(begin, end,
    //                                [](const Vector& a, const Vector& b) {
    //                                    return a[1] < b[1];
    //                                }))[1];
    //         auto v1_max =
    //             (*std::max_element(begin, end,
    //                                [](const Vector& a, const Vector& b) {
    //                                    return a[1] < b[1];
    //                                }))[1];

    //         auto v2_min =
    //             (*std::min_element(begin, end,
    //                                [](const Vector& a, const Vector& b) {
    //                                    return a[2] < b[2];
    //                                }))[2];
    //         auto v2_max =
    //             (*std::max_element(begin, end,
    //                                [](const Vector& a, const Vector& b) {
    //                                    return a[2] < b[2];
    //                                }))[2];

    //         Vector bbox_min(v0_min, v1_min, v2_min);
    //         Vector bbox_max(v0_max, v1_max, v2_max);
    //         return BBox{bbox_min, bbox_max};
    // }

    void compute_initial_bbox() {
        if (vertices.empty()) return;
        bbox = compute_bbox_simple(0, indices.size());
        computed_bbox = true;
    }

    void compute_bvh(int starting_triangle, int ending_triangle,
                     int node_index) {
        bvh[node_index].tri_start = starting_triangle;
        bvh[node_index].tri_end = ending_triangle;
        bvh[node_index].bbox = compute_bbox_simple(starting_triangle, ending_triangle);

        // std::cerr << "Node " << node_index << ": triangles ["
        //           << starting_triangle << ", " << ending_triangle
        //           << ") bbox min " << bvh[node_index].bbox.min[0] << ","
        //           << bvh[node_index].bbox.min[1] << ","
        //           << bvh[node_index].bbox.min[2] << " max "
        //           << bvh[node_index].bbox.max[0] << ","
        //           << bvh[node_index].bbox.max[1] << ","
        //           << bvh[node_index].bbox.max[2] << std::endl;

        if (ending_triangle - starting_triangle <= 5) {
            bvh[node_index].left = -1;
            bvh[node_index].right = -1;
            return;
        }

        Vector diag = bvh[node_index].bbox.max - bvh[node_index].bbox.min;
        int axis = 0;
        if (diag[1] > diag[0] && diag[1] > diag[2]) axis = 1;
        if (diag[2] > diag[0] && diag[2] > diag[1]) axis = 2;

        std::sort(
            indices.begin() + starting_triangle,
            indices.begin() + ending_triangle,
            [this, axis](const TriangleIndices& a, const TriangleIndices& b) {
                double center_a = 0, center_b = 0;
                for (int i = 0; i < 3; i++) {
                    center_a += vertices[a.vtx[i]][axis];
                    center_b += vertices[b.vtx[i]][axis];
                }
                return center_a < center_b;
            });

        int mid = (starting_triangle + ending_triangle) / 2;

        int left_idx = bvh.size();
        bvh.push_back(BVHNode{});
        int right_idx = bvh.size();
        bvh.push_back(BVHNode{});

        bvh[node_index].left = left_idx;
        bvh[node_index].right = right_idx;

        compute_bvh(starting_triangle, mid, left_idx);
        compute_bvh(mid, ending_triangle, right_idx);
    }

    void build_bvh() {
        bvh.clear();
        if (indices.empty()) return;
        bvh.push_back(BVHNode{});
        compute_bvh(0, indices.size(), 0);
    }

    bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {
        if (bvh.empty()) return false;
        bool intersected = false;
        t = 1e100;

        std::vector<int> stack = {0};
        while (!stack.empty()) {
            int idx = stack.back();
            stack.pop_back();

            if (!bvh[idx].is_hit(ray, 0, t)) continue;

            if (bvh[idx].is_leaf()) {
                for (int tri_index = bvh[idx].tri_start;
                     tri_index < bvh[idx].tri_end; tri_index++) {
                    auto tri = indices[tri_index];
                    auto v0 = vertices[tri.vtx[0]];
                    auto v1 = vertices[tri.vtx[1]];
                    auto v2 = vertices[tri.vtx[2]];
                    auto edge1 = v1 - v0;
                    auto edge2 = v2 - v0;
                    auto denom = dot(ray.u, cross(edge1, edge2));
                    if (std::abs(denom) < 1e-12) continue;

                    auto beta = dot(edge2, cross((v0 - ray.O), ray.u)) / denom;
                    auto gamma =
                        -1 * dot(edge1, cross((v0 - ray.O), ray.u)) / denom;
                    auto alpha = 1 - beta - gamma;
                    auto tri_t = dot((v0 - ray.O), cross(edge1, edge2)) / denom;

                    if (tri_t > 0 && alpha >= 0 && beta >= 0 && gamma >= 0 &&
                        tri_t < t) {
                        t = tri_t;
                        intersected = true;
                        P = ray.O + t * ray.u;
                        N = cross(edge1, edge2);
                        N.normalize();
                        if (dot(ray.u, N) > 0) N = -1.0 * N;
                    }
                }
            } else {
                stack.push_back(bvh[idx].left);
                stack.push_back(bvh[idx].right);
            }
        }

        return intersected;
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;

    BBox bbox;
    bool computed_bbox;

    std::vector<BVHNode> bvh;
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
                if (tmp_t < min_t) {
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
                auto new_ray = Ray(P + 1e-10 * N, ref_direction);

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
                double cos_i =
                    -dot(ray.u, normal);  // cos of incident angle (>0)
                double sin2_t =
                    sqr(n_ratio) * (1.0 - sqr(cos_i));  // sin²(θ_t) via Snell

                if (sin2_t > 1.0) {
                    // Total internal reflection — treat like a mirror
                    Vector ref_direction = ray.u - 2 * dot(ray.u, N) * N;
                    Ray reflected_ray(P + 1e-10 * N, ref_direction);
                    return getColor(reflected_ray, recursion_depth + 1);
                }

                // Refracted ray direction (Snell's law in vector form)
                double cos_t = sqrt(1.0 - sin2_t);
                Vector refracted_dir =
                    n_ratio * ray.u + (n_ratio * cos_i - cos_t) * normal;
                refracted_dir.normalize();

                // Offset origin inward (against the normal) to avoid
                // self-intersection
                Ray refracted_ray(P - 1e-10 * normal, refracted_dir);
                return getColor(refracted_ray, recursion_depth + 1);
            }

            // return Vector(255, 255, 255);

            // test if there is a shadow by sending a new ray
            Vector direction_to_light = light_position - P;
            auto distance_to_light = direction_to_light.norm();
            direction_to_light.normalize();

            auto shadow_ray = Ray(P + 1e-10 * N, direction_to_light);

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
                direct_light =
                    (light_intensity /
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
            auto new_ray = Ray(P + 1e-10 * N, new_direction);
            indirect_light = objects[object_id]->albedo *
                             getColor(new_ray, recursion_depth + 1);

            return direct_light + indirect_light;
        }

        return Vector(0, 0, 0);
    }

    std::vector<const Object*> objects;

    Vector camera_center, light_position;
    double fov, gamma, light_intensity;
    int max_light_bounce;
};

void boxMuller(double stdev, double& x, double& y) {
    // since MacOS issues with fopenmp, I will just use the first engine for all
    // threads double r1 = uniform ( engine[omp_get_thread_num()] ) ; double r2
    // = uniform ( engine[omp_get_thread_num()] ) ;
    double r1 = uniform(engine[0]);
    double r2 = uniform(engine[0]);
    x = sqrt(-2 * log(r1 + 1e-15)) * cos(2 * M_PI * r2) * stdev;
    y = sqrt(-2 * log(r1 + 1e-15)) * sin(2 * M_PI * r2) * stdev;
}

int main() {
    int W = 512;
    int H = 512;

    for (int i = 0; i < 32; i++) {
        engine[i].seed(i);
    }

    Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8), false,
                         true);
    Sphere center_sphere_2(Vector(-15, 0, 15), 10., Vector(0.8, 0.8, 0.8),
                           true, false);
    Sphere center_sphere_3(Vector(15, 0, -15), 10., Vector(0.8, 0.8, 0.8), false,
                           true);
    Sphere wall_left(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
    Sphere wall_right(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
    Sphere wall_front(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
    Sphere wall_behind(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
    Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3), false);
    Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7), false);

    TriangleMesh cat(Vector(0.6, 0.5, 0.4));
    cat.readOBJ("cat.obj");
    cat.scale_translate(0.6, Vector(0, -5, 0));
    cat.compute_initial_bbox();
    cat.build_bvh();
    // cat.computeBoundingBox();

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

    scene.addObject(&cat);

    int SAMPLES = 32;
    double FOCAL_LENGTH = 40.;
    double APERTURE_SIZE = 1;
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

                Vector focal_point =
                    scene.camera_center + FOCAL_LENGTH * ray_direction;
                boxMuller(APERTURE_SIZE, ax, ay);
                Vector new_origin =
                    scene.camera_center + APERTURE_SIZE * Vector(ax, ay, 0);
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
    stbi_write_png("image_lab5__1.png", W, H, 3, &image[0], 0);

    return 0;
}