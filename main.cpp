#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>
#include <map>
#include <string>
#include <fstream>

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
	double norm() const {
		return sqrt(norm2());
	}
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
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

class Ray {
public:
	Ray(const Vector& origin, const Vector& unit_direction) : O(origin), u(unit_direction) {};
	Vector O, u;
};

class Object {
public:
	Object(const Vector& albedo, bool mirror = false, bool transparent = false) : albedo(albedo), mirror(mirror), transparent(transparent) {};

	virtual bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const = 0;

	Vector albedo;
	bool mirror, transparent;
};

class Sphere : public Object {
public:
	Sphere(const Vector& center, double radius, const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent), C(center), R(radius) {};

	// returns true iif there is an intersection between the ray and the sphere
	// if there is an intersection, also computes the point of intersection P, 
	// t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
	// and the unit normal N
	bool intersect(const Ray& ray, Vector& P, double &t, Vector& N) const {
		// TODO (lab 1) : compute the intersection (just true/false at the begining of lab 1, then P, t and N as well)
        
        Vector L = ray.O - C;
        double a = dot(ray.u, ray.u);
        double b = 2.0 * dot(ray.u, L);
        double c = dot(L, L) - R * R;
        double delta = b*b - 4*a*c;
        if(delta < 0) return false;

        double sqrtd = sqrt(delta);
        double t0 = (-b - sqrtd) / (2*a);
        double t1 = (-b + sqrtd) / (2*a);

        const double EPS = 1e-6;
        if(t0 >= EPS) t = t0;
        else if(t1 >= EPS) t = t1;
        else return false;

        P = ray.O + ray.u * t;
        N = P - C;
        N.normalize();
        return true;

		// return false;
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
		// TODO (labs 3 and 4)
		
		// lab 3 : for each triangle, compute the ray-triangle intersection with Moller-Trumbore algorithm
		// lab 3 : once done, speed it up by first checking against the mesh bounding box
		// lab 4 : recursively apply the bounding-box test from a BVH datastructure

		// find bonding box
		Vector bound_mn, bound_mx;
		if(vertices.empty()){
			bound_mn = Vector(0, 0, 0);
			bound_mx = Vector(0, 0, 0);
		}
		bound_mn = vertices[0];
		bound_mx = vertices[0];
		for(auto v:vertices){
			for(int i=0;i<3;i++){
				bound_mn[i] = std::min(bound_mn[i], v[i]);
				bound_mx[i] = std::max(bound_mx[i], v[i]);
			}
		}

		// check if the ray intersect with the bounding box
		bool intersect_bx = true;
		double mn = -1e30, mx = 1e30;
		for(int i=0;i<3;i++){
			if(std::abs(ray.u[0]) < 1e-6){
				// ray parall to the plan
				if(ray.O[i] < bound_mn[i] || ray.O[i] > bound_mx[i]){
					intersect_bx = false;
					break;
				}
			}
			else{
				// B - O / u
				double t0 = (bound_mn[i] - ray.O[i]) / ray.u[i];
				double t1 = (bound_mx[i] - ray.O[i]) / ray.u[i];
				if(t0 > t1){std::swap(t0, t1);}
				mn = std::max(mn, t0);
				mx = std::min(mx, t1);
				if(mn > mx){
					intersect_bx = false;
					break;
				}
			}
		}
		
		if(!intersect_bx){return false;}


		// Moller-Trumbore, iterate through all vertex index to check if the ray intersect with the traingle
		bool intersect_tri = false;
		Vector PP, NN;
		double t_mn = 1e30;
		for(auto idx : indices){
			Vector P_i, N_i;
			double t_i;

			Vector A = vertices[idx.vtx[0]];
			Vector B = vertices[idx.vtx[1]];
			Vector C = vertices[idx.vtx[2]];
			Vector e1 = B - A;
			Vector e2 = C - A;

			Vector h = cross(ray.u, e2);
			double a = dot(e1, h);
			if(std::abs(a) < 1e-6){continue;}

			double f = 1.0 / a;
			Vector s = ray.O - A;
			double u = f * dot(s, h);
			if(u < 0.0 || u > 1.0){continue;}

			Vector q = cross(s, e1);
			double v = f * dot(ray.u, q);
			if(v < 0.0 || u + v > 1.0){continue;}

			t_i = f * dot(e2, q);
			if(t_i < 1e-6){continue;}

			P_i = ray.O + ray.u * t_i;

			if(idx.n[0] != -1 && idx.n[1] != -1 && idx.n[2] != -1){
				double w = 1.0 - u - v;
				N_i = w * normals[idx.n[0]] + u * normals[idx.n[1]] + v * normals[idx.n[2]];
			}
			else{
				N_i = cross(e1, e2);
			}
			N_i.normalize();

			if(t_i < t_mn){
				intersect_tri = true;
				t_mn = t_i;
				PP = P_i;
				NN = N_i;
			}
		}

		if(intersect_tri){
			P = PP;
			N = NN;
			N.normalize();
			t = t_mn;
		}

		return intersect_tri;
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
	void addObject(const Object* obj) {
		objects.push_back(obj);
	}

	// returns true iif there is an intersection between the ray and any object in the scene
    // if there is an intersection, also computes the point of the *nearest* intersection P, 
    // t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
    // and the unit normal N. 
	// Also returns the index of the object within the std::vector objects in object_id
	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N, int &object_id) const  {

		// TODO (lab 1): iterate through the objects and check the intersections with all of them, 
		// and keep the closest intersection, i.e., the one if smallest positive value of t

        bool hit = false;
        double t_mn = 1e30;
        Vector PP, NN;
        double tt;
        int id;

        // iterate through all obj to find the closest obj
        for (int i = 0; i < (int)objects.size(); i++) {
            if (objects[i]->intersect(ray, PP, tt, NN)) {
                if (tt > 1e-6 && tt < t_mn) {
                    hit = true;
                    t_mn = tt;
                    id = i;
                    P = PP;
                    N = NN;
                }
            }
        }

        if (hit) {
            t = t_mn;
            object_id = id;
            N.normalize();
            return true;
        }

		return false;
	}


	// return the radiance (color) along ray
	Vector getColor(const Ray& ray, int recursion_depth) {

		if (recursion_depth >= max_light_bounce) return Vector(0, 0, 0);

		// TODO (lab 1) : if intersect with ray, use the returned information to compute the color ; otherwise black 
		// in lab 1, the color only includes direct lighting with shadows

		Vector P, N;
		double t;
		int object_id;
		if (intersect(ray, P, t, N, object_id)) {

			if (objects[object_id]->mirror) {
                Vector r = ray.u - 2.0 * dot(ray.u, N) * N;
                r.normalize();
                return getColor(Ray(P + N * 1e-6, r), recursion_depth + 1);
				// return getColor in the reflected direction, with recursion_depth+1 (recursively)
			} // else

			if (objects[object_id]->transparent) { // optional

				// return getColor in the refraction direction, with recursion_depth+1 (recursively)
			} // else

			// test if there is a shadow by sending a new ray
			// if there is no shadow, compute the formula with dot products etc.
            
            Vector l = this->light_position - P;
            double dist2 = dot(l, l);
            double dist = sqrt(dist2);
            l = l / dist;   // normalize
            Ray shadow_ray(P + N * 1e-6, l);
            bool is_shadow = false;
            for (int k = 0; k < (int)objects.size(); k++) {
                if (k == object_id) continue;
                Vector P_shadow, N_shadow;
                double t_shadow;
                if (objects[k]->intersect(shadow_ray, P_shadow, t_shadow, N_shadow)) {
                    if (t_shadow > 1e-6 && t_shadow < dist) {
                        is_shadow = true;
                        break;
                    }
                }
            }

            // if(is_shadow){return Vector(0, 0, 0);}

            double n_dot_l = std::max(0.0, dot(N, l));
            double irr = light_intensity / (4.0 * M_PI * dist2);
            Vector color = is_shadow? Vector(0, 0, 0) : (objects[object_id]->albedo / M_PI) * irr * std::max(0.0, dot(N, l));
			
			// TODO (lab 2) : add indirect lighting component with a recursive call
			Vector indirect_light(0, 0, 0);
			int thread_id = omp_get_thread_num();

			double r1 = uniform(engine[thread_id]);
			double r2 = uniform(engine[thread_id]);
			double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
			double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
			double z = sqrt(r2);
			Vector t1 = Vector(-N[1], N[0], 0);
			t1.normalize();
			Vector t2 = cross(N, t1);
			Vector dir = x*t1 + y*t2 + z*N;
			dir.normalize();
			Ray indirect_ray(P + N * 1e-6, dir);
			Vector indirect_color = getColor(indirect_ray, recursion_depth + 1);
			// Vector albedo = objects[object_id]->albedo;
			indirect_light = indirect_light + (objects[object_id]->albedo / M_PI) * indirect_color;
			return color + indirect_light;
		}

		

		return Vector(0, 0, 0);
	}

	std::vector<const Object*> objects;

	Vector camera_center, light_position;
	double fov, gamma, light_intensity;
	int max_light_bounce;
};


int main() {
	int W = 512;
	int H = 512;

	for (int i = 0; i<32; i++) {
		engine[i].seed(i);
	}

	Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8), true);
	Sphere front_sphere(Vector(-15, -3, 20), 6., Vector(0.8, 0.8, 0.8), false);
	Sphere back_sphere(Vector(25, -2, -15), 8., Vector(0.8, 0.8, 0.8), true);
	Sphere wall_left(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
	Sphere wall_right(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
	Sphere wall_front(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
	Sphere wall_behind(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3));
	Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));

	Scene scene;
	scene.camera_center = Vector(0, 0, 55);
	scene.light_position = Vector(-10,20,40);
	// scene.light_intensity = 2E6;
	// scene.light_intensity = 1E5;
	scene.light_intensity = 3E7;
	scene.fov = 60 * M_PI / 180.;
	// scene.gamma = 1.0;    // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
	scene.gamma = 2.2;    // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
	scene.max_light_bounce = 5;

	// scene.addObject(&center_sphere);
	// scene.addObject(&front_sphere);
	// scene.addObject(&back_sphere);

	TriangleMesh* cat = new TriangleMesh(Vector(0.8, 0.8, 0.8), false, false);
	cat->readOBJ("cadnav.com_model/Models_F0202A090/cat.obj");
	cat->scale_translate(0.6, Vector(0, -10, 0));
	scene.addObject(cat);

	// /*
	scene.addObject(&wall_left);
	scene.addObject(&wall_right);
	scene.addObject(&wall_front);
	scene.addObject(&wall_behind);
	scene.addObject(&ceiling);
	scene.addObject(&floor);
	// */

	std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector color;

			// TODO (lab 1) : correct ray_direction so that it goes through each pixel (j, i)			
			// Vector ray_direction(0., 0., -1);
            double aspect_ratio = W / (double)H;
			double fov_rad = scene.fov;
			double tan_fov = scene.fov;
			double d = W / (2. * std::tan(scene.fov / 2.));
			double x = j + 0.5 - W / 2;
			double y = (H - i - 1) + 0.5 - H / 2.;
			Vector ray_direction(x, y, -d);
			ray_direction.normalize();

			Ray ray(scene.camera_center, ray_direction);

            // TODO (lab 2) : add Monte Carlo / averaging of random ray contributions here
            // TODO (lab 2) : add antialiasing by altering the ray_direction here
            // TODO (lab 2) : add depth of field effect by altering the ray origin (and direction) here
            
            int n_samples = 64;
            color = Vector(0, 0, 0);
			int thread_id = omp_get_thread_num();
            
            for (int s = 0; s < n_samples; s++) {
				// antialiasing by altering the ray_direction
				double mu_x = j + 0.5 - W / 2.0;
				double mu_y = (H - i - 1) + 0.5 - H / 2.0;
				// double r1 = uniform(engine[0]);
				// double r2 = uniform(engine[0]);
				double r1 = uniform(engine[thread_id]);
				double r2 = uniform(engine[thread_id]);
				double sigma = 0.5;
				double xx = sigma * sqrt(-2.0 * log(r1)) * cos(2*M_PI*r1);
				double yy = sigma * sqrt(-2.0 * log(r1)) * sin(2*M_PI*r1);

				Vector ray_dir(mu_x + xx, mu_y + yy, -d);
				ray_dir.normalize();

                // color = color + scene.getColor(Ray(scene.camera_center, ray_dir), 0);		// color without depth of field

				// depth of field
				double focal_dist = 55.0;
				double aperture = 1.5;
				Vector focus_point = scene.camera_center + ray_dir * (focal_dist / fabs(ray_dir[2]));

				double r_ap = aperture * sqrt(r1);
				double theta_ap = 2.0 * M_PI * r2;
				Vector origin_dof = scene.camera_center + Vector(r_ap * cos(theta_ap), r_ap * sin(theta_ap), 0.0);

				Vector dir_dof = focus_point - origin_dof;
				dir_dof.normalize();

				color = color + scene.getColor(Ray(origin_dof, dir_dof), 0);
            }
            
            // averaging of random ray
            color = color / n_samples;

			image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., 255. * std::pow(color[0] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., 255. * std::pow(color[1] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., 255. * std::pow(color[2] / 255., 1. / scene.gamma)));
		}
	}
	stbi_write_png("lab4.png", W, H, 3, &image[0], 0);
	delete cat;

	return 0;
}

// g++ -O3 -fopenmp main.cpp -o raytracer
// git fetch upstream
// git merge upstream/master
// git push origin main

/*
image
lab1: rendering
lab1_basic: show the sphere
lab1_mirror: with mirror effect

lab2: antialiasing + depth of field (DoF)
lab2_1: n_sample = 16, sigma = 0.5
lab2_2: n_sample = 32, sigma = 0.5
lab2_3: n_sample = 64, sigma = 0.5
lab2_4: n_sample = 64, sigma = 0.1
lab2_5: n_sample = 64, sigma = 3.0
lab2_threeSphere: three sphere in the image
*/