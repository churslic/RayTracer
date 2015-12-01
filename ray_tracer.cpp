/**
 * ray_tracer.cpp
 * CS230
 * -------------------------------
 * Implement ray tracer here.
 */

#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))

#include "ray_tracer.h"

using namespace std;

const double Object::small_t=1e-6;

//--------------------
// utility functions
//--------------------
double sqr(const double x)
{
    return x*x;
}

Pixel Pixel_Color(const Vector_3D<double>& color)
{
    Pixel pixel=0;
    SET_RED(pixel,(unsigned char)(min(color.x,1.0)*255));
    SET_GREEN(pixel,(unsigned char)(min(color.y,1.0)*255));
    SET_BLUE(pixel,(unsigned char)(min(color.z,1.0)*255));
    return pixel;
}
//-------------------
// Shader
//-------------------
Vector_3D<double> Phong_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,
              const Vector_3D<double>& intersection_point,
              const Vector_3D<double>& same_side_normal) const
{
    Vector_3D<double> color;

    //First iterate over the vector of lights
    for(unsigned int i = 0; i<world.lights.size(); ++i) {
        //Create direction vector from intersection point to light source
        Vector_3D<double> dirLight =
            intersection_point - world.lights[i]->position;
        //Create endpoint for the ray
        Vector_3D<double> end = intersection_point;

        //Create ray
        Ray shadow_ray(end, dirLight);

        //Check for intersections with other objects
        for(unsigned int j = 0; j<world.objects.size(); ++j) {
            if(world.objects[i]->Intersection(shadow_ray)) {
                world.enable_shadows = true;
                break;
            }
        }

        //If there was an intersection, change the value of the color
        if(world.enable_shadows) {
            Vector_3D<double> lightColor = world.lights[i]->Emitted_Light(ray);

            //Create diffuse component
            Vector_3D<double> diffuse;
            Vector_3D<double> toLight = world.lights[i]->position
                - intersection_point;
            toLight.Normalize();
            double lDotNormal =
                Vector_3D<double>::Dot_Product(same_side_normal,
                toLight);

            diffuse = lightColor * lDotNormal * color_diffuse;

            Vector_3D<double> phong;
            Vector_3D<double> ambient;
        }
    }

    return color;
}

Vector_3D<double> Reflective_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,
              const Vector_3D<double>& intersection_point,
              const Vector_3D<double>& same_side_normal) const
{
    Vector_3D<double> color;

    // TODO: determine the color

    return color;
}

Vector_3D<double> Flat_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,
              const Vector_3D<double>& intersection_point,
              const Vector_3D<double>& same_side_normal) const
{
    return color;
}

//------------------
// Objects
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and
// semi_infinite as appropriate and return true
//------------------
bool Sphere::
Intersection(Ray& ray) const
{
    // TODO
    // p(t) = s + td, s is start point, d is direction
    // (d dot d)t^2 + (2s dot d)t + (s dot s - r^2) = 0

    //Here is the startpoint vector a, center provided from
    //Sphere class
    Vector_3D<double> startpoint = ray.endpoint - center;

    //Find a, b, c
    double a = Vector_3D<double>::Dot_Product(ray.direction,
               ray.direction);
    double b = 2 * Vector_3D<double>::Dot_Product(startpoint,
               ray.direction);
    double c = Vector_3D<double>::Dot_Product(startpoint,
               startpoint) - pow(radius, 2);

    double delta = b * b - (4*a*c);

    if(delta < 0) return false;
    //One intersection
    else if(delta == 0) {
        double t = -b / (2*a);

        if(ray.semi_infinite == true && t > small_t) {
            ray.current_object = this;
            ray.semi_infinite = false;
            ray.t_max = t;
            return true;
        }
        else if(ray.semi_infinite == false && t > small_t
                && t < ray.t_max) {
            ray.current_object = this;
            ray.t_max = t;
            return true;
        }
    }
    //Two intersections
    else if(delta > 0) {
        double minus_t = (-b - sqrt(delta)) / 2.0*a;
        double plus_t = (-b + sqrt(delta)) / 2.0*a;
        double smallest_t = 0;

        //Take the smallest of the t's
        //if both are positive
        if(minus_t >= 0 && plus_t >= 0) smallest_t = minus_t;
        //if one is positive, it must be the plus_t since it
        //is always larger
        else if(plus_t >= 0) smallest_t = plus_t;

        if(ray.semi_infinite == true && smallest_t > small_t) {
            ray.current_object = this;
            ray.semi_infinite = false;
            ray.t_max = smallest_t;
            return true;
        }
        else if(ray.semi_infinite == false && smallest_t > small_t
                && smallest_t < ray.t_max) {
            ray.current_object = this;
            ray.t_max = smallest_t;
            return true;
        }
    }
    return false;
}

Vector_3D<double> Sphere::
Normal(const Vector_3D<double>& location) const
{
    Vector_3D<double> normal;
    // TODO: set the normal
    normal = location - center;
    normal.Normalize();
    return normal;
}

// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and
// semi_infinite as appropriate and return true
bool Plane::
Intersection(Ray& ray) const
{
    // TODO
    // p(t) = s + td
    //Using equation t = (N dot (p0 - s)) / (N dot d)
    Vector_3D<double> v = x1 - ray.endpoint;
    double a = Vector_3D<double>::Dot_Product(normal, v);
    double b = Vector_3D<double>::Dot_Product(normal, ray.direction);
    double t = a / b;
    if(t > small_t) {
        ray.current_object = this;
        ray.t_max = t;
        ray.semi_infinite = false;
        return true;
    }
    return false;
}

Vector_3D<double> Plane::
Normal(const Vector_3D<double>& location) const
{
    return normal;
}

//-----------------
// Camera
// Find the world position of the input pixel
//-----------------
Vector_3D<double> Camera::
World_Position(const Vector_2D<int>& pixel_index)
{
    Vector_3D<double> result;
    //pixel will already be incremented
    Vector_2D<double> pixel = film.pixel_grid.X(pixel_index);
    result = horizontal_vector*pixel.x + vertical_vector*pixel.y
             + focal_point;
    return result;
}

//-----------------
// Render_World
// Find the closest object of intersection and return a pointer to
// it if the ray intersects with an object, then ray.t_max,
// ray.current_object, and ray.semi_infinite will be set
// appropriately if there is no intersection do not modify the ray
// and return 0
//----------------
const Object* Render_World::
Closest_Intersection(Ray& ray)
{
    // TODO
    Ray temp = ray;
    for(unsigned int i = 0; i < objects.size(); ++i) {
        //We find an intersection
        if(objects[i]->Intersection(temp)) {
            //Want to make sure that the intersection has the
            //smallest t value, appropriately update
            if(temp.t_max < ray.t_max) {
                ray.t_max = temp.t_max;
                ray.current_object = temp.current_object;
                ray.semi_infinite = temp.semi_infinite;
            }
        }
    }
    return ray.current_object;
}

// set up the initial view ray and call
void Render_World::
Render_Pixel(const Vector_2D<int>& pixel_index)
{
    // TODO
    // Render_World has Camera camera
    // first parameter is the endpoint, so camera.position
    // second parameter is the direciton, so we need the
    // world position - camera.position
    Vector_3D<double> direction;
    direction = camera.World_Position(pixel_index)
                - camera.position;
    Ray ray(camera.position, direction);
    ray.t_max = FLT_MAX;

    Ray dummy_root;
    Vector_3D<double> color=Cast_Ray(ray,dummy_root);
    camera.film.Set_Pixel(pixel_index,Pixel_Color(color));
}

// cast ray and return the color of the closest intersected surface
// point, or the background color if there is no object intersection
Vector_3D<double> Render_World::
Cast_Ray(Ray& ray,const Ray& parent_ray)
{
    // TODO
    Vector_3D<double> color;
    const Object *closest = Closest_Intersection(ray);
    // determine the color here
    if(closest != NULL) {
        //Need intersection point and the normal
        Vector_3D<double> intersection;
        Vector_3D<double> obj_normal;

        //To find intersection, use the ray equation:
        //p(t) = a + tu
        intersection = ray.endpoint + ray.direction*ray.t_max;

        //To find the normal to the object, call the normal function
        //for that object
        obj_normal = closest->Normal(intersection);
        obj_normal.Normalize();

        //Now set color appropriately
        return closest->material_shader->Shade_Surface(ray,
                *closest, intersection, obj_normal);
    }

    return color;
}
