#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

class NormalIntegrator : public Integrator {
public:
    NormalIntegrator(const PropertyList &propList) {
        /* Parameters from the ajax-simple.xml method */
         fi = propList.getColor("energy", Color3f(3000, 3000, 3000));
         p = propList.getPoint("position", Point3f(-20, 40, 20));
        
    }

    //Computing the radiance
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */

        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        //If there is an intersection that is found, its gets updated, hence we have to check for intersection with the scene

        //x is the point at which you want to find the radiance
        Point3f x = its.p;

        //p is the position of the light source given to us
        Point3f p = this->p;

        //dir is the vector between light source and point at which we want to find the intersection
        Vector3f dir = p-x;

        //unit_dir is the magnitude of this dir vector
        Vector3f unit_dir = dir/dir.norm();

        //theta is the cos of the angle between incident ray from the point p to the light ray and the normal
        //our coordinates are in the world coordinate system and we want to use the local coordinate system 
        float theta = its.shFrame.cosTheta(its.shFrame.toLocal(unit_dir));
        //std::cout<<"theta: "<<theta;

        //tempray is the ray we make; (shadow ray) that goes from point at which we want to find the radiance in the direction of light source
        //this acts like a shadow ray
        Ray3f tempray = Ray3f(x, dir);

        Intersection its1;

        //whether point x and p are visible to each other by passing shadow ray and checking for intersection.
        //if intersection not found, they are visible to each other and hence v = 1
        int v = 0;
        if (!scene->rayIntersect(tempray,its1))
        {
            v = 1;
        }

        //calculate the radiance using the formula
        float rad = (this->fi.x() * INV_PI * INV_PI * v)/4;
        rad = rad * fmax(0, theta)/(dir.norm()*dir.norm());

        //return the radiance that is found
        return Color3f(rad, rad, rad);

    }

    std::string toString() const {
        return "NormalIntegrator[]";
    }
    private:
    Color3f fi;
    Point3f p;  
};

NORI_REGISTER_CLASS(NormalIntegrator, "simple");
NORI_NAMESPACE_END