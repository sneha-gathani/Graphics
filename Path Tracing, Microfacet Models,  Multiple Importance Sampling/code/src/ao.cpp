#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class NormalIntegrator : public Integrator {
public:
    NormalIntegrator(const PropertyList &propList) {
        /* Parameters from the ajax-simple.xml method */
         number_samples = propList.getInteger("number_samples", 4);
         //sampler = propList.getSampler("squareToCosineHemisphere", warp::squareToCosineHemisphere(Sampler));      
    }

    //Computing the radiance
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */

        //Initial radiance set to 0
        float rad = 0.0f;

        //Checking for initial intersection of ray and the scene to get the point on which we want to find radiance
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        //If there is an intersection that is found, its gets updated, hence we have to check for intersection with the scene

        //its is the point where intersection has happened
        Point3f origin = its.p; 

        //Samples a 2D random point using the 2D finction of the sampler.h
        Point2f sample = sampler->next2D();
        //std::cout << "sample" << sample << endl;
        //Manual random generator used to generate samples
        //float u[2];
        // u[0] = (float) rand()/RAND_MAX;
        // u[1] = (float) rand()/RAND_MAX;
        // Point2f sample;
        // sample.x() = u[0];
        // sample.y() = u[1];

        //w is the direction of all the diffused rays from the intersection point its
        Vector3f w = Warp::squareToCosineHemisphere(sample);
        //std::cout << "w" << w.x() << " " << w.y() << " " << w.z() << endl;
        //converting the directions to the world coordinates
        Vector3f worldw = its.shFrame.toWorld(w);

        //Making temporary ray from origin its to direction worldw
        Ray3f tempray = Ray3f(origin, worldw);

        Intersection its1;
        // whether ray from origin in the direction worldw; that is shadow ray is intersected by scene
        // if intersection not found, it reaches the source light and hence v = 1
        int v = 0;
        if(!scene->rayIntersect(tempray,its1))
        {
            v = 1;
        }

        //calculating the radiance using formula
        rad = rad + float(v);
        //return the radiance that is found
        return Color3f(rad, rad, rad);

    }
    std::string toString() const {
        return "NormalIntegrator[]";
    }
    private:
    int number_samples;
};

NORI_REGISTER_CLASS(NormalIntegrator, "ao");
NORI_NAMESPACE_END