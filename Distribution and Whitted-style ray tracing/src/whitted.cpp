#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/common.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
public:
    std::vector<Mesh *> allmeshes;

    WhittedIntegrator(const PropertyList &propList) {
        
    }

    //Computing the radiance
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        
        Color3f rad;
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        //Intersection point on scene where we want to calculate radiance
        Point3f x = its.p;
        const Mesh *intersectedmesh = its.mesh; //mesh of the scene
        const BSDF *bsdf = intersectedmesh->getBSDF(); //bsdf of the scene

        //checking if its a diffuse material
        if(bsdf->isDiffuse())
        {
            //get all the meshes of the emitter in a vector
            std::vector<Mesh *> allmeshes = scene->getEmittermeshes();

            //sample randomly to get one of the meshes to be emitted from
            int numberofmeshes = allmeshes.size();
            float randomemittersample = (sampler->next1D()) * numberofmeshes;
            int emittertosamplefrom = (int)(randomemittersample);

            //get the mesh of this emittertosamplefrom
            Mesh *emittertoemit = allmeshes[emittertosamplefrom];

            //get the emitter instance of the emitter we have chosen
            Emitter *emittermesh = emittertoemit->getEmitter();

            //Make a structure of emitter having origin at x  
            EmitterQueryRecord EQR(x);

            //sample a point on light source
            Point3f y = emittermesh->sample(EQR, sampler->next2D(), *emittertoemit);
            //get the radiance
            Color3f emitterradiance = emittermesh->eval(EQR);
            //get the pdf
            float pdf = emittermesh->pdf(EQR, *emittertoemit);
          
            //checking visibility of y and x
            Vector3f dir = y-x;
            Vector3f oppdir = x-y;
            Ray3f shadowRay = Ray3f(x, dir);
            Point3f lightnormal = EQR.n;

            //checks for emitter point sampled normal to always be in downward direction
            float check = lightnormal.dot(oppdir.normalized());
            //if it is in upward direction; then ignore the ray
            if(check < 0.0f)
                return Color3f(0.0f);
            else
            {
                Intersection its1;
                int v = 0;
                if(scene->rayIntersect(shadowRay, its1))
                {
                    if(its1.mesh == emittertoemit)
                        v = 1;
                }

               
                //finding costheta value; that is normal at x and ray to emitter
                float costheta = std::abs(its.shFrame.cosTheta(its.shFrame.toLocal(dir).normalized()));
                //finding cosalpha value; that is normal at y and ray to point x
                float cosalpha = std::abs(lightnormal.dot(oppdir.normalized()));
                //Calculating the bsdf

                //getting the origin point of ray from camera from w_o
                Vector3f w_o = its.toLocal((ray.o - x).normalized());
                BSDFQueryRecord bsdfquery(its.shFrame.toLocal(dir), w_o, EMeasure::ESolidAngle);
                //get bsdf for the scene at point x
                Color3f brdfofmesh = bsdf->eval(bsdfquery); 

                float r = (brdfofmesh.r() * costheta * cosalpha * emitterradiance.r() * v)/(pdf * dir.norm() * dir.norm());
                float g = (brdfofmesh.g() * costheta * cosalpha * emitterradiance.g() * v)/(pdf * dir.norm() * dir.norm());
                float b = (brdfofmesh.b() * costheta * cosalpha * emitterradiance.b() * v)/(pdf * dir.norm() * dir.norm());

                return Color3f(r, g, b);
            }
        }
        //if the material bsdf is dielectric or not diffuse
        else
        {        
            //getting the origin point of ray from camera from w_o
            Vector3f w_o = its.toLocal((ray.o-x).normalized());
            BSDFQueryRecord bsdfquery(w_o);
            //get bsdf for scene at point x
            //return back the bsdf and also the w_t of the transmitted ray into the medium
            Color3f brdfofmesh = bsdf->sample(bsdfquery, sampler->next2D());
            
            float addsampler = sampler->next1D();
            Ray3f newray = Ray3f(x, its.toWorld(bsdfquery.wo));
            if(addsampler < 0.95f)
            {
                rad = (1.0f/0.95f) * (brdfofmesh) * Li(scene, sampler, newray);
                return rad;
            }
            else
            {
                return Color3f(0.0f);
            }
        }
    }

    std::string toString() const {
        return "WhittedIntegrator[]";
    }
    
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END