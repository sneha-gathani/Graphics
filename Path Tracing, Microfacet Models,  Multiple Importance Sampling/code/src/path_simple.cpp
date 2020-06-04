#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/common.h>
#include <nori/bsdf.h>
#include <math.h>
#include <algorithm>

NORI_NAMESPACE_BEGIN

class path_simple : public Integrator {
public:
    std::vector<Mesh *> allmeshes;

    path_simple(const PropertyList &propList) {
        
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
    {

        //final radiance at pixel N
        Color3f rad;
        //intersection point of the ray with the scene
        Intersection its;
        //factor by which radiance will be affected
        Color3f alpha(1.0f);

        Ray3f recursiveray = ray;
        bool is_specular = false;
        float RR = -1.0f;

        while(true)
        {
            //if the ray intersects the scene
            if(scene->rayIntersect(recursiveray, its))
            {
                //get all the emitter meshes in a vector
                std::vector<Mesh *> all_emitter_meshes = scene->getEmittermeshes();
                //number of emitter meshes
                int number_of_emitter_meshes = all_emitter_meshes.size();
            
                //select a random emitter
                int emitter_index_tosamplefrom = (int)(sampler->next1D() * number_of_emitter_meshes);
                Mesh *meshof_emitter_to_emit = all_emitter_meshes[emitter_index_tosamplefrom];
                Emitter *emitter_to_emit = meshof_emitter_to_emit->getEmitter();
                EmitterQueryRecord E(recursiveray.o, its.p, its.shFrame.n);

                //get a sampled point on that emitter
                Point3f sampled_point = emitter_to_emit->sample(E, sampler->next2D(), *meshof_emitter_to_emit);

                //get the direction of the ray from scene intersection point to the new sampled point
                Vector3f diff_vector = sampled_point - its.p;

                //if direct illumination or last intersection is on emitter and one before last is a specular material
                if(its.mesh->isEmitter() && (RR == -1.0f || is_specular))
                {
                    rad = rad + (alpha * emitter_to_emit->eval(E));
                }

                //if intersected mesh is diffuse
                if(its.mesh->getBSDF()->isDiffuse())
                {
                    //make the ray in the direction that we got before
                    Ray3f newray(its.p, diff_vector);
                    Intersection next_hit;

                    //check if this new ray intersects the scene
                    bool is_hit = scene->rayIntersect(newray, next_hit);

                    //if it does not intersect the scene or if the same mesh is hit back (post recursion of the rays)
                    if(!is_hit || next_hit.mesh == meshof_emitter_to_emit)
                    {
                        BSDFQueryRecord bsdfquery(its.shFrame.toLocal(-recursiveray.d), its.shFrame.toLocal(diff_vector), EMeasure::ESolidAngle);

                        Color3f bsdf = its.mesh->getBSDF()->eval(bsdfquery);

                        //finding costheta value; that is normal at x and ray to emitter
                        float costheta = std::abs(its.shFrame.cosTheta(its.shFrame.toLocal(diff_vector).normalized()));
                        //finding cosalpha value; that is normal at y and ray to point x
                        float cosalpha = std::abs((E.n).dot(-diff_vector.normalized()));

                        float pdf = emitter_to_emit->pdf(sampled_point, *meshof_emitter_to_emit);

                        Color3f source_radiance_value = emitter_to_emit->eval(E);

                        //check if the emitter is emitting light in the downward direction
                        if((E.n).dot(-diff_vector) > 0.0f)
                        {
                            rad = rad + (alpha * source_radiance_value * costheta * cosalpha * bsdf * number_of_emitter_meshes) / (pdf * diff_vector.norm() * diff_vector.norm());
                        }
                    }
                }

                //sample the bsdf
                BSDFQueryRecord bsdfquery(its.shFrame.toLocal(-recursiveray.d));

                auto sampled_color = its.mesh->getBSDF()->sample(bsdfquery, sampler->next2D());
                
                //if RR is less than 0.95 then this ray becomes the new ray and alpha is factored by the sampled color
                if(RR < 0.95)
                {
                    Ray3f ray_in_new_direction(its.p, its.shFrame.toWorld(bsdfquery.wo));
                    alpha = (alpha * sampled_color) / 0.95;
                    recursiveray = ray_in_new_direction;
                }
                //termination condition using RR
                else
                    break;
            }
            //if ray does not intersect the scene, break out from the loop
            else
                break;
        //sample a new RR
        RR = sampler->next1D();
        //check if the intersected mesh is specular or not
        is_specular = !its.mesh->getBSDF()->isDiffuse();
    }
    //return the final radiance
    return rad;
}

    std::string toString() const {
        return "path_simple[]";
    }
    
};

NORI_REGISTER_CLASS(path_simple, "path_simple");
NORI_NAMESPACE_END