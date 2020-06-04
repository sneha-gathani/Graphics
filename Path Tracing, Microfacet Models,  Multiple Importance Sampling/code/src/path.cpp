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

class path : public Integrator {
public:
    std::vector<Mesh *> allmeshes;

    path(const PropertyList &propList) {
        
    }
    //sample the light
    Color3f sample_light(const Scene *scene, Sampler *sampler, const Ray3f &ray, Intersection its, Mesh *meshof_emitter_to_emit, Vector3f sampled_point, Vector3f diff_vector, EmitterQueryRecord E, int number_of_emitter_meshes) const
    {
        //make the ray in the direction that we got before
        Ray3f newray(its.p, diff_vector);
        Intersection next_hit;
        Color3f final_rad(0.0f);

        //check if this new ray intersects the scene
        bool is_hit = scene->rayIntersect(newray, next_hit);

        //if it does not intersect the scene or if the same mesh is hit back (post recursion of the rays)
        if(!is_hit || next_hit.mesh == meshof_emitter_to_emit)
        {
            BSDFQueryRecord bsdfquery(its.shFrame.toLocal(-ray.d), its.shFrame.toLocal(diff_vector), EMeasure::ESolidAngle);

            Color3f bsdf = its.mesh->getBSDF()->eval(bsdfquery);

            //finding costheta value; that is normal at x and ray to emitter
            float costheta = std::abs(its.shFrame.cosTheta(its.shFrame.toLocal(diff_vector).normalized()));
            //finding cosalpha value; that is normal at y and ray to point x
            float cosalpha = std::abs((E.n).dot(-diff_vector.normalized()));
            
            const Emitter *emitter_to_emit = meshof_emitter_to_emit->getEmitter();
            // if(emitter_to_emit == NULL)
            
            float pdf = emitter_to_emit->pdf(E, *meshof_emitter_to_emit);

            Color3f source_radiance_value = emitter_to_emit->eval(E);
            
            float bsdf_pdf = its.mesh->getBSDF()->pdf(bsdfquery);
            //factor of contribution of the light sampled
            float mis_weight;

            ////this is the pdf over angle; change from coordinate system to angles
            float light_pdf = (pdf * diff_vector.norm() * diff_vector.norm()) / cosalpha;
            //balance heuristic for sampling the light
            mis_weight = light_pdf / (light_pdf + bsdf_pdf);

            final_rad = (mis_weight * source_radiance_value * costheta * cosalpha * bsdf * number_of_emitter_meshes) / (pdf * diff_vector.norm() * diff_vector.norm());
        }
        return final_rad;
    }


    //sample the brdf
    Color3f sample_brdf(const Scene *scene, Sampler *sampler, const Ray3f &ray, Intersection its, Mesh *meshof_emitter_to_emit, Vector3f sampled_point, Vector3f diff_vector, EmitterQueryRecord E, int number_of_emitter_meshes) const
    {
        BSDFQueryRecord bsdfquery(its.shFrame.toLocal(-ray.d));

        Color3f sampled_color = its.mesh->getBSDF()->sample(bsdfquery, sampler->next2D());
        float bsdf_pdf = its.mesh->getBSDF()->pdf(bsdfquery);

        Intersection next_hit;
        Ray3f newray(its.p, its.shFrame.toLocal(bsdfquery.wo.normalized()));

        float mis_weight = 0.0f;
        Color3f final_rad(0.0f);
        Color3f source_radiance_value(0.0f);

        //check if this new ray intersects the scene
        bool is_hit = scene->rayIntersect(newray, next_hit);

        if(is_hit && next_hit.mesh->isEmitter())
        {
            //get the direction of the ray from scene intersection point to the new sampled point
            Vector3f diff_vector = next_hit.p - its.p;

            auto *emitter_instance = next_hit.mesh->getEmitter();
            EmitterQueryRecord E1(next_hit.p);
            //pdf of the light source
            float pdf_from_light = emitter_instance->pdf(E1, *next_hit.mesh); 

            source_radiance_value = emitter_instance->eval(E);

            //finding cosalpha value; that is normal at y and ray to point x
            float cosalpha = std::abs(((next_hit.shFrame.n).dot(-diff_vector)) / (diff_vector.norm() * next_hit.shFrame.n.norm()));

            //this is the pdf over angle; change from coordinate system to angles
            float light_pdf = (pdf_from_light * ((diff_vector.norm() * diff_vector.norm()))) / cosalpha;
            //balance heuristic for sampling the brdf
            mis_weight = bsdf_pdf / (bsdf_pdf + light_pdf);
        }
        
        final_rad = mis_weight * sampled_color * source_radiance_value;
        return final_rad;
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

                    //sample the light to get contribution of light in the total radiance
                    Color3f light_contribution = sample_light(scene, sampler, ray, its, meshof_emitter_to_emit, sampled_point, diff_vector, E, number_of_emitter_meshes);
                    //sample the brdf to get contribution of light in the total radiance
                    Color3f brdf_contribution = sample_brdf(scene, sampler, ray, its, meshof_emitter_to_emit, sampled_point, diff_vector, E, number_of_emitter_meshes);
            
                    //check if the emitter is emitting light in the downward direction
                    if((E.n).dot(-diff_vector) > 0.0f)
                    {
                        //total radiance is the radiance + alpha times both the contributions of light and brdf
                        rad = rad + (alpha * (light_contribution + brdf_contribution));
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
        return "path[]";
    }
    
};

NORI_REGISTER_CLASS(path, "path");
NORI_NAMESPACE_END