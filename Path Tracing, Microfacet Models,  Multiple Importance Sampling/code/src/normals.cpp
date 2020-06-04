// #include <nori/integrator.h>

// NORI_NAMESPACE_BEGIN

// class NormalIntegrator : public Integrator {
// public:
//     NormalIntegrator(const PropertyList &props) {
//         m_myProperty = props.getString("myProperty");
//         std::cout << "Parameter value was : " << m_myProperty << std::endl;
//     }

//     /// Compute the radiance value for a given ray. Just return green here
//     Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
//         return Color3f(0, 1, 0);
//     }

//     /// Return a human-readable description for debugging purposes
//     std::string toString() const {
//         return tfm::format(
//             "NormalIntegrator[\n"
//             "  myProperty = \"%s\"\n"
//             "]",
//             m_myProperty
//         );
//     }
// protected:
//     std::string m_myProperty;
// };

// NORI_REGISTER_CLASS(NormalIntegrator, "normals");
// NORI_NAMESPACE_END


#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class NormalIntegrator : public Integrator {
public:
    NormalIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        /* Return the component-wise absolute
           value of the shading normal as a color */
        Normal3f n = its.shFrame.n.cwiseAbs();
        return Color3f(n.x(), n.y(), n.z());
    }

    std::string toString() const {
        return "NormalIntegrator[]";
    }
};

NORI_REGISTER_CLASS(NormalIntegrator, "normals");
NORI_NAMESPACE_END