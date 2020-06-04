/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/common.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    float smithmasking(float b) const {

        if(b < 1.6)
            return (((3.535 * b) + (2.181 * b * b)) / (1.f + (2.276 * b) + (2.577 * b * b)));
        else
            return 1.0f;
    }

    // Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {

        Vector3f wi = bRec.wi.normalized();
        Vector3f wo = bRec.wo.normalized();

        if(Frame::cosTheta(wo) <= 0.0f)
            return Color3f(0.0f);

        Vector3f wh = (wi + wo).normalized();

        Vector3f n = Vector3f(0.0f, 0.0f, 0.1f);

        auto b = 0.0f, c = 0.0f;
        float masking_incident = 0.0f, masking_out = 0.0f;
        float cosine_theta_v, tan_theta_v;
        float multiplying_factor = 0.0f;
        
        float cosThetai = (wh).dot(wi);
        float fresnelterm = fresnel(cosThetai, m_extIOR, m_intIOR);
        float D = Warp::squareToBeckmannPdf(wh, m_alpha);

        //for incoming direction - masking in
        if((wi.dot(wh) / wi.dot(n)) > 0.0f)
            c = 1.0f;
        else
            c = 0.0f;

        cosine_theta_v = wi.dot(n) / (wi.norm() * n.norm());
        tan_theta_v = (1 - (cosine_theta_v * cosine_theta_v)) / (cosine_theta_v * cosine_theta_v);
        b = 1 / (m_alpha * tan_theta_v);

        multiplying_factor = smithmasking(b);

        masking_incident = c * multiplying_factor; //g1

        //for outgoing direction - masking out
        if((wo.dot(wh) / wo.dot(n)) > 0.0f)
            c = 1.f;
        else
            c = 0.f;

        cosine_theta_v = wo.dot(n) / (wo.norm() * n.norm());
        tan_theta_v = (1 - (cosine_theta_v * cosine_theta_v)) / (cosine_theta_v * cosine_theta_v);

        b = 1 / (m_alpha * tan_theta_v);

        multiplying_factor = smithmasking(b);

        masking_out = c * multiplying_factor; //g2

        //visibility = g1 * g2
        auto part1 = m_kd / M_PI;
        auto part2 = (m_ks * D * fresnelterm * masking_incident * masking_out) / (4 * Frame::cosTheta(wi) * Frame::cosTheta(wo) * Frame::cosTheta(wh));
        return (part1 + part2);
    }


    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        Vector3f wi = bRec.wi.normalized();
        Vector3f wo = bRec.wo.normalized();

        if(Frame::cosTheta(wo) <= 0.0f)
            return 0.0f;

        Vector3f wh = (wi + wo).normalized();

        float jh = 1.0f / (4.0f * (wh.dot(wo)));
        float D = Warp::squareToBeckmannPdf(wh, m_alpha);
        float pdf = (m_ks * D * jh) + (((1.0f - m_ks) * Frame::cosTheta(wo)) / M_PI);
        return pdf;
    }


    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        
        Color3f sampledbsdf;
        if (Frame::cosTheta(bRec.wi) <= 0.0f)
            return Color3f(0.0f);
        auto sample = _sample;
        //check if specular or diffuse
        if (_sample.y() < m_ks)
        {
            //specular

            //scaling and offsetting the sample; range of _sample.y() is [0, m_ks] and want it to [0,1]
            float scaledsample = _sample.y() / m_ks;
            Point2f scaledpoint(_sample.x(), scaledsample);
            //sample a normal using squareToBeckmann as specular
            Vector3f whpoint = Warp::squareToBeckmann(scaledpoint, m_alpha);
            Vector3f wh = whpoint.normalized();

            //now we reflect wi using n in order to generate the outgoing direction wo
            Vector3f wo = -bRec.wi + (2 * (wh.dot(bRec.wi)) * wh);
            bRec.wo = wo;
            if (Frame::cosTheta(wo) <= 0)
                return Color3f(0.f);
            sampledbsdf = (eval(bRec) * Frame::cosTheta(bRec.wo)) / pdf(bRec) ;
        }
        else
        {
            //diffuse

            //scaling and offsetting the sample; range of _sample.y() is [0, m_ks] and want it to [0,1]
            float scaledsample = (_sample.y() - m_ks) / (1 - m_ks);
            Point2f scaledpoint(_sample.x(), scaledsample);
            //sample a normal using squareToCosineHemisphere as specular
            Vector3f wopoint = Warp::squareToCosineHemisphere(scaledpoint);
            bRec.wo = wopoint;
            if (Frame::cosTheta(wopoint) <= 0)
                return Color3f(0.0f);
            sampledbsdf = (eval(bRec) * Frame::cosTheta(bRec.wo))/ pdf(bRec);
        }
        return sampledbsdf;
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
