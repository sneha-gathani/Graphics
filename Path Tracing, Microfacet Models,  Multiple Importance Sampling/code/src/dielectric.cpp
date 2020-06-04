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
#include <nori/dpdf.h>
#include <nori/common.h>
#include <nori/mesh.h>
#include <Eigen/Geometry>
#include <atomic>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }
        /**
     * \brief Sample the BSDF and return the importance weight (i.e. the
     * value of the BSDF * cos(theta_o) divided by the probability density
     * of the sample with respect to solid angles).
     *
     * \param bRec    A BSDF query record
     * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
     *
     * \return The BSDF value divided by the probability density of the sample
     *         sample. The returned value also includes the cosine
     *         foreshortening factor associated with the outgoing direction,
     *         when this is appropriate. A zero value means that sampling
     *         failed.
     */
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {

        //calculating the fresnel reflection coefficient
        //cos of incident direction and normal
        //m_extIOR is the refractive index of the side containing the surface normal
        //m_intIOR is the refractive index of the interior
        float cosThetai = Frame::cosTheta(bRec.wi);
        float fresnelterm = fresnel(cosThetai, m_extIOR, m_intIOR);

        //as fresnel coefficient gives amount of reflected and refracted light, it is used to check whether reflection or refraction is happening
        //if random sample is less than the fresnel term, then reflection will happen otherwise refraction
        if(sample.y() < fresnelterm)
        {
            //reflection; in local coordinate system
            bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
            bRec.measure = EDiscrete;

            /* Relative index of refraction: no change */
            bRec.eta = 1.0f;
            return Color3f(1.0f);
        }
        else
        {
            //refraction
            float eta_incidence = m_extIOR;
            float eta_outgoing = m_intIOR;
            Vector3f n = Vector3f(0.f, 0.f, 1.f);
            //checking sign of cos theta - outer surface pointing; within same hemisphere; no signs change
            //inner surface pointing; inside surface pointing; opposite side of the hemisphere; signs change
            //happens because cos theta value is positive in hemisphere and negative opposite of hemisphere, but we want eta_incidence to always be positive
            //need to check whether it is inside or outside refraction
            if(Frame::cosTheta(bRec.wi) < 0)
            {
                eta_incidence = m_intIOR;
                eta_outgoing = m_extIOR;
                n = -n;
            }
            //final eta value now; eta_incidence/eta_outgoing
            bRec.eta = eta_incidence/eta_outgoing;

            // WHY WOULD THIS BRUTE FORCE NOT WORK?
            // //finding costhetat using snell's law
            // //we already know costhetai = costheta(wi)
            // float sinsqthetai = std::max(0.0f, 1.0f - (cosThetai * cosThetai));
            // //snells law
            // float sinsqthetat = bRec.eta * bRec.eta * sinsqthetai;
            // //handle total internal reflection
            // if(sinsqthetat >= 1)
            //     return Color3f(1.0);
            // float cosThetat = std::sqrt(1 - sinsqthetat);
            // bRec.wo = bRec.eta * -bRec.wi + (bRec.eta * cosThetai + cosThetat) * n;

            //calculating direction of w_t using snells law
            auto wt_part1 = -eta_incidence / eta_outgoing * (bRec.wi - bRec.wi.dot(n)*n);
            auto wt_part2 = -std::sqrt(1 - (eta_incidence / eta_outgoing) * (eta_incidence / eta_outgoing) * (1 - std::pow(bRec.wi.dot(n),2))) * n;
            bRec.wo = wt_part1 + wt_part2;
         
            return Color3f(1.0f/(bRec.eta * bRec.eta));

        }
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
    DiscretePDF m_pdf;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
