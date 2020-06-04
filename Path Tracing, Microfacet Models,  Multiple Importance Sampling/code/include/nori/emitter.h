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

#pragma once

#include <nori/object.h>
#include <nori/dpdf.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Convenience data structure used to pass multiple
 * parameters to the evaluation and sampling routines in \ref emitter
 */
struct EmitterQueryRecord {

    //origin from where we sample the emitter
    Point3f origin;
    
    // Sampled point on emitter
    Point3f p;

    // Normal at emitter point
    Normal3f n;

    //direction between hit point and emitter point
    Vector3f wi;

    // probability of sampled point being an emitter
    float pdf;

    //Shadow Ray
    Ray3f shadowRay;

    //Empty constructor
    EmitterQueryRecord() {}

    /// Create a new record for sampling the emitter
    EmitterQueryRecord(const Point3f &origin): origin(origin) { }

    /**
     * \brief Create a query record that can be used to query the
     * sampling density after having intersected an area emitter
     */
    EmitterQueryRecord(const Point3f &origin, const Point3f &p, const Normal3f &n) :
        origin(origin), p(p), n(n) {
        wi = (p - origin).normalized();
    }

};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:

    virtual Color3f eval(const EmitterQueryRecord &lRec) const = 0;

    virtual Point3f sample(EmitterQueryRecord &lRec, const Point2f &sample, const Mesh &inst) const = 0;

    virtual float pdf(const EmitterQueryRecord &lRec, const Mesh &inst) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }
    

    protected:

    DiscretePDF m_pdf;


};

NORI_NAMESPACE_END
