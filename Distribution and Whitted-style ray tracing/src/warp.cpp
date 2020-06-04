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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    if(sample.x() <= 0.5f && sample.y() <= 0.5f && sample.x() >= 0.0f && sample.y() >= 0.0f)
        return Point2f((sqrt(2 * sample.x()) - 1.0f), (sqrt(2 * sample.y()) - 1.0f));
    else if(sample.x() <= 0.5f && sample.y() > 0.5f && sample.x() >= 0.0f && sample.y() <= 1.0f)
        return Point2f((sqrt(2 * sample.x())) - 1.0f, (1.0f + sqrt(2 * (1.0f - sample.y()))));
    else if(sample.x() > 0.5f && sample.y() <= 0.5f && sample.x() < 1.0f && sample.y() >= 0.0f)
        return Point2f((1.0f + sqrt(2 * (1.0f - sample.x()))), (sqrt(2 * sample.y())) - 1.0f);
     else if(sample.x() > 0.5f && sample.y() > 0.5f && sample.x() < 1.0f && sample.y() < 1.0f)
        return Point2f((1.0f + sqrt(2 * (1.0f - sample.x()))), (1.0f + sqrt(2 * (1.0f - sample.y()))));
}

float Warp::squareToTentPdf(const Point2f &p) {
    if (p.x()>=-1 && p.y()>=-1 && p.x()<=1 && p.y()<=1)
        return (1-abs(p.x()))*(1-abs(p.y()));
    else 
        return 0;
}

Vector3f Warp::squareToUniformTriangle(const Point2f &sample) {
    float u = 1.0f - sqrt(sample.x());
    float v = sample.y() * sqrt(sample.x());
    return Vector3f(u, v, 1.0f - u - v);
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float two_pi_y = 2 * M_PI * sample.y();
    float sqrt_x = sqrt(sample.x());
    return Point2f(cos(two_pi_y) * sqrt_x, sin(two_pi_y) * sqrt_x);
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return (p.x() * p.x() + p.y() * p.y() < 1.0f) ? 1.0f / M_PI : 0.0f;
    //throw NoriException("Warp::squareToUniformDiskPdf() is not yet implemented!");
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float z = 1.0f - 2.0f * sample.x();
    //float r = sqrt(std::max(0.0f, 1.0f - z * z));
    float r = sqrt(1.0f - z * z);
    float two_pi_y = 2.0f * M_PI * sample.y();
    float sin_phi = sin(two_pi_y);
    float cos_phi = cos(two_pi_y);
    return Vector3f(r * cos_phi, r * sin_phi, z);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return (abs(v.x() * v.x() + v.y() * v.y() + v.z() * v.z() - 1.0f) < 1e-6) ? INV_FOURPI : 0.0f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float phi = 2 * M_PI * sample.x();
    float costheta = sample.y();
    float x = cos(phi) * sqrt(1 - (costheta * costheta));
    float y = sin(phi) * sqrt(1 - (costheta * costheta));
    return Vector3f(x, y, costheta);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    float radius = v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
    float margin = 0.000001;
    if((radius <= 1 + margin) && (radius >= 1 - margin) && v.z() >= 0)
        return 1 / (2 * M_PI);
    else
        return 0;
    //return (abs(v.x() * v.x() + v.y() * v.y() + v.z() * v.z() - 1.0f) < 1e-6 && v.z() > 0.0f) ? 1.0f / 2 * M_PI : 0.0f;

}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    float x = std::sqrt(sample.x()) * cos(sample.y() * 2 * M_PI);
    float y = std::sqrt(sample.x()) * sin(sample.y() * 2 * M_PI);
    float z = std::sqrt(1.0f - sample.x());
    return Vector3f(x, y, z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    float radius = v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
    float margin = 0.000001;
    if((radius <= 1 + margin) && (radius >= 1 - margin) && v.z() >= 0)
        return v.z() / (M_PI);
    else
        return 0;
    //return (abs(v.x() * v.x() + v.y() * v.y() + v.z() * v.z() - 1.0f) < 1e-6) ? v.z() / M_PI: 0.0f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
}

NORI_NAMESPACE_END
