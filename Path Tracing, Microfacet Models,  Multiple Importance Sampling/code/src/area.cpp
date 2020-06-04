#include <nori/emitter.h>
#include <nori/dpdf.h>
#include <nori/scene.h>
#include <nori/frame.h>
#include <nori/mesh.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN


class AreaLight : public Emitter
{
public:
	AreaLight(const PropertyList &propList) {
        radiance = propList.getColor("radiance");
    }

    Color3f eval(const EmitterQueryRecord & lRec) const {
		return radiance;
	}

	//samples the point on the triangle
    Point3f sample(EmitterQueryRecord & lRec, const Point2f & sample, const Mesh &inst) const {

        //lRec.origin is the point that the ray from eye intersects on the scene; point x on the scene where we want to find the radiance
		EmitterQueryRecord SQR(lRec.origin);
		const MatrixXf &m_V  = inst.getVertexPositions();
        const MatrixXf &m_N  = inst.getVertexNormals();
        const MatrixXf &m_UV = inst.getVertexTexCoords();
        const MatrixXu &m_F  = inst.getIndices();

	    //take random sample of the index of the triangle
	    DiscretePDF m_pdf1 = inst.getTrianglePDF();

	    size_t index = m_pdf.sample(sample.x());

		//warp samples to get vector from triangle using warping function
		Vector3f sampledvector = Warp::squareToUniformTriangle(sample);

	    //emitted ray position is found by interpolating the warped directions onto the triangle using the index
	    Point3f interpolatedvertex = sampledvector.x() * m_V.col(m_F(0, index)) + sampledvector.y() * m_V.col(m_F(1, index)) + sampledvector.z() * m_V.col(m_F(2, index));
	    //Normal3f interpolatednormal = sampledvector.x() * m_N.col(m_F(0, index)) + sampledvector.y() * m_N.col(m_F(1, index)) + sampledvector.z() * m_N.col(m_F(2, index)).normalized();

	    lRec.p = interpolatedvertex;

	    //computing the normals
	    //so if the mesh of triangle gives normal of the point at position p then just interpolate the normal, otherwise find the face normal instead
	    if (m_N.size() > 0) {
	        lRec.n = sampledvector.x() * m_N.col(m_F(0, index)) + sampledvector.y() * m_N.col(m_F(1, index)) + sampledvector.z() * m_N.col(m_F(2, index)).normalized();
	    }
	    else {
	        Point3f p0 = m_V.col(m_F(0, index));
	        Point3f p1 = m_V.col(m_F(1, index));
	        Point3f p2 = m_V.col(m_F(2, index));
	        Normal3f n = (p1-p0).cross(p2-p0).normalized();
	        lRec.n = n;
	    }
    	return lRec.p;
    }

     float pdf(const EmitterQueryRecord & lRec, const Mesh &inst) const {
 		//pdf is just 1/total area of the emitter
		return 1.0f/inst.getArea();
	}


    std::string toString() const {
    return "Emitter[]";
    }


private:
    Color3f radiance;
};


NORI_REGISTER_CLASS(AreaLight, "area")
NORI_NAMESPACE_END
