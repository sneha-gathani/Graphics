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

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <tbb/tbb.h>


NORI_NAMESPACE_BEGIN
int num_leafnodes = 0;
int num_interiornodes = 0;
int num_triangles = 0;
void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build()
{
    std::cout << "Building mesh ..." << endl;
    std::cout << "Memory Requirement: " << sizeof(Nodetype) << "bytes" <<endl;
    //List of all the indices of the triangles in the mesh
    std::list<int> indextriangle;
    for(uint32_t idx = 0; idx < m_mesh->getTriangleCount(); idx++)
        indextriangle.push_back(idx);
    std::cout << "Built mesh ..." << endl;
    //Gets the start time before building the octree
    std::time_t start = std::time(nullptr);
    //depth of the octree is set to 6
    root = build_octree(m_bbox, indextriangle, 6);
    //Gets the end time after building the octree
    std::time_t end = std::time(nullptr);
    std::cout << "Time for making octree: " << end - start << " seconds." << endl;

    //std::cout << "Root" << m_bbox.toString() << endl;
    std::cout << "Octree made successfully ..." << endl;
    std::cout << "Number of leaf nodes in the octree are: " << num_leafnodes << endl;
    std::cout << "Number of interior nodes in the octree are: " << num_interiornodes - num_leafnodes - 1 << endl;
    std::cout << "Average number of triangles in leaf nodes of the octree are: " << num_triangles/num_leafnodes << endl;
}

Node *Accel::build_octree(BoundingBox3f boundingbox, std::list<int> triindices, int depth) const {
    //make a new node
    Node* node = new Node();
    num_interiornodes++;
    //setting the bounding box of the node
    node->bbox = boundingbox;
    //if no triangles are found in the node
    if(triindices.empty())
        return node;
    //threshold of the number of ttriangles in the leaf node is set to 10
    //if number of triangles are less than threshold; here 10 and also if the depth of the octree is reached; make the leaf node
    if(triindices.size()<10 || depth == 1)
    {
        num_leafnodes++;
        node->triangles = triindices;
        return node;
    }
    //array of 8 elements which will be the children of node; have the list of triangles belonging the bounding box enclosed in the node
    std::list<int> trichildindices[8];
    //parallel computing, but did not work, as OpenMP did not sync with nori
    //#pragma omp parallel num_threads(8)
    //{
        //#pragma omp parallel for
            //tried to use the existing tbb block
            //tbb::parallel_for(tbb::blocked_range<size_t>(0,8), 
            for(int i = 0; i<8; i++)
            {
                //Get the minimum point of the new bounding box which is one of the 8 boxes
                Vector3f minpoint(std::min(boundingbox.getCenter()[0], boundingbox.getCorner(i)[0]),
                    std::min(boundingbox.getCenter()[1], boundingbox.getCorner(i)[1]), 
                    std::min(boundingbox.getCenter()[2], boundingbox.getCorner(i)[2]));
                //Get the maximum point of the new bounding box which is one of the 8 boxes
                Vector3f maxpoint(std::max(boundingbox.getCenter()[0], boundingbox.getCorner(i)[0]),
                    std::max(boundingbox.getCenter()[1], boundingbox.getCorner(i)[1]), 
                    std::max(boundingbox.getCenter()[2], boundingbox.getCorner(i)[2]));
                //Make the new bounding box with these min and max points
                BoundingBox3f partbox = BoundingBox3f(minpoint, maxpoint);

                //std::cout << partbox.toString() << endl;

                //iterator to traverse over triangle indices in the bounding box defined by the node
                std::list<int>::iterator triangle = triindices.begin();

                //until all triangles belonging to that bounding box are not traversed
                while(triangle != triindices.end())
                {
                    //get bounding box of the independent triangle
                    BoundingBox3f trianglebox = m_mesh->getBoundingBox(*triangle);
                    //check if triangle box overlaps with the partitioned bounding box
                    bool overlap = trianglebox.overlaps(partbox);
                    //if yes, send the triangle to the respective child node
                    if(overlap)
                    {
                        trichildindices[i].push_back(*triangle);
                        num_triangles++;
                    }
                    triangle++;
                }
                node->triangles = trichildindices[i];
                //recursively call build_octree for each of the child nodes
                node->child[i] = build_octree(partbox, trichildindices[i], depth-1);
                //);
            }
    //}
    return node;
}

//Part 3 of where the sorted array of triangles is passed
bool Accel::traversefaster(const Ray3f &ray_, Intersection &its, bool shadowRay, Node *node, bool &foundIntersection, uint32_t &f) const{
    Ray3f ray(ray_);
    //std::cout << "NodeBox" << node->bbox.toString() << endl;
    if(node->bbox.rayIntersect(ray))
    {
        if(node->child[0] == nullptr)
        {
            if(!node->triangles.empty())
            {
                num_triangles += node->triangles.size();
                std::list<int>:: iterator triangle = node->triangles.begin();
                while(triangle != node->triangles.end())
                {
                    float u, v, t;
                    if (m_mesh->rayIntersect(*triangle, ray, u, v, t)) 
                    {
                        //std::cout << "triangle" << *triangle << endl;
                    //An intersection was found! Can terminate immediately if this is a shadow ray query 
                    if (shadowRay)
                        return true;
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = m_mesh;
                    f = *triangle;
                    foundIntersection = true;
                    }
                    triangle++;
                }
            }
        }
        else
        {
            //Part 3 - fast intersection

            //get closest intersection of the ray and bounding box
            float nearestIntersection;
            //get farthest intersection of the ray and bounding box
            float farthestIntersection;
            std::vector<std::pair<float, int> > indexmap;
            //make a vector of pair of values having point of intersection between ray and bounding box and child it belongs to
            #pragma omp parallel num_threads(8)
            {
                #pragma omp parallel for
                for(int i = 0; i<8; i++)
                {
                    node->child[i]->bbox.rayIntersect(ray, nearestIntersection, farthestIntersection);
                    indexmap.push_back(std::make_pair(nearestIntersection, i));
                }
                //sort this vector from closest to farthest point
                sort(indexmap.begin(), indexmap.end());

                //check intersection for each of the 8 child nodes
                #pragma omp parallel for
                for(int i = 0; i<8; i++)
                {
                    if(!foundIntersection)
                    {
                        foundIntersection = traversefaster(ray_, its, shadowRay, node->child[indexmap[i].second], foundIntersection, f);
                    }
                }
            }
        }
    }
    return foundIntersection;
}

//Part 2 - Plain traverse through the octree
bool Accel::traverse(const Ray3f &ray_, Intersection &its, bool shadowRay, Node *node, bool &foundIntersection, uint32_t &f) const{
    Ray3f ray(ray_);
    if(node->bbox.rayIntersect(ray) )
    {
        if(node->child[0]==nullptr)
        {
            if(!node->triangles.empty())
            {
                //std::cout << "Checking for leafnode" << endl;
                num_triangles += node->triangles.size();
                //std::cout << "Num triangles" << num_triangles << endl;
                std::list<int>:: iterator triangle = node->triangles.begin();
                while(triangle != node->triangles.end())
                {
                    float u, v, t;
                    if (m_mesh->rayIntersect(*triangle, ray, u, v, t)) 
                    {
                    //An intersection was found! Can terminate immediately if this is a shadow ray query 
                    if (shadowRay)
                        return true;
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = m_mesh;
                    f = *triangle;
                    foundIntersection = true;
                    }
                    triangle++;
                }
            }
        }
        else
        {
            //get closest intersection of the ray and bounding box
            float nearestIntersection;
            //get farthest intersection of the ray and bounding box
            float farthestIntersection;
            std::vector<std::pair<float, int> > indexmap;
            //make a vector of pair of values having point of intersection between ray and bounding box and child it belongs to
            #pragma omp parallel num_threads(8)
            {
                #pragma omp parallel for
                for(int i = 0; i<8; i++)
                {
                    node->child[i]->bbox.rayIntersect(ray, nearestIntersection, farthestIntersection);
                    indexmap.push_back(std::make_pair(nearestIntersection, i));
                }
            }
            for(int i = 0; i<8; i++)
            {
                foundIntersection = traverse(ray_, its, shadowRay, node->child[i], foundIntersection, f);
            }
        }
    }
    return foundIntersection;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection
    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)


// /* Brute force search through all triangles */
//     for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
//         float u, v, t;
//         if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
// //          An intersection was found! Can terminate immediately if this is a shadow ray query 
//         if (shadowRay)
//              return true;
//         ray.maxt = its.t = t;
//         its.uv = Point2f(u, v);
//         its.mesh = m_mesh;
//         f = idx;
//         foundIntersection = true;
//     }
// }
        //Part 2 - Ray Traversal
        //std::cout << "Rendering using only traveral through octree ..." << endl;
        foundIntersection = traverse(ray_, its, shadowRay, root, foundIntersection, f);
        

        //Part 3 - Faster Ray Traversal
        // std::cout << "Rendering using faster traveral through octree ..." << endl;
        //foundIntersection = traversefaster(ray_, its, shadowRay, root, foundIntersection, f);
   


    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }
    return foundIntersection;
}

NORI_NAMESPACE_END
