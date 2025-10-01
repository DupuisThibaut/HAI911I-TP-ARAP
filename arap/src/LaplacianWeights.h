#ifndef LAPLACIANWEIGHTS_H
#define LAPLACIANWEIGHTS_H

#include <vector>
#include <map>
#include "Mesh.h"



//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//
// This class describes the implementation of the well known cotangent weights
// Here, weight( v1 , v2 ) is half the sum of the cotangent angles of the opposite corners in the triangles (v1,v2,other)
//   Careful ! these weights can be negative !
//
// Additionally, we provide a scheme for non-symmetric (!!!) POSITIVE edge weights
//   Careful ! these weights are non-symetric ! e_ij != e_ji
//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//





// access to an edge is of complexity O( log(val) ), with val the average valence of the vertices

//---------------------------------   YOU DO NOT NEED TO CHANGE THE FOLLOWING CODE  --------------------------------//
class LaplacianWeights{
private:

public:
    unsigned int n_vertices;
    std::vector< std::map< unsigned int , double > > edge_weights;
    std::vector< double > vertex_weights;
    LaplacianWeights() : n_vertices(0) {}
    void clear()
    {
        n_vertices = 0;
        edge_weights.clear();
        vertex_weights.clear();
    }
    ~LaplacianWeights() { clear(); }

    double sumVertexWeights() const {
        double s = 0.0;
        for(unsigned int v = 0 ; v < n_vertices ; ++v) {
            s += vertex_weights[v];
        }
        return s;
    }

    void resize( unsigned int nVertices )
    {
        clear();
        if( nVertices > 0 )
        {
            n_vertices = nVertices;
            edge_weights.resize(nVertices);
            vertex_weights.resize(nVertices , 0.0);

            for( unsigned int v = 0 ; v < nVertices ; ++v )
            {
                edge_weights[v].clear();
                vertex_weights[v] = 0.0;
            }
        }
    }
    unsigned int get_n_adjacent_edges( unsigned int vertex_index ) const
    {
        return edge_weights[vertex_index].size();
    }
    double get_edge_weight( unsigned int v1 , unsigned int v2 ) const
    {
        std::map< unsigned int , double >::const_iterator it = edge_weights[v1].find(v2);
        if( it == edge_weights[v1].end() ) return 0.0;
        return it->second;
    }
    unsigned int get_n_vertices() const
    {
        return n_vertices;
    }


    std::map< unsigned int , double >::iterator get_weight_of_adjacent_edges_it_begin( unsigned int v1 )
    {
        return edge_weights[v1].begin();
    }
    std::map< unsigned int , double >::iterator get_weight_of_adjacent_edges_it_end( unsigned int v1 )
    {
        return edge_weights[v1].end();
    }

    std::map< unsigned int , double >::const_iterator get_weight_of_adjacent_edges_it_begin( unsigned int v1 ) const
    {
        return edge_weights[v1].begin();
    }
    std::map< unsigned int , double >::const_iterator get_weight_of_adjacent_edges_it_end( unsigned int v1 ) const
    {
        return edge_weights[v1].end();
    }

    double get_vertex_weight( unsigned int v ) const
    {
        return vertex_weights[v];
    }

    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //---------------------------------  CODE TO CHANGE  --------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//

    //--------------------------   Cotangent weights:   ------------------------//
    // Weight of edge eij : i<->j is the sum of cotangent of opposite angles divided by 2
    // wij = 1/2 * (cot(alpha_ij) + cot(beta_ij)) alpha_ij and beta_ij being the two opposite angles of the edge ij

    //Fonction test des deux triangles avant nouvelle explication
    // unsigned int testTriangle(unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4){
    //     if(v0==v2){
    //         if(v1==v3){
    //             return v4;
    //         }
    //         if(v1==v4){
    //             return v3;
    //         }
    //     }if(v0==v3){
    //         if(v1==v2){
    //             return v4;
    //         }
    //         if(v1==v4){
    //             return v2;
    //         }
    //     }if(v0==v4){
    //         if(v1==v2){
    //             return v3;
    //         }
    //         if(v1==v3){
    //             return v2;
    //         }
    //     }
    //     return UINT_MAX;
    // }
    
    void buildCotangentWeightsOfTriangleMesh( const Mesh& mesh){

        // For each triangle : 
            // Compute its edge p0, p1, p2
            // Compute opposite angle for each edge
            // Compute cotangent of this angle
            // Add to edge weight 

        resize(mesh.V.size());
        for(size_t i=0;i<mesh.T.size();i++){
            unsigned int v0=mesh.T[i][0];
            unsigned int v1=mesh.T[i][1];
            unsigned int v2=mesh.T[i][2];
            Vec3 s1 = mesh.V[v0].p;
            Vec3 s2 = mesh.V[v1].p;
            Vec3 s3 = mesh.V[v2].p;
            double e1=(s2-s1).squareLength();
            double e2=(s3-s2).squareLength();
            double e3=(s1-s3).squareLength();
            double w1=Vec3::dot(s2-s1,s3-s1);
            double w2=Vec3::dot(s1-s2,s3-s2);
            double w3=Vec3::dot(s1-s3,s2-s3);
            Vec3 m;
            double l1=0.0,l2=0.0,l3=0.0;
            double t=Vec3::cross(s2-s1,s3-s1).norm()/2.0;
            if(w1<0.0) {
                m=(s2+s3)/2.0;
                l2=sqrt(((s1+s3)/2.0-m).squareLength()/e3);
                l1=sqrt(((s1+s2)/2.0-m).squareLength()/e1);
                edge_weights[v0][v2]+=l2;
                edge_weights[v2][v0]+=l2;
                edge_weights[v0][v1]+=l1;
                edge_weights[v1][v0]+=l1;
                vertex_weights[v0]+=t/2.0;
                vertex_weights[v1]+=t/4.0;
                vertex_weights[v2]+=t/4.0;
            }else if(w2<0.0){
                m=(s1+s3)/2.0;
                l3=sqrt(((s3+s2)/2.0-m).squareLength()/e2);
                l1=sqrt(((s1+s2)/2.0-m).squareLength()/e1);
                edge_weights[v1][v2]+=l3;
                edge_weights[v2][v1]+=l3;
                edge_weights[v0][v1]+=l1;
                edge_weights[v1][v0]+=l1;
                vertex_weights[v0]+=t/4.0;
                vertex_weights[v1]+=t/2.0;
                vertex_weights[v2]+=t/4.0;
            }else if(w3<0.0){
                m=(s1+s2)/2.0;
                l3=sqrt(((s3+s2)/2.0-m).squareLength()/e2);
                l2=sqrt(((s1+s3)/2.0-m).squareLength()/e3);
                edge_weights[v1][v2]+=l3;
                edge_weights[v2][v1]+=l3;
                edge_weights[v0][v2]+=l2;
                edge_weights[v2][v0]+=l2;
                vertex_weights[v0]+=t/4.0;
                vertex_weights[v1]+=t/4.0;
                vertex_weights[v2]+=t/2.0;
            }else{
                double cot1=w1/(2.0*sqrt(e1*e3-w1*w1));
                double cot2=w2/(2.0* sqrt(e2*e1-w2*w2));
                double cot3=w3/(2.0*sqrt(e2*e3-w3*w3));
                edge_weights[v1][v2]+=cot1;
                edge_weights[v2][v1]+=cot1;
                edge_weights[v0][v2]+=cot2;
                edge_weights[v2][v0]+=cot2;
                edge_weights[v1][v0]+=cot3;
                edge_weights[v0][v1]+=cot3;
                vertex_weights[v1]+=cot1*e2/2.0;
                vertex_weights[v2]+=cot1*e2/2.0;
                vertex_weights[v0]+=cot2*e3/2.0;
                vertex_weights[v2]+=cot2*e3/2.0;
                vertex_weights[v0]+=cot3*e1/2.0;
                vertex_weights[v1]+=cot3*e1/2.0;
            }
        }
    }

    //---------------------------------   YOU DO NOT NEED TO CHANGE THE FOLLOWING CODE  --------------------------------//

    //------------------------------------   Barycentric weights:  --------------------------------//
    // Weight of vertex i is its barycentric area (sum of areas of connected triangles divided by 3)
    // Weight of directed edge i->j is its barycentric area divided by the weight of vertex i
    // CAREFUL: weight(i->j) is different from weight(j->i) !
    //          EDGE WEIGHTS ARE NOT SYMMETRIC !
    // These weights are all positive

    template< class vertex_t , class triangle_t >
    void buildBarycentricWeightsOfTriangleMesh( const std::vector< vertex_t > & vertices , const std::vector< triangle_t > & triangles )
    {
        resize(vertices.size());

        for( unsigned int t = 0 ; t < triangles.size() ; ++t )
        {
            unsigned int v0 = triangles[t][0];
            unsigned int v1 = triangles[t][1];
            unsigned int v2 = triangles[t][2];

            Vec3 p0( vertices[v0][0] , vertices[v0][1] , vertices[v0][2] );
            Vec3 p1( vertices[v1][0] , vertices[v1][1] , vertices[v1][2] );
            Vec3 p2( vertices[v2][0] , vertices[v2][1] , vertices[v2][2] );

            double t_area_by_3 = Vec3::cross( p1 - p0 , p2 - p0 ).norm() / 6.0;

            // Barycentric weights:
            edge_weights[v1][v2] += t_area_by_3;
            edge_weights[v2][v1] += t_area_by_3;

            edge_weights[v0][v2] += t_area_by_3;
            edge_weights[v2][v0] += t_area_by_3;

            edge_weights[v1][v0] += t_area_by_3;
            edge_weights[v0][v1] += t_area_by_3;

            vertex_weights[v0] += t_area_by_3;
            vertex_weights[v1] += t_area_by_3;
            vertex_weights[v2] += t_area_by_3;
        }
        for( unsigned int v = 0 ; v < n_vertices ; ++v )
        {
            double v_area = vertex_weights[v];
            for( std::map< unsigned int , double >::iterator it = edge_weights[v].begin() ; it != edge_weights[v].end() ; ++it )
            {
                // todo:: make sure it does not mess up the iterator (it should not, they are supposed to be ordered solely based on it->first, but, not sure)
                edge_weights[v][it->first] = it->second / v_area;
            }
        }
    }

};





#endif // LAPLACIANWEIGHTS_H
