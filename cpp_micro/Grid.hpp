#ifndef CGRID_HPP
#define CGRID_HPP

#include "eigen3/Eigen/Dense"
#include <vector>

/*
template< typename TValue >
struct CRectGrid : public Eigen::MatrixXX<TValue>
{
public:
    CRectGrid( int dim ):
        m_dim( dim )
    {}


    void setValue( Eigen::Vector3d pos , TValue val )
    {
        this->at( indexOf(pos) ) = val;
    }

    void setValue( size_t i , size_t j , size_t k , TValue val )
    {
        this->at( indexOf(i,j,k) ) = val;
    }

    TValue& getValue( Eigen::Vector3d pos )
    {
        this->at( indexOf(pos) );
    }

    const TValue& getValue( size_t i , size_t j , size_t k )
    {
        return this->at( indexOf(i,j,k) );
    }

    size_t indexOf( size_t i , size_t j , size_t k )
    {
        return i + m_size[0]*( j + m_size[1]*k );
    }

    size_t indexOf( Eigen::Vector3d pos )
    {
        pos = pos - m_min;
        pos.array() = pos.array() / (m_max-m_min).array();
        size_t i,j,k;
        i = pos[0] * m_size[0];
        j = pos[1] * m_size[1];
        k = pos[2] * m_size[2];

        return indexOf( i,j,k );
    }

    Eigen::Vector3d positionOf( size_t i , size_t j , size_t k )
    {
        Eigen::Vector3d v;
        v << i, j, k;
        return m_min + (v.array()*m_stepSize.array()).matrix();
    }

    Eigen::Vector3i multiIndexOf( typename std::vector< TValue >::iterator it )
    {
        size_t n = (size_t)(it - this->begin());
        size_t i = n % m_size[0];

        n = (n-i)/m_size[0];
        size_t j = n % m_size[1];

        n = (n-j) / m_size[1];
        size_t k = n % m_size[2];

        return Eigen::Vector3i() << i,j,k;
    }

    Eigen::Vector3d positionOf( typename std::vector< TValue >::iterator it )
    {
        return m_min + (m_stepSize.array() * multiIndexOf(it).array()).matrix();
    }


private:
    const int m_dim;

    Eigen::Vector3d m_min;
    Eigen::Vector3d m_max;

    Eigen::Vector3d             m_stepSize;
    Eigen::Matrix<size_t,3,1>   m_size;


};
*/

/*
#include <unordered_map>
#include <vector>
#include <array>

#include "boost/functional/hash.hpp"
#include "eigen3/Eigen/Dense"

class CGrid;

enum EPrimitive
{
    EPoint,
    ELine,
    ETriangle,
    EQuad,
    ETetraeda,
    EPolygon
};

class CPrimitive
{

public:

    CPrimitive( CGrid* grid );
    CGrid* getGridPointer();


private:
    CGrid*  m_pGrid;


};



using TVec3 = Eigen::Vector3d;

namespace std {
  template <>
  struct hash<TVec3>
  {
    std::size_t operator()(const TVec3& k) const
    {
        size_t seed = 0;
        boost::hash_combine(seed , k[0]);
        boost::hash_combine(seed , k[1]);
        boost::hash_combine(seed , k[2]);
        return seed;
    }
  };
}

template <typename TValue>
class CGridIterator
{
public:

    std::vector< TValue >::iterator operator()

private:
};

class CGrid
{
public:

    using TIndex = size_t;

public:
    CGrid( int dim );

    std::vector< TVec3 >& vertices()
    {
        return m_vertices;
    }

    std::vector<TIndex> getNearestPoints( int count )
    {
        auto vec = std::vector<TIndex>(count);


    }

    void createMeshgrid( TVec3 start , TVec3 end , int nx , int ny , int nz )
    {
        TVec3 diag = end - start;
        diag[0] /= nx;     diag[1] /= ny;     diag[2] /= nz;

        TVec3 dx = TVec3(diag[0]/nx , 0          , 0 );
        TVec3 dy = TVec3( 0         , diag[0]/nx , 0 );
        TVec3 dz = TVec3( 0         , 0          , diag[0]/nx );

        m_vertices.resize( (nx+1) * (ny+1) * (nz+1) );
        m_primitives.resize( nx*ny*nz );

        size_t  i_vertices = 0;
        auto it_primitves = m_vertices.begin();
        TVec3 pos;

        for( int x = 0; x < nx+1 ; ++x )
            for( int y = 0; y < ny+1 ; ++y )
                for( int z = 0; z < ny+1 ; ++z )
                {
                    pos = start + x*dx + y*dy + z*dz;

                    // generate vertices
                    m_vertices[ i_vertices ] = pos;
                    ++i_vertices;

                    // add a quad primitive
                    // *it_primitves = std::vector<TIndex>({ i_vertices , i_vertices - 1 , i_vertices - (nx+1) , i_vertices - (nx+1) - 1 ,
                    //                                      i_vertices , i_vertices - 1 , i_vertices - (nx+1) , i_vertices - (nx+1) - 1 });
                }

        // generate primitives

    }


private:

    const int dim;


    std::vector< TIndex >   m_primitives;
    std::vector< TVec3 >        m_vertices;



    std::unordered_map< TVec3 , std::vector< CPrimitive* > > m_point_to_neighbours;

};
*/

#endif // CGRID_HPP
