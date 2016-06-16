#ifndef CFORCEFIELD_HPP
#define CFORCEFIELD_HPP

// std includes
#include <unordered_map>
#include <map>
#include <iostream>

// eigen includes
#include <Eigen/Dense>

#include "EigenMatrix_hash.hpp"


template< typename T >
T findValue( int             index , std::unordered_map<int,T>&             mp );

template< typename T >
T findValue( Eigen::Vector2i index , std::unordered_map<Eigen::Vector2i,T>& mp );

template< typename T >
T findValue( Eigen::Vector3i index , std::unordered_map<Eigen::Vector3i,T>& mp );

template< typename T >
T findValue( Eigen::Vector4i index , std::unordered_map<Eigen::Vector4i,T>& mp );

class CForceField
{

    using TVec3  = Eigen::Vector3d;

public:

    void setAtomicNumber( const std::string& symbol , const std::string& name , int id );
    int symbolToNumber( const std::string& symbol );


    bool loadFromFile( const std::string& filename );

    Eigen::Vector3d force_bond( const Eigen::Vector3d& a1 , const Eigen::Vector3d& a2 , const Eigen::Vector2i& atomicNumbers );
    Eigen::Vector3d force_lenard_jones( const Eigen::Vector3d& a1 , const Eigen::Vector3d& a2 , const Eigen::Vector2i& atomicNumbers );



#define getMacro( Name , MapName , IndexType , ValueType ) \
    ValueType get##Name( const IndexType& index ) \
    {   return findValue( index , MapName );  }


    getMacro( Mass , m_mass , int , double )
    getMacro( Epsilon , m_epsilon , int , double )
    getMacro( RMinHalf , m_r_min_half , int , double )
    getMacro( Kb , m_kb , Eigen::Vector2i , double )
    getMacro( B0 , m_b0 ,  Eigen::Vector2i  , double )
    getMacro( KTheta , m_k_theta , Eigen::Vector3i  , double )
    getMacro( Theta0 , m_theta0 , Eigen::Vector3i , double )
    getMacro( KChi , m_k_chi , Eigen::Vector4i , double )
    getMacro( N , m_n , Eigen::Vector4i , int )
    getMacro( Delta , m_delta , Eigen::Vector4i , double )
    getMacro( KPsi , m_k_psi , Eigen::Vector4i , double )
    getMacro( Psi0 , m_psi0 , Eigen::Vector4i , double )

private:

    inline int atomID( const std::string& ref );



    std::map<std::string,int>   m_atomicNumbers;

    // parameters...

    std::unordered_map<int,double>  m_mass;
    std::map<int,std::string>       m_names;
    std::map<int,std::string>       m_symbols;

    std::unordered_map<int,double> m_epsilon;
    std::unordered_map<int,double> m_r_min_half;

    std::unordered_map<Eigen::Vector2i,double> m_kb;
    std::unordered_map<Eigen::Vector2i,double> m_b0;

    std::unordered_map<Eigen::Vector3i,double> m_k_theta;
    std::unordered_map<Eigen::Vector3i,double> m_theta0;

    std::unordered_map<Eigen::Vector4i,double> m_k_chi;
    std::unordered_map<Eigen::Vector4i,int>    m_n;
    std::unordered_map<Eigen::Vector4i,double> m_delta;

    std::unordered_map<Eigen::Vector4i,double> m_k_psi;
    std::unordered_map<Eigen::Vector4i,double> m_psi0;

};


template< typename T >
T findValue( int index , std::unordered_map<int,T>& mp )
{
    // try to find it directly
    auto it = mp.find( index );
    if( it != mp.end() )
        return it->second;

    std::cout << "Missing parameter for " << index << "\n" << std::flush;

    mp[index] = 0;

    return 0;
}

template< typename T >
T findValue( Eigen::Vector2i index , std::unordered_map<Eigen::Matrix<int,2,1>,T>& mp )
{
    // try to find it directly
    auto it = mp.find( index );
    if( it != mp.end() )
        return it->second;

    // swap a0 and a1
    std::swap( index[0] , index[1] );

    it = mp.find( index );
    if( it != mp.end() )
    {
        T val = it->second;
        // we don't want to do this step every frame, hence insert it into the map
        mp[index] = val;
        return val;
    }

    // report the missing parameter:

    std::cout << "Missing parameter for " << index.transpose() << "\n" << std::flush;

    mp[index] = 0;
    return 0;
}

template< typename T >
T findValue( Eigen::Vector3i index , std::unordered_map<Eigen::Matrix<int,2,1>,T>& mp )
{
    // try to find it directly
    auto it = mp.find( index );
    if( it != mp.end() )
        return it->second;

    // swap a0 and a2
    std::swap( index[0] , index[2] );

    it = mp.find( index );
    if( it != mp.end() )
    {
        T val = it->second;
        // we don't want to do this step every frame, hence insert it into the map
        mp[index] = val;
        return val;
    }
    // report the missing parameter:

    std::cout << "Missing parameter for " << index.transpose() << "\n" << std::flush;

    mp[index] = 0;
    // T is expected to be numeric
    return 0;
}


template< typename T >
T findValue( Eigen::Vector4i index , std::unordered_map<Eigen::Matrix<int,2,1>,T>& mp )
{
    // try to find it directly
    auto it = mp.find( index );
    if( it != mp.end() )
        return *it;

    // swap a0 and a3
    std::swap( index[0] , index[3] );

    it = mp.find( index );
    if( it != mp.end() )
    {
        T val = *it;
        // we don't want to do this step every frame, hence insert it into the map
        mp[index] = val;
        return val;
    }

    // swap a1 and a2
    std::swap( index[1] , index[2] );

    it = mp.find( index );
    if( it != mp.end() )
    {
        T val = *it;
        // we don't want to do this step every frame, hence insert it into the map
        mp[index] = val;
        return val;
    }

    // report the missing parameter:

    std::cout << "Missing parameter for " << index.transpose() << "\n" << std::flush;

    mp[index] = 0;
    return 0;
}


#endif // CFORCEFIELD_HPP
