#ifndef EIGENMATRIX_HASH_HPP
#define EIGENMATRIX_HASH_HPP


// boost includes
#include <boost/functional/hash.hpp>

// eigen includes
#include <Eigen/Dense>

namespace std {
  template < typename T, int N >
  struct hash<Eigen::Matrix<T,N,1>>
  {
    std::size_t operator()(const Eigen::Matrix<T,N,1>& k) const
    {
        size_t seed = 0;
        for( int i = 0 ; i < N ; ++i )
        boost::hash_combine(seed , k[i]);
        return seed;
    }
  };
}


#endif // EIGENMATRIX_HASH_HPP
