#ifndef CONFIGMANAGER
#define CONFIGMANAGER

#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>

class CConfigManager : public std::map<std::string,std::string>
{
public:

    void parse( int argc , char** argv , char delimiter = '=' );
    void parse( std::ifstream& file , char delimiter = '=' );
    void parse( const std::string& str , char delimiter = '=' );
    void parse( std::stringstream& ss , char delimiter = '=' );

    template< typename T>
    bool hasKey( const std::string& key ) const;

    bool hasKey( const std::string& key ) const;

    template<typename T>
    T   getValue( const std::string& key ) const;

    std::string   getValue( const std::string& key ) const;

};

std::ostream& operator<<( std::ostream& os , const CConfigManager& config );

template< typename T>
bool CConfigManager::hasKey( const std::string& key ) const
{
    auto it = this->find(key);

    // check if key exists
    if( it != this->end() )
    {
        // check type
        std::stringstream ss;
        ss << it->second;
        T value;
        if( ss >> value )
            return true;
    }

    return false;
}

template<typename T>
T   CConfigManager::getValue( const std::string& key ) const
{
    T value;
    auto it = this->find(key);

    if( it == this->end() )
        return value;   // undefined behavior...

    std::stringstream ss( it->second );
    ss >> value;
    return value;
}

#endif // CONFIGMANAGER

