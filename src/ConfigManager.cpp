#include "ConfigManager.hpp"

std::ostream &operator<<(std::ostream &os, const CConfigManager &config)
{
    for( const auto& item : config )
    {
        os << item.first;
        os << " = ";
        os << item.second;
        os << std::endl;
    }
    return os;
}

void CConfigManager::parse(int argc, char **argv, char delimiter)
{
    for( int i = 1 ; i < argc ; ++i )
        parse( std::string(argv[i]) , delimiter );
}

void CConfigManager::parse(std::ifstream &file, char delimiter)
{
    std::stringstream ss;
    std::copy( std::istreambuf_iterator<char>(file)
               , std::istreambuf_iterator<char>()
               , std::ostreambuf_iterator<char>(ss) );
    parse( ss );
}

void CConfigManager::parse(const std::string &str, char delimiter)
{
    // find first key
    std::stringstream ss(str);
    parse( ss , delimiter );
}

void CConfigManager::parse(std::stringstream &ss, char delimiter)
{
    if( ss.good() == 0 )
        return;

    // find first key
    std::string key;
    std::getline(ss, key , delimiter );


    std::remove_if( key.begin() , key.end() , [](char c){return std::isspace(c);} );

    if( !key.empty() )
    {
        // get the value
        std::string value;
        std::getline( ss, value , '\n' );

        // insert the value into the map
        std::remove_if( value.begin() , value.end() , [](char c){return std::isspace(c);} );
        (*this)[key] = value;
    }

    parse( ss );
}

std::string CConfigManager::getValue(const std::string &key) const
{
    return this->at(key);
}

bool CConfigManager::hasKey(const std::string &key) const
{
    return ( find(key) != end() );
}
