#include "ForceField.hpp"

// std includes
#include <fstream>
#include <string>
#include <cmath>

using std::pow;


void CForceField::setAtomicNumber(const std::string &symbol, const std::string &name, int id)
{
    m_atomicNumbers[symbol] = id;
    m_symbols[id] = symbol;
    m_names[id]   = name;
}

int CForceField::symbolToNumber(const std::string &symbol)
{
    return m_atomicNumbers[symbol];
}

bool CForceField::loadFromFile(const std::string &filename)
{
    setAtomicNumber("X","placeholder", -1);

    std::ifstream file( filename );

    enum EParameterType{
        EUnknown,
        EMass,
        EBonds,
        EAngles,
        EDihedrals,
        EImproper,
        ENonBonded
    } state;

    state = EUnknown;

    std::string line;
    std::string value;

    while( file.good() )
    {
        std::getline( file , line );

        auto it = line.find('!');
        if( it != std::string::npos )
            line = line.substr( 0 , it );

        it = line.find('*');
        if( it != std::string::npos )
            line = line.substr( 0 , it );

        std::stringstream ss;
        ss << line;
        ss >> value;


        bool newCategory = true;

        // there are different cases

        // is the current value a new category?

        if( value == "BONDS" )
            state = EBonds;
        else if( value == "ANGLES" )
            state = EAngles;
        else if( value == "DIHEDRALS" )
            state = EDihedrals;
        else if( value == "IMPROPER" )
            state = EImproper;
        else if( value == "NONBONDED")
            state = ENonBonded;
        else if( value == "MASS")
            state = EMass;
        else if( value == "CMAP")
            state = EUnknown;
        else
            newCategory = false;

        // ok, if the current line is not the title of a new category,
        // then we exect it to contain parameters or a comment

        if( !newCategory )
        {
            switch(state)
            {
            case EBonds:
            {
                std::string a1,a2;
                double kb,b0;

                a1 = value;
                ss >> a2 >> kb >> b0;

                Eigen::Vector2i v;
                v[0] = atomID(a1);
                v[1] = atomID(a2);

                m_kb[v] = kb;
                m_b0[v] = b0;
            }
                break;


            case EAngles:
            {
                std::string a1,a2,a3;
                double k_theta,theta0;

                a1 = value;
                ss >> a2 >> a3 >> k_theta >> theta0;

                Eigen::Vector3i v;
                v[0] = atomID(a1);
                v[1] = atomID(a2);
                v[2] = atomID(a3);

                m_k_theta[v] = k_theta;
                m_theta0[v] = theta0;
            }
                break;


            case EDihedrals:
            {
                std::string a1,a2,a3,a4;
                double k_chi,delta;
                int n;

                a1 = value;
                ss >> a2 >> a3 >> a4 >> k_chi >> n >> delta;

                Eigen::Vector4i v;
                v[0] = atomID(a1);
                v[1] = atomID(a2);
                v[2] = atomID(a3);
                v[3] = atomID(a4);

                m_k_chi[v] = k_chi;
                m_n[v] = n;
                m_delta[v] = delta;
            }
                break;


            case EImproper:
            {
                std::string a1,a2,a3,a4;
                double k_psi,psi0;

                a1 = value;
                ss >> a2 >> a3 >> a4 >> k_psi >> psi0;

                Eigen::Vector4i v;
                v[0] = atomID(a1);
                v[1] = atomID(a2);
                v[2] = atomID(a3);
                v[3] = atomID(a4);

                m_k_psi[v] = k_psi;
                m_psi0[v] = psi0;
            }
                break;

            case ENonBonded:
            {
                std::string a1;
                double epsilon,r_min_half;

                a1 = value;
                ss >> epsilon >> r_min_half;

                int v = atomID(a1);

                m_epsilon[v] = epsilon;
                m_r_min_half[v] = r_min_half;
            }
                break;

            case EMass:
            {
                std::string a1;
                double mass;
                int id;

                ss >> id >> a1 >> mass;

                int v = atomID(a1);

                m_mass[v] = mass;
            }
                break;


            default: // do nothing
                break;
            }
        }
    }


    return true;
}

Eigen::Vector3d CForceField::force_bond(const Eigen::Vector3d &a1, const Eigen::Vector3d &a2, const Eigen::Vector2i &atomicNumbers)
{
    // V_bond = kb(r-b0)**2
    // grad V = 2*kb*(r-b0) * vec r / norm r;

    TVec3 r = a2 - a1;

    return 2. * getKb(atomicNumbers) * ( r.norm() - getB0(atomicNumbers) ) * r.normalized();
}

Eigen::Vector3d CForceField::force_lenard_jones(const Eigen::Vector3d &a1, const Eigen::Vector3d &a2, const Eigen::Vector2i &atomicNumbers)
{
    //V = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
    //
    //epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
    //Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
    //
    // grad V = epsij * ...

    double eps  = sqrt( getEpsilon(atomicNumbers[0]) + getEpsilon(atomicNumbers[1]) );
    double Rmin = getRMinHalf(atomicNumbers[0]) + getRMinHalf(atomicNumbers[1]);

    TVec3 r = a2 - a1;

    double d = Rmin / r.squaredNorm();

    return eps * ( -12. * pow(d,12) + 12 * pow(d,6) ) * r / r.squaredNorm();
}

int CForceField::atomID(const std::string &ref)
{
    return m_atomicNumbers[ref];
}
