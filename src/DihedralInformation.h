#pragma once
#include <boost/concept_check.hpp>

struct Residue;
const float DEFAULT = -999;

class DihedralInformation
{
public :
    DihedralInformation()
    {
    }

    float phi    = DEFAULT;
    float psi    = DEFAULT;
    float omega  = DEFAULT;
    float tco    = DEFAULT;
    float alpha  = DEFAULT;
    float kappa  = DEFAULT;
    float d      = DEFAULT;
    float theta  = DEFAULT;
    float r      = DEFAULT;
    DihedralInformation(const Residue*res_prev_prev, const Residue*res_prev, const Residue*res, const Residue*res_next, const Residue*res_next_next);

    inline friend std::ostream& operator<<(std::ostream &output, const DihedralInformation &dhi)
    {
        output << "Dihedral " << dhi.phi << ", " <<  dhi.psi;
        return output;
    }
};