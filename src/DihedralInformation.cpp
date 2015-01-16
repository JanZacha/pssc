#include <limits>
#include <cmath>

#include <boost/math/special_functions/round.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/tr1/tuple.hpp>

#include "standardlibs.h"
#include "DihedralInformation.h"

#include "Residue.h"

#define PI 3.14159265358979323846

namespace tr1 = std::tr1;
using std::abs; // important. for some reason an integer version is used otherwise.

tr1::tuple<double, double, double> ThetaD(double phi, double psi)
{

    phi = phi*PI/180;
    psi = psi*PI/180;
    const double su = .5 * (phi + psi);
    const double di = .5 * (psi - phi);

    const double R = -0.8235 * sin(su) - 0.0222 * sin(di);
    const double S =  2.999  * cos(su) - 0.657  * cos(di);

    int handedness = 1;
    handedness *= boost::math::sign(R);

    const double theta_half = acos(abs(R));
    handedness *= boost::math::sign(S);

    const double d = abs((S / sin(theta_half)));
    const double theta = 2 * theta_half;

    const double r_sq = (3.8*3.8 - d*d) / (2 - 2 * cos(theta));

    return tr1::tuple<double, double, double>(handedness*theta * 180 / PI, d, sqrt(r_sq));
}

float calc_dihedral_angle(const vector3d& p1, const vector3d& p2, const vector3d& p3, const vector3d& p4)
{
    vector3d v12 = p1 - p2;   // vector from p2 to p1
    vector3d v43 = p4 - p3;   // vector from p3 to p4

    vector3d z = p2 - p3;             // vector from p3 to p2

    vector3d p = z.crossProduct(v12);
    vector3d x = z.crossProduct(v43);
    vector3d y = z.crossProduct(x);

    float u = x.dotProduct(x);
    float v = y.dotProduct(y);

    float result = 360;
    if (u > 0 and v > 0)
    {
        u = p.dotProduct(x) / sqrt(u);
        v = p.dotProduct(y) / sqrt(v);
        if (u != 0 or v != 0)
            result = atan2(v, u) * 180 / PI;
    }

    return result;
}


float calc_phi(const Residue*res_prev, const Residue*res)
{
    if (res_prev == nullptr) return NAN;
    if (res == nullptr) return NAN;

    if (!res_prev->has_atom("C")) return NAN;
    if (!res->has_atom("N"))  return NAN;
    if (!res->has_atom("CA")) return NAN;
    if (!res->has_atom("C"))  return NAN;

    assert(res_prev->_pp_no == res->_pp_no);
    return calc_dihedral_angle(res_prev->get("C").mLoc, res->get("N").mLoc, res->get("CA").mLoc, res->get("C").mLoc);
}

float calc_psi(const Residue*res, const Residue*res_next)
{
    if (res_next == nullptr) return NAN ;
    if (res == nullptr) return NAN ;

    if (!res->has_atom("N"))  return NAN;
    if (!res->has_atom("CA")) return NAN;
    if (!res->has_atom("C"))  return NAN;
    if (!res_next->has_atom("N"))return NAN;

    assert(res_next->_pp_no == res->_pp_no);
    return calc_dihedral_angle(res->get("N").mLoc, res->get("CA").mLoc, res->get("C").mLoc, res_next->get("N").mLoc);
}



float calc_alpha(const Residue*res_prev, const Residue*res, const Residue*res_next, const Residue*res_next_next)
{
    if (res_prev == nullptr      || !res_prev->has_atom("CA")) return NAN;
    if (res == nullptr           || !res->has_atom("CA")) return NAN;
    if (res_next == nullptr      || !res_next->has_atom("CA")) return NAN;
    if (res_next_next == nullptr || !res_next_next->has_atom("CA")) return NAN;

    return calc_dihedral_angle(res_prev->get("CA").mLoc, res->get("CA").mLoc, res_next->get("CA").mLoc, res_next_next->get("CA").mLoc);
}

float calc_cosinus_angle(const vector3d& pa, const vector3d& pb, const vector3d& pc  )
{
    vector3d vba = pa - pb;
    vector3d vbc = pc - pb;
    vba.normalize();
    vbc.normalize();
    return vba.dotProduct(vbc);
}

float calc_cosinus_angle(const vector3d& p1, const vector3d& p2, const vector3d& p3, const vector3d& p4  )
{
    vector3d v12 = p1 - p2;
    vector3d v34 = p3 - p4;

    v12.normalize();
    v34.normalize();
    return v12.dotProduct(v34);
}

float calc_kappa(const Residue*res_prev_prev, const Residue*res, const Residue*res_next_next)
{

    if (res_prev_prev == nullptr || !res_prev_prev->has_atom("CA")) return NAN;
    if (res == nullptr           || !res->has_atom("CA")) return NAN;
    if (res_next_next == nullptr || !res_next_next->has_atom("CA")) return NAN;

    assert(res_prev_prev->has_atom("CA"));
    assert(res->has_atom("CA"));
    assert(res_next_next->has_atom("CA"));

    float c = calc_cosinus_angle(res_prev_prev->get("CA").mLoc, res->get("CA").mLoc, res_next_next->get("CA").mLoc);

    return acos(c) * 180 / PI;;
}

float calc_tco(const Residue*res_prev, const Residue*res)
{
    if (res_prev == nullptr) return NAN;
    if (res == nullptr) return NAN;
    return calc_cosinus_angle(res->get("C").mLoc, res->get("O").mLoc, res_prev->get("C").mLoc, res_prev->get("O").mLoc);
}

DihedralInformation::DihedralInformation(const Residue*res_prev_prev, const Residue*res_prev, const Residue*res, const Residue*res_next, const Residue*res_next_next)
{
    phi = calc_phi(res_prev, res);
    psi = calc_psi(res, res_next);


    tco = calc_tco(res_prev, res);

    alpha = calc_alpha(res_prev, res, res_next, res_next_next);
    kappa = calc_kappa(res_prev_prev, res, res_next_next);
    tr1::tie(theta, d, r) = ThetaD(phi, psi);
}
