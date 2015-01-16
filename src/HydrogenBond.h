#pragma once

#include "standardlibs.h"
#include <boost/format.hpp>

struct Residue;

typedef double (*hbond_energy_callback_function)(const Residue& donor, const Residue& acceptor);

struct HBond
{
    HBond(Residue& inDonor, Residue& inAcceptor):
        mDonor(inDonor),
        mAcceptor(inAcceptor)
    {
        mEnergy = _calculate_hbond_energy_callback(mDonor, mAcceptor);
    }

    static hbond_energy_callback_function _calculate_hbond_energy_callback;
    static double calculate_energy_estat(const Residue& donor, const Residue& acceptor);
    static double calculate_energy_dssp(const Residue& donor, const Residue& acceptor);
    static double calculate_energy_test(const Residue& donor, const Residue& acceptor);

    double mEnergy;
    Residue& mDonor;
    Residue& mAcceptor;
    double ca_dist;

private:
    friend ostream &operator<<(ostream &output, const HBond &hb);

};

inline bool operator< ( HBond lhs,  HBond rhs){return lhs.mEnergy < rhs.mEnergy;}
inline bool operator> (const HBond& lhs, const HBond& rhs){return rhs < lhs;}
inline bool operator<=(const HBond& lhs, const HBond& rhs){return !(lhs > rhs);}
inline bool operator>=(const HBond& lhs, const HBond& rhs){return !(lhs < rhs);}

bool HBond_from_to(const Residue* donor, const Residue* acceptor);
bool HBond_from_to(const Residue* donor, const Residue* acceptor, HBond*& outHBond);
