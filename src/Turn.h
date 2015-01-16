#pragma once
#include "standardlibs.h"
#include "ResidueID.h"
#include "HydrogenBond.h"
#include "Residue.h"


struct Turn
{
    HBond *mHBond;
    int n;
    bool valid = false;

    Turn(HBond* hb):
        mHBond(hb)
    {
        assert(hb->mDonor._pp_no == hb->mDonor._pp_no);
        n = hb->mDonor._res_no - hb->mAcceptor._res_no;
    }

    static bool check_for_turn(const HBond& hb, int& outN)
    {
        if (hb.mDonor._pp_no != hb.mAcceptor._pp_no)
            return false;
        outN = hb.mDonor._res_no - hb.mAcceptor._res_no;
        return (outN <= 6 && outN >= -6);
    }

private:
    friend ostream &operator<< (ostream &output, const Turn &t)
    {
            output << t.n << "-Turn " << "(" << " " << *t.mHBond << ")";
            return output;
    }
};