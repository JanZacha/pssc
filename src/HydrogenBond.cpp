#include "HydrogenBond.h"
#include "Residue.h"
hbond_energy_callback_function HBond::_calculate_hbond_energy_callback = nullptr;

bool HBond_from_to(const Residue* donor, const Residue* acceptor)
{
   //cerr <<  "HBond_from_to " <<  donor->HBond_donating.size() <<  endl;
    for (const HBond* hb: donor->HBond_donating)
    {
        assert(hb);
        //cerr <<  ">>>>" <<  (*hb) <<  endl;
       // if (!(hb->mAcceptor.ID != acceptor->ID))
        if (hb->mAcceptor._total_res_no == acceptor->_total_res_no)
            return true;
    }
    return false;
}

bool HBond_from_to(const Residue* donor, const Residue* acceptor, HBond*& outHBond)
{
    for (HBond* hb: donor->HBond_donating)
    {
        assert(hb);
        //if (!(hb->mAcceptor.ID != acceptor->ID))
        if (hb->mAcceptor._total_res_no == acceptor->_total_res_no)
        {
            outHBond = hb;
            return true;
        }
    }
    outHBond = nullptr;
    return false;
}


double HBond::calculate_energy_estat(const Residue& donor, const Residue& acceptor)
{
    double distanceHC = donor.get("H").mLoc.getDistanceFrom(acceptor.get("C").mLoc);
    double distanceHO = donor.get("H").mLoc.getDistanceFrom(acceptor.get("O").mLoc);

    double distanceNC = donor.get("N").mLoc.getDistanceFrom(acceptor.get("C").mLoc);
    double distanceNO = donor.get("N").mLoc.getDistanceFrom(acceptor.get("O").mLoc);

    double distanceCCA = donor.get("CA").mLoc.getDistanceFrom(acceptor.get("C").mLoc);
    double distanceOCA = donor.get("CA").mLoc.getDistanceFrom(acceptor.get("O").mLoc);

    // coulomb energies
    double qh  =   0.31;
    double qn  =  -0.47;
    double qc  =   0.51;
    double qo  =  -0.51;
    double qca =   0.16;
    double f   = 332.0;
    double qNO = qo*qn;
    double qHO = qo*qh;
    double qNC = qc*qn;
    double qHC = qc*qh;
    double qCCA = qc*qca;
    double qOCA = qo*qca;
    return f*(qHO/distanceHO + qHC/distanceHC + qNC/distanceNC + qNO/distanceNO + qCCA/distanceCCA + qOCA/distanceOCA);
}

double HBond::calculate_energy_dssp(const Residue& donor, const Residue& acceptor)
{
        double distanceHC = donor.get("H").mLoc.getDistanceFrom(acceptor.get("C").mLoc);
        double distanceHO = donor.get("H").mLoc.getDistanceFrom(acceptor.get("O").mLoc);
        double distanceNC = donor.get("N").mLoc.getDistanceFrom(acceptor.get("C").mLoc);
        double distanceNO = donor.get("N").mLoc.getDistanceFrom(acceptor.get("O").mLoc);
        double kCouplingConstant = -27.888;
        return kCouplingConstant / distanceHO - kCouplingConstant / distanceHC + kCouplingConstant / distanceNC - kCouplingConstant / distanceNO;
}

double HBond::calculate_energy_test(const Residue& donor, const Residue& acceptor)
{
        return -1;
}

ostream &operator<<(ostream &output, const HBond &hb)
{
    output << "HB " << hb.mDonor.ID <<  "->" <<  hb.mAcceptor.ID << "="
        << (boost::format("%.1f") % hb.mEnergy ) <<  " ";
    return output;
}