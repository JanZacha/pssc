#pragma once

#include "standardlibs.h"
#include "ResidueID.h"
#include "ChemicalData.h"

struct Atom;
struct Residue;
typedef std::pair<Atom, ResidueID> AtomResidPair;

const std::map<string, double> k_radii_for_asa = {
    {"N",  1.65},
    {"CA",  1.87},
    {"C",  1.76},
    {"O",  1.4}
};
const double k_radius_for_asa_sidechain = 1.8;
const double k_radius_for_asa_water = 1.4;
const double k_radius_for_asa_max = 1.87;

enum class AtomClass : int {backbone, sidechain, hydrogen};

struct Atom
{
    static AtomResidPair create_from_pdb_line(string record_data);
    static Atom create_hydrogen_from_N(const Atom& iN);
    int        mSerial = -1;
    string     mName   = "";
    string     altLoc = "";

    MAtomType  mType   = MAtomType();

    vector3d   mLoc        = vector3d(0, 0, 0);
    double     radius_for_asa = 0;
    double     mOccupancy  = 0;
    double     mTempFactor = 0;
    string     mElement    = "";
    int        mCharge     = 0;
    double     asa = -999;
    Residue*   parent      = nullptr;
    AtomClass  atomclass;
private:
    friend ostream &operator<<  ( ostream &output, const Atom &c )
    {
            //output << c.mSerial << " " <<  c.mName<<" "<< c.mLoc.X<< ", "   << c.mLoc.Y<< ", "   << c.mLoc.Z     ;
            output << c.mSerial << " " <<  c.mName << " " << (int) c.asa;
            return output;
    }
};
