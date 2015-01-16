#include "Atom.h"

Atom Atom::create_hydrogen_from_N(const Atom& iN)
{
    Atom atom;
    atom.mName = "H";
    atom.mElement = "H";
    atom.mType = MAtomType::from_element_string(atom.mElement);
    atom.parent = iN.parent;
    atom.mLoc = iN.mLoc;
    atom.atomclass = AtomClass::hydrogen;
    return atom;
}

AtomResidPair Atom::create_from_pdb_line(string line)
{
        Atom atom;
        //      7 - 11  Integer serial Atom serial number.
        atom.mSerial = boost::lexical_cast<uint32_t>(ba::trim_copy(line.substr(6, 5)));
        //      13 - 16 Atom name Atom name.
        atom.mName = ba::trim_copy(line.substr(12, 4));
        atom.altLoc = ba::trim_copy(line.substr(16, 1));

        //      18 - 20 Residue name resName Residue name.
        //      22              Character chainID Chain identifier.
        //      23 - 26 Integer resSeq Residue sequence number.
        //      27              AChar iCode Code for insertion of residues.
        auto resid = ResidueID::from_PDB_line(line);

        //      31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms.
        atom.mLoc.X = std::stod(line.substr(30, 8));
        //      39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms.
        atom.mLoc.Y = std::stod(line.substr(38, 8));
        //      47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms.
        atom.mLoc.Z = std::stod(line.substr(46, 8));
        //      55 - 60 Real(6.2) occupancy Occupancy.
        atom.mOccupancy = std::stod(line.substr(54, 6));
        //      61 - 66 Real(6.2) tempFactor Temperature factor.
        atom.mTempFactor = std::stod(line.substr(60, 6));
        //      77 - 78 LString(2) element Element symbol, right-justified.
        if (line.length() > 76) {
            atom.mElement = ba::trim_copy(line.substr(76, 3));
            atom.mType = MAtomType::from_element_string(atom.mElement);
        }
        //      79 - 80 LString(2) charge Charge on the atom.
        atom.mCharge = 0;

        if (k_radii_for_asa.count(atom.mName))
        {
            atom.radius_for_asa = k_radii_for_asa.at(atom.mName);
            // TODO: for pdbs containing hydrogens this is wrong!
            atom.atomclass = AtomClass::backbone;
        } else {
            atom.radius_for_asa = k_radius_for_asa_sidechain;
            atom.atomclass = AtomClass::sidechain;
        }
        return std::make_pair(atom, resid);
}