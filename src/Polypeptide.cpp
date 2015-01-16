#include "Polypeptide.h"

    /**
     * Return true if residue is covalently bonded to this Polypeptide or if this Polypeptide is empty.
     * If Polypeptide is closed return false.
     * If residue is covalently bonded, this Polypeptide is not the right place for it.
     */

unsigned Polypeptide::pp_count_ = 1;

bool Polypeptide::is_covalently_bonded(const Residue& residue)
{
    //cerr <<  "Polypeptide::is_covalently_bonded " <<  _residues.size() <<"\n";
    if (!size())
        return true;

    if (_closed)
        return false;

    const Residue& previous = *(_residues.back());
    //cerr <<  "previous " << _residues.back()->ID<<  previous.ID <<  previous.size()<< endl;
    assert (previous.has_atom("C"));
    assert (residue.has_atom("N"));
    return (previous.get("C").mLoc - residue.get("N").mLoc).getLength() <= kMaxPeptideBondLength;
}

void Polypeptide::add_missing_hydrogens()
{
    if (_residues.size() <= 2) return;
    Residue* previous;
    auto it = _residues.begin() + 1;
    std::for_each(it, std::end(_residues), [&it, &previous](Residue* residue)
    {
        previous = *std::prev(it++);
        if (residue->_resinfo->mRestype != ResidueType::withPolarHydrogen) return;
        //bool hydrogen_missing = false;
        if (!residue->has_atom("H"))
        {
            residue->add_atom(Atom::create_hydrogen_from_N( residue->get("N")) );
           // hydrogen_missing = true;
        }
        vector3d vec_oc = previous->get("C").mLoc - previous->get("O").mLoc;
        vector3d interpolated_pos = residue->get("N").mLoc + vec_oc.normalize();
        /*if (!hydrogen_missing)
        {
            double diff = (interpolated_pos - residue->get("H").mLoc).getLength();
            cerr <<  "diff: "<< diff <<  " ";
        }*/
        residue->get("H").mLoc = interpolated_pos;

    });
}


void Polypeptide::add_dihedrals()
{
    if (_residues.size() <= 2) return;

    std::size_t i = 0;
    Residue* res_prev_prev;
    Residue* res_prev;
    Residue* res_next;
    Residue* res_next_next;
    std::for_each(std::begin(_residues), std::end(_residues), [&](Residue* residue)
    {

        res_prev_prev = i > 1 ? _residues[i-2] : nullptr;
        res_prev      = i > 0 ? _residues[i-1] : nullptr;
        res_next      = i < -1 + _residues.size() ? _residues[i+1] : nullptr;
        res_next_next = i < -2 + _residues.size() ? _residues[i+2] : nullptr;
        residue->dhi = DihedralInformation(res_prev_prev, res_prev, residue, res_next, res_next_next);
        ++i;
    });
/*
    it = std::begin(_residues);
    (*it)->dhi = DihedralInformation(nullptr, *it, *std::next(it));

    it = std::end(_residues) - 1;
    (*it)->dhi = DihedralInformation(*std::prev(it), *it, nullptr);*/
}

