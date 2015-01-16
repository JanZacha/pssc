#pragma once

#include "standardlibs.h"
#include "ResidueID.h"
#include "Residue.h"
const double kMaxPeptideBondLength = 2.5;

struct Polypeptide
{
    static unsigned pp_count_ ;
    Polypeptide():
    _number(pp_count_++),
    _residues({}),
    _residues_map({}),
    _closed(false)
    {
    }

    /**
     * Add Residue.
     */
    Residue* add_residue(Residue&& residue)
    {
        //cerr <<  "Polypeptide::add_residue " <<  residue <<  endl;
        string id = residue.ID.shortstr();
        residue.set_index(_number, _residues.size()+1);
        _residues_map[id] = std::move(residue);
        _residues.push_back(&_residues_map[id]);
        return &_residues_map[id];
    }

    /**
     * Return true if residue is covalently bonded to this Polypeptide or if this Polypeptide is empty.
     * If Polypeptide is closed return false.
     * If residue is covalently bonded, this Polypeptide is not the right place for it.
     */
    bool is_covalently_bonded(const Residue& residue);
    void close() {_closed = true;}
    std::size_t size() const
    {
        return _residues.size();
    }

    void add_missing_hydrogens();
    void add_dihedrals();
    int _number;
    std::vector<Residue*> _residues;
    std::unordered_map<string, Residue> _residues_map;
    bool _closed;

private:

    friend ostream &operator<<(ostream &output, const Polypeptide &self)
    {
        output << "Polypeptide " << self._number << "(";
        for (const auto residue: self._residues)
            //residue->mID.shortstr()
            output << *residue <<" ";
        output <<  ")";
        return output;
    }
};