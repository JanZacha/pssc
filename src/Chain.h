#pragma once

#include "standardlibs.h"
#include "ResidueID.h"
#include "Residue.h"
#include "Polypeptide.h"

struct Chain
{
    string _id;

    Chain(string id = ""):
    _id(id),
    _polypeptides({Polypeptide()}),
    _terminated(false)
    {}

    /**
     * Add Residue.
     * Add polypeptide if needed.
     */
    Residue* add_residue(Residue&& residue);

    /**
     * handle gap. close last Polypeptide
     */
    void recognize_gap()
    {
        _polypeptides.back().close();
    }

    void add_missing_hydrogens()
    {
        for (auto& pp: _polypeptides)
        {
            pp.add_missing_hydrogens();
        }
    }

    void add_dihedrals()
    {
        for (auto& pp: _polypeptides)
        {
            pp.add_dihedrals();
        }
    }

    void terminate()
    {
        _terminated = true;
    }
    bool has_been_terminated()
    {
        return _terminated;
    }

    size_t get_num_residues() const
    {
        size_t ret = 0;
        for (auto p: _polypeptides)
            ret += p.size();
        return ret;
    }

    string get_fasta_raw(int row_length=50, bool indicate_gaps=true) const;

    deque<Polypeptide> _polypeptides;

private:

    bool _terminated;
    friend ostream &operator<<(ostream &output, const Chain &self)
    {

        output << "Chain " << self._id << ", size " << self.get_num_residues() << " residues in " << self._polypeptides.size() << " polypeptide(s)";
        //for (const auto pp: self._polypeptides)
        //    output << pp <<  " ";
        //output <<  ")";
        return output;
    }
};