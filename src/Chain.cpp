#include "Chain.h"

/**
* Add Residue.
* Add polypeptide if needed.
*/
Residue* Chain::add_residue(Residue&& residue)
{
    assert (residue.has_complete_backbone());
    assert (_polypeptides.size());

    //cerr <<  "Chain::add_residue" << residue << " \n";
    if (!_polypeptides.back().is_covalently_bonded(residue))
    {
        _polypeptides.push_back(Polypeptide());
    }

    residue._polypeptide = &_polypeptides.back();
    _polypeptides.back().add_residue(std::move(residue));
    return _polypeptides.back()._residues.back();
}



string Chain::get_fasta_raw(int row_length, bool indicate_gaps) const
{
    string ret = "";
    bool first = true;

    for (const auto& pp :  _polypeptides)
    {
        string pp_oneletters = "";
        int col_pos = 0;

        if (indicate_gaps && !first)
            ret += "!\n";

        std::for_each(std::begin(pp._residues), std::end(pp._residues), [&pp_oneletters, &row_length, &col_pos](const Residue* res) {
            pp_oneletters += res->_resinfo->mSingleLetter;
            if (++col_pos >= row_length)
            {
                pp_oneletters += "\n";
                col_pos = 0;
            }
        });
        ret += pp_oneletters;
        first = false;
    }
    return ret;
}