#pragma once

#include "standardlibs.h"
#include "ResidueID.h"
#include "Atom.h"
#include "SSInformation.h"
#include "DihedralInformation.h"
#include "pssc_config.h"

struct Polypeptide;
struct HBond;
struct Turn;
struct BetaBridge;

struct Residue
{
public :
    Residue()
    {
    }
    bool has_betabridge_assignend = false;
    Residue(ResidueID id):
        ID(id),
        _resinfo(ChemicalPropertiesResidue::from_three_letter(id.resname)),
        _polypeptide(nullptr)
    {
        //cerr << "Residue (" <<  id << ")" <<  (resinfo->mName) <<  ":" <<  int(resinfo->mRestype) << endl;
    }

    /*
     * should be like this. but then copy constructor needs to be defined too
     *
    Residue(Residue && tmp):
    ID(tmp.ID),
    _resinfo(tmp._resinfo),
    _polypeptide(tmp._polypeptide)
    {
    cerr << "MOVE";
    for (auto& kv : _atoms_map)
    {
    // moving a residue invalidates internal parent reference
    kv.second.parent = this;
}
}*/

    std::vector<Residue*>::iterator get_iter();

    void set_index(int pp_no, int res_no)
    {
        _pp_no = pp_no;
        _res_no = res_no;
    }

    void add_atom(Atom&& atom)
    {
        assert(_resinfo);
        // moving a residue invalidates internal parent reference
        // all parent links are reset in the end, externally (in class Protein)!
        atom.parent = this;
        _atoms_map[atom.mName] = std::move(atom);
    }

    void finalize()
    {
        if (_atoms_map.count("N"))
            N = &_atoms_map["N"];
        if (_atoms_map.count("CA"))
            CA = &_atoms_map["CA"];
        if (_atoms_map.count("C"))
            C = &_atoms_map["C"];
        if (_atoms_map.count("O"))
            O = &_atoms_map["O"];
    }

    const Atom& operator[](string iAtomName)
    {
        return _atoms_map[iAtomName];
    }

    const Atom& get(string iAtomName) const
    {
        return _atoms_map.at(iAtomName);
    }

    Atom& get(string iAtomName)
    {
        return _atoms_map.at(iAtomName);
    }

    bool has_atom(string iAtomName) const
    {
        return _atoms_map.count(iAtomName);
    }

    bool has_complete_backbone() const
    {
        return has_atom("N") && has_atom("CA") && has_atom("C") && has_atom("O");
    }

    std::size_t size() const
    {
        return _atoms_map.size();
    }

    std::deque<HBond*> get_sorted_HBond_donating();
    std::deque<HBond*> get_sorted_HBond_accepting();

    float get_total_dist_to(const Residue& other) const
    {
       return other._total_res_no - _total_res_no;
    }

    std::size_t get_sequence_dist(const Residue& other) const
    {
        if (_pp_no == other._pp_no)
            return abs(_res_no - other._res_no);
        return 9999999; //TODO: std maxint oder wie es heisst.;
    }

    std::size_t get_sequence_dist(const Residue* other) const
    {
        assert(_res_no);
        assert(other->_res_no);
        if (_pp_no == other->_pp_no)
            return abs(_res_no - other->_res_no);
        return 9999999; //TODO: std maxint oder wie es heisst.;
    }

    bool operator<(const Residue& other) const
    {
        if (_pp_no == other._pp_no)
            return _res_no < other._res_no;
        return _pp_no < other._pp_no;
    }

    bool operator>=(const Residue& other) const
    {
       return !(operator<(other));
    }

    string get_vmd_sel() const
    {
        return "chain " + ID.chain + " and " + "resid " + std::to_string(ID.seqNumber) + ID.insertionCode;;
    }

    ResidueID ID;
    ChemicalPropertiesResidue* _resinfo = nullptr;
    double asa = -999;
    string get_ssi_output_dssp();
    string get_ssi_output_col();
    string get_ssi_output_tab();
    string get_ssi_output_json();

    SSInformation ssi;
    DihedralInformation dhi;

    int _pp_no;
    int _res_no;
    int _total_res_no;
    Polypeptide *_polypeptide;

    std::deque<HBond*> HBond_donating  = {};
    std::deque<HBond*> HBond_accepting = {};

    std::map<int, Turn*> _turns_starting = {{3, nullptr},  {4,  nullptr}, {5, nullptr}};

    bool betabridged_b = false;
    bool betabridged_r = false;

    Atom* N = nullptr;
    Atom* CA = nullptr;
    Atom* C = nullptr;
    Atom* O = nullptr;

private:
    std::map<string, Atom> _atoms_map = {};
    friend struct Protein;
    inline friend std::ostream& operator<<(std::ostream &output, const Residue &res)
    {
        assert(res._resinfo);
        output << res.ID.chain  << ":";
        if (res._resinfo)
            output << res._resinfo->mThreeLetter;
        else
            output << "res???";
        output << res.ID.seqNumber << res.ID.insertionCode;

        output << "(" <<  res._pp_no <<  ", " << res._res_no <<  ")";
        /*
        output << "(" << res.size() << ")" << "[";
        for (const auto atom: res._atoms_map)
            output << atom.second.mName <<  ",";
        output << "]";
        */
        return output;
    }
    friend class ASACalculator;
};



