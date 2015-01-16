#pragma once

#include "standardlibs.h"
#include "ResidueID.h"
#include "Residue.h"
#include "Chain.h"
#include "HydrogenBond.h"
#include "Turn.h"
#include "BetaBridge.h"
#include "BetaLadder.h"
#include "tools_string.h"

typedef std::pair<Residue*, Residue*> ResPair;
typedef std::vector<std::pair<Residue*, Residue*>> VecResPair;

struct ProteinInfo
{
    string mCompound, mSource, mAuthor;
    string mID = "UNDF";
    std::vector<string> mDbRef;
};

struct hash_vector3i
{
    size_t operator()(const vector3i &x) const{
        return (73856093*x.X) ^ (19349669*x.Y) ^ (83492791*x.Z);
      }
};

struct Protein
{
public:
    static Protein create_from_pdb_file(istream& is, bool cAlphaOnly);

    ~Protein()
    {
    }

    /*Chain& operator[](string chainid) const
    {
    }*/
    void calculate_secondary_structure();
    /**
     * Add an residue
     * Add new Chain if needed
     */
    void add_residue(Residue&& residue);
    void add_missing_hydrogens()
    {
        cerr <<  "add_missing_hydrogens" << endl;
        for (auto chain: _chains)
            chain->add_missing_hydrogens();
    }

    void add_dihedrals()
    {
        cerr <<  "add_dihedrals" << endl;
        for (auto chain: _chains)
            chain->add_dihedrals();
    }

    VecResPair find_neighs_naive();
    VecResPair find_neighs_hash();
    bool assign_single_HB(Residue *res_donor, Residue *res_acceptor);
    void assign_hydrogen_bonds();
    void assign_hydrogen_bonds_from_pairs(VecResPair& pairs);
    void assign_turns();
    void assign_beta_bridges();

    void add_turn(Turn&& turn);
    ProteinInfo mProteinInfo = {};


    string get_fasta(int row_length=50, bool indicate_gaps=true) const;

    void dump_ss_col();
    void dump_ss_dssp();
    void dump_ss_tab();
    void dump_ss_json();
    void dump_ladders();
    void dump_hbonds();
    void dump_hbonds2();
    int residue_counter = 0;
    double asa = -999;

    std::deque<Chain*> _chains = {};
    std::unordered_map<string, Chain> _chains_map = {{}, {}};
    std::deque<HBond> mHbonds = {};
private:
    std::deque<BetaLadder> _ladders = {};
    std::vector<Residue*> _residues = {};
    std::deque<BetaBridge> bridges;
    std::map<int, std::deque<Turn>> _turns = {{3, {}},  {4,  {}}, {5, {}}};


    void collect_secondary_structure();

    void collect_for_dssp_turns();
    void collect_for_dssp_ladders();
    void collect_for_dssp_bulges();

    void collect_for_pssc_turns();
    void collect_for_pssc_strands();

    void collect_for_dihedrals();

    friend ostream &operator<<(ostream &output, const Protein &self)
    {
        output << "#\tmID: "       << self.mProteinInfo.mID << endl;
        output << "#\tmCompound: " << self.mProteinInfo.mCompound << endl;
        output << "#\tmSource: "   << self.mProteinInfo.mSource << endl;
        output << "#\t" << self._chains.size() << " chain(s)";
        for (const auto chain: self._chains)
            output << endl << "#\t\t" <<  *chain;

        return output;
    }
    friend class ASACalculator;
};
