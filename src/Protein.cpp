#include "Protein.h"

#include <fstream>
#include <algorithm>
#include <boost/tr1/tuple.hpp>

const double max_HB_energy = -.5;
const double kMinimalCADistance = 9.0;
namespace tr1 = std::tr1;

/**
* Add an residue
* Add new Chain if needed
*/
void Protein::add_residue(Residue&& residue)
{
    ;
    //cerr <<  "Protein::add_residue " << residue << endl;

    if (residue.has_complete_backbone())
    {
        residue._total_res_no = ++residue_counter;
        if (!_chains_map.count(residue.ID.chain))
        {
            _chains_map[residue.ID.chain] = Chain(residue.ID.chain);
            _chains.push_back(&_chains_map[residue.ID.chain]);
        }
        if (_chains_map[residue.ID.chain].has_been_terminated())
        {
            //TODO: raise some Exception
            cerr << "chain was already terminated";
        }

        Residue *res = _chains_map[residue.ID.chain].add_residue(std::move(residue));
        //cerr <<  "Protein::add_residue " << res << " " <<(int)res->has_complete_backbone() ;
        for (auto& kv : res->_atoms_map)
        {
            //      This is extremly unelegant. Should be implemented in the move constructor I guess
            // moving a residue invalidates internal parent reference
            kv.second.parent = res;
        }
        _residues.push_back(res);
    }
    else
    {
        residue._total_res_no = -1;
        //TODO: do some logging
        //cerr << "incomplete residue";
        if (_chains_map.count(residue.ID.chain))
        {
            _chains_map[residue.ID.chain].recognize_gap();
        }
    }

}

void Protein::assign_hydrogen_bonds_from_pairs(VecResPair& pairs)
{
    cerr << "assign_hydrogen_bonds_from_pairs (" << pairs.size() << " candidates)" << endl;
    int added = 0;
    for (auto pair: pairs)
    {
        //cerr << *pair.first << *pair.second << endl;
        //cerr << (int)pair.first->_resinfo->mRestype  << " - " <<(int)pair.second->_resinfo->mRestype << endl;
        if (pair.first->_resinfo->mRestype == ResidueType::withPolarHydrogen)
        added += assign_single_HB(pair.first, pair.second);
        if (pair.second->_resinfo->mRestype == ResidueType::withPolarHydrogen)
        added += assign_single_HB(pair.second, pair.first);
    }
    cerr << "assigned " << added << " hydrogen bonds" << endl;
}


void Protein::assign_hydrogen_bonds()
{
    cerr << "assign_hydrogen_bonds" << endl;

    auto pairs = find_neighs_hash();
   // auto pairs = find_neighs_naive();
    assign_hydrogen_bonds_from_pairs(pairs);

    return;
/*
 * can be removed now. was just left here to check if neighbor via hash algorithm does not make any mistakes.

    typedef std::pair<vector3d, Residue*> posResPair;
    std::vector<posResPair> posPairs = {};
    posPairs.resize(_residues.size());

    auto makePosResPair = [](Residue *res) {
        return std::make_pair(res->get("CA").mLoc, res);
    };

    std::transform(std::begin(_residues), std::end(_residues), std::begin(posPairs), makePosResPair);

    int i = 0;
    int j = 0;
    for (posResPair pair_i: posPairs)
    {
        ++i;
        j = 0;
        for (posResPair pair_j: posPairs)
        {
            ++j;
            if (j > i)
                break;

            if (pair_i.first.getDistanceFrom(pair_j.first) < kMinimalCADistance)
            {
                if (pair_i.second->_resinfo->mRestype == ResidueType::withPolarHydrogen)
                    assign_single_HB(pair_i.second, pair_j.second);
                if (pair_j.second->_resinfo->mRestype == ResidueType::withPolarHydrogen)
                    assign_single_HB(pair_j.second, pair_i.second);
            }
        }
    };*/
}


const double MAXDIST = 3.67;

VecResPair Protein::find_neighs_naive()
{
    cerr << "find_neighs_naive " << _residues.size() <<endl;
    int tries = 0;
    VecResPair ret;
    for (auto res1: _residues)
    for (auto res2: _residues)
    {
        tries += 1;
        if (res1 == res2) continue;

        if (!res1->has_atom("CA") || !res2->has_atom("CA") )
            continue;
        const Atom ca1 = (*res1)["CA"];
        const Atom ca2 = (*res2)["CA"];
        double dist = ca1.mLoc.getDistanceFrom(ca2.mLoc);
        if (dist < kMinimalCADistance)
        {
            //cerr << ca1 << " -- " << ca2 << " " << dist << endl;
            ret.push_back(std::make_pair(res1, res2));
        }

    }
    cerr << "found " << ret.size() << endl;
    cerr  << "in " << tries << " tries " << endl;
    return ret;
}

vector3i grid_vec(const vector3d &v)
{
    double fact = kMinimalCADistance;
    return vector3i(v.X/fact, v.Y/fact, v.Z/fact);
}



/**
 * see "Optimized Spatial Hashing for Collision Detection of Deformable Objects"
 * Matthias Teschner Bruno Heidelberger Matthias M¨uller Danat Pomeranets Markus Gross
Computer Graphics Laboratory
ETH Zurich
http://www.beosil.com/download/CollisionDetectionHashing_VMV03.pdf
 */
VecResPair Protein::find_neighs_hash()
{
    cerr << "find_neighs_hash"<<endl;
    int tries = 0;
    VecResPair ret;
    std::unordered_map<vector3i,std::vector<const Atom*>,hash_vector3i> my_map;
    my_map.reserve(5000);
    for (auto res: _residues)
    {
        if (!res->has_atom("CA"))
            continue;
        const auto& ca1 = ((*res)["CA"]);
        vector3i ca1_pos = grid_vec(ca1.mLoc);

        my_map[ ca1_pos ].push_back(&((*res)["CA"]));
    }
    std::vector<vector3i> connected_fields = {
        {-1, -1, 0}, {0, -1, 0}, {1,-1, 0},
        {-1,  0, 0}, {0,  0, 0}, {1, 0, 0},
        {-1,  1, 0}, {0,  1, 0}, {1, 1, 0},

        {-1, -1, 1}, {0, -1, 1}, {1,-1, 1},
        {-1,  0, 1}, {0,  0, 1}, {1, 0, 1},
        {-1,  1, 1}, {0,  1, 1}, {1, 1, 1},

        {-1, -1,-1}, {0, -1,-1}, {1,-1,-1},
        {-1,  0,-1}, {0,  0,-1}, {1, 0,-1},
        {-1,  1,-1}, {0,  1,-1}, {1, 1,-1}
      };

    for (auto res: _residues)
    {
        if (!res->has_atom("CA"))
            continue;
        const auto& ca1 = ((*res)["CA"]);
        vector3i ca1_pos = grid_vec(ca1.mLoc);
        for (const auto& con : connected_fields)
        {
            vector3i grid_pos = ca1_pos + con;
            tries += my_map[grid_pos].size();

            for (auto ca2p :  my_map[grid_pos])
            {
                if (ca2p == &(*res)["CA"])
                    continue;
                double dist = ca1.mLoc.getDistanceFrom(ca2p->mLoc);
                if (dist < kMinimalCADistance)
                {
                    ret.push_back(std::make_pair(res, ca2p->parent));
                }
            }
        }
    }
    cerr  << "needed " << tries << " tries " << endl;
    return ret;
}

/**
 * return true if Hb was added
 */
bool Protein::assign_single_HB(Residue *res_donor, Residue *res_acceptor)
{
    if (!res_donor->has_atom("H")) return false;

    const std::size_t diff = res_donor->get_sequence_dist(*res_acceptor);
    assert (diff >= 0);
    if (res_donor == res_acceptor)
        assert(diff == 0);

    if (diff < 2)
        return false;

    HBond hb = HBond(*res_donor, *res_acceptor);
    if (hb.mEnergy > max_HB_energy)
        return false;
    mHbonds.push_back(hb);
    res_donor->HBond_donating.push_back(&(mHbonds.back()));
    res_acceptor->HBond_accepting.push_back(&(mHbonds.back()));
    return true;
}

void Protein::assign_turns()
{
    cerr << "assign_turns" << endl;
    for (Chain* chain: _chains)
    {
        for (Polypeptide pp: chain->_polypeptides)
        {
            for (Residue* res: pp._residues)
            {
                int n;
                for (HBond* hb : res->HBond_accepting)
                {
                    bool val = Turn::check_for_turn(*hb, n);

                    if (val) {
                        add_turn(Turn(hb));
                    }
                }
            }
        }
    }
}

void Protein::add_turn(Turn&& turn)
{
    _turns[turn.n].push_back(std::move(turn));

    turn.mHBond->mAcceptor._turns_starting[turn.n] = &_turns[turn.n].back();
}

// TODO: this needs to be revised. a quadratic loop is not necessary.
// why not use the already assigned hydrogen bonds directly?
void Protein::assign_beta_bridges()
{
    cerr <<  "assign_beta_bridges " << _residues.size() <<  endl;
    if (_residues.size() < 4) return;
    auto it_outer = std::begin(_residues) + 1;
    if(!_residues.size()) return;
    std::for_each(it_outer, _residues.end() - 1, [&](Residue* element)
    {
        //cerr << **it_outer <<  " ";
        auto it_inner = std::begin(_residues) + 3;
        std::for_each(it_inner, _residues.end() - 1, [&](Residue* element)
        {
           // cerr << **it_outer << " to "<< **it_inner << endl;
            auto res_a = *(it_outer-1);
            auto res_b = *(it_outer);
            auto res_c = *(it_outer+1);
            auto res_q = *(it_inner-1);
            auto res_r = *(it_inner);
            auto res_s = *(it_inner+1);
            ++it_inner;

            if (!(*res_b < *res_r)) return;
            BetaBridge bridge = BetaBridge::create_from_residues(res_a, res_b, res_c, res_q, res_r, res_s);
            if (bridge.mType != BetaBridgeType::none ) {
                bridges.push_back(bridge);
                bridge.res_b->betabridged_b = true;
                bridge.res_r->betabridged_r = true;
                //bridge.res_b->has_betabridge_assignend = true;
            }
        });

      it_outer++;
    });

    cerr << "found " << bridges.size() << " bridges.\n";
    /*for (BetaBridge bridge: bridges)
    {
        cout <<  bridge <<  endl;
    }*/

    deque<BetaLadder>& ladders = _ladders;

    for (auto &bridge: bridges)
    {
        //cerr <<  "!,  " <<  bridge<< endl;
        bool found = false;
        for (auto &ladder: ladders)
        {
            if (ladder.can_be_extended(&bridge))
            {
                //cerr <<  ladder << " can_be_extended by " << bridge << "\n";
                ladder.append_bridge(&bridge);
                found = true;
                break;
            }
        }
        if (!found)
        {
            BetaLadder ladder = BetaLadder(&bridge, ladders.size() + 1);
            ladders.push_back(ladder);
        }
    }

    /**
    std::ofstream outfile("/home/zacha/ladders.tcl", std::ios_base::out | std::ios_base::binary);
    cerr <<  "######\nLadders:"<< endl;
    for (BetaLadder &ladder: ladders)
    {
        cerr << "- " <<  ladder << " "  << ladder.get_marker() <<  endl;
        outfile <<  "mol selection backbone and (";
        bool first = true;
        for (BetaBridge* b: ladder.mBridges)
        {
            if (!first)
                outfile <<  " or ";
            first = false;
            cerr <<  "\t" << *b->res_b << "  " <<   *b->res_r <<  endl;
            outfile << b->res_b->get_vmd_sel() << " ";
        }
        outfile << ")" <<  endl;
        outfile << "mol addrep [molinfo top]" <<  endl;


        outfile <<  "mol selection backbone and (";
        first = true;
        for (BetaBridge* b: ladder.mBridges)
        {
            if (!first)
                outfile <<  " or ";
            first = false;
            outfile << b->res_r->get_vmd_sel() << " ";
        }
        outfile << ")" <<  endl;
        outfile << "mol addrep [molinfo top]" <<  endl;

    }
    **/
}

Protein Protein::create_from_pdb_file(istream& instream, bool cAlphaOnly)
{
        cerr << "create_from_pdb_file" << endl;
        Protein protein;

        string line;

        bool model = false;

        bool break_demanded = false;
        bool altlocs_present = false;
        string selected_altloc;

        Residue current_residue = Residue(ResidueID());

        auto atom_function = [&](string record_data)
        {
            auto atom_resid = Atom::create_from_pdb_line(line);
            if ((!altlocs_present) && (atom_resid.first.altLoc != ""))
            {
                altlocs_present = true;
                if (selected_altloc == "")
                {
                    selected_altloc = atom_resid.first.altLoc;
                    cerr <<  "Found altLoc '" <<  selected_altloc <<  "' will only use this (no selection was given)";
                }
            }

            if (altlocs_present && atom_resid.first.altLoc != "" and atom_resid.first.altLoc != selected_altloc)
                {
                    cerr << "skipping alternate atom record '"<<  atom_resid.first.altLoc  << "' for " << atom_resid.first << endl;
                    return;
                }

            if (atom_resid.second != current_residue.ID)
            {
                if (current_residue.size())
                {
                    current_residue.finalize();
                    protein.add_residue(std::move(current_residue));
                }
                current_residue = Residue(atom_resid.second);
            }
            current_residue.finalize();
            current_residue.add_atom(std::move(atom_resid.first));
        };

        std::map<std::string, std::function<void(string)>> record_functions = {

                {"HEADER", [&](string record_data) {
                    if (record_data.length() >= 56)
                        protein.mProteinInfo.mID = record_data.substr(52, 4);
                }},

                {"COMPND", [&](string record_data) {
                    protein.mProteinInfo.mCompound += record_data;
                }},

                {"SOURCE", [&](string record_data) {
                                protein.mProteinInfo.mSource += record_data;
                }},

                {"AUTHOR", [&](string record_data) {
                                protein.mProteinInfo.mAuthor += line.substr(10);
                }},

                {"DBREF", [&](string record_data) {
                                protein.mProteinInfo.mDbRef.push_back(line);
                }},

                {"MODEL", [&](string record_data) {
                        /* brain dead support for only the first model in the file (NMR) */
                        model = true;
                }},

                {"ENDMDL", [&](string record_data) {
                        assert(model);
                        break_demanded = true;
                }},

                {"SSBOND", [&](string record_data) {
                        /**
                         * don't need them actually
                         **/
                }},

                {"ATOM", atom_function},
                {"HETATM", atom_function},
                /*
                // add ATOMs only if the chain isn't terminated
                if (terminatedChains.count(line[21]))
                        continue;*/

                {"TER", [&](string record_data) {

                }},

        };

        while (std::getline(instream, line) and not break_demanded)
        {
            if (line.size() < 10)
                continue;
            string record_type = line.substr(0,6);
            string record_data = line.substr(10);
            trim(record_type);
            trim(record_data);
            if (!record_functions.count(record_type) )
                continue;
            record_functions[record_type](record_data);
        }

        if (current_residue.size())
            protein.add_residue(std::move(current_residue));


        return protein;
}

void Protein::calculate_secondary_structure()
{
    add_dihedrals();
    add_missing_hydrogens();
    assign_hydrogen_bonds();
    assign_turns();
    assign_beta_bridges();
    collect_secondary_structure();
}

void Protein::collect_for_dssp_turns()
{
  // This is only for pseudo-DSSP output. For easy comparison only.
    if (_residues.size()<4) return;
    for (std::size_t i = 0; i < _residues.size() - 2; ++i)
    {

        Residue* res = _residues[i];
        Residue* res_next = _residues[i+1];

        const std::vector<std::size_t> ns = {3, 4, 5};
        for (std::size_t n: ns)
        {
            if (res->_turns_starting[n] && res_next->_turns_starting[n])
            {
                res_next->_turns_starting[n]->valid = true;
            }

            if (res->_turns_starting[n])
            {
                res->_turns_starting[n]->mHBond->mDonor.ssi._pseudo_dssp_full[n-3] = '<';

                if (res->ssi._pseudo_dssp_full[n-3] == '<')
                {
                    res->ssi._pseudo_dssp_full[n-3] = 'X';
                }
                else
                    res->ssi._pseudo_dssp_full[n-3] = '>';

               for (std::size_t j = 0; j <= n; ++j)
               {
                    if (i + j >= _residues.size())
                        break;
                    Residue *res_j = _residues[i + j];

                    if (j < n)
                    {
                        if (res->_turns_starting[n]->valid)
                        {
                            char GHI = 'G' + (n-3);
                            res_j->ssi.dssp_one_csi = GHI;
                        }
                    }
                    if (j && j < n)
                    {
                        if (res_j->ssi.dssp_one_csi == "-")
                            res_j->ssi.dssp_one_csi = 'T';
                    }
                    if (
                        res_j->ssi._pseudo_dssp_full[n-3] != '>' &&
                        res_j->ssi._pseudo_dssp_full[n-3] != 'X' &&
                        res_j->ssi._pseudo_dssp_full[n-3] != '<'
                    )
                        res_j->ssi._pseudo_dssp_full[n-3] = '0' + n;
                }
            }
        }
    }
}

void Protein::collect_for_pssc_turns()
{
    // TODO: res_next->_turns_starting[n]->valid is set by collect_for_dssp_turns

    const std::vector<std::size_t>   ns = {3, 4, 5};
    std::unordered_map<int, string> GHI = {{3, "G"}, {4, "H"}, {5, "I"}};
    std::unordered_map<int, string> ghi = {{3, "g"}, {4, "h"}, {5, "i"}};

    for (auto res: _residues)
    {
        res->ssi.values["G"] = ".";
        res->ssi.values["H"] = ".";
        res->ssi.values["I"] = ".";
    }
    for (std::size_t n: ns)
    for (Chain* c :  _chains)
    for (Polypeptide& pp: c->_polypeptides)
    {
        int helices_to_fill = 0;
        int turns_to_fill = 0;

        for (std::size_t i = 1; i < pp._residues.size() ; ++i)
        {
            Residue* res = pp._residues[i];
            Residue* res_prev = pp._residues[i-1];
            if (turns_to_fill)
            {
                res->ssi.values[GHI[n]] = ghi[n];
                turns_to_fill--;
            }
            if (helices_to_fill)
            {
                res_prev->ssi.values[GHI[n]] = GHI[n];
                helices_to_fill--;
            }
            if (res->_turns_starting[n])
            {
                if (helices_to_fill || res_prev->_turns_starting[n])
                {
                    helices_to_fill = n;
                }
                turns_to_fill = n-1;
            }

            assert(helices_to_fill>=0);
            assert(turns_to_fill>=0);

        }
    }
}

void  Protein::collect_for_dssp_ladders()
{
    for (BetaLadder &ladder: _ladders)
    {
        //cerr << ladder << endl;

        std::deque<BetaBridge*>::iterator it_begin = ladder.mBridges.begin();
        std::deque<BetaBridge*>::iterator it_end = ladder.mBridges.end();
        std::deque<BetaBridge*>::value_type& first_bridge = ladder.mBridges.front();
        std::deque<BetaBridge*>::value_type& last_bridge = ladder.mBridges.back();
        if (ladder.mType == BetaLadderType::antiparallel)
        {
            int dist_last =  last_bridge->res_b->get_sequence_dist(last_bridge->res_r);
            int dist_first = first_bridge->res_b->get_sequence_dist(first_bridge->res_r) ;
            if (dist_last > dist_first)
            {
                cerr <<  "will fail:" <<  dist_last <<  "\t" << dist_first << endl;
                cerr <<  ladder<< endl;
                cerr <<  "last_bridge  b:" <<  *last_bridge->res_b <<   "\tlast_bridge r "  << *last_bridge->res_r << endl;
                cerr <<  "first_bridge b:" <<  *first_bridge->res_b <<  "\tfirst_bridge b"  << *first_bridge->res_r << endl;
            }
            assert (dist_last <= dist_first);
          //  cerr << "AP-dist: " << *last_bridge  <<  ":" <<   dist_last << endl;
           // cerr << "         " << *first_bridge <<  ":" <<  dist_first << endl;
            BetaBridge* bridge = last_bridge;
            unsigned int start_i;
            unsigned int end_i;
            assert(*bridge->res_b < *bridge->res_r);
            if (*bridge->res_b < *bridge->res_r)
            {
                start_i = bridge->res_b->_total_res_no;
                end_i   = bridge->res_r->_total_res_no;
            } else {
                end_i   = bridge->res_b->_total_res_no;
                start_i = bridge->res_r->_total_res_no;
            }
            assert(start_i < end_i);
            assert(end_i);
            if (dist_last <= 5)
            {
                assert(end_i-1 <_residues.size());
                for (unsigned int i = start_i; i < end_i-1; ++i)
                {
                    _residues[i]->ssi.values["beta-turn"] = 'T';
                }
            }
        }

        bool b5_clear = std::all_of(it_begin, it_end, [](BetaBridge*& bridge){return bridge->res_b->ssi._pseudo_dssp_full[5] == ' ';});
       // bool b6_clear = std::all_of(it_begin, it_end, [](BetaBridge*& bridge){return bridge->res_b->ssi._pseudo_dssp_full[6] == ' ';});
        //bool r5_clear = std::all_of(it_begin, it_end, [](BetaBridge*& bridge){return bridge->res_r->ssi._pseudo_dssp_full[5] == ' ';});
        bool r6_clear = std::all_of(it_begin, it_end, [](BetaBridge*& bridge){return bridge->res_r->ssi._pseudo_dssp_full[6] == ' ';});

        string marker = ladder.get_marker();
        int index_b = 5;
        int index_r = 6;
        if ( !b5_clear )
        {
            //assert(b6_clear); sadly this is not always true
            index_b = 6;
            index_r = 5;
        }

        if ( !r6_clear )
        {
            // assert(r5_clear); sadly this is not always true
            index_b = 5;
            index_r = 6;
        }
        assert (index_b != index_r);
        std::for_each( it_begin, it_end, [&](BetaBridge*& bridge) {
            if (bridge->res_b->ssi._pseudo_dssp_full[index_b] == ' ')
                bridge->res_b->ssi._pseudo_dssp_full[index_b] = marker[0] ;
            else
                bridge->res_b->ssi._pseudo_dssp_full[index_r] = marker[0] ;
            if (bridge->res_r->ssi._pseudo_dssp_full[index_r] == ' ')
                bridge->res_r->ssi._pseudo_dssp_full[index_r] = marker[0] ;
            else
                bridge->res_r->ssi._pseudo_dssp_full[index_b] = marker[0] ;
        });

        /*if (ladder.mBridges.front()->res_b->ssi._pseudo_dssp_full[index_b] == ' ')
            ladder.mBridges.front()->res_b->ssi._pseudo_dssp_full[index_b] = '[';
        if (ladder.mBridges.back()->res_b->ssi._pseudo_dssp_full[index_b] == ' ')
            ladder.mBridges.back()->res_b->ssi._pseudo_dssp_full[index_b] = ']';

        if (ladder.mBridges.front()->res_r->ssi._pseudo_dssp_full[index_r] == ' ')
            ladder.mBridges.front()->res_r->ssi._pseudo_dssp_full[index_r] = '^';
        if (ladder.mBridges.back()->res_r->ssi._pseudo_dssp_full[index_r] == ' ')
            ladder.mBridges.back()->res_r->ssi._pseudo_dssp_full[index_r] = 'v';*/
    }
}


void Protein::collect_for_dssp_bulges()
{
    cerr << "collect_for_dssp_bulges" << endl;
    for (BetaLadder &ladder: _ladders)
    {
        //cerr << ladder << endl;
        string ap_marker;
        if (ladder.mType == BetaLadderType::antiparallel)
            if (ladder.mBridges.size() > 1)
                ap_marker = "A";
            else
                ap_marker = "a";
        else
            if (ladder.mBridges.size() > 1)
                ap_marker = "P";
            else
                ap_marker = "p";                            // std::deque<BetaBridge*>::iterator it_begin = ladder.mBridges.begin();
       // std::deque<BetaBridge*>::iterator it_end = ladder.mBridges.end();
        //std::deque<BetaBridge*>::value_type& first_bridge = ladder.mBridges.front();
        //std::deque<BetaBridge*>::value_type& last_bridge = ladder.mBridges.back();
        //std::vector<Residue*>::iterator it_begin(ladder.mBridges.front()->res_b); //Polypeptide::_residues::iterator
        auto it_begin = ladder.mBridges.front()->res_b->get_iter();
        auto it_end   = ladder.mBridges.back()->res_b->get_iter();
        std::for_each(it_begin, it_end+1, [&ladder, &ap_marker]( Residue* it ){
            //it->ssi._pseudo_dssp_full[0] = 'b';
            if (it->ssi.values["bulge1"] == "")
                it->ssi.values["bulge1"] = ladder.get_marker();
            else
                it->ssi.values["bulge2"] = ladder.get_marker();

            it->ssi.values["AP"] += ap_marker;
            ladder.append_bulged_residue(&(*it));
        });

        auto it_beginr = ladder.mBridges.back()->res_r->get_iter();
        auto it_endr   = ladder.mBridges.front()->res_r->get_iter();
      //  cerr << **it_beginr <<  " to " <<  **it_endr <<  endl;
        if ((**it_beginr)._res_no > (**it_endr)._res_no)
            std::swap(it_beginr, it_endr);
        std::for_each(it_beginr, it_endr+1, [&ladder, &ap_marker]( Residue* it ){
            //cerr <<  *it <<  endl;
           // it->ssi._pseudo_dssp_full[1] = 'r';
            if (it->ssi.values["bulge1"] == "")
                it->ssi.values["bulge1"] = ladder.get_marker();
            else
                it->ssi.values["bulge2"] = ladder.get_marker();
            it->ssi.values["AP"] += ap_marker;
            ladder.append_bulged_residue(&(*it));
        });
        //cerr <<  endl;
        // std::vector<Residue*> _residues;
       // cerr <<  **it_begin <<  endl;

    }
}

void Protein::collect_for_pssc_strands()
{
// important: collect_for_dssp_bulges needs to be called first!
  //  cerr << "collect_for_pssc_strands" << endl;
/*  for (BetaLadder &ladder: _ladders)
    {
        for (Residue* res : ladder.mBulgeResidues)
        {

        }
    }*/
}

double to360(double x)
{
    if (x < 0)
        x += 360;

   /*
   if (x > 360)
       cerr <<  "assertion failed: " <<  x << endl;

    if (x < 0)
        cerr <<  "assertion failed: " <<  x << endl;
        */

    return x;
}

bool is_extended(double theta, double d)                    //  region F
{
    theta = to360(theta);
    if (d < 2.5)
    return false;

    return (149.1 < theta && theta < 231.9);
}

bool is_extended_strict(double theta, double d)             //  region E
{
    theta = to360(theta);
    if (d < 2.5)
        return false;
    return (158.2 < theta && theta < 221.3);
}

bool is_pii_helical(double theta, double d)                 //  region P
{
    theta = to360(theta);
    if (d < 2.5)
        return false;
    // if (-135 < theta and theta < -80)
    return (212.4 < theta && theta < 264.6);
}

bool is_helical(double theta, double d)                     //  region H
{
    theta = to360(theta);
    if (d < 0.5)
        return false;
    return (30 < theta && theta < 135);
}

bool is_small_helical(double theta, double d)               //  region D
{
    return (d < 1.0);
}

bool is_left_handed(double theta, double d)                 //  region L
{
    theta = to360(theta);

    if (d < .5 || d > 2.8)
        return false;
    return (210 < theta && theta < 300);
}

void Protein::collect_for_dihedrals()
{
    cerr <<  "collect_for_dihedrals" << endl;

    for (Chain* c :  _chains)
    for (Polypeptide& pp: c->_polypeptides)
    {
        for (std::size_t i = 1; i < pp._residues.size() - 4; ++i)
        {

            double theta, d;

            bool all_helical = true;
            bool all_extended = true;
            bool all_extended_strict = true;
            bool all_left_handed = true;
            bool all_pii_helical = true;
            bool all_small_helical = true;

            for (int32 j = -1; j <= 1; ++j)
            {
                theta = pp._residues[i+j]->dhi.theta;
                d = pp._residues[i+j]->dhi.d;
                all_extended        &= is_extended(theta, d);
                all_extended_strict &= is_extended_strict(theta, d);
                all_pii_helical     &= is_pii_helical(theta, d);
                all_helical         &= (is_helical(theta, d) or is_small_helical(theta, d));
                all_small_helical   &= is_small_helical(theta, d);
                all_left_handed     &= is_left_handed(theta, d);
            }

            for (int32 j = -1; j <= 1; ++j)
            {
                theta = pp._residues[i+j]->dhi.theta;
                d = pp._residues[i+j]->dhi.d;
                if (all_helical)
                {
                    pp._residues[i+j]->ssi.values["dh_H"] = "H";
                }
                else
                {
                    if (is_helical(theta, d))
                        if (pp._residues[i+j]->ssi.values["dh_H"] != "H")
                            pp._residues[i+j]->ssi.values["dh_H"] = "h";
                }

                if (all_extended)
                {
                    pp._residues[i+j]->ssi.values["dh_F"] = "F";
                    if (is_extended_strict(theta, d))
                        pp._residues[i+j]->ssi.values["dh_E"] = "E";
                }
                else
                {
                    if (is_extended(theta, d))
                        if (pp._residues[i+j]->ssi.values["dh_F"] != "F")
                            pp._residues[i+j]->ssi.values["dh_F"] = "f";
                    if (is_extended_strict(theta, d))
                        if (pp._residues[i+j]->ssi.values["dh_E"] != "E")
                            pp._residues[i+j]->ssi.values["dh_E"] = "e";

                }
                if (all_left_handed)
                {
                    pp._residues[i+j]->ssi.values["dh_L"] = "L";
                }
                else
                {
                    if (is_left_handed(theta, d))
                        if (pp._residues[i+j]->ssi.values["dh_L"] != "L" )
                            pp._residues[i+j]->ssi.values["dh_L"] = "l";
                }

                if (all_pii_helical)
                    pp._residues[i+j]->ssi.values["dh_P"] = "P";
                else
                {
                    if (is_pii_helical(theta, d))
                        if (pp._residues[i+j]->ssi.values["dh_P"] != "P")
                            pp._residues[i+j]->ssi.values["dh_P"] = "p";
                }

                if (all_small_helical)
                {
                    pp._residues[i+j]->ssi.values["dh_D"] = "D";
                }
                else
                {
                    if (is_small_helical(theta, d))
                        if (pp._residues[i+j]->ssi.values["dh_D"] != "D")
                            pp._residues[i+j]->ssi.values["dh_D"] = "d";
                }
            }
        }
    }
}

void Protein::collect_secondary_structure()
{
    std::cerr << "collect_secondary_structure " << _residues.size() << std::endl;

    for (Residue *res :  _residues)
    {
       res->ssi.values["S"]   = (res->dhi.kappa < 110 ?  "S" :  ".");
       res->ssi.values["+/-"] = (res->dhi.alpha < 0   ?  "-" :  "+");
    }

    collect_for_dssp_turns();
    collect_for_dssp_ladders();
    collect_for_dssp_bulges();
    collect_for_pssc_turns();
    collect_for_pssc_strands();
    collect_for_dihedrals();

    for (Residue *res :  _residues)
    {
        if (res->ssi.dssp_one_csi != "-")
            continue;
        if (res->ssi._pseudo_dssp_full[5] != '.' ||
            res->ssi._pseudo_dssp_full[6] != '.' ||
            res->ssi.values["bulge1"] != "" ||
            res->ssi.values["bulge2"] != ""
          )
            res->ssi.dssp_one_csi = "E";
        if (res->ssi.dssp_one_csi == "-" && res->ssi.values["S"] == "S")
            res->ssi.dssp_one_csi = "S";
    }
}

string Protein::get_fasta(int row_length, bool indicate_gaps) const
{
    /**
    >1F34:A|PDBID|CHAIN|SEQUENCE
    IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCSSLACSDHNQFNPDDSSTFEATSQELSITYGTGSM
    TGILGYDTVQVGGISDTNQIFGLSETEPGSFLYYAPFDGILGLAYPSISASGATPVFDNLWDQGLVSQDLFSVYLSSNDD
    SGSVVLLGGIDSSYYTGSLNWVPVSVEGYWQITLDSITMDGETIACSGGCQAIVDTGTSLLTGPTSAIANIQSDIGASEN
    SDGEMVISCSSIDSLPDIVFTINGVQYPLSPSAYILQDDDSCTSGFEGMDVPTSSGELWILGDVFIRQYYTVFDRANNKV
    GLAPVA
    >1F34:B|PDBID|CHAIN|SEQUENCE
    QFLFSMSTGPFICTVKDNQVFVANLPWTMLEGDDIQVGKEFAARVEDCTNVKHDMAPTCTKPPPFCGPQDMKMFNFVGCS
    VLGNKLFIDQKYVRDLTAKDHAEVQTFREKIAAFEEQQENQPPSSGMPHGAVPAGGLSPPPPPSFCTVQ
    **/
    string ret = "";
    for (const auto& chain :  _chains)
    {
        ret += ">" + mProteinInfo.mID + ":" + chain->_id + "|PDBID|CHAIN|SEQUENCE\n";
        ret += chain->get_fasta_raw(row_length, indicate_gaps) + "\n";
    }
    return ret;
}

void Protein::dump_ss_col()
{
    for (Chain* chain: _chains)
    {
        cout << "# Chain " <<  chain->_id << endl;
        for (Polypeptide pp: chain->_polypeptides)
        {
            for (Residue* res: pp._residues)
            {
                cout << res->get_ssi_output_col()<< endl;
            }
            cout << "# Ter" << endl;
        }
    }
    cerr <<  "## Total ASA: " << asa <<   endl;
}

void Protein::dump_ss_dssp()
{
    for (Chain* chain: _chains)
    {
        cout << "# Chain " <<  chain->_id << endl;
        for (Polypeptide pp: chain->_polypeptides)
        {
            for (Residue* res: pp._residues)
            {
                cout << res->get_ssi_output_dssp()<< endl;
            }
            cout << "# Ter" << endl;
        }
    }
}

void Protein::dump_ss_tab()
{
    for (Chain* chain: _chains)
    {
        cout << "# Chain " <<  chain->_id << endl;
        for (Polypeptide pp: chain->_polypeptides)
        {
            for (Residue* res: pp._residues)
            {
                cout << res->get_ssi_output_tab()<< endl;
            }
            cout << "# Ter" << endl;
        }
    }
}


void Protein::dump_ss_json()
{
    cout << "{\"version\":\"" << VERSION << "\"," << endl;
    string pdbinfo = "{\"title\":\"" + string("TITLE") + " \"}";

    cout << "\"pdbinfo\":" << pdbinfo << "," << endl;
    cout << "\"protein\":{\"chains\": [" << endl;

    uint32 chain_ind = 0;
    uint32 res_ind;
    uint32 pp_ind;

    for (Chain* chain: _chains)
    {
        if (chain_ind++)
            cout << "," << endl;

        cout << "{\"chainid\":\"" << chain->_id << "\",\"residues\":[" << endl;

        res_ind = 0;
        pp_ind = 0;
        for (Polypeptide pp: chain->_polypeptides)
        {
            if (pp_ind++) {
                if (res_ind)
                    cout <<  ",";
                cout << "{\"resid\":\"GAP\"}";
            }
            for (Residue* res: pp._residues)
            {
                if (res_ind++)
                    cout << ",";

                cout << res->get_ssi_output_json();
            }
        }
        cout << "]";          //  residues
        cout << "}" << endl;  //  chain
    } // chains
    cout << (boost::format("], \"acc\": %.0f}") % asa).str() << endl;


    cout << "}" << endl;                                       //  protein
}


void Protein::dump_ladders()
{
    for (BetaLadder &ladder1 : _ladders)
    for (BetaLadder &ladder2 : _ladders)
    {
        if (ladder1._no == ladder2._no)
            continue;

        if (ladder1.mType == BetaLadderType::antiparallel and
            ladder2.mType == BetaLadderType::antiparallel)
        {
            Residue& res_b1 = *ladder1.get_first_bridge()->res_b;
            Residue& res_b2 = *ladder2.get_first_bridge()->res_b;
            Residue& res_r1 = *ladder1.get_first_bridge()->res_r;
            Residue& res_r2 = *ladder2.get_first_bridge()->res_r;
            int diff_b = res_b1.get_total_dist_to(res_b2);
            int diff_r = res_r2.get_total_dist_to(res_r1);
            if (abs(diff_b) > 8) continue;
            if (abs(diff_r) > 8) continue;
            std::pair<int, int> bounds = std::minmax(diff_b, diff_r);
            int max_len_smaller_gap = 2;
            int max_len_larger_gap = 5;
            if (ladder1.get_first_bridge()->mType == BetaBridgeType::AP_short)
            {
                max_len_smaller_gap += 1;
                max_len_larger_gap += 1;
            }
            if (ladder2.get_first_bridge()->mType == BetaBridgeType::AP_short)
            {
                max_len_smaller_gap += 1;
                max_len_larger_gap += 1;
            }
            if (bounds.first > 0 and bounds.first <= max_len_smaller_gap and bounds.second > 0 and bounds.second <= max_len_larger_gap)
            {
                cerr <<  "bulge" << "\n";
                cerr << res_b1 << " to " << res_b2 <<  endl;
                cerr << res_r1 << " to " << res_r2 <<  endl;
            }
            else
                cerr <<  "     " << "\t";
            cerr << ladder1 << "\tvs\t" << ladder2 << "\t"  ;
            cerr << diff_b << "\t" << diff_r << endl;
        }

    }
}

void Protein::dump_hbonds()
{
    for (Chain* chain: _chains)
    {
        cout << "Chain " <<  chain->_id << endl;
        for (Polypeptide pp: chain->_polypeptides)
        {
            for (Residue* res: pp._residues)
            {
                cout << *res << "\t";
                int i = 0;
                for (HBond* hb: res->get_sorted_HBond_donating())
                {
                    cout << *hb;
                    cout << "\n(d_HO=" << hb->mDonor["H"].mLoc.getDistanceFrom(hb->mAcceptor["O"].mLoc)<< "\n";
                    cout << " d_HC=" << hb->mDonor["H"].mLoc.getDistanceFrom(hb->mAcceptor["C"].mLoc) << "\n";
                    cout << " d_NO=" << hb->mDonor["N"].mLoc.getDistanceFrom(hb->mAcceptor["O"].mLoc)<< "\n";
                    cout << " d_NC=" << hb->mDonor["N"].mLoc.getDistanceFrom(hb->mAcceptor["C"].mLoc);
                    cout << "\nH =" << hb->mDonor["H"].mLoc.X <<  ", "
                    <<  hb->mDonor["H"].mLoc.Y<<  ", "
                    <<  hb->mDonor["H"].mLoc.Z<<  "\n";

                    cout <<  ")" ;
                    if (++i >= 2) break;
                }
                //cout << endl;
                cout << " | " ;
                i = 0;
                for (HBond* hb: res->get_sorted_HBond_accepting())
                {
                    cout << *hb ;
                    if (++i >= 2) break;
                }
                cout << endl;
            }
            cout << "Ter\n" << endl;
        }
    }
}


void Protein::dump_hbonds2()
{
   for (HBond hb :  mHbonds)
   {
       cerr << hb << "\t" << hb.ca_dist <<  endl;
   }
}
