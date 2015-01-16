#pragma once

#include "standardlibs.h"
#include "BetaBridge.h"

enum class BetaLadderType : int {parallel, antiparallel};

struct BetaLadder
{
    BetaLadderType mType;
    deque<BetaBridge*> mBridges = {};
    std::vector<Residue*> mBulgeResidues = {};

    int _no;
    BetaLadder(BetaLadderType type):
        mType(type)
    {
    }

    BetaLadder(BetaBridge *first_bridge, int no):
    _no(no)
    {
        if (first_bridge->mType == BetaBridgeType::AP_long || first_bridge->mType == BetaBridgeType::AP_short)
            mType = BetaLadderType::antiparallel;
        else
            mType = BetaLadderType::parallel;
        //cerr <<  "New BetaLadder first_bridge=" << *first_bridge <<  endl;
        mBridges.push_back(first_bridge);
    }

    bool can_be_extended(BetaBridge* bridge);

    void append_bridge(BetaBridge* bridge)
    {
        //cerr << "append_bridge "<< *bridge << endl;
        mBridges.push_back(bridge);
    }

    void append_bulged_residue(Residue* res)
    {
        mBulgeResidues.push_back(res);
    }

    BetaBridge* get_last_bridge()
    {
        assert(mBridges.size());
        BetaBridge* a = mBridges.front();
        BetaBridge* b = mBridges.back();
        if (a->res_b->get_total_dist_to(*b->res_b) < 0)
            return a;
        return b;
    }

    BetaBridge* get_first_bridge()
    {
        assert(mBridges.size());
        BetaBridge* a = mBridges.front();
        BetaBridge* b = mBridges.back();
        if (a->res_b->get_total_dist_to(*b->res_b) < 0)
            return b;
        return a;
    }

    string get_marker()
    {
        string ret = " ";
        if (mType == BetaLadderType::antiparallel)
            ret[0] = 'A' + (_no - 1) % 26;
        else
            ret[0] = 'a' + (_no - 1) % 26;
        return ret;
    }

private:
    friend ostream &operator<<  ( ostream &output, const BetaLadder &self )
    {
        string typestring = "?";
        if (self.mType == BetaLadderType::parallel)
            typestring = "    parallel";
        else if (self.mType == BetaLadderType::antiparallel)
            typestring = "antiparallel";
        output << "Ladder " << self._no << " (" << self.mBridges.size() << ") " << typestring;
        //if (self.mBridges.size())
        //    output << " " <<  *self.mBridges.front() <<  " to " <<  *self.mBridges.back();
        return output;
    }
};

