#include "BetaLadder.h"


bool BetaLadder::can_be_extended(BetaBridge* bridge)
{
    assert(mBridges.size());
    const BetaBridge* first = mBridges.front();
    const BetaBridge* last  = mBridges.back();
    bool seq_dist_ok = false;
    int dista, distb,  distc,  distd;
    dista = first->res_b->get_sequence_dist(bridge->res_b) ;
    distb = first->res_r->get_sequence_dist(bridge->res_r);
    distc = last->res_b->get_sequence_dist(bridge->res_b) ;
    distd = last->res_r->get_sequence_dist(bridge->res_r) ;
    //cerr <<  dista <<  " " <<  distb <<  " " <<  distc <<  " " <<  distd <<  endl;
    if (dista <= 3)
    if (distb <= 3) seq_dist_ok = true;

    if (distc <= 3)
    if (distd <= 3) seq_dist_ok = true;

    if (!seq_dist_ok) return false;
    if (mType == BetaLadderType::parallel)
        return (bridge->mType == BetaBridgeType::P_backward || bridge->mType == BetaBridgeType::P_forward);
    return (bridge->mType == BetaBridgeType::AP_long || bridge->mType == BetaBridgeType::AP_short);
}