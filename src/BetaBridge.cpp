#include "BetaBridge.h"

BetaBridge BetaBridge::create_from_residues(
        Residue* res_a, Residue* res_b, Residue* res_c,
        Residue* res_q, Residue* res_r, Residue* res_s
      )
{
    HBond* hBondA = nullptr;
    HBond* hBondB = nullptr;

    /*if (res_b->has_betabridge_assignend)
    {
        cerr <<  "SWITCHING due to "<< *res_b <<  endl;
        Residue* res_1; Residue* res_2; Residue* res_3;
        res_1 = res_q; res_2 = res_r; res_3 = res_s;
        res_q = res_a; res_r = res_b; res_a = res_c;
        res_a = res_1; res_b = res_2; res_s = res_3;
    }*/
    int is_AP_short, is_AP_long, is_P_forward, is_P_backward;

    #define ARGS res_a, res_b, res_c, res_q, res_r, res_s

    if ((is_AP_short   = test_AP_short_cycle(ARGS, hBondA, hBondB)))
    {
        assert(hBondB);
        assert(hBondA);
    }
    if ((is_AP_long    = test_AP_long_cycle(ARGS, hBondA, hBondB)))
    {
        assert(hBondB);
        assert(hBondA);
    }
    if ((is_P_forward  = test_P_forward_cycle(ARGS, hBondA, hBondB)))
    {
        assert(hBondB);
        assert(hBondA);
    }
    if ((is_P_backward = test_P_backward_cycle(ARGS, hBondA, hBondB)))
    {
        assert(hBondB);
        assert(hBondA);
    }


    int sum = (is_AP_short+is_AP_long+is_P_forward+is_P_backward);
    if (sum > 1) {
        std::cerr << "Warning: more than one strand cycle possible at " << *res_b << "-" << *res_r << std::endl;
        std::cerr << "Will consider this to be AP" << std::endl;
    }
    BetaBridgeType type = BetaBridgeType::none;

    if (is_P_forward)
        type = BetaBridgeType::P_forward;
    if (is_P_backward)
        type = BetaBridgeType::P_backward;

    if (is_AP_short)
        type = BetaBridgeType::AP_short;
    if (is_AP_long)
        type = BetaBridgeType::AP_long;

    BetaBridge bridge(type);
    if (*res_b >= *res_r)
    {
        std::cerr << "Assertion error " << *res_b << " ... " << *res_r << std::endl;
    }
    assert(*res_b < *res_r);
    bridge.mHBondA = hBondA;
    bridge.mHBondB = hBondB;
    bridge.res_b = res_b;
    bridge.res_r = res_r;

    return bridge;
}

bool BetaBridge::test_AP_long_cycle (
        const Residue* res_a, const Residue* res_b, const Residue* res_c,
        const Residue* res_q, const Residue* res_r, const Residue* res_s,
        HBond*& outHBondA, HBond*& outHBondB)
    {
        /**
        * ..=(a)==(b)==(c)=>.
        *     O         H
        *     |         |
        *     H         O
        * .<=(s)==(r)==(q)=..
        */
            assert(res_a);
            assert(res_c);
            assert(res_q);
            assert(res_s);
            outHBondA = outHBondB = nullptr;
            return (HBond_from_to(res_s, res_a, outHBondA) and HBond_from_to(res_c, res_q, outHBondB));
    };

bool BetaBridge::test_AP_short_cycle(
        const Residue* res_b, const Residue* res_c, const Residue* res_d,
        const Residue* res_q, const Residue* res_r, const Residue* res_s,
        HBond*& outHBondA, HBond*& outHBondB)
    {
        /**
        * ..=(b)==(c)==(d)=>.
        *         H O
        *         | |
        *         O H
        * .<=(s)==(r)==(q)=..
        *
        */
        assert(res_c);
        assert(res_r);
        outHBondA = outHBondB = nullptr;
        return (HBond_from_to(res_c, res_r, outHBondA) && HBond_from_to(res_r, res_c, outHBondB));
    };

bool BetaBridge::test_P_backward_cycle(
        const Residue* res_a, const Residue* res_b, const Residue* res_c,
        const Residue* res_q, const Residue* res_r, const Residue* res_s,
        HBond*& outHBondA, HBond*& outHBondB)
    {
        /**
        * ..=(a)==(b)==(c)=>..
        *     O_       _H
        *       \_   _/
        *         H O
        * ..=(q)==(r)==(s)=>..
        */
        assert(res_b);
        assert(res_q);
        assert(res_s);
        outHBondA = outHBondB = nullptr;
        return (HBond_from_to(res_r, res_a, outHBondA) and HBond_from_to(res_c, res_r, outHBondB));
    };

bool BetaBridge::test_P_forward_cycle(
        const Residue* res_b, const Residue* res_c, const Residue* res_d,
        const Residue* res_q, const Residue* res_r, const Residue* res_s,
        HBond*& outHBondA, HBond*& outHBondB)
    {
        /**
        * ..=(b)==(c)==(d)=>..
        *        _H O_
        *      _/     \_
        *     O         H
        * ..=(q)==(r)==(s)=>..
        */
        assert(res_c);
        assert(res_q);
        assert(res_s);
        outHBondA = outHBondB = nullptr;
        return (HBond_from_to(res_c, res_q, outHBondA) and HBond_from_to(res_s, res_c, outHBondB));
    };

std::ostream& operator<<  ( std::ostream &output, const BetaBridge &self )
{
    string typestring = "? ";
    if (self.mType == BetaBridgeType::AP_long)
        typestring = "AP_long";
    else if (self.mType == BetaBridgeType::AP_short)
        typestring = "AP_short";
    else if (self.mType == BetaBridgeType::P_backward)
        typestring = "P_backward";
    else if (self.mType == BetaBridgeType::P_forward)
        typestring = "P_forward";
    output << "Bridge " << typestring<< " "<< *self.res_b <<  " : " <<  *self.res_r;
    return output;
}
