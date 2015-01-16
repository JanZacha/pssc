#pragma once

#include "standardlibs.h"
#include "ResidueID.h"
#include "HydrogenBond.h"
#include "Residue.h"

struct Atom;
struct Residue;
typedef std::pair<Atom, ResidueID> AtomResidPair;

enum class BetaBridgeType : int {none=0, AP_short, AP_long, P_forward, P_backward};

struct BetaBridge
{
    BetaBridgeType mType;
    HBond *mHBondA = nullptr;
    HBond *mHBondB = nullptr;

    Residue* res_b;
    Residue* res_r;

    BetaBridge( BetaBridgeType type = BetaBridgeType::none):
        mType(type)
    {
    }

    static BetaBridge create_from_residues(
        Residue* res_a, Residue* res_b, Residue* res_c,
        Residue* res_q, Residue* res_r, Residue* res_s
      );

    static bool test_AP_long_cycle (
        const Residue* res_a, const Residue* res_b, const Residue* res_c,
        const Residue* res_s, const Residue* res_t, const Residue* res_u,
        HBond*& outHBondA, HBond*& outHBondB);

    static bool test_AP_short_cycle(
        const Residue* res_b, const Residue* res_c, const Residue* res_d,
        const Residue* res_r, const Residue* res_s, const Residue* res_t,
        HBond*& outHBondA, HBond*& outHBondB);

    static bool test_P_backward_cycle(
        const Residue* res_a, const Residue* res_b, const Residue* res_c,
        const Residue* res_q, const Residue* res_r, const Residue* res_s,
        HBond*& outHBondA, HBond*& outHBondB);

    static bool test_P_forward_cycle(
        const Residue* res_b, const Residue* res_c, const Residue* res_d,
        const Residue* res_r, const Residue* res_s, const Residue* res_t,
        HBond*& outHBondA, HBond*& outHBondB);

private:
   friend std::ostream& operator<<  ( std::ostream &output, const BetaBridge &self );
};
