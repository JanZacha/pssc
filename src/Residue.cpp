
#include "NamedFormat.h"
#include "Residue.h"
#include "HydrogenBond.h"
#include "Polypeptide.h"

std::vector<Residue*>::iterator Residue::get_iter()
{
    return _polypeptide->_residues.begin() + (_res_no - 1);
}

bool sort_HBond(HBond* a,  HBond*b)
{
    return *a < *b;
}

std::deque<HBond*> Residue::get_sorted_HBond_donating()
{
    auto ret (HBond_donating);
    std::sort (ret.begin(), ret.end(), sort_HBond);
    return ret;
}

std::deque<HBond*> Residue::get_sorted_HBond_accepting()
{
    auto ret (HBond_accepting);
    std::sort (ret.begin(), ret.end(), sort_HBond);
    return ret;
}


string Residue::get_ssi_output_col()
{
    NamedFormat fmt("%resno$5.5d%seqno$5.5d%icode$1.1s%chain$1.1s %AA$c %1CSSI_dssp$c  |%7CSSI_dssp$7s| [%bulge1$1s%bulge2$1s] |%pssc$s| %ASA$3i %GARBAGE$12s "
        "| %tco$6.3f %kappa$6.1f %alpha$6.1f | %phi$6.1f %psi$6.1f |%d$4.1f  %theta$6.1f %r$4.1f | %CAx$6.1f %CAy$6.1f %CAz$6.1f");

    fmt ("resno",  _total_res_no)
        ("seqno",  ID.seqNumber)
        ("icode",  ID.insertionCode)
        ("chain",  ID.chain)
        ("AA",     _resinfo->mSingleLetter)
        ("1CSSI_dssp",  ssi.dssp_one_csi)
        ("7CSSI_dssp",  ssi.get_pseudo_dssp_full())
        ("bulge1",  ssi.values["bulge1"])
        ("bulge2",  ssi.values["bulge2"])
      //  ("beta-turn", ssi.values["beta-turn"])
        ("pssc",  ssi.get_pssc_full())
        ("ASA",  floor(asa + 0.5))
        ("GARBAGE",  "...")
        ("tco",    dhi.tco)
        ("kappa",  dhi.kappa)
        ("alpha",  dhi.alpha)
        ("phi",    dhi.phi)
        ("psi",  dhi.psi)
        ("d",    dhi.d)
        ("theta",  dhi.theta)
        ("r",  dhi.r)
        ("CAx",  get("CA").mLoc.X)
        ("CAy",  get("CA").mLoc.Y)
        ("CAz",  get("CA").mLoc.Z) ;

    return fmt.str();
}

string Residue::get_ssi_output_dssp()
{
    //"%5.5d%5.5d%1.1s%1.1s %c  %c %c%c%c%c%c%c%c%4.4d%4.4d%c%4.4d %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f");
    NamedFormat fmt("%resno$5.5d%seqno$5.5d%icode$1.1s%chain$1.1s %AA$c %1CSSI_dssp$c  |%7CSSI_dssp$7s| [%bulge1$1s%bulge2$1s] |%pssc$s| %ASA$3i %GARBAGE$12s "
        "| %tco$6.3f %kappa$6.1f %alpha$6.1f | %phi$6.1f %psi$6.1f |%d$4.1f  %theta$6.1f %r$4.1f | %CAx$6.1f %CAy$6.1f %CAz$6.1f");

    fmt ("resno",  _total_res_no)
    ("seqno",  ID.seqNumber)
    ("icode",  ID.insertionCode)
    ("chain",  ID.chain)
    ("AA",     _resinfo->mSingleLetter)
    ("1CSSI_dssp",  ssi.dssp_one_csi)
    ("7CSSI_dssp",  ssi.get_pseudo_dssp_full())
    ("bulge1",  ssi.values["bulge1"])
    ("bulge2",  ssi.values["bulge2"])
    //  ("beta-turn", ssi.values["beta-turn"])
    ("pssc",  ssi.get_pssc_full())
    ("ASA",  floor(asa + 0.5))
    ("GARBAGE",  "...")
    ("tco",    dhi.tco)
    ("kappa",  dhi.kappa)
    ("alpha",  dhi.alpha)
    ("phi",    dhi.phi)
    ("psi",  dhi.psi)
    ("d",    dhi.d)
    ("theta",  dhi.theta)
    ("r",  dhi.r)
    ("CAx",  get("CA").mLoc.X)
    ("CAy",  get("CA").mLoc.Y)
    ("CAz",  get("CA").mLoc.Z) ;
    /*%
    residue.GetNumber() % ca.mResSeq % ca.mICode % ca.mChainID % ssi.GetCode() %
    ssi.GetSS() % helix[0] % helix[1] % helix[2] % ssi.GetBend() % ssi.GetChirality() % bridgelabel[0] % bridgelabel[1] %
    bp[0] % bp[1] % sheet % floor(residue.Accessibility() + 0.5) %
    NHO[0] % ONH[0] % NHO[1] % ONH[1] %
    residue.TCO() % residue.Kappa() % ssi.GetAlpha() % residue.Phi() % residue.Psi() %
    ca.mLoc.mX % ca.mLoc.mY % ca.mLoc.mZ).str();*/
    return fmt.str();
}

string Residue::get_ssi_output_tab()
{
    //"%5.5d%5.5d%1.1s%1.1s %c  %c %c%c%c%c%c%c%c%4.4d%4.4d%c%4.4d %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f");
    NamedFormat fmt("%resno$5.5d%seqno$5.5d%icode$1.1s%chain$1.1s %AA$c %1CSSI_dssp$c  |%7CSSI_dssp$7s| [%bulge1$1s%bulge2$1s] < %1CSSI_pssc$1s > |%pssc$s| %ASA$3i %GARBAGE$12s "
        " %phi$6.1f %psi$6.1f |%d$4.1f  %theta$6.1f %r$4.1f ");

    fmt ("resno",  _total_res_no)
    ("seqno",  ID.seqNumber)
    ("icode",  ID.insertionCode)
    ("chain",  ID.chain)
    ("AA",     _resinfo->mSingleLetter)
    ("1CSSI_dssp",  ssi.dssp_one_csi)
    ("7CSSI_dssp",  ssi.get_pseudo_dssp_full())
    ("bulge1",  ssi.values["bulge1"])
    ("bulge2",  ssi.values["bulge2"])
    //  ("beta-turn", ssi.values["beta-turn"])
    ("1CSSI_pssc",  ssi.get_pssc_one_csi())
    ("pssc",  ssi.get_pssc_full())
    ("ASA",  floor(asa + 0.5))
    ("GARBAGE",  "...")
    ("phi",    dhi.phi)
    ("psi",  dhi.psi)
    ("d",    dhi.d)
    ("theta",  dhi.theta)
    ("r",  dhi.r)
 ;
    return fmt.str();
}

string Residue::get_ssi_output_json()
{
    NamedFormat fmt_resid("%AA3$s%seqno$d%icode$1.1s");

    fmt_resid("AA3",  _resinfo->mThreeLetter)
             ("seqno",  ID.seqNumber)
             ("icode",  ID.insertionCode);
    string resid = fmt_resid.str();
    ba::trim(resid);
    NamedFormat fmt_resline("{"
            "\"resid\":\"%resid$s\","
            "\"AA\":\"%AA$s\","
            "\"pssc1\":\"%1CSSI_dssp$1.1s\","
            "\"pssc\":\"%pssc$s\","
            "\"phi\":%phi$.1f,"
            "\"psi\":%psi$.1f,"
            "\"d\":%d$.2f,"
            "\"theta\":%theta$.1f,"
            "\"r\":%r$.2f,"
            "\"ASA\":%ASA$.1f,"
            "\"tco\":%tco$.1f,"
            "\"kappa\":%kappa$.1f"
        "}");
    fmt_resline ("resid", resid)
                ("AA",     _resinfo->mSingleLetter)
                ("1CSSI_dssp",  ssi.get_pssc_one_csi())
                ("7CSSI_dssp",  ssi.get_pseudo_dssp_full())
                ("pssc",  ssi.get_pssc_full())
                ("ASA",  floor(asa + 0.5))
                ("tco",    dhi.tco)
                ("kappa",  dhi.kappa)
                ("alpha",  dhi.alpha)
                ("phi",    dhi.phi)
                ("psi",  dhi.psi)
                ("d",    dhi.d)
                ("theta",  dhi.theta)
                ("r",  dhi.r)
    ;
    fmt_resline.setDefaultForNonFinite(-999);
    return fmt_resline.str();
    /*


        const MAtom& ca = residue.GetCAlpha();

        MSecondaryStructureInformation ssi = MSecondaryStructureInformation(residue);
        const MPoint helical_center = residue.GetHelicalCenter();
        //const MAtom& ca = residue.GetCAlpha();
        //const MAtom& h = residue.GetH();
        //const MPoint helical_center = residue.GetHelicalCenter();
        const uint32 acc = floor(residue.Accessibility() + 0.5);
        string helix = ssi.GetHelix();
        //char code = kResidueInfo[residue.GetType()].code;
        string name = string(kResidueInfo[residue.GetType()].name);

        //if (residue.GetType() == kCysteine and residue.GetSSBridgeNr() != 0)
        //      code = 'a' + ((residue.GetSSBridgeNr() - 1) % 26);

        string centerstr = (boost::format("\"center\":[%.2f, %.2f, %.2f]") % helical_center.mX % helical_center.mY % helical_center.mZ).str();

        double theta, d, r;
        tr1::tie(theta, d, r) = residue.ThetaD();
        string icode = ca.mICode;
        if (icode == " ")
        icode = "";
        return (kDSSPResidueLine % name % ca.mResSeq % icode % residue.GetPSSC1() %
            residue.GetPSSCode().GetCode() % residue.Phi() % residue.Psi() % d % theta % r % acc % centerstr).str();

 */
}