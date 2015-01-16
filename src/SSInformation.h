#pragma once
#include "standardlibs.h"


struct SSInformation {

    //string single_letter = "?";
    string dssp_one_csi = "-";
    string helix = "___";
    string _pseudo_dssp_full = "___.....";
    string get_pseudo_dssp_full()
    {
        string ap_marker = values["AP"];
        if (ap_marker.size() == 1)
            ap_marker += " ";
        if (ap_marker.size() == 0)
            ap_marker += "  ";
        return _pseudo_dssp_full + values["S"] + values["+/-"];// + ap_marker + values["G"] + values["H"] + values["I"];
    }

    string get_pssc_full()
    {
        string ap_marker = values["AP"];
        if (ap_marker.size() == 1)
            ap_marker += ".";
        if (ap_marker.size() == 0)
            ap_marker += "..";
        return values["G"] + values["H"] + values["I"] + ap_marker + values["S"] + values["+/-"] + values["beta-turn"] + values["dtheta"]
        + values["dh_E"]+ values["dh_F"]+ values["dh_P"]+ values["dh_H"]+ values["dh_D"]+ values["dh_L"];
    }


    string get_pssc_one_csi()
    {
        if (values["dh_L"] == "L") {
                if (values["G"] == "G") return "L";
                if (values["H"] == "H") return "L";
                if (values["I"] == "I") return "L";
        }

        if (values["I"] == "I") return "I";
        if (values["H"] == "H") return "H";
        if (values["G"] == "G") return "G";

        for (auto c: values["AP"])
            if (c == 'A' || c == 'P')
                return "E";

        if (values["dh_E"] == "E")
        {
            for (auto c: values["AP"])
                if (c == 'a' || c == 'p')
                    return "E";
        }

        //if (values["beta-turn"] == "T")
        //    return "U";

        if (values["I"] == "i") return "T";
        if (values["H"] == "h") return "T";
        if (values["G"] == "g") return "T";

        if (values["dh_P"] == "P") return "P";
        if (values["dh_F"] == "F") return "F";

        if (values["S"] == "S") return "S";
        return "C";
    }

    std::unordered_map<string, string> values = {
        {"S",  "."} ,
        {"+/-",  ""} ,
        {"AP",  ""} ,
        {"G",  ""} ,
        {"H",  ""} ,
        {"I",  ""} ,
        {"bulge1",  ""} ,
        {"bulge2",  ""} ,
        {"bridge1",  ""} ,
        {"bridge2",  ""} ,
        {"beta-turn",  "."},
        {"dtheta",  ""},

        {"dh_E",  " "},             // dihedral_pssc_extended_strict
        {"dh_F",  " "},             // dihedral_pssc_extended
        {"dh_P",  " "},             // dihedral_pssc_pii
        {"dh_H",  " "},                                      // dihedral_pssc_helix
        {"dh_D",  " "},                                      // dihedral_pssc_helix_small
        {"dh_L",  " "}                                       // dihedral_pssc_left

    };
};

