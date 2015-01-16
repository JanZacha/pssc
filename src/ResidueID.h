#pragma once

#include "standardlibs.h"
#include "tools_string.h"

struct ResidueID
{
    string resname;
    string chain;
    int seqNumber;
    string insertionCode;

    string shortstr()
    {
        return chain + std::to_string(seqNumber) + insertionCode;
    }

    ResidueID():
        resname("empty"),
        chain(""),
        seqNumber(999),
        insertionCode("")
    {};
    ResidueID(string iResname, string iChain, int iSeqNumber, string iInsertionCode):
        resname(iResname),
        chain(iChain),
        seqNumber(iSeqNumber),
        insertionCode(iInsertionCode)
    {
    }

    bool operator<(const ResidueID& o) const
    {
        return
        chain < o.chain or
        (chain == o.chain and seqNumber < o.seqNumber) or
        (chain == o.chain and seqNumber == o.seqNumber and insertionCode < o.insertionCode);
    }

    /**
     * Comparison operator does not take resname into accout.
     */
    bool operator!=(const ResidueID& o) const
    {
        return chain != o.chain or seqNumber != o.seqNumber or insertionCode != o.insertionCode;
    }

    static ResidueID from_PDB_line(string iLine)
    {
        //      18 - 20 Residue name resName Residue name.
        string resname = ba::trim_copy(iLine.substr(17, 4));
        //      22              Character chainID Chain identifier.
        string chain{iLine[21]};
        //      23 - 26 Integer resSeq Residue sequence number.
        uint16 seqNumber = boost::lexical_cast<int>(ba::trim_copy(iLine.substr(22, 4)));
        //      27              AChar iCode Code for insertion of residues.
        string insertionCode = iLine.substr(26, 1);
        trim(insertionCode);
        return ResidueID{resname,  chain, seqNumber, insertionCode};
    }

private:
    friend ostream &operator<<(ostream &output, const ResidueID &id)
    {
        output << id.chain << ":" << id.resname << id.seqNumber<< id.insertionCode;
        return output;
    }

};