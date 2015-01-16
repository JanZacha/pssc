#include "standardlibs.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <stdexcept>

#include <irrlicht/irrlicht.h>
#include <irrlicht/quaternion.h>
#include <irrlicht/vector3d.h>

#include "Atom.h"
#include "Residue.h"
#include "Chain.h"
#include "Protein.h"
#include "pssc_config.h"
#include "ASACalculator.h"
using std::cerr;
using std::endl;

int main(int argc, char* argv[])
{

    PSSC_config config(argc, argv) ;

    if (config._exit)
    {
        cerr << config._what;
        return config._status;
    }
    cerr << config;

    ifstream infile(config.input.c_str(), ios_base::in | ios_base::binary);
    io::filtering_stream<io::input> in;

    if (ba::ends_with(config.input, ".bz2"))
    {
            in.push(io::bzip2_decompressor());
    }
    else if (ba::ends_with(config.input, ".gz"))
    {
            in.push(io::gzip_decompressor());
    }
    in.push(infile);
    Protein my_prot = Protein::create_from_pdb_file(in, false);
    my_prot.calculate_secondary_structure();
    //cout << my_prot << endl;
    //cout <<  "#" << endl;
    //cout << my_prot.get_fasta(100, true) << endl;


    ASACalculator c = ASACalculator(my_prot, 1.4, 200);
    cerr <<  "c.calculate();\n";
    c.calculate();
    cerr <<  "done with asa;\n";
    switch (config.output_format)
    {
        case  pssc_output_format::COL:
            my_prot.dump_ss_col();
            break;
        case  pssc_output_format::TAB:
            my_prot.dump_ss_tab();
            break;
        case  pssc_output_format::DSSP:
            my_prot.dump_ss_dssp();
            break;
        case  pssc_output_format::JSON:
            my_prot.dump_ss_json();
            break;
    }

    return 0;
}
