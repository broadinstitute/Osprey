
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <set>
#include <string>

#include "omp.h"

#include "OspreyParams.hpp"
#include "PhaseImpMissing.hpp"
#include "Variant.hpp"
#include "VCFReader.hpp"
#include "VCFWriter.hpp"
#include "format.hpp"
#include "timestamp.hpp"

using namespace std;
using namespace Osprey;

int main(int argc, const char* argv[]) {

    OspreyParams params;
    int status = params.processCommandLineArgs(argc, argv);
    if (status == 0) {
        exit(0);
    }
    if (status < 0) {
        // cerr << "Aborting due to error processing command line arguments" << endl;
        exit(-status);
    }

    if (params.iterations <= 20 || params.iterations > 1000) {
        cerr << "Invalid number of iterations: " << params.iterations << endl;
        exit(1);
    }
    if (params.threads > 64) {
        cerr << "Invalid number of threads: " << params.threads << endl;
        exit(1);
    }

    if (params.verbose > 0) {
        cout << timestamp() << " Osprey version " << OspreyParams::getVersionString() << endl;
        cout << timestamp() << " Number of threads: " << params.threads << endl;
        cout << timestamp() << " Number of iterations of simulated annealing: " << params.iterations << endl;
    }

    const int seed = 1;
    srand(seed);
    omp_set_num_threads(params.threads);

    set<string> siteSet(params.siteList.begin(), params.siteList.end());

    PhaseImpMissing algorithm;
    algorithm.setDebug(params.debug);
    algorithm.setVerbose(params.verbose);
    algorithm.setIterations(params.iterations);

    if (params.verbose > 0) {
        cout << timestamp() << " Reading input file " << params.ibsMatrixFile << " ..." << endl;
    }
    algorithm.readIBSMatrix(params.ibsMatrixFile);

    if (params.verbose > 0) {
        cout << timestamp() << " Reading input file " << params.inputFile << " ..." << endl;
    }
    VCFReader vcfReader(params.inputFile);
    VCFWriter vcfWriter(params.outputFile);
    algorithm.updateVCFHeader(vcfReader, vcfWriter);
    while (true) {
        Variant* variant = vcfReader.nextVariant();
        if (variant == NULL) {
            break;
        }
        if (!siteSet.empty() && siteSet.find(variant->getId()) == siteSet.end()) {
            continue;
        }
        if (params.verbose > 1) {
            cout << timestamp() << " Processing site " << variant->getId() << " ..." << endl;
        }
        srand(seed);
        variant = algorithm.processVariant(variant);
        vcfWriter.writeVariant(variant);
    }
    vcfReader.close();
    vcfWriter.close();
    if (params.verbose > 0) {
        cout << timestamp() << " Run complete." << endl;
    }

    return(0);
}

