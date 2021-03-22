
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <algorithm>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include "omp.h"
#include "timestamp.hpp"

using namespace std;

#include "FileUtils.cpp"
#include "StringUtils.cpp"
#include "NumericUtils.cpp"


inline double sq(double x) { return x*x; }

double r2(vector<double> x, vector<double> y) {
    int N = x.size();
    double mu = NumericUtils::mean(x);
    for (int i = 0; i < N; i++) {
        x[i] -= mu;
    }
    mu = NumericUtils::mean(y);
    for (int i = 0; i < N; i++) {
        y[i] -= mu;
    }
    return sq(NumericUtils::dot(&x[0], &y[0], N)) / (NumericUtils::norm2(&x[0], N) * NumericUtils::norm2(&y[0], N));
}

void computeProbs(vector<double> &pVec, int jMax, const vector< pair<double, int> > &ibslVec,
		  const vector<bool> &masked, double mul) {

    const double cMoffset = 0.01;

    double pSum = 0;
    int ctr = 0;
    for (int j = 0; j < jMax; j++) {
        int hx = ibslVec[j].second;
        if (masked[hx/2]) {
            continue;
        }
        ctr++;
        double x = ibslVec[j].first + cMoffset;
        pVec[j] = max(exp(-mul / x), 1e-100);
        pSum += pVec[j];
    }
    for (int j = 0; j < jMax; j++) {
        int hx = ibslVec[j].second;
        if (masked[hx/2]) {
            pVec[j] = 0;
        } else {
            pVec[j] = pVec[j] / pSum;
        }
    }
}

vector< vector<double> > phase(const double params[],
                               const vector< vector< pair<double, int> > > &ibsl,
                               const vector< vector<double> > &dipCNPs,
                               vector< vector<double> > hapCNPs,
                               const vector<bool> &masked,
                               const vector<int> &order,
                               const int debug) {
    int H = ibsl.size();
    int Nchop = (int) params[0];
    double resampleP = params[1];
    double mul = params[2];

    vector< vector<double> > pVecs(2, vector<double> (Nchop));

    vector<double> popCNP;
    int numHaps = 0;
    for (int h = 0; h < H; h++) {
        if (!masked[h/2]) {
            numHaps++;
            if (hapCNPs[h].size() > popCNP.size()) {
                popCNP.resize(hapCNPs[h].size());
            }
            for (int c = 0; c < (int) hapCNPs[h].size(); c++) {
                popCNP[c] += hapCNPs[h][c];
            }
        }
        for (int c = 0; c < (int) popCNP.size(); c++) {
            popCNP[c] /= numHaps;
        }
    }

    for (int j = 0; j < H/2; j++) {

        int i = order[j];
        if (masked[i]) {
            continue;
        }
        int C = dipCNPs[i].size();
        vector< vector<double> > CNPs(2, vector<double> (C));
        vector< vector<double> > CNPpairs(C, vector<double> (C));

#pragma omp parallel for
        for (int k = 0; k < 2; k++) {
            vector<double> &pVec = pVecs[k];
            vector<double> &CNP = CNPs[k];
            int h = 2*i + k;
            int jMax = min(Nchop, (int) ibsl[h].size());
            computeProbs(pVec, jMax, ibsl[h], masked, mul);
            for (int j = 0; j < jMax; j++) {
                int hRef = ibsl[h][j].second;
                for (int c = 0; c < min(C, (int) hapCNPs[hRef].size()); c++) {
                    CNP[c] += pVec[j] * hapCNPs[hRef][c];
                }
            }
            for (int c = 0; c < C; c++) {
                CNP[c] += resampleP * popCNP[c];
            }
        }

        // create cross product
        for (int c1 = 0; c1 < C; c1++) {
            for (int c2 = 0; c2 < C; c2++) {
                CNPpairs[c1][c2] = CNPs[0][c1] * CNPs[1][c2];
            }
        }

        // subtract double-IBD
        for (int j1 = 0; j1 < min(Nchop, (int) ibsl[2*i].size()); j1++) {
            for (int j2 = 0; j2 < min(Nchop, (int) ibsl[2*i+1].size()); j2++) {
                int h1x = ibsl[2*i][j1].second;
                int h2x = ibsl[2*i+1][j2].second;
                if (h1x/2 == h2x/2) {
                    for (int c1 = 0; c1 < min(C, (int) hapCNPs[h1x].size()); c1++) {
                        for (int c2 = 0; c2 < min(C, (int) hapCNPs[h2x].size()); c2++) {
                            CNPpairs[c1][c2] -= (pVecs[0][j1]*hapCNPs[h1x][c1])*(pVecs[1][j2]*hapCNPs[h2x][c2]);
                        }
                    }
                }
            }
        }

        vector< vector<double> > hapCNPsIn(2);
        for (int k = 0; k < 2; k++) {
            hapCNPsIn[k] = hapCNPs[2*i+k];
            for (int c = 0; c < C; c++) {
                hapCNPs[2*i+k][c] = 0;
            }
        }

        // multiply in diploid CN probabilities
        double totP = 0;
        for (int c1 = 0; c1 < C; c1++) {
            for (int c2 = 0; c2 < C; c2++) {
                assert(CNPpairs[c1][c2]>=-1e-9);
                if (CNPpairs[c1][c2] < 0) {
                    CNPpairs[c1][c2] = 0;
                }
                CNPpairs[c1][c2] *= (c1+c2>=C) ? 0 : dipCNPs[i][c1+c2];
                totP += CNPpairs[c1][c2];
                hapCNPs[2*i][c1] += CNPpairs[c1][c2];
                hapCNPs[2*i+1][c2] += CNPpairs[c1][c2];
            }
        }

        if (totP < 0.1) {
            if (debug > 0) {
                cout << "Warning: normalization failed: i=" << i << " j=" << j << " C=" << C << " totP=" << totP << endl;
            }
            for (int k = 0; k < 2; k++) {
                hapCNPs[2*i+k] = hapCNPsIn[k];
            }
            continue;
        }

        // normalize hapCNPs
        for (int k = 0; k < 2; k++) {
            double check = 0;
            for (int c = 0; c < C; c++) {
                hapCNPs[2*i+k][c] /= totP;
                check += hapCNPs[2*i+k][c];
            }
            assert(fabs(1-check)<1e-6);
        }
    }
    return hapCNPs;
}

double meanCN(const vector <double> &CNPs) {
    double meanCN = 0;
    for (int c = 0; c < (int) CNPs.size(); c++) {
        meanCN += CNPs[c] * c;
    }
    return meanCN;
}

double r2imp(const vector< vector<double> > &hapCNPs,
             const vector<double> &impDipCNs,
	     const vector<int> &impInds,
             const vector<bool> &masked,
	     const vector< vector< pair<double, int> > > &ibsl,
	     int Nchop=100,
             double mul=2) {

    int impN = impInds.size();
    vector<double> estDipCNs(impN);

    vector<double> pVec(Nchop);

    for (int i = 0; i < impN; i++) {
        for (int h = 2*impInds[i]; h < 2*impInds[i]+2; h++) {
            int jMax = min(Nchop, (int) ibsl[h].size());
            computeProbs(pVec, jMax, ibsl[h], masked, mul);
            double pCheck = 0;
            for (int j = 0; j < jMax; j++) {
                estDipCNs[i] += pVec[j] * meanCN(hapCNPs[ibsl[h][j].second]);
                pCheck += pVec[j];
            }
            assert(fabs(pCheck-1) < 1e-6);
        }
    }
    return r2(estDipCNs, impDipCNs);
}

void impMissing(vector< vector<double> > &hapCNPs,
                const vector<bool> &masked,
		const vector< vector< pair<double, int> > > &ibsl,
                int Nchop,
                double mul) {

    vector<double> pVec(Nchop);

    for (int h = 0; h < (int) hapCNPs.size(); h++) {
        if (hapCNPs[h].empty()) {
            int jMax = min(Nchop, (int) ibsl[h].size());
            computeProbs(pVec, jMax, ibsl[h], masked, mul);
            double pCheck = 0;
            for (int j = 0; j < jMax; j++) {
                const vector <double> &CNPs = hapCNPs[ibsl[h][j].second];
                if (hapCNPs[h].size() < CNPs.size()) {
                    hapCNPs[h].resize(CNPs.size());
                }
                for (int c = 0; c < (int) CNPs.size(); c++) {
                    hapCNPs[h][c] += pVec[j] * CNPs[c];
                }
                pCheck += pVec[j];
            }
            assert(fabs(pCheck-1) < 1e-6);
        }
    }
}

void printIterOutput(int t,
                     const vector< vector<double> > &hapCNPs,
                     const double params[],
		     const vector<double> &impDipCNs,
                     const vector<int> &impInds,
		     const vector<bool> &masked,
		     const vector< vector< pair<double, int> > > &ibsl) {

    printf("iter %d:  %.3f  Nchop: %3d  resampleP: %.3f  invIBSmul: %.3f",
           t, r2imp(hapCNPs, impDipCNs, impInds, masked, ibsl), (int) params[0], params[1], params[2]);
    cout << endl;
}

// Algorithm entry point

extern vector< vector<double> > phaseImpCore(const vector< vector< pair<double, int> > >& ibsMatrix,
                                             const vector< vector<double> >& dipCNPs,
                                             const int nIterations,
                                             const int debug) {

    // Fixed parameters (currently)
    const int tBurn = 20;
    const int tFinal = 10;

    // Fitted per-site parameters
    const int N_PARAMS = 3;
    double params[N_PARAMS] = {50, 0.05, 1};
    const int N = dipCNPs.size();
    const int H = ibsMatrix.size();

    assert(H == 2*N);

    if (debug > 0) {
        cout << "Initial parameter values:" << endl;
        cout << "  Nchop = " << params[0] << endl;
        cout << "  resampleP = " << params[1] << endl;
        cout << "  invIBSmul = " << params[2] << endl;
        cout << "Iterations of simulated annealing: " << nIterations << endl;
    }

    // Initialize hapCNPs

    vector< vector<double> > hapCNPs(H);
    vector<bool> masked(H/2);
    int Nimpute = 0;
    for (int i = 0; i < H/2; i++) {
        if (dipCNPs[i].empty()) {
            masked[i] = 1;
            Nimpute++;
        } else {
            int maxCN = max_element(dipCNPs[i].begin(), dipCNPs[i].end()) - dipCNPs[i].begin();
            hapCNPs[2*i] = hapCNPs[2*i+1] = vector <double> ((maxCN+1), 1.0/(maxCN+1));
            hapCNPs[2*i].resize(dipCNPs[i].size());
            hapCNPs[2*i+1].resize(dipCNPs[i].size());
        }
    }
    if (Nimpute > 0 && debug > 0) {
        cout << "Missing CNLs for " << Nimpute << " individuals; will impute" << endl;
    }

    // mask individuals to impute
    vector<bool> maskedJustMiss = masked;
    vector<double> impDipCNs;
    vector<int> impInds;
    for (int i = 0; i < H/2; i += 5) {
        if (!masked[i]) {
            masked[i] = true;
            impDipCNs.push_back(meanCN(dipCNPs[i]));
            impInds.push_back(i);
        }
    }

    if (debug > 0) {
        cout << "Using " << impInds.size() << " individuals for imputation benchmark" << endl;
        cout << endl;
    }

    if (debug > 0) {
        cout << "Beginning phasing (" << tBurn << " burn-in iters; then parameter-optimization)" << endl;
        cout << "ITER #: IMP_R2 params" << endl;
        printIterOutput(0, hapCNPs, params, impDipCNs, impInds, masked, ibsMatrix);
    }

    vector<int> order(N);
    for (int i = 0; i < N; i++) {
        order[i] = i;
    }

    double r2best = 0;
    vector< vector<double> > hapCNPsBest;
    double paramsBest[N_PARAMS];
    int tBest = 0;

    // run phasing benchmark + parameter optimization
    for (int t = 0; t < nIterations; t++) {
        random_shuffle(order.begin(), order.end());

        vector< vector<double> > hapCNPsSame = phase(params, ibsMatrix, dipCNPs, hapCNPs, masked, order, debug);
        double r2same = r2imp(hapCNPsSame, impDipCNs, impInds, masked, ibsMatrix);
        double r2iter = r2same;
    
        if (t >= tBurn) {
            int pInd = t % N_PARAMS;
            double pMult = exp(2 * ((rand() % 1000) / 1000.0 - 0.5) * (nIterations-t) / nIterations);
            if (pInd == 0) {
                pMult = min(pMult, 100/params[pInd]); // Nchop <= 100
                pMult = max(pMult, 25/params[pInd]); // Nchop >= 25
            }
            if (pInd == 1) {
                pMult = min(pMult, 0.5/params[pInd]); // resampleP <= 0.5
            }
            params[pInd] *= pMult;
            vector< vector<double> > hapCNPsDiff = phase(params, ibsMatrix, dipCNPs, hapCNPs, masked, order, debug);
            double r2diff = r2imp(hapCNPsDiff, impDipCNs, impInds, masked, ibsMatrix);
            if (1e4*(r2same-r2diff)*max(0.1, t/(double) nIterations) < (rand()%1000)/1000.0) {
                hapCNPs = hapCNPsDiff;
                r2iter = r2diff;
            } else {
                hapCNPs = hapCNPsSame;
                params[pInd] /= pMult;
            }
        } else {
            hapCNPs = hapCNPsSame;
        }

        if (debug > 0) {
            printIterOutput(t+1, hapCNPs, params, impDipCNs, impInds, masked, ibsMatrix);
        }

        if (r2iter > r2best) {
            r2best = r2iter;
            hapCNPsBest = hapCNPs;
            memcpy(paramsBest, params, sizeof(params));
            tBest = t;
        }
    }

    hapCNPs = hapCNPsBest;
    memcpy(params, paramsBest, sizeof(params));
    if (debug > 0) {
        cout << "BEST ITER: " << endl;
        printIterOutput(tBest+1, hapCNPs, params, impDipCNs, impInds, masked, ibsMatrix);
        cout << endl;
    }

    int opt_Nchop = 0;
    double opt_mul = 0;
    double opt_impR2 = 0;
    for (int Nchop = 20; Nchop <= 200; Nchop += 20) {
        for (double mul = 0.1; mul < 20; mul *= 1.414) {
            double r2 = r2imp(hapCNPs, impDipCNs, impInds, masked, ibsMatrix, Nchop, mul);
            if (r2 > opt_impR2) {
                opt_impR2 = r2;
                opt_Nchop = Nchop;
                opt_mul = mul;
            }
        }
    }

    if (debug > 0) {
        for (int imp_version = 1; imp_version <= 3; imp_version++) {
            string ver_str;
            int Nchop;
            double mul;
            if (imp_version == 1) {
                ver_str = "OPT";
                Nchop = opt_Nchop;
                mul = opt_mul;
            } else if (imp_version == 2) {
                ver_str = "PHASE";
                Nchop = (int) params[0];
                mul = params[2];
            } else {
                ver_str = "FIXED";
                Nchop = 100;
                mul = 1;
            }
            printf("IMP_%-5s:  R2: %.3f  Nchop: %3d  mul: %5.2f\n",
                   ver_str.c_str(), r2imp(hapCNPs, impDipCNs, impInds, masked, ibsMatrix, Nchop, mul), Nchop, mul);
        }
    }

    // phase all individuals using optimized parameters
    if (debug > 0) {
        cout << endl;
        cout << "Beginning phasing iterations including previously-masked individuals" << endl;
        cout << endl;
        cout << "NOTE: test data is now used during phasing; benchmarks are only sanity-checks" << endl;
    }

    masked = maskedJustMiss;
    for (int t = 0; t < tFinal; t++) {
        random_shuffle(order.begin(), order.end());
        hapCNPs = phase(params, ibsMatrix, dipCNPs, hapCNPs, masked, order, debug);
        if (debug > 0) {
            printIterOutput(-(t+1), hapCNPs, params, impDipCNs, impInds, masked, ibsMatrix);
        }
    }

    // impute any individuals in IBS file who had missing CNLs
    if (debug > 0) {
        cout << endl;
        cout << "Imputing " << Nimpute << " individuals in IBS file with missing CNLs" << endl;
    }
    impMissing(hapCNPs, masked, ibsMatrix, opt_Nchop, opt_mul);

    return(hapCNPs);
}





