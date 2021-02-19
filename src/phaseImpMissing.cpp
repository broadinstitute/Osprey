// g++ -fopenmp -O3 -Wall phaseImpMissing.cpp -Iinclude -I/home/pl88/boost_1_58_0/install/include -Wl,-Bstatic -lboost_iostreams -lz -Wl,-Bdynamic -o phaseImpMissing -L/home/pl88/boost_1_58_0/install/lib -L/n/groups/price/poru/external_software/zlib/zlib-1.2.11

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

using namespace std;

#include "FileUtils.cpp"
#include "StringUtils.cpp"
#include "NumericUtils.cpp"

const int N_PARAMS = 3;

inline double sq(double x) { return x*x; }

double r2(vector <double> x, vector <double> y) {
  int N = x.size();
  double mu = NumericUtils::mean(x);
  for (int i = 0; i < N; i++) x[i] -= mu;
  mu = NumericUtils::mean(y);
  for (int i = 0; i < N; i++) y[i] -= mu;
  return sq(NumericUtils::dot(&x[0], &y[0], N)) / (NumericUtils::norm2(&x[0], N) * NumericUtils::norm2(&y[0], N));
}

void computeProbs(vector <double> &pVec, int jMax, const vector < pair <double, int> > &ibslVec,
		  const vector <bool> &masked, double mul) {

  const double cMoffset = 0.01;

  double pSum = 0; int ctr = 0;
  for (int j = 0; j < jMax; j++) {
    int hx = ibslVec[j].second; if (masked[hx/2]) continue;
    ctr++;
    double x = ibslVec[j].first + cMoffset;
    pVec[j] = max(exp(-mul / x), 1e-100);
    pSum += pVec[j];
  }
  for (int j = 0; j < jMax; j++) {
    int hx = ibslVec[j].second;
    if (masked[hx/2])
      pVec[j] = 0;
    else
      pVec[j] = pVec[j] / pSum;
  }
}

vector < vector <double> > phase(double params[],
				 const vector < vector < pair <double, int> > > &ibsl,
				 const vector < vector <double> > &dipCNPs,
				 vector < vector <double> > hapCNPs, const vector <bool> &masked,
				 const vector <int> &order) {
  int H = ibsl.size();
  int Nchop = (int) params[0]; double resampleP = params[1], mul = params[2];

  vector < vector <double> > pVecs(2, vector <double> ((int) Nchop));

  vector <double> popCNP;
  int numHaps = 0;
  for (int h = 0; h < H; h++)
    if (!masked[h/2]) {
      numHaps++;
      if (hapCNPs[h].size() > popCNP.size())
	popCNP.resize(hapCNPs[h].size());
      for (int c = 0; c < (int) hapCNPs[h].size(); c++)
	popCNP[c] += hapCNPs[h][c];
    }
  for (int c = 0; c < (int) popCNP.size(); c++)
    popCNP[c] /= numHaps;

  for (int j = 0; j < H/2; j++) {
    int i = order[j];
    if (masked[i]) continue;

    int C = dipCNPs[i].size();
    vector < vector <double> > CNPs(2, vector <double> (C));
    vector < vector <double> > CNPpairs(C, vector <double> (C));

#pragma omp parallel for
    for (int k = 0; k < 2; k++) {
      vector <double> &pVec = pVecs[k];
      vector <double> &CNP = CNPs[k];
      int h = 2*i + k;
      int jMax = min(Nchop, (int) ibsl[h].size());
      computeProbs(pVec, jMax, ibsl[h], masked, mul);
      for (int j = 0; j < jMax; j++) {
	int hRef = ibsl[h][j].second;
	for (int c = 0; c < min(C, (int) hapCNPs[hRef].size()); c++)
	  CNP[c] += pVec[j] * hapCNPs[hRef][c];
      }
      for (int c = 0; c < C; c++) {
	CNP[c] += resampleP * popCNP[c];
      }
    }

    // create cross product
    for (int c1 = 0; c1 < C; c1++)
      for (int c2 = 0; c2 < C; c2++)
	CNPpairs[c1][c2] = CNPs[0][c1] * CNPs[1][c2];

    // subtract double-IBD
    for (int j1 = 0; j1 < min(Nchop, (int) ibsl[2*i].size()); j1++)
      for (int j2 = 0; j2 < min(Nchop, (int) ibsl[2*i+1].size()); j2++) {
	int h1x = ibsl[2*i][j1].second;
	int h2x = ibsl[2*i+1][j2].second;
	if (h1x/2 == h2x/2)
	  for (int c1 = 0; c1 < min(C, (int) hapCNPs[h1x].size()); c1++)
	    for (int c2 = 0; c2 < min(C, (int) hapCNPs[h2x].size()); c2++)
	      CNPpairs[c1][c2] -= (pVecs[0][j1]*hapCNPs[h1x][c1])*(pVecs[1][j2]*hapCNPs[h2x][c2]);
      }
    
    // reset hapCNPs for output
    for (int k = 0; k < 2; k++)
      for (int c = 0; c < C; c++)
	hapCNPs[2*i+k][c] = 0;

    // multiply in diploid CN probabilities
    double totP = 0;
    for (int c1 = 0; c1 < C; c1++)
      for (int c2 = 0; c2 < C; c2++) {
	assert(CNPpairs[c1][c2]>=-1e-9);
	if (CNPpairs[c1][c2] < 0) CNPpairs[c1][c2] = 0;
	CNPpairs[c1][c2] *= c1+c2>=C ? 0 : dipCNPs[i][c1+c2];
	totP += CNPpairs[c1][c2];
	hapCNPs[2*i][c1] += CNPpairs[c1][c2];
	hapCNPs[2*i+1][c2] += CNPpairs[c1][c2];
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
  for (int c = 0; c < (int) CNPs.size(); c++)
    meanCN += CNPs[c] * c;
  return meanCN;
}

double r2imp(const vector < vector <double> > &hapCNPs, const vector <double> &impDipCNs,
	     const vector <int> &impInds, const vector <bool> &masked,
	     const vector < vector < pair <double, int> > > &ibsl,
	     int Nchop=100, double mul=2) {

  int impN = impInds.size();
  vector <double> estDipCNs(impN);

  vector <double> pVec(Nchop);

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

void impMissing(vector < vector <double> > &hapCNPs, const vector <bool> &masked,
		const vector < vector < pair <double, int> > > &ibsl, int Nchop, double mul) {

  vector <double> pVec(Nchop);

  for (int h = 0; h < (int) hapCNPs.size(); h++) {
    if (hapCNPs[h].empty()) {
      int jMax = min(Nchop, (int) ibsl[h].size());
      computeProbs(pVec, jMax, ibsl[h], masked, mul);
      double pCheck = 0;
      for (int j = 0; j < jMax; j++) {
	const vector <double> &CNPs = hapCNPs[ibsl[h][j].second];
	if (hapCNPs[h].size() < CNPs.size())
	  hapCNPs[h].resize(CNPs.size());
	for (int c = 0; c < (int) CNPs.size(); c++)
	  hapCNPs[h][c] += pVec[j] * CNPs[c];
	pCheck += pVec[j];
      }
      assert(fabs(pCheck-1) < 1e-6);
    }
  }
}

void printIterOutput(int t, const vector < vector <double> > &hapCNPs, const double params[],
		     const vector <double> &impDipCNs, const vector <int> &impInds,
		     const vector <bool> &masked,
		     const vector < vector < pair <double, int> > > &ibsl) {

  printf("iter %d:  %.3f  Nchop: %3d  resampleP: %.3f  invIBSmul: %.3f",
	 t, r2imp(hapCNPs, impDipCNs, impInds, masked, ibsl),
	 (int) params[0], params[1], params[2]);
  cout << endl;
}

int main(int argc, char *argv[]) {
  if (argc != 6) {
    cerr << "Usage:" << endl;
    cerr << "- arg1: IBS file (from computeIBS)" << endl;
    cerr << "- arg2: CNL file (ID CNL0 CNL1 ...)" << endl;
    cerr << "- arg3: iters" << endl;
    cerr << "- arg4: threads" << endl;
    cerr << "- arg5: output hapCNP file" << endl;
    exit(1);
    // cerr << "- arg3: Nchop (e.g., 50)" << endl;
    // cerr << "- arg4: resampleP (e.g., 0.05)" << endl;
    // cerr << "- arg5: inverse IBS multiplier (e.g., 2)" << endl;
    // cerr << "- arg7: random seed" << endl;
  }

  const double weight0err = 0.5, weightTotalLen = 0.25;
  int iters, threads;
  double params[N_PARAMS] = {50, 0.05, 1};

  const char *IBSfile = argv[1];
  const char *CNLfile = argv[2];
  sscanf(argv[3], "%d", &iters);
  sscanf(argv[4], "%d", &threads);
  const char *outFile = argv[5];

  const int tBurn = 20, tFinal = 10;

  assert(iters>=20);
  assert(iters<=1000);
  assert(threads>=1);
  assert(threads<=64);

  /*
  sscanf(argv[3], "%lf", &params[0]);
  sscanf(argv[4], "%lf", &params[1]);
  sscanf(argv[5], "%lf", &params[2]);
  sscanf(argv[7], "%d", &seed);
  cout << "Setting weight0err: " << weight0err << endl;
  cout << "Setting weightTotalLen: " << weightTotalLen << endl;
  */
  cout << "Initial parameter values:" << endl;
  cout << "  Nchop = " << params[0] << endl;
  cout << "  resampleP = " << params[1] << endl;
  cout << "  invIBSmul = " << params[2] << endl;
  cout << "Iterations of simulated annealing: " << iters << endl;
  //cout << "Random seed: " << seed << endl;
  cout << "Setting number of threads to " << threads << endl;
  omp_set_num_threads(threads);

  srand(1); //srand((int) (seed + 1000*weight0err + 1000000*weightTotalLen));

  FileUtils::AutoGzIfstream fin;

  // read IDs
  vector <string> IDs; map <string, int> IDtoInd;
  fin.openOrExit(IBSfile);
  string ID; string line;
  while (fin >> ID) {
    IDtoInd[ID] = IDs.size();
    IDs.push_back(ID);
    getline(fin, line);
    getline(fin, line); // ignore line for second haplotype
  }
  fin.close();
  int N = IDs.size(); int H = 2*N;
  cout << "Read " << N << " IDs" << endl;

  istringstream iss(line);
  int Htop = 0;
  while (iss >> line) Htop++;
  Htop /= 6;
  cout << "Reading " << Htop << " IBS partners per line" << endl;

  // read IBS lengths
  vector < vector < pair <double, int> > > ibsl(H, vector < pair <double, int> > (Htop));
  fin.openOrExit(IBSfile);
  for (int h = 0; h < H; h++) {
    int hap12;
    fin >> ID >> hap12;
    assert(ID == IDs[h/2]);
    for (int i = 0; i < Htop; i++) {
      double cMs[4];
      fin >> ID >> hap12 >> cMs[0] >> cMs[1] >> cMs[2] >> cMs[3];
      double lenLeft = weight0err * -cMs[1] + (1-weight0err) * -cMs[0];
      double lenRight = weight0err * cMs[2] + (1-weight0err) * cMs[3];
      double lenTot = lenLeft + lenRight;
      double lenMin = min(lenLeft, lenRight);
      ibsl[h][i].first = weightTotalLen * lenTot + (1-weightTotalLen) * lenMin;
      ibsl[h][i].second = 2*IDtoInd[ID] + hap12-1;
      if (ibsl[h][i].second == -1) {
	cout << "ERROR: " << h << " " << i << " " << ID << " " << hap12 << endl;
	exit(1);
      }
    }
    sort(ibsl[h].begin(), ibsl[h].end(), greater < pair <double, int> > ());
  }
  fin.close();

  // read dipCNLs
  vector < vector <double> > dipCNPs(H/2);
  vector < vector <double> > hapCNPs(H);
  fin.openOrExit(CNLfile);
  int numFound = 0, numConf = 0;
  while (fin >> ID) {
    getline(fin, line);
    if (IDtoInd.find(ID) != IDtoInd.end()) {
      istringstream iss(line);
      vector <double> CNPs;
      double CNL, totP = 0;
      while (iss >> CNL) {
	CNPs.push_back(pow(10.0, CNL));
	totP += CNPs.back();
      }
      double maxP = 0;
      for (int c = 0; c < (int) CNPs.size(); c++) {
	CNPs[c] /= totP;
	if (CNPs[c] > maxP)
	  maxP = CNPs[c];
      }
      dipCNPs[IDtoInd[ID]] = CNPs;
      numFound++;
      if (maxP > 0.95) numConf++;
    }
  }
  fin.close();
  cout << "Read diploid CNLs for " << numFound << " individuals" << endl;
  cout << "(" << numConf << " individuals have 95% confident dipCNs)" << endl;
  vector <bool> masked(H/2);
  int Nimpute = 0;
  for (int i = 0; i < H/2; i++) {
    if (dipCNPs[i].empty()) {
      //cout << "Missing CNLs for indiv " << IDs[i] << "; will impute" << endl;
      masked[i] = 1;
      Nimpute++;
    }
    else {
      int maxCN = max_element(dipCNPs[i].begin(), dipCNPs[i].end()) - dipCNPs[i].begin();
      hapCNPs[2*i] = hapCNPs[2*i+1] = vector <double> ((maxCN+1), 1.0/(maxCN+1));
      hapCNPs[2*i].resize(dipCNPs[i].size());
      hapCNPs[2*i+1].resize(dipCNPs[i].size());
    }
  }
  if (Nimpute)
    cout << "Missing CNLs for " << Nimpute << " individuals; will impute" << endl;
  vector <bool> maskedJustMiss = masked;

  // mask indivs to impute
  vector <double> impDipCNs; vector <int> impInds;
  for (int i = 0; i < H/2; i += 5)
    if (!masked[i]) {
      masked[i] = true;
      impDipCNs.push_back(meanCN(dipCNPs[i]));
      impInds.push_back(i);
    }
  cout << "Using " << impInds.size() << " individuals for imputation benchmark" << endl;
  cout << endl;



  cout << "Beginning phasing (" << tBurn << " burn-in iters; then parameter-optimization)" << endl;
  cout << "ITER #: IMP_R2 params" << endl;

  printIterOutput(0, hapCNPs, params, impDipCNs, impInds, masked, ibsl);

  vector <int> order(N);
  for (int i = 0; i < N; i++) order[i] = i;

  double r2best = 0; vector < vector <double> > hapCNPsBest; double paramsBest[N_PARAMS];
  int tBest = 0;

  // run phasing benchmark + parameter optimization
  for (int t = 0; t < iters; t++) {
    random_shuffle(order.begin(), order.end());

    vector < vector <double> > hapCNPsSame = phase(params, ibsl, dipCNPs, hapCNPs, masked, order);
    double r2same = r2imp(hapCNPsSame, impDipCNs, impInds, masked, ibsl), r2iter = r2same;
    
    if (t >= tBurn) {
      int pInd = t % N_PARAMS;
      double pMult = exp(2 * ((rand() % 1000) / 1000.0 - 0.5) * (iters-t) / iters);
      if (pInd == 0) {
	pMult = min(pMult, 100/params[pInd]); // Nchop <= 100
	pMult = max(pMult, 25/params[pInd]); // Nchop >= 25
      }
      if (pInd == 1) {
	pMult = min(pMult, 0.5/params[pInd]); // resampleP <= 0.5
      }
      params[pInd] *= pMult;
      vector < vector <double> > hapCNPsDiff = phase(params, ibsl, dipCNPs, hapCNPs, masked, order);
      double r2diff = r2imp(hapCNPsDiff, impDipCNs, impInds, masked, ibsl);
      if (1e4*(r2same-r2diff)*max(0.1, t/(double) iters) < (rand()%1000)/1000.0) {
	hapCNPs = hapCNPsDiff;
	r2iter = r2diff;
      }
      else {
	hapCNPs = hapCNPsSame;
	params[pInd] /= pMult;
      }
    }
    else
      hapCNPs = hapCNPsSame;

    printIterOutput(t+1, hapCNPs, params, impDipCNs, impInds, masked, ibsl);
    
    if (r2iter > r2best) {
      r2best = r2iter;
      hapCNPsBest = hapCNPs;
      memcpy(paramsBest, params, N_PARAMS*sizeof(params[0]));
      tBest = t;
    }
  }

  cout << "BEST ITER: " << endl;
  hapCNPs = hapCNPsBest;
  memcpy(params, paramsBest, N_PARAMS*sizeof(params[0]));
  printIterOutput(tBest+1, hapCNPs, params, impDipCNs, impInds, masked, ibsl);
  cout << endl;

  int opt_Nchop = 0; double opt_mul = 0, opt_impR2 = 0;
  for (int Nchop = 20; Nchop <= 200; Nchop += 20)
    for (double mul = 0.1; mul < 20; mul *= 1.414) {
      double r2 = r2imp(hapCNPs, impDipCNs, impInds, masked, ibsl, Nchop, mul);
      if (r2 > opt_impR2) {
	opt_impR2 = r2;
	opt_Nchop = Nchop;
	opt_mul = mul;
      }
    }
  for (int imp_version = 1; imp_version <= 3; imp_version++) {
    string ver_str;
    int Nchop; double mul;
    if (imp_version == 1) {
      ver_str = "OPT"; Nchop = opt_Nchop; mul = opt_mul;
    }
    else if (imp_version == 2) {
      ver_str = "PHASE"; Nchop = (int) params[0]; mul = params[2];
    }
    else {
      ver_str = "FIXED"; Nchop = 100; mul = 1;
    }
    printf("IMP_%-5s:  R2: %.3f  Nchop: %3d  mul: %5.2f\n",
	   ver_str.c_str(), r2imp(hapCNPs, impDipCNs, impInds, masked, ibsl, Nchop, mul),
	   Nchop, mul);
  }

  // phase all individuals using optimized parameters
  cout << endl;
  cout << "Beginning phasing iterations including previously-masked individuals" << endl;
  cout << endl;

  cout << "NOTE: test data is now used during phasing; benchmarks are only sanity-checks" << endl;
  masked = maskedJustMiss;
  for (int t = 0; t < tFinal; t++) {
    random_shuffle(order.begin(), order.end());
    hapCNPs = phase(params, ibsl, dipCNPs, hapCNPs, masked, order);
    printIterOutput(-(t+1), hapCNPs, params, impDipCNs, impInds, masked, ibsl);
  }

  // impute any individuals in IBS file who had missing CNLs
  cout << endl;
  cout << "Imputing " << Nimpute << " individuals in IBS file with missing CNLs" << endl;
  impMissing(hapCNPs, masked, ibsl, opt_Nchop, opt_mul);

  FileUtils::AutoGzOfstream fout; fout.openOrExit(outFile);
  fout << std::setprecision(4) << std::fixed;
  fout << "ID\tHAP\tSTATUS\tTOP_CN\tTOP_P\tALL_P" << endl;
  for (int h = 0; h < H; h++) {
    fout << IDs[h/2] << "\t" << h%2+1 << "\t" << (maskedJustMiss[h/2] ? "IMPUTED" : "PHASED");
    int cMax = max_element(hapCNPs[h].begin(), hapCNPs[h].end()) - hapCNPs[h].begin();
    fout << "\t" << cMax << "\t" << hapCNPs[h][cMax] << "\t";
    for (int c = 0; c < (int) hapCNPs[h].size(); c++) {
      if (c) fout << ",";
      fout << hapCNPs[h][c];
    }
    fout << endl;
  }
  fout.close();

  cout << endl << "Successfully completed phasing" << endl;

  return 0;
}
