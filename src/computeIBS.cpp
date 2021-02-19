// g++ -fopenmp -O3 -Wall computeIBS.cpp -Iinclude -I/home/pl88/boost_1_58_0/install/include -Wl,-Bstatic -lboost_iostreams -lz -Wl,-Bdynamic -o computeIBS -L/home/pl88/boost_1_58_0/install/lib -L/n/groups/price/poru/external_software/zlib/zlib-1.2.11

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <cstdio>
#include <cmath>

#include "omp.h"

#include "FileUtils.cpp"
#include "StringUtils.cpp"
#include "Timer.cpp"
#include "MapInterpolater.cpp"

using namespace std;

typedef unsigned char uchar;

const int MAX_ERR = 2;

struct IBSmatch {
  int hap;
  int bpErrs[2][MAX_ERR];
  double cMerrs[2][MAX_ERR];
  double len;
  bool operator < (const IBSmatch &match2) const {
    /*
    double len1 = 0, len2 = 0;
    for (int e = 0; e < MAX_ERR; e++) {
      len1 += cMerrs[1][e] - cMerrs[0][e];
      len2 += match2.cMerrs[1][e] - match2.cMerrs[0][e];
    }
    if (len1 != len2) return len1 > len2;
    */
    if (len != match2.len) return len > match2.len;
    else return hap < match2.hap;
  }
};

void updateGenos(uchar *genos, int N, int M, int m, const string &line) {
  for (int n = 0; n < N; n++) {
    genos[(2LL*n)*M + m] = line[4*n+1]-'0';
    genos[(2LL*n+1)*M + m] = line[4*n+3]-'0';
  }
}

string findBestIBS(int h1, uchar *genoPtrs[2], const vector <string> &sampleIDs,
		   const vector <int> bps[2], const vector <double> cMs[2], double weight0err,
		   double weightTotalLen, const vector <bool> &isRefSample) {
  ostringstream oss;
  //cout << bps[0].back() << " " << bps[0][0] << " " << bps[1][0] << " " << bps[1].back() << endl;
  const int N = sampleIDs.size();
  const int H = 2*N;
  const int Mflanks[2] = {(int) bps[0].size(), (int) bps[1].size()};
  //cout << "N = " << N << " Mflanks = " << Mflanks[0] << " " << Mflanks[1] << endl;

  set <IBSmatch> bestMatches;
  double num_checks = 0;
  const int MAX_MATCHES = 200;
  for (int h2 = 0; h2 < H; h2++) {
    if (h2 == h1) continue;
    if (!isRefSample[h2/2]) continue;
    IBSmatch matchInfo; matchInfo.hap = h2;
    for (int k = 0; k < 2; k++) {
      if (Mflanks[k] == 0) { // at telomere; no flanking region
	for (int errs = 0; errs < MAX_ERR; errs++) {
	  matchInfo.bpErrs[k][errs] = 0;
	  matchInfo.cMerrs[k][errs] = 0;
	}
	continue;
      }
      int errs = 0;
      const uchar *genos_h1 = genoPtrs[k] + h1 * (long long) Mflanks[k];
      const uchar *genos_h2 = genoPtrs[k] + h2 * (long long) Mflanks[k];
      for (int m = 0; m < Mflanks[k]; m++) {
	num_checks++;
	if (genos_h1[m] != genos_h2[m]) {
	  matchInfo.bpErrs[k][errs] = bps[k][m];
	  matchInfo.cMerrs[k][errs] = cMs[k][m];
	  errs++;
	  if (errs == MAX_ERR)
	    break;
	}
      }
      while (errs < MAX_ERR) {
	matchInfo.bpErrs[k][errs] = bps[k].back();
	matchInfo.cMerrs[k][errs] = cMs[k].back();
	errs++;
      }
      double lenLeft = weight0err*-matchInfo.cMerrs[0][0] + (1-weight0err)*-matchInfo.cMerrs[0][1];
      double lenRight = weight0err*matchInfo.cMerrs[1][0] + (1-weight0err)*matchInfo.cMerrs[1][1];
      double lenTot = lenLeft + lenRight;
      double lenMin = min(lenLeft, lenRight);
      // randomize for sorting
      matchInfo.len = weightTotalLen*lenTot + (1-weightTotalLen)*lenMin + (h1*h2%1000003)*1e-12;
    }
    bestMatches.insert(matchInfo);
    if (bestMatches.size() > MAX_MATCHES)
      bestMatches.erase(--bestMatches.end());
  }
  //cout << "avg_checks: " << num_checks / H << endl;
  oss << sampleIDs[h1/2] << " " << h1%2+1;
  oss << std::setprecision(3) << std::fixed;
  int ctr = 0;
  for (set <IBSmatch>::iterator it = bestMatches.begin(); it != bestMatches.end(); it++) {
    ctr++;
    //cout << it->hap << " " << it->bpErrs[0][1] << " " << it->bpErrs[0][0] << " " << it->bpErrs[1][0] << " " << it->bpErrs[1][1] << " " << it->cMerrs[0][1] << " " << it->cMerrs[0][0] << " " << it->cMerrs[1][0] << " " << it->cMerrs[1][1] << endl;
    oss << " " << sampleIDs[it->hap/2] << " " << it->hap%2+1 << " " << it->cMerrs[0][1] << " " << it->cMerrs[0][0] << " " << it->cMerrs[1][0] << " " << it->cMerrs[1][1];
    //if (ctr == 5) break;
  }

  return oss.str();
}

int main(int argc, char *argv[]) {

  if (argc != 7 && argc != 8) {
    cerr << "ERROR: 6 or 7 args required" << endl;
    cerr << "- arg1: input vcf[.gz] (including VNTR + ~2cM on each side)" << endl;
    cerr << "- arg2: genetic map file (genetic_map_hg##_withX.txt.gz)" << endl;
    cerr << "- arg3: VNTR start (to exclude SNPs within VNTR)" << endl;
    cerr << "- arg4: VNTR end (to exclude SNPs within VNTR)" << endl;
    cerr << "- arg5: number of threads" << endl;
    cerr << "- arg6: output IBS file" << endl;
    cerr << "- (optional) arg7: file containing subset of samples to use as ref" << endl;
    exit(1);
  }
  const char *vcfFile = argv[1];
  const char *geneticMapFile = argv[2];
  int vntrStart; sscanf(argv[3], "%d", &vntrStart);
  int vntrEnd; sscanf(argv[4], "%d", &vntrEnd);
  int threads; sscanf(argv[5], "%d", &threads);
  const char *outFile = argv[6];
  const char *refFile = argc==8 ? argv[7] : NULL;

  assert(vntrStart>=1);
  assert(vntrEnd>=1);
  assert(vntrStart<=vntrEnd);
  assert(threads>=1);
  assert(threads<=64);

  const double weight0err = 0.5;
  const double weightTotalLen = 0.25;

  omp_set_num_threads(threads);
  cout << "Set number of threads to " << threads << endl;

  Timer timer;
  
  FileUtils::AutoGzIfstream fin;
  fin.openOrExit(vcfFile);

  /***** read VCF file: first pass (count variants before and after VNTR *****/

  int Mflanks[2] = {0, 0};
  string chrStr, line; int pos; string POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT;
  bool headerFound = false;
  vector <string> sampleIDs;
  vector <int> bps[2];
  while (fin >> chrStr) {
    if (headerFound) {
      fin >> pos;
      if (pos < vntrStart) { Mflanks[0]++; bps[0].push_back(pos); }
      if (pos > vntrEnd) { Mflanks[1]++; bps[1].push_back(pos); }
      getline(fin, line);
    }
    else {
      if (chrStr == "#CHROM") {
	headerFound = true;
	fin >> POS >> ID >> REF >> ALT >> QUAL >> FILTER >> INFO >> FORMAT;
	getline(fin, line);
	istringstream iss(line);
	string sampleID;
	while (iss >> sampleID) {
	  sampleIDs.push_back(sampleID);
	}
	cout << "Reading data for " << sampleIDs.size() << " samples" << endl;
      }
    }
  }
  reverse(bps[0].begin(), bps[0].end());
  int chr;
  if (chrStr == "chrX" || chrStr == "X")
    chr = 23;
  else if (chrStr[0]=='c')
    sscanf(chrStr.c_str(), "chr%d", &chr);
  else
    sscanf(chrStr.c_str(), "%d", &chr);

  cout << "Found " << Mflanks[0] << " variants before " << chr << ":" << vntrStart << endl;
  cout << "Found " << Mflanks[1] << " variants after " << chr << ":" << vntrEnd << endl;
  fin.close();

  cout << "Finished first pass VCF read: " << timer.update_time() << " sec" << endl;

  const int N = sampleIDs.size();
  const int H = 2*N;

  /***** read ref sample subset (if provided) *****/
  vector <bool> isRefSample;
  if (refFile != NULL) {
    isRefSample = vector <bool> (N, false);
    map <string, int> IDtoInd;
    for (int i = 0; i < N; i++)
      IDtoInd[sampleIDs[i]] = i;
    fin.openOrExit(refFile);
    while (fin >> ID) {
      if (IDtoInd.find(ID) == IDtoInd.end()) {
	cerr << "ERROR: VCF file does not contain sample " << ID << " listed in ref file" << endl;
	exit(1);
      }
      isRefSample[IDtoInd[ID]] = true;
    }
    fin.close();
  }
  else
    isRefSample = vector <bool> (N, true);

  /***** interpolate genetic map coordinates *****/

  cout << "Filling in genetic map coordinates using reference file:" << endl;
  cout << "  " << geneticMapFile << endl;
  Genetics::MapInterpolater mapInterpolater(geneticMapFile);
  double cMmid = mapInterpolater.interp(chr, (vntrStart+vntrEnd)/2);
  vector <double> cMs[2];
  for (int k = 0; k < 2; k++) {
    cMs[k].resize(Mflanks[k]);
    for (int m = 0; m < Mflanks[k]; m++)
      cMs[k][m] = 100*(mapInterpolater.interp(chr, bps[k][m])-cMmid);
  }

  cout << "Finished interpolating coordinates: " << timer.update_time() << " sec" << endl;

  /***** read VCF file: second pass (store genotypes) *****/

  uchar *genoPtrs[2];
  for (int k = 0; k < 2; k++) {
    genoPtrs[k] = new uchar[H * (long long) Mflanks[k]];
    memset(genoPtrs[k], 0, H * (long long) Mflanks[k]);
  }
  int Mcurs[2] = {Mflanks[0], 0};

  fin.openOrExit(vcfFile);
  headerFound = false;
  while (fin >> chrStr) {
    if (headerFound) {
      fin >> pos >> ID >> REF >> ALT >> QUAL >> FILTER >> INFO >> FORMAT;
      getline(fin, line);
      if (pos < vntrStart) {
	Mcurs[0]--;
	updateGenos(genoPtrs[0], N, Mflanks[0], Mcurs[0], line);
      }
      if (pos > vntrEnd) {
	updateGenos(genoPtrs[1], N, Mflanks[1], Mcurs[1], line);
	Mcurs[1]++;
      }	
    }
    else {
      if (chrStr == "#CHROM") {
	headerFound = true;
	getline(fin, line);
      }
    }
  }
  fin.close();

  cout << "Finished second pass VCF read: " << timer.update_time() << " sec" << endl;

  /***** compute IBS *****/

  vector <string> outStrs(H);
#pragma omp parallel for
  for (int h1 = 0; h1 < H; h1++) {
    outStrs[h1] = findBestIBS(h1, genoPtrs, sampleIDs, bps, cMs, weight0err, weightTotalLen,
			      isRefSample);
    //cout << "." << flush;
  }
  //cout << endl;

  cout << "Finished IBS computations: " << timer.update_time() << " sec" << endl;

  FileUtils::AutoGzOfstream fout; fout.openOrExit(outFile);
  for (int h1 = 0; h1 < H; h1++)
    fout << outStrs[h1] << endl;
  fout.close();

  cout << "Finished writing output: " << timer.update_time() << " sec" << endl;

  for (int k = 0; k < 2; k++)
    delete[] genoPtrs[k];

  return 0;
}
