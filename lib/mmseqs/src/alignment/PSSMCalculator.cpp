//
// Created by mad on 3/24/15.
//
#include "PSSMCalculator.h"
#include "simd.h"
#include "MathUtil.h"
#include "BaseMatrix.h"
#include "Sequence.h"
#include "Util.h"
#include "Debug.h"
#include "MultipleAlignment.h"

PSSMCalculator::PSSMCalculator(SubstitutionMatrix *subMat, size_t maxSeqLength, size_t maxSetSize, int pcmode,
                               MultiParam<PseudoCounts> pca, MultiParam<PseudoCounts> pcb
#ifdef GAP_POS_SCORING
                               , int gapOpen
                               , int gapPseudoCount
#endif
    ) : subMat(subMat),  ps(NULL), maxSeqLength(maxSeqLength)
#ifdef GAP_POS_SCORING
    , gapOpen(gapOpen)
    , gapPseudoCount(gapPseudoCount)
#endif
    , pca(pca), pcb(pcb) {

    if (pcmode == Parameters::PCMODE_CONTEXT_SPECIFIC) {
        ps = new CSProfile(maxSeqLength + 1);
    }
    this->maxSeqLength = maxSeqLength;
    this->maxSetSize = maxSetSize;
    this->profile            = new float[(maxSeqLength + 1) * Sequence::PROFILE_AA_SIZE];
    this->Neff_M             = new float[(maxSeqLength + 1)];
    this->seqWeight          = (float*)malloc(maxSetSize * sizeof(float));
    this->pssm               = new char[(maxSeqLength + 1) * Sequence::PROFILE_AA_SIZE];
    this->consensusSequence  = new unsigned char[maxSeqLength + 1];
    this->matchWeight        = (float *) malloc_simd_float(Sequence::PROFILE_AA_SIZE * (maxSeqLength + 1) * sizeof(float));
    this->pseudocountsWeight = (float *) malloc_simd_float(Sequence::PROFILE_AA_SIZE * (maxSeqLength + 1) * sizeof(float));
    this->counts             = (float *) mem_align(ALIGN_FLOAT, (Sequence::PROFILE_AA_SIZE + 4) * (maxSeqLength + 1) * sizeof(float));
    memset(this->counts, 0, (Sequence::PROFILE_AA_SIZE + 4) * (maxSeqLength + 1) * sizeof(float));
    this->nseqs              = new int[maxSeqLength + 1];
    unsigned int NAA_ALIGNSIZE = ((((MultipleAlignment::NAA + 3) + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT) * sizeof(float);
    NAA_ALIGNSIZE = ((NAA_ALIGNSIZE + ALIGN_FLOAT - 1) / ALIGN_FLOAT) * ALIGN_FLOAT;
    w_contrib         = new float*[maxSeqLength + 1];
    w_contrib_backing = (unsigned char*)mem_align(ALIGN_FLOAT, NAA_ALIGNSIZE * (maxSeqLength + 1));
    for (size_t j = 0; j < (maxSeqLength + 1); j++) {
        w_contrib[j] = (float*)(w_contrib_backing + (NAA_ALIGNSIZE * j));
    }
    wi = (float*)malloc(maxSetSize * sizeof(float));
    naa = new int[maxSeqLength + 1];
    f = malloc_matrix<float>(maxSeqLength + 1, MultipleAlignment::NAA + 3);
    n = new int*[maxSeqLength + 2];
    n_backing = (unsigned char*)mem_align(ALIGN_INT, NAA_ALIGNSIZE * (maxSeqLength + 2));
    for (size_t j = 0; j < (maxSeqLength + 2); j++) {
        n[j] = (int*)(n_backing + (NAA_ALIGNSIZE * j));
    }
    this->pcmode = pcmode;
    this->pca = pca;
    this->pcb = pcb;
#ifdef GAP_POS_SCORING
    gDel = new uint8_t[(maxSeqLength + 1)];
    gIns = new uint8_t[(maxSeqLength + 1)];
    gapWeightsIns.reserve(maxSeqLength + 1);
#endif
}

PSSMCalculator::~PSSMCalculator() {
    if(ps != NULL){
        delete ps;
    }
    delete[] profile;
    delete[] Neff_M;
    free(seqWeight);
    delete[] pssm;
    delete[] nseqs;
    delete[] consensusSequence;
    free(matchWeight);
    free(counts);
    free(pseudocountsWeight);
    free(w_contrib_backing);
    delete[] w_contrib;
    free(wi);
    delete[] naa;
    free(n_backing);
    delete[] n;
    free(f);
#ifdef GAP_POS_SCORING
    delete[] gDel;
    delete[] gIns;
#endif
}

//    PSSMCalculator::Profile PSSMCalculator::computePSSMFromMSA(size_t setSize,
//                                                               size_t queryLength,
//                                                               const char **msaSeqs,
//                                                               bool wg, 0.0) {
//        increaseSetSize(setSize);
//        // Quick and dirty calculation of the weight per sequence wg[k]
//        computeSequenceWeights(seqWeight, queryLength, setSize, msaSeqs);
//        MathUtil::NormalizeTo1(seqWeight, setSize);
//        if (wg == false) {
//            // compute context specific counts and Neff
//            computeContextSpecificWeights(matchWeight, seqWeight, Neff_M, queryLength, setSize, msaSeqs);
//        } else {
//            // compute matchWeight based on sequence weight
//            computeMatchWeights(matchWeight, seqWeight, setSize, queryLength, msaSeqs);
//            // compute NEFF_M
//            computeNeff_M(matchWeight, seqWeight, Neff_M, queryLength, setSize, msaSeqs);
//        }
//        // compute consensus sequence
//        computeConsensusSequence(consensusSequence, matchWeight, queryLength, subMat->pBack, subMat->num2aa);
//        if(pcmode == Parameters::PCMODE_SUBSTITUTION_SCORE && pca.values.normal() > 0.0){
//            // add pseudocounts (compute the scalar product between matchWeight and substitution matrix with pseudo counts)
//            preparePseudoCounts(matchWeight, pseudocountsWeight, Sequence::PROFILE_AA_SIZE, queryLength, (const float **) subMat->subMatrixPseudoCounts);
//            //    SubstitutionMatrix::print(subMat->subMatrixPseudoCounts, subMat->num2aa, 20 );
//            computePseudoCounts(profile, matchWeight, pseudocountsWeight, Sequence::PROFILE_AA_SIZE, Neff_M, queryLength, pca.values.normal(), pcb.values.normal());
//        }else if (pcmode == Parameters::PCMODE_CONTEXT_SPECIFIC && pca.values.cs() > 0.0){
//            fillCounteProfile(counts, matchWeight, Neff_M, queryLength);
//            float * csprofile = ps.computeProfileCs(static_cast<int>(queryLength), counts, Neff_M, pca.values.cs(), pcb.values.cs());
//            for (size_t pos = 0; pos < queryLength; pos++) {
//                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
//                    this->profile[pos * Sequence::PROFILE_AA_SIZE + aa] = csprofile[(pos * (Sequence::PROFILE_AA_SIZE + 4)) + aa];
//                }
//            }
//        } else{
//            for (size_t pos = 0; pos < queryLength; pos++) {
//                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
//                    profile[pos * Sequence::PROFILE_AA_SIZE + aa] = matchWeight[pos * Sequence::PROFILE_AA_SIZE + aa];;
//                }
//            }
//        }
//
//        computeLogPSSM(subMat, pssm, profile, 8.0, queryLength, 0.0);
//        return Profile(pssm, profile, Neff_M, consensusSequence);
////    PSSMCalculator::printProfile(queryLength);
////    PSSMCalculator::printPSSM(queryLength);
//    }
// this overload is only used in TestProfileAlignment.cpp and TestPSSM.cpp
PSSMCalculator::Profile PSSMCalculator::computePSSMFromMSA(size_t setSize, size_t queryLength, const char **msaSeqs, bool wg, float scoreBias) {
#ifdef GAP_POS_SCORING
    std::vector<Matcher::result_t> dummy;
    return computePSSMFromMSA(setSize, queryLength, msaSeqs, dummy, wg, scoreBias);
}
#endif

#ifdef GAP_POS_SCORING
PSSMCalculator::Profile PSSMCalculator::computePSSMFromMSA(
    size_t setSize, size_t queryLength, const char **msaSeqs,
    const std::vector<Matcher::result_t> &alnResults,
    bool wg,
    float scoreBias
) {
#endif
    increaseSetSize(setSize);
    // Quick and dirty calculation of the weight per sequence wg[k]
    computeSequenceWeights(seqWeight, queryLength, setSize, msaSeqs);
    seqWeightTotal = 0.0; // TODO: MathUtil::NormalizeTo1 computes the same sum again, could be optimized?
    for (size_t i = 0; i < setSize; ++i) {
        seqWeightTotal += seqWeight[i];
    }
    MathUtil::NormalizeTo1(seqWeight, setSize);
    if (wg == false) {
        // compute context specific counts and Neff
        computeContextSpecificWeights(matchWeight, seqWeight, Neff_M, queryLength, setSize, msaSeqs);
    } else {
        // compute matchWeight based on sequence weight
        computeMatchWeights(matchWeight, seqWeight, setSize, queryLength, msaSeqs);
        // compute NEFF_M
        computeNeff_M(matchWeight, seqWeight, Neff_M, queryLength, setSize, msaSeqs);
    }
    // compute consensus sequence
    computeConsensusSequence(consensusSequence, matchWeight, queryLength, subMat->pBack, subMat->num2aa);
    if(pcmode == Parameters::PCMODE_SUBSTITUTION_SCORE && pca.values.normal() > 0.0){
        // add pseudocounts (compute the scalar product between matchWeight and substitution matrix with pseudo counts)
        preparePseudoCounts(matchWeight, pseudocountsWeight, Sequence::PROFILE_AA_SIZE, queryLength, (const float **) subMat->subMatrixPseudoCounts);
        //    SubstitutionMatrix::print(subMat->subMatrixPseudoCounts, subMat->num2aa, 20 );
        computePseudoCounts(profile, matchWeight, pseudocountsWeight, Sequence::PROFILE_AA_SIZE, Neff_M, queryLength, pca.values.normal(), pcb.values.normal());
    }else if (pcmode == Parameters::PCMODE_CONTEXT_SPECIFIC && pca.values.cs() > 0.0){
        fillCounteProfile(counts, matchWeight, Neff_M, queryLength);
        float * csprofile = ps->computeProfileCs(static_cast<int>(queryLength), counts, Neff_M, pca.values.cs(), pcb.values.cs());
        for (size_t pos = 0; pos < queryLength; pos++) {
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
                this->profile[pos * Sequence::PROFILE_AA_SIZE + aa] = csprofile[(pos * (Sequence::PROFILE_AA_SIZE + 4)) + aa];
            }
        }
    } else{
        for (size_t pos = 0; pos < queryLength; pos++) {
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
                profile[pos * Sequence::PROFILE_AA_SIZE + aa] = matchWeight[pos * Sequence::PROFILE_AA_SIZE + aa];;
            }
        }
    }

//    PSSMCalculator::printProfile(queryLength);
//    PSSMCalculator::printPSSM(queryLength);

    // create final Matrix
    computeLogPSSM(subMat, pssm, profile, 8.0, queryLength, scoreBias);
//    PSSMCalculator::printProfile(queryLength);
//    PSSMCalculator::printPSSM(queryLength);
#ifdef GAP_POS_SCORING
    computeGapPenalties(queryLength, setSize, msaSeqs, alnResults);
    return Profile(pssm, profile, Neff_M, gDel, gIns, consensusSequence);
#else
    return Profile(pssm, profile, Neff_M, consensusSequence);
#endif
}

void PSSMCalculator::printProfile(size_t queryLength) {
    printf("Pos");
    for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
        printf(" %6c", subMat->num2aa[aa]);
    }
    for (size_t i = 0; i < queryLength; i++) {
        printf("%3zu", i);
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
            printf(" %.4f", profile[i * Sequence::PROFILE_AA_SIZE + aa]);
        }
#ifdef GAP_POS_SCORING
        printf(" %7d %7d %7d\n", gDel[i] & 0xF, gDel[i] >> 4, gIns[i]);
#endif
    }
}

void PSSMCalculator::printPSSM(size_t queryLength){
    printf("Pos ");
    for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
        printf("%3c ", subMat->num2aa[aa]);
    }
    printf("\n");
    for(size_t i = 0; i <  queryLength; i++) {
        printf("%3zu ", i);
        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++){
//            char pssmVal = (pssm[i * Sequence::PROFILE_AA_SIZE + aa] == -128) ? 0 : pssm[i * Sequence::PROFILE_AA_SIZE + aa]  ;
            char pssmVal = pssm[i * Sequence::PROFILE_AA_SIZE + aa];

            printf("%3d ",  pssmVal);
        }
        printf("\n");
    }
}

void PSSMCalculator::computeLogPSSM(BaseMatrix *subMat, char *pssm, const float *profile, float bitFactor, size_t queryLength, float scoreBias) {
    for(size_t pos = 0; pos < queryLength; pos++) {
        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
            const float aaProb = profile[pos * Sequence::PROFILE_AA_SIZE + aa];
            const unsigned int idx = pos * Sequence::PROFILE_AA_SIZE + aa;
            float logProb = MathUtil::flog2(aaProb / subMat->pBack[aa]);
            float pssmVal = bitFactor * logProb + bitFactor * scoreBias;
            pssmVal = static_cast<char>((pssmVal < 0.0) ? pssmVal - 0.5 : pssmVal + 0.5);
            float truncPssmVal =  std::min(pssmVal, 127.0f);
            truncPssmVal       =  std::max(-128.0f, truncPssmVal);
            pssm[idx] = truncPssmVal;
        }
    }
}

void PSSMCalculator::preparePseudoCounts(float *frequency, float *frequency_with_pseudocounts, size_t entrySize,
                                         size_t queryLength, float const ** R) {
    for (size_t pos = 0; pos < queryLength; pos++) {
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
            frequency_with_pseudocounts[pos * entrySize + aa] = ScalarProd20(R[aa],
                                                                             &frequency[pos * entrySize]);
        }
    }
}
void PSSMCalculator::computeNeff_M(float *frequency, float *seqWeight, float *Neff_M,
                                   size_t queryLength, size_t setSize, char const **msaSeqs) {
    float Neff_HMM = 0.0f;
    for (size_t pos = 0; pos < queryLength; pos++) {
        float sum = 0.0f;
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa){
            float freq_pos_aa = frequency[pos * Sequence::PROFILE_AA_SIZE + aa];
            if (freq_pos_aa > 1E-10) {
                sum -= freq_pos_aa * MathUtil::flog2(freq_pos_aa);
            }
        }
        Neff_HMM += MathUtil::fpow2(sum);
    }
    Neff_HMM /= queryLength;
    float Nlim = fmax(10.0, Neff_HMM + 1.0);    // limiting Neff
    float scale = MathUtil::flog2((Nlim - Neff_HMM) / (Nlim - 1.0));  // for calculating Neff for those seqs with inserts at specific pos
    for (size_t pos = 0; pos < queryLength; pos++) {
        float w_M = -1.0 / setSize;
        for (size_t k = 0; k < setSize; ++k){
            if ( msaSeqs[k][pos] != MultipleAlignment::GAP) {
                w_M += seqWeight[k];
            }
        }
        Neff_M[pos] = (w_M < 0) ? 1.0 : Nlim - (Nlim - 1.0) * MathUtil::fpow2(scale * w_M);
//        fprintf(stderr,"M  i=%3i  ncol=---  Neff_M=%5.2f  Nlim=%5.2f  w_M=%5.3f  Neff_M=%5.2f\n",pos,Neff_HMM,Nlim,w_M,Neff_M[pos]);
    }
}

void PSSMCalculator::computeSequenceWeights(float *seqWeight, size_t queryLength,
                                            size_t setSize, const char **msaSeqs) {
    unsigned int *number_res = new unsigned int[setSize];
    // initialized wg[k] with tiny pseudo counts
    std::fill(seqWeight, seqWeight + setSize,  1e-6);
    // count number of residues per sequence
    for (size_t k = 0; k < setSize; ++k)
    {
        unsigned int nr = 0;
        for (size_t pos = 0; pos < queryLength; pos++) {
            if (msaSeqs[k][pos] != MultipleAlignment::GAP) {
                nr++;
            }
        }
        number_res[k] = nr;
    }
    for (size_t pos = 0; pos < queryLength; pos++) {
        int nl[ Sequence::PROFILE_AA_SIZE ];  //nl[a] = number of seq's with amino acid a at position l
        //number of different amino acids (ignore X)
        std::fill(nl, nl + Sequence::PROFILE_AA_SIZE,  0);
        for (size_t k = 0; k < setSize; ++k) {
            if (msaSeqs[k][pos] != MultipleAlignment::GAP) {
                const unsigned int aa_pos = msaSeqs[k][pos];
                if(aa_pos < Sequence::PROFILE_AA_SIZE){
                    nl[aa_pos]++;
                }
            }
        }
        //count distinct amino acids (ignore X)
        int distinct_aa_count = 0;
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
            if (nl[aa]) {
                ++distinct_aa_count;
            }
        }
//        if(distinct_aa_count == 0){
//            Debug(Debug::ERROR) << "Error in computeSequenceWeights. Distinct amino acid count is 0.\n";
//            EXIT(EXIT_FAILURE);
//        }
        // Compute sequence Weight
        // "Position-based Sequence Weights", Henikoff (1994)
        for (size_t k = 0; k < setSize; ++k) {
            if (msaSeqs[k][pos] != MultipleAlignment::GAP) {
                if (distinct_aa_count != 0) {
                    const unsigned int aa_pos = msaSeqs[k][pos];
//                    std::cout << "k="<< k << "\t";
                    if(aa_pos < Sequence::PROFILE_AA_SIZE){ // Treat score of X with other amino acid as 0.0
                        seqWeight[k] += 1.0f / (float(nl[aa_pos]) * float(distinct_aa_count) * (float(number_res[k]) + 30.0f));
//                        std::cout << number_res[k] << "\t" << distinct_aa_count << "\t" <<  nl[aa_pos] << "\t";
                    }
                }
                // ensure that each residue of a short sequence contributes as much as a residue of a long sequence:
                // contribution is proportional to one over sequence length nres[k] plus 30.
            }
        }
    }
//    std::cout << setSize << std::endl;
//    std::cout << " Seq. Weight: " << std::endl;
//    for (size_t k = 0; k < setSize; ++k) {
//        std::cout << " k="<< k << "\t" << seqWeight[k] << std::endl;
//    }
    delete [] number_res;
}

void PSSMCalculator::computePseudoCounts(float *profile, float *frequency,
                                         float *frequency_with_pseudocounts, size_t entrySize,
                                         float * Neff_M, size_t queryLength,
                                         float pca, float pcb) {
    for (size_t pos = 0; pos < queryLength; pos++) {
        float tau = fmin(1.0, pca / (1.0 + Neff_M[pos] / pcb));
        //float tau = fmin(1.0, pca * (1.0 + pcb)/ (Neff_M[pos] + pcb));

        //std::cout<< "Tau: "<< tau << ", Neff: " << Neff_M[pos] <<std::endl;
//        printf("%.6f\n", tau);

        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
            // compute proportion of pseudo counts and signal
            float pseudoCounts    = tau * frequency_with_pseudocounts[pos * entrySize + aa];
            float frequencySignal = (1.0 - tau) * frequency[pos * entrySize + aa];
            profile[pos * entrySize + aa] = frequencySignal + pseudoCounts;
//            printf("%f %f %f %f\n", tau, frequencySignal, pseudoCounts,  profile[pos * Sequence::PROFILE_AA_SIZE + aa]);
        }
    }
}

void PSSMCalculator::computeMatchWeights(float * matchWeight, float * seqWeight, size_t setSize, size_t queryLength, const char **msaSeqs) {
    for (size_t pos = 0; pos < queryLength; pos++) {
        memset(matchWeight + pos * Sequence::PROFILE_AA_SIZE, 0,
               Sequence::PROFILE_AA_SIZE * sizeof(float));
        for (size_t k = 0; k < setSize; ++k){
            if(msaSeqs[k][pos] != MultipleAlignment::GAP){
                unsigned int aa_pos = msaSeqs[k][pos];
                if(aa_pos < Sequence::PROFILE_AA_SIZE) { // Treat score of X with other amino acid as 0.0
                    matchWeight[pos * Sequence::PROFILE_AA_SIZE + aa_pos] += seqWeight[k];
                }
            }
        }
        MathUtil::NormalizeTo1(&matchWeight[pos * Sequence::PROFILE_AA_SIZE], Sequence::PROFILE_AA_SIZE, subMat->pBack);
    }
}

void PSSMCalculator::computeContextSpecificWeights(float * matchWeight, float *wg, float * Neff_M, size_t queryLength, size_t setSize,
                                                   const char **X) {
    //For weighting: include only columns into subalignment i that have a max fraction of seqs with endgap
    const float MAXENDGAPFRAC=0.1;
    const int NCOLMIN=20;   //min number of cols in subalignment for calculating pos-specific weights w[k][i]
    const int ENDGAP=22;    //Important to distinguish because end gaps do not contribute to tansition counts

    int nseqi = 0;
//    unsigned int NAA_VECSIZE = ((MultipleAlignment::NAA+ 3 + VECSIZE_INT - 1) / VECSIZE_INT) * VECSIZE_INT; // round NAA+3 up to next multiple of VECSIZE_INT
//    for (size_t j = 0; j < queryLength; j++) {
//        memset(n[j], 0, NAA_VECSIZE * sizeof(int));
//        memset(w_contrib[j], 0, NAA_VECSIZE * sizeof(int));
//    }
    unsigned int NAA_ALIGNSIZE = ((((MultipleAlignment::NAA + 3) + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT) * sizeof(float);
    NAA_ALIGNSIZE = ((NAA_ALIGNSIZE + ALIGN_FLOAT - 1) / ALIGN_FLOAT) * ALIGN_FLOAT;
    memset(n_backing, 0, NAA_ALIGNSIZE * queryLength);
    memset(w_contrib_backing, 0, NAA_ALIGNSIZE * queryLength);
    // insert endgaps
    for (size_t k = 0; k < setSize; ++k) {
        for (size_t i = 0; i < queryLength && X[k][i] == MultipleAlignment::GAP; ++i)
            ((char**)X)[k][i] = ENDGAP;
        for (int i = queryLength - 1; i >= 0 && X[k][i] == MultipleAlignment::GAP; i--)
            ((char**)X)[k][i] = ENDGAP;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Main loop through alignment columns
    for (size_t i = 0; i < queryLength; i++)  // Calculate wi[k] at position i as well as Neff[i]
    {
        bool change = false;
        // Check all sequences k and update n[j][a] and ri[j] if necessary
        for (size_t k = 0; k < setSize; ++k) {
            // Update amino acid and GAP / ENDGAP counts for sequences with AA in i-1 and GAP/ENDGAP in i or vice versa
//            printf("%d %d %d\n", k, i, (int) X[k][i - 1]);
            if ((i == 0  && X[k][i] < MultipleAlignment::ANY) ||
                (i != 0  && X[k][i - 1] >= MultipleAlignment::ANY && X[k][i] < MultipleAlignment::ANY)) {  // ... if sequence k was NOT included in i-1 and has to be included for column i
                change = true;
                nseqi++;
                for (size_t j = 0; j < queryLength; ++j){
                    n[j][(int) X[k][j]]++;
                }
            } else if ( i != 0 && X[k][i - 1] < MultipleAlignment::ANY && X[k][i] >= MultipleAlignment::ANY) {  // ... if sequence k WAS included in i-1 and has to be thrown out for column i
                change = true;
                nseqi--;
                for (size_t j = 0; j < queryLength; ++j)
                    n[j][(int) X[k][j]]--;
            }

        }  //end for (k)
        nseqs[i] = nseqi;

//        printf("%d\n", nseqi);
        // Only if subalignment changed we need to update weights wi[k] and Neff[i]
        if (change) {

            // We gained a factor ~8.0 for the following computation of weights
            // and profile by exchanging the inner two loops (j, k => k, j)
            // and precomputing the weight contributions w_contrib[j][a].
            // M. Steinegger and J. Soeding (29 July 2014)

            // Initialize weights and numbers of residues for subalignment i
            int ncol = 0;
            for (size_t k = 0; k < setSize; ++k)
                wi[k] = 1E-8;  // for pathological alignments all wi[k] can get 0;

            // Find min and max borders between which > fraction MAXENDGAPFRAC of sequences in subalignment contain an aa
            int jmin;
            int jmax;
            for (jmin = 0; jmin < static_cast<int>(queryLength) && n[jmin][ENDGAP] > MAXENDGAPFRAC * nseqi;
                 ++jmin) {
            };
            //TODO maybe wrong jmax >= 0
            for (jmax = queryLength - 1; jmax >= 0 && n[jmax][ENDGAP] > MAXENDGAPFRAC * nseqi;
                 --jmax) {
            };
            ncol = jmax - jmin + 1;
//            printf("%d %d %d\n", ncol, jmax, jmin);

            // Check whether number of columns in subalignment is sufficient
            if (ncol < NCOLMIN) {
                // Take global weights
                for (size_t k = 0; k < setSize; ++k){
                    wi[k] = (X[k][i] < MultipleAlignment::ANY)? wg[k] : 0.0f;
                }
            } else {
                // Count number of different amino acids in column j
                for (int j = jmin; j <= jmax; ++j){
                    naa[j] = 0;
                    for (int a = 0; a < MultipleAlignment::ANY; ++a){
                        naa[j] += (n[j][a] ? 1 : 0);
                    }
                }
                // Compute the contribution of amino acid a to the weight
                //for (a = 0; a < ANY; ++a)
                //      w_contrib[j][a] = (n[j][a] > 0) ? 1.0/ float(naa[j]*n[j][a]): 0.0f;
                for (int j = jmin; j <= jmax; ++j) {
                    simd_float naa_j = simdi32_i2f(simdi32_set(naa[j]));
                    const simd_int *nj = (const simd_int *) n[j];
                    const int aa_size = (MultipleAlignment::ANY + VECSIZE_INT - 1) / VECSIZE_INT;
                    for (int a = 0; a < aa_size; ++a) {
                        simd_float nja = simdi32_i2f(simdi_load(nj + a));
                        simd_float res = simdf32_mul(nja, naa_j);
                        simd_float rcp = simdf32_rcp(res);
                        // Add one iteration Newton-Raphson to improve approximate rcp
                        // https://stackoverflow.com/questions/31555260/fast-vectorized-rsqrt-and-reciprocal-with-sse-avx-depending-on-precision
                        simd_float mul = simdf32_mul(res, simdf32_mul(rcp, rcp));
                        simdf32_store(w_contrib[j] + (a * VECSIZE_INT), simdf32_sub(simdf32_add(rcp, rcp), mul));
                    }
                    for (int a = MultipleAlignment::ANY; a < MultipleAlignment::NAA + 3; ++a)
                        w_contrib[j][a] = 0.0f;  // set non-amino acid values to 0 to avoid checking in next loop for X[k][j]<ANY
                }

                // Compute pos-specific weights wi[k]
                for (size_t k = 0; k < setSize; ++k) {
                    if (X[k][i] >= MultipleAlignment::ANY)
                        continue;
                    for (int j = jmin; j <= jmax; ++j)  // innermost, time-critical loop; O(L*setSize*L)
                        wi[k] += w_contrib[j][(int) X[k][j]];
                }
            }


            // Calculate Neff[i]
            Neff_M[i] = 0.0;

            // Allocate and reset amino acid frequencies
            for (int j = jmin; j <= jmax; ++j)
                memset(f[j], 0, MultipleAlignment::ANY * sizeof(float));

            // Update f[j][a]
            for (size_t k = 0; k < setSize; ++k) {
                if (X[k][i] >= MultipleAlignment::ANY)
                    continue;
                for (int j = jmin; j <= jmax; ++j)  // innermost loop; O(L*setSize*L)
                    f[j][(int) X[k][j]] += wi[k];
            }

            // Add contributions to Neff[i]
            for (int j = jmin; j <= jmax; ++j) {
                MathUtil::NormalizeTo1(f[j], MultipleAlignment::NAA);
                for (int a = 0; a < 20; ++a)
                    if (f[j][a] > 1E-10)
                        Neff_M[i] -= f[j][a]
                                     * MathUtil::flog2(f[j][a]);
            }

            if (ncol > 0)
                Neff_M[i] = MathUtil::fpow2(Neff_M[i] / ncol);
            else
                Neff_M[i] = 1.0;

//            printf("%f\n", Neff_M[i]);
//
        }
        else  //no update was necessary; copy values for i-1
        {
            if(i==0){
                Neff_M[i] = 0.0f;
            }else{
                Neff_M[i] = Neff_M[i - 1];
            }
        }

        // Calculate amino acid frequencies q->f[i][a] from weights wi[k]
        for (int a = 0; a < 20; ++a)
            matchWeight[i * Sequence::PROFILE_AA_SIZE + a] = 0.0;
        for (size_t k = 0; k < setSize; ++k)
            matchWeight[i * Sequence::PROFILE_AA_SIZE + (int) X[k][i]] += wi[k];
        MathUtil::NormalizeTo1((matchWeight+ i * Sequence::PROFILE_AA_SIZE), MultipleAlignment::NAA, subMat->pBack);
    }
    // remove end gaps
    for (size_t k = 0; k < setSize; ++k) {
        for (size_t i = 0; i < queryLength && X[k][i] == ENDGAP; ++i)
            ((char**)X)[k][i] = MultipleAlignment::GAP;
        for (int i = queryLength - 1; i >= 0 && X[k][i] == ENDGAP; i--)
            ((char**)X)[k][i] = MultipleAlignment::GAP;
    }
}

#ifdef GAP_POS_SCORING
void PSSMCalculator::computeGapPenalties(size_t queryLength, size_t setSize, const char **msaSeqs, const std::vector<Matcher::result_t> &alnResults) {
    gapWeightsIns.clear();
    const float pseudoCounts = gapPseudoCount / seqWeightTotal;
    const float gapWeightStart = pseudoCounts * MathUtil::fpow2(-gapOpen);
    const float gapG = 0.8; // TODO: make parameter or move somewhere else?
    // compute penalties for insertions
    gapWeightsIns.resize(queryLength, gapWeightStart);
    for (size_t i = 0; i < alnResults.size(); ++i) {
        size_t queryPos = alnResults[i].qStartPos; // current position in the reference sequence
        for (size_t k = 0; k < alnResults[i].backtrace.length(); ++k) {
            char c = alnResults[i].backtrace[k];
            if (c == 'M' || c == 'I') {
                ++queryPos;
            } else if (c == 'D') {
                if (k > 0 && alnResults[i].backtrace[k - 1] == 'M') {
                    // seqWeight includes the query, alnResults does not
                    gapWeightsIns[queryPos] += seqWeight[i + 1];
                }
            }
        }
    }
    // compute penalties for deletions
    // we need the seqWeigthSum of two consecutive columns, precalculate for the first column
    float seqWeightSumPrev = 0.0;
    for (size_t i = 0; i < setSize; ++i) {
        if (msaSeqs[i][0] != MultipleAlignment::GAP) {
            seqWeightSumPrev += seqWeight[i];
        }
    }
    size_t k = queryLength - 1;
    gDel[0] = 0;
    gIns[k] = 0;
    for (size_t pos = 1; pos < queryLength; ++pos) {
        float gapWeightDelOpen = gapWeightStart;
        float gapWeightDelClose = gapWeightStart;
        float seqWeightSum = 0.0;
        for (size_t i = 0; i < setSize; ++i) {
            if (msaSeqs[i][pos] == MultipleAlignment::GAP) {
                if (msaSeqs[i][pos - 1] != MultipleAlignment::GAP) {
                    gapWeightDelOpen += seqWeight[i];
                }
            } else {
                seqWeightSum += seqWeight[i];
                if (msaSeqs[i][pos - 1] == MultipleAlignment::GAP) {
                    gapWeightDelClose += seqWeight[i];
                }
            }
        }
        float gDelOpenFwd = 0.5 * -gapG * MathUtil::flog2(gapWeightDelOpen / (seqWeightSumPrev + pseudoCounts)) + 0.5;
        float gDelCloseFwd = 0.5 * -gapG * MathUtil::flog2(gapWeightDelClose / (1 - seqWeightSumPrev + pseudoCounts)) + 0.5;
        // TODO: warning if any of the penalties is > 15 (overflow condition)
        gDel[pos] = static_cast<uint8_t>(gDelOpenFwd) | (static_cast<uint8_t>(gDelCloseFwd) << 4);
        gIns[pos - 1] = -gapG * MathUtil::flog2(gapWeightsIns[pos - 1] / (seqWeightSumPrev + pseudoCounts));
        seqWeightSumPrev = seqWeightSum;
        --k;
    }
}
#endif


void PSSMCalculator::computeConsensusSequence(unsigned char * consensusSeq, float *frequency, size_t queryLength, double *pBack, char *num2aa) {
    // TODO: consensusSeq std::vect?
    for (size_t pos = 0; pos < queryLength; pos++) {
        float maxw = 1E-8;
        int maxa = MultipleAlignment::ANY;
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
            float prob = frequency[pos * Sequence::PROFILE_AA_SIZE + aa];
            if (prob - pBack[aa] > maxw) {
                maxw = prob - pBack[aa];
                maxa = aa;
            }
        }
        consensusSeq[pos]=num2aa[maxa];
    }
}

void PSSMCalculator::Profile::toBuffer(Sequence &centerSequence, BaseMatrix& subMat, std::string &result) {
    toBuffer(centerSequence.numSequence, centerSequence.L, subMat, result);
}

void PSSMCalculator::Profile::toBuffer(const unsigned char* centerSequence, size_t centerSeqLen, BaseMatrix& subMat, std::string& result) {
    for (size_t pos = 0; pos < centerSeqLen; pos++) {
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
            result.push_back(pssm[pos * Sequence::PROFILE_AA_SIZE + aa]);
        }
        result.push_back(static_cast<unsigned char>(centerSequence[pos]));
        result.push_back(static_cast<unsigned char>(subMat.aa2num[static_cast<int>(consensus[pos])]));
        result.push_back(static_cast<unsigned char>(MathUtil::convertNeffToChar(neffM[pos])));
#ifdef GAP_POS_SCORING
        result.push_back(gDel[pos]);
        result.push_back(gIns[pos]);
#else
        result.push_back(static_cast<unsigned char>(0));
        result.push_back(static_cast<unsigned char>(0));
#endif
    }
}

void PSSMCalculator::increaseSetSize(size_t newSetSize) {
    if (newSetSize > maxSetSize) {
        maxSetSize = newSetSize * 1.5;
        seqWeight = (float*)realloc(seqWeight, maxSetSize * sizeof(float));
        wi = (float*)realloc(wi, maxSetSize * sizeof(float));
    }
}

void PSSMCalculator::fillCounteProfile(float *counts, float *matchWeight, float *Neff_M, size_t queryLength) {
    for (size_t pos = 0; pos < queryLength; ++pos) {
        //csProfile->neff[i] = Neff_M[i];
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++){
            counts[(pos * (Sequence::PROFILE_AA_SIZE + 4)) + aa] = matchWeight[pos * Sequence::PROFILE_AA_SIZE + aa] * Neff_M[pos];
        }
    }
}
