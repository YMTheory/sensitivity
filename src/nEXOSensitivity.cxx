#include <algorithm>

#include "nEXOSensitivity.hh"
#include "TObjectTable.h"
//#include "RooFit.h"
//

ClassImp(nEXOSensitivity)

nEXOSensitivity::nEXOSensitivity(int seed, const char* treeFileName) : fExcelTree(0), fWsp(0) {
    fVerboseLevel = 10;
    std::cout << "Creating nEXOSensitivity object...\n";

    SetSeed(seed);

    if (not treeFileName) {
        std::cout << "Must give tree file to load pdfs ...\n";
        exit(1);
    }
    fTreeFileName = treeFileName;

    fWriteWsp = false;
    fWspFileName = "";

    SetBinning(270, 800, 3500, 21, 10, 640);
    //fNbinsX = 270; fNbinsY = 21; //10; //40; // 56; //65; //9; //10;
    //fXmin = 800; fYmin = 10; //250; //90; //0.; //90; //0.;
    //fXmax = 3500; fYmax = 640;//650; //650.;

    //Float_t* fYbins = new Float_t();
    //fYbins[0] = fYmin;

    //fNbinsY = 7;
    //fYbins = new Double_t[fNbinsY+1];
    //fYbins[0] = 0; fYbins[1] = 90; fYbins[2] = 256; fYbins[3] = 650;
    //fYbins[0] = 0; fYbins[1] = 90; fYbins[2] = 122; fYbins[3] = 159; fYbins[4] = 202; fYbins[5] = 256; fYbins[6] = 333; fYbins[7] = 650;
    //fYbins[0] = 0; fYbins[1] = 90; fYbins[2] = 120; fYbins[3] = 160; fYbins[4] = 200; fYbins[5] = 250; fYbins[6] = 330; fYbins[7] = 650;
    //fYbins[7] = 650;

    //Suppress RooFit messages
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    fWithEff = false; //true;	//If true include the efficiency variable in the fit
    fInterpOrder = 0; //Interpolation order of the RooHistPdfs

    //fMeanSignalEff = 0.800;//900;	//mean efficiency for bb0n taken from MC
    //fSignalEffError = 0.007;	//relative efficiency error
    fFidVol = 3740.; //345; //3740.;
    fMeanSignalEff = -1;
    fSignalEffError = -1;

    //These are the assumed uncertainties for the systematics
    fFracError = 0.059; //relative ss fraction error
    fRnRateError = 0.10; //relative rate error (for Rn222 specifically)

    fSignalName = "LXeBb0n"; //"bb0n";  //The name of the signal pdf
    fWsp = new RooWorkspace("wsp"); //Global workspace

    fBaTag = false;
    fTurnedOffGroups.clear();
    fScaleBkgds = 1.;
    fScale2nu = false;

    fRandomizeMeanCounts = false;
    fRunTruthValFit = false;
    fWithStandoff = true;

    fNcpu = 1; //number of cpus to parallelize over
    fNFitPdfs = 0; //number of fitting pdfs
    fErrorLevel = 4.29778/2.; //minuit error level, to get the 90% UL
//    fErrorLevel = 1.35; //minuit error level, to get the 90% UL
    fPrintLevel = -1; //minuit print level

    fExcelTree = 0;

    fExpectCountMethod = kRdmCV;
    fUserMeanCounts.clear();

    fSSFracImprovement = 1.;
    fRn222RateCorrection = 1.;
}

nEXOSensitivity::~nEXOSensitivity() {
    std::cout << "Killing nEXOSensitivity object...\n";

    if (fWsp)  fWsp->Delete();
}

void nEXOSensitivity::SetSeed(int seed) {
    fRandom.SetSeed(seed);
    RooRandom::randomGenerator()->SetSeed(seed);
}

void nEXOSensitivity::SetBinning(int nBinsX, double xMin, double xMax, int nBinsY, double yMin, double yMax) {
    fNbinsX = nBinsX;
    fNbinsY = nBinsY;
    fXmin = xMin;
    fYmin = yMin;
    fXmax = xMax;
    fYmax = yMax;
    std::cout << "Binning set to: " << fNbinsX << " x bins = ( " << xMin << " , " << xMax << " ); " << fNbinsY << " y bins = ( " << yMin << " , " << yMax << " )." << std::endl;
}

void nEXOSensitivity::LoadExcelTree(const char* filename, const char* treename) {
    TFile treeFile(filename, "read");
    TTree* tree = dynamic_cast<TTree*> (treeFile.Get(treename));

    fExcelTree = tree->CloneTree();
    fExcelTree->SetDirectory(0);

    treeFile.Close();
}

void nEXOSensitivity::ReadExcelTree() {
    // Set vectors of groups: group fractions and mean counts - fGroups, fGroupFractions and fGroupMeanCounts
    fComponentNames.clear();
    fGroups.clear();
    fGroupFractions.clear();
    fGroupMeanCounts.clear();

    std::map<TString, Double_t> groupCV;
    std::map<TString, Double_t> groupError;
    std::map<TString, Double_t> specActivities;

    std::map<TString, Double_t> groupRatio0p5tfwhm;
    std::map<TString, Double_t> groupRatio1p5tfwhm;
    std::map<TString, Double_t> groupRatio2p5tfwhm;
    std::map<TString, Double_t> groupRatio3p5tfwhm;
    std::map<TString, Double_t> groupRatio1tfwhm;
    std::map<TString, Double_t> groupRatio2tfwhm;
    std::map<TString, Double_t> groupRatio3tfwhm;
    std::map<TString, Double_t> groupRatioFVfwhm;

    std::map<TString, Double_t> groupRatio0p5t1sig;
    std::map<TString, Double_t> groupRatio1p5t1sig;
    std::map<TString, Double_t> groupRatio2p5t1sig;
    std::map<TString, Double_t> groupRatio3p5t1sig;
    std::map<TString, Double_t> groupRatio1t1sig;
    std::map<TString, Double_t> groupRatio2t1sig;
    std::map<TString, Double_t> groupRatio3t1sig;
    std::map<TString, Double_t> groupRatioFV1sig;

    std::map<TString, Double_t> groupRatio0p5t2sig;
    std::map<TString, Double_t> groupRatio1p5t2sig;
    std::map<TString, Double_t> groupRatio2p5t2sig;
    std::map<TString, Double_t> groupRatio3p5t2sig;
    std::map<TString, Double_t> groupRatio1t2sig;
    std::map<TString, Double_t> groupRatio2t2sig;
    std::map<TString, Double_t> groupRatio3t2sig;
    std::map<TString, Double_t> groupRatioFV2sig;

    std::map<TString, Double_t> groupRatio0p5t3sig;
    std::map<TString, Double_t> groupRatio1p5t3sig;
    std::map<TString, Double_t> groupRatio2p5t3sig;
    std::map<TString, Double_t> groupRatio3p5t3sig;
    std::map<TString, Double_t> groupRatio1t3sig;
    std::map<TString, Double_t> groupRatio2t3sig;
    std::map<TString, Double_t> groupRatio3t3sig;
    std::map<TString, Double_t> groupRatioFV3sig;

    ExcelTableValues* table = new ExcelTableValues();
    fExcelTree->SetBranchAddress("table", &table);
    for (int i = 0; i < fExcelTree->GetEntriesFast(); i++) {
        fExcelTree->GetEntry(i);

        //std::cout << "Working on " << table->fPdf << std::endl;
        if (fVerboseLevel > 0)
            table->Print();

        fComponentNames.push_back(table->fPdf.Data());

        if (fWithEff && table->fGroup == fSignalName) // set efficiency
        {
            Double_t factor = (table->fSpecActivError > 0) ? table->fActivError / table->fSpecActivError : table->fActivCV / table->fSpecActivCV;
            factor *= 1e3; // from mBq to Bq in table
            fMeanSignalEff = table->GetEfficiency() * factor / fFidVol; // correct for fiducial volume;
            fSignalEffError = table->GetEfficiencyError() * factor / fFidVol; // correct for fiducial volume;

            fMeanSignalEff *= fNormHistoBins.at(Form("%s_effcor", fSignalName.Data()));
            fSignalEffError *= fNormHistoBins.at(Form("%s_effcor", fSignalName.Data()));
        }

        if (fGroups.count(table->fGroup.Data()) <= 0)
            fGroups.insert(std::make_pair(table->fGroup.Data(), std::vector<TString>()));
        fGroups[table->fGroup.Data()].push_back(table->fPdf.Data());

        // Fill group fractions and group mean counts
        for (size_t s = 0; s < table->fSuffixes.size(); s++) {

            // First - fill group fraction
            TString groupFractionName(Form("%s_%s", table->fGroup.Data(), table->fSuffixes[s].Data()));

            if (fGroupFractions.count(groupFractionName.Data()) <= 0)
                fGroupFractions.insert(std::make_pair(groupFractionName.Data(), std::vector<Double_t>()));

            if (fGroupSeparation.count(groupFractionName.Data()) <= 0)
                fGroupSeparation.insert(std::make_pair(groupFractionName.Data(), std::make_pair(0., 0.)));


            Double_t counts = 0.;
            switch (fExpectCountMethod) {
                case kUL:
                case kPosUL:
                    counts = table->fCountsUL[s];
                    break;

                case kPosCV:
                    counts = (table->fCountsCV[s] > 0) ? table->fCountsCV[s] : 0;
                    break;

                case kRdmCV:
                    if (groupFractionName == "LXeBb2n_SS" or groupFractionName == "LXeBb2n_MS" or groupFractionName == "LXeBb0n_SS" or groupFractionName == "LXeBb0n_MS")
                        counts = table->fCountsCV[s];
                    else {
                        Double_t specActivity = table->fSpecActivCV;
                        TString activID = Form("%s_%s", table->fActivID.Data(), table->fIsotope.Data());
                        if (specActivities.count(activID.Data()) > 0) {
                            specActivity = specActivities[activID.Data()];
                        } else {
                            if (fRandomizeMeanCounts) {
                                specActivity = -1;
                                int nDraws = 1000;
                                while (specActivity < 0 && nDraws > 0) {
                                    specActivity = fRandom.Gaus(table->fSpecActivCV, table->fSpecActivError);
                                    //std::cout << "Seed " << fRandom.GetSeed() << " randomized activity " << specActivity << " from " << table->fSpecActivCV << " " << table->fSpecActivError << std::endl;
                                    //double unif = fRandom.Uniform();
                                    //double alternative = GenTruncGaus(table->fSpecActivCV,table->fSpecActivError,unif);
                                    //std::cout << "Uniform " << unif << " alternative " << alternative << std::endl;
                                    //std::cin.get();
                                    nDraws--;
                                }
                                if (specActivity < 0 && nDraws == 0)
                                    specActivity = 1e-16;
                            }
                            specActivities[activID.Data()] = specActivity;
                        }

                        Double_t factor = (table->fSpecActivError > 0) ? table->fActivError / table->fSpecActivError : table->fActivCV / table->fSpecActivCV;
                        Double_t activity = factor * specActivity;

                        if (fVerboseLevel > 0)
                            std::cout << Form("Using activity ID: %s , spec = %g +- %g , full = %g +- %g , eval = %g (from %g and %g)", table->fActivID.Data(), table->fSpecActivCV, table->fSpecActivError, table->fActivCV, table->fActivError, activity, factor, specActivity) << std::endl;
                        counts = EvalCounts(table->fHitEffK[s] / table->fHitEffN[s], activity, 1., table->fHalflife);
                    }
                    break;
                default:
                    counts = 0.;
            }
            counts = (counts > 0) ? counts : 0.;
            TString suffix = table->fSuffixes[s];
            suffix.ToLower();
            TString h_id_name = Form("%s_%s", table->fPdf.Data(), suffix.Data());
            if (fVerboseLevel > 0)
                std::cout << "Reading bin fraction with name " << h_id_name << std::endl;
            counts *= fNormHistoBins.at(h_id_name);
            fGroupFractions[groupFractionName.Data()].push_back(counts);

            // Second - fill group mean counts
            if (fGroupMeanCounts.count(groupFractionName.Data()) <= 0) {
                fGroupMeanCounts.insert(std::make_pair(groupFractionName.Data(), 0));
                groupCV.insert(std::make_pair(groupFractionName.Data(), 0));
                groupError.insert(std::make_pair(groupFractionName.Data(), 0));

                groupRatioFVfwhm.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio3tfwhm.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio2tfwhm.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio1tfwhm.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio3p5tfwhm.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio2p5tfwhm.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio1p5tfwhm.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio0p5tfwhm.insert(std::make_pair(groupFractionName.Data(), 0.));

                groupRatioFV1sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio3t1sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio2t1sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio1t1sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio3p5t1sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio2p5t1sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio1p5t1sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio0p5t1sig.insert(std::make_pair(groupFractionName.Data(), 0.));

                groupRatioFV2sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio3t2sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio2t2sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio1t2sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio3p5t2sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio2p5t2sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio1p5t2sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio0p5t2sig.insert(std::make_pair(groupFractionName.Data(), 0.));

                groupRatioFV3sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio3t3sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio2t3sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio1t3sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio3p5t3sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio2p5t3sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio1p5t3sig.insert(std::make_pair(groupFractionName.Data(), 0.));
                groupRatio0p5t3sig.insert(std::make_pair(groupFractionName.Data(), 0.));
            }

            Double_t mean = 0.;
            Double_t error = pow(table->fCountsError[s], 2);

            Double_t ratioFVfwhm(table->fRatioFWHMfv[s]), ratio3tfwhm(table->fRatioFWHM3t[s]), ratio2tfwhm(table->fRatioFWHM2t[s]), ratio1tfwhm(table->fRatioFWHM1t[s]), ratio3p5tfwhm(table->fRatioFWHM3p5t[s]), ratio2p5tfwhm(table->fRatioFWHM2p5t[s]), ratio1p5tfwhm(table->fRatioFWHM1p5t[s]), ratio0p5tfwhm(table->fRatioFWHM0p5t[s]);
            Double_t ratioFV1s(table->fRatio1sfv[s]), ratio3t1s(table->fRatio1s3t[s]), ratio2t1s(table->fRatio1s2t[s]), ratio1t1s(table->fRatio1s1t[s]), ratio3p5t1s(table->fRatio1s3p5t[s]), ratio2p5t1s(table->fRatio1s2p5t[s]), ratio1p5t1s(table->fRatio1s1p5t[s]), ratio0p5t1s(table->fRatio1s0p5t[s]);
            Double_t ratioFV2s(table->fRatio2sfv[s]), ratio3t2s(table->fRatio2s3t[s]), ratio2t2s(table->fRatio2s2t[s]), ratio1t2s(table->fRatio2s1t[s]), ratio3p5t2s(table->fRatio2s3p5t[s]), ratio2p5t2s(table->fRatio2s2p5t[s]), ratio1p5t2s(table->fRatio2s1p5t[s]), ratio0p5t2s(table->fRatio2s0p5t[s]);
            Double_t ratioFV3s(table->fRatio3sfv[s]), ratio3t3s(table->fRatio3s3t[s]), ratio2t3s(table->fRatio3s2t[s]), ratio1t3s(table->fRatio3s1t[s]), ratio3p5t3s(table->fRatio3s3p5t[s]), ratio2p5t3s(table->fRatio3s2p5t[s]), ratio1p5t3s(table->fRatio3s1p5t[s]), ratio0p5t3s(table->fRatio3s0p5t[s]);

            switch (fExpectCountMethod) {
                case kUL:
                    mean = table->fCountsCV[s];
                    break;
                case kPosUL:
                    mean = (table->fCountsCV[s] > 0) ? table->fCountsCV[s] : 0.;
                    break;
                case kPosCV:
                case kRdmCV:
                    mean = counts;
                    break;
                default:
                    mean = 0.;
            }

            groupCV[groupFractionName.Data()] += mean;
            groupError[groupFractionName.Data()] += error;

            groupRatioFVfwhm[groupFractionName.Data()] += counts*ratioFVfwhm;
            groupRatio3tfwhm[groupFractionName.Data()] += counts*ratio3tfwhm;
            groupRatio2tfwhm[groupFractionName.Data()] += counts*ratio2tfwhm;
            groupRatio1tfwhm[groupFractionName.Data()] += counts*ratio1tfwhm;
            groupRatio3p5tfwhm[groupFractionName.Data()] += counts*ratio3p5tfwhm;
            groupRatio2p5tfwhm[groupFractionName.Data()] += counts*ratio2p5tfwhm;
            groupRatio1p5tfwhm[groupFractionName.Data()] += counts*ratio1p5tfwhm;
            groupRatio0p5tfwhm[groupFractionName.Data()] += counts*ratio0p5tfwhm;
            groupRatioFV1sig[groupFractionName.Data()] += counts*ratioFV1s;
            groupRatio3t1sig[groupFractionName.Data()] += counts*ratio3t1s;
            groupRatio2t1sig[groupFractionName.Data()] += counts*ratio2t1s;
            groupRatio1t1sig[groupFractionName.Data()] += counts*ratio1t1s;
            groupRatio3p5t1sig[groupFractionName.Data()] += counts*ratio3p5t1s;
            groupRatio2p5t1sig[groupFractionName.Data()] += counts*ratio2p5t1s;
            groupRatio1p5t1sig[groupFractionName.Data()] += counts*ratio1p5t1s;
            groupRatio0p5t1sig[groupFractionName.Data()] += counts*ratio0p5t1s;
            groupRatioFV2sig[groupFractionName.Data()] += counts*ratioFV2s;
            groupRatio3t2sig[groupFractionName.Data()] += counts*ratio3t2s;
            groupRatio2t2sig[groupFractionName.Data()] += counts*ratio2t2s;
            groupRatio1t2sig[groupFractionName.Data()] += counts*ratio1t2s;
            groupRatio3p5t2sig[groupFractionName.Data()] += counts*ratio3p5t2s;
            groupRatio2p5t2sig[groupFractionName.Data()] += counts*ratio2p5t2s;
            groupRatio1p5t2sig[groupFractionName.Data()] += counts*ratio1p5t2s;
            groupRatio0p5t2sig[groupFractionName.Data()] += counts*ratio0p5t2s;
            groupRatioFV3sig[groupFractionName.Data()] += counts*ratioFV3s;
            groupRatio3t3sig[groupFractionName.Data()] += counts*ratio3t3s;
            groupRatio2t3sig[groupFractionName.Data()] += counts*ratio2t3s;
            groupRatio1t3sig[groupFractionName.Data()] += counts*ratio1t3s;
            groupRatio3p5t3sig[groupFractionName.Data()] += counts*ratio3p5t3s;
            groupRatio2p5t3sig[groupFractionName.Data()] += counts*ratio2p5t3s;
            groupRatio1p5t3sig[groupFractionName.Data()] += counts*ratio1p5t3s;
            groupRatio0p5t3sig[groupFractionName.Data()] += counts*ratio0p5t3s;
        }
    }

    for (std::map<TString, Double_t>::iterator group = fGroupMeanCounts.begin(); group != fGroupMeanCounts.end(); group++) {
        TString groupName = group->first;
        Double_t poscv = (groupCV[groupName.Data()] > 0) ? groupCV[groupName.Data()] : 0;
        Double_t error = sqrt(groupError[groupName.Data()]);

        switch (fExpectCountMethod) {
            case kUL:
            case kPosUL:
                fGroupMeanCounts[groupName.Data()] = poscv + 1.64 * error;
                break;

            case kPosCV:
            case kRdmCV:
                fGroupMeanCounts[groupName.Data()] = poscv;
                break;

            default:
                fGroupMeanCounts[groupName.Data()] = groupCV[groupName.Data()];
        }

        Double_t sumCounts = std::accumulate(fGroupFractions[groupName.Data()].begin(), fGroupFractions[groupName.Data()].end(), 0.);

        fGroupMeanRatioFVfwhm[groupName.Data()] = sumCounts > 0 ? groupRatioFVfwhm[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio3tfwhm[groupName.Data()] = sumCounts > 0 ? groupRatio3tfwhm[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio2tfwhm[groupName.Data()] = sumCounts > 0 ? groupRatio2tfwhm[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio1tfwhm[groupName.Data()] = sumCounts > 0 ? groupRatio1tfwhm[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio3p5tfwhm[groupName.Data()] = sumCounts > 0 ? groupRatio3p5tfwhm[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio2p5tfwhm[groupName.Data()] = sumCounts > 0 ? groupRatio2p5tfwhm[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio1p5tfwhm[groupName.Data()] = sumCounts > 0 ? groupRatio1p5tfwhm[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio0p5tfwhm[groupName.Data()] = sumCounts > 0 ? groupRatio0p5tfwhm[groupName.Data()] / sumCounts : 1.;

        fGroupMeanRatioFV1sig[groupName.Data()] = sumCounts > 0 ? groupRatioFV1sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio3t1sig[groupName.Data()] = sumCounts > 0 ? groupRatio3t1sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio2t1sig[groupName.Data()] = sumCounts > 0 ? groupRatio2t1sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio1t1sig[groupName.Data()] = sumCounts > 0 ? groupRatio1t1sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio3p5t1sig[groupName.Data()] = sumCounts > 0 ? groupRatio3p5t1sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio2p5t1sig[groupName.Data()] = sumCounts > 0 ? groupRatio2p5t1sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio1p5t1sig[groupName.Data()] = sumCounts > 0 ? groupRatio1p5t1sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio0p5t1sig[groupName.Data()] = sumCounts > 0 ? groupRatio0p5t1sig[groupName.Data()] / sumCounts : 1.;

        fGroupMeanRatioFV2sig[groupName.Data()] = sumCounts > 0 ? groupRatioFV2sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio3t2sig[groupName.Data()] = sumCounts > 0 ? groupRatio3t2sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio2t2sig[groupName.Data()] = sumCounts > 0 ? groupRatio2t2sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio1t2sig[groupName.Data()] = sumCounts > 0 ? groupRatio1t2sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio3p5t2sig[groupName.Data()] = sumCounts > 0 ? groupRatio3p5t2sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio2p5t2sig[groupName.Data()] = sumCounts > 0 ? groupRatio2p5t2sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio1p5t2sig[groupName.Data()] = sumCounts > 0 ? groupRatio1p5t2sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio0p5t2sig[groupName.Data()] = sumCounts > 0 ? groupRatio0p5t2sig[groupName.Data()] / sumCounts : 1.;

        fGroupMeanRatioFV3sig[groupName.Data()] = sumCounts > 0 ? groupRatioFV3sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio3t3sig[groupName.Data()] = sumCounts > 0 ? groupRatio3t3sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio2t3sig[groupName.Data()] = sumCounts > 0 ? groupRatio2t3sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio1t3sig[groupName.Data()] = sumCounts > 0 ? groupRatio1t3sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio3p5t3sig[groupName.Data()] = sumCounts > 0 ? groupRatio3p5t3sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio2p5t3sig[groupName.Data()] = sumCounts > 0 ? groupRatio2p5t3sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio1p5t3sig[groupName.Data()] = sumCounts > 0 ? groupRatio1p5t3sig[groupName.Data()] / sumCounts : 1.;
        fGroupMeanRatio0p5t3sig[groupName.Data()] = sumCounts > 0 ? groupRatio0p5t3sig[groupName.Data()] / sumCounts : 1.;

        if (fVerboseLevel > 0)
            std::cout << Form("Group name : %s , mean counts = %g , ratio FWHM-3t = %g , ratio FWHM-2t = %g , ratio FWHM-1t = %g ", groupName.Data(), fGroupMeanCounts[groupName.Data()], fGroupMeanRatio3tfwhm[groupName.Data()], fGroupMeanRatio2t2sig[groupName.Data()], fGroupMeanRatio1t1sig[groupName.Data()]) << std::endl;
    }
    delete table;
}

void nEXOSensitivity::SetUserMeanCounts(TString groupName, Double_t value) {
    // User value should correspond to inner ROI (SS) in FWHM-3t

    TString groupNameSS = Form("%s_SS", groupName.Data());
    TString groupNameMS = Form("%s_MS", groupName.Data());

    Double_t total = fGroupMeanCounts[groupNameSS.Data()] + fGroupMeanCounts[groupNameMS.Data()];
    Double_t frac = fGroupMeanCounts[groupNameSS.Data()] / total;

    if (fVerboseLevel > 0)
        std::cout << Form("Manually setting group : %s to  %g - Originals : SS = %g | MS = %g | Total = %g | Fraction = %g", groupName.Data(), value, fGroupMeanCounts[groupNameSS.Data()], fGroupMeanCounts[groupNameMS.Data()], total, frac) << std::endl;

    Double_t ratio = -1;
    if (fGroupMeanRatio3tfwhm[groupNameSS.Data()] > 0) // first try to scale by SS ratio in of FWHM-3t
    {
        ratio = fGroupMeanRatio3tfwhm[groupNameSS.Data()];
        fGroupMeanCounts[groupNameSS.Data()] = value / ratio;
        fGroupMeanCounts[groupNameMS.Data()] = fGroupMeanCounts[groupNameSS.Data()] * (1. - frac) / frac; // preserve SS fraction
    } else if (fGroupMeanRatio3tfwhm[groupNameSS.Data()] > 0) // second try to scale by MS ratio in of FWHM-3t
    {
        ratio = fGroupMeanRatio3tfwhm[groupNameMS.Data()];
        fGroupMeanCounts[groupNameMS.Data()] = value / ratio;
        fGroupMeanCounts[groupNameSS.Data()] = fGroupMeanCounts[groupNameMS.Data()] * frac / (1. - frac); // preserve SS fraction
    } else // that means there should be no events in FWHM-3t of this group, so just set to given value
    {
        fGroupMeanCounts[groupNameSS.Data()] = value;
        fGroupMeanCounts[groupNameMS.Data()] = fGroupMeanCounts[groupNameSS.Data()] * (1. - frac) / frac; // preserve SS fraction    
    }

    total = fGroupMeanCounts[groupNameSS.Data()] + fGroupMeanCounts[groupNameMS.Data()];
    frac = fGroupMeanCounts[groupNameSS.Data()] / total;
    if (fVerboseLevel > 0)
        std::cout << Form("Manually setting group : %s to %g (ratio %g) - Finals : SS = %g | MS = %g | Total = %g | Fraction = %g", groupName.Data(), value, ratio, fGroupMeanCounts[groupNameSS.Data()], fGroupMeanCounts[groupNameMS.Data()], total, frac) << std::endl;

}

void nEXOSensitivity::SetAllGroupMeanCounts(Double_t value, TString except) {
    // Set the mean counts of all groups to value
    // Value must correspond top inner 3t scaled to full by ratio3t

    for (std::map<TString, std::vector<TString> >::iterator group = fGroups.begin(); group != fGroups.end(); group++) {
        if (group->first == except)
            continue;
        SetUserMeanCounts(group->first, value);
    }
}

void nEXOSensitivity::MakeFittingHistogramFile() {
}

void nEXOSensitivity::LoadComponentHistograms() {
    TFile* fIn = 0;

    TH2D* h_ss = 0;
    TH2D* h_ms = 0;

    ExcelTableValues* table = new ExcelTableValues();
    fExcelTree->SetBranchAddress("table", &table);
    for (int i = 0; i < fExcelTree->GetEntriesFast(); i++) {
        fExcelTree->GetEntry(i);

        //std::cout << "Working on " << table->fPdf << std::endl;
        if (fVerboseLevel)
            table->Print();

        fIn = new TFile(table->fFileName.Data()); // new TFile(Form("%s/nEXO_Histos_%s_R.root", fHistoPathIn.Data(), pdfNames[i].Data()));

        h_ss = (TH2D*) fIn->Get("h_StandoffVsEnergySS_Smear");
        h_ms = (TH2D*) fIn->Get("h_StandoffVsEnergyMS_Smear");
        h_ss->SetName(Form("h_%s_ss", table->fPdf.Data()));
        h_ss->SetTitle("SS: Stand Off Vs. Smeared Energy Histogram");
        h_ms->SetName(Form("h_%s_ms", table->fPdf.Data()));
        h_ms->SetTitle("MS: Stand Off Vs. Smeared Energy Histogram");

        //outFile->cd();

        //TH2D* hh_ss = (TH2D*)h_ss->Clone();
        //TH2D* hh_ms = (TH2D*)h_ms->Clone();

        //int nBinsXcomp = hh_ss->GetNbinsX();
        //int xRebin = (int) nBinsXcomp / fNbinsX;
        //int nBinsYcomp = hh_ss->GetNbinsY();
        //int yRebin = (int) nBinsYcomp / fNbinsY;    
        //hh_ss->Rebin2D(xRebin, yRebin);
        //hh_ms->Rebin2D(xRebin, yRebin);

        //TAxis* xAxis = h_ss->GetXaxis();
        //TAxis* yAxis = h_ss->GetYaxis();

        //double histFull = h_ss->Integral();
        //double histFWHMfv = h_ss->Integral(xAxis->FindBin(2428),xAxis->FindBin(2488),1,genHist_ss->GetNbinsY());
        //double histFWHM3t = h_ss->Integral(xAxis->FindBin(2428),xAxis->FindBin(2488),yAxis->FindBin(90),genHist_ss->GetNbinsY());
        //double histFWHM1t = h_ss->Integral(xAxis->FindBin(2428),xAxis->FindBin(2488),yAxis->FindBin(256),genHist_ss->GetNbinsY());

        TH2D* hh_ss = AdjustedBinHist(*h_ss);
        TH2D* hh_ms = AdjustedBinHist(*h_ms);

        TString h_ss_name = Form("%s_ss", table->fPdf.Data());
        TString h_ms_name = Form("%s_ms", table->fPdf.Data());

        std::cout<<"size "<<fComponentHistos.size()<<std::endl;
        fComponentHistos.insert(std::make_pair(h_ss_name, hh_ss));
        fComponentHistos.insert(std::make_pair(h_ms_name, hh_ms));

        fComponentHistos.at(h_ss_name)->SetDirectory(0);
        fComponentHistos.at(h_ms_name)->SetDirectory(0);

        if (fVerboseLevel > 0)
            std::cout << "Adding " << h_ss_name << " bin fraction = " << hh_ss->Integral() / h_ss->Integral() << std::endl;

        fNormHistoBins.insert(std::make_pair(h_ss_name, hh_ss->Integral() / h_ss->Integral()));
        fNormHistoBins.insert(std::make_pair(h_ms_name, hh_ms->Integral() / h_ms->Integral()));

        //fNormHistoFull.insert(std::make_pair(h_ss_name,histFull));
        //fRatioHistoFV.insert(std::make_pair(h_ss_name,histFWHMfv/histFull));
        //fRatioHisto3t.insert(std::make_pair(h_ss_name,histFWHM3t/histFull));
        //fRatioHisto1t.insert(std::make_pair(h_ss_name,histFWHM1t/histFull));

        if (fWithEff && table->fGroup == fSignalName) // set efficiency
        {
            fNormHistoBins.insert(std::make_pair(Form("%s_effcor", fSignalName.Data()), (hh_ss->Integral() + hh_ms->Integral()) / (h_ss->Integral() + h_ms->Integral())));
        }

        fIn->Close();
        delete fIn;
    }
    delete table;
}

TH2D* nEXOSensitivity::AdjustedBinHist(TH2D& inHist) {
    if (fXmin < inHist.GetXaxis()->GetXmin() or fXmax > inHist.GetXaxis()->GetXmax() or fYmin < inHist.GetYaxis()->GetXmin() or fYmax > inHist.GetYaxis()->GetXmax()) {
        std::cout << "Requested limits outside of given boundaries, please fix this issue. Will quit now...\n";
        exit(1);
    }
    TH2D* resHist = new TH2D(Form("%s_binned", inHist.GetName()), inHist.GetTitle(), fNbinsX, fXmin, fXmax, fNbinsY, fYmin, fYmax); //fNbinsY, fYbins); //fNbinsY, fYmin, fYmax);
    if(fbinAdjustMap.size()!=inHist.GetNbinsX()*inHist.GetNbinsY()){
        fbinAdjustMap.clear();
//        std::cout<<"Making histogram adjustment map: "<<inHist.GetNbinsX()*inHist.GetNbinsY()<<std::endl;
        //make adjustment map for fast adjustment
        for (int bx = 1; bx <= inHist.GetNbinsX(); bx++) {
            double centerx = inHist.GetXaxis()->GetBinCenter(bx);
            for (int by = 1; by <= inHist.GetNbinsY(); by++) {
                double centery = inHist.GetYaxis()->GetBinCenter(by);
                fbinAdjustMap[inHist.GetBin(bx,by)]=resHist->FindBin(centerx,centery);
            }
        }
//        std::cout<<"Made histogram adjustment map, "<<fbinAdjustMap.size()<<" entries"<<std::endl;
    }
    
    for (auto ia:fbinAdjustMap){
        double prevContent = resHist->GetBinContent(ia.second);
        double prevError = resHist->GetBinError(ia.second);
        resHist->SetBinContent(ia.second,prevContent + inHist.GetBinContent(ia.first));
        resHist->SetBinError(ia.second,sqrt(prevError*prevError + inHist.GetBinError(ia.first)*inHist.GetBinError(ia.first)));
        resHist->SetEntries(resHist->GetEntries()+inHist.GetBinContent(ia.first));
    }
//    std::cout<<"Adjusted histogram, with "<<resHist->GetEntries()<<" entries"<<std::endl;
    
    return resHist;
}

void nEXOSensitivity::MakeGroupHistograms() {
    // Make the histograms (pdfs) of the combined groups - fGroupHistos
    
//    gObjectTable->Print();
    for (std::map<TString, TH1*>::iterator groupHisto = fGroupHistos.begin(); groupHisto != fGroupHistos.end(); groupHisto++)
        groupHisto->second->Delete();
    fGroupHistos.clear();

    double tot_ss = 0.;
    double tot_ms = 0.;
    double roi = 0.;
    double roi_3t = 0.;

    fBkgdTotal = 0.;
    fBkgdFwhmFV = 0.;
    fBkgdFwhm3t = 0.;
    fBkgdFwhm2t = 0.;
    fBkgdFwhm1t = 0.;
    fBkgdFwhm3p5t = 0.;
    fBkgdFwhm2p5t = 0.;
    fBkgdFwhm1p5t = 0.;
    fBkgdFwhm0p5t = 0.;
    fBkgd1sigmaFV = 0.;
    fBkgd1sigma3t = 0.;
    fBkgd1sigma2t = 0.;
    fBkgd1sigma1t = 0.;
    fBkgd1sigma3p5t = 0.;
    fBkgd1sigma2p5t = 0.;
    fBkgd1sigma1p5t = 0.;
    fBkgd1sigma0p5t = 0.;
    fBkgd2sigmaFV = 0.;
    fBkgd2sigma3t = 0.;
    fBkgd2sigma2t = 0.;
    fBkgd2sigma1t = 0.;
    fBkgd2sigma3p5t = 0.;
    fBkgd2sigma2p5t = 0.;
    fBkgd2sigma1p5t = 0.;
    fBkgd2sigma0p5t = 0.;
    fBkgd3sigmaFV = 0.;
    fBkgd3sigma3t = 0.;
    fBkgd3sigma2t = 0.;
    fBkgd3sigma1t = 0.;
    fBkgd3sigma3p5t = 0.;
    fBkgd3sigma2p5t = 0.;
    fBkgd3sigma1p5t = 0.;
    fBkgd3sigma0p5t = 0.;

    //double roi_r3t = 0.;
    //double sum = 0.;
    for (std::map<TString, std::vector<TString> >::iterator group = fGroups.begin(); group != fGroups.end(); group++) {
        TString groupName = group->first;
        std::vector<TString>& groupComponents = group->second;

        std::vector<Double_t> groupFractionsSS;
        std::vector<Double_t> groupFractionsMS;

        for (size_t i = 0; i < groupComponents.size(); i++) {
            groupFractionsSS.push_back(fGroupFractions[Form("%s_%s", groupName.Data(), "SS")][i]);
            groupFractionsMS.push_back(fGroupFractions[Form("%s_%s", groupName.Data(), "MS")][i]);
        }
        TH2D* h_ss = MakeCombinedHisto(Form("h_%s_ss", groupName.Data()), groupComponents.size(), &groupComponents[0], &groupFractionsSS[0], "", true);
        TH2D* h_ms = MakeCombinedHisto(Form("h_%s_ms", groupName.Data()), groupComponents.size(), &groupComponents[0], &groupFractionsMS[0], "", false);

        //tot_ss += h_ss->Integral()*fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        //tot_ms += h_ms->Integral()*fGroupMeanCounts[Form("%s_%s",groupName.Data(),"MS")];
        //roi += h_ss->Integral(h_ss->GetXaxis()->FindBin(2428),h_ss->GetXaxis()->FindBin(2488),1,h_ss->GetNbinsY())*fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        //roi_3t += h_ss->Integral(h_ss->GetXaxis()->FindBin(2428),h_ss->GetXaxis()->FindBin(2488),h_ss->GetYaxis()->FindBin(90),h_ss->GetNbinsY())*fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double full_ss = fGroupMeanCounts[Form("%s_%s", groupName.Data(), "SS")];

        double fwhm_fv = fGroupMeanRatioFVfwhm[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double fwhm_3t = fGroupMeanRatio3tfwhm[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double fwhm_2t = fGroupMeanRatio2tfwhm[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double fwhm_1t = fGroupMeanRatio1tfwhm[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double fwhm_3p5t = fGroupMeanRatio3p5tfwhm[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double fwhm_2p5t = fGroupMeanRatio2p5tfwhm[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double fwhm_1p5t = fGroupMeanRatio1p5tfwhm[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double fwhm_0p5t = fGroupMeanRatio0p5tfwhm[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];

        double onesig_fv = fGroupMeanRatioFV1sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double onesig_3t = fGroupMeanRatio3t1sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double onesig_2t = fGroupMeanRatio2t1sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double onesig_1t = fGroupMeanRatio1t1sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double onesig_3p5t = fGroupMeanRatio3p5t1sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double onesig_2p5t = fGroupMeanRatio2p5t1sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double onesig_1p5t = fGroupMeanRatio1p5t1sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double onesig_0p5t = fGroupMeanRatio0p5t1sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];

        double twosig_fv = fGroupMeanRatioFV2sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double twosig_3t = fGroupMeanRatio3t2sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double twosig_2t = fGroupMeanRatio2t2sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double twosig_1t = fGroupMeanRatio1t2sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double twosig_3p5t = fGroupMeanRatio3p5t2sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double twosig_2p5t = fGroupMeanRatio2p5t2sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double twosig_1p5t = fGroupMeanRatio1p5t2sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double twosig_0p5t = fGroupMeanRatio0p5t2sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];

        double threesig_fv = fGroupMeanRatioFV3sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double threesig_3t = fGroupMeanRatio3t3sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double threesig_2t = fGroupMeanRatio2t3sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double threesig_1t = fGroupMeanRatio1t3sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double threesig_3p5t = fGroupMeanRatio3p5t3sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double threesig_2p5t = fGroupMeanRatio2p5t3sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double threesig_1p5t = fGroupMeanRatio1p5t3sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        double threesig_0p5t = fGroupMeanRatio0p5t3sig[Form("%s_%s", groupName.Data(), "SS")] * full_ss; //fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];

        fBkgdTotal += full_ss;
        fBkgdFwhmFV += fwhm_fv;
        fBkgdFwhm3t += fwhm_3t;
        fBkgdFwhm2t += fwhm_2t;
        fBkgdFwhm1t += fwhm_1t;
        fBkgdFwhm3p5t += fwhm_3p5t;
        fBkgdFwhm2p5t += fwhm_2p5t;
        fBkgdFwhm1p5t += fwhm_1p5t;
        fBkgdFwhm0p5t += fwhm_0p5t;
        fBkgd1sigmaFV += onesig_fv;
        fBkgd1sigma3t += onesig_3t;
        fBkgd1sigma2t += onesig_2t;
        fBkgd1sigma1t += onesig_1t;
        fBkgd1sigma3p5t += onesig_3p5t;
        fBkgd1sigma2p5t += onesig_2p5t;
        fBkgd1sigma1p5t += onesig_1p5t;
        fBkgd1sigma0p5t += onesig_0p5t;
        fBkgd2sigmaFV += twosig_fv;
        fBkgd2sigma3t += twosig_3t;
        fBkgd2sigma2t += twosig_2t;
        fBkgd2sigma1t += twosig_1t;
        fBkgd2sigma3p5t += twosig_3p5t;
        fBkgd2sigma2p5t += twosig_2p5t;
        fBkgd2sigma1p5t += twosig_1p5t;
        fBkgd2sigma0p5t += twosig_0p5t;
        fBkgd3sigmaFV += threesig_fv;
        fBkgd3sigma3t += threesig_3t;
        fBkgd3sigma2t += threesig_2t;
        fBkgd3sigma1t += threesig_1t;
        fBkgd3sigma3p5t += threesig_3p5t;
        fBkgd3sigma2p5t += threesig_2p5t;
        fBkgd3sigma1p5t += threesig_1p5t;
        fBkgd3sigma0p5t += threesig_0p5t;

        //if(fVerboseLevel)
        //  std::cout << Form("Group %s SS : Total = %g , FWHM-FV = %g , FWHM-3t = %g, FWHM-1t = %g",Form("%s_%s",groupName.Data(),"SS"),tot_ss,fwhm_fv,fwhm_3t,fwhm_1t) << std::endl;
        //sum += fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];
        //roi_r3t += fGroupMeanRatio3t[Form("%s_%s",groupName.Data(),"SS")]*fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")];

        //std::cout << "Group mean ratio SS " << fGroupMeanRatio3t[Form("%s_%s",groupName.Data(),"SS")] << " mean counts " << fGroupMeanCounts[Form("%s_%s",groupName.Data(),"SS")] << std::endl;
        fGroupHistos.insert(std::make_pair(h_ss->GetName(), h_ss));
        fGroupHistos.insert(std::make_pair(h_ms->GetName(), h_ms));
    }
    if (fVerboseLevel > 0)
        std::cout << Form("Total integral : SS = %g , MS = %g , ROI = %g , 3t = %g , a3t = %g, a2t = %g, aFV = %g", tot_ss, tot_ms, roi, roi_3t, fBkgdFwhm3t, fBkgdFwhm2t, fBkgdFwhmFV) << std::endl;
}

void nEXOSensitivity::BuildWorkspace(Double_t yrs, Double_t signalCounts) {
    // Build the workspace for the fit
    if (fWsp){
        std::list<RooAbsData*> alldata = fWsp->allEmbeddedData();
        for(auto d:alldata){
            d->Delete();
        }
        
        
//        RooArgSet thepdfs = fWsp->allPdfs();
//        TIterator* pdfit = thepdfs.createIterator();
//        RooAbsPdf* thepdf = 0;
//        while((thepdf = dynamic_cast<RooAbsPdf*>(pdfit->Next()))){
//            if(thepdf->IsA()->InheritsFrom(RooHistPdf::Class())){
//                RooHistPdf* thehistpdf = dynamic_cast<RooHistPdf*>(thepdf);
//                thehistpdf->dataHist().Delete();
//            }
//        }
//        fWsp->Print();
        fWsp->Delete();
    }
    

    fWsp = new RooWorkspace("wsp");

    std::vector<TString> pdfNames;
    std::vector<TString> fitPdfNames;

    std::vector<Double_t> meanPerYear_ss;
    std::vector<Double_t> meanPerYear_ms;
    for (std::map<TString, std::vector<TString> >::iterator group = fGroups.begin(); group != fGroups.end(); group++) {
        TString groupName = group->first;

        if (fBaTag and groupName != "LXeBb2n" and groupName != "LXeBb0n")
            continue;

        if (fTurnedOffGroups.count(groupName) > 0)
            continue;

        pdfNames.push_back(groupName.Data());
        fitPdfNames.push_back(groupName.Data());

        meanPerYear_ss.push_back(fGroupMeanCounts.at(Form("%s_%s", groupName.Data(), "SS")));
        if (fVerboseLevel > 0)
            std::cout << "Mean counts group " << Form("%s_%s", groupName.Data(), "SS") << " = " << fGroupMeanCounts.at(Form("%s_%s", groupName.Data(), "SS")) << std::endl;
        meanPerYear_ms.push_back(fGroupMeanCounts.at(Form("%s_%s", groupName.Data(), "MS")));
        if (fVerboseLevel > 0)
            std::cout << "Mean counts group " << Form("%s_%s", groupName.Data(), "MS") << " = " << fGroupMeanCounts.at(Form("%s_%s", groupName.Data(), "MS")) << std::endl;
    }
    if (fVerboseLevel > 0) {
        std::cout << "Total counts groups in SS = " << std::accumulate(meanPerYear_ss.begin(), meanPerYear_ss.end(), 0.) << std::endl;
        std::cout << "Total counts groups in MS = " << std::accumulate(meanPerYear_ms.begin(), meanPerYear_ms.end(), 0.) << std::endl;
    }

    //Add the values of the systematic uncertainties to the workspace
    RooRealVar effErr("effError", "", fSignalEffError);
    RooRealVar fracErr("fracError", "", fFracError);
    RooRealVar rateErr("rateError", "", fRnRateError);
    fWsp->import(effErr);
    fWsp->import(fracErr);
    fWsp->import(rateErr);
    //Build the generating histos, and the fit pdfs and add them to the workspace
    BuildGenHistos(&pdfNames[0], pdfNames.size(), &meanPerYear_ss[0], &meanPerYear_ms[0], yrs, signalCounts);
    BuildFitPdfs(&fitPdfNames[0], fitPdfNames.size());

    //Write the workspace to file, and print the fit setup
    if (fWriteWsp) {
        std::cout << "Writing wsp into file " << fWspFileName << " ..." << std::endl;
        fWsp->writeToFile(fWspFileName.Data());
    }

    if (fVerboseLevel > 0) {
        std::cout << "////////////////////////////////////////////////////////////" << std::endl;
        std::cout << "Years of exposure = " << yrs << std::endl;
        std::cout << "Signal name = " << fSignalName << std::endl;
        std::cout << "Signal Counts = " << ((RooRealVar*) fWsp->var(Form("mean_num_%s", fSignalName.Data())))->getVal() << std::endl;
        std::cout << "Energy bins = " << ((RooRealVar*) fWsp->var("energy"))->getBins();
        std::cout << ", Energy min  = " << ((RooRealVar*) fWsp->var("energy"))->getMin();
        std::cout << ", Energy max  = " << ((RooRealVar*) fWsp->var("energy"))->getMax() << std::endl;
        std::cout << "Standoff bins = " << ((RooRealVar*) fWsp->var("standoff"))->getBins();
        std::cout << ", Standoff min  = " << ((RooRealVar*) fWsp->var("standoff"))->getMin();
        std::cout << ", Standoff max  = " << ((RooRealVar*) fWsp->var("standoff"))->getMax() << std::endl;
        std::cout << "Using signal efficiency variable in the fit: " << (fWithEff ? "yes" : "false") << std::endl;
        std::cout << "Mean signal efficiency = " << fMeanSignalEff << std::endl;
        std::cout << "SS fraction uncertainty = " << fFracError << std::endl;
        std::cout << "Rn222 rate uncertainty = " << fRnRateError << std::endl;
        std::cout << "Signal efficiency uncertainty = " << fSignalEffError << std::endl;
        std::cout << "////////////////////////////////////////////////////////////" << std::endl;
        std::cout << "Fitting pdf names: " << std::endl;

        for (size_t i = 0; i < fitPdfNames.size(); i++) {
            std::cout << fitPdfNames[i] << std::endl;
        }


        std::cout << "////////////////////////////////////////////////////////////" << std::endl;
        std::cout << Form("%18s|%9s |%9s", "Names", "Counts", "Fraction") << std::endl;
        for (size_t i = 0; i < pdfNames.size(); i++) {
            std::cout << Form("%18s|", pdfNames[i].Data());
            std::cout << Form("% 1.2e |", ((RooRealVar*) fWsp->var(Form("mean_num_%s", pdfNames[i].Data())))->getVal());
            std::cout << Form("% 1.2e", ((RooRealVar*) fWsp->var(Form("mean_frac_%s", pdfNames[i].Data())))->getVal()) << std::endl;
        }
        std::cout << "////////////////////////////////////////////////////////////" << std::endl;
    }
}

void nEXOSensitivity::AddMagicNumber(double signal, double magicN){//wrapper for python
    fMagic_numbers[signal]=magicN;
}

void nEXOSensitivity::GenAndFitData(Int_t nRuns, Double_t yrs, Double_t signalCounts, Int_t rdmRate) {
    std::cout << "Generate and fit data...\n";

    // check on mean counts randomization
    if (rdmRate <= 0)
        rdmRate = nRuns;
    if (rdmRate < nRuns)
        fRandomizeMeanCounts = true;

    LoadExcelTree(fTreeFileName.Data());
    LoadComponentHistograms();

    TString outName = fResultFileName;
    TFile outFile(outName.Data(), "recreate");
    TString suffix = fWithStandoff ? "2D" : "1D";
    //Set up a TTree to store the results of the fits
    TTree* tree = new TTree("tree", "tree");
    nEXOFitResult* fitResult = new nEXOFitResult(fSignalName.Data());

    tree->Branch("nEXOFitResult", &fitResult);

    TStopwatch clock;
    int iRun = 0;
    while (iRun < nRuns) {
        fitResult->Reset();

        if (iRun % rdmRate == 0) {
            ReadExcelTree();
            if (not fUserMeanCounts.empty()) {
                SetAllGroupMeanCounts(1e-16);
                for (std::map<TString, Double_t>::iterator userGroupCount = fUserMeanCounts.begin(); userGroupCount != fUserMeanCounts.end(); userGroupCount++)
                    SetUserMeanCounts(userGroupCount->first, userGroupCount->second);
            }
            MakeGroupHistograms();
            BuildWorkspace(yrs, signalCounts);
        }

        fitResult->bkg_tot = fBkgdTotal*fScaleBkgds;
        fitResult->bkg_fwhm_fv = fBkgdFwhmFV*fScaleBkgds;
        fitResult->bkg_fwhm_3t = fBkgdFwhm3t*fScaleBkgds;
        fitResult->bkg_fwhm_2t = fBkgdFwhm2t*fScaleBkgds;
        fitResult->bkg_fwhm_1t = fBkgdFwhm1t*fScaleBkgds;
        fitResult->bkg_fwhm_3p5t = fBkgdFwhm3p5t*fScaleBkgds;
        fitResult->bkg_fwhm_2p5t = fBkgdFwhm2p5t*fScaleBkgds;
        fitResult->bkg_fwhm_1p5t = fBkgdFwhm1p5t*fScaleBkgds;
        fitResult->bkg_fwhm_0p5t = fBkgdFwhm0p5t*fScaleBkgds;
        fitResult->bkg_1sigma_fv = fBkgd1sigmaFV*fScaleBkgds;
        fitResult->bkg_1sigma_3t = fBkgd1sigma3t*fScaleBkgds;
        fitResult->bkg_1sigma_2t = fBkgd1sigma2t*fScaleBkgds;
        fitResult->bkg_1sigma_1t = fBkgd1sigma1t*fScaleBkgds;
        fitResult->bkg_1sigma_3p5t = fBkgd1sigma3p5t*fScaleBkgds;
        fitResult->bkg_1sigma_2p5t = fBkgd1sigma2p5t*fScaleBkgds;
        fitResult->bkg_1sigma_1p5t = fBkgd1sigma1p5t*fScaleBkgds;
        fitResult->bkg_1sigma_0p5t = fBkgd1sigma0p5t*fScaleBkgds;
        fitResult->bkg_2sigma_fv = fBkgd2sigmaFV*fScaleBkgds;
        fitResult->bkg_2sigma_3t = fBkgd2sigma3t*fScaleBkgds;
        fitResult->bkg_2sigma_2t = fBkgd2sigma2t*fScaleBkgds;
        fitResult->bkg_2sigma_1t = fBkgd2sigma1t*fScaleBkgds;
        fitResult->bkg_2sigma_3p5t = fBkgd2sigma3p5t*fScaleBkgds;
        fitResult->bkg_2sigma_2p5t = fBkgd2sigma2p5t*fScaleBkgds;
        fitResult->bkg_2sigma_1p5t = fBkgd2sigma1p5t*fScaleBkgds;
        fitResult->bkg_2sigma_0p5t = fBkgd2sigma0p5t*fScaleBkgds;
        fitResult->bkg_3sigma_fv = fBkgd3sigmaFV*fScaleBkgds;
        fitResult->bkg_3sigma_3t = fBkgd3sigma3t*fScaleBkgds;
        fitResult->bkg_3sigma_2t = fBkgd3sigma2t*fScaleBkgds;
        fitResult->bkg_3sigma_1t = fBkgd3sigma1t*fScaleBkgds;
        fitResult->bkg_3sigma_3p5t = fBkgd3sigma3p5t*fScaleBkgds;
        fitResult->bkg_3sigma_2p5t = fBkgd3sigma2p5t*fScaleBkgds;
        fitResult->bkg_3sigma_1p5t = fBkgd3sigma1p5t*fScaleBkgds;
        fitResult->bkg_3sigma_0p5t = fBkgd3sigma0p5t*fScaleBkgds;

        //Get the constants from the workspace
        Double_t effError = ((RooRealVar*) fWsp->var("effError"))->getVal();
        Double_t fracError = ((RooRealVar*) fWsp->var("fracError"))->getVal();
        Double_t rateError = ((RooRealVar*) fWsp->var("rateError"))->getVal();

        //Get the observables from the workspace
        RooRealVar* energy = (RooRealVar*) fWsp->var("energy");

        RooRealVar* standoff = (RooRealVar*) fWsp->var("standoff");
        RooArgSet obs(*energy);

        //Double_t val_energy = ((RooRealVar*)fWsp->var("energy"))->getVal();

        if (fWithStandoff) obs.add(*standoff);

        //Get whether or not the signal efficiency variable is included in the fit
        Bool_t withEff = false;
        if (fWsp->var(Form("mean_eff_%s", fSignalName.Data()))) withEff = true;

        //Get the fitting pdfs
        RooAddPdf* fitPdf_ss = 0;
        RooAddPdf* fitPdf_ms = (RooAddPdf*) fWsp->pdf("pdf_ms");
        if (fWithStandoff) {
            fitPdf_ss = (RooAddPdf*) fWsp->pdf("pdf_ss");
            fitPdf_ms = (RooAddPdf*) fWsp->pdf("pdf_ms");
        } else {
            fitPdf_ss = (RooAddPdf*) fWsp->pdf("pdf1D_ss");
            fitPdf_ms = (RooAddPdf*) fWsp->pdf("pdf1D_ms");
        }

        //Fill the list of fit pdf names
        GetFitPdfNames(fitPdf_ss);

        //Get the generating pdfs
        RooAbsPdf* genPdf_ss = 0;
        RooAbsPdf* genPdf_ms = 0;
        if (fWithStandoff) {
            genPdf_ss = fWsp->pdf("genPdf_ss");
            genPdf_ms = fWsp->pdf("genPdf_ms");
        } else {
            genPdf_ss = fWsp->pdf("genPdf1D_ss");
            genPdf_ms = fWsp->pdf("genPdf1D_ms");
        }

        //TH1* genHist_ss = genPdf_ss->createHistogram(energy->GetName());
        //roi_bkg = genHist_ss->Integral(genHist_ss->FindBin(2428),genHist_ss->FindBin(2488),"width");
        //double tot_ss = genHist_ss->Integral("width");
        //std::cout << Form("Integral of gen pdf = %g and ROI = %g (%.0f - %.0f) ", tot_ss, roi_bkg, genHist_ss->GetBinLowEdge(genHist_ss->FindBin(2428)), genHist_ss->GetBinLowEdge(genHist_ss->FindBin(2488)) + genHist_ss->GetBinWidth(genHist_ss->FindBin(2488))) << std::endl;
        //delete genHist_ss;


        double expected_ss = genPdf_ss->expectedEvents(obs);
        TH2* genHist_ss = (TH2*) genPdf_ss->createHistogram(Form("%s_hist", genPdf_ss->GetName()), *energy, RooFit::YVar(*standoff)); //,RooFit::IntrinsicBinning());//,RooFit::Extended(false),RooFit::Scaling(true));
        genPdf_ss->fillHistogram(genHist_ss, RooArgList(*energy, *standoff), 1, 0, true);
        //TFile rout("test7.root","recreate");
        //genHist_ss->Write();
        //genHist_ss->ProjectionX()->Write();
        //genHist_ss->ProjectionY()->Write();
        //rout.Close();
        //TH2* genHist_ss = (TH2*)genPdf_ss->createHistogram("energy,standoff",energy->getMax()-energy->getMin(),standoff->getMax()-standoff->getMin());
        //genPdf_ss->fillHistogram(genHist_ss,RooArgList(*energy,*standoff),expected_ss,0,true);

        //expected_ss = genHist_ss->Integral();

        TAxis* xAxis = genHist_ss->GetXaxis();
        TAxis* yAxis = genHist_ss->GetYaxis();
        double pdf_roi_bkg = genHist_ss->Integral(xAxis->FindBin(2428), xAxis->FindBin(2487), 1, genHist_ss->GetNbinsY()); ///genHist_ss->Integral();
        double pdf_roi_bkg_3t = genHist_ss->Integral(xAxis->FindBin(2428), xAxis->FindBin(2487), yAxis->FindBin(90), genHist_ss->GetNbinsY()); ///genHist_ss->Integral();
        double pdf_roi_bkg_1t = genHist_ss->Integral(xAxis->FindBin(2428), xAxis->FindBin(2487), yAxis->FindBin(256), genHist_ss->GetNbinsY()); ///genHist_ss->Integral();

        //energy->setRange("fwhm",2428,2488);
        //standoff->setRange("fv",standoff->getMin(),standoff->getMax());
        //standoff->setRange("3t",90,standoff->getMax());    

        fitResult->all_bkg = expected_ss;
        fitResult->roi_bkg = pdf_roi_bkg * expected_ss;
        fitResult->roi_bkg_3t = pdf_roi_bkg_3t * expected_ss;
        fitResult->roi_bkg_1t = pdf_roi_bkg_1t * expected_ss;
        if (fVerboseLevel > 0)
            std::cout << Form("Expected events in gen pdf SS = %g and ROI = %g (%.0f - %.0f), 3t = %g and 1t = %g", expected_ss, fitResult->roi_bkg, xAxis->GetBinLowEdge(xAxis->FindBin(2428)), xAxis->GetBinLowEdge(xAxis->FindBin(2488)) + xAxis->GetBinWidth(xAxis->FindBin(2488)), fitResult->roi_bkg_3t, fitResult->roi_bkg_1t) << std::endl;
        genHist_ss->Delete();


        //Start of main fitting loop//////////////////////////////////////////////////////////////////////////
        clock.Start();
        if (iRun % 10 == 0) std::cout << "Fit number: " << iRun << std::endl;

        //Generate the ss and ms datasets
        RooDataSet* d_ss = (RooDataSet*) GenerateData(genPdf_ss, obs, false);
        RooDataSet* d_ms = (RooDataSet*) GenerateData(genPdf_ms, obs, false);
        RooDataHist* data_ss = d_ss->binnedClone();
        RooDataHist* data_ms = d_ms->binnedClone();
        
        delete d_ss; delete d_ms;
        //if (iRun%10 == 0) cout << "Fit number: " << iRun << endl;

        //Make the category to define SS and MS
        RooCategory sites("sites", "sites");
        sites.defineType("SS");
        sites.defineType("MS");
        RooDataHist combData("combData", "combined data", obs, RooFit::Index(sites), RooFit::Import("SS", *data_ss), RooFit::Import("MS", *data_ms));
        //if (iRun%10 == 0) cout << "Fit number: " << iRun << endl;

        //Make the simultaneous pdf
        RooSimultaneous simPdf("simPdf", "simPdf", sites);
        simPdf.addPdf(*fitPdf_ss, "SS");
        simPdf.addPdf(*fitPdf_ms, "MS");
        //if (iRun%10 == 0) cout << "Fit number: " << iRun << endl;

        //Make the constraints
        RooArgSet constraints;

        RooMultiVarGaussian* fracConstPdf = GetFracConstraint(fracError, true);
        constraints.add(*fracConstPdf);
        //if(not fBaTag)
        //{
        RooGaussian* rn222Const = GetRn222Constraint(rateError, true);
        if (rn222Const)
            constraints.add(*rn222Const);
        //}    

        RooGaussian* effConst = 0;
        if (withEff) {
            effConst = GetEfficiencyConstraint(effError, true);
            constraints.add(*effConst);
        }

        //Set the floating pars to random starting values or just the mean values
        bool co60Flag = false;
        for (int i = 0; i < fNFitPdfs; i++) {
            TString name = fFitPdfNames.at(i);
            double meanNum = fWsp->var(Form("mean_num_%s", name.Data()))->getVal();
            meanNum = fRandom.Gaus(meanNum, meanNum * 0.01);
            fWsp->var(Form("num_%s", name.Data()))->setVal(meanNum);
            //Pay attention to the internal Co60 ss fraction
            if (co60Flag) {//name.Contains("Co60") && !co60Flag) {
                double meanFrac = fWsp->var("mean_frac_Internal_Co60")->getVal();
                meanFrac = fRandom.Gaus(meanFrac, meanFrac * fracError);
                fWsp->var("mean_frac_Internal_Co60")->setVal(meanFrac);
                co60Flag = true;
                //} else if (name.Contains("Co60")) {
                //continue;
            } else {
                double meanFrac = fWsp->var(Form("mean_frac_%s", name.Data()))->getVal();
                meanFrac = fRandom.Gaus(meanFrac, meanFrac * 0.01);
                fWsp->var(Form("frac_%s", name.Data()))->setVal(meanFrac);
            }
        }

        if (withEff) {
            double meanEff = fWsp->var(Form("mean_eff_%s", fSignalName.Data()))->getVal();
            meanEff = fRandom.Gaus(meanEff, meanEff * 0.01);
            fWsp->var(Form("eff_%s", fSignalName.Data()))->setVal(meanEff);
        }

        if (signalCounts == 0.) fWsp->var(Form("num_%s", fSignalName.Data()))->setVal(50.);
        //Create the negative log-likelihood function to be minimized - include all external constraints
        RooAbsReal* nll = simPdf.createNLL(combData, RooFit::Extended(true), RooFit::CloneData(false), RooFit::Verbose(false), RooFit::ExternalConstraints(constraints), RooFit::NumCPU(fNcpu, 0));
        fitResult->nll_offset = nll->getVal();
        RooConstVar offset("offset", "", -fitResult->nll_offset);
        RooAddition offsetNll("offsetNll", "", RooArgSet(offset, *nll));
        //Minimize
        RooMinuit m(offsetNll);
        m.setPrintLevel(fPrintLevel);
        m.optimizeConst(true);
        m.setErrorLevel(fErrorLevel);
        m.migrad();
        //Get the best fit results for the bkgd + signal fit
        fitResult->num_signal = fWsp->var(Form("num_%s", fSignalName.Data()))->getVal();
//        fitResult->num_signal_eHi = fWsp->var(Form("num_%s", fSignalName.Data()))->getErrorHi();
//        fitResult->num_signal_eLo = fWsp->var(Form("num_%s", fSignalName.Data()))->getErrorLo();
        fitResult->fitres_sig = m.save();
        if (fVerboseLevel > 0) {
            std::cout << "Fit signal results: \n";
            fitResult->fitres_sig->Print();
        }
        fitResult->nll_sig = fitResult->fitres_sig->minNll();
        fitResult->stat_sig = fitResult->fitres_sig->status();
        fitResult->covQual_sig = fitResult->fitres_sig->covQual();
        if (fRunTruthValFit) {
            //Set the floating pars to random starting values or just the mean values
            co60Flag = false;
            for (int i = 0; i < fNFitPdfs; i++) {
                TString name = fFitPdfNames.at(i);
                if (name == fSignalName) continue;
                double meanNum = fWsp->var(Form("mean_num_%s", name.Data()))->getVal();
                meanNum = fRandom.Gaus(meanNum, meanNum * 0.001);
                fWsp->var(Form("num_%s", name.Data()))->setVal(meanNum);
                //Pay attention to the internal Co60 ss fraction
                if (co60Flag) {//name.Contains("Co60") && !co60Flag) {
                    double meanFrac = fWsp->var("mean_frac_Internal_Co60")->getVal();
                    meanFrac = fRandom.Gaus(meanFrac, meanFrac * fracError);
                    fWsp->var("mean_frac_Internal_Co60")->setVal(meanFrac);
                    co60Flag = true;
                    //} else if (name.Contains("Co60")) {
                    //continue;
                } else {
                    double meanFrac = fWsp->var(Form("mean_frac_%s", name.Data()))->getVal();
                    meanFrac = fRandom.Gaus(meanFrac, meanFrac * 0.01);
                    fWsp->var(Form("frac_%s", name.Data()))->setVal(meanFrac);
                }
            }
            //Set the signal variables constant
            fWsp->var(Form("num_%s", fSignalName.Data()))->setVal(signalCounts);
            fWsp->var(Form("num_%s", fSignalName.Data()))->setConstant(true);

            double meanFrac = fWsp->var(Form("mean_frac_%s", fSignalName.Data()))->getVal();
            fWsp->var(Form("frac_%s", fSignalName.Data()))->setVal(meanFrac);
            fWsp->var(Form("frac_%s", fSignalName.Data()))->setConstant(true);

            if (withEff) {
                double meanEff = fWsp->var(Form("mean_eff_%s", fSignalName.Data()))->getVal();
                fWsp->var(Form("eff_%s", fSignalName.Data()))->setVal(meanEff);
                fWsp->var(Form("eff_%s", fSignalName.Data()))->setConstant(true);
            }

            //Do the bkgd only fit
            
            m.setPrintLevel(fPrintLevel);
            m.optimizeConst(true);
            m.migrad();
            fitResult->fitres_bkg = m.save();
            if (fVerboseLevel > 0) {
                std::cout << "Fit background results: \n";
                fitResult->fitres_bkg->Print();
            }
            //Get the best fit results for the truth-value fit
            fitResult->fitres_bkg = m.save();
            fitResult->nll_bkg = fitResult->fitres_bkg->minNll();
            fitResult->nll_ratio = 2. * (fitResult->nll_bkg - fitResult->nll_sig);
            fitResult->stat_bkg = fitResult->fitres_bkg->status();
            fitResult->covQual_bkg = fitResult->fitres_bkg->covQual();

            //Remove constant status of floating pars
            fWsp->var(Form("num_%s", fSignalName.Data()))->setConstant(false);
            fWsp->var(Form("frac_%s", fSignalName.Data()))->setConstant(false);
            if (withEff) fWsp->var(Form("eff_%s", fSignalName.Data()))->setConstant(false);
        }
        else{
//            std::cout<<"a"<<std::endl;
            for (std::map<double,double>::reverse_iterator rit=fMagic_numbers.rbegin(); rit!=fMagic_numbers.rend(); ++rit){
                if(rit->first<fitResult->num_signal ){
                    rit--;
                    fitResult->num_signal_eHi = rit->first - fitResult->num_signal;
                    break;
                }
                if(std::distance(rit,fMagic_numbers.rend())==1){
                    fitResult->num_signal_eHi = rit->first - fitResult->num_signal;
                    break;
                }
          
                
                //Set the floating pars to random starting values or just the mean values
                co60Flag = false;
                for (int i = 0; i < fNFitPdfs; i++) {
                    TString name = fFitPdfNames.at(i);
                    if (name == fSignalName) continue;
                    double meanNum = fWsp->var(Form("mean_num_%s", name.Data()))->getVal();
                    meanNum = fRandom.Gaus(meanNum, meanNum * 0.001);
                    fWsp->var(Form("num_%s", name.Data()))->setVal(meanNum);
                    //Pay attention to the internal Co60 ss fraction
                    if (co60Flag) {//name.Contains("Co60") && !co60Flag) {
                        double meanFrac = fWsp->var("mean_frac_Internal_Co60")->getVal();
                        meanFrac = fRandom.Gaus(meanFrac, meanFrac * fracError);
                        fWsp->var("mean_frac_Internal_Co60")->setVal(meanFrac);
                        co60Flag = true;
                        //} else if (name.Contains("Co60")) {
                        //continue;
                    } else {
                        double meanFrac = fWsp->var(Form("mean_frac_%s", name.Data()))->getVal();
                        meanFrac = fRandom.Gaus(meanFrac, meanFrac * 0.01);
                        fWsp->var(Form("frac_%s", name.Data()))->setVal(meanFrac);
                    }
                }

                //Set the signal variables constant
                fWsp->var(Form("num_%s", fSignalName.Data()))->setVal(rit->first);
                fWsp->var(Form("num_%s", fSignalName.Data()))->setConstant(true);

                double meanFrac = fWsp->var(Form("mean_frac_%s", fSignalName.Data()))->getVal();
                fWsp->var(Form("frac_%s", fSignalName.Data()))->setVal(meanFrac);
                fWsp->var(Form("frac_%s", fSignalName.Data()))->setConstant(true);

                if (withEff) {
                    double meanEff = fWsp->var(Form("mean_eff_%s", fSignalName.Data()))->getVal();
                    fWsp->var(Form("eff_%s", fSignalName.Data()))->setVal(meanEff);
                    fWsp->var(Form("eff_%s", fSignalName.Data()))->setConstant(true);
                }

                //Do the bkgd only fit

                m.setPrintLevel(fPrintLevel);
                m.optimizeConst(true);
                m.migrad();
                RooFitResult* fitres_hyp = m.save();
                if (fVerboseLevel > 0) {
                    std::cout << "Fit "<<rit->first<<" hypothesis results: \n";
                    fitres_hyp->Print();
                }

                //Get the best fit results for the truth-value fit
                
                double nll_ratio = 2. * (fitres_hyp->minNll() - fitResult->nll_sig);
                                
                //Remove constant status of floating pars
                fWsp->var(Form("num_%s", fSignalName.Data()))->setConstant(false);
                fWsp->var(Form("frac_%s", fSignalName.Data()))->setConstant(false);
                if (withEff) fWsp->var(Form("eff_%s", fSignalName.Data()))->setConstant(false);
                
                if(nll_ratio<rit->second ){//experiment is consistent with this hypothesis
                    rit--;
                    fitResult->num_signal_eHi = rit->first-fitResult->num_signal;
                    break;
                }
                
            }
//            if(fitResult->num_signal_eHi+fitResult->num_signal <0.15) {
//                std::cout<<fitResult->num_signal<<"  "<<fitResult->num_signal_eHi<<std::endl;
//                throw;
//            }
            
            for (auto it=fMagic_numbers.begin(); it!=fMagic_numbers.end(); ++it) {
                if (it->first > fitResult->num_signal) {
                    if(it==fMagic_numbers.begin()){fitResult->num_signal_eLo =  fitResult->num_signal;}
                    else{
                        it--;
                        fitResult->num_signal_eLo = -(it->first - fitResult->num_signal);
                    }
                    break;
                }
                if (std::distance(it, fMagic_numbers.end()) == 1) {
                    fitResult->num_signal_eLo = -(it->first - fitResult->num_signal);
                    break;
                }
                //Set the floating pars to random starting values or just the mean values
                co60Flag = false;
                for (int i = 0; i < fNFitPdfs; i++) {
                    TString name = fFitPdfNames.at(i);
                    if (name == fSignalName) continue;
                    double meanNum = fWsp->var(Form("mean_num_%s", name.Data()))->getVal();
                    meanNum = fRandom.Gaus(meanNum, meanNum * 0.001);
                    fWsp->var(Form("num_%s", name.Data()))->setVal(meanNum);
                    //Pay attention to the internal Co60 ss fraction
                    if (co60Flag) {//name.Contains("Co60") && !co60Flag) {
                        double meanFrac = fWsp->var("mean_frac_Internal_Co60")->getVal();
                        meanFrac = fRandom.Gaus(meanFrac, meanFrac * fracError);
                        fWsp->var("mean_frac_Internal_Co60")->setVal(meanFrac);
                        co60Flag = true;
                        //} else if (name.Contains("Co60")) {
                        //continue;
                    } else {
                        double meanFrac = fWsp->var(Form("mean_frac_%s", name.Data()))->getVal();
                        meanFrac = fRandom.Gaus(meanFrac, meanFrac * 0.01);
                        fWsp->var(Form("frac_%s", name.Data()))->setVal(meanFrac);
                    }
                }

                //Set the signal variables constant
                fWsp->var(Form("num_%s", fSignalName.Data()))->setVal(it->first);
                fWsp->var(Form("num_%s", fSignalName.Data()))->setConstant(true);

                double meanFrac = fWsp->var(Form("mean_frac_%s", fSignalName.Data()))->getVal();
                fWsp->var(Form("frac_%s", fSignalName.Data()))->setVal(meanFrac);
                fWsp->var(Form("frac_%s", fSignalName.Data()))->setConstant(true);

                if (withEff) {
                    double meanEff = fWsp->var(Form("mean_eff_%s", fSignalName.Data()))->getVal();
                    fWsp->var(Form("eff_%s", fSignalName.Data()))->setVal(meanEff);
                    fWsp->var(Form("eff_%s", fSignalName.Data()))->setConstant(true);
                }

                //Do the bkgd only fit

                m.setPrintLevel(fPrintLevel);
                m.optimizeConst(true);
                m.migrad();
                RooFitResult* fitres_hyp = m.save();
                if (fVerboseLevel > 0) {
                    std::cout << "Fit "<<it->first<<" hypothesis results: \n";
                    fitres_hyp->Print();
                }

                //Get the best fit results for the truth-value fit
                
                double nll_ratio = 2. * (fitres_hyp->minNll() - fitResult->nll_sig);
                                
                //Remove constant status of floating pars
                fWsp->var(Form("num_%s", fSignalName.Data()))->setConstant(false);
                fWsp->var(Form("frac_%s", fSignalName.Data()))->setConstant(false);
                if (withEff) fWsp->var(Form("eff_%s", fSignalName.Data()))->setConstant(false);
                
                if(nll_ratio<it->second){//experiment is consistent with this hypothesis
                    it--;
                    fitResult->num_signal_eLo = -it->first+fitResult->num_signal;
                    break;
                }
                
            }
            std::cout<<"ll: "<<fitResult->num_signal-fitResult->num_signal_eLo<<std::endl;
            
        }

        fitResult->real_time = clock.RealTime();

        //Fill the tree
        tree->Fill();

        //RooPlot* frame = energy->frame();
        //data_ss->plotOn(frame);
        //fitPdf_ss->plotOn(frame);
        //TCanvas* cc = new TCanvas("cc","cc");
        //frame->Draw();
        //cc->SaveAs("cc.C");

        //Cleanup
        delete nll;
        delete data_ss;
        delete data_ms;

        delete fracConstPdf;
        if (withEff) delete effConst;
        //delete rn222Const;

        iRun++;
    }
    for (auto pr : fComponentHistos) {
        delete pr.second;
    }
    fComponentHistos.clear();
    hdata1D_ss.clear();
    hdata1D_ms.clear();
    hdata_ss.clear();
    hdata_ms.clear();
    //Write the tree to file
    outFile.cd();
//    outFile.Print();
    tree->Write();
    outFile.Close();
    

}


TH2D* nEXOSensitivity::MakeCombinedHisto(TString histName, Int_t nComp, TString* compNames, Double_t* fractions, TString , Bool_t isSS) {
    if (fVerboseLevel > 0)
        std::cout << "Working on " << histName << std::endl;
    TString sites = (isSS ? "SS" : "MS");
    TH2D* retHist = new TH2D(histName, "", fNbinsX, fXmin, fXmax, fNbinsY, fYmin, fYmax); //fYbins);//fYmin, fYmax); 

    TH2* compHist = 0;

    Double_t overallNorm = 0.;

    //double combRatioFV(0.), combRatio3t(0.), combRatio1t(0.);
    //double combFull(0.);
    for (int i = 0; i < nComp; i++) {
        TString compName = Form("%s_%s", compNames[i].Data(), (isSS ? "ss" : "ms"));
        compHist = dynamic_cast<TH2*> (fComponentHistos.at(compName.Data()));
        //int nBinsXcomp = compHist->GetNbinsX();
        //int xRebin = (int) nBinsXcomp / fNbinsX;
        //int nBinsYcomp = compHist->GetNbinsY();
        //int yRebin = (int) nBinsYcomp / fNbinsY;    
        //compHist->Rebin2D(xRebin, yRebin);

        overallNorm += fractions[i]; //compHist->Integral()*fractions[i];

        double histIntegral = compHist->Integral();
        compHist->Scale(fractions[i] / histIntegral);

        if (fVerboseLevel > 0)
            std::cout << Form("Group fractions histName : %s , comp : %s , fraction : %g , integral : %g", histName.Data(), compNames[i].Data(), fractions[i], histIntegral) << std::endl;

        retHist->Add(compHist);

        //combFull += fractions[i]*fNormHistoFull.at(compName);
        //combRatioFV += fractions[i]*fRatioHistoFV.at(compName);
        //combRatio3t += fractions[i]*fRatioHisto3t.at(compName);
        //combRatio1t += fractions[i]*fRatioHisto1t.at(compName);
    }

    //if(isSS)
    //{
    //  fRatioHistoFV.insert(std::make_pair(histName,combRatioFV/combFull));
    //  fRatioHisto3t.insert(std::make_pair(histName,combRatio3t/combFull));
    //  fRatioHisto1t.insert(std::make_pair(histName,combRatio1t/combFull));
    //}

    retHist->Scale(1. / overallNorm);
    //if(fVerboseLevel > 0)
    //  std::cout << Form("Group integral : %g , ROI = %g", retHist->Integral(), retHist->Integral(retHist->GetXaxis()->FindBin(2428),retHist->GetXaxis()->FindBin(2488),1,fNbinsY)) << std::endl;

    return retHist;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//This function does two necessary things and must be called first.
//1. It loads all component pdfs into the workspace.
//2. It builds the generating pdfs

void nEXOSensitivity::BuildGenHistos(TString* pdfNames, const int nPdfs, Double_t* meanPerYear_ss, Double_t* meanPerYear_ms, Double_t yrs, Double_t signalCounts) {
    //Get the file containing the simulated histograms

    //TFile fIn(fHistoFileName.Data());

    //Get the bins, and ranges for the histograms (they all have to have same binning)
    Int_t nBinsEnergy = fNbinsX; //((TH2D*)fIn.Get("h_bb0n_ss"))->GetNbinsX();
    Double_t energyLo = fXmin; //((TH2D*)fIn.Get("h_bb0n_ss"))->GetXaxis()->GetXmin();
    Double_t energyHi = fXmax; //((TH2D*)fIn.Get("h_bb0n_ss"))->GetXaxis()->GetXmax();	

    Int_t nBinsStandoff = fNbinsY; //((TH2D*)fIn.Get("h_bb0n_ss"))->GetNbinsY();
    Double_t standoffLo = fYmin; //((TH2D*)fIn.Get("h_bb0n_ss"))->GetYaxis()->GetXmin();
    Double_t standoffHi = fYmax; //((TH2D*)fIn.Get("h_bb0n_ss"))->GetYaxis()->GetXmax();	

    //Create the fitting variables
    RooRealVar energy("energy", "energy", (energyLo + energyHi) / 2., energyLo, energyHi);
    energy.setBins(nBinsEnergy);

    RooRealVar standoff("standoff", "standoff", (standoffLo + standoffHi) / 2., standoffLo, standoffHi);
    standoff.setBins(nBinsStandoff);
    //RooBinning yBins(fNbinsY,fYbins);
    //standoff.setBinning(yBins);

    RooArgSet obs;
    obs.add(energy);
    obs.add(standoff);

    //To make the individual 2D and 1D hist pdfs
    TH1 *histo_ss = 0;
    TH1 *histo_ms = 0;
    hdata_ss.clear();
    hdata_ms.clear();
    RooHistPdf *pdf_ss = 0;
    RooHistPdf *pdf_ms = 0;

    hdata1D_ss.clear();
    hdata1D_ms.clear();
    RooHistPdf *pdf1D_ss = 0;
    RooHistPdf *pdf1D_ms = 0;

    //Making the sum generating pdf
    RooRealVar* mean_num_ss = 0;
    RooRealVar* mean_num_ms = 0;
    RooArgList pdfList_ss;
    RooArgList pdfList_ms;
    RooArgList pdfList1D_ss;
    RooArgList pdfList1D_ms;
    RooArgList coefList_ss;
    RooArgList coefList_ms;

    //Lists and variables to store the mean values
    RooRealVar* mean_num = 0;
    RooRealVar* mean_frac = 0;

    //Loop through the pdf components and save the pdfs, mean fractions, mean counts, etc.
    for (int i = 0; i < nPdfs; i++) {
        TString name = pdfNames[i];

        //Calculated the mean expected events for each of the components - if the component is Co60, scale appropriately
        Double_t mean_ss = meanPerYear_ss[i] * yrs;
        Double_t mean_ms = meanPerYear_ms[i] * yrs;

        if (not name.Contains("LXe")) // check if beta-like PDF
        {
            //std::cout << "improving fraction " << name << " " << frac << " " << fSSFracImprovement << std::endl;
            //frac *= fSSFracImprovement;
            mean_ss *= fSSFracImprovement;
            mean_ms *= 1 + (1 - fSSFracImprovement) * mean_ss / mean_ms;
        }

        if (name.Contains("Rn222")) {
            mean_ss *= fRn222RateCorrection;
            mean_ms *= fRn222RateCorrection;
        }


        Double_t mean = mean_ss + mean_ms;
        Double_t frac = mean_ss / mean;


        if (name.Contains("Co60")) {
            // N0 = Ndecayed(1 year) / (1. - exp(-tau))
            // Ndecayed(years) = N0*(1.- exp(-tau*yrs)) = Ndecayed(1 year)*(1.- exp(-tau*yrs)) / (1. - exp(-tau))
            // tau = ln(2)/t_1/2, t_1/2 = 5.27yrs
            Double_t tau = log(2.) / 5.27;
            Double_t scale = (1. - exp(-tau * yrs)) / (1. - exp(-tau));
            mean_ss = mean_ss * scale / yrs;
            mean_ms = mean_ms * scale / yrs;
            mean = mean_ss + mean_ms;
        }
        if (name.Contains("LXeBb0n")) {
            mean = signalCounts;
            mean_ss = signalCounts*frac;
            mean_ms = signalCounts * (1. - frac);
        } else {
            if (fScale2nu or (not name.Contains("LXeBb2n")))
                mean *= fScaleBkgds;
        }
        mean_num = new RooRealVar(Form("mean_num_%s", name.Data()), "", mean);
        mean_frac = new RooRealVar(Form("mean_frac_%s", name.Data()), "", frac);
        mean_num_ss = new RooRealVar(Form("mean_num_ss_%s", name.Data()), "", mean * frac);
        mean_num_ms = new RooRealVar(Form("mean_num_ms_%s", name.Data()), "", mean * (1. - frac));

        //Get the MC histograms, and create RooHistPdfs from them
        histo_ss = fGroupHistos.at(Form("h_%s_ss", name.Data())); //(TH2D*)fIn.Get(Form("h_%s_ss", name.Data()));
        hdata_ss.emplace_back(new RooDataHist(Form("dh_%s_ss", name.Data()), "", obs, histo_ss));
        pdf_ss = new RooHistPdf(Form("pdf_%s_ss", name.Data()), "", obs, *hdata_ss.back().get(), fInterpOrder);

        
        histo_ms = fGroupHistos.at(Form("h_%s_ms", name.Data())); //(TH2D*)fIn.Get(Form("h_%s_ms", name.Data()));
        hdata_ms.emplace_back(new RooDataHist(Form("dh_%s_ms", name.Data()), "", obs, histo_ms));
        pdf_ms = new RooHistPdf(Form("pdf_%s_ms", name.Data()), "", obs, *hdata_ms.back().get(), fInterpOrder);

        TH2* histo_ss_2d = dynamic_cast<TH2*> (histo_ss);
        TH2* histo_ms_2d = dynamic_cast<TH2*> (histo_ms);
        if (histo_ss_2d && histo_ms_2d) {
            hdata1D_ss.emplace_back(new RooDataHist(Form("dh1D_%s_ss", name.Data()), "", energy, histo_ss_2d->ProjectionX()));
            pdf1D_ss = new RooHistPdf(Form("pdf1D_%s_ss", name.Data()), "", energy, *hdata1D_ss.back().get(), fInterpOrder);
  
            hdata1D_ms.emplace_back(new RooDataHist(Form("dh1D_%s_ms", name.Data()), "", energy, histo_ms_2d->ProjectionX()));
            pdf1D_ms = new RooHistPdf(Form("pdf1D_%s_ms", name.Data()), "", energy, *hdata1D_ms.back().get(), fInterpOrder);
  
        }

        coefList_ss.addOwned(*mean_num_ss);
        coefList_ms.addOwned(*mean_num_ms);
        pdfList_ss.addOwned(*pdf_ss);
        pdfList_ms.addOwned(*pdf_ms);
        pdfList1D_ss.addOwned(*pdf1D_ss);
        pdfList1D_ms.addOwned(*pdf1D_ms);
        
        //Add the pdfs and constants to the workspace
//        fWsp->import(*pdf_ss);
//        fWsp->import(*pdf_ms);
//        fWsp->import(*pdf1D_ss);
//        fWsp->import(*pdf1D_ms);
        fWsp->import(*mean_num);
        fWsp->import(*mean_frac);
         mean_num->Delete();
         mean_frac->Delete();
         
    }
    //Create the sum generating hist pdf and import into the workspace
    RooAddPdf genPdf_ss("genPdf_ss", "", pdfList_ss, coefList_ss);
    RooAddPdf genPdf_ms("genPdf_ms", "", pdfList_ms, coefList_ms);
    fWsp->import(genPdf_ss, RooFit::RecycleConflictNodes(true));
    fWsp->import(genPdf_ms, RooFit::RecycleConflictNodes(true));

    //energy.setRange("FULL_RANGE",700,3500);
    //standoff.setRange("FULL_RANGE",0,650);

    //double alt_gen_full = genPdf_ss.createIntegral(RooArgSet(standoff,energy),"FULL_RANGE")->getVal();
    //double alt_gen_roi = genPdf_ss.createIntegral(RooArgSet(standoff,energy),"ROI_RANGE")->getVal();

    ////TH2* genpdf_hist_ss_2d = (TH2*)genPdf_ss.createHistogram("energy,standoff");
    //TH2* genpdf_hist_ss_2d = (TH2*)genPdf_ss.createHistogram("genpdf_ss",energy,RooFit::YVar(standoff),RooFit::IntrinsicBinning(),RooFit::Extended(false));
    //double genprod_ss = genpdf_hist_ss_2d->Integral();
    //double genprod_roi = genpdf_hist_ss_2d->Integral(genpdf_hist_ss_2d->GetXaxis()->FindBin(2428),genpdf_hist_ss_2d->GetXaxis()->FindBin(2488),1,genpdf_hist_ss_2d->GetNbinsY());
    //double expected = genPdf_ss.expectedEvents(obs);
    //std::cout << Form("Pdf SS: alt full = %g, roi = %g, prod = %g , prod roi = %g, exp = %g",alt_gen_full,alt_gen_roi,genprod_ss,genprod_roi,expected) << std::endl;
    RooAddPdf genPdf1D_ss("genPdf1D_ss", "", pdfList1D_ss, coefList_ss);
    RooAddPdf genPdf1D_ms("genPdf1D_ms", "", pdfList1D_ms, coefList_ms);
    fWsp->import(genPdf1D_ss, RooFit::RecycleConflictNodes(true));
    fWsp->import(genPdf1D_ms, RooFit::RecycleConflictNodes(true));
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Build the fit pdf in SS and MS

void nEXOSensitivity::BuildFitPdfs(TString* pdfNames, const int nFitPdfs) {
    BuildFitPdf(pdfNames, nFitPdfs, true);
    BuildFitPdf(pdfNames, nFitPdfs, false);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Build the fit pdf in SS or MS

void nEXOSensitivity::BuildFitPdf(TString* pdfNames, const int nFitPdfs, Bool_t isSS) {
    //Suffix and formula for ss or ms
    TString suffix = "ms";
    TString formula = "@0*(1.-@1)";
    if (isSS) {
        suffix = "ss";
        formula = "@0*@1";
    }

    RooHistPdf* pdf = 0; //Individual fit pdf
    RooHistPdf* pdf1D = 0; //Energy only fit pdf
    num_tomakefit.clear(); //number of events for component (floating)
    frac_tomakefit.clear(); //ss fraction for component (floating)
    coef_tomakefit.clear(); //coefficient in sum pdf for component

    RooArgList pdfList; //List of component pdfs
    RooArgList pdfList1D; //List of 1D component pdfs
    RooArgList coefList; //List of component coeficients

    for (Int_t i = 0; i < nFitPdfs; i++) {
        TString name = pdfNames[i];

        if (fVerboseLevel > 0)
            std::cout << "Adding group " << name << " to pdfs... " << std::endl;

        //We'll handle the signal pdf separately
        if (pdfNames[i] == fSignalName) continue;

        //Get the mean number of events and ss fraction for the component
        Double_t meanNum = ((RooRealVar*) fWsp->var(Form("mean_num_%s", name.Data())))->getVal();
        Double_t meanFrac = ((RooRealVar*) fWsp->var(Form("mean_frac_%s", name.Data())))->getVal();

        //Create the component's coefficient
        num_tomakefit.emplace_back(new RooRealVar(Form("num_%s", name.Data()), "", meanNum, 0., meanNum * 10.)); //meanNum*(1.-0.3), meanNum*(1.+0.3));
        frac_tomakefit.emplace_back(new RooRealVar(Form("frac_%s", name.Data()), "", meanFrac, 0., 1.)); //meanFrac*(1.-3.*fFracError), meanFrac*(1.+3.*fFracError));
        coef_tomakefit.emplace_back(new RooFormulaVar(Form("num_%s_%s", name.Data(), suffix.Data()), "", formula.Data(), RooArgList(*num_tomakefit.back().get(), *frac_tomakefit.back().get())));

        //Get the pdf from the workspace
        pdf = (RooHistPdf*) fWsp->pdf(Form("pdf_%s_%s", name.Data(), suffix.Data()));
        pdf1D = (RooHistPdf*) fWsp->pdf(Form("pdf1D_%s_%s", name.Data(), suffix.Data()));

        //Add the component pdf and coefficient to the list
        pdfList.add(*pdf);
        pdfList1D.add(*pdf1D);
        coefList.add(*coef_tomakefit.back().get());
    }

    //Create the signal component
    Double_t meanNum = ((RooRealVar*) fWsp->var(Form("mean_num_%s", fSignalName.Data())))->getVal();
    Double_t meanFrac = ((RooRealVar*) fWsp->var(Form("mean_frac_%s", fSignalName.Data())))->getVal();

    num_tomakefit.emplace_back(new RooRealVar(Form("num_%s", fSignalName.Data()), "", meanNum, 0., 500.));
    frac_tomakefit.emplace_back(new RooRealVar(Form("frac_%s", fSignalName.Data()), "", meanFrac, 0., 1.)); //meanFrac*(1.-3.*fFracError), meanFrac*(1.+3.*fFracError));
    RooRealVar* efficiency = 0;
    if (fWithEff) {
        RooRealVar meanSignalEfficiency(Form("mean_eff_%s", fSignalName.Data()), "", fMeanSignalEff);
        fWsp->import(meanSignalEfficiency);
        efficiency = new RooRealVar(Form("eff_%s", fSignalName.Data()), "", fMeanSignalEff, fMeanSignalEff * (1. - fSignalEffError * 4.), fMeanSignalEff * (1. + fSignalEffError * 4.));
        formula += "*@2";
        //efficiency->setConstant(true);
        coef_tomakefit.emplace_back(new RooFormulaVar(Form("num_%s_%s", fSignalName.Data(), suffix.Data()), "", formula.Data(), RooArgList(*num_tomakefit.back().get(), *frac_tomakefit.back().get(), *efficiency)));
    } else {
        coef_tomakefit.emplace_back(new RooFormulaVar(Form("num_%s_%s", fSignalName.Data(), suffix.Data()), "", formula.Data(), RooArgList(*num_tomakefit.back().get(), *frac_tomakefit.back().get())));
    }
    pdf = (RooHistPdf*) fWsp->pdf(Form("pdf_%s_%s", fSignalName.Data(), suffix.Data()));
    pdf1D = (RooHistPdf*) fWsp->pdf(Form("pdf1D_%s_%s", fSignalName.Data(), suffix.Data()));
    pdfList.add(*pdf);
    pdfList1D.add(*pdf1D);
    coefList.add(*coef_tomakefit.back().get());

    //Create the sum pdf and import it into the workspace - reusing any common functions/pdfs already in the workspace
    RooAddPdf fitPdf(Form("pdf_%s", suffix.Data()), "", pdfList, coefList);
    fWsp->import(fitPdf, RooFit::RecycleConflictNodes(true));
    RooAddPdf fitPdf1D(Form("pdf1D_%s", suffix.Data()), "", pdfList1D, coefList);
    fWsp->import(fitPdf1D, RooFit::RecycleConflictNodes(true));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Make the Rn222 rate constraint

RooGaussian* nEXOSensitivity::GetRn222Constraint(Double_t rateError, Bool_t randomize) {

    RooRealVar* num = (RooRealVar*) fWsp->var("num_LXeRn222");
    if (not num)
        return 0;

    Double_t meanNum = ((RooRealVar*) fWsp->var("mean_num_LXeRn222"))->getVal();
    Double_t meanErr = meanNum*rateError;

    //If we want to draw a random constraint value
    if (randomize) {
        meanNum = fRandom.Gaus(meanNum, meanErr);
    }
    RooGaussian* rn222Const = new RooGaussian("rn222Const", "", *num, RooFit::RooConst(meanNum), RooFit::RooConst(meanErr));

    return rn222Const;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Make the bb0n efficiency constraint

RooGaussian* nEXOSensitivity::GetEfficiencyConstraint(Double_t effError, Bool_t randomize) {

    RooRealVar* eff = (RooRealVar*) fWsp->var("eff_LXeBb0n");
    Double_t meanEff = ((RooRealVar*) fWsp->var("mean_eff_LXeBb0n"))->getVal();
    Double_t meanErr = meanEff*effError;

    //If we want to draw a random constraint value
    if (randomize) {
        meanEff = fRandom.Gaus(meanEff, meanErr);
    }
    RooGaussian* effConst = new RooGaussian("effConst", "", *eff, RooFit::RooConst(meanEff), RooFit::RooConst(meanErr));

    return effConst;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Make the SS fraction constraints

RooMultiVarGaussian* nEXOSensitivity::GetFracConstraint(Double_t fracError, Bool_t randomize) {
    RooArgList varList; //list of ss frac parameters
    RooArgList meanList; //list of ss frac mean values
    std::vector<Double_t> sigmas; //list of errors on ss frac parameters
    RooRealVar* frac = 0;

    //Bool_t co60Flag = false; //to keep track of TPC Co60, as these components share their ss frac parameter
    for (int i = 0; i < fNFitPdfs; i++) {
        TString name = fFitPdfNames.at(i); //get the fit pdf name, and replace the prefix, and suffix to get base component name
        name.ReplaceAll("pdf_", "");
        name.ReplaceAll("pdf1D_", "");
        name.ReplaceAll("_ss", "");
        name.ReplaceAll("_ms", "");

        //if (name.Contains("Co60")) { //check to see if this is Co60
        //	name = "Internal_Co60";
        //}

        frac = (RooRealVar*) fWsp->var(Form("frac_%s", name.Data())); //get the ss frac par from the workspace

        Double_t meanFrac = ((RooRealVar*) fWsp->var(Form("mean_frac_%s", name.Data())))->getVal(); //get the mean value
        Double_t meanErr = fracError*meanFrac; //get the error
        frac->setVal(meanFrac);
        sigmas.push_back(meanErr);

        //If we want to draw random constraints - renormalize the gaussian over the allowed range of 0 to 1
        if (randomize) {
            Double_t randMean = 0.;
            do {
                randMean = fRandom.Gaus(meanFrac, meanErr);
            } while (randMean > 1.);
            meanFrac = randMean;
        }

        varList.add(*frac);
        meanList.add(RooFit::RooConst(meanFrac));
    }

    //Create the covariance matrix - assume uncorrelated parameters for now
    TMatrixDSym covMat(varList.getSize());
    for (int i = 0; i < varList.getSize(); i++) {
        for (int j = 0; j < varList.getSize(); j++) {
            covMat[i][j] = 0.;
            if (i == j) covMat[i][j] = sigmas.at(i) * sigmas.at(i);
        }
    }

    //Create the ss fraction constraint
    RooMultiVarGaussian* fracConst = new RooMultiVarGaussian("fracConst", "", varList, meanList, covMat);
    return fracConst;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Fill the vector list of the fitting pdf names

void nEXOSensitivity::GetFitPdfNames(RooAddPdf* fitPdf) {
    RooArgList pdfList = fitPdf->pdfList();
    fNFitPdfs = (int) pdfList.getSize();
    for (int i = 0; i < fNFitPdfs; i++) {
        TString name = pdfList.at(i)->GetName();
        name.ReplaceAll("pdf_", "");
        name.ReplaceAll("pdf1D_", "");
        name.ReplaceAll("_ss", "");
        name.ReplaceAll("_ms", "");
        fFitPdfNames.push_back(name);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generate a fake dataset

RooAbsData* nEXOSensitivity::GenerateData(RooAbsPdf* genPdf, RooArgSet obs, Bool_t isBinned) {

    //generate the dataset
    RooDataSet* data = genPdf->generate(obs, RooFit::Extended(true));
    if (isBinned) { //Return binned dataset if wanted
        RooDataHist* bdata = data->binnedClone("data");
        delete data;
        return bdata;
    } else { //Return unbinned dataset
        data->SetName("data");
        return data;
    }
}

Double_t nEXOSensitivity::EvalCounts(Double_t hitEfficiency, Double_t activity, Double_t time, Double_t halflife) {
    if (time / halflife > 0.01) {
        double lhl = TMath::Log2(halflife);
        time = lhl * (1 - exp(-time / lhl));
    }

    double counts = time * hitEfficiency * activity * 31556736; // seconds per year conversion

    return counts;
}

Double_t nEXOSensitivity::GenTruncGaus(Double_t mean, Double_t sigma, Double_t unif, Int_t low) {
    // generate truncated Gauss(mean,sigma) sample from unif[0,1] number
    // low truncation at low (default = 0)
    // code from https://people.sc.fsu.edu/~jburkardt/cpp_src/truncated_normal/truncated_normal.html

    double alpha = (low - mean) / sigma;
    double alpha_cdf = ROOT::Math::normal_cdf(alpha);
    double xi_cdf = alpha_cdf + unif * (1.0 - alpha_cdf);
    double xi = ROOT::Math::gaussian_quantile(xi_cdf, 1);

    double x = mean + sigma * xi;

    //alpha = ( low - mean ) / sigma;
    //alpha_cdf = normal_01_cdf ( alpha );
    //u = r8_uniform_01 ( seed );
    //xi_cdf = alpha_cdf + u * ( 1.0 - alpha_cdf );
    // xi = normal_01_cdf_inv ( xi_cdf );
    //x = mu + sigma * xi;

    return x;
}
