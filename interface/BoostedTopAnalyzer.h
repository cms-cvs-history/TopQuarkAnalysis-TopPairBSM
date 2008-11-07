#ifndef BoostedTopAnalyzer_H
#define BoostedTopAnalyzer_H

/** \class BoostedTopAnalyzer
 *
 *  Author: Francisco Yumiceva
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h" 
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h" 

//#include "Analyzers/TQAFAnalyzer/test/TtEventDummyAnalysisHisto.h"
#include "TopQuarkAnalysis/TopPairBSM/interface/TPBHistograms.h"
#include "TopQuarkAnalysis/TopPairBSM/interface/JetCombinatorics.h"

#include "TFile.h"
#include "TLorentzVector.h"
#include <vector>
#include <map>
#include <string>
#include <utility> 
#include <fstream>

class BoostedTopAnalyzer : public edm::EDAnalyzer {

  public:


    /// Constructor
    BoostedTopAnalyzer(const edm::ParameterSet& pset);

    /// Destructor
    virtual ~BoostedTopAnalyzer();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);

    /// Jet to parton matching
    bool IsTruthMatch( Combo acombo, const std::vector<pat::Jet> jets );

    /// 3d angles
	double Psi(TLorentzVector p1, TLorentzVector p2, double mass);
	double dij(TLorentzVector p1, TLorentzVector p2, double mass, bool min = true);
	double PtRel(TLorentzVector p, TLorentzVector paxis);
	
	typedef math::XYZTLorentzVector LorentzVector;
	
  private:


    // Histogram containers
	TPBHistograms *hcounter;
	TPBHistograms *hmuons_;
	TPBHistograms *hmet_;
	TPBHistograms *hjets_;
	TPBHistograms *hgen_;
	TPBHistograms *hmass_;
	TPBHistograms *hdisp_;
	
	std::map<TString, TString> cut_map;
	    
    // The file which will store the histos
    TFile *theFile;
    // ascii outpur filename
	std::ofstream fasciiFile;

	JetCombinatorics myCombi_;
        JetCombinatorics myCombi0_;
	JetCombinatorics myCombi2_;
	JetCombinatorics myCombi3_;
	
    // configuration
	bool fwriteAscii;// flag to dump ASCII file
    std::string fasciiFileName; // ASCII filename
    // csa07 weights
    bool fApplyWeights;
    // verbose
    bool debug;
	bool fdisplayJets; // make lego plots
    int feventToProcess;
	
    std::string rootFileName;
    int leptonFlavor;
	edm::InputTag genEvnSrc;
    edm::InputTag leptonSrc;
    //edm::InputTag electronSrc;
    edm::InputTag metSrc;
    edm::InputTag jetSrc;
	edm::InputTag jetSrc1;
    edm::InputTag jetSrc2;

    edm::InputTag evtsols;

	int nevents;
	int nbadmuons;
	int nWcomplex;
	
	double fMinLeptonPt;
	double fMinLeptonEta;
	double fTrackIso;
	double fCaloIso;
	double fMinLeadingJetEt;
	double fMinJetEt;
	double fMinJetEta;
	double fMinHt;
	double fMinMET;
	
};


#endif
