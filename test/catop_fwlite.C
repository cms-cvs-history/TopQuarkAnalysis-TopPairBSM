#include "DataFormats/FWLite/interface/Handle.h"


#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"
#endif

#include <iostream>
#include <map>
#include <string>

using namespace std;

void catop_fwlite()
{
   
   
  TFile  * file = new TFile("/uscms_data/d1/rappocc/zprime_ca_214_output.root");

  using namespace std;
  using namespace reco;

  TH1D * hist_top_jetPt = new TH1D("hist_top_jetPt", "Jet p_{T}", 50, 0, 2000 );
  TH1D * hist_nontop_jetPt = new TH1D("hist_nontop_jetPt", "Jet p_{T}", 50, 0, 2000 );

  TH1D * hist_top_jetEta = new TH1D("hist_top_jetEta", "Jet #eta", 50, -3.0, 3.0 );
  TH1D * hist_nontop_jetEta = new TH1D("hist_nontop_jetEta", "Jet #eta", 50, -3.0, 3.0);

  TH1D * hist_top_jetMass = new TH1D("hist_top_jetMass", "Jet Mass", 50, 0, 500 );
  TH1D * hist_nontop_jetMass = new TH1D("hist_nontop_jetMass", "Jet Mass", 50, 0, 500 );

  TH1D * hist_top_jetMinMass = new TH1D("hist_top_jetMinMass", "Jet Min Mass", 50, 0, 200 );
  TH1D * hist_nontop_jetMinMass = new TH1D("hist_nontop_jetMinMass", "Jet Min Mass", 50, 0, 200 );

  TH1D * hist_top_jetWMass = new TH1D("hist_top_jetWMass", "Jet p_{T}", 50, 0, 200 );
  TH1D * hist_nontop_jetWMass = new TH1D("hist_nontop_jetWMass", "Jet p_{T}", 50, 0, 200 );

  fwlite::Event ev(file);

  map<string, int> top_counts;
  map<string, int> nontop_counts;

  const char * cutnames [] = {
    "Jet Pt Cuts" ,
    ">= 3 Subjets",
    "Jet Mass Cut",
    "W Mass Cut"  ,
    "Min Mass Cut"
  };

  top_counts[cutnames[0]] = 0;
  top_counts[cutnames[1]] = 0;
  top_counts[cutnames[2]] = 0;
  top_counts[cutnames[3]] = 0;
  top_counts[cutnames[4]] = 0;
  

  nontop_counts[cutnames[0]] = 0;
  nontop_counts[cutnames[1]] = 0;
  nontop_counts[cutnames[2]] = 0;
  nontop_counts[cutnames[3]] = 0;
  nontop_counts[cutnames[4]] = 0;
  
  
  for( ev.toBegin();
          ! ev.atEnd();
          ++ev) {


    fwlite::Handle<std::vector<pat::Jet> > h_jet;

    h_jet   .getByLabel(ev,"selectedLayer1Jets");

    vector<pat::Jet> const & jets = *h_jet;

     for ( int i = 0; i < jets.size();  ++i ) {

       const reco::CATopJetTagInfo * catopTag = dynamic_cast<CATopJetTagInfo*>(jets[i].tagInfo("CATopJetTagger"));
       
       if ( abs(jets[i].partonFlavour()) == 6) {

	 hist_top_jetPt->Fill( jets[i].pt() );
	 hist_top_jetEta->Fill( jets[i].eta() );
	 hist_top_jetMass->Fill( catopTag->properties().topMass );
	 hist_top_jetMinMass->Fill( catopTag->properties().minMass );
	 hist_top_jetWMass->Fill( catopTag->properties().wMass );

	 top_counts["Jet Pt Cuts"]++;

	 if  ( jets[i].getJetConstituents().size() < 3 )  continue;
	 top_counts[">= 3 Subjets"]++;

	 if ( catopTag->properties().topMass < 100 || catopTag->properties().topMass > 250 ) continue;
	 top_counts["Jet Mass Cut"]++;

	 if ( catopTag->properties().wMass < 0 ) continue;
	 top_counts["W Mass Cut"]++;

	 if ( catopTag->properties().minMass < 50 ) continue;
	 top_counts["Min Mass Cut"]++;
       }
       else {

	 hist_nontop_jetPt->Fill( jets[i].pt() );
	 hist_nontop_jetEta->Fill( jets[i].eta() );
	 hist_nontop_jetMass->Fill( catopTag->properties().topMass );
	 hist_nontop_jetMinMass->Fill( catopTag->properties().minMass );
	 hist_nontop_jetWMass->Fill( catopTag->properties().wMass );

	 nontop_counts["Jet Pt Cuts"]++;

	 if  ( jets[i].getJetConstituents().size() < 3 )  continue;
	 nontop_counts[">= 3 Subjets"]++;

	 if ( catopTag->properties().topMass < 100 || catopTag->properties().topMass > 250 ) continue;
	 nontop_counts["Jet Mass Cut"]++;

	 if ( catopTag->properties().wMass < 0 ) continue;
	 nontop_counts["W Mass Cut"]++;

	 if ( catopTag->properties().minMass < 50 ) continue;
	 nontop_counts["Min Mass Cut"]++;
       }
     }

   }

  cout << "About to print summaries" << endl;

  for ( int i = 0; i < top_counts.size(); ++i ) {
    char buff[1000];
    sprintf(buff, "%20s  : %6d   %6d", cutnames[i], nontop_counts[cutnames[i]], top_counts[cutnames[i]]);
    cout << buff << endl;
  }


  TFile * f = new TFile("histograms_catop_fwlite.root", "RECREATE");
  f->cd();

  hist_nontop_jetPt->Write();
  hist_top_jetPt->Write();

  hist_nontop_jetEta->Write();
  hist_top_jetEta->Write();

  hist_nontop_jetMinMass->Write();
  hist_top_jetMinMass->Write();

  hist_nontop_jetMass->Write();
  hist_top_jetMass->Write();

  hist_nontop_jetWMass->Write();
  hist_top_jetWMass->Write();  
}
