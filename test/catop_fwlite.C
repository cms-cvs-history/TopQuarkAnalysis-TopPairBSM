#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include <TH1.h>
#include <TFile.h>


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
   
   
  TFile  * file = new TFile("PATLayer1_Output.fromAOD_full.root");

  fwlite::Event ev(file);
 
  for( ev.toBegin();
          ! ev.atEnd();
          ++ev) {


    fwlite::Handle<std::vector<pat::Jet> > h_jet;

    h_jet   .getByLabel(ev,"selectedLayer1Jets");

    vector<pat::Jet> const & jets = *h_jet;

     for ( int i = 0; i < jets.size();  ++i ) {

       const reco::CATopJetTagInfo * catopTag = dynamic_cast<reco::CATopJetTagInfo const *>(jets[i].tagInfo("CATopCaloJet"));

       if ( catopTag !=0 ) {
	 char buff[1000];
	 sprintf(buff, "Jet %6d, pt = %6.2f, mass = %6.2f, minMass = %6.2f, wMass = %6.2f",
		 i,
		 jets[i].pt(), 
		 catopTag->properties().topMass,
		 catopTag->properties().minMass,
		 catopTag->properties().wMass );
	 cout << buff << endl;
       }

     }

   }

}
