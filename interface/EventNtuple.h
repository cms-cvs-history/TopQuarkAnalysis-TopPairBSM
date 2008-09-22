#ifndef EventNtuple_h
#define EventNtuple_h

/**_________________________________________________________________
   class:   EventNtuple.h
   package: TopPairBSM


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: EventNtuple.h,v 1.1 2008/08/11 06:06:21 bazterra Exp $

________________________________________________________________**/

// ROOT
#include "TObject.h"
#include "TMatrixDSym.h"

class EventNtuple : public TObject
{

public:

    EventNtuple();
    ~EventNtuple();

    void                 Reset();

	//_______ event ID_______________________________
    Int_t event;       // event number
    Int_t run;         // run number
	Int_t dataType;    // type of data: MC, cosmics, colisions,
	
    //_______ event ID_______________________________
    Int_t njets;       // number of jets
    Int_t nmuons;      // number of muons
    Int_t nvertices;   // number of vertices
  
   	//_____ total number of MC objects ______________
    Int_t ngenjets;    // number of generated jets
  
	//_____ jets ____________________________________
    std::vector< float > jet_pt;
    std::vector< float > jet_eta;
    std::vector< float > jet_phi;
    std::vector< float > jet_e;
    std::vector< float > jet_et;
    std::vector< int > jet_ntrks;
	
    std::vector< int > jet_flavour;
    std::vector< float > jetcorrection;
    
    std::vector< float > genjet_pt;
    std::vector< float > genjet_eta;
    std::vector< float > genjet_phi;
    std::vector< float > genjet_e;

	//_____ muons ____________________________________
	std::vector< double > muon_px;
	std::vector< double > muon_py;
	std::vector< double > muon_pz;
	std::vector< double > muon_e;
	std::vector< double > muon_normchi2;
	std::vector< double > muon_d0;
	std::vector< double > muon_d0Error;

    std::vector< double > genmuon_px;
	std::vector< double > genmuon_py;
	std::vector< double > genmuon_pz;
	std::vector< double > genmuon_e;
	std::vector< int > genmuon_pdg;
	std::vector< int > genmoun_motherpdg;

	std::vector< double > gentop_px;
	std::vector< double > gentop_py;
	std::vector< double > gentop_pz;
	std::vector< double > gentop_e;
	std::vector< double > gentop_charge;
	std::vector< int > gentop_hadronic;
	std::vector< double > gennu_px;
	std::vector< double > gennu_py;
	std::vector< double > gennu_pz;
	std::vector< double > gennu_e;
	std::vector< int > gennu_pdg;

    ClassDef(EventNtuple,1);

};

#endif
