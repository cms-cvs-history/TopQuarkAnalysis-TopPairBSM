#ifndef JetCombinatorics_h
#define JetCombinatorics_h

/**_________________________________________________________________
   class:   JetCombinatorics.h
   package: 


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: JetCombinatorics.h,v 1.1.4.6 2009/03/13 20:51:23 yumiceva Exp $

________________________________________________________________**/

#include "TLorentzVector.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"
#include <map>
#include <vector>
#include <iostream>

class Combo {


  public:

	Combo() {
		
		MW = 79.8;
		Mtop = 175.;
		SumEt_ = 0.;
		usebtag_ = false;
	}
	~Combo(){};

	void SetWp(TLorentzVector Wp) { Wp_ = Wp; }
	void SetWq(TLorentzVector Wq) { Wq_ = Wq; }
	void SetHadb(TLorentzVector Hadb) { Hadb_ = Hadb; }
	void SetLepW(TLorentzVector LepW) { LepW_ = LepW; }
	void SetLepb(TLorentzVector Lepb) { Lepb_ = Lepb; }
	void SetWp_disc(double disc) { Wp_disc_ = disc;}
	void SetWq_disc(double disc) { Wq_disc_= disc;}
	void SetHadb_disc(double disc) { Hadb_disc_= disc;}
	void SetLepb_disc(double disc) { Lepb_disc_= disc;}
	void SetbDiscPdf(TString filename) { 
	  pdffile_ = TFile::Open(filename);
	  hdisc_b_ = (TH1F*) gDirectory->Get("hdiscNorm_b");
	  hdisc_cl_ = (TH1F*) gDirectory->Get("hdiscNorm_cl");
	}
	void Usebtagging(bool option = true) { usebtag_ = option;}
	void SetMinMassLepW( double mass ) { minMassLepW_ = mass; }
	void SetMaxMassLepW( double mass ) { maxMassLepW_ = mass; }
	void SetMinMassHadW( double mass ) { minMassHadW_ = mass; }
	void SetMaxMassHadW( double mass ) { maxMassHadW_ = mass; }
	void SetMinMassLepTop( double mass ) { minMassLepTop_ = mass; }
	void SetMaxMassLepTop( double mass ) { maxMassLepTop_ = mass; }
	
	void analyze() {

		HadW_ = Wp_ + Wq_;
		HadTop_ = HadW_ + Hadb_;
		LepTop_ = LepW_ + Lepb_;
		TopPair_ = HadTop_ + LepTop_;

		double sigmaHadW = 2.*7.6;
		double sigmaHadt = 2.*12.5;
		double sigmaLept = 2.*15.6;
		
		double chiHadW = (HadW_.M() - MW)/sigmaHadW;
		double chiHadt = (HadTop_.M() - Mtop)/sigmaHadt;
		double chiLept = (LepTop_.M() - Mtop)/sigmaLept;

		chi2_ = chiHadW*chiHadW + chiHadt*chiHadt + chiLept*chiLept;

		SumEt_ = HadTop_.Pt();

		if ( usebtag_ ) {

			double gauss_norm = (2.)*TMath::Log(sigmaHadW*TMath::Sqrt(2*TMath::Pi())) +
				(2.)*TMath::Log(sigmaHadt*TMath::Sqrt(2*TMath::Pi())) + (2.)*TMath::Log(sigmaLept*TMath::Sqrt(2*TMath::Pi()));

			double LR_Wp; double LR_Wq;
			double LR_Hadb; double LR_Lepb;

			double LR_den = 0;
			LR_den = ( getPdfValue("cl", Wp_disc_) + getPdfValue("b", Wp_disc_));
			if (LR_den == 0 ) LR_Wp = 1e-5;
			else LR_Wp = getPdfValue( "cl", Wp_disc_ )/ LR_den;

			LR_den = ( getPdfValue("cl", Wq_disc_) + getPdfValue("b", Wq_disc_));
			if (LR_den == 0 ) LR_Wq = 1e-5;
			else LR_Wq = getPdfValue( "cl", Wq_disc_ )/ LR_den;

			LR_den = ( getPdfValue("cl", Hadb_disc_) + getPdfValue("b", Hadb_disc_));
			if (LR_den == 0 ) LR_Hadb = 1e-5;
			else LR_Hadb = getPdfValue( "b", Hadb_disc_ )/ LR_den;

			LR_den = ( getPdfValue("cl", Lepb_disc_) + getPdfValue("b", Lepb_disc_));
			if (LR_den == 0 ) LR_Lepb = 1e-5;
			else LR_Lepb = getPdfValue( "b", Lepb_disc_ )/ LR_den;

			double btag_norm = (-0.25-TMath::Log(4)/2);
			double btag_N2LL = btag_norm*4.*( LR_Wp * TMath::Log(LR_Wp/4) + LR_Wq*TMath::Log(LR_Wq/4) + LR_Hadb*TMath::Log(LR_Hadb/4) + LR_Lepb*TMath::Log(LR_Lepb/4) );
		  
			chi2_ += btag_N2LL + gauss_norm;
			
			pdffile_->Close();
		}
	}

	TLorentzVector GetWp() { return Wp_; }
	TLorentzVector GetWq() { return Wq_; }
	TLorentzVector GetHadW() { return HadW_; }
	TLorentzVector GetLepW() { return LepW_; }
	TLorentzVector GetHadb() { return Hadb_; }
	TLorentzVector GetLepb() { return Lepb_; }
	TLorentzVector GetHadTop() { return HadTop_; }
	TLorentzVector GetLepTop() { return LepTop_; }
	TLorentzVector GetTopPair() { return TopPair_; }
	double GetChi2() { return chi2_; }
	double GetSumEt() { return SumEt_; }
	int GetIdHadb() { return IdHadb_;}
	int GetIdWp() { return IdWp_; }
	int GetIdWq() { return IdWq_; }
	int GetIdLepb() { return IdLepb_;}
	void SetIdHadb(int id) { IdHadb_ = id;}
	void SetIdWp(int id) { IdWp_ = id; }
	void SetIdWq(int id) { IdWq_ = id; }
	void SetIdLepb(int id) { IdLepb_ = id;}
	void Print() {
	  std::cout << " jet Wp  : px = " << Wp_.Px() << " py = " <<  Wp_.Py() << " pz = " << Wp_.Pz() << " e = " << Wp_.E() << std::endl;
	  std::cout << " jet Wq  : px = " << Wq_.Px() << " py = " <<  Wq_.Py() << " pz = " << Wq_.Pz() << " e = "<< Wq_.E() << std::endl;
	  std::cout << " jet Hadb: px = " << Hadb_.Px() << " py = " <<  Hadb_.Py() <<" pz = " << Hadb_.Pz() <<" e = "<< Hadb_.E() << std::endl;
	  std::cout << " jet Lepb: px = " << Lepb_.Px() << " py = " <<  Lepb_.Py() <<" pz = " << Lepb_.Pz() <<" e = "<< Lepb_.E() << std::endl;
	  std::cout << " chi-squared = " << chi2_ << " sumEt = " << SumEt_ << std::endl;
	}
	double getPdfValue(std::string flavor, double disc) {
	  double pdf= 0;
	  TH1F *hpdf;
	  if ( flavor == "b" ) hpdf = hdisc_b_;
	  else hpdf = hdisc_cl_;
	  int bin = hpdf->GetXaxis()->FindBin( disc );
	  pdf = hpdf->GetBinContent( bin );
	  if ( disc < -10 || disc >50 ) return 0;
	  //if ( pdf == 0 ) return 1.e-7;
	  return pdf;
	}
	
  private:
	
	TLorentzVector Wp_;
	TLorentzVector Wq_;
	TLorentzVector HadW_;
	TLorentzVector Hadb_;
	TLorentzVector HadTop_;
	TLorentzVector LepW_;
	TLorentzVector Lepb_;	
	TLorentzVector LepTop_;
	TLorentzVector TopPair_;
	
	bool usebtag_;
	double Wp_disc_;
	double Wq_disc_;
	double Hadb_disc_;
	double Lepb_disc_;
	TFile *pdffile_;
	TH1F *hdisc_b_;
	TH1F *hdisc_cl_;

	double chi2_;
	double SumEt_;
	double minMassLepW_;
	double maxMassLepW_;
	double minMassHadW_;
	double maxMassHadW_;
	
	double minMassLepTop_;
	double maxMassLepTop_;

	double MW;
	double Mtop;

	int IdHadb_;
	int IdWp_;
	int IdWq_;
	int IdLepb_;
	
};

struct minChi2
{
  bool operator()(Combo s1, Combo s2) const
  {
    return s1.GetChi2() <= s2.GetChi2();
  }
};

struct maxSumEt
{
  bool operator()(Combo s1, Combo s2) const
  {
    return s1.GetSumEt() >= s2.GetSumEt();
  }
};


class JetCombinatorics {

  public:

	JetCombinatorics();
	~JetCombinatorics();

	void Verbose() {
	  verbosef = true;
	}

	std::map< int, std::string > Combinatorics(int k, int max = 6);
	std::map< int, std::string > NestedCombinatorics();

	void FourJetsCombinations(std::vector<TLorentzVector> jets, std::vector<double> bdiscriminators = 0);
	void SetMaxNJets(int n) { maxNJets_ = n; }
	Combo GetCombination(int n=0);
	Combo GetCombinationSumEt(int n=0);
	int GetNumberOfCombos() { return ( (int)allCombos_.size() ); } 
	//void SetCandidate( std::vector< TLorentzVector > JetCandidates );
	
	void SetLeptonicW( TLorentzVector LepW ) { theLepW_ = LepW; }

	void SetMinMassLepW( double mass ) { minMassLepW_ = mass; }
	void SetMaxMassLepW( double mass ) { maxMassLepW_ = mass; }
	void SetMinMassHadW( double mass ) { minMassHadW_ = mass; }
	void SetMaxMassHadW( double mass ) { maxMassHadW_ = mass; }
	void SetMinMassLepTop( double mass ) { minMassLepTop_ = mass; }
	void SetMaxMassLepTop( double mass ) { maxMassLepTop_ = mass; }

	void UsebTagging( bool option = true ) { UsebTagging_ = option; }
	void SetbTagPdf( TString name ) { bTagPdffilename_ = name; }
	void Clear();

	std::vector< TLorentzVector > TwoCombos();
	std::vector< TLorentzVector > ThreeCombos();

	void RemoveDuplicates( bool option) { removeDuplicates_ = option; }

	std::vector< TLorentzVector > GetComposites();
	void AnalyzeCombos();


  private:

	//int kcombos_;
	//int maxcombos_;
	bool verbosef;
	std::map< int, std::string > Template4jCombos_;
	std::map< int, std::string > Template5jCombos_;
	std::map< int, std::string > Template6jCombos_;
	std::map< int, std::string > Template7jCombos_;

	int maxNJets_;
	bool UsebTagging_;
	TString bTagPdffilename_;
	
	TLorentzVector theLepW_;

	double minMassLepW_;
	double maxMassLepW_;
	double minMassHadW_;
	double maxMassHadW_;
	double minMassLepTop_;
	double maxMassLepTop_;
	
	std::map< Combo, int, minChi2 > allCombos_;
	std::map< Combo, int, maxSumEt > allCombosSumEt_;

	Double_t minPhi_;
	double chi2_;
	int ndf_;
	bool removeDuplicates_;
	
	std::vector< TLorentzVector > cand1_;
	std::vector< TLorentzVector > cand2_;
	std::vector< TLorentzVector > cand3_;

	//int nLists_;
	
	//std::vector< TLorentzVector > composites_;
	
};

#endif
