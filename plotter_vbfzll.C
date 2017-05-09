#include <ctype.h>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm> 
#include <TROOT.h>
#include <TString.h>




#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TF1.h>
#include <TProfile.h>
#include "math.h"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/RoccoR.cc"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/RoccoR.h"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/rochcor2016.cc"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/rochcor2016.h"
//#include "EWcorr.C"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/rochcor2016.h"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/rochcor2016.cc"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/RoccoR.cc"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/RoccoR.h"
#include "/afs/cern.ch/user/g/gimandor/private/CMSSW_8_0_24/src/giulioMandorli/2016/RoccoR.cc"

Double_t erf( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]));
}
Double_t erf2( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]))+ (1.-par[0]);
}





#define SWAP2(A, B) { TLorentzVector t = A; A = B; B = t; }
void SortByEta(std::vector<TLorentzVector> &jets){
  int i, j;
	int n=jets.size();
  for (i = n - 1; i >= 0; i--){
    for (j = 0; j < i; j++){
      if (jets[j].Eta() < jets[j + 1].Eta() ){
        SWAP2( jets[j], jets[j + 1] );
		}
    }
	}
}

float getScaleFactor(TH2F *scaleMap, double pt, double eta, float sf_err, bool abs) {
//    std::cout<<"called getScaleFactor"<<std::endl;
  //  std::cout<<pt<<":, "<<eta<<std::endl;
    float sfactor = 1.0;
    int binx = scaleMap->GetXaxis()->FindBin(pt);
	 int biny;
    if (abs==0) biny = scaleMap->GetYaxis()->FindBin(eta);
    else biny = scaleMap->GetYaxis()->FindBin(TMath::Abs(eta));
  //  std::cout<<binx<<": ,"<<biny<<std::endl;
    if ( (binx != 0) && (binx != scaleMap->GetNbinsX()+1) && (biny != 0) && (biny != scaleMap->GetNbinsY()+1)) {
        sfactor = scaleMap->GetBinContent(binx, biny);
        sf_err = scaleMap->GetBinError(binx, biny);
	//		cout<<sfactor<<endl;
        if (sfactor == 0.0) {
            // bin was not filled for w/e reason, assume we don't have value in this 2D bin from the json
            sfactor = 1.0;
            sf_err = 0.0;
        }
    }
    //std::cout<<sfactor<<std::endl;
    return sfactor;
}
float getScaleFactor1D(TH1F *scaleMap, double eta, float sf_err, bool abs) {
//    std::cout<<"called getScaleFactor"<<std::endl;
  //  std::cout<<pt<<":, "<<eta<<std::endl;
    float sfactor = 1.0;
	 int binx;
    if (abs==0) binx = scaleMap->GetXaxis()->FindBin(eta);
    else binx = scaleMap->GetXaxis()->FindBin(TMath::Abs(eta));
    if ( (binx != 0) && (binx != scaleMap->GetNbinsX()+1) ) {
        sfactor = scaleMap->GetBinContent(binx);
        sf_err = scaleMap->GetBinError(binx);
	//		cout<<sfactor<<endl;
        if (sfactor == 0.0) {
            // bin was not filled for w/e reason, assume we don't have value in this 2D bin from the json
            sfactor = 1.0;
            sf_err = 0.0;
        }
    }
    //std::cout<<sfactor<<std::endl;
    return sfactor;
}







typedef std::map<double, int> JetList;
const int njets = 30;




typedef struct {
   Float_t eta[njets];
   Float_t pt[njets];
   Float_t JEC_corr[njets];
   Float_t JEC_corr_up[njets];
   Float_t JEC_corr_down[njets];
   Float_t JER_corr[njets];
   Float_t JER_corr_up[njets];
   Float_t JER_corr_down[njets];
   Float_t phi[njets];
	Float_t mass[njets];
	Float_t btag[njets];
	Int_t nsoft;
	Float_t soft_pt[njets];
	Float_t soft_eta[njets];
	Float_t soft_mass[njets];
	Float_t qgl[njets];
	Int_t nsoft2;
	Int_t nsoft5;
	Int_t nsoft10;
	Int_t EWKnsoft2;
	Int_t EWKnsoft5;
	Int_t EWKnsoft10;
	Int_t id[njets];
	Int_t puId[njets];
	Float_t HTsoft;
	Int_t partonFlavour[njets];
	Float_t EWKHTsoft;
	Float_t EWKsoft_pt[njets];
	Float_t EWKsoft_eta[njets];
	Float_t pt_regVBF[njets];	
	Float_t ptd[njets];
	Float_t axis2[njets];
	Int_t mult[njets];
	Float_t leadTrackPt[njets];
	Float_t blike_VBF[njets];
} Jets;

using namespace std;


int main(int argc, char* argv[]){

//gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.cc++");
//gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.cc++");



TString file_name = std::string(argv[1]);
TString file_tag = std::string(argv[2]);
TString region = std::string(argv[3]); 
int data = atoi(argv[4]);
int applyQCDScaleWeight = atoi(argv[5]);
TString QCDScaleWeight_str = std::string(argv[6]);
int applyJESWeight = atoi(argv[7]);
TString JESWeight_str = std::string(argv[8]);
TString heppyVersion = std::string(argv[9]);
TString postfix = std::string(argv[10]);
TString output = std::string(argv[11]);


std::map <TString, float> xsec;
std::map <TString, float> qgl_norm;
xsec["SingleMuon"] = 1.;
xsec["SingleMuonB"] = 1.;
xsec["SingleMuonC"] = 1.;
xsec["SingleMuonD"] = 1.;
xsec["SingleMuonE"] = 1.;
xsec["SingleMuonF"] = 1.;
xsec["SingleMuonG"] = 1.;
xsec["SingleElectron"] =  1.;
xsec["SingleElectronB"] =  1.;
xsec["SingleElectronC"] =  1.;
xsec["SingleElectronD"] =  1.;
xsec["SingleElectronE"] =  1.;
xsec["SingleElectronF"] =  1.;
xsec["SingleElectronG"] =  1.;
xsec["DYJetstoLL"] =  5765.4;
xsec["DYJetstoLL_amc"] =  5765.4;
xsec["DYJetstoLL_HT100"] =  5765.4;
//xsec["DYJetstoLL_HT100_200"] = 173.96106;
//xsec["DYJetstoLL_HT200_400"] = 48.27802 ;
//xsec["DYJetstoLL_HT400_600"] =6.68755 ;
//xsec["DYJetstoLL_HT600_Inf"] = 2.588804;
//xsec["DYJetstoLL_HT100_200"] = 147.40 ; 
//xsec["DYJetstoLL_HT200_400"] = 40.99 ; 
//xsec["DYJetstoLL_HT400_600"] = 5.678 ; 
//xsec["DYJetstoLL_HT600_Inf"] = 2.198; 
xsec["DYJetstoLL_HT100_200"] = 181.302; 
xsec["DYJetstoLL_HT200_400"] =50.4177  ; 
xsec["DYJetstoLL_HT400_600"] =6.98394; 
xsec["DYJetstoLL_HT600_Inf"] =2.70354 ;
xsec["DYJetstoLL_HT600_800"] = 1.6814;
xsec["DYJetstoLL_HT800_1200"] = 0.7754;
xsec["DYJetstoLL_HT1200_2500"] = 0.186;
xsec["DYJetstoLL_HT2500_Inf"] = 0.00438495;
 
xsec["DYJetstoLL_Pt-100_amc"] = 5765.4; 
xsec["DYJetstoLL_Pt-100To250_amc"] = 83.12; 
xsec["DYJetstoLL_Pt-250To400_amc"] =3.047 ; 
xsec["DYJetstoLL_Pt-400To650_amc"] = 0.3921 ; 
xsec["DYJetstoLL_Pt-650ToInf_amc"] = 0.03636 ;
 
xsec["DYJetstoLL_amc_0J"] = 4585.27; //4732.;  normalize to 5765.4 pb, 1.032 
xsec["DYJetstoLL_amc_1J"] = 853.198;//880.5; 
xsec["DYJetstoLL_amc_2J"] = 325.194;//335.6; 



xsec["VBF_HToMuMu"] = 1.;//------------------------------------------------------------------------------------------------
xsec["GluGlu_HToMuMu"] = 1.;//------------------------------------------------------------------------------------------------

xsec["TT"] =809.;
xsec["WW"] =118.7;
xsec["WZ"] = 47.13;
xsec["ZZ"] =16.523;
xsec["EWK_LLJJ"]=1.664;
xsec["EWK_LLJJ_herwig"]=1.664;
xsec["interference"]=1.664;

xsec["QCD_HT100to200"] = 27990000;
xsec["QCD_HT200to300"] = 1712000 ;
xsec["QCD_HT300to500"] = 347700;
xsec["QCD_HT500to700"] = 32100 ;
xsec["QCD_HT700to1000"] = 6831;
xsec["QCD_HT1000to1500"] = 1207 ;
xsec["QCD_HT1500to2000"] =  119.9;
xsec["QCD_HT2000toInf"] = 25.24;

xsec["ST_tW"] = 71.7 ;			//inclusive decays
xsec["ST_tW_top"] = 35.85  ;			//inclusive decays
xsec["ST_tW_antitop"] = 35.85  ;			//inclusive decays
xsec["ST_s-channel"] = 3.36; //leptonic decays
xsec["ST_t-channel_top_4f_inclusiveDecays"] = 136.02;
xsec["ST_t-channel_antitop_4f_inclusiveDecays"] = 80.95;
xsec["ST_t-channel_top_4f_leptonDecays"] = 44.33;  //leptonDecays  , multiplied with BR 0.325
xsec["ST_t-channel_antitop_4f_leptonDecays"] = 26.38;//leptonDecays ,multiplied with BR 0.325

xsec["WJetsToLNu_amc"]  = 61526.7; //not going to use these ones
xsec["WJetsToLNu"]  = 61526.7;

//k factors 1.21 are not included// 
/*xsec["WJetsToLNu_HT100To200"] = 1345 ;
xsec["WJetsToLNu_HT200To400"] = 359.7  ;
xsec["WJetsToLNu_HT400To600"] = 48.91;
xsec["WJetsToLNu_HT600To800"] =12.05;
xsec["WJetsToLNu_HT800To1200"] = 5.501;
xsec["WJetsToLNu_HT1200To2500"] = 1.329;
xsec["WJetsToLNu_HT2500ToInf"] = 0.03216;
*/
xsec["WJetsToLNu_HT100"]  = 61526.7;
xsec["WJetsToLNu_HT100To200"] = 1627.45 ;
xsec["WJetsToLNu_HT200To400"] = 435.236  ;
xsec["WJetsToLNu_HT400To600"] = 59.18109;
xsec["WJetsToLNu_HT600To800"] =14.5805;
xsec["WJetsToLNu_HT800To1200"] = 6.656210;
xsec["WJetsToLNu_HT1200To2500"] = 1.608089;
xsec["WJetsToLNu_HT2500ToInf"] = 0.0389135;

xsec["TTZToLLNuNu"] = 0.2529;
xsec["tZq_ll"]=0.0758;



if (region.CompareTo("el")==0) {
qgl_norm["EWK_LLJJ"]=0.938595977;
qgl_norm["EWK_LLJJ_herwig"]=1;
qgl_norm["TT"]=1.05998369;
qgl_norm["WW"]=0.981059169;
qgl_norm["WZ"]=0.956107275;
qgl_norm["ZZ"]=0.970928133;
qgl_norm["ST_tW_antitop"]=0.999659674;
qgl_norm["ST_tW_top"]=0.980208626;
qgl_norm["ST_s-channel"]=0.99449528;
qgl_norm["ST_t-channel_top_4f_inclusiveDecays"]=0.998267176;
qgl_norm["ST_t-channel_antitop_4f_inclusiveDecays"]=0.856539545;
qgl_norm["WJetsToLNu"]=1;
qgl_norm["DYJetstoLL_amc_0J"]=0.996136594;
qgl_norm["DYJetstoLL_amc_1J"]=0.956949008;
qgl_norm["DYJetstoLL_amc_2J"]=0.952277759;

}

if (region.CompareTo("mu")==0) {
qgl_norm["EWK_LLJJ"]=0.939774091;
qgl_norm["EWK_LLJJ_herwig"]=1;
qgl_norm["TT"]=1.069615388;
qgl_norm["WW"]=0.927930632;
qgl_norm["WZ"]=0.967820083;
qgl_norm["ZZ"]=0.964110717;
qgl_norm["ST_tW_antitop"]=1.012466766;
qgl_norm["ST_tW_top"]=0.990673372;
qgl_norm["ST_s-channel"]=0.911006075;
qgl_norm["ST_t-channel_top_4f_inclusiveDecays"]=0.986731357;
qgl_norm["ST_t-channel_antitop_4f_inclusiveDecays"]=0.994925429;
qgl_norm["WJetsToLNu"]=1;
qgl_norm["DYJetstoLL_amc_0J"]=1.002031915;
qgl_norm["DYJetstoLL_amc_1J"]=0.966710372;
qgl_norm["DYJetstoLL_amc_2J"]=0.954117783;

}


 int counter=0;


int whichQCDScaleWeight;
if ((QCDScaleWeight_str.CompareTo("none")==0)||(QCDScaleWeight_str.CompareTo("nom")==0)) whichQCDScaleWeight=0;
if (QCDScaleWeight_str.CompareTo("up")==0) whichQCDScaleWeight=1;
if (QCDScaleWeight_str.CompareTo("down")==0) whichQCDScaleWeight=2;
int whichJESWeight;
if ((JESWeight_str.CompareTo("none")==0)||(JESWeight_str.CompareTo("nom")==0)) whichJESWeight=0;
if (JESWeight_str.CompareTo("up")==0) whichJESWeight=1;
if (JESWeight_str.CompareTo("down")==0) whichJESWeight=2;

    
float gen_pos=0; 
float gen_neg=0; 
float gen_pos_weight=0; 
float gen_neg_weight=0; 


	
	Float_t presel=0;
	Float_t presel_vtype[10] = {0,0,0,0,0,0,0,0,0};
	Float_t pos_puweight=0;
	Float_t all_puweight=0.;
	Float_t puweight;
	Float_t PU=1.;
	Float_t genweight;
	Float_t bTagWeight;
	Float_t genweight0;
	float  trigWeight_tree;
	Int_t global_counter = 0;
	Int_t HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200;
	Int_t HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460;
	Int_t HLT_IsoMu22;
	Int_t HLT_IsoTkMu22;
	Int_t HLT_IsoMu27;
	Int_t HLT_IsoTkMu27;
	Int_t HLT_IsoTkMu24;
	Int_t HLT_IsoMu24;
	Int_t HLT_Ele27_eta2p1;
	TFile *file_initial;
	TChain *tree_initial;

	Int_t nvLeptons, nselLeptons;
	const int brLeptons=13;
	Float_t vLeptons_pt[30], vLeptons_eta[30], vLeptons_phi[30], vLeptons_mass[30], vLeptons_SF_IdCutLoose[30], vLeptons_SF_IdCutTight[30], vLeptons_SF_IsoLoose[30], vLeptons_SF_IsoTight[30],vLeptons_SF_trk_eta[30], vLeptons_SF_HLT_RunD4p2[30],vLeptons_SF_HLT_RunD4p3[30], vLeptons_relIso03[30], vLeptons_eleSieie[30], vLeptons_eleHoE[30], vLeptons_eleDEta[30],vLeptons_eleDPhi[30], vLeptons_eleEcalClusterIso[30], vLeptons_eleHcalClusterIso[30],vLeptons_dr03TkSumPt[30]  ;
	Int_t vLeptons_charge[30], vLeptons_pdgId[30],vLeptons_trackerLayers[30] ; 

	Float_t selLeptons_pt[30], selLeptons_eta[30], selLeptons_phi[30], selLeptons_mass[30], selLeptons_SF_IdCutLoose[30], selLeptons_SF_IdCutTight[30], selLeptons_SF_IsoLoose[30], selLeptons_SF_IsoTight[30],selLeptons_SF_trk_eta[30], selLeptons_SF_HLT_RunD4p2[30],selLeptons_SF_HLT_RunD4p3[30], selLeptons_relIso04[30], selLeptons_relIso03[30], selLeptons_eleSieie[30], selLeptons_eleHoE[30], selLeptons_eleDEta[30],selLeptons_eleDPhi[30], selLeptons_eleEcalClusterIso[30], selLeptons_eleHcalClusterIso[30],selLeptons_dr03TkSumPt[30] ;
	Int_t selLeptons_charge[30], selLeptons_pdgId[30], selLeptons_looseIdPOG[30], selLeptons_trackerLayers[30],  selLeptons_eleMVAIdSppring16GenPurp[30]; 

	TString str_leptons[brLeptons] = {"vLeptons_pt", "vLeptons_eta", "vLeptons_phi", "vLeptons_mass", "vLeptons_charge", "vLeptons_pdgId", "vLeptons_SF_IdCutLoose", "vLeptons_SF_IdCutTight", "vLeptons_SF_IsoLoose","vLeptons_SF_IsoTight","vLeptons_SF_trk_eta","vLeptons_SF_HLT_RunD4p2","vLeptons_SF_HLT_RunD4p3"};


	
////////////////////////////
//	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF.root");
//	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF");
//	TFile* file_trig_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix.root");
//	TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix");
//	TFile* file_trig_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix.root");
///	TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix");
//	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF.root");
//	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF");
/*	TFile* file_id_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunBCDEF.root");
	TH2F* id_mu_bf = (TH2F*)file_id_mu_bf->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");
	TFile* file_id_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunGH.root");
	TH2F* id_mu_aft = (TH2F*)file_id_mu_aft->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");

	TFile* file_trig_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunBCDEF.root");
	TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");
	TFile* file_trig_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunGH.root");
	TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");

	TFile* file_iso_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_LooseISO_LooseID_pt_eta_RunBCDEF.root");
	TH2F* iso_mu_bf = (TH2F*)file_iso_mu_bf->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");
	TFile* file_iso_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_LooseISO_LooseID_pt_eta_RunGH.root");
	TH2F* iso_mu_aft = (TH2F*)file_iso_mu_aft->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");

//	TFile* file_id_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_ScaleFactor_MVAIDWP80_80x.root");
//	TH2F* id_el = (TH2F*)file_id_el->Get("TriggerEffMap_ScaleFactor_MVAIDWP80_80x");
	TFile* file_tracker_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_ScaleFactor_tracker_80x.root");
	TH2F* tracker_el = (TH2F*)file_tracker_el->Get("TriggerEffMap_ScaleFactor_tracker_80x");
	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_Tight27AfterIDISO.root");
	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_Tight27AfterIDISO");
	TFile* file_id_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_EIDISO_ZH.root");
	TH2F* id_el = (TH2F*)file_id_el->Get("TriggerEffMap_EIDISO_ZH");

	TFile* file_track_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_Muons_trk_SF_RunBCDEF.root");
	TH1F* track_mu_bf = (TH1F*)file_track_mu_bf->Get("TriggerEffMap_Graph");
	TFile* file_track_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_Muons_trk_SF_RunGH.root");
	TH1F* track_mu_aft = (TH1F*)file_track_mu_aft->Get("TriggerEffMap_Graph");
*/

	RoccoR  *rc = new RoccoR("/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/mucorr/2016/rcdata.2016.v3/");


	
	file_initial = TFile::Open(file_name);
	
	tree_initial = (TChain*)file_initial->Get("tree");
	Int_t events_generated;
	TH1F *countPos;
	TH1F *countNeg;
	TH1F *countLHEScale;
	TH1F *countLHEPdf;
	TH1F *countWeighted;
	if ((data!=1)){
		countPos = (TH1F*)file_initial->Get("CountPosWeight");
		countNeg = (TH1F*)file_initial->Get("CountNegWeight");
 		countWeighted = (TH1F*)file_initial->Get("CountWeighted");
 		countLHEScale = (TH1F*)file_initial->Get("CountWeightedLHEWeightScale");
		countLHEPdf=	(TH1F*)file_initial->Get("CountWeightedLHEWeightPdf");
 	//	events_generated = countWeighted->GetBinContent(1);
 	//	if (whichQCDScaleWeight==0) events_generated = countPos->GetBinContent(1) - countNeg->GetBinContent(1);
 		if (whichQCDScaleWeight==0) 	events_generated = countWeighted->GetBinContent(1);
 		else {
			events_generated = countLHEScale->GetBinContent( countLHEScale->FindBin( 3 + whichQCDScaleWeight) );
			if (events_generated==0) events_generated =  countPos->GetBinContent(1) - countNeg->GetBinContent(1);
		}
	} else events_generated = 1;
//	if (file_tag.CompareTo("EWK_LLJJ")==0)  events_generated=events_generated/2.00703;
    Jets Jet;
    Float_t v_type;
    Float_t wrong_type=0.;
    Int_t nJets;
	Float_t JSON;	
	Float_t nPVs;	
	Float_t rho;	
	Float_t lheHT;	
	Float_t lheV_pt;	

	Float_t bdt;
	Float_t met_pt;
	Float_t met_phi;

	Jets GenHiggsSisters;

	int pos_weight_presel=0;
 	Int_t selLeptons_tightId[20];
	Float_t  selLeptons_chargedHadRelIso03[20], selLeptons_pfRelIso03[20];
	Float_t vLeptons_dz[20], vLeptons_edz[20];

	Int_t nGenVbosons;
	Float_t GenVbosons_pt[1];
	Int_t GenVbosons_pdgId[1];
	Float_t VtypeSim; 

Float_t LHE_weights_pdf_wgt[103];
Float_t LHE_weights_scale_wgt[10];
	
	float V_mass;
	ULong64_t evt;


    tree_initial->SetBranchAddress("Vtype",&v_type);
    tree_initial->SetBranchAddress("V_mass",&V_mass);
    tree_initial->SetBranchAddress("rho",&rho);
    tree_initial->SetBranchAddress("nJet",&nJets);
    tree_initial->SetBranchAddress("Jet_pt",Jet.pt);
    tree_initial->SetBranchAddress("Jet_corr_JECUp",Jet.JEC_corr_up);
    tree_initial->SetBranchAddress("Jet_corr_JECDown",Jet.JEC_corr_down);
    tree_initial->SetBranchAddress("Jet_corr",Jet.JEC_corr);
    tree_initial->SetBranchAddress("Jet_corr_JERUp",Jet.JER_corr_up);
    tree_initial->SetBranchAddress("Jet_corr_JERDown",Jet.JER_corr_down);
    tree_initial->SetBranchAddress("Jet_corr_JER",Jet.JER_corr);
    tree_initial->SetBranchAddress("Jet_eta",Jet.eta);
    tree_initial->SetBranchAddress("Jet_phi",Jet.phi);
	tree_initial->SetBranchAddress("Jet_mass",Jet.mass);
	tree_initial->SetBranchAddress("Jet_btagCSV",Jet.btag);
	tree_initial->SetBranchAddress("Jet_blike_VBF",Jet.blike_VBF);
	tree_initial->SetBranchAddress("Jet_id",Jet.id);	
	tree_initial->SetBranchAddress("Jet_puId",Jet.puId);
 	tree_initial->SetBranchAddress("Jet_leadTrackPt",Jet.leadTrackPt);
 	tree_initial->SetBranchAddress("Jet_partonFlavour",Jet.partonFlavour);
	
	tree_initial->SetBranchAddress("met_pt",&met_pt);
	tree_initial->SetBranchAddress("met_phi",&met_phi);
	
	tree_initial->SetBranchAddress("softActivityJets_pt",Jet.soft_pt);
	tree_initial->SetBranchAddress("softActivityJets_eta",Jet.soft_eta);
	tree_initial->SetBranchAddress("softActivityJets_mass",Jet.soft_mass);
	tree_initial->SetBranchAddress("softActivity_HT",&Jet.HTsoft);
	tree_initial->SetBranchAddress("softActivity_njets2",&Jet.nsoft2);
	tree_initial->SetBranchAddress("softActivity_njets5",&Jet.nsoft5);
	tree_initial->SetBranchAddress("softActivity_njets10",&Jet.nsoft10);
	tree_initial->SetBranchAddress("softActivityEWK_HT",&Jet.EWKHTsoft);
	tree_initial->SetBranchAddress("softActivityEWK_njets2",&Jet.EWKnsoft2);
	tree_initial->SetBranchAddress("softActivityEWK_njets5",&Jet.EWKnsoft5);
	tree_initial->SetBranchAddress("softActivityEWK_njets10",&Jet.EWKnsoft10);
	tree_initial->SetBranchAddress("softActivityEWKJets_pt",Jet.EWKsoft_pt);
	tree_initial->SetBranchAddress("softActivityEWKJets_eta",Jet.EWKsoft_eta);
	tree_initial->SetBranchAddress("Jet_qgl",Jet.qgl);
	tree_initial->SetBranchAddress("genWeight",&genweight);
	tree_initial->SetBranchAddress("bTagWeight",&bTagWeight);
	tree_initial->SetBranchAddress("puWeight",&puweight);
	tree_initial->SetBranchAddress("nPVs",&nPVs);
	tree_initial->SetBranchAddress("Jet_ptd",Jet.ptd);
	tree_initial->SetBranchAddress("Jet_axis2",Jet.axis2);
	tree_initial->SetBranchAddress("Jet_mult",Jet.mult);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v",&HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v",&HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu22_v",&HLT_IsoMu22);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu22_v",&HLT_IsoTkMu22);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu27_v",&HLT_IsoMu27);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu24_v",&HLT_IsoMu24);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu24_v",&HLT_IsoTkMu24);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu27_v",&HLT_IsoTkMu27);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WPTight_Gsf_v",&HLT_Ele27_eta2p1);
	tree_initial->SetBranchAddress("Jet_pt_regVBF",Jet.pt_regVBF);
	tree_initial->SetBranchAddress("json",&JSON);
    
	tree_initial->SetBranchAddress("GenHiggsSisters_pt",GenHiggsSisters.pt);
    tree_initial->SetBranchAddress("GenHiggsSisters_eta",GenHiggsSisters.eta);
    tree_initial->SetBranchAddress("GenHiggsSisters_phi",GenHiggsSisters.phi);
	tree_initial->SetBranchAddress("GenHiggsSisters_mass",GenHiggsSisters.mass);
	tree_initial->SetBranchAddress("nGenVbosons",&nGenVbosons);
	tree_initial->SetBranchAddress("GenVbosons_pt",GenVbosons_pt);
	tree_initial->SetBranchAddress("GenVbosons_pdgId",GenVbosons_pdgId);
	tree_initial->SetBranchAddress("VtypeSim",&VtypeSim);

	tree_initial->SetBranchAddress("nvLeptons",&nvLeptons);
	tree_initial->SetBranchAddress("vLeptons_pt",vLeptons_pt);
	tree_initial->SetBranchAddress("vLeptons_eta",vLeptons_eta);
	tree_initial->SetBranchAddress("vLeptons_phi",vLeptons_phi);
	tree_initial->SetBranchAddress("vLeptons_mass",vLeptons_mass);
	tree_initial->SetBranchAddress("vLeptons_charge",vLeptons_charge);
	tree_initial->SetBranchAddress("vLeptons_pdgId",vLeptons_pdgId);
	tree_initial->SetBranchAddress("vLeptons_SF_IdCutLoose",vLeptons_SF_IdCutLoose);
	tree_initial->SetBranchAddress("vLeptons_SF_IdCutTight",vLeptons_SF_IdCutTight);
	tree_initial->SetBranchAddress("vLeptons_SF_IsoLoose",vLeptons_SF_IsoLoose);
	tree_initial->SetBranchAddress("vLeptons_SF_IsoTight",vLeptons_SF_IsoTight);
	tree_initial->SetBranchAddress("vLeptons_SF_trk_eta",vLeptons_SF_trk_eta);
	tree_initial->SetBranchAddress("vLeptons_SF_HLT_RunD4p2",vLeptons_SF_HLT_RunD4p2);
	tree_initial->SetBranchAddress("vLeptons_SF_HLT_RunD4p3",vLeptons_SF_HLT_RunD4p3);
	tree_initial->SetBranchAddress("vLeptons_trackerLayers", vLeptons_trackerLayers);





	tree_initial->SetBranchAddress("nselLeptons",&nselLeptons);
	tree_initial->SetBranchAddress("selLeptons_pt",selLeptons_pt);
	tree_initial->SetBranchAddress("selLeptons_eta",selLeptons_eta);
	tree_initial->SetBranchAddress("selLeptons_phi",selLeptons_phi);
	tree_initial->SetBranchAddress("selLeptons_mass",selLeptons_mass);
	tree_initial->SetBranchAddress("selLeptons_charge",selLeptons_charge);
	tree_initial->SetBranchAddress("selLeptons_pdgId",selLeptons_pdgId);
	tree_initial->SetBranchAddress("selLeptons_looseIdPOG",selLeptons_looseIdPOG);
	tree_initial->SetBranchAddress("selLeptons_relIso04",selLeptons_relIso04);
	tree_initial->SetBranchAddress("selLeptons_relIso03",selLeptons_relIso03);
   tree_initial->SetBranchAddress("selLeptons_eleMVAIdSppring16GenPurp",selLeptons_eleMVAIdSppring16GenPurp); 
	tree_initial->SetBranchAddress("selLeptons_trackerLayers", selLeptons_trackerLayers);
	tree_initial->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	tree_initial->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	tree_initial->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	tree_initial->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	tree_initial->SetBranchAddress("selLeptons_eleHoE",selLeptons_eleHoE);
	tree_initial->SetBranchAddress("selLeptons_eleDEta",selLeptons_eleDEta);
	tree_initial->SetBranchAddress("selLeptons_eleDPhi",selLeptons_eleDPhi);
	tree_initial->SetBranchAddress("selLeptons_eleEcalClusterIso",selLeptons_eleEcalClusterIso);
	tree_initial->SetBranchAddress("selLeptons_eleHcalClusterIso",selLeptons_eleHcalClusterIso);
	tree_initial->SetBranchAddress("selLeptons_dr03TkSumPt",selLeptons_dr03TkSumPt);
	
	tree_initial->SetBranchAddress("evt",&evt);



//	for (int i=0;i<brLeptons;i++){
//		tree_initial->SetBranchAddress(str_leptons[i],);
//	}	

	tree_initial->SetBranchAddress("LHE_weights_pdf_wgt",LHE_weights_pdf_wgt);
	tree_initial->SetBranchAddress("LHE_weights_scale_wgt",LHE_weights_scale_wgt);
	tree_initial->SetBranchAddress("lheHT",&lheHT);
	tree_initial->SetBranchAddress("lheV_pt",&lheV_pt);
	tree_initial->SetBranchAddress("BDT_VBF",&bdt);



	if (data==1){
		genweight = 1.;
		bTagWeight = 1.;
		puweight=1.;
	}

 	
    TH1F *hVtype = new TH1F("hVtype","", 6,-1.,6.);
    hVtype->GetXaxis()->SetTitle("vtype");

	TH1F *hMqq = new TH1F("hMqq","",100.,0.,3000.);
	hMqq->GetXaxis()->SetTitle("m(qq) (GeV)");
	TH1F *hMqq_log = new TH1F("hMqq_log","",150.,0.,15.);
	hMqq_log->GetXaxis()->SetTitle("ln(m(qq)) (GeV)");
	TH1F *hqq_pt = new TH1F("hqq_pt","",80.,0.,800.);
	hqq_pt->GetXaxis()->SetTitle("p_{T}(qq) (GeV)");
	TH1F *hlepton1_pt = new TH1F("hlepton1_pt","",40.,0.,400.);
	hlepton1_pt->GetXaxis()->SetTitle("leading lepton p_{T} (GeV)");
	TH1F *hlepton2_pt = new TH1F("hlepton2_pt","",30.,0.,300.);
	hlepton2_pt->GetXaxis()->SetTitle("subleading lepton p_{T} (GeV)");
	TH1F *hlepton1_eta = new TH1F("hlepton1_eta","",80.,-4.,4.);
	hlepton1_eta->GetXaxis()->SetTitle("leading lepton #eta");
	TH1F *hlepton2_eta = new TH1F("hlepton2_eta","",80.,-4.,4.);
	hlepton2_eta->GetXaxis()->SetTitle("subleading lepton #eta");
  

	TH1F *hlepton1_iso03 = new TH1F("hlepton1_iso03","",40.,0.,0.4);
	hlepton1_iso03->GetXaxis()->SetTitle("leading lepton iso");
	
	TH1F *hlepton2_iso03= new TH1F("hlepton2_iso03","",40.,0.,0.4);
	hlepton2_iso03->GetXaxis()->SetTitle("subleading lepton iso");



 
	TH1F *hEtaQQ = new TH1F("hEtaQQ","",90,0.,9.);
	hEtaQQ->GetXaxis()->SetTitle("|#Delta#eta_{qq}|");
	
	TH1F *hbdt = new TH1F("hbdt","",100,-1.,1.);
	hbdt->GetXaxis()->SetTitle("BDT output");
	TH1F *hbdt_atanh = new TH1F("hbdt_atanh","",500,0.,5.);
	hbdt_atanh->GetXaxis()->SetTitle("AThanH((BDT+1)/2)");

	float bining[50];
	bining[0]=0.;
	for (int i=1;i<31;i++)
		bining[i]=bining[i-1]+0.1;
	bining[31] = 3.5;
			
	TH1F *hbdt_atanh2 = new TH1F("hbdt_atanh2","",31,bining);
	hbdt_atanh2->GetXaxis()->SetTitle("AThanH((BDT+1)/2)");

	TH1F *hPhiQQ = new TH1F("hPhiQQ","",32,0.,3.2);
	hPhiQQ->GetXaxis()->SetTitle("|#Delta#phi_{qq}|");
    
	
	TH1F *hEtaSoftJets = new TH1F("hEtaSoftJets","",12,-3.,3.);
	hEtaSoftJets->GetXaxis()->SetTitle("|#eta^{soft}|");
    
	
	TH1F *hMassSoftJets = new TH1F("hMassSoftJets","",10,0.,100.);
	hMassSoftJets->GetXaxis()->SetTitle("m^{soft}");

	TH1F *hHTsoft = new TH1F("hHTsoft","",30,0.,300.);
	hHTsoft->GetXaxis()->SetTitle("H_{T}^{soft} (GeV)" );
	TH1F *hSoft_n2 = new TH1F("hSoft_n2","",25,0.,25.);
	hSoft_n2->GetXaxis()->SetTitle("N soft jets, p_{T} > 2 GeV");
	TH1F *hSoft_n5 = new TH1F("hSoft_n5","",10,0.,10.);
	hSoft_n5->GetXaxis()->SetTitle("N soft jets, p_{T} > 5 GeV");
	TH1F *hSoft_n10 = new TH1F("hSoft_n10","",6,0.,6.);
	hSoft_n10->GetXaxis()->SetTitle("N soft jets, p_{T} > 10 GeV");

	TH1F *hHTsoftEWK = new TH1F("hHTsoftEWK","",30,0.,300.);
	hHTsoftEWK->GetXaxis()->SetTitle("EWK H_{T}^{soft} (GeV)" );
	TH1F *hSoft_n2EWK = new TH1F("hSoft_n2EWK","",25,0.,25.);
	hSoft_n2EWK->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV");
	TH1F *hSoft_n5EWK = new TH1F("hSoft_n5EWK","",10,0.,10.);
	hSoft_n5EWK->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV");
	TH1F *hSoft_n10EWK = new TH1F("hSoft_n10EWK","",6,0.,6.);
	hSoft_n10EWK->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV");
	
	TH1F *hHTsoftEWK_bdt = new TH1F("hHTsoftEWK_bdt","",30,0.,300.);
	hHTsoftEWK_bdt->GetXaxis()->SetTitle("EWK H_{T}^{soft} , BDT>0.92 (GeV)" );
	TH1F *hSoft_n2EWK_bdt = new TH1F("hSoft_n2EWK_bdt","",25,0.,25.);
	hSoft_n2EWK_bdt->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , BDT>0.92");
	TH1F *hSoft_n5EWK_bdt = new TH1F("hSoft_n5EWK_bdt","",10,0.,10.);
	hSoft_n5EWK_bdt->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , BDT>0.92");
	TH1F *hSoft_n10EWK_bdt = new TH1F("hSoft_n10EWK_bdt","",6,0.,6.);
	hSoft_n10EWK_bdt->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , BDT>0.92");


	TH1F *hHTsoftEWK_bdt2 = new TH1F("hHTsoftEWK_bdt2","",30,0.,300.);
	hHTsoftEWK_bdt2->GetXaxis()->SetTitle("EWK H_{T}^{soft} , BDT>0.84 (GeV)" );
	TH1F *hSoft_n2EWK_bdt2 = new TH1F("hSoft_n2EWK_bdt2","",25,0.,25.);
	hSoft_n2EWK_bdt2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , BDT>0.84");
	TH1F *hSoft_n5EWK_bdt2 = new TH1F("hSoft_n5EWK_bdt2","",10,0.,10.);
	hSoft_n5EWK_bdt2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , BDT>0.84");
	TH1F *hSoft_n10EWK_bdt2 = new TH1F("hSoft_n10EWK_bdt2","",6,0.,6.);
	hSoft_n10EWK_bdt2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , BDT>0.84");

	TH1F *hHTsoftEWK_mjj1 = new TH1F("hHTsoftEWK_mjj1","",30,0.,300.);
	hHTsoftEWK_mjj1->GetXaxis()->SetTitle("EWK H_{T}^{soft} , m(qq) > 1500 (GeV)" );
	TH1F *hSoft_n2EWK_mjj1 = new TH1F("hSoft_n2EWK_mjj1","",25,0.,25.);
	hSoft_n2EWK_mjj1->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , m(qq) > 1500");
	TH1F *hSoft_n5EWK_mjj1 = new TH1F("hSoft_n5EWK_mjj1","",10,0.,10.);
	hSoft_n5EWK_mjj1->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , m(qq) > 1500");
	TH1F *hSoft_n10EWK_mjj1 = new TH1F("hSoft_n10EWK_mjj1","",6,0.,6.);
	hSoft_n10EWK_mjj1->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , m(qq) > 1500");
	
	TH1F *hHTsoftEWK_mjj2 = new TH1F("hHTsoftEWK_mjj2","",30,0.,300.);
	hHTsoftEWK_mjj2->GetXaxis()->SetTitle("EWK H_{T}^{soft} , m(qq) > 2500 (GeV)" );
	TH1F *hSoft_n2EWK_mjj2 = new TH1F("hSoft_n2EWK_mjj2","",25,0.,25.);
	hSoft_n2EWK_mjj2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , m(qq) > 2500");
	TH1F *hSoft_n5EWK_mjj2 = new TH1F("hSoft_n5EWK_mjj2","",10,0.,10.);
	hSoft_n5EWK_mjj2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , m(qq) > 2500");
	TH1F *hSoft_n10EWK_mjj2 = new TH1F("hSoft_n10EWK_mjj2","",6,0.,6.);
	hSoft_n10EWK_mjj2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , m(qq) > 2500");



	TH1F* hqgl = new TH1F("hqgl","",20.,0.,1.);
	hqgl->GetXaxis()->SetTitle("QGL 1^{st} q-jet");
	
	TH1F* hqgl2 = new TH1F("hqgl2","",20.,0.,1.);
	hqgl2->GetXaxis()->SetTitle("QGL 2^{nd} q-jet");

	
	TH1F *hPtSoftJets = new TH1F("hPtSoftJets","",30,0.,300);
	hPtSoftJets->GetXaxis()->SetTitle("p_{T}^{soft} (GeV)");
    TH1F *hPtSoftJets2 = new TH1F("hPtSoftJets2", "", 20, 0., 200.);
    hPtSoftJets2->GetXaxis()->SetTitle("2nd Soft Jet p_{T} (GeV)");
    TH1F *hPtSoftJets3 = new TH1F("hPtSoftJets3", "", 20, 0., 200.);
    hPtSoftJets3->GetXaxis()->SetTitle("3rd Soft Jet p_{T} (GeV)");
	
	TH1F *hcosOqqbb = new TH1F("hcosOqqbb","",100,-1.,1.);
	hcosOqqbb->GetXaxis()->SetTitle("cos(#theta_{bb_qq})");
	TH1F *hEtaQB1 = new TH1F("hEtaQB1","",160.,-8.,8.);
	hEtaQB1->GetXaxis()->SetTitle("#Delta#eta_{qb}^{forward}");
	TH1F *hEtaQB2 = new TH1F("hEtaQB2","",160.,-8.,8.);
	hEtaQB2->GetXaxis()->SetTitle("#Delta#eta_{qb}^{backward}");
	TH1F *hPhiQB1 = new TH1F("hPhiQB1","",32,0.,3.2);
	hPhiQB1->GetXaxis()->SetTitle("#Delta#phi_{qb}^{forward}");
	TH1F *hPhiQB2 = new TH1F("hPhiQB2","",32,0.,3.2);
	hPhiQB2->GetXaxis()->SetTitle("#Delta#phi_{qb}^{backward}");
	TH1F *hx1 = new TH1F("hx1","",100.,0.,1.);
	hx1->GetXaxis()->SetTitle("x_{1}");
	TH1F *hx2 = new TH1F("hx2","",100.,0.,1.);
	hx2->GetXaxis()->SetTitle("x_{2}");
	TH1F *hVB1_mass = new TH1F("hVB1_mass","",100,0.,1000.);
	hVB1_mass->GetXaxis()->SetTitle("M_{W'_{1}} (GeV)");
	TH1F *hVB2_mass = new TH1F("hVB2_mass","",100.,0.,1000.);
	hVB2_mass->GetXaxis()->SetTitle("M_{W'_{2}} (GeV)");

	TH1F* hEtot = new TH1F("hEtot","",150.,0.,6000.);
	hEtot->GetXaxis()->SetTitle("E^{tot} (GeV)");
	TH1F* hPxtot= new TH1F("hPxtot","",100,-500.,500.);
	hPxtot->GetXaxis()->SetTitle("P_{x}^{tot} (GeV)");
	TH1F* hPytot= new TH1F("hPytot","",100,-500.,500.);
	hPytot->GetXaxis()->SetTitle("P_{y}^{tot} (GeV)");
	TH1F* hPztot= new TH1F("hPztot","",100,-5000.,5000);
	hPztot->GetXaxis()->SetTitle("P_{z}^{tot} (GeV)");

	
	TH1F *hPtqqll = new TH1F("hPtqqll","",50.,0.,500.);
	hPtqqll->GetXaxis()->SetTitle("p_{T} of qqll system (GeV)");
	TH1F *hPhiqqll = new TH1F("hPhiqqll","",32,-3.2,3.2);
	hPhiqqll->GetXaxis()->SetTitle("-#phi of qqll system");
	TH1F *hEtaqqll = new TH1F("hEtaqqll","",160,0,8);
	hEtaqqll->GetXaxis()->SetTitle("#eta of qqll system");

	TH1F *hnPVs = new TH1F("hPVs","",50,0,50);
	hnPVs->GetXaxis()->SetTitle("nPVs");
	

	TH1F* hV_mass = new TH1F("hV_mass","",20,70.,110.);
	hV_mass->GetXaxis()->SetTitle("V_mass (GeV)");

	TH1F* hZll_mass = new TH1F("hZll_mass","",20,70.,110.);
	hZll_mass->GetXaxis()->SetTitle("m(Z) (GeV)");
	TH1F* hZll_pt = new TH1F("hZll_pt","",40,0.,400.);
	hZll_pt->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
	TH1F* hZll_eta = new TH1F("hZll_eta","",20,-5,5.);
	hZll_eta->GetXaxis()->SetTitle("#eta(Z)");
	TH1F* hZll_phi = new TH1F("hZll_phi","",32,-3.2,3.2);
	hZll_phi->GetXaxis()->SetTitle("#phi(Z)");

	TH1F *hJet1q_pt = new TH1F("hJet1q_pt","",235,30,500.);
	hJet1q_pt->GetXaxis()->SetTitle("p_{T} 1^{st} q-jet");
	TH1F *hJet1q_eta = new TH1F("hJet1q_eta","",20,-5,5);
	hJet1q_eta->GetXaxis()->SetTitle("#eta 1^{st} q-jet");
	TH1F *hJet1q_ptd = new TH1F("hJet1q_ptd","",100,0,1);
	hJet1q_ptd->GetXaxis()->SetTitle("ptd 1^{st} q-jet");
	TH1F *hJet1q_axis2= new TH1F("hJet1q_axis2","",80,0.,0.16);
	hJet1q_axis2->GetXaxis()->SetTitle("#sigma_{2} 1^{st} q-jet");
	TH1F *hJet1q_mult= new TH1F("hJet1q_mult","",30,0,30.);
	hJet1q_mult->GetXaxis()->SetTitle("N 1^{st} q-jet");
	TH1F *hJet1q_leadTrackPt= new TH1F("hJet1q_leadTrackPt","",20,0,100.);
	hJet1q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 1^{st} q-jet");
	TH1F *hJet1q_leadTrackEta= new TH1F("hJet1q_leadTrackEta","",20,-5,5.);
	hJet1q_leadTrackEta->GetXaxis()->SetTitle("leading track #eta 1^{st} q-jet");
	TH1F *hJets12_pt = new TH1F("hJets12_pt","",335,30,600.);
	hJets12_pt->GetXaxis()->SetTitle("|p_{T}_{1q} + p_{T}_{2q}| (GeV)");
	TH1F *hJets12_pt_log = new TH1F("hJets12_pt_log","",100,0,10);
	hJets12_pt_log->GetXaxis()->SetTitle("ln|p_{T}_{1q} + p_{T}_{2q}| (GeV)");
	TH1F *hJet1q_phi = new TH1F("hJet1q_phi","",32,-3.2,3.2);
	hJet1q_phi->GetXaxis()->SetTitle("#phi 1^{st} q-jet");
	TH1F *hJet2q_phi = new TH1F("hJet2q_phi","",32,-3.2,3.2);
	hJet2q_phi->GetXaxis()->SetTitle("#phi 2^{nd} q-jet");


	TH1F *hJet2q_pt = new TH1F("hJet2q_pt","",235,30,500);
	hJet2q_pt->GetXaxis()->SetTitle("p_{T} 2^{nd} q-jet");
	TH1F *hJet2q_eta = new TH1F("hJet2q_eta","",20,-5,5);
	hJet2q_eta->GetXaxis()->SetTitle("#eta 2^{nd} q-jet");
	TH1F *hJet2q_ptd = new TH1F("hJet2q_ptd","",100,0,1);
	hJet2q_ptd->GetXaxis()->SetTitle("ptd 2^{nd} q-jet");
	TH1F *hJet2q_axis2= new TH1F("hJet2q_axis2","",80,0.,0.16);
	hJet2q_axis2->GetXaxis()->SetTitle("#sigma_{2} 2^{nd} q-jet");

	TH1F *hist_bins = new TH1F("bins","",80,0.,0.16);
	
	TH1F *hJet2q_pt_log = new TH1F("hJet2q_pt_log","",80,0,8);
	hJet2q_pt_log->GetXaxis()->SetTitle("log p_{T} 2^{nd} q-jet");
	TH1F *hJet1q_pt_log = new TH1F("hJet1q_pt_log","",80,0,8);
	hJet1q_pt_log->GetXaxis()->SetTitle("log p_{T} 1^{st} q-jet");
	


	TH1F *hJet2q_mult= new TH1F("hJet2q_mult","",30,0,30.);
	hJet2q_mult->GetXaxis()->SetTitle("N 2^{nd} q-jet");
	TH1F *hJet2q_leadTrackPt= new TH1F("hJet2q_leadTrackPt","",20,0,100.);
	hJet2q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 2^{nd} q-jet");
	
	TH1F *hJet3_pt = new TH1F("hJet3_pt","",18,20,200);
	hJet3_pt->GetXaxis()->SetTitle("p_{T} 3^{rd} jet");
	TH1F *hJet3_pt_new = new TH1F("hJet3_pt_new","",13,0,195);
	hJet3_pt_new->GetXaxis()->SetTitle("p_{T} 3^{rd} jet");
	TH1F *hJet3_eta = new TH1F("hJet3_eta","",20,-5,5);
	hJet3_eta->GetXaxis()->SetTitle("#eta 3^{rd} jet");
	TH1F *hJet3_eta_bdt = new TH1F("hJet3_eta_bdt","",20,-5,5);
	hJet3_eta_bdt->GetXaxis()->SetTitle("#eta 3^{rd} jet, BDT > 0.92");
	TH1F *hJet3_eta_bdt2 = new TH1F("hJet3_eta_bdt2","",20,-5,5);
	hJet3_eta_bdt2->GetXaxis()->SetTitle("#eta 3^{rd} jet, BDT > 0.84");
	
	TH1F *hsoftleadTrackPt= new TH1F("hsoftleadTrackPt","",40,0,200.);
	hsoftleadTrackPt->GetXaxis()->SetTitle("leading track p_{T}");
	TH1F *hsoftleadTrackEta= new TH1F("hsoftleadTrackEta","",20,-5,5.);
	hsoftleadTrackEta->GetXaxis()->SetTitle("leading s track #eta ");
	
	TH1F *hAdJetHT = new TH1F("hAdJetHT","",62,0,930);
	hAdJetHT->GetXaxis()->SetTitle("additional jets H_{T} (GeV)");


	TH1F *hmet = new TH1F("hmet","",40,0.,400.);
	hmet->GetXaxis()->SetTitle("MET p_{T} (GeV)");
	TH1F *hrho = new TH1F("hrho","",60,0.,30.);
	hrho->GetXaxis()->SetTitle("rho");
	TH1F *hHT = new TH1F("hHT","",50,0.,1000.);
	hHT->GetXaxis()->SetTitle("lhe H_{T} (GeV)" );
	TH1F *hlheHT_log = new TH1F("hlheHT_log","",100,0.,10.);
	hlheHT_log->GetXaxis()->SetTitle("ln(lhe H_{T}) (GeV)" );


	TH1F *hDeltaRelQQ = new TH1F("hDeltaRelQQ","",25.,0.,1.);
	hDeltaRelQQ->GetXaxis()->SetTitle("#Delta_{rel}(qq)");
	TH1F *hRptHard = new TH1F("hRptHard","",25.,0.,1.);
	hRptHard->GetXaxis()->SetTitle("R(p_{T}^{hard})");
	TH1F *hEtaQQSum = new TH1F("hEtaQQSum","",90,0.,9.);
	hEtaQQSum->GetXaxis()->SetTitle("|#Delta_{q_{1}}| + |#Delta_{q_{2}}| ");
	TH1F *hPhiZQ1 = new TH1F("hPhiZQ1","",32,0.,3.2);
	hPhiZQ1->GetXaxis()->SetTitle("|#Delta#phi(Z,q_{1})|");
	TH1F* hZll_y = new TH1F("hZll_y","",20,-5,5.);
	hZll_y->GetXaxis()->SetTitle("y(Z)");
	TH1F* hZll_ystar = new TH1F("hZll_ystar","",20,-5,5.);
	hZll_ystar->GetXaxis()->SetTitle("y*(Z)");
	TH1F* hZll_zstar = new TH1F("hZll_zstar","",15,0,3.);
	hZll_zstar->GetXaxis()->SetTitle("z*(Z)");
	TH1F* hlheV_pt = new TH1F("hlheV_pt","",40,0.,400.);
	hlheV_pt->GetXaxis()->SetTitle("lheV_pt (GeV)");


	TH1F *hJet3_pt_bdt = new TH1F("hJet3_pt_bdt","",13,0,195);
	hJet3_pt_bdt->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, BDT>0.92 (GeV)");
	TH1F *hAdJetHT_bdt = new TH1F("hAdJetHT_bdt","",30,0,450);
	hAdJetHT_bdt->GetXaxis()->SetTitle("additional jets HT, BDT>0.92 (GeV)");
	TH1F *hNAdJets_bdt = new TH1F("hNAdJets_bdt","",10,0,10);
	hNAdJets_bdt->GetXaxis()->SetTitle("N of jets, BDT>0.92");
	TH1F *hNAdJets = new TH1F("hNAdJets","",10,0,10);
	hNAdJets->GetXaxis()->SetTitle("N of jets");

	TH1F *hJet3_pt_bdt2 = new TH1F("hJet3_pt_bdt2","",13,0,195);
	hJet3_pt_bdt2->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, BDT>0.84 (GeV)");
	TH1F *hAdJetHT_bdt2 = new TH1F("hAdJetHT_bdt2","",30,0,450);
	hAdJetHT_bdt2->GetXaxis()->SetTitle("additional jets HT, BDT>0.84 (GeV)");
	TH1F *hNAdJets_bdt2 = new TH1F("hNAdJets_bdt2","",10,0,10);
	hNAdJets_bdt2->GetXaxis()->SetTitle("N of jets, BDT>0.84");
	TH1F *hJet3_pt_mjj1 = new TH1F("hJet3_pt_mjj1","",13,0,195);
	hJet3_pt_mjj1->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, m(qq) > 1500 (GeV)");
	TH1F *hAdJetHT_mjj1 = new TH1F("hAdJetHT_mjj1","",30,0,450);
	hAdJetHT_mjj1->GetXaxis()->SetTitle("additional jets HT, m(qq) > 1500 (GeV)");
	TH1F *hNAdJets_mjj1 = new TH1F("hNAdJets_mjj1","",10,0,10);
	hNAdJets_mjj1->GetXaxis()->SetTitle("N of jets, m(qq) > 1500");
	TH1F *hJet3_pt_mjj2 = new TH1F("hJet3_pt_mjj2","",13,0,195);
	hJet3_pt_mjj2->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, m(qq) > 2500 (GeV)");
	TH1F *hAdJetHT_mjj2 = new TH1F("hAdJetHT_mjj2","",30,0,450);
	hAdJetHT_mjj2->GetXaxis()->SetTitle("additional jets HT, m(qq) > 2500 (GeV)");
	TH1F *hNAdJets_mjj2 = new TH1F("hNAdJets_mjj2","",10,0,10);
	hNAdJets_mjj2->GetXaxis()->SetTitle("N of jets, m(qq) > 2500 ");


	TH1F *hJet1q_eta_bdt = new TH1F("hJet1q_eta_bdt","",20,-5,5);
	hJet1q_eta_bdt->GetXaxis()->SetTitle("#eta 1^{st} q-jet, BDT > 0.92");
	TH1F *hJet2q_eta_bdt = new TH1F("hJet2q_eta_bdt","",20,-5,5);
	hJet2q_eta_bdt->GetXaxis()->SetTitle("#eta 2^{nd} q-jet, BDT > 0.92");
	
	TH1F *hJet1q_eta_bdt2 = new TH1F("hJet1q_eta_bdt2","",20,-5,5);
	hJet1q_eta_bdt2->GetXaxis()->SetTitle("#eta 1^{st} q-jet, BDT > 0.84");
	TH1F *hJet2q_eta_bdt2 = new TH1F("hJet2q_eta_bdt2","",20,-5,5);
	hJet2q_eta_bdt2->GetXaxis()->SetTitle("#eta 2^{nd} q-jet, BDT > 0.84");


	
	TH1F *hveto_jet3pt_nom = new TH1F("hveto_jet3pt_nom","",9,0,270);
	TH1F *hveto_jet3pt_denom = new TH1F("hveto_jet3pt_denom","",9,0,270);
	TH1F *hveto_ht_nom = new TH1F("hveto_ht_nom","",14,0,420);
	TH1F *hveto_ht_denom = new TH1F("hveto_ht_denom","",14,0,420);
	TH1F *hveto_softht_nom = new TH1F("hveto_softht_nom","",8,0,320);
	TH1F *hveto_softht_denom = new TH1F("hveto_softht_denom","",8,0,320);
	TH1F *hveto_softpt_nom = new TH1F("hveto_softpt_nom","",6,0,180);
	TH1F *hveto_softpt_denom = new TH1F("hveto_softpt_denom","",6,0,180);
	

	TProfile *hprof_htsoft_pu  = new TProfile("hprof_htsoft_pu","",50,0.,50,0.,300.);
	hprof_htsoft_pu->GetXaxis()->SetTitle("number of PVs");
	hprof_htsoft_pu->GetYaxis()->SetTitle("<EWK H_{T}^{soft}> (GeV)");

	TProfile *hprof_htsoft_pu_bdt  = new TProfile("hprof_htsoft_pu_bdt","",50,0.,50,0.,300.);
	hprof_htsoft_pu_bdt->GetXaxis()->SetTitle("number of PVs");
	hprof_htsoft_pu_bdt->GetYaxis()->SetTitle("<EWK H_{T}^{soft}>, BDT > 0.92 (GeV)");

	TProfile *hprof_htsoft_pu_rms  = new TProfile("hprof_htsoft_pu_rms","",50,0.,50,0.,300.,"s");
	hprof_htsoft_pu_rms->GetXaxis()->SetTitle("number of PVs");
	hprof_htsoft_pu_rms->GetYaxis()->SetTitle("<EWK H_{T}^{soft}> (GeV)");

	TProfile *hprof_htsoft_pu_rms_bdt  = new TProfile("hprof_htsoft_pu_rms_bdt","",50,0.,50,0.,300.,"s");
	hprof_htsoft_pu_rms_bdt->GetXaxis()->SetTitle("number of PVs");
	hprof_htsoft_pu_rms_bdt->GetYaxis()->SetTitle("<EWK H_{T}^{soft}>, BDT > 0.92 (GeV)");



   		const int numArray= 109;  //64+8 
   		TH1F* histArray[numArray] = { hMqq, hEtaQQ,hHTsoft,hSoft_n2,hSoft_n5,hSoft_n10,hHTsoftEWK,hSoft_n2EWK,hSoft_n5EWK,hSoft_n10EWK,hHTsoftEWK_bdt,hSoft_n2EWK_bdt,hSoft_n5EWK_bdt,hSoft_n10EWK_bdt,hnPVs, hJet1q_pt, hJet1q_eta, hJet1q_ptd, hJet1q_axis2, hJet1q_mult, hJet2q_pt, hJet2q_eta, hJet2q_ptd, hJet2q_axis2, hJet2q_mult, hmet,   hJet1q_leadTrackPt, hJet2q_leadTrackPt, hqq_pt,hV_mass, hqgl, hqgl2, hZll_mass, hZll_pt, hZll_phi, hZll_eta, hrho, hlepton1_pt, hlepton2_pt, hlepton1_eta, hlepton2_eta, hHT, hDeltaRelQQ, hRptHard, hEtaQQSum, hPhiZQ1, hZll_y, hZll_ystar, hZll_zstar, hMqq_log, hlheV_pt, hJet3_pt, hlheHT_log, hPhiQQ, hJets12_pt_log, hJets12_pt, hJet1q_pt_log, hJet2q_pt_log, hbdt, hbdt_atanh,hbdt_atanh2 , hlepton1_iso03, hlepton2_iso03, hveto_jet3pt_nom, hveto_jet3pt_denom, hveto_ht_nom, hveto_ht_denom, hveto_softht_nom, hveto_softht_denom, hveto_softpt_nom, hveto_softpt_denom, hJet2q_phi, hJet1q_phi, hNAdJets, hNAdJets_bdt, hJet3_pt_bdt, hAdJetHT_bdt, hNAdJets_bdt2, hJet3_pt_bdt2, hAdJetHT_bdt2,hNAdJets_mjj1, hJet3_pt_mjj1, hAdJetHT_mjj1,hNAdJets_mjj2, hJet3_pt_mjj2, hAdJetHT_mjj2, hHTsoftEWK_bdt2,hSoft_n2EWK_bdt2,hSoft_n5EWK_bdt2,hSoft_n10EWK_bdt2,hHTsoftEWK_mjj1,hSoft_n2EWK_mjj1,hSoft_n5EWK_mjj1,hSoft_n10EWK_mjj1 ,hHTsoftEWK_mjj2,hSoft_n2EWK_mjj2,hSoft_n5EWK_mjj2,hSoft_n10EWK_mjj2 ,hJet1q_eta_bdt, hJet1q_eta_bdt2, hJet2q_eta_bdt, hJet2q_eta_bdt2,
	hsoftleadTrackPt, hsoftleadTrackEta, hAdJetHT, hJet3_eta , hJet3_pt_new , hJet3_eta_bdt, hJet3_eta_bdt2};
			for (int i=0;i<numArray;i++){
				histArray[i]->Sumw2();
			}
			hprof_htsoft_pu->Sumw2();
			hprof_htsoft_pu_bdt->Sumw2();
			hprof_htsoft_pu_rms->Sumw2();
			hprof_htsoft_pu_rms_bdt->Sumw2();
	

		
		TString cut_flow_names[30] = {"triggers","2jets events","q1_pt>50","q2_pt>30","Mqq>200","leptons_pt<20","(mll-mz)<15"};
		Float_t cut_flow[30] = {0,0,0,0,0,0,0};

	float qq_matching = 0;
	float qq_matching_all = 0;
	

	int nentries = tree_initial->GetEntries() ;
	

//	TF1 *func_lheHT = new TF1("func_lheHT","([0]+[1]*x+[2]*x*x+[3]*x*x*x)*TMath::Exp(-1*[4]*x)",60,4000);
//	func_lheHT->FixParameter(0,  -1.71063e+00);
//	func_lheHT->FixParameter(1, 6.90159e-02 );
///	func_lheHT->FixParameter(2, -2.83168e-04);
//	func_lheHT->FixParameter(3, 4.69007e-07);
//	func_lheHT->FixParameter(4, 7.07950e-03 );

	TF1* func_lheHT = new TF1("func_lheHT","pol6",4.2,7.8);
	func_lheHT->FixParameter(0,89.0139);
	func_lheHT->FixParameter(1,-275.535);
	func_lheHT->FixParameter(2,195.308);
	func_lheHT->FixParameter(3,-61.2467);
	func_lheHT->FixParameter(4,9.8217);
	func_lheHT->FixParameter(5,-0.791744);
	func_lheHT->FixParameter(6,0.0255211);
//	TF1* func_EtaQQ = new TF1("func_EtaQQ","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.,16.);
//	func_EtaQQ->FixParameter(0,1.33558e+00);
//	func_EtaQQ->FixParameter(1,-2.60070e-01 );
//	func_EtaQQ->FixParameter(2,5.87759e-02);
//	func_EtaQQ->FixParameter(3,-4.64109e-03 );

	TF1* func_JetsPt = new TF1("func_JetsPt","pol5",4.3,10);
	TF1* func_EtaQQ = new TF1("func_EtaQQ","pol7",0,10);
	TF1* func_Mqq = new TF1("func_Mqq","pol6",0,10);

	func_Mqq->FixParameter(0,-1378.361758);
	func_Mqq->FixParameter(1,1327.383178);
	func_Mqq->FixParameter(2,-529.261759);
	func_Mqq->FixParameter(3,111.854413);
	func_Mqq->FixParameter(4,-13.208941);
	func_Mqq->FixParameter(5,0.826140);
	func_Mqq->FixParameter(6,-0.021376);


	TF1* func_qgl_q = new TF1("func_qgl_q","pol3",0.,1.);
	func_qgl_q->FixParameter(0,0.981581);
	func_qgl_q->FixParameter(1,-0.255505);
	func_qgl_q->FixParameter(2,0.929524);
	func_qgl_q->FixParameter(3,-0.666978);
	TF1* func_qgl_g = new TF1("func_qgl_g","pol7",0.,1.);
	func_qgl_g->FixParameter(0,0.612992);
	func_qgl_g->FixParameter(1,6.27);
	func_qgl_g->FixParameter(2,-34.3663);
	func_qgl_g->FixParameter(3,92.8668);
	func_qgl_g->FixParameter(4,-99.927);
	func_qgl_g->FixParameter(5,-21.1421);
	func_qgl_g->FixParameter(6, 113.218);
	func_qgl_g->FixParameter(7,-55.7067);
//pythia8 quark (|eta|<2.0, pT inclusive, pythia ): -0.666978*x*x*x + 0.929524*x*x -0.255505*x + 0.981581
//pythia8 gluon (|eta|<2.0, pT inclusive, pythia ): -55.7067*x^7 + 113.218*x^6 -21.1421*x^5 -99.927*x^4 + 92.8668*x^3 -34.3663*x^2 + 6.27*x + 0.612992

	TF1* interference_func = new TF1("interference_func","pol7",0,10);
	interference_func->FixParameter(0,-3236.73);
	interference_func->FixParameter(1,3158.71);
	interference_func->FixParameter(2,-1314.93);
	interference_func->FixParameter(3,302.849);
	interference_func->FixParameter(4,-41.6913);
	interference_func->FixParameter(5,3.4312);
	interference_func->FixParameter(6,-0.156337);
	interference_func->FixParameter(7,0.00304253);

//	rochcor2016 *rmcor = new rochcor2016();


//	cout<<nentries<<endl;
	for (int entry=0; entry<nentries;++entry){
//	for (int entry=0; entry<2000;++entry){
//        if (entry%1000 == 0) std::cout << "Writing " << entry << "th event" << std::endl;
        tree_initial->GetEntry(entry);
	

		for (int i=0;i<nJets;i++){
			if (whichJESWeight==1) Jet.pt[i]*= Jet.JEC_corr_up[i]/Jet.JEC_corr[i]; 	
			if (whichJESWeight==2) Jet.pt[i]*= Jet.JEC_corr_down[i]/Jet.JEC_corr[i];	
		}

		if (JSON!=1) {
			continue;
		}

	//	if ((file_tag.CompareTo("EWK_LLJJ")==0) &&(evt%2!=0)) continue;

	//	if (region.CompareTo("mu")==0) if (!(v_type==0)) continue;
	//	if (region.CompareTo("el")==0) if (!(v_type==1)) continue;

			

		if (data==1) PU=1.;
		else PU=puweight;
		genweight0 = genweight/TMath::Abs(genweight);
		genweight=genweight/TMath::Abs(genweight)*PU;   
		genweight/=events_generated/xsec[file_tag]; 
		if  ((data!=1 ) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)&& (file_tag.CompareTo("ST_tW_top")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)) if (whichQCDScaleWeight==1)  genweight*=LHE_weights_scale_wgt[4];
		if  ((data!=1 ) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)) if (whichQCDScaleWeight==2)  genweight*=LHE_weights_scale_wgt[5];

		
		double GenVbosons_pt_first = GenVbosons_pt[0];
		int GenVbosons_pdgId_first = GenVbosons_pdgId[0];


 
		if  ((file_tag.CompareTo("DYJetstoLL_HT100")==0)) if (lheHT>100) continue;  
		if  ((file_tag.CompareTo("DYJetstoLL_Pt-100_amc")==0)) if (lheV_pt>100) continue;  
		if  ((file_tag.CompareTo("WJetsToLNu_HT100")==0)) if (lheHT>100) continue;  

		int pt_num1 = -1;
		int pt_num2 = -1;
		TLorentzVector Qjet1;
		TLorentzVector Qjet2;
		TLorentzVector qq;
		int good_jets = 0;
		vector<TLorentzVector> jets_pv;
		vector<int> jets_indices;
		///////////////////////
		//preselection/////
		//////////////////////
		cut_flow[0]+=genweight;

		for (int i=0;i<nJets;i++){
			TLorentzVector jet0;
		//	if (!((Jet.id[i]>2)&&(Jet.puId[i]>0)&&(Jet.pt[i]>20))) continue;
			if (!((Jet.id[i]>2)&&(Jet.puId[i]>0))) continue;
			jet0.SetPtEtaPhiM(Jet.pt[i],Jet.eta[i],Jet.phi[i],Jet.mass[i]);
			jets_pv.push_back(jet0);
			jets_indices.push_back(i);
			good_jets++;
		}
		if (good_jets<2) continue;
		
		Qjet1 = jets_pv[0];
		Qjet2 = jets_pv[1];
		float jet3_pt = 0;
		float jet3_eta;
		if (good_jets>=3) {
			jet3_pt=jets_pv[2].Pt();
			jet3_eta=jets_pv[2].Eta();
		}
		qq=Qjet1+Qjet2;
		Float_t Mqq = qq.M();
		Float_t qq_pt = qq.Pt();
		Float_t qqDeltaEta = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
		Float_t qqDeltaPhi = TMath::Abs(Qjet1.DeltaPhi(Qjet2));
//////////////////leptons////////////////
		TLorentzVector lepton1;
		TLorentzVector lepton2;
		TLorentzVector Zll;
		int idx_1stLepton = 0;
		int idx_2ndLepton = 0;
		int count_l=0;
		if (region.CompareTo("el")==0) {
		for (int i=0; i<nselLeptons;i++ ){
			if (!((selLeptons_eleMVAIdSppring16GenPurp[i]>=2)&& (selLeptons_relIso03[i]<0.15)&& (TMath::Abs(selLeptons_pdgId[i])==11))) continue;
			
		//	if  (!(   (selLeptons_pt[i]>15) && 
//		( ( (TMath::Abs(selLeptons_eta[i]) < 1.4442) &&  (selLeptons_eleSieie[i] < 0.012) && (selLeptons_eleHoE[i] < 0.09) && ((selLeptons_eleEcalClusterIso[i]/selLeptons_pt[i] ) < 0.37) && ((selLeptons_eleHcalClusterIso[i]/selLeptons_pt[i]) < 0.25) && ((selLeptons_dr03TkSumPt[i]/selLeptons_pt[i]) < 0.18) && (TMath::Abs(selLeptons_eleDEta[i]) < 0.0095)  && (TMath::Abs(selLeptons_eleDPhi[i]) < 0.065 ) )  || 
//		  ( TMath::Abs(selLeptons_eta[i]) > 1.5660) && (selLeptons_eleSieie[i]< 0.033) && (selLeptons_eleHoE[i] <0.09)&& ((selLeptons_eleEcalClusterIso[i]/selLeptons_pt[i] ) < 0.45) && ((selLeptons_eleHcalClusterIso[i]/selLeptons_pt[i]) < 0.28) && ((selLeptons_dr03TkSumPt[i]/selLeptons_pt[i]) < 0.18) ) )  )   continue;


	//		if ((count_l==1) && (selLeptons_charge[idx_1stLepton]*selLeptons_charge[i] > 0)) continue;
			if (count_l==1) {
				idx_2ndLepton=i;
				lepton2.SetPtEtaPhiM(selLeptons_pt[idx_2ndLepton], selLeptons_eta[idx_2ndLepton], selLeptons_phi[idx_2ndLepton], selLeptons_mass[idx_2ndLepton]);
				count_l++;
				break;
			}
			if (count_l==0) {
				idx_1stLepton=i;
				lepton1.SetPtEtaPhiM(selLeptons_pt[idx_1stLepton], selLeptons_eta[idx_1stLepton], selLeptons_phi[idx_1stLepton], selLeptons_mass[idx_1stLepton]);
				count_l++;
			}
		}
		}
		if (region.CompareTo("mu")==0) {
		count_l=0;
		idx_1stLepton = 0;
		idx_2ndLepton = 0;
		for (int i=0; i<nselLeptons;i++ ){
			if (!((selLeptons_looseIdPOG[i]>0) && (selLeptons_relIso04[i]<0.25) && (TMath::Abs(selLeptons_pdgId[i])==13 )) ) continue;
		//	if ((count_l==1) && (selLeptons_charge[idx_1stLepton]*selLeptons_charge[i] > 0)) continue;
			if (count_l==1) {
				idx_2ndLepton=i;
				lepton2.SetPtEtaPhiM(selLeptons_pt[idx_2ndLepton], selLeptons_eta[idx_2ndLepton], selLeptons_phi[idx_2ndLepton], selLeptons_mass[idx_2ndLepton]);
				count_l++;
				break;
			}
			if (count_l==0) {
				idx_1stLepton=i;
				lepton1.SetPtEtaPhiM(selLeptons_pt[idx_1stLepton], selLeptons_eta[idx_1stLepton], selLeptons_phi[idx_1stLepton], selLeptons_mass[idx_1stLepton]);
				count_l++;
			}
		}
		}
		if (count_l<2)  continue;
		if ((selLeptons_charge[idx_1stLepton]*selLeptons_charge[idx_2ndLepton]) >0) continue;

/*
		lepton1.SetPtEtaPhiM(vLeptons_pt[0], vLeptons_eta[0], vLeptons_phi[0], vLeptons_mass[0]);	
		int idx_2ndLepton = 0;
		for (int i=1; i<nvLeptons;i++ ){
			if (vLeptons_charge[0]*vLeptons_charge[i] < 0) {
				idx_2ndLepton=i;
				break;
			}
		}
		lepton2.SetPtEtaPhiM(vLeptons_pt[idx_2ndLepton], vLeptons_eta[idx_2ndLepton], vLeptons_phi[idx_2ndLepton], vLeptons_mass[idx_2ndLepton]);

		float qter1 = 1.0;
		float qter2 = 1.0;
		float mu_correction1 = 1.0;
		float mu_correction2 = 1.0;
		if (data!=1) 	if (region.CompareTo("mu")==0) {
			rmcor->momcor_mc(lepton1, vLeptons_charge[0], vLeptons_trackerLayers[0], qter1);
			rmcor->momcor_mc(lepton2, vLeptons_charge[idx_2ndLepton], vLeptons_trackerLayers[idx_2ndLepton], qter2);
			}
		if (data==1) 	if (region.CompareTo("mu")==0){
			rmcor->momcor_data(lepton1, vLeptons_charge[0],  0, qter1);
			rmcor->momcor_data(lepton2, vLeptons_charge[idx_2ndLepton], 0, qter2);
			}
*/
///////////////muon corrections 2016 calibration////////////////
/*		if (region.CompareTo("mu")==0) {
			double dataSF1, dataSF2 ;
			double mcSF1, mcSF2;
			if (data==1) {
				dataSF1 = rc->kScaleDT(selLeptons_charge[idx_1stLepton], lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), 0, 0);
				dataSF2 = rc->kScaleDT(selLeptons_charge[idx_2ndLepton], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), 0, 0);
				lepton1.SetPtEtaPhiM(lepton1.Pt()*dataSF1,lepton1.Eta(), lepton1.Phi(), lepton1.M() );
				lepton2.SetPtEtaPhiM(lepton2.Pt()*dataSF2,lepton2.Eta(), lepton2.Phi(), lepton2.M() );
			}
			if (data!=1) {
				double u1 = gRandom->Rndm();
				double u2 = gRandom->Rndm();
				mcSF1 = rc->kScaleAndSmearMC(selLeptons_charge[idx_1stLepton], lepton1.Pt(), lepton1.Eta(), lepton1.Phi(),  selLeptons_trackerLayers[idx_1stLepton], u1, u2, 0, 0);
				mcSF2 = rc->kScaleAndSmearMC(selLeptons_charge[idx_2ndLepton], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(),  selLeptons_trackerLayers[idx_2ndLepton], u1, u2, 0, 0);
				lepton1.SetPtEtaPhiM(lepton1.Pt()*mcSF1,lepton1.Eta(), lepton1.Phi(), lepton1.M() );
				lepton2.SetPtEtaPhiM(lepton2.Pt()*mcSF2,lepton2.Eta(), lepton2.Phi(), lepton2.M() );
			}
		}
		
*/
////////////////////////////////////////////////////////////////
			if  (region.CompareTo("mu")==0) if (!((HLT_IsoMu24==1) || (HLT_IsoTkMu24==1)  )) continue; 
			if  (region.CompareTo("el")==0) if (!(HLT_Ele27_eta2p1 == 1)) continue;

/*		if (data!=1) {
			if (region.CompareTo("mu")==0) {
				float SF_mu_bf_err1 = 0.;
				float SF_mu_bf_err2 = 0.;
				float SF_mu_aft_err1 = 0.;
				float SF_mu_aft_err2 = 0.;
				bool abs=1;
				float eff1 =20.1/36.4*getScaleFactor(trig_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
	
				float eff1_id =20.1/36.4*getScaleFactor(id_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(id_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2_id =20.1/36.4*getScaleFactor(id_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ;  
				float eff1_iso =20.1/36.4*getScaleFactor(iso_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(id_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2_iso =20.1/36.4*getScaleFactor(iso_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ; 
				abs=0; 
				float eff1_tracker =20.1/36.4*getScaleFactor1D(track_mu_bf, lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor1D(track_mu_aft,  lepton1.Eta(), SF_mu_bf_err1,abs );  	
				float eff2_tracker =20.1/36.4*getScaleFactor1D(track_mu_bf, lepton2.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor1D(track_mu_aft,  lepton2.Eta(), SF_mu_bf_err1,abs );  	

				genweight*= eff1*eff1_id*eff2_id*eff1_iso*eff2_iso*eff1_tracker*eff2_tracker; 	
			}
			if (region.CompareTo("el")==0) {
				float SF_el_err1 = 0.;
				float SF_el_err2 = 0.;
				bool abs=0;
				float eff1 = getScaleFactor(trig_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs );  	
				float eff1_id =getScaleFactor(id_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs ) ;  	
				float eff2_id =getScaleFactor(id_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs  ) ;  
				float eff1_tracker =getScaleFactor(tracker_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs ) ;  	
				float eff2_tracker =getScaleFactor(tracker_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs  ) ;  
				genweight*=eff1* eff1_id*eff2_id*eff1_tracker*eff2_tracker; 	
			//	genweight*= eff1_id*eff2_id*eff1_tracker*eff2_tracker; 	
			}
		}
*/
		string file_tag_str = file_tag.Data();


		Float_t Mqq_log = TMath::Log(Mqq);	

	//	if ((Mqq_log< 8.3227 )&&(Mqq_log>5.101)) if ( (file_tag_str.find("DYJetstoLL_HT")!=std::string::npos) || (file_tag.CompareTo("DYJetstoLL")==0))  genweight*=func_Mqq->Eval(Mqq_log);		
		if ( (file_tag_str.find("DYJetstoLL_HT")!=std::string::npos) || (file_tag.CompareTo("DYJetstoLL")==0))  genweight*=func_Mqq->Eval(Mqq_log);		


		float qgl_weight=1.;
		int apply_qgl=0;
		if (!( (data==1)|| (Jet.partonFlavour[jets_indices[0]] ==0 ) || (TMath::Abs(Jet.eta[jets_indices[0]])>=2) || (Jet.qgl[jets_indices[0]] < 0) ) ) {
			if (TMath::Abs(Jet.partonFlavour[jets_indices[0]]) < 4 ) qgl_weight=func_qgl_q->Eval(Jet.qgl[jets_indices[0]]);
			if (TMath::Abs(Jet.partonFlavour[jets_indices[0]]) ==21 ) qgl_weight=func_qgl_g->Eval(Jet.qgl[jets_indices[0]]);
		}
		if (qgl_weight!=1.) apply_qgl+=1;
	//	cout<<qgl_weight<<endl;
		genweight*=qgl_weight;
		qgl_weight=1.;
		if (!( (data==1)|| (Jet.partonFlavour[jets_indices[1]] ==0 ) || (TMath::Abs(Jet.eta[jets_indices[1]])>=2) || (Jet.qgl[jets_indices[1]] < 0)) ) {
			if (TMath::Abs(Jet.partonFlavour[jets_indices[1]]) < 4 ) qgl_weight=func_qgl_q->Eval(Jet.qgl[jets_indices[1]]);
			if (TMath::Abs(Jet.partonFlavour[jets_indices[1]]) ==21 ) qgl_weight=func_qgl_g->Eval(Jet.qgl[jets_indices[1]]);
		}
		if (qgl_weight!=1.) apply_qgl+=1;
		if (data!=1) genweight*=qgl_norm[file_tag];
	//	cout<<qgl_weight<<endl;
		genweight*=qgl_weight;
			


	
		Zll = lepton1+lepton2;	
		Float_t Zll_mass = Zll.M();
		Float_t Zll_pt = Zll.Pt();
		Float_t Zll_eta = Zll.Eta();
		Float_t Zll_phi = Zll.Phi();


		cut_flow[1]+=genweight;
		if (Qjet1.Pt() < 45) continue;
		cut_flow[2]+=genweight;
		if (Qjet2.Pt() < 25) continue;
		cut_flow[3]+=genweight;
		if (Mqq<200) continue;
		cut_flow[4]+=genweight;
	
	 	if (region.CompareTo("mu")==0) {
			if (lepton1.Pt()<30) continue;	
			if (lepton2.Pt()<20) continue;
		}	
	 	if (region.CompareTo("el")==0) {
			if (lepton1.Pt()<30) continue;	
			if (lepton2.Pt()<20) continue;
		}	
		cut_flow[5]+=genweight;
	 	if (region.CompareTo("el")==0) {
			if (TMath::Abs(lepton1.Eta())>2.1) continue;	
			if (TMath::Abs(lepton2.Eta())>2.1) continue;
		}	
	 	if (region.CompareTo("mu")==0) {
			if (TMath::Abs(lepton1.Eta())>2.4) continue;	
			if (TMath::Abs(lepton2.Eta())>2.4) continue;
		}	
		if (Zll_mass < 105 ) continue;
		cut_flow[6]+=genweight;

//		cout<<genweight<<endl;

		presel+=genweight;

		presel_vtype[(int)(v_type+1)]+=genweight;

		float Zll_ystar = Zll.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity()) ;
		float Zll_zstar = TMath::Abs( Zll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() )) ;
		Float_t DeltaEtaQQSum = TMath::Abs(Qjet1.Eta()) +  TMath::Abs(Qjet2.Eta());
		Float_t PhiZQ1 = TMath::Abs(Zll.DeltaPhi(Qjet1));
		Float_t DeltaRelQQ = (Qjet1+Qjet2).Pt()/( Qjet1.Pt()+Qjet2.Pt()); 
		Float_t RptHard = (Qjet1+Qjet2+ Zll).Pt()/( Qjet1.Pt()+Qjet2.Pt() + Zll.Pt()); 



/////////////////////////////////////////////////////////


			counter++;

            hVtype->Fill(v_type,genweight);
		   	hMqq->Fill(Mqq,genweight);
		   	hMqq_log->Fill(TMath::Log(Mqq),genweight);
			   hqq_pt->Fill(qq_pt,genweight);
			   hEtaQQ->Fill(qqDeltaEta,genweight);
		 	   hPhiQQ->Fill(qqDeltaPhi,genweight);
		   	hZll_mass->Fill(Zll_mass,genweight);
		   	hZll_pt->Fill(Zll_pt,genweight);
		   	hZll_phi->Fill(Zll_phi,genweight);
		   	hZll_eta->Fill(Zll_eta,genweight);
			   hHTsoft->Fill(Jet.HTsoft,genweight);
			   hSoft_n2->Fill(Jet.nsoft2, genweight);
			   hSoft_n5->Fill(Jet.nsoft5, genweight);
			   hSoft_n10->Fill(Jet.nsoft10, genweight);
			   hHTsoftEWK->Fill(Jet.EWKHTsoft,genweight);
			   hSoft_n2EWK->Fill(Jet.EWKnsoft2, genweight);
			   hSoft_n5EWK->Fill(Jet.EWKnsoft5, genweight);
			   hSoft_n10EWK->Fill(Jet.EWKnsoft10, genweight);
				hnPVs->Fill(nPVs,genweight);
				hqgl->Fill(Jet.qgl[jets_indices[0]],genweight);
				hqgl2->Fill(Jet.qgl[jets_indices[1]],genweight);
			
				hJet1q_pt->Fill(jets_pv[0].Pt(),genweight);
				hJet1q_eta->Fill(jets_pv[0].Eta(),genweight);
				hJet1q_ptd->Fill(Jet.ptd[jets_indices[0]],genweight);
				hJet1q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[0]]),genweight);
				hJet1q_mult->Fill(Jet.mult[jets_indices[0]],genweight);
				hJet1q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[0]],genweight);
				hJet1q_phi->Fill(jets_pv[0].Phi(),genweight);
				hJet2q_pt->Fill(jets_pv[1].Pt(),genweight);
				hJet2q_eta->Fill(jets_pv[1].Eta(),genweight);
				hJet2q_ptd->Fill(Jet.ptd[jets_indices[1]],genweight);
				hJet2q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[1]]),genweight);
				hJet2q_mult->Fill(Jet.mult[jets_indices[1]],genweight);
				hJet2q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[1]],genweight);
				hJet2q_phi->Fill(jets_pv[1].Phi(),genweight);
				hmet->Fill(met_pt,genweight);
				hV_mass->Fill(V_mass,genweight);
				hrho->Fill(rho,genweight);
				hlepton1_pt->Fill(lepton1.Pt(),genweight);
				hlepton2_pt->Fill(lepton2.Pt(),genweight);
				hlepton1_eta->Fill(lepton1.Eta(),genweight);
				hlepton2_eta->Fill(lepton2.Eta(),genweight);
				hlepton1_iso03->Fill(selLeptons_relIso03[idx_1stLepton],genweight);
				hlepton2_iso03->Fill(selLeptons_relIso03[idx_2ndLepton],genweight);
				hHT->Fill(lheHT ,genweight);
				hlheHT_log->Fill(TMath::Log(lheHT) ,genweight);
				
				hDeltaRelQQ->Fill(DeltaRelQQ,genweight);
				hRptHard->Fill(RptHard,genweight);
				hEtaQQSum->Fill(DeltaEtaQQSum,genweight);
				hPhiZQ1->Fill(PhiZQ1,genweight);
				hZll_y->Fill(Zll.Rapidity(),genweight);
				hZll_ystar->Fill(Zll_ystar   ,genweight);
				hZll_zstar->Fill(Zll_zstar,genweight);
				hlheV_pt->Fill(lheV_pt  ,genweight);
				hJet3_pt->Fill(jet3_pt ,genweight);	
				if (good_jets>=3) hJet3_eta->Fill(jet3_eta ,genweight);	
				if (good_jets>=3) hJet3_pt_new->Fill(jets_pv[2].Pt(),genweight);
				if (good_jets==2) hJet3_pt_new->Fill(10.,genweight);
				float AdJetHT = 0;
				if (good_jets>=3)
					for (int i=2;i<good_jets;i++)
						if (jets_pv[i].Pt() > 15 ) AdJetHT+=jets_pv[i].Pt();
				if (good_jets==2) hAdJetHT->Fill(10.,genweight);
				if (good_jets>=3) hAdJetHT->Fill(AdJetHT,genweight);
				hsoftleadTrackPt->Fill(Jet.EWKsoft_pt[0],genweight);	
				hsoftleadTrackEta->Fill(Jet.EWKsoft_eta[0],genweight);	
			
				hJets12_pt->Fill((jets_pv[0].Pt() + jets_pv[1].Pt()),genweight);
				hJets12_pt_log->Fill(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()),genweight);
				hJet1q_pt_log->Fill(TMath::Log(jets_pv[0].Pt()),genweight);
				hJet2q_pt_log->Fill(TMath::Log(jets_pv[1].Pt()),genweight);
				hNAdJets->Fill(good_jets, genweight);	

				if  (file_tag.CompareTo("interference")!=0) {
					hbdt->Fill(bdt,genweight);
					hbdt_atanh->Fill(TMath::ATanH((bdt+1)/2),genweight);
					hbdt_atanh2->Fill(TMath::ATanH((bdt+1)/2),genweight);
				}
				else {
					float interference_weight= interference_func->Eval(TMath::Log(Mqq))  ;
					hbdt->Fill(bdt,genweight*interference_weight);
					hbdt_atanh->Fill(TMath::ATanH((bdt+1)/2),genweight*interference_weight);
					hbdt_atanh2->Fill(TMath::ATanH((bdt+1)/2),genweight*interference_weight);
				}
		 

				hprof_htsoft_pu->Fill(nPVs,Jet.EWKHTsoft,genweight);
				hprof_htsoft_pu_rms->Fill(nPVs,Jet.EWKHTsoft,genweight);

				if (bdt>0.92) {
					hprof_htsoft_pu_bdt->Fill(nPVs,Jet.EWKHTsoft,genweight);
					hprof_htsoft_pu_rms_bdt->Fill(nPVs,Jet.EWKHTsoft,genweight);
				
				//	hveto_jet3pt_denom->Fill(,genweight)
					for (int i=0;i<hveto_jet3pt_denom->GetNbinsX();i++)
						hveto_jet3pt_denom->Fill(hveto_jet3pt_denom->GetBinCenter(i+1),genweight);
					for (int i=0;i<hveto_ht_denom->GetNbinsX();i++)
						hveto_ht_denom->Fill(hveto_ht_denom->GetBinCenter(i+1),genweight);
					for (int i=0;i<hveto_softht_denom->GetNbinsX();i++)
						hveto_softht_denom->Fill(hveto_softht_denom->GetBinCenter(i+1),genweight);
					for (int i=0;i<hveto_softpt_denom->GetNbinsX();i++)
						hveto_softpt_denom->Fill(hveto_softpt_denom->GetBinCenter(i+1),genweight);
			//		cout<<hveto_jet3pt_denom->GetBinCenter(1)<<" , "<<hveto_jet3pt_denom->GetBinCenter(2)<<endl;
			
					hNAdJets_bdt->Fill(good_jets, genweight);	
			   	hHTsoftEWK_bdt->Fill(Jet.EWKHTsoft,genweight);
			 	 	hSoft_n2EWK_bdt->Fill(Jet.EWKnsoft2, genweight);
			  		hSoft_n5EWK_bdt->Fill(Jet.EWKnsoft5, genweight);
			  		hSoft_n10EWK_bdt->Fill(Jet.EWKnsoft10, genweight);
					if (good_jets>=3) {
						hJet3_pt_bdt->Fill(jets_pv[2].Pt(),genweight);
						hJet3_eta_bdt->Fill(jets_pv[2].Eta(),genweight);
					}
					if (good_jets==2) hJet3_pt_bdt->Fill(10.,genweight);
					if (good_jets==2) hAdJetHT_bdt->Fill(10.,genweight);
					if (good_jets>=3) hAdJetHT_bdt->Fill(AdJetHT,genweight);	

					hJet1q_eta_bdt->Fill(Qjet1.Eta(),genweight);
					hJet2q_eta_bdt->Fill(Qjet2.Eta(),genweight);
			
					if (good_jets>2) {
						for (int i=0;i<hveto_jet3pt_nom->GetNbinsX();i++)
							if (jets_pv[2].Pt()>hveto_jet3pt_nom->GetBinCenter(i+1)) hveto_jet3pt_nom->Fill(hveto_jet3pt_nom->GetBinCenter(i+1),genweight);
						for (int i=0;i<hveto_ht_nom->GetNbinsX();i++)
							if (AdJetHT>hveto_ht_nom->GetBinCenter(i+1)) hveto_ht_nom->Fill(hveto_ht_nom->GetBinCenter(i+1),genweight);
					}
					for (int i=0;i<hveto_softht_nom->GetNbinsX();i++)
						if (Jet.EWKHTsoft>hveto_softht_nom->GetBinCenter(i+1)) hveto_softht_nom->Fill(hveto_softht_nom->GetBinCenter(i+1),genweight);
					for (int i=0;i<hveto_softpt_nom->GetNbinsX();i++)
						if (Jet.EWKsoft_pt[0] > hveto_softpt_nom->GetBinCenter(i+1)) hveto_softpt_nom->Fill(hveto_softpt_nom->GetBinCenter(i+1),genweight);
				}
				if (bdt>0.84) {
					hNAdJets_bdt2->Fill(good_jets, genweight);	
			   	hHTsoftEWK_bdt2->Fill(Jet.EWKHTsoft,genweight);
			 	 	hSoft_n2EWK_bdt2->Fill(Jet.EWKnsoft2, genweight);
			  		hSoft_n5EWK_bdt2->Fill(Jet.EWKnsoft5, genweight);
			  		hSoft_n10EWK_bdt2->Fill(Jet.EWKnsoft10, genweight);
					if (good_jets>=3) {
						hJet3_pt_bdt2->Fill(jets_pv[2].Pt(),genweight);
						hJet3_eta_bdt2->Fill(jets_pv[2].Eta(),genweight);
					}
					if (good_jets==2) hJet3_pt_bdt2->Fill(10.,genweight);
					if (good_jets==2) hAdJetHT_bdt2->Fill(10.,genweight);
					if (good_jets>=3) hAdJetHT_bdt2->Fill(AdJetHT,genweight);	
					hJet1q_eta_bdt2->Fill(Qjet1.Eta(),genweight);
					hJet2q_eta_bdt2->Fill(Qjet2.Eta(),genweight);
				}
				if (Mqq > 1500) {
					hNAdJets_mjj1->Fill(good_jets, genweight);	
			   	hHTsoftEWK_mjj1->Fill(Jet.EWKHTsoft,genweight);
			 	 	hSoft_n2EWK_mjj1->Fill(Jet.EWKnsoft2, genweight);
			  		hSoft_n5EWK_mjj1->Fill(Jet.EWKnsoft5, genweight);
			  		hSoft_n10EWK_mjj1->Fill(Jet.EWKnsoft10, genweight);
					if (good_jets>=3) hJet3_pt_mjj1->Fill(jets_pv[2].Pt(),genweight);
					if (good_jets==2) hJet3_pt_mjj1->Fill(10.,genweight);
					if (good_jets==2) hAdJetHT_mjj1->Fill(10.,genweight);
					if (good_jets>=3) hAdJetHT_mjj1->Fill(AdJetHT,genweight);	
				}
				if (Mqq > 2500) {
					hNAdJets_mjj2->Fill(good_jets, genweight);	
			   	hHTsoftEWK_mjj2->Fill(Jet.EWKHTsoft,genweight);
			 	 	hSoft_n2EWK_mjj2->Fill(Jet.EWKnsoft2, genweight);
			  		hSoft_n5EWK_mjj2->Fill(Jet.EWKnsoft5, genweight);
			  		hSoft_n10EWK_mjj2->Fill(Jet.EWKnsoft10, genweight);
					if (good_jets>=3) hJet3_pt_mjj2->Fill(jets_pv[2].Pt(),genweight);
					if (good_jets==2) hJet3_pt_mjj2->Fill(10.,genweight);
					if (good_jets==2) hAdJetHT_mjj2->Fill(10.,genweight);
					if (good_jets>=3) hAdJetHT_mjj2->Fill(AdJetHT,genweight);	
				}

		
		if (genweight>0) pos_weight_presel++;
		float mcweight=genweight*events_generated/xsec[file_tag]; 

		if (genweight>0) gen_pos_weight+=mcweight;
		if (genweight<0) gen_neg_weight+=mcweight;
		if (genweight>0) gen_pos+=genweight0;
		if (genweight<0) gen_neg+=genweight0;

				
			global_counter++;
        }

		cout<<counter<<endl;
		TFile file(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_"+heppyVersion+"_"+postfix+".root","recreate");
    
		for (int i=0;i<numArray;++i){
    	    	histArray[i]->SetLineWidth(2);
    	   	histArray[i]->GetYaxis()->SetTitle("N_{events}");
       		histArray[i]->GetYaxis()->SetTitleFont(42);
       		histArray[i]->GetYaxis()->SetTitleSize(0.060);
        		histArray[i]->GetYaxis()->SetTitleOffset(0.8);
        		histArray[i]->SetLineColor(kBlue);
        		histArray[i]->Draw();
        		histArray[i]->Write();
   		}
			hprof_htsoft_pu->SetLineWidth(2); 
			hprof_htsoft_pu->SetLineColor(kBlue); 
			hprof_htsoft_pu->Draw();
			hprof_htsoft_pu->Write();
			hprof_htsoft_pu_rms->Draw();
			hprof_htsoft_pu_rms->Write();
			hprof_htsoft_pu_bdt->SetLineWidth(2); 
			hprof_htsoft_pu_bdt->SetLineColor(kBlue); 
			hprof_htsoft_pu_bdt->Draw();
			hprof_htsoft_pu_bdt->Write();
			hprof_htsoft_pu_rms_bdt->Draw();
			hprof_htsoft_pu_rms_bdt->Write();
    		file.Write();
    		file.Close();
	 ofstream out(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_"+heppyVersion+"_"+postfix+".txt");
	out<< "positive pure selected = "<<gen_pos<<"  , positive weighted selected =  "<<gen_pos_weight<<" , negative pure selected = "<<gen_neg<< ", negative weighted selected = "<<gen_neg_weight<< ", all evetns in the begining = "<<events_generated<<" , xsec = "<<xsec[file_tag]<<endl;
	out<<"positive weight in so many events : "<<  pos_weight_presel<<endl;
	for (int i=0;i<7;i++)
		out<<cut_flow_names[i]<<"\t";
	out<<endl;
	for (int i=0;i<7;i++){
		cut_flow[i]=cut_flow[i];
		out<<cut_flow[i]<<"\t";
	}
	out<<endl;

return 0;
    
}
