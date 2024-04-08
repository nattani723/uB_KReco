#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/MCBase/MCShower.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"


#include "ubana/ParticleID/Algorithms/FiducialVolume.h"
#include "ubana/ParticleID/Algorithms/dQdxSeparatorMarco.h"
#include "ubana/ParticleID/Algorithms/Bragg_Likelihood_Estimator.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "ubobj/UBXSec/SelectionResult.h"
#include "ubobj/UBXSec/TPCObject.h"

#include "ubana/UBXSec/Algorithms/FiducialVolume.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

using namespace std;

// ====================================================================== Local Function Definition to get the reco origin ======================
art::Ptr<simb::MCTruth>TrackIDToMCTruth(art::Event const & evt, std::string _geant_producer, int geant_track_id)
{
	lar_pandora::MCTruthToMCParticles truthToParticles;
	lar_pandora::MCParticlesToMCTruth particlesToTruth;

	lar_pandora::LArPandoraHelper::CollectMCParticles(evt, _geant_producer, truthToParticles, particlesToTruth);

	for (auto iter : particlesToTruth) {
		if (iter.first->TrackId() == geant_track_id) {
			return iter.second;
		}
	}

	art::Ptr<simb::MCTruth> null_ptr;
	return null_ptr;
}	
//================================================================================================================================================

namespace microboone{

class LetsFindKaonsinMCv1 : public art::EDAnalyzer {
public:

    explicit LetsFindKaonsinMCv1(fhicl::ParameterSet const& pset);
    virtual ~LetsFindKaonsinMCv1();

    void beginJob();
    void endSubRun(const art::SubRun& sr);
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    double length(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);
    void reset();
    
private:
    TTree* fEventTree;
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Float_t  reco_numuCC_vtx_xyz[3];
    Float_t  k_mu_min_dis;
    Int_t    is_kmu; //=1 when event has kaon-muon pair
    Int_t    is_mu_sec_k; // =1 when muon is from secondary kaon
    Int_t    k_reco_trkID;
    Float_t  k_reco_tlen;
    Float_t  k_reco_stxyz[3];
    Float_t  k_reco_enxyz[3];
    Float_t  k_reco_theta; // in deg;
    Float_t  k_reco_phi; // in deg;
    Float_t  k_reco_thetaXZ; // in deg
    Float_t  k_reco_thetaYZ; // in deg
    Float_t  k_reco_stcos[3];
    Float_t  k_reco_encos[3];
    Float_t  k_reco_nhits[3];
    Float_t  k_reco_KE[3];
    Float_t  k_reco_rr[3][5000];
    Float_t  k_reco_dqdx[3][5000];
    Float_t  k_reco_dedx[3][5000];
    Float_t  k_reco_x[3][5000];
    Float_t  k_reco_y[3][5000];
    Float_t  k_reco_z[3][5000];
    Float_t  k_reco_chi_p[3];
    Float_t  k_reco_chi_k[3];
    Float_t  k_reco_chi_pi[3];
    Float_t  k_reco_chi_mu[3];
    Float_t  k_reco_pln_pur[3];
    Float_t  k_reco_tot_pur;
    Float_t  k_reco_tot_cmp;
    Float_t  k_reco_pln_cmp[3];
    Int_t    k_tr_pdg;
    Int_t    k_tr_g4ID;
    Int_t    k_tr_origin;
    std::string k_tr_born_process;
    std::string k_tr_end_process;
    Float_t  k_tr_st[3]; 
    Float_t  k_tr_en[3]; 
    Float_t  k_tr_len;  
    Float_t  k_tr_mom; 
    Float_t  k_tr_Pxyz[3]; // momentum components in X,Y and Z directions in MeV/c
    Float_t  k_tr_phi; // Azimuthal angle w.r.t beam direction in degrees
    Float_t  k_tr_theta; // Polar angle w.r.t beam direction in degrees
    Float_t  k_tr_endE; // in MeV
    Float_t  k_tr_stE; // in MeV
    Float_t  k_tr_mass; // in MeV
    Float_t  k_tr_ke;
    Int_t    is_tr_k_cont;
    Int_t    mu_reco_trkID;
    Float_t  mu_reco_tlen;
    Float_t  mu_reco_stxyz[3];
    Float_t  mu_reco_enxyz[3];
    Float_t  mu_reco_theta; // in deg;
    Float_t  mu_reco_phi; // in deg;
    Float_t  mu_reco_thetaXZ; // in deg
    Float_t  mu_reco_thetaYZ; // in deg
    Float_t  mu_reco_stcos[3];
    Float_t  mu_reco_encos[3];
    Float_t  mu_reco_nhits[3];
    Float_t  mu_reco_KE[3];
    Float_t  mu_reco_rr[3][5000];
    Float_t  mu_reco_dqdx[3][5000];
    Float_t  mu_reco_dedx[3][5000];
    Float_t  mu_reco_x[3][5000];
    Float_t  mu_reco_y[3][5000];
    Float_t  mu_reco_z[3][5000];
    Float_t  mu_reco_chi_p[3];
    Float_t  mu_reco_chi_k[3];
    Float_t  mu_reco_chi_pi[3];
    Float_t  mu_reco_chi_mu[3];
    Float_t  mu_reco_pln_pur[3];
    Float_t  mu_reco_tot_pur;
    Float_t  mu_reco_tot_cmp;
    Float_t  mu_reco_pln_cmp[3];
    Int_t    mu_tr_pdg;
    Int_t    mu_tr_g4ID;
    Int_t    mu_tr_origin;
    std::string mu_tr_born_process;
    std::string mu_tr_end_process;
    Float_t  mu_tr_st[3]; 
    Float_t  mu_tr_en[3]; 
    Float_t  mu_tr_len;  
    Float_t  mu_tr_mom; 
    Float_t  mu_tr_Pxyz[3]; 
    Float_t  mu_tr_phi; 
    Float_t  mu_tr_theta; 
    Float_t  mu_tr_endE; 
    Float_t  mu_tr_stE; 
    Float_t  mu_tr_mass; 
    Float_t  mu_tr_ke;
    Int_t    is_tr_mu_cont;
    std::string sec_k_born_process;
    Int_t   is_eve_5cm;
    Int_t   is_eve_10cm;
    Float_t k_mu_reco_op_ang;
    Int_t    n_nu;
    Float_t  nu_eng[2]; // in MeV
    Int_t    nu_pdg[2];
    Int_t    nu_ccnc[2]; // neutrino interaction type: 0=Charged current (CC), 1=Neutral current (NC)
    Int_t    nu_mode[2]; // neutrino nucleus 0=Quasi-elastic or Elastic, 1=Resonant (RES), 2=DIS, 3=Coherent production, 10=MEC
    Float_t  nu_st[2][3];
    Int_t    true_k_ext;
    Int_t    true_k_cont;
    Int_t    n_true_dauts;
    Int_t    n_true_mu;
    Int_t    n_true_pi;
    Int_t    true_mu_cont[2];
    Int_t    true_pi_cont[2];
    
    TTree* fPOT;
    Double_t potbnb;
    
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fShowerModuleLabel;
    std::string fPOTModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fParticleIDModuleLabel;
    std::string fHitTruthAssns;
    std::string fHitTrackAssns;
    std::string fHitShowerAssns;
    std::string fMCShowerModuleLabel;
    std::string m_pfp_producer;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
    bool  fSaveGenieInfo;
    float  fklen;
    int    fkhit;
    float  fk_vtx_dis;
    float  fk_low_k;
    float  fk_up_k;
    float  fk_low_p;
    float  fk_up_p;
    float  fk_mu_dis;
    float  fmu_tlen_up;
    float  fmu_tlen_low;
    int    fmuhit;
    float  fmu_low_p;		    
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
}; 

//========================================================================
LetsFindKaonsinMCv1::LetsFindKaonsinMCv1(fhicl::ParameterSet const& pset) :
EDAnalyzer(pset),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","gaushit")),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel","largeant")),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","generator")),  
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
  fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel","")        ),		
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel","")          ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
  fParticleIDModuleLabel    (pset.get< std::string >("ParticleIDModuleLabel","pandoracalipidSCE")),
  fHitTruthAssns            (pset.get< std::string >("HitTruthAssn","gaushitTruthMatch")), 
  fHitTrackAssns            (pset.get< std::string >("HitTrackAssns","pandora")), 
  fHitShowerAssns           (pset.get< std::string >("HitShowerAssns","pandora")), 
  fMCShowerModuleLabel      (pset.get< std::string >("MCShowerModuleLabel","mcreco")), 
  m_pfp_producer            (pset.get< std::string >("pfp_producer", "pandora")),
  fSaveCaloInfo             (pset.get< bool >("SaveCaloInfo",false)),
  fSaveTrackInfo            (pset.get< bool >("SaveTrackInfo",false)),  
  fSaveGenieInfo            (pset.get< bool >("SaveGenieInfo",false)),
  fklen                     (pset.get< float >("klen",5)),
  fkhit                     (pset.get< int >("khit",5)),
  fk_vtx_dis                (pset.get< float >("k_vtx_dis",5)),
  fk_low_k                  (pset.get< float >("k_low_k",4)),
  fk_up_k                   (pset.get< float >("k_up_k",15)),
  fk_low_p                  (pset.get< float >("k_low_p",10)),
  fk_up_p                   (pset.get< float >("k_up_p",45)),
  fk_mu_dis                 (pset.get< float >("k_mu_dis",5)),
  fmu_tlen_up               (pset.get< float >("mu_tlen_up",70)),
  fmu_tlen_low              (pset.get< float >("mu_tlen_low",30)),
  fmuhit                    (pset.get< int >("muhit",10)),
  fmu_low_p                 (pset.get< float >("mu_low_p",50))
{
  if (fSaveTrackInfo == false) fSaveCaloInfo = false;
}
 
//========================================================================
LetsFindKaonsinMCv1::~LetsFindKaonsinMCv1(){
}
//========================================================================

//========================================================================
void LetsFindKaonsinMCv1::beginJob(){
  std::cout<<"job begin..."<<std::endl;
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
  fEventTree->Branch("event", &event,"event/I");
  fEventTree->Branch("run", &run,"run/I");
  fEventTree->Branch("subrun", &subrun,"surbrun/I");
  fEventTree->Branch("reco_numuCC_vtx_xyz", reco_numuCC_vtx_xyz,"reco_numuCC_vtx_xyz[3]/F");
  fEventTree->Branch("k_mu_min_dis", &k_mu_min_dis,"k_mu_min_dis/F");
  fEventTree->Branch("is_kmu", &is_kmu,"is_kmu/I");
  fEventTree->Branch("is_mu_sec_k", &is_mu_sec_k,"is_mu_sec_k/I");
  fEventTree->Branch("k_reco_trkID", &k_reco_trkID,"k_reco_trkID/I");
  fEventTree->Branch("k_reco_tlen", &k_reco_tlen,"k_reco_tlen/F");
  fEventTree->Branch("k_reco_stxyz", k_reco_stxyz,"k_reco_stxyz[3]/F");
  fEventTree->Branch("k_reco_enxyz", k_reco_enxyz,"k_reco_enxyz[3]/F");
  fEventTree->Branch("k_reco_theta", &k_reco_theta,"k_reco_theta/F");
  fEventTree->Branch("k_reco_phi", &k_reco_phi,"k_reco_phi/F");
  fEventTree->Branch("k_reco_thetaXZ", &k_reco_thetaXZ,"k_reco_thetaXZ/F");
  fEventTree->Branch("k_reco_thetaYZ", &k_reco_thetaYZ,"k_reco_thetaYZ/F");
  fEventTree->Branch("k_reco_stcos", k_reco_stcos,"k_reco_stcos[3]/F");
  fEventTree->Branch("k_reco_encos", k_reco_encos,"k_reco_encos[3]/F");
  fEventTree->Branch("k_reco_nhits", k_reco_nhits,"k_reco_nhits[3]/F");
  fEventTree->Branch("k_reco_KE", k_reco_KE,"k_reco_KE[3]/F");
  fEventTree->Branch("k_reco_rr", k_reco_rr,"k_reco_rr[3][5000]/F");
  fEventTree->Branch("k_reco_dqdx", k_reco_dqdx,"k_reco_dqdx[3][5000]/F");
  fEventTree->Branch("k_reco_dedx", k_reco_dedx,"k_reco_dedx[3][5000]/F");
  fEventTree->Branch("k_reco_x", k_reco_x,"k_reco_x[3][5000]/F");
  fEventTree->Branch("k_reco_y", k_reco_y,"k_reco_y[3][5000]/F");
  fEventTree->Branch("k_reco_z", k_reco_z,"k_reco_z[3][5000]/F");
  fEventTree->Branch("k_reco_chi_p", k_reco_chi_p,"k_reco_chi_p[3]/F");
  fEventTree->Branch("k_reco_chi_k", k_reco_chi_k,"k_reco_chi_k[3]/F");
  fEventTree->Branch("k_reco_chi_pi", k_reco_chi_pi,"k_reco_chi_pi[3]/F");
  fEventTree->Branch("k_reco_chi_mu", k_reco_chi_mu,"k_reco_chi_mu[3]/F");
  fEventTree->Branch("k_reco_pln_pur", k_reco_pln_pur,"k_reco_pln_pur[3][3]/F");
  fEventTree->Branch("k_reco_tot_pur", &k_reco_tot_pur,"k_reco_tot_pur/F");
  fEventTree->Branch("k_reco_pln_cmp", k_reco_pln_cmp,"k_reco_pln_cmp[3][3]/F");
  fEventTree->Branch("k_reco_tot_cmp", &k_reco_tot_cmp,"k_reco_tot_cmp/F");
  fEventTree->Branch("k_tr_pdg", &k_tr_pdg,"k_tr_pdg/I");
  fEventTree->Branch("k_tr_g4ID", &k_tr_g4ID,"k_tr_g4ID/I");
  fEventTree->Branch("k_tr_origin", &k_tr_origin,"k_tr_origin/I");
  fEventTree->Branch("k_tr_born_process", &k_tr_born_process);
  fEventTree->Branch("k_tr_end_process", &k_tr_end_process);
  fEventTree->Branch("k_tr_st", &k_tr_st,"k_tr_st[3]/F");
  fEventTree->Branch("k_tr_en", &k_tr_en,"k_tr_en[3]/F");
  fEventTree->Branch("k_tr_len", &k_tr_len,"k_tr_len/F");
  fEventTree->Branch("k_tr_mom", &k_tr_mom,"k_tr_mom/F");
  fEventTree->Branch("k_tr_Pxyz", k_tr_Pxyz,"k_tr_Pxyz[3]/F");
  fEventTree->Branch("k_tr_phi", &k_tr_phi,"k_tr_phi/F");
  fEventTree->Branch("k_tr_theta", &k_tr_theta,"k_tr_theta/F");
  fEventTree->Branch("k_tr_endE", &k_tr_endE,"k_tr_endE/F");
  fEventTree->Branch("k_tr_stE", &k_tr_stE,"k_tr_stE/F");
  fEventTree->Branch("k_tr_mass", &k_tr_mass,"k_tr_mass/F");
  fEventTree->Branch("k_tr_ke", &k_tr_ke,"k_tr_ke/F");
  fEventTree->Branch("is_tr_k_cont", &is_tr_k_cont,"is_tr_k_cont/I");
  fEventTree->Branch("mu_reco_trkID", &mu_reco_trkID,"mu_reco_trkID/I");
  fEventTree->Branch("mu_reco_tlen", &mu_reco_tlen,"mu_reco_tlen/F");
  fEventTree->Branch("mu_reco_stxyz", mu_reco_stxyz,"mu_reco_stxyz[3]/F");
  fEventTree->Branch("mu_reco_enxyz", mu_reco_enxyz,"mu_reco_enxyz[3]/F");
  fEventTree->Branch("mu_reco_theta", &mu_reco_theta,"mu_reco_theta/F");
  fEventTree->Branch("mu_reco_phi", &mu_reco_phi,"mu_reco_phi/F");
  fEventTree->Branch("mu_reco_thetaXZ", &mu_reco_thetaXZ,"mu_reco_thetaXZ/F");
  fEventTree->Branch("mu_reco_thetaYZ", &mu_reco_thetaYZ,"mu_reco_thetaYZ/F");
  fEventTree->Branch("mu_reco_stcos", mu_reco_stcos,"mu_reco_stcos[3]/F");
  fEventTree->Branch("mu_reco_encos", mu_reco_encos,"mu_reco_encos[3]/F");
  fEventTree->Branch("mu_reco_nhits", mu_reco_nhits,"mu_reco_nhits[3]/F");
  fEventTree->Branch("mu_reco_KE", mu_reco_KE,"mu_reco_KE[3]/F");
  fEventTree->Branch("mu_reco_rr", mu_reco_rr,"mu_reco_rr[3][5000]/F");
  fEventTree->Branch("mu_reco_dqdx", mu_reco_dqdx,"mu_reco_dqdx[3][5000]/F");
  fEventTree->Branch("mu_reco_dedx", mu_reco_dedx,"mu_reco_dedx[3][5000]/F");
  fEventTree->Branch("mu_reco_x", mu_reco_x,"mu_reco_x[3][5000]/F");
  fEventTree->Branch("mu_reco_y", mu_reco_y,"mu_reco_y[3][5000]/F");
  fEventTree->Branch("mu_reco_z", mu_reco_z,"mu_reco_z[3][5000]/F");
  fEventTree->Branch("mu_reco_chi_p", mu_reco_chi_p,"mu_reco_chi_p[3]/F");
  fEventTree->Branch("mu_reco_chi_k", mu_reco_chi_k,"mu_reco_chi_k[3]/F");
  fEventTree->Branch("mu_reco_chi_pi", mu_reco_chi_pi,"mu_reco_chi_pi[3]/F");
  fEventTree->Branch("mu_reco_chi_mu", mu_reco_chi_mu,"mu_reco_chi_mu[3]/F");
  fEventTree->Branch("mu_reco_pln_pur", mu_reco_pln_pur,"mu_reco_pln_pur[3][3]/F");
  fEventTree->Branch("mu_reco_tot_pur", &mu_reco_tot_pur,"mu_reco_tot_pur/F");
  fEventTree->Branch("mu_reco_pln_cmp", mu_reco_pln_cmp,"mu_reco_pln_cmp[3][3]/F");
  fEventTree->Branch("mu_reco_tot_cmp", &mu_reco_tot_cmp,"mu_reco_tot_cmp/F");
  fEventTree->Branch("mu_tr_pdg", &mu_tr_pdg,"mu_tr_pdg/I");
  fEventTree->Branch("mu_tr_g4ID", &mu_tr_g4ID,"mu_tr_g4ID/I");
  fEventTree->Branch("mu_tr_origin", &mu_tr_origin,"mu_tr_origin/I");
  fEventTree->Branch("mu_tr_born_process", &mu_tr_born_process);
  fEventTree->Branch("mu_tr_end_process", &mu_tr_end_process);
  fEventTree->Branch("mu_tr_st", &mu_tr_st,"mu_tr_st[3]/F");
  fEventTree->Branch("mu_tr_en", &mu_tr_en,"mu_tr_en[3]/F");
  fEventTree->Branch("mu_tr_len", &mu_tr_len,"mu_tr_len/F");
  fEventTree->Branch("mu_tr_mom", &mu_tr_mom,"mu_tr_mom/F");
  fEventTree->Branch("mu_tr_Pxyz", mu_tr_Pxyz,"mu_tr_Pxyz[3]/F");
  fEventTree->Branch("mu_tr_phi", &mu_tr_phi,"mu_tr_phi/F");
  fEventTree->Branch("mu_tr_theta", &mu_tr_theta,"mu_tr_theta/F");
  fEventTree->Branch("mu_tr_endE", &mu_tr_endE,"mu_tr_endE/F");
  fEventTree->Branch("mu_tr_stE", &mu_tr_stE,"mu_tr_stE/F");
  fEventTree->Branch("mu_tr_mass", &mu_tr_mass,"mu_tr_mass/F");
  fEventTree->Branch("mu_tr_ke", &mu_tr_ke,"mu_tr_ke/F");
  fEventTree->Branch("is_tr_mu_cont", &is_tr_mu_cont,"is_tr_mu_cont/I");
  fEventTree->Branch("sec_k_born_process", &sec_k_born_process);
  fEventTree->Branch("is_eve_5cm", &is_eve_5cm,"is_eve_5cm/I");
  fEventTree->Branch("is_eve_10cm", &is_eve_10cm,"is_eve_10cm/I");
  fEventTree->Branch("k_mu_reco_op_ang", &k_mu_reco_op_ang,"k_mu_reco_op_ang/F");
  fEventTree->Branch("n_nu", &n_nu,"n_nu/I");
  fEventTree->Branch("nu_eng", nu_eng,"nu_eng[2]/F");
  fEventTree->Branch("nu_pdg", nu_pdg,"nu_pdg[2]/I");
  fEventTree->Branch("nu_ccnc", nu_ccnc,"nu_ccnc[2]/I");
  fEventTree->Branch("nu_mode", nu_mode,"nu_mode[2]/I");
  fEventTree->Branch("nu_st", nu_st,"nu_st[2][3]/F");
  fEventTree->Branch("true_k_ext", &true_k_ext,"true_k_ext/I");
  fEventTree->Branch("true_k_cont", &true_k_cont,"true_k_cont/I");
  fEventTree->Branch("n_true_dauts", &n_true_dauts,"n_true_dauts/I");
  fEventTree->Branch("n_true_mu", &n_true_mu,"n_true_mu/I");
  fEventTree->Branch("n_true_pi", &n_true_pi,"n_true_pi/I");
  fEventTree->Branch("true_mu_cont", true_mu_cont,"true_mu_cont[2]/I");
  fEventTree->Branch("true_pi_cont", true_pi_cont,"true_pi_cont[2]/I");
  
  fPOT = tfs->make<TTree>("POT", "Pot Tree from Reco");
  fPOT->Branch("potbnb",&potbnb,"potbnb/D");
}

//========================================================================
//========================================================================
void LetsFindKaonsinMCv1::beginRun(const art::Run&){
  mf::LogInfo("LetsFindKaonsinMCv1")<<"begin run..."<<std::endl;
}
//========================================================================

//========================================================================

//========================================================================

void LetsFindKaonsinMCv1::analyze( const art::Event& evt){
     reset();  
     art::Handle< std::vector<recob::Track> > trackListHandle;
     std::vector<art::Ptr<recob::Track> > tracklist;
     if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
     
     art::Handle< std::vector<recob::Shower> > showerListHandle;
     std::vector<art::Ptr<recob::Shower> > showerlist;
     if(evt.getByLabel(fShowerModuleLabel,showerListHandle)) art::fill_ptr_vector(showerlist, showerListHandle);
     
     art::Handle< std::vector<recob::Hit> > hitListHandle;
     std::vector<art::Ptr<recob::Hit> > hitlist;
     if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
        art::fill_ptr_vector(hitlist, hitListHandle);
     
     art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
     std::vector<art::Ptr<simb::MCTruth> > mclist;
     if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);

     art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
     std::vector< art::Ptr<simb::MCParticle> > ptList;
     if (evt.getByLabel(fLArG4ModuleLabel, mcParticleHandle))
         art::fill_ptr_vector(ptList, mcParticleHandle); 
     
     art::Handle< std::vector<sim::MCShower> > mcshowerh;
     evt.getByLabel(fMCShowerModuleLabel, mcshowerh);
     
     art::Handle<std::vector<recob::PFParticle>> pfparticles;
     evt.getByLabel(m_pfp_producer, pfparticles);
     
     art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
     art::FindManyP<anab::ParticleID> trackPIDAssn(trackListHandle, evt, fParticleIDModuleLabel);
     art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, evt, fHitTrackAssns);
     art::FindManyP<recob::Hit> hits_from_showers(showerListHandle, evt, fHitShowerAssns);
     art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle,evt,fHitTruthAssns);
     
     if(pfparticles->size()>0) {
	art::FindManyP<recob::Vertex> pfparticleVertexAssn(pfparticles, evt, "pandora");
	if(pfparticleVertexAssn.isValid()){
           lar_pandora::PFParticleVector pfneutrinos(0);
	   for(unsigned int i=0; i<pfparticles->size(); ++i) {
               art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);
               if(pfparticle->IsPrimary() && pfparticle->PdgCode()==14) {
                  pfneutrinos.push_back(pfparticle);
               }
           }
	   
	   if(pfneutrinos.size()==1){
	      art::Ptr<recob::PFParticle> pfnu=pfneutrinos.front();		     		     
              const recob::Vertex::Point_t &neutrino_vtx=pfparticleVertexAssn.at(pfnu.key()).front()->position();
              reco_numuCC_vtx_xyz[0]=neutrino_vtx.X();reco_numuCC_vtx_xyz[1]=neutrino_vtx.Y();reco_numuCC_vtx_xyz[2]=neutrino_vtx.Z();
	   }
	}
     }
     
     run = evt.run();
     subrun = evt.subRun();
     event = evt.id().event();
     
     std::vector<int> kaon_can_trkID_vec;
     std::vector<int> muon_can_trkID_vec;
     
     int k_can_trkid=-9999;
     int mu_can_trkid=-9999;
     
     int is_kaon=0;
     int k_g4_trkID=-9999;
     int is_muon=0;
     int is_sec_k_muon=0;
     
     if(reco_numuCC_vtx_xyz[0]!=-9999 && reco_numuCC_vtx_xyz[1]!=-9999 && reco_numuCC_vtx_xyz[2]!=-9999){
        int NTracks=tracklist.size();
        int trk_mult=0;
        
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	for(int i=0; i<NTracks;++i){
	    art::Ptr<recob::Track> ptrack(trackListHandle, i);
	    const recob::Track& track = *ptrack;
	    TVector3 pos;
	    pos=track.Vertex<TVector3>();
	    float st_vtx=TMath::Sqrt((reco_numuCC_vtx_xyz[0]-pos.X())*(reco_numuCC_vtx_xyz[0]-pos.X()) + (reco_numuCC_vtx_xyz[1]-pos.Y())*(reco_numuCC_vtx_xyz[1]-pos.Y()) + (reco_numuCC_vtx_xyz[2]-pos.Z())*(reco_numuCC_vtx_xyz[2]-pos.Z()));
	    if(st_vtx<fk_vtx_dis || st_vtx==fk_vtx_dis)	trk_mult++;
	    if(trk_mult>2 || trk_mult==2) break;
	} // loop over reco tracks
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	if((trk_mult==2 || trk_mult>2) && (NTracks==3 || NTracks>3)){
		
	    ///////////////////////////////////// Loop over tracklist to find kaon candidates ///////////////////////////////////////////////////////////////	
		
	    for(int i=0; i<NTracks;++i){
	        art::Ptr<recob::Track> ptrack(trackListHandle,i);
	        const recob::Track& track = *ptrack;
	        TVector3 pos,end;
		pos=track.Vertex<TVector3>();
                end=track.End<TVector3>();
		if((pos.X()>0 && pos.X()<256.35) && (end.X()>0 && end.X()<256.35) && (pos.Y()>-116.5 && pos.Y()<116.5) && (end.Y()>-116.5 && end.Y()<116.5) && (pos.Z()>0 && pos.Z()<1036.8) && (end.Z()>0 && end.Z()<1036.8)){
		    float st_vtx=TMath::Sqrt((reco_numuCC_vtx_xyz[0]-pos.X())*(reco_numuCC_vtx_xyz[0]-pos.X()) + (reco_numuCC_vtx_xyz[1]-pos.Y())*(reco_numuCC_vtx_xyz[1]-pos.Y()) + (reco_numuCC_vtx_xyz[2]-pos.Z())*(reco_numuCC_vtx_xyz[2]-pos.Z()));
		    if(st_vtx<fk_vtx_dis || st_vtx==fk_vtx_dis){
		       float tlen=track.Length();
		       if(tlen>fklen || tlen==fklen){
		          int p0_nhit=0;int p1_nhit=0;int p2_nhit=0;
			  
			  std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);
			  for(unsigned int ical=0; ical<calos.size(); ++ical){
		              if(!calos[ical]) continue;
		              if(!calos[ical]->PlaneID().isValid) continue;
		              int planenum=calos[ical]->PlaneID().Plane;
	                      if(planenum<0||planenum>2) continue; 
		              const size_t NHits=calos[ical]->dEdx().size();
		              if(planenum==0) p0_nhit=int(NHits);
			      if(planenum==1) p1_nhit=int(NHits);
			      if(planenum==2) p2_nhit=int(NHits);
		          }
			  
			  if((p0_nhit>fkhit || p0_nhit==fkhit) || (p1_nhit>fkhit || p1_nhit==fkhit) || (p2_nhit>fkhit || p2_nhit==fkhit)){
			      float p0_chi_k=-9999;float p0_chi_p=-9999;
			      float p1_chi_k=-9999;float p1_chi_p=-9999;
			      float p2_chi_k=-9999;float p2_chi_p=-9999;
			      float p0_low_bound=-9999;float p0_up_bound=-9999;
			      float p1_low_bound=-9999;float p1_up_bound=-9999;
			      float p2_low_bound=-9999;float p2_up_bound=-9999;
			      if(trackPIDAssn.isValid()){
		              std::vector<art::Ptr<anab::ParticleID>> pids=trackPIDAssn.at(i);
		              for(size_t ipid=0; ipid<pids.size(); ipid++){
		                  std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[ipid]->ParticleIDAlgScores();
			          for(size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
			              anab::sParticleIDAlgScores AlgScore=AlgScoresVec.at(i_algscore);
			              int planenum=UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
	                              if(planenum<0 || planenum>2) continue;
			              if(AlgScore.fAlgName == "Chi2"){
			                 if(planenum==0){
					    if(TMath::Abs(AlgScore.fAssumedPdg)==321) p0_chi_k=AlgScore.fValue;
					    if(TMath::Abs(AlgScore.fAssumedPdg)==2212) p0_chi_p=AlgScore.fValue;
					 }
					 
					 if(planenum==1){
					    if(TMath::Abs(AlgScore.fAssumedPdg)==321) p1_chi_k=AlgScore.fValue;
					    if(TMath::Abs(AlgScore.fAssumedPdg)==2212) p1_chi_p=AlgScore.fValue;
					 }
					 
					 if(planenum==2){
					    if(TMath::Abs(AlgScore.fAssumedPdg)==321) p2_chi_k=AlgScore.fValue;
					    if(TMath::Abs(AlgScore.fAssumedPdg)==2212) p2_chi_p=AlgScore.fValue;
					 }
				      }
			          }  
		               }
		           }
			   
			   if(p0_chi_k!=-9999 && p0_chi_p!=-9999){
			      p0_low_bound=float(p0_chi_k*p0_chi_k)/(fk_low_k*fk_low_k) + float(p0_chi_p*p0_chi_k)/(fk_low_p*fk_low_p );
			      p0_up_bound=float(p0_chi_k*p0_chi_k)/(fk_up_k*fk_up_k) + float(p0_chi_p*p0_chi_k)/(fk_up_p*fk_up_p );
			   }	  
			   
			   if(p1_chi_k!=-9999 && p1_chi_p!=-9999){
			      p1_low_bound=float(p1_chi_k*p1_chi_k)/(fk_low_k*fk_low_k) + float(p1_chi_p*p1_chi_k)/(fk_low_p*fk_low_p );
			      p1_up_bound=float(p1_chi_k*p1_chi_k)/(fk_up_k*fk_up_k) + float(p1_chi_p*p1_chi_k)/(fk_up_p*fk_up_p );
			   }
			   
			   if(p2_chi_k!=-9999 && p2_chi_p!=-9999){
			      p2_low_bound=float(p2_chi_k*p2_chi_k)/(fk_low_k*fk_low_k) + float(p2_chi_p*p2_chi_k)/(fk_low_p*fk_low_p );
			      p2_up_bound=float(p2_chi_k*p2_chi_k)/(fk_up_k*fk_up_k) + float(p2_chi_p*p2_chi_k)/(fk_up_p*fk_up_p );
			   }
			   
			   int chi_pass=0;
			   
			   if(p0_low_bound!=-9999 && p0_up_bound!=-9999){
			      if(p0_low_bound>1 && p0_up_bound<1) chi_pass=1;
			   }	
			   
			   if(p1_low_bound!=-9999 && p1_up_bound!=-9999){
			      if(p1_low_bound>1 && p1_up_bound<1) chi_pass=1;
			   }
			   
			   if(p2_low_bound!=-9999 && p2_up_bound!=-9999){
			      if(p2_low_bound>1 && p2_up_bound<1) chi_pass=1;
			   }
			    
			   if(chi_pass==1){
			      kaon_can_trkID_vec.push_back(i);
			   } // chi 2 cut passed
			 } // hit cut passed
		       } // pass tlen cut
		    } // start of the track is close to nu vtx
	        } // found a contained track
	    } // loop over track list to find kaon candidates
	    
	    ////////////////////////////////////////////////////////////////// Loop over tracklist to find muon candidates //////////////////////////////////////////////////////////////////////
	    
	    for(int i=0; i<NTracks;++i){
	        art::Ptr<recob::Track> ptrack(trackListHandle,i);
	        const recob::Track& track = *ptrack;
	        TVector3 pos,end;
		pos=track.Vertex<TVector3>();
                end=track.End<TVector3>();
		if((pos.X()>0 && pos.X()<256.35) && (end.X()>0 && end.X()<256.35) && (pos.Y()>-116.5 && pos.Y()<116.5) && (end.Y()>-116.5 && end.Y()<116.5) && (pos.Z()>0 && pos.Z()<1036.8) && (end.Z()>0 && end.Z()<1036.8)){
		    float start_dis=TMath::Sqrt((reco_numuCC_vtx_xyz[0]-pos.X())*(reco_numuCC_vtx_xyz[0]-pos.X()) + (reco_numuCC_vtx_xyz[1]-pos.Y())*(reco_numuCC_vtx_xyz[1]-pos.Y()) + (reco_numuCC_vtx_xyz[2]-pos.Z())*(reco_numuCC_vtx_xyz[2]-pos.Z()));
		    if(start_dis>fklen){
		       float tlen=track.Length();
		       if((tlen>fmu_tlen_low || tlen==fmu_tlen_low) && (tlen<fmu_tlen_up || tlen==fmu_tlen_up)){
		           int p0_nhit=0;int p1_nhit=0;int p2_nhit=0;
			  
			   std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);
			   for(unsigned int ical=0; ical<calos.size(); ++ical){
		               if(!calos[ical]) continue;
		               if(!calos[ical]->PlaneID().isValid) continue;
		               int planenum=calos[ical]->PlaneID().Plane;
	                       if(planenum<0||planenum>2) continue; 
		               const size_t NHits=calos[ical]->dEdx().size();
		               if(planenum==0) p0_nhit=int(NHits);
			       if(planenum==1) p1_nhit=int(NHits);
			       if(planenum==2) p2_nhit=int(NHits);
		          }
			  
			  if((p0_nhit>fmuhit || p0_nhit==fmuhit) || (p1_nhit>fmuhit || p1_nhit==fmuhit) || (p2_nhit>fmuhit || p2_nhit==fmuhit)){
			      float p0_chi_p=-9999;float p1_chi_p=-9999;float p2_chi_p=-9999;
			      if(trackPIDAssn.isValid()){
		              std::vector<art::Ptr<anab::ParticleID>> pids=trackPIDAssn.at(i);
		              for(size_t ipid=0; ipid<pids.size(); ipid++){
		                  std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[ipid]->ParticleIDAlgScores();
			          for(size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
			              anab::sParticleIDAlgScores AlgScore=AlgScoresVec.at(i_algscore);
			              int planenum=UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
	                              if(planenum<0 || planenum>2) continue;
			              if(AlgScore.fAlgName == "Chi2"){
			                 if(planenum==0){
					    if(TMath::Abs(AlgScore.fAssumedPdg)==2212) p0_chi_p=AlgScore.fValue;
				         }
					 
					 if(planenum==1){
					    if(TMath::Abs(AlgScore.fAssumedPdg)==2212) p1_chi_p=AlgScore.fValue;
					 }
					 
					 if(planenum==2){
					    if(TMath::Abs(AlgScore.fAssumedPdg)==2212) p2_chi_p=AlgScore.fValue;
					 }
				      }
			           }  
		                }
		             }
			     
			     if(p0_chi_p>fmu_low_p || p1_chi_p>fmu_low_p || p2_chi_p>fmu_low_p){
			        muon_can_trkID_vec.push_back(i);
			     } // chi p cut passed
			  } // trk hit cut passed
		       } // tlen cut
		    } // start point is away from nu vtx
	        } // track is contained
	    } // loop over track list to find muon candidates
	    
	    /////////////////////////////////////////// Matching track pairs ////////////////////////////////////////////////////////////////////////////////
	    
	    if(kaon_can_trkID_vec.size() && muon_can_trkID_vec.size()){
	       std::vector<int>mu_mindis_can_vec;
	       std::vector<int>k_mindis_can_vec;
	       std::vector<float>close_dis_vec;
	       for(unsigned int i=0; i<kaon_can_trkID_vec.size(); i++){
	           for(unsigned int j=0; j<muon_can_trkID_vec.size(); j++){
		       if(muon_can_trkID_vec[j]!=kaon_can_trkID_vec[i]){
		          
			  art::Ptr<recob::Track> ptrack_k(trackListHandle,kaon_can_trkID_vec[i]);
		          const recob::Track& track_k = *ptrack_k;
	                  TVector3 end_k;
		          end_k=track_k.End<TVector3>();
		      
	                  art::Ptr<recob::Track> ptrack_mu(trackListHandle,muon_can_trkID_vec[j]);
		          const recob::Track& track_mu = *ptrack_mu;
	                  TVector3 pos_mu;
	                  pos_mu=track_mu.Vertex<TVector3>(); 
			  
			  float k_en_mu_st=TMath::Sqrt((end_k.X()-pos_mu.X())*(end_k.X()-pos_mu.X()) + (end_k.Y()-pos_mu.Y())*(end_k.Y()-pos_mu.Y()) + (end_k.Z()-pos_mu.Z())*(end_k.Z()-pos_mu.Z()));
			  
			  if(k_en_mu_st<fk_mu_dis || k_en_mu_st==fk_mu_dis){
			     k_mindis_can_vec.push_back(kaon_can_trkID_vec[i]);
			     mu_mindis_can_vec.push_back(muon_can_trkID_vec[j]);
			     close_dis_vec.push_back(k_en_mu_st);
			  } // start and end of kaon-muon less than given value
		       } // kaon and muon does not have some track ID
		    } // loop over muon ID vector
	        } // looping over kaon iD  vector
		
		if(k_mindis_can_vec.size() && mu_mindis_can_vec.size()){
		   float min2_dis=5e20;
		   int track_index=-1;
		   for(unsigned int k=0; k<close_dis_vec.size(); k++){
		       if(close_dis_vec[k]<min2_dis){
		          min2_dis=close_dis_vec[k];
		          track_index=k;
		       } 
		   }
		   k_can_trkid=k_mindis_can_vec[track_index];
		   mu_can_trkid=mu_mindis_can_vec[track_index];
		   k_mu_min_dis=min2_dis;
		} // have kaon/muon pairs which satisfy minimum-distance criteria
	     } // have candidates for muons and kaons
	 } // at least two tracks coming out of nu vtx and more thatn 3 tracks in the event
     } // numuCC vtx is real
     
     //////////////////////////////////////////////////////////// Filter found Kaon-Muon pair ////////////////////////////////////////////////////////////
     
     if(k_can_trkid!=-9999 && mu_can_trkid!=-9999){
        
	///////////////////////////////////////////////////// Reco-truth matching for kaon /////////////////////////////////////////////  
	
	if(k_can_trkid!=-9999){
	   simb::MCParticle const* matched_mcparticle = NULL;
           std::unordered_map<int,float> trkide;
	   float maxe=-9999, tote=0;
	   std::vector<simb::MCParticle const*> particle_vec;
           std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
           std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(k_can_trkid);
	   
	   for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
               particle_vec.clear(); match_vec.clear();
               particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
	       for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                   trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
		   tote += match_vec[i_p]->energy;
		   if(trkide[ particle_vec[i_p]->TrackId() ]>maxe){
                      maxe=trkide[ particle_vec[i_p]->TrackId() ];
		      matched_mcparticle = particle_vec[i_p];
		   }
	       }
           }
	   
	   if(matched_mcparticle){
	      const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
	      int origin=int(mc_truth->Origin());
	      if(matched_mcparticle->PdgCode()==321 && matched_mcparticle->Process()=="primary" && origin==1){
	         float stx=matched_mcparticle->Vx();float sty=matched_mcparticle->Vy();float stz=matched_mcparticle->Vz();
		 if((stx>0 && stx<256.35) && (sty>-116.5 && sty<116.5) && (stz>0 && stz<1036.8)){
		     is_kaon=1;
	             k_g4_trkID=matched_mcparticle->TrackId();
	         }
	      }
           }
	} // has kaon candiate  
	
	////////////////////////////////////////////// Reco-truth matching for muon //////////////////////////////////////////////////
	
	if(mu_can_trkid!=-9999){
	   simb::MCParticle const* matched_mcparticle = NULL;
           std::unordered_map<int,float> trkide;
	   float maxe=-9999, tote=0;
	   std::vector<simb::MCParticle const*> particle_vec;
           std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
           std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(mu_can_trkid);
	   
	   for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
               particle_vec.clear(); match_vec.clear();
               particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
	       for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                   trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
		   tote += match_vec[i_p]->energy;
		   if(trkide[ particle_vec[i_p]->TrackId() ]>maxe){
                      maxe=trkide[ particle_vec[i_p]->TrackId() ];
		      matched_mcparticle = particle_vec[i_p];
		   }
	       }
           }
	   
	   if(matched_mcparticle){
	      const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
	      int origin=int(mc_truth->Origin());
	      if((matched_mcparticle->PdgCode()==-13 || matched_mcparticle->PdgCode()==211) && matched_mcparticle->Process()=="Decay" && matched_mcparticle->Mother()==k_g4_trkID && origin==1){
	          is_muon=1;
	      }
	      if((matched_mcparticle->PdgCode()==-13 || matched_mcparticle->PdgCode()==211) && matched_mcparticle->Process()=="Decay" && origin==1 && matched_mcparticle->Mother()!=k_g4_trkID){
	          for(auto const& pPart : ptList){
		      const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",pPart->TrackId());
	              int origin=int(mc_truth->Origin());	
		      if(pPart->TrackId()==matched_mcparticle->Mother()){
		         if(pPart->PdgCode()==321 && pPart->Process()!="primary" && pPart->Mother()==k_g4_trkID && origin==1){
			    is_muon=1;
			    is_sec_k_muon=1;
			    break;
			  }
		      }
	          }
	       }
            }
	 } // has muon candidate     
	 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
     } // kaon and muon candidates found
     /////////////////////////////////////////////////////////// Filter found Kaon-Muon pair /////////////////////////////////////////////////////////////
     
     if(is_kaon==1 && is_muon==1){
        is_kmu=1; 
	if(is_sec_k_muon==1) is_mu_sec_k=1;
     }
	
     ///////////////////////////////////////// Accessing primary kaon reconstructed information ////////////////////////////////
	
     if(k_can_trkid!=-9999){
	art::Ptr<recob::Track> ptrack(trackListHandle, k_can_trkid);
	const recob::Track& track = *ptrack;
	TVector3 pos,end;
        pos=track.Vertex<TVector3>();
        end=track.End<TVector3>();
	TVector3 dir_start,dir_end;	
        dir_start=track.VertexDirection<TVector3>();
	dir_end=track.EndDirection<TVector3>();
	   
	k_reco_trkID=track.ID();
	k_reco_tlen=track.Length();
	k_reco_stxyz[0]=pos.X();k_reco_stxyz[1]=pos.Y();k_reco_stxyz[2]=pos.Z();
        k_reco_enxyz[0]=end.X();k_reco_enxyz[1]=end.Y();k_reco_enxyz[2]=end.Z();
	k_reco_theta=dir_start.Theta()*57.3;
	k_reco_phi=dir_start.Phi()*57.3;
	k_reco_thetaXZ=(std::atan2(dir_start.X(), dir_start.Z()))*57.3;
	k_reco_thetaYZ=(std::atan2(dir_start.Y(), dir_start.Z()))*57.3;
	k_reco_stcos[0]=dir_start.X();k_reco_stcos[1]=dir_start.Y();k_reco_stcos[2]=dir_start.Z();
	k_reco_encos[0]=dir_end.X();k_reco_encos[1]=dir_end.Y();k_reco_encos[2]=dir_end.Z();
	   
	std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(k_can_trkid);
        for(unsigned int ical=0; ical<calos.size(); ++ical){
	    if(!calos[ical]) continue;
	    if(!calos[ical]->PlaneID().isValid) continue;
	    int planenum=calos[ical]->PlaneID().Plane;
	    if(planenum<0||planenum>2) continue; 
	    k_reco_nhits[planenum]=int(calos[ical]->dEdx().size());
	    k_reco_KE[planenum]=calos[ical]->KineticEnergy();
	    const size_t NHits=calos[ical]->dEdx().size();
	    for(size_t iHit=0; iHit<NHits; ++iHit){
		k_reco_rr[planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
	        k_reco_dqdx[planenum][iHit]=(calos[ical]->dQdx())[iHit];
	        k_reco_dedx[planenum][iHit]=(calos[ical]->dEdx())[iHit];
		const auto& TrkPos=(calos[ical] -> XYZ())[iHit];
	        k_reco_x[planenum][iHit]=TrkPos.X();
                k_reco_y[planenum][iHit]=TrkPos.Y();
                k_reco_z[planenum][iHit]=TrkPos.Z();
	    }
	 }
	    
	 if(trackPIDAssn.isValid()){
	    std::vector<art::Ptr<anab::ParticleID>> pids=trackPIDAssn.at(k_can_trkid);
	    for(size_t ipid=0; ipid<pids.size(); ipid++){
		std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[ipid]->ParticleIDAlgScores();
	        for(size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
		    anab::sParticleIDAlgScores AlgScore=AlgScoresVec.at(i_algscore);
		    int planenum=UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
	            if(planenum<0 || planenum>2) continue;
		    if(AlgScore.fAlgName == "Chi2"){
		       if(TMath::Abs(AlgScore.fAssumedPdg)==13) k_reco_chi_mu[planenum]=AlgScore.fValue; 
		       else if(TMath::Abs(AlgScore.fAssumedPdg)==211) k_reco_chi_pi[planenum]=AlgScore.fValue;
		       else if(TMath::Abs(AlgScore.fAssumedPdg)==321) k_reco_chi_k[planenum]=AlgScore.fValue;
		       else if(TMath::Abs(AlgScore.fAssumedPdg)==2212) k_reco_chi_p[planenum]=AlgScore.fValue;
		     }
	          }  
	       }
	   }
	     
	   simb::MCParticle const* matched_mcparticle = NULL;
           std::unordered_map<int,float> trkide;
	   std::vector<std::unordered_map<int, float>> trkide_planes(3);
           float maxe=-9999, tote=0;
	   float maxe_p0=0;float maxe_p1=0; float maxe_p2=0;
	   std::vector<simb::MCParticle const*> particle_vec;
           std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
           std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(k_can_trkid);
         
	   for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
               particle_vec.clear(); match_vec.clear();
               particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
	       for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                   trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
		   trkide_planes[hits_from_track[i_h]->WireID().Plane][particle_vec[i_p]->TrackId()] += match_vec[i_p]->energy;
		   tote += match_vec[i_p]->energy;
		   if(trkide[ particle_vec[i_p]->TrackId() ]>maxe){
                      maxe=trkide[ particle_vec[i_p]->TrackId() ];
		      matched_mcparticle = particle_vec[i_p];
		   }
	        }
            }
	      
	    if(matched_mcparticle){
	       for(int p=0; p<3; p++){
	           float plane_tot_e=0;
	           for(auto const & p_id :trkide_planes[p]){
		       plane_tot_e += p_id.second;    
	           }
	           if(p==0) maxe_p0=trkide_planes[p][matched_mcparticle->TrackId()]; 
	           if(p==1) maxe_p1=trkide_planes[p][matched_mcparticle->TrackId()];	
	           if(p==2) maxe_p2=trkide_planes[p][matched_mcparticle->TrackId()];  
	           if(p==0 && plane_tot_e!=0) k_reco_pln_pur[0]=(maxe_p0/plane_tot_e)*100;
	           if(p==1 && plane_tot_e!=0) k_reco_pln_pur[1]=(maxe_p1/plane_tot_e)*100;
	           if(p==2 && plane_tot_e!=0) k_reco_pln_pur[2]=(maxe_p2/plane_tot_e)*100;
	       }
	       if(tote!=0) k_reco_tot_pur=(maxe/tote)*100;     
		   
	       std::vector<simb::MCParticle const *>allhit_particle_vec;
	       std::vector<anab::BackTrackerHitMatchingData const *> allhit_match_vec;
	       std::vector<std::unordered_map<int, float>> allhit_plane_trkide(3);
	       std::unordered_map<int, float> allhit_trkide;	
		 
	       for(size_t i_h=0; i_h<hitlist.size(); i_h++){
	           allhit_particle_vec.clear(); allhit_match_vec.clear();
                   particles_per_hit.get(hitlist[i_h].key(),allhit_particle_vec,allhit_match_vec);
                   for(size_t i_p=0; i_p<allhit_particle_vec.size(); ++i_p){
		       allhit_trkide[allhit_particle_vec[i_p]->TrackId()]+= allhit_match_vec[i_p]->energy;
		       allhit_plane_trkide[hitlist[i_h]->WireID().Plane][allhit_particle_vec[i_p]->TrackId()]+= allhit_match_vec[i_p]->energy;
	           }
	       }
		 
	       float particle_total_e=allhit_trkide[matched_mcparticle->TrackId()]; 	
	       if(particle_total_e!=0)k_reco_tot_cmp=(maxe/particle_total_e)*100;
		 
	       for(int p=0; p<3; p++){
	           float particle_total_pln_e=allhit_plane_trkide[p][matched_mcparticle->TrackId()];
                   if(p==0 && particle_total_pln_e!=0){ 
		      k_reco_pln_cmp[0]=(maxe_p0/particle_total_pln_e)*100;
	           }
	           if(p==1 && particle_total_pln_e!=0){ 
		      k_reco_pln_cmp[1]=(maxe_p1/particle_total_pln_e)*100;
	           }
	           if(p==2 && particle_total_pln_e!=0){ 
		      k_reco_pln_cmp[2]=(maxe_p2/particle_total_pln_e)*100;
	           }
	       }
		 
	       k_tr_pdg=matched_mcparticle->PdgCode();
	       k_tr_g4ID=matched_mcparticle->TrackId();  
	       const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
	       k_tr_origin=int(mc_truth->Origin());
	       k_tr_born_process=matched_mcparticle->Process();
	       k_tr_end_process=matched_mcparticle->EndProcess();
	       k_tr_st[0]=matched_mcparticle->Vx();k_tr_st[1]=matched_mcparticle->Vy();k_tr_st[2]=matched_mcparticle->Vz();
	       k_tr_en[0]=matched_mcparticle->EndX();k_tr_en[1]=matched_mcparticle->EndY();k_tr_en[2]=matched_mcparticle->EndZ();
	       TLorentzVector mcstart, mcend;
	       unsigned int pstarti, pendi;
	       k_tr_len=length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
	       k_tr_mom=matched_mcparticle->Momentum().Vect().Mag()*1000;
	       k_tr_Pxyz[0]=matched_mcparticle->Px()*1000;k_tr_Pxyz[1]=matched_mcparticle->Py()*1000;k_tr_Pxyz[2]=matched_mcparticle->Pz()*1000; 
	       k_tr_phi=matched_mcparticle->Momentum().Phi()*57.3;; 
               k_tr_theta=matched_mcparticle->Momentum().Theta()*57.3;    
	       k_tr_endE=matched_mcparticle->EndE()*1000;
	       k_tr_stE=matched_mcparticle->E()*1000;
	       k_tr_mass=matched_mcparticle->Mass()*1000;
	       k_tr_ke=(matched_mcparticle->E()-matched_mcparticle->Mass())*1000;
	       if((k_tr_st[0]>0 && k_tr_st[0]<256.35) && (k_tr_en[0]>0 && k_tr_en[0]<256.35) && (k_tr_st[1]>-116.5 && k_tr_st[1]<116.5) && (k_tr_en[1]>-116.5 && k_tr_en[1]<116.5) && (k_tr_st[2]>0 && k_tr_st[2]<1036.8) && (k_tr_en[2]>0 && k_tr_en[2]<1036.8)){
		   is_tr_k_cont=1;
	       }
	    }
         } // kaon candidate is there
	
	//////////////////////////////////////// Accessing primaary kaon reconstructed information ////////////////////////////////
	
	/////////////////////////////////////// Accessing decay muon/pion reconstructed information ////////////////////////////// 
	 
	 if(mu_can_trkid!=-9999){
	    art::Ptr<recob::Track> ptrack(trackListHandle, mu_can_trkid);
	    const recob::Track& track = *ptrack;
	    TVector3 pos,end;
            pos=track.Vertex<TVector3>();
            end=track.End<TVector3>();
	    TVector3 dir_start,dir_end;	
            dir_start=track.VertexDirection<TVector3>();
	    dir_end=track.EndDirection<TVector3>();
	   
	    mu_reco_trkID=track.ID();
	    mu_reco_tlen=track.Length();
	    mu_reco_stxyz[0]=pos.X();mu_reco_stxyz[1]=pos.Y();mu_reco_stxyz[2]=pos.Z();
	    mu_reco_enxyz[0]=end.X();mu_reco_enxyz[1]=end.Y();mu_reco_enxyz[2]=end.Z();
	    mu_reco_theta=dir_start.Theta()*57.3;
	    mu_reco_phi=dir_start.Phi()*57.3;
	    mu_reco_thetaXZ=(std::atan2(dir_start.X(), dir_start.Z()))*57.3;
	    mu_reco_thetaYZ=(std::atan2(dir_start.Y(), dir_start.Z()))*57.3;
	    mu_reco_stcos[0]=dir_start.X();mu_reco_stcos[1]=dir_start.Y();mu_reco_stcos[2]=dir_start.Z();
	    mu_reco_encos[0]=dir_end.X();mu_reco_encos[1]=dir_end.Y();mu_reco_encos[2]=dir_end.Z();
	    
	    std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(mu_can_trkid);
            for(unsigned int ical=0; ical<calos.size(); ++ical){
	        if(!calos[ical]) continue;
	        if(!calos[ical]->PlaneID().isValid) continue;
	        int planenum=calos[ical]->PlaneID().Plane;
	        if(planenum<0||planenum>2) continue; 
	        mu_reco_nhits[planenum]=int(calos[ical]->dEdx().size());
	        mu_reco_KE[planenum]=calos[ical]->KineticEnergy();
	        const size_t NHits=calos[ical]->dEdx().size();
	        for(size_t iHit=0; iHit<NHits; ++iHit){
		    mu_reco_rr[planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
	            mu_reco_dqdx[planenum][iHit]=(calos[ical]->dQdx())[iHit];
	            mu_reco_dedx[planenum][iHit]=(calos[ical]->dEdx())[iHit];
		    const auto& TrkPos=(calos[ical] -> XYZ())[iHit];
	            mu_reco_x[planenum][iHit]=TrkPos.X();
                    mu_reco_y[planenum][iHit]=TrkPos.Y();
                    mu_reco_z[planenum][iHit]=TrkPos.Z();
	       }
	    }
	    
	    if(trackPIDAssn.isValid()){
	       std::vector<art::Ptr<anab::ParticleID>> pids=trackPIDAssn.at(mu_can_trkid);
	       for(size_t ipid=0; ipid<pids.size(); ipid++){
		   std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[ipid]->ParticleIDAlgScores();
	           for(size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
		       anab::sParticleIDAlgScores AlgScore=AlgScoresVec.at(i_algscore);
		       int planenum=UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
	               if(planenum<0 || planenum>2) continue;
		       if(AlgScore.fAlgName == "Chi2"){
			  if(TMath::Abs(AlgScore.fAssumedPdg)==13) mu_reco_chi_mu[planenum]=AlgScore.fValue; 
		          else if(TMath::Abs(AlgScore.fAssumedPdg)==211) mu_reco_chi_pi[planenum]=AlgScore.fValue;
			  else if(TMath::Abs(AlgScore.fAssumedPdg)==321) mu_reco_chi_k[planenum]=AlgScore.fValue;
			  else if(TMath::Abs(AlgScore.fAssumedPdg)==2212) mu_reco_chi_p[planenum]=AlgScore.fValue;
		       }
	            }  
		 }
	     }
	     
	     simb::MCParticle const* matched_mcparticle = NULL;
             std::unordered_map<int,float> trkide;
	     std::vector<std::unordered_map<int, float>> trkide_planes(3);
             float maxe=-9999, tote=0;
	     float maxe_p0=0;float maxe_p1=0; float maxe_p2=0;
	     std::vector<simb::MCParticle const*> particle_vec;
             std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
             std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(mu_can_trkid);
         
	     for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
                 particle_vec.clear(); match_vec.clear();
                 particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
	         for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                     trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
		     trkide_planes[hits_from_track[i_h]->WireID().Plane][particle_vec[i_p]->TrackId()] += match_vec[i_p]->energy;
		     tote += match_vec[i_p]->energy;
		     if(trkide[ particle_vec[i_p]->TrackId() ]>maxe){
                        maxe=trkide[ particle_vec[i_p]->TrackId() ];
		        matched_mcparticle = particle_vec[i_p];
		     }
	          }
              }
	      
	      if(matched_mcparticle){
	         for(int p=0; p<3; p++) {
	             float plane_tot_e=0;
	             for(auto const & p_id :trkide_planes[p]){
		         plane_tot_e += p_id.second;    
	             }
	             if(p==0) maxe_p0=trkide_planes[p][matched_mcparticle->TrackId()]; 
	             if(p==1) maxe_p1=trkide_planes[p][matched_mcparticle->TrackId()];	
	             if(p==2) maxe_p2=trkide_planes[p][matched_mcparticle->TrackId()];  
	             if(p==0 && plane_tot_e!=0) mu_reco_pln_pur[0]=(maxe_p0/plane_tot_e)*100;
	             if(p==1 && plane_tot_e!=0) mu_reco_pln_pur[1]=(maxe_p1/plane_tot_e)*100;
	             if(p==2 && plane_tot_e!=0) mu_reco_pln_pur[2]=(maxe_p2/plane_tot_e)*100;
	         }
		 if(tote!=0) mu_reco_tot_pur=(maxe/tote)*100;     
		   
		 std::vector<simb::MCParticle const *>allhit_particle_vec;
	         std::vector<anab::BackTrackerHitMatchingData const *> allhit_match_vec;
	         std::vector<std::unordered_map<int, float>> allhit_plane_trkide(3);
	         std::unordered_map<int, float> allhit_trkide;	
		 
		 for(size_t i_h=0; i_h<hitlist.size(); i_h++){
	             allhit_particle_vec.clear(); allhit_match_vec.clear();
                     particles_per_hit.get(hitlist[i_h].key(),allhit_particle_vec,allhit_match_vec);
                     for(size_t i_p=0; i_p<allhit_particle_vec.size(); ++i_p){
		         allhit_trkide[allhit_particle_vec[i_p]->TrackId()]+= allhit_match_vec[i_p]->energy;
		         allhit_plane_trkide[hitlist[i_h]->WireID().Plane][allhit_particle_vec[i_p]->TrackId()]+= allhit_match_vec[i_p]->energy;
	             }
	         }
		 
		 float particle_total_e=allhit_trkide[matched_mcparticle->TrackId()]; 	
		 if(particle_total_e!=0)mu_reco_tot_cmp=(maxe/particle_total_e)*100;
		 
		 for(int p=0; p<3; p++){
	             float particle_total_pln_e=allhit_plane_trkide[p][matched_mcparticle->TrackId()];
                     if(p==0 && particle_total_pln_e!=0){ 
		        mu_reco_pln_cmp[0]=(maxe_p0/particle_total_pln_e)*100;
	             }
	             if(p==1 && particle_total_pln_e!=0){ 
		        mu_reco_pln_cmp[1]=(maxe_p1/particle_total_pln_e)*100;
	             }
	             if(p==2 && particle_total_pln_e!=0){ 
		        mu_reco_pln_cmp[2]=(maxe_p2/particle_total_pln_e)*100;
	             }
	         }
		 
		 mu_tr_pdg=matched_mcparticle->PdgCode();
		 mu_tr_g4ID=matched_mcparticle->TrackId();  
		 const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
	         mu_tr_origin=int(mc_truth->Origin());
		 mu_tr_born_process=matched_mcparticle->Process();
		 mu_tr_end_process=matched_mcparticle->EndProcess();
		 mu_tr_st[0]=matched_mcparticle->Vx();mu_tr_st[1]=matched_mcparticle->Vy();mu_tr_st[2]=matched_mcparticle->Vz();
		 mu_tr_en[0]=matched_mcparticle->EndX();mu_tr_en[1]=matched_mcparticle->EndY();mu_tr_en[2]=matched_mcparticle->EndZ();
		 TLorentzVector mcstart, mcend;
	         unsigned int pstarti, pendi;
	         mu_tr_len=length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
		 mu_tr_mom=matched_mcparticle->Momentum().Vect().Mag()*1000;
		 mu_tr_Pxyz[0]=matched_mcparticle->Px()*1000;mu_tr_Pxyz[1]=matched_mcparticle->Py()*1000;mu_tr_Pxyz[2]=matched_mcparticle->Pz()*1000; 
		 mu_tr_phi=matched_mcparticle->Momentum().Phi()*57.3;; 
                 mu_tr_theta=matched_mcparticle->Momentum().Theta()*57.3;    
	         mu_tr_endE=matched_mcparticle->EndE()*1000;
	         mu_tr_stE=matched_mcparticle->E()*1000;
	         mu_tr_mass=matched_mcparticle->Mass()*1000;
		 mu_tr_ke=(matched_mcparticle->E()-matched_mcparticle->Mass())*1000;
		 if((mu_tr_st[0]>0 && mu_tr_st[0]<256.35) && (mu_tr_en[0]>0 && mu_tr_en[0]<256.35) && (mu_tr_st[1]>-116.5 && mu_tr_st[1]<116.5) && (mu_tr_en[1]>-116.5 && mu_tr_en[1]<116.5) && (mu_tr_st[2]>0 && mu_tr_st[2]<1036.8) && (mu_tr_en[2]>0 && mu_tr_en[2]<1036.8)){
		     is_tr_mu_cont=1;
	         }
		 
		 if(is_mu_sec_k==1){
		    for(auto const& pPart : ptList){
		        if(pPart->TrackId()==matched_mcparticle->Mother()){
			   sec_k_born_process=pPart->Process();
			   break;
			}
		    } 
	         }
	      }
	   } // muon candidate is there
	
	   /////////////////////////////////////// Accessing decay muon/pion reconstructed information ////////////////////////////  
	   
	   /////////////////////////////////////////// Accessing common reco-variables ////////////////////////////////////////////////
	   
	   if(k_can_trkid!=-9999 && mu_can_trkid!=-9999){
	      float stkx=k_reco_stxyz[0];float stky=k_reco_stxyz[1];float stkz=k_reco_stxyz[2];
	      float enkx=k_reco_enxyz[0];float enky=k_reco_enxyz[1];float enkz=k_reco_enxyz[2];	
	      int is_k_reco_5cm=0; 
	      int is_k_reco_10cm=0;
	      if((stkx>5 && stkx<251.35) && (enkx>5 && enkx<251.35) && (stky>-111.5 && stky<111.5) && (enky>-111.5 && enky<111.5) && (stkz>5 && stkz<1031.8) && (enkz>5 && enkz<1031.8)) is_k_reco_5cm=1;
	      if((stkx>10 && stkx<246.35) && (enkx>10 && enkx<246.35) && (stky>-106.5 && stky<106.5) && (enky>-106.5 && enky<106.5) && (stkz>10 && stkz<1026.8) && (enkz>10 && enkz<1026.8)) is_k_reco_10cm=1;
	      float stmx=mu_reco_stxyz[0];float stmy=mu_reco_stxyz[1];float stmz=mu_reco_stxyz[2];
	      float enmx=mu_reco_enxyz[0];float enmy=mu_reco_enxyz[1];float enmz=mu_reco_enxyz[2];	
	      int is_mu_reco_5cm=0; 
	      int is_mu_reco_10cm=0;
	      if((stmx>5 && stmx<251.35) && (enmx>5 && enmx<251.35) && (stmy>-111.5 && stmy<111.5) && (enmy>-111.5 && enmy<111.5) && (stmz>5 && stmz<1031.8) && (enmz>5 && enmz<1031.8)) is_mu_reco_5cm=1;
	      if((stmx>10 && stmx<246.35) && (enmx>10 && enmx<246.35) && (stmy>-106.5 && stmy<106.5) && (enmy>-106.5 && enmy<106.5) && (stmz>10 && stmz<1026.8) && (enmz>10 && enmz<1026.8)) is_mu_reco_10cm=1;
	      if(is_k_reco_5cm==1 && is_mu_reco_5cm==1) is_eve_5cm=1;
	      if(is_k_reco_10cm==1 && is_mu_reco_10cm==1) is_eve_10cm=1;
	      k_mu_reco_op_ang=(TMath::ACos(k_reco_encos[0]*mu_reco_stcos[0] + k_reco_encos[1]*mu_reco_stcos[1] + k_reco_encos[2]*mu_reco_stcos[2]))*57.3;
	   } // kaon-muon candidates found
	   
	   ///////////////////////////////// Accessing common reco-variables ////////////////////////////////////////////////////////
	   
	   ////////////////////////////// Storing neutrino interaction type information ////////////////////////////////////////////
	   int nu_int=0;
	   for(unsigned int iList = 0; iList<mclist.size(); ++iList){
               nu_int++;
	       nu_eng[iList]=mclist[iList]->GetNeutrino().Nu().E()*1000;
               nu_pdg[iList]=mclist[iList]->GetNeutrino().Nu().PdgCode();
               nu_ccnc[iList]=mclist[iList]->GetNeutrino().CCNC();
               nu_mode[iList]=mclist[iList]->GetNeutrino().Mode();
	       nu_st[iList][0]=mclist[iList]->GetNeutrino().Nu().Vx();nu_st[iList][1]=mclist[iList]->GetNeutrino().Nu().Vy();nu_st[iList][2]=mclist[iList]->GetNeutrino().Nu().Vz();
	       if(nu_int==2) break;
          }
	  n_nu=nu_int;
	  ////////////////////////////// neutrino interaction type information ////////////////////////////////////////////////////
	  
	  /////////////////////////////////////// Storing kaon/muon truth information /////////////////////////////////////////////
	  
	  for(auto const& pPart : ptList){
	      std::string pri("primary");
	      bool isPrimary=0;
	      isPrimary=pPart->Process()==pri;
	      if(isPrimary && pPart->PdgCode()==321){
	         const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",pPart->TrackId());
	         int origin=int(mc_truth->Origin());
		 if(origin==1){
		    float stx=pPart->Vx();float sty=pPart->Vy();float stz=pPart->Vz();
		    if((stx>0 && stx<256.35) && (sty>-116.5 && sty<116.5) && (stz>0 && stz<1036.8)){
		        true_k_ext=1;
			float enx=pPart->EndX();float eny=pPart->EndY();float enz=pPart->EndZ();
			if((enx>0 && enx<256.35) && (eny>-116.5 && eny<116.5) && (enz>0 && enz<1036.8)) true_k_cont=1;
			int g4id=pPart->TrackId();
			for(auto const& sPart : ptList){
			    if(sPart->Process()=="Decay" && sPart->Mother()==g4id){ 
			       n_true_dauts++;
			       float stx=sPart->Vx();float sty=sPart->Vy();float stz=sPart->Vz();
			       float enx=sPart->EndX();float eny=sPart->EndY();float enz=sPart->EndZ();
			       if(sPart->PdgCode()==211){ 
				  n_true_pi++;
				  if((stx>0 && stx<256.35) && (sty>-116.5 && sty<116.5) && (stz>0 && stz<1036.8) && (enx>0 && enx<256.35) && (eny>-116.5 && eny<116.5) && (enz>0 && enz<1036.8) && (n_true_pi<3)){
				      true_pi_cont[n_true_pi-1]=1;
				  }
			       }
			       if(sPart->PdgCode()==-13){ 
				  n_true_mu++;
				  if((stx>0 && stx<256.35) && (sty>-116.5 && sty<116.5) && (stz>0 && stz<1036.8) && (enx>0 && enx<256.35) && (eny>-116.5 && eny<116.5) && (enz>0 && enz<1036.8) && (n_true_mu<3)){
				      true_mu_cont[n_true_mu-1]=1;
				  }
			       }
			    }
			}
			break;
		    }
		 }
	      }
	  }
	  
	  ///////////////////////////////////// Storing kaon/muon truth informaiton //////////////////////////////////////////////
	   
	  fEventTree->Fill();
} // end of analyze function
	   
 /////////////////// Defintion of reset function //////////////////////////////////////////

void LetsFindKaonsinMCv1::reset(){
     run = -9999;
     subrun = -9999;
     event = -9999;
     k_mu_min_dis=-9999;
     is_kmu=-9999;
     is_mu_sec_k=-9999;
     k_reco_trkID=-9999;
     k_reco_tlen=-9999;
     k_reco_theta=-9999;
     k_reco_phi=-9999;
     k_reco_thetaXZ=-9999;
     k_reco_thetaYZ=-9999; 
     k_reco_tot_pur=-9999;
     k_reco_tot_cmp=-9999;
     k_tr_pdg=-9999;
     k_tr_g4ID=-9999;
     k_tr_origin=-9999;
     k_tr_born_process="N/A";
     k_tr_end_process="N/A";
     k_tr_len=-9999;  
     k_tr_mom=-9999; 
     k_tr_phi=-9999; 
     k_tr_theta=-9999; 
     k_tr_endE=-9999; 
     k_tr_stE=-9999; 
     k_tr_mass=-9999; 
     k_tr_ke=-9999;
     is_tr_k_cont=-9999;
     mu_reco_trkID=-9999;
     mu_reco_tlen=-9999;
     mu_reco_theta=-9999; 
     mu_reco_phi=-9999; 
     mu_reco_thetaXZ=-9999; 
     mu_reco_thetaYZ=-9999;
     mu_reco_tot_pur=-9999;
     mu_reco_tot_cmp=-9999;
     mu_tr_pdg=-9999;
     mu_tr_g4ID=-9999;
     mu_tr_origin=-9999;
     mu_tr_born_process="N/A";
     mu_tr_end_process="N/A";
     mu_tr_len=-9999;  
     mu_tr_mom=-99999; 
     mu_tr_phi=-9999; 
     mu_tr_theta=-9999; 
     mu_tr_endE=-9999; 
     mu_tr_stE=-9999; 
     mu_tr_mass=-9999; 
     mu_tr_ke=-9999;
     is_tr_mu_cont=-9999;
     sec_k_born_process="N/A";
     is_eve_5cm=-9999;
     is_eve_10cm=-9999;
     k_mu_reco_op_ang=-9999;
     n_nu=0;
     true_k_ext=-9999;
     true_k_cont=-9999;
     n_true_dauts=0;
     n_true_mu=0;
     n_true_pi=0;
     for(int i=0; i<2; i++){
         nu_eng[i]=-9999;
         nu_pdg[i]=-9999;
         nu_ccnc[i]=-9999; 
         nu_mode[i]=-9999;
	 true_mu_cont[i]=-9999;
	 true_pi_cont[i]=-9999; 
	 for(int j=0; j<3; j++){
             nu_st[i][j]=-9999;
         }
     }
     for(int i=0; i<3; i++){
         reco_numuCC_vtx_xyz[i]=-9999;
	 k_reco_stxyz[i]=-9999;
	 k_reco_enxyz[i]=-9999;
	 k_reco_stcos[i]=-9999;
	 k_reco_encos[i]=-9999;
	 k_reco_nhits[i]=-9999;
	 k_reco_KE[i]=-9999;
	 k_reco_chi_p[i]=-9999;
         k_reco_chi_k[i]=-9999;
         k_reco_chi_pi[i]=-9999;
         k_reco_chi_mu[i]=-9999;
	 k_reco_pln_pur[i]=-9999;
	 k_reco_pln_cmp[i]=-9999;
	 k_tr_st[i]=-9999; 
         k_tr_en[i]=-9999; 
	 k_tr_Pxyz[i]=-9999;
	 mu_reco_stxyz[i]=-9999;
         mu_reco_enxyz[i]=-9999;
	 mu_reco_stcos[i]=-9999;
         mu_reco_encos[i]=-9999;
         mu_reco_nhits[i]=-9999;
	 mu_reco_KE[i]=-9999;
	 mu_reco_chi_p[i]=-9999;
         mu_reco_chi_k[i]=-9999;
         mu_reco_chi_pi[i]=-9999;
         mu_reco_chi_mu[i]=-9999;
	 mu_reco_pln_pur[i]=-9999;
	 mu_reco_pln_cmp[i]=-9999;
	 mu_tr_st[i]=-9999; 
         mu_tr_en[i]=-9999; 
	 mu_tr_Pxyz[i]=-9999; 
	 for(int j=0; j<5000; j++){
	     k_reco_rr[i][j]=-9999;
	     k_reco_dqdx[i][j]=-9999;
	     k_reco_dedx[i][j]=-9999;
	     k_reco_x[i][j]=-9999;
	     k_reco_y[i][j]=-9999;
	     k_reco_z[i][j]=-9999;
	     mu_reco_rr[i][j]=-9999;
             mu_reco_dqdx[i][j]=-9999;
             mu_reco_dedx[i][j]=-9999;
             mu_reco_x[i][j]=-9999;
             mu_reco_y[i][j]=-9999;
             mu_reco_z[i][j]=-9999;
	 }
     }
 }

 //////////////////////// End of definition //////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
 
 double LetsFindKaonsinMCv1::length(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
{
  art::ServiceHandle<geo::Geometry> geom;
  double bnd[6] = {0.,2.*geom->DetHalfWidth(),-geom->DetHalfHeight(),geom->DetHalfHeight(),0.,geom->DetLength()};
  double result = 0.;
  TVector3 disp;
  bool first = true;

  for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) {
    // check if the particle is inside a TPC
    if (p.Vx(i) >= bnd[0] && p.Vx(i) <= bnd[1] && p.Vy(i) >= bnd[2] && p.Vy(i) <= bnd[3] && p.Vz(i) >= bnd[4] && p.Vz(i) <= bnd[5]){
      if(first){
	start = p.Position(i);
	first = false;
	starti = i;
      }else{
	disp -= p.Position(i).Vect();
	result += disp.Mag();
      }
      disp = p.Position(i).Vect();
      end = p.Position(i);
      endi = i;
    }
  }
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////

void LetsFindKaonsinMCv1::endSubRun(const art::SubRun& sr){ 
     art::Handle< sumdata::POTSummary > potListHandle;
     if(sr.getByLabel(fPOTModuleLabel,potListHandle))
        potbnb=potListHandle->totpot;
     fPOT->Fill(); 
} 

//////////////////////////////////////////////////////////////////////////////////////////
	
DEFINE_ART_MODULE(LetsFindKaonsinMCv1)
}


