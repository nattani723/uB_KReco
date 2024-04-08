#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "ubana/ParticleID/Algorithms/FiducialVolume.h"
#include "ubana/ParticleID/Algorithms/dQdxSeparatorMarco.h"
#include "ubana/ParticleID/Algorithms/Bragg_Likelihood_Estimator.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "ubobj/UBXSec/SelectionResult.h"
#include "ubobj/UBXSec/TPCObject.h"

#include "ubana/UBXSec/Algorithms/FiducialVolume.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

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
#include <algorithm>

const int kMaxTracks=20;

using namespace std;

namespace microboone{
	
// ======================== Local Function Definition to get the reco origin ======================
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

//==============================================================================================	
class CCKaonFilter : public art::EDAnalyzer {
  public:

    explicit CCKaonFilter(fhicl::ParameterSet const& pset);
    virtual ~CCKaonFilter();

    void beginJob();
    void analyze(const art::Event& evt);
    void reset();

    double length(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);

    bool isInsideVolume(string volume, double x, double y, double z);

  private:

    TTree* fEventTree;

    Int_t   run;                  
    Int_t   subrun;               
    Int_t   event;

    Float_t true_nu_energy;
    Int_t   true_nu_pdg;
    Int_t   true_nu_mode;
    Int_t   true_nu_ccnc;
    Float_t true_nu_vtx_x;
    Float_t true_nu_vtx_y;
    Float_t true_nu_vtx_z;
    Bool_t  true_nu_vtx_inTPC;
    Bool_t  true_nu_vtx_in5cmTPC;
    Bool_t  true_nu_vtx_inCCInclusiveTPC;
    Bool_t  true_nu_vtx_inOldCCInclusiveTPC;

    Int_t   true_lepton_pdg;
    Float_t true_lepton_p;
    Float_t true_lepton_theta;
    Float_t true_lepton_phi;

    Float_t true_kaon_length;
    Float_t true_kaon_ke;
    Float_t true_kaon_theta;
    Float_t true_kaon_phi;
    Float_t true_kaon_lepton_angle;
    Int_t   true_kaon_end_process;
    Float_t true_kaon_end_ke;
    Float_t true_kaon_end_x;
    Float_t true_kaon_end_y;
    Float_t true_kaon_end_z;
    Bool_t  true_kaon_end_inTPC;
    Bool_t  true_kaon_end_in5cmTPC;
    Bool_t  true_kaon_end_inCCInclusiveTPC;
    Bool_t  true_kaon_end_inOldCCInclusiveTPC;

    Float_t true_kaon_daughter_length;
    Float_t true_kaon_daughter_ke;
    Float_t true_kaon_daughter_angle;
    Int_t   true_kaon_daughter_pdg;
    Float_t true_kaon_daughter_end_x;
    Float_t true_kaon_daughter_end_y;
    Float_t true_kaon_daughter_end_z;
    Bool_t  true_kaon_daughter_end_inTPC;
    Bool_t  true_kaon_daughter_end_in5cmTPC;
    Bool_t  true_kaon_daughter_end_inCCInclusiveTPC;
    Bool_t  true_kaon_daughter_end_inOldCCInclusiveTPC;

    Bool_t reco_nu_cc_filter;

    Float_t reco_nu_vtx_x;
    Float_t reco_nu_vtx_y;
    Float_t reco_nu_vtx_z;
    Bool_t  reco_nu_vtx_inTPC;
    Bool_t  reco_nu_vtx_in5cmTPC;
    Bool_t  reco_nu_vtx_inCCInclusiveTPC;
    Int_t   reco_nu_ndaughters;

    Float_t reco_ccmu_vtx_x;
    Float_t reco_ccmu_vtx_y;
    Float_t reco_ccmu_vtx_z;
    Bool_t  reco_ccmu_vtx_inTPC;
    Bool_t  reco_ccmu_vtx_in5cmTPC;
    Bool_t  reco_ccmu_vtx_inCCInclusiveTPC;

    Int_t   reco_ntracks;
    Float_t reco_track_distance[kMaxTracks];
    Int_t   reco_track_nhits0[kMaxTracks];
    Int_t   reco_track_nhits1[kMaxTracks];
    Int_t   reco_track_nhits2[kMaxTracks];
    Float_t reco_track_length[kMaxTracks];
    Bool_t  reco_track_dir[kMaxTracks];
    Float_t reco_track_chi2k[kMaxTracks];
    Float_t reco_track_chi2p[kMaxTracks];
    Float_t reco_track_chi2pi[kMaxTracks];
    Float_t reco_track_chi2mu[kMaxTracks];
    Bool_t  reco_track_vtx_inTPC[kMaxTracks];
    Bool_t  reco_track_vtx_in5cmTPC[kMaxTracks];
    Bool_t  reco_track_vtx_inCCInclusiveTPC[kMaxTracks];
    Bool_t  reco_track_end_inTPC[kMaxTracks];
    Bool_t  reco_track_end_in5cmTPC[kMaxTracks];
    Bool_t  reco_track_end_inCCInclusiveTPC[kMaxTracks];
    Int_t   reco_track_pdg[kMaxTracks];
    Int_t   reco_track_origin[kMaxTracks];
    Bool_t  reco_track_primary[kMaxTracks];
    Int_t   reco_track_ndaughters2;

    Int_t   reco_track_ndaughters[kMaxTracks];
    Float_t reco_track_daughter_distance[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_vtx_distance[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits0[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits1[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_length[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2k[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2p[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_inTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_in5cmTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_inCCInclusiveTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_inTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_in5cmTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_inCCInclusiveTPC[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_pdg[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_origin[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_primary[kMaxTracks][kMaxTracks];

    Int_t   k_can_trkid; //track id of reco kaon
    Int_t   mu_can_trkid; //track id of reco muon
    Float_t k_mu_can_dis; //distance between kaon and muon end
    Float_t k_mu_open_angle; //angle between kaon and muon end
    Float_t k_vtx_dis; //distance between kaon and muon start

    //Int_t mtc_k_5cm; //kaon track and matched kaon true track end or start is within 5 cm
    //Int_t mtc_k_10cm; //kaon track and matched kaon true track end or start is within 10 cm
    //Float_t mtc_k_endE; //matched true kaon track end energy
    //Int_t mtc_mu_5cm; //muon track and matched muon true track end or start is within 5 cm
    //Int_t mtc_mu_10cm; //muon track and matched muon true track end or start is within 10 cm
    //Int_t mtc_mu_pid; //matched muon true track (should be 13)

    Int_t k_geant_ID; //kaon track matched true particle geant ID
    Int_t k_origin; //kaon track matched true particle origin
    Int_t k_pdg; //kaon track matched true particle pdg
    Int_t k_isPri; //kaon track matched true particle is primary
    Float_t k_endE; //kaon track matched true particle end energy
    Int_t k_ness; //kaon track matched true particle is (pdg==321,primary,endE<510)

    Float_t kaon_vtx_dis; //kaon track and vertex distance (same as k_vtx_dis?)
    Float_t k_plen; //kaon track length
    Float_t k_phi; //kaon track phi
    Float_t k_theta; //kaon track theta
    Int_t k_in_5_TPC; //kaon track is within 5 cm from the TPC edges
    Int_t k_in_CC_TPC; //kaon track is within CC fiducial volume
    Int_t k_hit; //number of kaon dEdx hits from calorimetry
    Float_t k_range; //kaon range from calorimetry
    Float_t k_KE; //kaon kinectic energy from calorimetry
    Float_t k_large_dedx; //dedx median?
    Float_t k_small_dedx; //dedx median?
    Float_t k_chi_p; //kaon track proton chi2 score
    Float_t k_chi_k; //kaon track kaon chi2 score
    Float_t k_chi_pi; //kaon track pion chi2 score
    Float_t k_chi_mu; //kaon track muon chi2 score
    Float_t k_p_max; //kaon track proton bragg peak likelihood maximum forward or backward
    Float_t k_k_max; //kaon track kaon bragg peak likelihood maximum forward or backward
    Float_t k_pi_max; //kaon track pion bragg peak likelihood maximum forward or backward
    Float_t k_mu_max; //kaon track muon bragg peak likelihood maximum forward or backward
    Float_t k_mip_max; //kaon track mip bragg peak likelihood maximum forward or backward
    Float_t k_L1_ratio; //kaon track mip / (proton+kaon)
    Float_t k_LL1; //kaon track log( mip /(proton+kaon) )
    Float_t k_L2_ratio; //kaon track mip+muon+pion / (mip+muon+pion+kaon)
    Float_t k_LL2; //kaon track log( mip+muon+pion / (mip+muon+pion+kaon) )
    Float_t k_Lp_ratio; //kaon track mip+muon / (mip+muon+proton)
    Float_t k_LLp;
    Float_t k_Lk_ratio; //kaon track mip+muon / (mip+muon+kaon)
    Float_t k_LLk; //kaon track log ( (mip+muon) / (mip+muon+kaon) )
    Float_t k_pida_mean; //PIDA
    Float_t k_pida_med; //PIDA
    Float_t k_kde; //PIDA
    Float_t k_trm_dqdx; //truncated mean
    Float_t k_trm_dedx; //truncated mean
    Float_t k_dedx[3][3000]; //kaon track dedx per plane per hit
    Float_t k_rr[3][3000]; //kaon track residual range per plane per hit
    Float_t mu_dedx[3][3000]; //muon track dedx per plane per hit
    Float_t mu_rr[3][3000]; //muon track residual range per plane per hit

    Int_t  mu_pdg; //muon track true matched particle pdg
    Int_t  mu_isDec; //muon track true matched particle is true decay?
    Int_t  mu_origin; //muon track true matched particle event origin (unknown=0,beam=1,cosmic=2,SN=3,single=4)
    Int_t  mu_k_is_Mother; //muon track true matched particle is daugther of true kaon
    Int_t  mu_mom_k_inelas; //muon track true matched particle from inelastic K interaction
    Int_t  mu_ness; //muon track true matched particle is an interesting true muon (mu_origin==1 && mu_pdg==-13 && mu_isDec==1 && mu_k_is_Mother==1)
    Float_t mu_plen; //muon track track length
    Float_t mu_phi; //muon track phi
    Float_t mu_theta; //muon track theta
    Int_t   mu_in_5_TPC; //muon track is within 5 cm from TPC edges
    Int_t   mu_in_CC_TPC; //muon track is within CC fiducial volume
    Float_t mu_KE; //muon track reconstructed kinetic energy
    Int_t   mu_hit; //muon track number of hits
    Float_t mu_range; //muon range from calorimetry
    Float_t mu_large_dedx;
    Float_t mu_small_dedx;
    Float_t mu_chi_p;
    Float_t mu_chi_k;
    Float_t mu_chi_pi;
    Float_t mu_chi_mu;
    Float_t mu_p_max;
    Float_t mu_k_max;
    Float_t mu_pi_max;
    Float_t mu_mu_max;
    Float_t mu_mip_max;
    Float_t mu_L1_ratio;
    Float_t mu_LL1;
    Float_t mu_L2_ratio;
    Float_t mu_LL2;
    Float_t mu_Lp_ratio;
    Float_t mu_LLp;
    Float_t mu_Lk_ratio;
    Float_t mu_LLk;
    Float_t mu_pida_mean;
    Float_t mu_pida_med;
    Float_t mu_kde;
    Float_t mu_trm_dqdx;
    Float_t mu_trm_dedx;
    std::string mu_mom_process; //muon track true matched particle creation process

    Int_t cc_mu_trkid;
    Float_t cc_mu_tlen;
    Float_t cc_mu_phi;
    Float_t cc_mu_theta;
    Float_t cc_mu_range;
    Float_t cc_mu_KE;
    Int_t   cc_mu_hit;
    Float_t cc_mu_large_dedx;
    Float_t cc_mu_small_dedx;
    Float_t cc_dis_vtx;
    Int_t cc_mu_pdg;
    Float_t cc_mu_chi_p;
    Float_t cc_mu_chi_k;
    Float_t cc_mu_chi_pi;
    Float_t cc_mu_chi_mu;
    Float_t cc_mu_p_max;
    Float_t cc_mu_k_max;
    Float_t cc_mu_pi_max;
    Float_t cc_mu_mu_max;
    Float_t cc_mu_mip_max;
    Float_t cc_mu_L1_ratio;
    Float_t cc_mu_LL1;
    Float_t cc_mu_L2_ratio;
    Float_t cc_mu_LL2;
    Float_t cc_mu_Lp_ratio;
    Float_t cc_mu_LLp;
    Float_t cc_mu_Lk_ratio;
    Float_t cc_mu_LLk;
    Float_t cc_mu_pida_mean;
    Float_t cc_mu_pida_med;
    Float_t cc_mu_kde;
    Float_t cc_mu_trm_dqdx;
    Float_t cc_mu_trm_dedx;
    Int_t pri_Mu_is;

    Int_t Tr_pri_mu_pdg;
    Int_t Tr_pri_mu_is;
    Int_t Tr_pri_st_k_is;
    Int_t Tr_K_Inelas;
    Float_t Tr_k_plen;
    Float_t Tr_k_endE;
    Float_t Tr_k_theta;
    Float_t Tr_k_phi;
    Int_t Tr_dec_mu_is;
    Int_t Tr_dec_mu_pi_pdg;
    Float_t Tr_mu_plen;
    Float_t Tr_mu_theta;
    Float_t Tr_mu_phi;
    Int_t Tr_k_inTPC;
    Int_t Tr_mu_inTPC;
    Int_t Tr_k_in_5_TPC;
    Int_t Tr_k_in_CC_TPC;
    Int_t Tr_mu_in_5_TPC;
    Int_t Tr_mu_in_CC_TPC;
    Float_t Tr_kmu_open_ang;

    Int_t longest_trkid;
    Float_t longest_trklen;
    Int_t vtx_5cm_mult;
    Float_t k_start_dedx;
    Float_t k_end_dedx;
    Float_t mu_start_dedx;
    Float_t mu_end_dedx;

    Int_t  kinelas_has_traks;
    Int_t  kinelas_reco_trkID;
    Float_t kinelas_tlen;
    Float_t True_kinelas_KE;
    Float_t True_kinelas_tlen;

    Int_t cut_1;
    Int_t cut_2;
    Int_t cut_3;
    Int_t cut_4;
    Int_t cut_5;
    Int_t cut_6;
    Int_t cut_7;
    Int_t cut_8;
    Int_t cut_9;
    Int_t cut_10;
    Int_t cut_11;
    Int_t cut_12;

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fPIDLabel;
    std::string fHitTruthAssns;
    std::string fHitTrackAssns;
    std::string m_pfp_producer;

    int fplane;

}; // class CCKaonFilter

//========================================================================
CCKaonFilter::CCKaonFilter(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","gaushit")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArG4ModuleLabel","largeant")        ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","generator")     ),  
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","pandora") ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","pandoracaliSCE")  ),
  fPIDLabel                 (pset.get< std::string >("PIDLabel","pandorapid") ),
  fHitTruthAssns            (pset.get< std::string >("HitTruthAssn","gaushitTruthMatch")), 
  fHitTrackAssns            (pset.get< std::string >("HitTrackAssn","pandora")), 
  m_pfp_producer            (pset.get< std::string >("pfp_producer", "pandora")),
  fplane                    (pset.get< int >("plane",2))
{
  fhicl::ParameterSet const p_labels = (pset.get<fhicl::ParameterSet>("ProducerLabels"));
  fhicl::ParameterSet const p_bragg  = (pset.get<fhicl::ParameterSet>("BraggAlgo"));
}
 
//========================================================================
CCKaonFilter::~CCKaonFilter(){
  //destructor
}
//========================================================================

//========================================================================
void CCKaonFilter::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");

  fEventTree->Branch("event", &event, "event/I");
  fEventTree->Branch("run", &run, "run/I");
  fEventTree->Branch("subrun", &subrun, "surbrun/I");

  fEventTree->Branch("true_nu_energy", &true_nu_energy, "true_nu_energy/F");
  fEventTree->Branch("true_nu_pdg", &true_nu_pdg, "true_nu_pdg/I");
  fEventTree->Branch("true_nu_mode", &true_nu_mode, "true_nu_mode/I");
  fEventTree->Branch("true_nu_ccnc", &true_nu_ccnc, "true_nu_ccnc/I");
  fEventTree->Branch("true_nu_vtx_x", &true_nu_vtx_x, "true_nu_vtx_x/F");
  fEventTree->Branch("true_nu_vtx_y", &true_nu_vtx_y, "true_nu_vtx_y/F");
  fEventTree->Branch("true_nu_vtx_z", &true_nu_vtx_z, "true_nu_vtx_z/F");
  fEventTree->Branch("true_nu_vtx_inTPC", &true_nu_vtx_inTPC, "true_nu_vtx_inTPC/O");
  fEventTree->Branch("true_nu_vtx_in5cmTPC", &true_nu_vtx_in5cmTPC, "true_nu_vtx_in5cmTPC/O");
  fEventTree->Branch("true_nu_vtx_inCCInclusiveTPC", &true_nu_vtx_inCCInclusiveTPC, "true_nu_vtx_inCCInclusiveTPC/O");
  fEventTree->Branch("true_nu_vtx_inOldCCInclusiveTPC", &true_nu_vtx_inOldCCInclusiveTPC, "true_nu_vtx_inOlcCCInclusiveTPC/O");

  fEventTree->Branch("true_lepton_pdg", &true_lepton_pdg, "true_lepton_pdg/I");
  fEventTree->Branch("true_lepton_p", &true_lepton_p, "true_lepton_p/F");
  fEventTree->Branch("true_lepton_theta", &true_lepton_theta, "true_lepton_theta/F");
  fEventTree->Branch("true_lepton_phi", &true_lepton_phi, "true_lepton_phi/F");

  fEventTree->Branch("true_kaon_length", &true_kaon_length, "true_kaon_length/F");
  fEventTree->Branch("true_kaon_ke", &true_kaon_ke, "true_kaon_ke/F");
  fEventTree->Branch("true_kaon_theta", &true_kaon_theta, "true_kaon_theta/F");
  fEventTree->Branch("true_kaon_phi", &true_kaon_phi, "true_kaon_phi/F");
  fEventTree->Branch("true_kaon_lepton_angle", &true_kaon_lepton_angle, "true_kaon_lepton_angle/F");
  fEventTree->Branch("true_kaon_end_process", &true_kaon_end_process, "true_kaon_end_process/I");
  fEventTree->Branch("true_kaon_end_ke", &true_kaon_end_ke, "true_kaon_end_ke/F");
  fEventTree->Branch("true_kaon_end_x", &true_kaon_end_x, "true_kaon_end_x/F");
  fEventTree->Branch("true_kaon_end_y", &true_kaon_end_y, "true_kaon_end_y/F");
  fEventTree->Branch("true_kaon_end_z", &true_kaon_end_z, "true_kaon_end_z/F");
  fEventTree->Branch("true_kaon_end_inTPC", &true_kaon_end_inTPC, "true_kaon_end_inTPC/O");
  fEventTree->Branch("true_kaon_end_in5cmTPC", &true_kaon_end_in5cmTPC, "true_kaon_end_in5cmTPC/O");
  fEventTree->Branch("true_kaon_end_inCCInclusiveTPC", &true_kaon_end_inCCInclusiveTPC, "true_kaon_end_inCCInclusiveTPC/O");
  fEventTree->Branch("true_kaon_end_inOldCCInclusiveTPC", &true_kaon_end_inOldCCInclusiveTPC, "true_kaon_end_inOldCCInclusiveTPC/O");

  fEventTree->Branch("true_kaon_daughter_length", &true_kaon_daughter_length, "true_kaon_daughter_length/F");
  fEventTree->Branch("true_kaon_daughter_ke", &true_kaon_daughter_ke, "true_kaon_daughter_ke/F");
  fEventTree->Branch("true_kaon_daughter_angle", &true_kaon_daughter_angle, "true_kaon_daughter_angle/F");
  fEventTree->Branch("true_kaon_daughter_pdg", &true_kaon_daughter_pdg, "true_kaon_daughter_pdg/I");
  fEventTree->Branch("true_kaon_daughter_end_x", &true_kaon_daughter_end_x, "true_kaon_daughter_end_x/F");
  fEventTree->Branch("true_kaon_daughter_end_y", &true_kaon_daughter_end_y, "true_kaon_daughter_end_y/F");
  fEventTree->Branch("true_kaon_daughter_end_z", &true_kaon_daughter_end_z, "true_kaon_daughter_end_z/F");
  fEventTree->Branch("true_kaon_daughter_end_inTPC", &true_kaon_daughter_end_inTPC, "true_kaon_daughter_endInTPC/O");
  fEventTree->Branch("true_kaon_daughter_end_in5cmTPC", &true_kaon_daughter_end_in5cmTPC, "true_kaon_daughter_end_in5cmTPC/O");
  fEventTree->Branch("true_kaon_daughter_end_inCCInclusiveTPC", &true_kaon_daughter_end_inCCInclusiveTPC, "true_kaon_daughter_end_inCCInclusiveTPC/O");
  fEventTree->Branch("true_kaon_daughter_end_inOldCCInclusiveTPC", &true_kaon_daughter_end_inOldCCInclusiveTPC, "true_kaon_daughter_end_inOldCCInclusiveTPC/O");

  fEventTree->Branch("reco_nu_cc_filter", &reco_nu_cc_filter, "reco_nu_cc_filter/O");

  fEventTree->Branch("reco_nu_vtx_x", &reco_nu_vtx_x, "reco_nu_vtx_x/F");
  fEventTree->Branch("reco_nu_vtx_y", &reco_nu_vtx_y, "reco_nu_vtx_y/F");
  fEventTree->Branch("reco_nu_vtx_z", &reco_nu_vtx_z, "reco_nu_vtx_z/F");
  fEventTree->Branch("reco_nu_vtx_inTPC", &reco_nu_vtx_inTPC, "reco_nu_vtx_inTPC/O");
  fEventTree->Branch("reco_nu_vtx_in5cmTPC", &reco_nu_vtx_in5cmTPC, "reco_nu_vtx_in5cmTPC/O");
  fEventTree->Branch("reco_nu_vtx_inCCInclusiveTPC", &reco_nu_vtx_inCCInclusiveTPC, "reco_nu_vtx_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_nu_ndaughters", &reco_nu_ndaughters, "reco_nu_ndaughters/I"          );

  fEventTree->Branch("reco_ccmu_vtx_x", &reco_ccmu_vtx_x, "reco_ccmu_vtx_x/F");
  fEventTree->Branch("reco_ccmu_vtx_y", &reco_ccmu_vtx_y, "reco_ccmu_vtx_y/F");
  fEventTree->Branch("reco_ccmu_vtx_z", &reco_ccmu_vtx_z, "reco_ccmu_vtx_z/F");
  fEventTree->Branch("reco_ccmu_vtx_inTPC",&reco_ccmu_vtx_inTPC, "reco_ccmu_vtx_inTPC/O");
  fEventTree->Branch("reco_ccmu_vtx_in5cmTPC", &reco_ccmu_vtx_in5cmTPC, "reco_ccmu_vtx_in5cmTPC/O");
  fEventTree->Branch("reco_ccmu_vtx_inCCInclusiveTPC", &reco_ccmu_vtx_inCCInclusiveTPC, "reco_ccmu_vtx_inCCInclusiveTPC/O");

  fEventTree->Branch("reco_ntracks", &reco_ntracks, "reco_ntracks/I");
  fEventTree->Branch("reco_track_distance", &reco_track_distance, "reco_track_distance[20]/F");
  fEventTree->Branch("reco_track_nhits0", &reco_track_nhits0, "reco_track_nhits0[20]/I");
  fEventTree->Branch("reco_track_nhits1", &reco_track_nhits1, "reco_track_nhits1[20]/I");
  fEventTree->Branch("reco_track_nhits2", &reco_track_nhits2, "reco_track_nhits2[20]/I");
  fEventTree->Branch("reco_track_length", &reco_track_length, "reco_track_length[20]/F");
  fEventTree->Branch("reco_track_dir", &reco_track_dir, "reco_track_dir[20]/O");
  fEventTree->Branch("reco_track_chi2k", &reco_track_chi2k, "reco_track_chi2k[20]/F");
  fEventTree->Branch("reco_track_chi2p", &reco_track_chi2p, "reco_track_chi2p[20]/F");
  fEventTree->Branch("reco_track_chi2pi", &reco_track_chi2pi, "reco_track_chi2pi[20]/F");
  fEventTree->Branch("reco_track_chi2mu", &reco_track_chi2mu, "reco_track_chi2mu[20]/F");
  fEventTree->Branch("reco_track_vtx_inTPC", &reco_track_vtx_inTPC, "reco_track_vtx_inTPC[20]/O");
  fEventTree->Branch("reco_track_vtx_in5cmTPC", &reco_track_vtx_in5cmTPC, "reco_track_vtx_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_vtx_inCCInclusiveTPC", &reco_track_vtx_inCCInclusiveTPC, "reco_track_vtx_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_end_inTPC", &reco_track_end_inTPC, "reco_track_end_inTPC[20]/O");
  fEventTree->Branch("reco_track_end_in5cmTPC", &reco_track_end_in5cmTPC, "reco_track_end_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_end_inCCInclusiveTPC", &reco_track_end_inCCInclusiveTPC, "reco_track_end_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_pdg", &reco_track_pdg, "reco_track_pdg[20]/I");
  fEventTree->Branch("reco_track_origin", &reco_track_origin, "reco_track_origin[20]/I");
  fEventTree->Branch("reco_track_primary", &reco_track_primary, "reco_track_primary[20]/O");
  fEventTree->Branch("reco_track_ndaughters2", &reco_track_ndaughters2, "reco_track_ndaughters2/I");

  fEventTree->Branch("reco_track_ndaughters", &reco_track_ndaughters, "reco_track_ndaughters[20]/I");
  fEventTree->Branch("reco_track_daughter_distance", &reco_track_daughter_distance, "reco_track_daughter_distance[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_distance", &reco_track_daughter_vtx_distance, "reco_track_daughter_vtx_distance[20][20]/F");
  fEventTree->Branch("reco_track_daughter_nhits0", &reco_track_daughter_nhits0, "reco_track_daughter_nhits0[20][20]/I");
  fEventTree->Branch("reco_track_daughter_nhits1", &reco_track_daughter_nhits1, "reco_track_daughter_nhits1[20][20]/I");
  fEventTree->Branch("reco_track_daughter_nhits2", &reco_track_daughter_nhits2, "reco_track_daughter_nhits2[20][20]/I");
  fEventTree->Branch("reco_track_daughter_length", &reco_track_daughter_length, "reco_track_daughter_length[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2k", &reco_track_daughter_chi2k, "reco_track_daughter_chi2k[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2p", &reco_track_daughter_chi2p, "reco_track_daughter_chi2p[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi", &reco_track_daughter_chi2pi, "reco_track_daughter_chi2pi[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu", &reco_track_daughter_chi2mu, "reco_track_daughter_chi2mu[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_inTPC", &reco_track_daughter_vtx_inTPC, "reco_track_daughter_vtx_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_in5cmTPC", &reco_track_daughter_vtx_in5cmTPC, "reco_track_daughter_vtx_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_inCCInclusiveTPC", &reco_track_daughter_vtx_inCCInclusiveTPC, "reco_track_daughter_vtx_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inTPC", &reco_track_daughter_end_inTPC, "reco_track_daughter_end_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_in5cmTPC", &reco_track_daughter_end_in5cmTPC, "reco_track_daughter_end_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inCCInclusiveTPC", &reco_track_daughter_end_inCCInclusiveTPC, "reco_track_daughter_end_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_pdg", &reco_track_daughter_pdg, "reco_track_daughter_pdg[20][20]/I");
  fEventTree->Branch("reco_track_daughter_origin", &reco_track_daughter_origin, "reco_track_daughter_origin[20][20]/I");
  fEventTree->Branch("reco_track_daughter_primary", &reco_track_daughter_primary, "reco_track_daughter_primary[20][20]/O");

  fEventTree->Branch("k_can_trkid", &k_can_trkid,"k_can_trkid/I");
  fEventTree->Branch("mu_can_trkid", &mu_can_trkid,"mu_can_trkid/I");
  fEventTree->Branch("k_mu_can_dis", &k_mu_can_dis,"k_mu_can_dis/F");
  fEventTree->Branch("k_mu_open_angle", &k_mu_open_angle,"k_mu_open_angle/F");
  fEventTree->Branch("k_vtx_dis", &k_vtx_dis,"k_vtx_dis/F");
  //fEventTree->Branch("mtc_k_5cm", &mtc_k_5cm,"mtc_k_5cm/I");
  //fEventTree->Branch("mtc_k_10cm", &mtc_k_10cm,"mtc_k_10cm/I");
  //fEventTree->Branch("mtc_k_endE", &mtc_k_endE,"mtc_k_endE/F");
  //fEventTree->Branch("mtc_mu_5cm", &mtc_mu_5cm,"mtc_mu_5cm/I");
  //fEventTree->Branch("mtc_mu_10cm", &mtc_mu_10cm,"mtc_mu_10cm/I");
  //fEventTree->Branch("mtc_mu_pid", &mtc_mu_pid,"mtc_mu_pid/I");
  fEventTree->Branch("k_geant_ID", &k_geant_ID,"k_geant_ID/I");
  fEventTree->Branch("k_origin", &k_origin,"k_origin/I");
  fEventTree->Branch("k_pdg", &k_pdg,"k_pdg/I");
  fEventTree->Branch("k_isPri", &k_isPri,"k_isPri/I");
  fEventTree->Branch("k_endE", &k_endE,"k_endE/F");
  fEventTree->Branch("k_ness", &k_ness,"k_ness/I");
  fEventTree->Branch("kaon_vtx_dis", &kaon_vtx_dis,"kaon_vtx_dis/F");
  fEventTree->Branch("k_plen", &k_plen,"k_plen/F");
  fEventTree->Branch("k_phi", &k_phi,"k_phi/F");
  fEventTree->Branch("k_theta", &k_theta,"k_theta/F");
  fEventTree->Branch("k_in_5_TPC", &k_in_5_TPC,"k_in_5_TPC/I");
  fEventTree->Branch("k_in_CC_TPC", &k_in_CC_TPC,"k_in_CC_TPC/I");
  fEventTree->Branch("k_hit", &k_hit,"k_hit/I");
  fEventTree->Branch("k_range", &k_range,"k_range/F");
  fEventTree->Branch("k_KE", &k_KE,"k_KE/F");
  fEventTree->Branch("k_large_dedx", &k_large_dedx,"k_large_dedx/F");
  fEventTree->Branch("k_small_dedx", &k_small_dedx,"k_small_dedx/F");
  fEventTree->Branch("k_chi_p", &k_chi_p,"k_chi_p/F");
  fEventTree->Branch("k_chi_k", &k_chi_k,"k_chi_k/F");
  fEventTree->Branch("k_chi_pi", &k_chi_pi,"k_chi_pi/F");
  fEventTree->Branch("k_chi_mu", &k_chi_mu,"k_chi_mu/F");
  fEventTree->Branch("k_p_max", &k_p_max,"k_p_max/F");
  fEventTree->Branch("k_k_max", &k_k_max,"k_k_max/F");
  fEventTree->Branch("k_pi_max", &k_pi_max,"k_pi_max/F");
  fEventTree->Branch("k_mu_max", &k_mu_max,"k_mu_max/F");
  fEventTree->Branch("k_mip_max", &k_mip_max,"k_mip_max/F");
  fEventTree->Branch("k_L1_ratio", &k_L1_ratio,"k_L1_ratio/F");
  fEventTree->Branch("k_LL1", &k_LL1,"k_LL1/F");
  fEventTree->Branch("k_L2_ratio", &k_L2_ratio,"k_L2_ratio/F");
  fEventTree->Branch("k_LL2", &k_LL2,"k_LL2/F");
  fEventTree->Branch("k_Lp_ratio", &k_Lp_ratio,"k_Lp_ratio/F");
  fEventTree->Branch("k_LLp", &k_LLp,"k_LLp/F");
  fEventTree->Branch("k_Lk_ratio", &k_Lk_ratio,"k_Lk_ratio/F");
  fEventTree->Branch("k_LLk", &k_LLk,"k_LLk/F");
  fEventTree->Branch("k_pida_mean", &k_pida_mean,"k_pida_mean/F");
  fEventTree->Branch("k_pida_med", &k_pida_med,"k_pida_med/F");
  fEventTree->Branch("k_kde", &k_kde,"k_kde/F");
  fEventTree->Branch("k_trm_dqdx", &k_trm_dqdx,"k_trm_dqdx/F");
  fEventTree->Branch("k_trm_dedx", &k_trm_dedx,"k_trm_dedx/F");
  fEventTree->Branch("k_dedx",k_dedx,"k_dedx[3][3000]/F");
  fEventTree->Branch("k_rr",k_rr,"k_rr[3][3000]/F");
  fEventTree->Branch("mu_pdg", &mu_pdg,"mu_pdg/I");
  fEventTree->Branch("mu_isDec", &mu_isDec,"mu_isDec/I");
  fEventTree->Branch("mu_origin", &mu_origin,"mu_origin/I");
  fEventTree->Branch("mu_k_is_Mother", &mu_k_is_Mother,"mu_k_is_Mother/I");
  fEventTree->Branch("mu_mom_k_inelas", &mu_mom_k_inelas,"mu_mom_k_inelas/I");
  fEventTree->Branch("mu_ness", &mu_ness,"mu_ness/I");
  fEventTree->Branch("mu_plen", &mu_plen,"mu_plen/F");
  fEventTree->Branch("mu_phi", &mu_phi,"mu_phi/F");
  fEventTree->Branch("mu_theta", &mu_theta,"mu_theta/F");
  fEventTree->Branch("mu_in_5_TPC", &mu_in_5_TPC,"mu_in_5_TPC/I");
  fEventTree->Branch("mu_in_CC_TPC", &mu_in_CC_TPC,"mu_in_CC_TPC/I");
  fEventTree->Branch("mu_KE", &mu_KE,"mu_KE/F");
  fEventTree->Branch("mu_hit", &mu_hit,"mu_hit/I");
  fEventTree->Branch("mu_range", &mu_range,"mu_range/F");
  fEventTree->Branch("mu_large_dedx", &mu_large_dedx,"mu_large_dedx/F");
  fEventTree->Branch("mu_small_dedx", &mu_small_dedx,"mu_small_dedx/F");
  fEventTree->Branch("mu_chi_p", &mu_chi_p,"mu_chi_p/F");
  fEventTree->Branch("mu_chi_k", &mu_chi_k,"mu_chi_k/F");
  fEventTree->Branch("mu_chi_pi", &mu_chi_pi,"mu_chi_pi/F");
  fEventTree->Branch("mu_chi_mu", &mu_chi_mu,"mu_chi_mu/F");
  fEventTree->Branch("mu_p_max", &mu_p_max,"mu_p_max/F");
  fEventTree->Branch("mu_k_max", &mu_k_max,"mu_k_max/F");
  fEventTree->Branch("mu_pi_max", &mu_pi_max,"mu_pi_max/F");
  fEventTree->Branch("mu_mu_max", &mu_mu_max,"mu_mu_max/F");
  fEventTree->Branch("mu_mip_max", &mu_mip_max,"mu_mip_max/F");
  fEventTree->Branch("mu_L1_ratio", &mu_L1_ratio,"mu_L1_ratio/F");
  fEventTree->Branch("mu_LL1", &mu_LL1,"mu_LL1/F");
  fEventTree->Branch("mu_L2_ratio", &mu_L2_ratio,"mu_L2_ratio/F");
  fEventTree->Branch("mu_LL2", &mu_LL2,"mu_LL2/F");
  fEventTree->Branch("mu_Lp_ratio", &mu_Lp_ratio,"mu_Lp_ratio/F");
  fEventTree->Branch("mu_LLp", &mu_LLp,"mu_LLp/F");
  fEventTree->Branch("mu_Lk_ratio", &mu_Lk_ratio,"mu_Lk_ratio/F");
  fEventTree->Branch("mu_LLk", &mu_LLk,"mu_LLk/F");
  fEventTree->Branch("mu_pida_mean", &mu_pida_mean,"mu_pida_mean/F");
  fEventTree->Branch("mu_pida_med", &mu_pida_med,"mu_pida_med/F");
  fEventTree->Branch("mu_kde", &mu_kde,"mu_kde/F");
  fEventTree->Branch("mu_trm_dqdx", &mu_trm_dqdx,"mu_trm_dqdx/F");
  fEventTree->Branch("mu_trm_dedx", &mu_trm_dedx,"mu_trm_dedx/F");
  fEventTree->Branch("mu_dedx",mu_dedx,"mu_dedx[3][3000]/F");
  fEventTree->Branch("mu_rr",mu_rr,"mu_rr[3][3000]/F");
  fEventTree->Branch("mu_mom_process",&mu_mom_process);
  fEventTree->Branch("cc_mu_trkid", &cc_mu_trkid,"cc_mu_trkid/I");
  fEventTree->Branch("cc_mu_tlen", &cc_mu_tlen,"cc_mu_tlen/F");
  fEventTree->Branch("cc_mu_phi", &cc_mu_phi,"cc_mu_phi/F");
  fEventTree->Branch("cc_mu_theta", &cc_mu_theta,"cc_mu_theta/F");
  fEventTree->Branch("cc_mu_range", &cc_mu_range,"cc_mu_range/F");
  fEventTree->Branch("cc_mu_KE", &cc_mu_KE,"cc_mu_KE/F");
  fEventTree->Branch("cc_mu_hit", &cc_mu_hit,"cc_mu_hit/I");
  fEventTree->Branch("cc_mu_large_dedx", &cc_mu_large_dedx,"cc_mu_large_dedx/F");
  fEventTree->Branch("cc_mu_small_dedx", &cc_mu_small_dedx,"cc_mu_small_dedx/F");
  fEventTree->Branch("cc_dis_vtx", &cc_dis_vtx,"cc_dis_vtx/F");
  fEventTree->Branch("cc_mu_pdg", &cc_mu_pdg,"cc_mu_pdg/I");
  fEventTree->Branch("cc_mu_chi_p", &cc_mu_chi_p,"cc_mu_chi_p/F");
  fEventTree->Branch("cc_mu_chi_k", &cc_mu_chi_k,"cc_mu_chi_k/F");
  fEventTree->Branch("cc_mu_chi_pi", &cc_mu_chi_pi,"cc_mu_chi_pi/F");
  fEventTree->Branch("cc_mu_chi_mu", &cc_mu_chi_mu,"cc_mu_chi_mu/F");
  fEventTree->Branch("cc_mu_p_max", &cc_mu_p_max,"cc_mu_p_max/F");
  fEventTree->Branch("cc_mu_k_max", &cc_mu_k_max,"cc_mu_k_max/F");
  fEventTree->Branch("cc_mu_pi_max", &cc_mu_pi_max,"cc_mu_pi_max/F");
  fEventTree->Branch("cc_mu_mu_max", &cc_mu_mu_max,"cc_mu_mu_max/F");
  fEventTree->Branch("cc_mu_mip_max", &cc_mu_mip_max,"cc_mu_mip_max/F");
  fEventTree->Branch("cc_mu_L1_ratio", &cc_mu_L1_ratio,"cc_mu_L1_ratio/F");
  fEventTree->Branch("cc_mu_LL1", &cc_mu_LL1,"cc_mu_LL1/F");
  fEventTree->Branch("cc_mu_L2_ratio", &cc_mu_L2_ratio,"cc_mu_L2_ratio/F");
  fEventTree->Branch("cc_mu_LL2", &cc_mu_LL2,"cc_mu_LL2/F");
  fEventTree->Branch("cc_mu_Lp_ratio", &cc_mu_Lp_ratio,"cc_mu_Lp_ratio/F");
  fEventTree->Branch("cc_mu_LLp", &cc_mu_LLp,"cc_mu_LLp/F");
  fEventTree->Branch("cc_mu_Lk_ratio", &cc_mu_Lk_ratio,"cc_mu_Lk_ratio/F");
  fEventTree->Branch("cc_mu_LLk", &cc_mu_LLk,"cc_mu_LLk/F");
  fEventTree->Branch("cc_mu_pida_mean", &cc_mu_pida_mean,"cc_mu_pida_mean/F");
  fEventTree->Branch("cc_mu_pida_med", &cc_mu_pida_med,"cc_mu_pida_med/F");
  fEventTree->Branch("cc_mu_kde", &cc_mu_kde,"cc_mu_kde/F");
  fEventTree->Branch("cc_mu_trm_dqdx", &cc_mu_trm_dqdx,"cc_mu_trm_dqdx/F");
  fEventTree->Branch("cc_mu_trm_dedx", &cc_mu_trm_dedx,"cc_mu_trm_dedx/F");
  fEventTree->Branch("longest_trkid", &longest_trkid,"longest_trkid/I");
  fEventTree->Branch("longest_trklen", &longest_trklen,"longest_trklen/F");
  fEventTree->Branch("Tr_pri_mu_pdg", &Tr_pri_mu_pdg,"Tr_pri_mu_pdg/I");
  fEventTree->Branch("pri_Mu_is", &pri_Mu_is,"pri_Mu_is/I");
  fEventTree->Branch("Tr_pri_st_k_is", &Tr_pri_st_k_is,"Tr_pri_st_k_is/I");
  fEventTree->Branch("Tr_K_Inelas", &Tr_K_Inelas,"Tr_K_Inelas/I");
  fEventTree->Branch("Tr_k_plen", &Tr_k_plen,"Tr_k_plen/F");
  fEventTree->Branch("Tr_k_theta", &Tr_k_theta,"Tr_k_theta/F");
  fEventTree->Branch("Tr_k_phi", &Tr_k_phi,"Tr_k_phi/F");
  fEventTree->Branch("Tr_dec_mu_is", &Tr_dec_mu_is,"Tr_dec_mu_is/I");
  fEventTree->Branch("Tr_dec_mu_pi_pdg", &Tr_dec_mu_pi_pdg,"Tr_dec_mu_pi_pdg/I");
  fEventTree->Branch("Tr_mu_plen", &Tr_mu_plen,"Tr_mu_plen/F");
  fEventTree->Branch("Tr_k_endE", &Tr_k_endE,"Tr_k_endE/F");
  fEventTree->Branch("Tr_mu_theta", &Tr_mu_theta,"Tr_mu_theta/F");
  fEventTree->Branch("Tr_mu_phi", &Tr_mu_phi,"Tr_mu_phi/F");
  fEventTree->Branch("Tr_k_inTPC", &Tr_k_inTPC,"Tr_k_inTPC/I");
  fEventTree->Branch("Tr_mu_inTPC", &Tr_mu_inTPC,"Tr_mu_inTPC/I");
  fEventTree->Branch("Tr_k_in_5_TPC", &Tr_k_in_5_TPC,"Tr_k_in_5_TPC/I");
  fEventTree->Branch("Tr_k_in_CC_TPC", &Tr_k_in_CC_TPC,"Tr_k_in_CC_TPC/I");
  fEventTree->Branch("Tr_mu_in_5_TPC", &Tr_mu_in_5_TPC,"Tr_mu_in_5_TPC/I");
  fEventTree->Branch("Tr_mu_in_CC_TPC", &Tr_mu_in_CC_TPC,"Tr_mu_in_CC_TPC/I");
  fEventTree->Branch("Tr_kmu_open_ang", &Tr_kmu_open_ang,"Tr_kmu_open_ang/F");
  fEventTree->Branch("vtx_5cm_mult", &vtx_5cm_mult,"vtx_5cm_mult/I");
  fEventTree->Branch("k_start_dedx", &k_start_dedx,"k_start_dedx/F");
  fEventTree->Branch("k_end_dedx", &k_end_dedx,"k_end_dedx/F");
  fEventTree->Branch("mu_start_dedx", &mu_start_dedx,"mu_start_dedx/F");
  fEventTree->Branch("mu_end_dedx", &mu_end_dedx,"mu_end_dedx/F");
  fEventTree->Branch("cut_1", &cut_1,"cut_1/I");
  fEventTree->Branch("cut_2", &cut_2,"cut_2/I");
  fEventTree->Branch("cut_3", &cut_3,"cut_3/I");
  fEventTree->Branch("cut_4", &cut_4,"cut_4/I");
  fEventTree->Branch("cut_5", &cut_5,"cut_5/I");
  fEventTree->Branch("cut_6", &cut_6,"cut_6/I");
  fEventTree->Branch("cut_7", &cut_7,"cut_7/I");
  fEventTree->Branch("cut_8", &cut_8,"cut_8/I");
  fEventTree->Branch("kinelas_has_traks", &kinelas_has_traks,"kinelas_has_traks/I");
  fEventTree->Branch("kinelas_reco_trkID", &kinelas_reco_trkID,"kinelas_reco_trkID/I");
  fEventTree->Branch("kinelas_tlen", &kinelas_tlen,"kinelas_tlen/F");
  fEventTree->Branch("True_kinelas_KE", &True_kinelas_KE,"True_kinelas_KE/F");
  fEventTree->Branch("True_kinelas_tlen", &True_kinelas_tlen,"True_kinelas_tlen/F");
}

void CCKaonFilter::analyze( const art::Event& evt){

  reset();  

  bool isMC = true;//!evt.isRealData();

  run=evt.run();
  subrun=evt.subRun();
  event=evt.id().event(); 
  cout << "Run " << run;
  cout << " Subrun " << subrun;
  cout << " Event " << event << endl;

  art::Handle< std::vector<simb::MCTruth> > mctruths;
  std::vector< art::Ptr<simb::MCParticle> > ptList;
  if (isMC) {

    // get MCTruth
    evt.getByLabel("generator", mctruths);
    if (mctruths->size()!=1) {
      std::cout << "Number of MCTruths objects in event " << mctruths->size() << std::endl;
      return;
    }
    simb::MCTruth mctruth = mctruths->at(0);

    // get MCParticles
    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
    if (evt.getByLabel(fLArG4ModuleLabel, mcParticleHandle)){
      art::fill_ptr_vector(ptList, mcParticleHandle); 
    }

    // true neutrino information
    true_nu_energy = mctruth.GetNeutrino().Nu().E();
    true_nu_pdg = mctruth.GetNeutrino().Nu().PdgCode();
    true_nu_ccnc = mctruth.GetNeutrino().CCNC();//0=CC,1=NC
    true_nu_mode = mctruth.GetNeutrino().Mode();//0=QE,1=RES,2=DIS
    true_nu_vtx_x = mctruth.GetNeutrino().Nu().Vx();
    true_nu_vtx_y = mctruth.GetNeutrino().Nu().Vy();
    true_nu_vtx_z = mctruth.GetNeutrino().Nu().Vz();
    true_nu_vtx_inTPC = isInsideVolume("TPC",true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);
    true_nu_vtx_in5cmTPC = isInsideVolume("5cmTPC",true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);
    true_nu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);
    true_nu_vtx_inOldCCInclusiveTPC = isInsideVolume("OldCCInclusiveTPC",true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);

    true_lepton_pdg = mctruth.GetNeutrino().Lepton().PdgCode();
    true_lepton_p = mctruth.GetNeutrino().Lepton().P();
    true_lepton_theta = mctruth.GetNeutrino().Lepton().Momentum().Theta();
    true_lepton_phi = mctruth.GetNeutrino().Lepton().Momentum().Phi();
    TVector3 true_lepton_pvector = mctruth.GetNeutrino().Lepton().Momentum().Vect();

    // find the the highest momentum k+
    double true_kaon_maxp = -1;
    int true_kaon_maxp_id = -1;
    int nkplus = 0;

    for (auto const& pPart : ptList) {
      if (pPart->Process()=="primary" && pPart->PdgCode()==321) {
        const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",pPart->TrackId());
        if(int(mc_truth->Origin())==1) { // make sure origin is neutrino beam
          if (true_kaon_maxp<pPart->P()) {
            true_kaon_maxp = pPart->P();
            true_kaon_maxp_id = pPart->TrackId();
            nkplus++;
          }
        }  
      }
    }

    // true signal interaction: nu CC with at least one k+
    if (true_nu_pdg==14 && true_nu_ccnc==0 && nkplus>0) {

      vector<int> decay_muplus_id;
      vector<int> decay_piplus_id;
      vector<int> decay_daughters;
      vector<int> inelastic_daughters;
      TVector3 true_kaon_pvector;

      for (auto const& pPart : ptList) {

        // primary kaon+
        if (pPart->TrackId()==true_kaon_maxp_id) {
          TLorentzVector mcstart, mcend;
          unsigned int pstarti, pendi;
          true_kaon_length = length(*pPart, mcstart, mcend, pstarti, pendi);
          true_kaon_ke = pPart->E()-pPart->Mass();
          true_kaon_theta = pPart->Momentum().Theta();
          true_kaon_phi = pPart->Momentum().Phi();
          true_kaon_lepton_angle = pPart->Momentum().Angle(true_lepton_pvector);
          true_kaon_pvector = pPart->Momentum().Vect();
          true_kaon_end_ke = pPart->EndE()-pPart->Mass();
          true_kaon_end_x = pPart->EndX();
          true_kaon_end_y = pPart->EndY();
          true_kaon_end_z = pPart->EndZ();
          true_kaon_end_inTPC = isInsideVolume("TPC",pPart->EndX(),pPart->EndY(),pPart->EndZ());
          true_kaon_end_in5cmTPC = isInsideVolume("5cmTPC",pPart->EndX(),pPart->EndY(),pPart->EndZ());
          true_kaon_end_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",pPart->EndX(),pPart->EndY(),pPart->EndZ());
          true_kaon_end_inOldCCInclusiveTPC = isInsideVolume("OldCCInclusiveTPC",pPart->EndX(),pPart->EndY(),pPart->EndZ());
        }

        // find kaon+ daughters
        else if (pPart->Mother()==true_kaon_maxp_id) {

          if (pPart->Process()=="Decay") {
            decay_daughters.push_back(pPart->PdgCode());
            if (pPart->PdgCode()==-13) {
              decay_muplus_id.push_back(pPart->TrackId());
            }
            else if (pPart->PdgCode()==211) {
              decay_piplus_id.push_back(pPart->TrackId());
            }
            cout << "Decay " << pPart->PdgCode() << endl;
          }

          if (pPart->Process()=="kaon+Inelastic") {
            inelastic_daughters.push_back(pPart->PdgCode());
            cout << "Inelastic " << pPart->PdgCode() << endl;
          }
        }

      }//MC particles loop

      int n_decay_muplus = decay_muplus_id.size();
      int n_decay_piplus = decay_piplus_id.size();
      int n_decay_numu   = count(decay_daughters.begin(),decay_daughters.end(),14);
      int n_decay_pi0    = count(decay_daughters.begin(),decay_daughters.end(),111);
      int n_inelastic_kplus = count(inelastic_daughters.begin(),inelastic_daughters.end(),321);
      int daughter_id = -1;

      if (decay_daughters.size()) {
        //kaon+ -> muon+ numu
        if (decay_daughters.size()==2 && n_decay_muplus==1 && n_decay_numu==1) {
          true_kaon_end_process = 0;
          daughter_id = decay_muplus_id[0];
        }
        //kaon+ -> muon+ nopion+
        else if (n_decay_muplus==1 && n_decay_piplus==0) {
          true_kaon_end_process = 1;
          daughter_id = decay_muplus_id[0];
        }
        //kaon+ -> pion+ + pi0
        else if (decay_daughters.size()==2 && n_decay_piplus==1 && n_decay_pi0==1) {
          true_kaon_end_process = 2;
          daughter_id = decay_piplus_id[0];
        }
        //kaon+ -> pion+ and nomuon+
        else if (n_decay_muplus==0 && n_decay_piplus==1) {
          true_kaon_end_process = 3;
          daughter_id = decay_piplus_id[0];
        }
        //other decays
        else {
          true_kaon_end_process = 4;
        }
      }

      if (inelastic_daughters.size()) {
        //inelasitic with kaon+
        if (n_inelastic_kplus>0) {
          true_kaon_end_process = 5;
        }
        //other inelastic
        else {
          true_kaon_end_process = 6;
        }
      }

      //this shouldn't happen but it does
      if (decay_daughters.size() && inelastic_daughters.size()) {
        cout << "Primary kaon+ decay and inelastic at the same time!!!" << endl;
        true_kaon_end_process = 7;
      }

      cout << "True kaon end process " << true_kaon_end_process << endl;

      //check kaon daughter if it exists
      if (daughter_id!=-1) {
        for (auto const& pPart : ptList) {
      
          // kaon+ daughter
          if (pPart->TrackId()==daughter_id) {
            TLorentzVector mcstart, mcend;
            unsigned int pstarti, pendi;
            true_kaon_daughter_length = length(*pPart, mcstart, mcend, pstarti, pendi);
            true_kaon_daughter_ke = pPart->E()-pPart->Mass();
            true_kaon_daughter_angle = pPart->Momentum().Angle(true_kaon_pvector);
            true_kaon_daughter_pdg = pPart->PdgCode();
            true_kaon_daughter_end_x = pPart->EndX();
            true_kaon_daughter_end_y = pPart->EndY();
            true_kaon_daughter_end_z = pPart->EndZ();
            true_kaon_daughter_end_inTPC = isInsideVolume("TPC",pPart->EndX(),pPart->EndY(),pPart->EndZ());
            true_kaon_daughter_end_in5cmTPC = isInsideVolume("5cmTPC",pPart->EndX(),pPart->EndY(),pPart->EndZ());
            true_kaon_daughter_end_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",pPart->EndX(),pPart->EndY(),pPart->EndZ());
            true_kaon_daughter_end_inOldCCInclusiveTPC = isInsideVolume("OldCCInclusiveTPC",pPart->EndX(),pPart->EndY(),pPart->EndZ());
            break;
          }
      
        }//MC particles loop
      }//kaon daughter exists

    }//is true signal

  }//isMC

  // Check if event passed the NuCC inclusive filter
  reco_nu_cc_filter = false;
  art::InputTag trigResInputTag("TriggerResults","","OverlayFiltersPostStage2"); // the last is the name of process where the filters were run
  art::ValidHandle<art::TriggerResults> trigRes = evt.getValidHandle<art::TriggerResults>(trigResInputTag);
  fhicl::ParameterSet pset;
  if (!fhicl::ParameterSetRegistry::get(trigRes->parameterSetID(), pset)) { throw cet::exception("PSet Not Found???"); }
  std::vector<std::string> trigger_path_names = pset.get<std::vector<std::string> >("trigger_paths", {});
  if (trigger_path_names.size()!=trigRes->size()) { throw cet::exception("Size mismatch???"); }
  for (size_t itp=0;itp<trigRes->size();itp++) {
    //cout << "Filter name " << trigger_path_names.at(itp) << " decision=" << trigRes->at(itp).accept() << endl;
    if (trigger_path_names.at(itp)=="NuCC") {
      reco_nu_cc_filter = trigRes->at(itp).accept();
    }
  }

  if (!reco_nu_cc_filter) {
    std::cout << "Event didn't pass NuCC filter" << std::endl;
    fEventTree->Fill();
    return;
  }

  std::cout << "Event passed NuCC filter" << std::endl;

  // Collect all recontructed particles
  art::Handle<std::vector<recob::PFParticle>> pfparticles;
  evt.getByLabel(m_pfp_producer, pfparticles);
  if (pfparticles->size()==0) {
    std::cout << "No PFParticles found" << std::endl;
    fEventTree->Fill();
    return;
  }
  std::cout << "Number of PFParticles " << pfparticles->size() << std::endl;

  // Get PFParticle associations
  art::FindManyP<anab::T0> pfp_muon_assn(pfparticles, evt, "NuCCproducer");
  if(!pfp_muon_assn.isValid()){
    cout << "PFParticle-T0 associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  art::FindManyP<recob::Track> pfparticleTrackAssn(pfparticles, evt, "pandora");
  if(!pfparticleTrackAssn.isValid()){
    cout << "PFParticle-Track associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  art::FindManyP<recob::Vertex> pfparticleVertexAssn(pfparticles, evt, "pandora");
  if(!pfparticleVertexAssn.isValid()){
    cout << "PFParticle-Vertex associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  // Find recontructed neutrino (there should be one)
  lar_pandora::PFParticleVector pfneutrinos(0);
  for (unsigned int i=0; i<pfparticles->size(); ++i) {

    art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);

    if (pfparticle->IsPrimary() && pfparticle->PdgCode()==14) {
      pfneutrinos.push_back(pfparticle);
    }

  }

  if (pfneutrinos.size() != 1) {
    cout << "Number of neutrinos is not one" << endl;
    cout << "Number of neutrinos " << pfneutrinos.size() << endl;
    fEventTree->Fill();
    return;
  }

  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
  cout << "Found one neutrino";
  cout << " ID " << pfnu->Self();
  cout << " PDG " << pfnu->PdgCode();
  cout << " Number of daughters " << pfnu->Daughters().size() << endl;

  reco_nu_vtx_x = pfparticleVertexAssn.at(pfnu.key()).front()->position().X();
  reco_nu_vtx_y = pfparticleVertexAssn.at(pfnu.key()).front()->position().Y();
  reco_nu_vtx_z = pfparticleVertexAssn.at(pfnu.key()).front()->position().Z();
  cout << "Neutrino vertex (x,y,z) = " << reco_nu_vtx_x << ", " << reco_nu_vtx_y << ", " << reco_nu_vtx_z << endl;
  reco_nu_vtx_inTPC = isInsideVolume("TPC",reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z);
  reco_nu_vtx_in5cmTPC = isInsideVolume("5cmTPC",reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z);
  reco_nu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z);

  // Find CC muon and daughters
  lar_pandora::PFParticleVector pfmuons(0);
  vector<int> reco_nu_daughters_id(0);

  for (unsigned int i=0; i<pfparticles->size(); ++i) {

    art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);

    // look at particles with neutrino parent and one associated track
    if (pfparticle->Parent()==pfnu->Self() && pfparticleTrackAssn.at(i).size()==1) {

      art::Ptr<recob::Track> track = pfparticleTrackAssn.at(i).front();
      cout << "Neutrino daughter";
      cout << " ID: " << pfparticle->Self();
      cout << " PDG: " << pfparticle->PdgCode();
      cout << " Track key: " << track.key() << endl;

      reco_nu_daughters_id.push_back(track.key());

      // muon has a T0 associated
      if (pfp_muon_assn.at(i).size()==1) {
        pfmuons.push_back(pfparticle);
      }

    }

  }

  reco_nu_ndaughters = reco_nu_daughters_id.size();
  std::cout << "Number of neutrino daughters with one associated track " << reco_nu_ndaughters << std::endl;

  if (pfmuons.size()!=1) {
    std::cout << "Number of CC inclusive muons is not 1: " << pfmuons.size() << std::endl;
    fEventTree->Fill();
    return;
  }

  art::Ptr<recob::PFParticle> pfmuon = pfmuons.front();
  art::Ptr<recob::Track> trkmuon = pfparticleTrackAssn.at(pfmuon.key()).front();
  cout << "Found CC muon";
  cout << " ID " << pfmuon->Self();
  cout << " PDG " << pfmuon->PdgCode();
  cout << " Parent " << pfmuon->Parent();
  cout << " Track key " << trkmuon.key() << endl;

  reco_ccmu_vtx_x = pfparticleVertexAssn.at(pfmuon.key()).front()->position().X();
  reco_ccmu_vtx_y = pfparticleVertexAssn.at(pfmuon.key()).front()->position().Y();
  reco_ccmu_vtx_z = pfparticleVertexAssn.at(pfmuon.key()).front()->position().Z();
  cout << "CC muon start (x,y,z) = " << reco_ccmu_vtx_x << ", " << reco_ccmu_vtx_y << ", " << reco_ccmu_vtx_z << endl;
  reco_ccmu_vtx_inTPC = isInsideVolume("TPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);
  reco_ccmu_vtx_in5cmTPC = isInsideVolume("5cmTPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);
  reco_ccmu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);

  // get hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitlist, hitListHandle);
  }
  cout << "Found " << hitlist.size() << " hits" << endl;

  // get tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) {
    art::fill_ptr_vector(tracklist, trackListHandle);
  }
  cout << "Found " << tracklist.size() << " tracks" << endl;

  // get track associations
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  if(!fmcal.isValid()){
    cout << "Track-Calorimetry associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, evt, fHitTrackAssns);
  if(!hits_from_tracks.isValid()){
    cout << "Track-Hit associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  art::FindManyP<anab::ParticleID> trackPIDAssn(trackListHandle, evt, fPIDLabel);
  if(!trackPIDAssn.isValid()){
    cout << "Track PID associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  // find track multiplicity around vertex
  int NTracks=tracklist.size();
  int ntracks = 0;

  //selection cuts
  std::vector<int> kaon_can_trkID;
  std::vector<int> muon_can_trkID;

  // loop over tracks again to look for kaon track
  cout << "Looking for kaon tracks from " << NTracks << " tracks available" << endl;

  // loop over nu daughters
  for (int i=0; i<reco_nu_ndaughters; i++) {

    art::Ptr<recob::Track> ptrack(trackListHandle,reco_nu_daughters_id[i]);
    const recob::Track& track = *ptrack;

    // skip cc muon track
    if (ptrack.key()==trkmuon.key()) continue;

    TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
    //double st_vtx=TMath::Sqrt((reco_ccmu_vtx_x-pos.X())*(reco_ccmu_vtx_x-pos.X()) +
    //                          (reco_ccmu_vtx_y-pos.Y())*(reco_ccmu_vtx_y-pos.Y()) +
    //                          (reco_ccmu_vtx_z-pos.Z())*(reco_ccmu_vtx_z-pos.Z()));
    double st_vtx=TMath::Sqrt((reco_nu_vtx_x-pos.X())*(reco_nu_vtx_x-pos.X()) +
                              (reco_nu_vtx_y-pos.Y())*(reco_nu_vtx_y-pos.Y()) +
                              (reco_nu_vtx_z-pos.Z())*(reco_nu_vtx_z-pos.Z()));

    // check distance to vertex
    //if (st_vtx<=10) { //7cm

    cut_1=1;	  

    reco_track_distance[ntracks] = st_vtx;

    // check track is contained?
    TVector3 end(track.End().X(),track.End().Y(),track.End().Z());

    reco_track_vtx_inTPC[ntracks] = isInsideVolume("TPC",pos.X(),pos.Y(),pos.Z());
    reco_track_vtx_in5cmTPC[ntracks] = isInsideVolume("5cmTPC",pos.X(),pos.Y(),pos.Z());
    reco_track_vtx_inCCInclusiveTPC[ntracks] = isInsideVolume("CCInclusiveTPC",pos.X(),pos.Y(),pos.Z());

    reco_track_end_inTPC[ntracks] = isInsideVolume("TPC",end.X(),end.Y(),end.Z());
    reco_track_end_in5cmTPC[ntracks] = isInsideVolume("5cmTPC",end.X(),end.Y(),end.Z());
    reco_track_end_inCCInclusiveTPC[ntracks] = isInsideVolume("CCInclusiveTPC",end.X(),end.Y(),end.Z());

    // check track length
    double trklen=track.Length();
    reco_track_length[ntracks] = trklen;

    // check kaon track start and end distance from vertex
    //double end_dis=TMath::Sqrt((reco_ccmu_vtx_x-end.X())*(reco_ccmu_vtx_x-end.X()) +
    //                           (reco_ccmu_vtx_y-end.Y())*(reco_ccmu_vtx_y-end.Y()) +
    //                           (reco_ccmu_vtx_z-end.Z())*(reco_ccmu_vtx_z-end.Z()));
    double end_dis=TMath::Sqrt((reco_nu_vtx_x-end.X())*(reco_nu_vtx_x-end.X()) +
                               (reco_nu_vtx_y-end.Y())*(reco_nu_vtx_y-end.Y()) +
                               (reco_nu_vtx_z-end.Z())*(reco_nu_vtx_z-end.Z()));
    
    reco_track_dir[ntracks] = (st_vtx<end_dis);

    // check calorimetry
    int hits_p0=0;
    int hits_p1=0;
    int hits_p2=0;

    std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(ptrack.key());
    
    for(unsigned int ical=0; ical<calos.size(); ++ical){
      if(!calos[ical]) continue;
      if(!calos[ical]->PlaneID().isValid) continue;
      int planenum = calos[ical]->PlaneID().Plane;
      if(planenum<0||planenum>2) continue; 
      if(planenum==0) hits_p0=calos[ical]->dEdx().size();
      if(planenum==1) hits_p1=calos[ical]->dEdx().size();
      if(planenum==2) hits_p2=calos[ical]->dEdx().size();
    }
    
    reco_track_nhits0[ntracks] = hits_p0;
    reco_track_nhits1[ntracks] = hits_p1;
    reco_track_nhits2[ntracks] = hits_p2;

    // check PID
    if(trackPIDAssn.isValid()){
      std::vector<art::Ptr<anab::ParticleID>> trackPID=trackPIDAssn.at(ptrack.key());
      if(trackPID.size()!=0){
        std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
        for(size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
          anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
          if(AlgScore.fAlgName == "Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF){
            if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
              if(AlgScore.fAssumedPdg==321)  reco_track_chi2k[ntracks]=AlgScore.fValue;
              if(AlgScore.fAssumedPdg==2212) reco_track_chi2p[ntracks]=AlgScore.fValue;
              if(AlgScore.fAssumedPdg==211)  reco_track_chi2pi[ntracks]=AlgScore.fValue;
              if(AlgScore.fAssumedPdg==13)   reco_track_chi2mu[ntracks]=AlgScore.fValue;
            }
          } 
        }
      }
    }

    // find true matched particle
    simb::MCParticle const* matched_mcparticle = NULL;
    std::unordered_map<int,double> trkide;
    double maxe=-1, tote=0;
    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
    std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
    art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle,evt,fHitTruthAssns);

    for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
      particle_vec.clear(); match_vec.clear();
      particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
      for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
        trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
        tote += match_vec[i_p]->energy;
        if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
          maxe = trkide[ particle_vec[i_p]->TrackId() ];
          matched_mcparticle = particle_vec[i_p];
        }
      }
    }

    if(matched_mcparticle){
      reco_track_pdg[ntracks] = matched_mcparticle->PdgCode();
      const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
      reco_track_origin[ntracks]=int(mc_truth->Origin());
      reco_track_primary[ntracks]=matched_mcparticle->Process()=="primary";
    }

    // check if there is a track at the end
    int ndaughters = 0;

    for(int j=0; j<NTracks; j++){
    
      art::Ptr<recob::Track> ptrack2(trackListHandle,j);
      const recob::Track& track2 = *ptrack2;
    
      // skip cc muon track
      if (ptrack2.key()==trkmuon.key()) continue;

      // skip primary track
      if (ptrack2.key()==ptrack.key()) continue;

      TVector3 pos2(track2.Vertex().X(),track2.Vertex().Y(),track2.Vertex().Z());

      double track2_distance=TMath::Sqrt((end.X()-pos2.X())*(end.X()-pos2.X()) +
                                         (end.Y()-pos2.Y())*(end.Y()-pos2.Y()) +
                                         (end.Z()-pos2.Z())*(end.Z()-pos2.Z()));
      // check distance to vertex
      if (track2_distance<10) { //7cm

        reco_track_daughter_distance[ntracks][ndaughters] = track2_distance;

        TVector3 end2(track2.End().X(),track2.End().Y(),track2.End().Z());
      
        reco_track_daughter_vtx_inTPC[ntracks][ndaughters] = isInsideVolume("TPC",pos2.X(),pos2.Y(),pos2.Z());
        reco_track_daughter_vtx_in5cmTPC[ntracks][ndaughters] = isInsideVolume("5cmTPC",pos2.X(),pos2.Y(),pos2.Z());
        reco_track_daughter_vtx_inCCInclusiveTPC[ntracks][ndaughters] = isInsideVolume("CCInclusiveTPC",pos2.X(),pos2.Y(),pos2.Z());
      
        reco_track_daughter_end_inTPC[ntracks][ndaughters] = isInsideVolume("TPC",end2.X(),end2.Y(),end2.Z());
        reco_track_daughter_end_in5cmTPC[ntracks][ndaughters] = isInsideVolume("5cmTPC",end2.X(),end2.Y(),end2.Z());
        reco_track_daughter_end_inCCInclusiveTPC[ntracks][ndaughters] = isInsideVolume("CCInclusiveTPC",end2.X(),end2.Y(),end2.Z());

        // check track length
        reco_track_daughter_length[ntracks][ndaughters] = track2.Length();

        double start_dis=TMath::Sqrt((reco_nu_vtx_x-pos2.X())*(reco_nu_vtx_x-pos2.X()) +
                                     (reco_nu_vtx_y-pos2.Y())*(reco_nu_vtx_y-pos2.Y()) +
                                     (reco_nu_vtx_z-pos2.Z())*(reco_nu_vtx_z-pos2.Z()));    

        reco_track_daughter_vtx_distance[ntracks][ndaughters] = start_dis;

        // check calorimetry
        int hits2_p0=0;
        int hits2_p1=0;
        int hits2_p2=0;
      
        std::vector<art::Ptr<anab::Calorimetry>> calos2=fmcal.at(ptrack2.key());
        
        for(unsigned int ical=0; ical<calos2.size(); ++ical){
          if(!calos2[ical]) continue;
          if(!calos2[ical]->PlaneID().isValid) continue;
          int planenum2 = calos2[ical]->PlaneID().Plane;
          if(planenum2<0||planenum2>2) continue; 
          if(planenum2==0) hits2_p0=calos2[ical]->dEdx().size();
          if(planenum2==1) hits2_p1=calos2[ical]->dEdx().size();
          if(planenum2==2) hits2_p2=calos2[ical]->dEdx().size();
        }
        
        reco_track_daughter_nhits0[ntracks][ndaughters] = hits2_p0;
        reco_track_daughter_nhits1[ntracks][ndaughters] = hits2_p1;
        reco_track_daughter_nhits2[ntracks][ndaughters] = hits2_p2;

        // check PID
        if(trackPIDAssn.isValid()){
          std::vector<art::Ptr<anab::ParticleID>> trackPID=trackPIDAssn.at(ptrack2.key());
          if(trackPID.size()!=0){
            std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
            for(size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
              anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
              if(AlgScore.fAlgName == "Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF){
                if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                  if(AlgScore.fAssumedPdg==321)  reco_track_daughter_chi2k[ntracks][ndaughters]=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==2212) reco_track_daughter_chi2p[ntracks][ndaughters]=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==211)  reco_track_daughter_chi2pi[ntracks][ndaughters]=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==13)   reco_track_daughter_chi2mu[ntracks][ndaughters]=AlgScore.fValue;
                }
              } 
            }
          }
        }
      
        // find true matched particle
        simb::MCParticle const* matched_mcparticle2 = NULL;
        std::unordered_map<int,double> trkide2;
        double maxe2=-1, tote2=0;
        std::vector<simb::MCParticle const*> particle_vec2;
        std::vector<anab::BackTrackerHitMatchingData const*> match_vec2;
        std::vector<art::Ptr<recob::Hit>> hits_from_track2 = hits_from_tracks.at(ptrack2.key());
        art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit2(hitListHandle,evt,fHitTruthAssns);
      
        for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
          particle_vec2.clear(); match_vec2.clear();
          particles_per_hit2.get(hits_from_track[i_h].key(),particle_vec2,match_vec2);
          for(size_t i_p=0; i_p<particle_vec2.size(); ++i_p){
            trkide[ particle_vec2[i_p]->TrackId() ] += match_vec2[i_p]->energy;
            tote2 += match_vec2[i_p]->energy;
            if( trkide2[ particle_vec2[i_p]->TrackId() ] > maxe2 ){
              maxe2 = trkide2[ particle_vec2[i_p]->TrackId() ];
              matched_mcparticle2 = particle_vec2[i_p];
            }
          }
        }
      
        if(matched_mcparticle2){
          reco_track_daughter_pdg[ntracks][ndaughters] = matched_mcparticle2->PdgCode();
          const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle2->TrackId());
          reco_track_daughter_origin[ntracks][ndaughters]=int(mc_truth->Origin());
          reco_track_daughter_primary[ntracks][ndaughters]=matched_mcparticle2->Process()=="primary";
        }

        ndaughters++;

      }

      reco_track_ndaughters[ntracks] = ndaughters;

    }

    ntracks++;
    
    if (reco_track_end_inTPC && reco_track_vtx_inTPC) {
      cut_2=1;
      if(hits_p2>=5){
        cut_3=1;     
        if(trklen>=5){ // initiall this was 5 cm
          cut_4=1;
          cut_5=1;   
          if (st_vtx<=7 && (st_vtx<end_dis)) { // 2.5 // 5
            cut_6=1;
            //save track ID for kaon candidate
            kaon_can_trkID.push_back(ptrack.key());
            cout << "track candidate " << ptrack.key() << endl;
          } // vtx distance to the K track < 2.5
        } // trklen > 7 cm
      } // more than 5 hits in the collection plane
    }// track is contained

    //} // vertex distance
  } // loop over K  reco trks	

  reco_ntracks = ntracks;

  cout << "Number of kaon candidate tracks " << reco_nu_ndaughters << endl;
  cout << "Number of kaon candidate tracks " << reco_ntracks << endl;
  cout << "Number of kaon candidate tracks " << kaon_can_trkID.size() << endl;

  for (int i=0;i<reco_ntracks;i++) {
    for (int j=0;j<reco_track_ndaughters[i];j++) {
      cout << "kaon " << i << " daughter " << j << " distance " << reco_track_daughter_distance[i][j] << endl;
    }
  }

  if (reco_ntracks>0) {

    // loop over tracks again to look for muon track
    cout << "Looking for muon tracks from " << NTracks << " tracks available" << endl;
    for(int j=0; j<NTracks; j++){

      art::Ptr<recob::Track> ptrack(trackListHandle,j);
      const recob::Track& track = *ptrack;

      // skip cc muon track
      if (ptrack.key()==trkmuon.key()) continue;

      TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
      TVector3 end(track.End().X(),track.End().Y(),track.End().Z());

      // check track is contained?
      if (isInsideVolume("TPC",pos.X(),pos.Y(),pos.Z()) &&
          isInsideVolume("TPC",end.X(),end.Y(),end.Z())) {

        cut_7=1;

        // check track start is far away from vertex
        //float start_dis=TMath::Sqrt((reco_ccmu_vtx_x-pos.X())*(reco_ccmu_vtx_x-pos.X()) +
        //                            (reco_ccmu_vtx_y-pos.Y())*(reco_ccmu_vtx_y-pos.Y()) +
        //                            (reco_ccmu_vtx_z-pos.Z())*(reco_ccmu_vtx_z-pos.Z()));    
        float start_dis=TMath::Sqrt((reco_nu_vtx_x-pos.X())*(reco_nu_vtx_x-pos.X()) +
                                    (reco_nu_vtx_y-pos.Y())*(reco_nu_vtx_y-pos.Y()) +
                                    (reco_nu_vtx_z-pos.Z())*(reco_nu_vtx_z-pos.Z()));    
        if(start_dis>5){   // 5

          cut_8=1;   

          // check hits in collection plane (plane 2)
          int hits_p2=0;
          std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(ptrack.key());

          for(unsigned int ical=0; ical<calos.size(); ++ical){
            if(!calos[ical]) continue;
            if(!calos[ical]->PlaneID().isValid) continue;
            int planenum = calos[ical]->PlaneID().Plane;
            if(planenum<0||planenum>2) continue; 
            if(planenum==2) hits_p2=calos[ical]->dEdx().size();
          }

          if(hits_p2>=10){

            cut_9=1;    

            // check track length
            float trklen=track.Length();
            if(trklen>5){ // trklen>40 && trklen<60 // trklen>30 && trklen<70 // trklen>20 && trklen<70

              cut_10=1;	
              cut_11=1;   

              //save track ID for muon candidate
              muon_can_trkID.push_back(ptrack.key());

            } // muon track length cut

          } // hits > 10
        } // track is 5 cm awary from the vertex 
      } // track is contained in the full TPC
    } // loop over mu reco trks

    //cout << "Number of candidate muon tracks found: " << muon_can_trkID.size() << endl;

    //find muon tracks that are 5cm from kaon track end
    std::vector<int> mu_mindis_can;
    std::vector<int> k_mindis_can;
    std::vector<float> close_dis;

    for (unsigned int i=0; i<kaon_can_trkID.size(); i++) {

      art::Ptr<recob::Track> ptrack_k(trackListHandle,kaon_can_trkID[i]);
      const recob::Track& track_k = *ptrack_k;
      TVector3 end_k(track_k.End().X(),track_k.End().Y(),track_k.End().Z());

      for (unsigned int j=0; j<muon_can_trkID.size(); j++) {

        if (kaon_can_trkID[i]==muon_can_trkID[j]) {
          continue;
        } // two track ID's are not identical

        art::Ptr<recob::Track> ptrack_mu(trackListHandle,muon_can_trkID[j]);
        const recob::Track& track_mu = *ptrack_mu;
        TVector3 pos_mu(track_mu.Vertex().X(),track_mu.Vertex().Y(),track_mu.Vertex().Z());

        double k_en_mu_st=TMath::Sqrt((end_k.X()-pos_mu.X())*(end_k.X()-pos_mu.X()) +
                                      (end_k.Y()-pos_mu.Y())*(end_k.Y()-pos_mu.Y()) +
                                      (end_k.Z()-pos_mu.Z())*(end_k.Z()-pos_mu.Z()));

        if (k_en_mu_st<=5) { // 2.5 // 5
          mu_mindis_can.push_back(muon_can_trkID[j]);
          k_mindis_can.push_back(kaon_can_trkID[i]);
          close_dis.push_back(k_en_mu_st);
        } // k/mu vtx < 2.5 cm

      } // loop over mu candidates
    } // loop over k candidates

    cout << "Number of candidate muon tracks found: " << close_dis.size() << endl;
    reco_track_ndaughters2 = close_dis.size();

    //found muon track 5cm away from kaon track end
    if (mu_mindis_can.size() && k_mindis_can.size()) {

      cut_12=1;

      //find muon track closest to a kaon track end
      float min2_dis=5e20;
      int track_index=-1;
      for(unsigned int k=0; k<close_dis.size(); k++){
        if(close_dis[k]<min2_dis){
          min2_dis=close_dis[k];
          track_index=k;
        } 
      }

      //save angle and distances between kaon and muon tracks
      k_can_trkid=k_mindis_can[track_index];
      mu_can_trkid=mu_mindis_can[track_index];
      k_mu_can_dis=min2_dis;

      art::Ptr<recob::Track> ptrack_k(trackListHandle,k_can_trkid);
      const recob::Track& track_k = *ptrack_k;

      TVector3 dir_end_k(track_k.EndDirection().X(),track_k.EndDirection().Y(),track_k.EndDirection().Z());
      TVector3 pos_k(track_k.Vertex().X(),track_k.Vertex().Y(),track_k.Vertex().Z());
      TVector3 end_k(track_k.End().X(),track_k.End().Y(),track_k.End().Z());

      art::Ptr<recob::Track> ptrack_mu(trackListHandle,mu_can_trkid);
      const recob::Track& track_mu = *ptrack_mu;

      TVector3 dir_start_mu(track_mu.VertexDirection().X(),track_mu.VertexDirection().Y(),track_mu.VertexDirection().Z());
      TVector3 pos_mu(track_mu.Vertex().X(),track_mu.Vertex().Y(),track_mu.Vertex().Z());
      TVector3 end_mu(track_mu.End().X(),track_mu.End().Y(),track_mu.End().Z());

      float k_cosA=dir_end_k.X();float k_cosB=dir_end_k.Y();float k_cosC=dir_end_k.Z();
      float mu_cosA=dir_start_mu.X();float mu_cosB=dir_start_mu.Y();float mu_cosC=dir_start_mu.Z();
      k_mu_open_angle=(TMath::ACos(k_cosA*mu_cosA + k_cosB*mu_cosB + k_cosC*mu_cosC));

      //k_vtx_dis=TMath::Sqrt((pos_k.X()-reco_ccmu_vtx_x)*(pos_k.X()-reco_ccmu_vtx_x) +
      //                      (pos_k.Y()-reco_ccmu_vtx_y)*(pos_k.Y()-reco_ccmu_vtx_y) +
      //                      (pos_k.Z()-reco_ccmu_vtx_z)*(pos_k.Z()-reco_ccmu_vtx_z));
      k_vtx_dis=TMath::Sqrt((pos_k.X()-reco_nu_vtx_x)*(pos_k.X()-reco_nu_vtx_x) +
                            (pos_k.Y()-reco_nu_vtx_y)*(pos_k.Y()-reco_nu_vtx_y) +
                            (pos_k.Z()-reco_nu_vtx_z)*(pos_k.Z()-reco_nu_vtx_z));

      cout << "Found muon track closest to kaon track end" << endl;
      cout << "kaon track key " << k_can_trkid << endl;
      cout << "muon track key " << mu_can_trkid << endl;
      cout << "distance between kaon and vertex " << k_vtx_dis << " cm" << endl;
      cout << "distance between muon and kaon " << k_mu_can_dis << " cm" << endl;

      ////////////////////////////////////////////////////////////////////// Doing End point matching for reco and true tracks /////////////////////////////////////
/*
      cout << "Look at true kaon and muon information" << endl;
      int New_True_Kaon_ID=-9999;
      for(auto const& pPart : ptList){
        std::string pri("primary");
        bool isPrimary=pPart->Process()==pri;
        if(isPrimary && pPart->PdgCode()==321){
          const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",pPart->TrackId());
          int origin=int(mc_truth->Origin());
          if(origin==1){ 
            New_True_Kaon_ID=pPart->TrackId();
            //distance to kaon track
            float ksx=pPart->Vx();float ksy=pPart->Vy();float ksz=pPart->Vz();
            float kex=pPart->EndX();float key=pPart->EndY();float kez=pPart->EndZ();
            float st_st=TMath::Sqrt((pos_k.X()-ksx)*(pos_k.X()-ksx) + (pos_k.Y()-ksy)*(pos_k.Y()-ksy) + (pos_k.Z()-ksz)*(pos_k.Z()-ksz));
            float st_en=TMath::Sqrt((pos_k.X()-kex)*(pos_k.X()-kex) + (pos_k.Y()-key)*(pos_k.Y()-key) + (pos_k.Z()-kez)*(pos_k.Z()-kez));
            float en_st=TMath::Sqrt((end_k.X()-ksx)*(end_k.X()-ksx) + (end_k.Y()-ksy)*(end_k.Y()-ksy) + (end_k.Z()-ksz)*(end_k.Z()-ksz));
            float en_en=TMath::Sqrt((end_k.X()-kex)*(end_k.X()-kex) + (end_k.Y()-key)*(end_k.Y()-key) + (end_k.Z()-kez)*(end_k.Z()-kez));
            if((st_st<=5 || st_en<=5) && (en_st<=5 || en_en<=5)) mtc_k_5cm=1;
            if((st_st<=10 || st_en<=10) && (en_st<=10 || en_en<=10)) mtc_k_10cm=1;
            mtc_k_endE=pPart->EndE();	    
          }  
        }
        if(pPart->Process()=="Decay" && pPart->PdgCode()==13 && pPart->Mother()==New_True_Kaon_ID){
          const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",pPart->TrackId());
          int origin=int(mc_truth->Origin());
          if(origin==1){ 
            float musx=pPart->Vx();float musy=pPart->Vy();float musz=pPart->Vz();
            float muex=pPart->EndX();float muey=pPart->EndY();float muez=pPart->EndZ();
            float st_st=TMath::Sqrt((pos_mu.X()-musx)*(pos_mu.X()-musx) + (pos_mu.Y()-musy)*(pos_mu.Y()-musy) + (pos_mu.Z()-musz)*(pos_mu.Z()-musz));
            float st_en=TMath::Sqrt((pos_mu.X()-muex)*(pos_mu.X()-muex) + (pos_mu.Y()-muey)*(pos_mu.Y()-muey) + (pos_mu.Z()-muez)*(pos_mu.Z()-muez));
            float en_st=TMath::Sqrt((end_mu.X()-musx)*(end_mu.X()-musx) + (end_mu.Y()-musy)*(end_mu.Y()-musy) + (end_mu.Z()-musz)*(end_mu.Z()-musz));
            float en_en=TMath::Sqrt((end_mu.X()-muex)*(end_mu.X()-muex) + (end_mu.Y()-muey)*(end_mu.Y()-muey) + (end_mu.Z()-muez)*(end_mu.Z()-muez));
            mtc_mu_pid=pPart->PdgCode();
            if((st_st<=5 || st_en<=5) && (en_st<=5 || en_en<=5)) mtc_mu_5cm=1;
            if((st_st<=10 || st_en<=10) && (en_st<=10 || en_en<=10)) mtc_mu_10cm=1;
          }
        }
      }// loop over plist
*/
      ////////////////////////////////////////////////////////////////// End of end point matching //////////////////////////////
    } // minimum distance vectors are filled
    //} // kaon and muon ID vectors are filled

    //////////////////////////////////////////////////////// RECO-TRUTH MATCHING /////////////////////////////////////////////////////

    if(k_can_trkid!=-9999){

      cout << "Find true particle for kaon track " << k_can_trkid << endl;
      simb::MCParticle const* matched_mcparticle = NULL;
      std::unordered_map<int,double> trkide;
      double maxe=-1, tote=0;
      std::vector<simb::MCParticle const*> particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(k_can_trkid);
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle,evt,fHitTruthAssns);

      for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
        particle_vec.clear(); match_vec.clear();
        particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
        for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
          trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
          tote += match_vec[i_p]->energy;
          if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
            maxe = trkide[ particle_vec[i_p]->TrackId() ];
            matched_mcparticle = particle_vec[i_p];
          }
        }
      }

      std::string pri("primary");
      if(matched_mcparticle){
        int k_geant_id=matched_mcparticle->TrackId();
        k_geant_ID=k_geant_id;
        const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",k_geant_id);
        k_origin=int(mc_truth->Origin());
        k_pdg=matched_mcparticle->PdgCode();
        cout << "Found truth particle " << k_pdg << endl;
        if(matched_mcparticle->Process()==pri) k_isPri=1;
        k_endE=matched_mcparticle->EndE();
        if(k_origin==1 && k_pdg==321 && k_endE<=510 && k_isPri==1) k_ness=1;
      }

      cout << "Check kaon track position" << endl;
      art::Ptr<recob::Track> ptrack(trackListHandle,k_can_trkid);
      const recob::Track& track = *ptrack;
      TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
      TVector3 end(track.End().X(),track.End().Y(),track.End().Z());
      TVector3 dir_start(track.VertexDirection().X(),track.VertexDirection().Y(),track.VertexDirection().Z());	
      k_plen=track.Length();
      k_phi=dir_start.Phi();
      k_theta=dir_start.Theta();
      if (isInsideVolume("5cmTPC",pos.X(),pos.Y(),pos.Z()) &&
          isInsideVolume("5cmTPC",end.X(),end.Y(),end.Z())) {
        k_in_5_TPC=1;
      }
      if (isInsideVolume("CCInclusiveTPC",pos.X(),pos.Y(),pos.Z()) &&
          isInsideVolume("CCInclusiveTPC",end.X(),end.Y(),end.Z())) {
        k_in_CC_TPC=1; 
      }
      //kaon_vtx_dis=TMath::Sqrt((reco_ccmu_vtx_x-pos.X())*(reco_ccmu_vtx_x-pos.X()) +
      //                         (reco_ccmu_vtx_y-pos.Y())*(reco_ccmu_vtx_y-pos.Y()) +
      //                         (reco_ccmu_vtx_z-pos.Z())*(reco_ccmu_vtx_z-pos.Z()));
      kaon_vtx_dis=TMath::Sqrt((reco_nu_vtx_x-pos.X())*(reco_nu_vtx_x-pos.X()) +
                               (reco_nu_vtx_y-pos.Y())*(reco_nu_vtx_y-pos.Y()) +
                               (reco_nu_vtx_z-pos.Z())*(reco_nu_vtx_z-pos.Z()));
      cout << "kaon and vertex distance " << kaon_vtx_dis << endl;

      cout << "Check kaon track calorimetry" << endl;
      float KE_K_p2=0;

      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(k_can_trkid);

      for(unsigned int ical=0; ical<calos.size(); ++ical){
        if(!calos[ical]) continue;
        if(!calos[ical]->PlaneID().isValid) continue;
        int planenum = calos[ical]->PlaneID().Plane;
        if(planenum<0||planenum>2) continue; 
        if(planenum==fplane){
          const size_t NHits = calos[ical] -> dEdx().size();
          k_hit=int(NHits);
          k_range=calos[ical]->Range();
          float bin_size=0;
          if(k_range<10) bin_size=float(k_range)/2;
          if(k_range>=10) bin_size=5.0;
          std::vector<float> st_vec,en_vec;
          float rr_x0=-1;float rr_y0=-1;float rr_z0=-1;
          for(size_t iHit = 0; iHit < NHits; ++iHit){
            const auto& TrkPos=(calos[ical] -> XYZ())[iHit];
            if(iHit==0){
              rr_x0=TrkPos.X();rr_y0=TrkPos.Y();rr_z0=TrkPos.Z();
            }
            if((calos[ical]->ResidualRange())[iHit]<=bin_size) en_vec.push_back((calos[ical] -> dEdx())[iHit]);
            if((k_range-(calos[ical]->ResidualRange())[iHit])<=bin_size) st_vec.push_back((calos[ical] -> dEdx())[iHit]);
            if(iHit>=1){
              if((calos[ical]->ResidualRange())[iHit]>=0 && (calos[ical]->ResidualRange())[iHit-1]>=0)
                KE_K_p2=KE_K_p2+TMath::Abs(((calos[ical]->ResidualRange())[iHit]-(calos[ical]->ResidualRange())[iHit-1])*(calos[ical] -> dEdx())[iHit]);
            }
            if(iHit==0){
              if((calos[ical]->ResidualRange())[iHit]>=0 && (calos[ical]->ResidualRange())[iHit+1]>=0)
                KE_K_p2=KE_K_p2+TMath::Abs(((calos[ical]->ResidualRange())[iHit]-(calos[ical]->ResidualRange())[iHit+1])*(calos[ical] -> dEdx())[iHit]);
            }
          } 

          if(st_vec.size() && en_vec.size()){
            float en_median=TMath::Median(en_vec.size(),&en_vec[0]);
            float st_median=TMath::Median(st_vec.size(),&st_vec[0]);
            float st_hit0=TMath::Sqrt((pos.X()-rr_x0)*(pos.X()-rr_x0) + (pos.Y()-rr_y0)*(pos.Y()-rr_y0) + (pos.Z()-rr_z0)*(pos.Z()-rr_z0));
            float en_hit0=TMath::Sqrt((end.X()-rr_x0)*(end.X()-rr_x0) + (end.Y()-rr_y0)*(end.Y()-rr_y0) + (end.Z()-rr_z0)*(end.Z()-rr_z0));

            if(st_hit0<en_hit0){
              k_start_dedx=en_median;
              k_end_dedx=st_median;
            }

            if(en_hit0<st_hit0){
              k_end_dedx=en_median;
              k_start_dedx=st_median;
            }

            if(en_hit0==st_hit0){
              k_end_dedx=en_median;
              k_start_dedx=st_median;
            }

            if(en_median>st_median){ 
              k_large_dedx=en_median;
              k_small_dedx=st_median;
            }
            if(en_median<st_median){ 
              k_large_dedx=st_median;
              k_small_dedx=en_median;
            }
            if(en_median==st_median){ 
              k_large_dedx=en_median;
              k_small_dedx=st_median;
            }
          }
        } // grab plane 2
      } // loop over calorimetry objects

      k_KE=KE_K_p2;

      cout << "Check kaon track PID" << endl;
      if(trackPIDAssn.isValid()){
        std::vector<art::Ptr<anab::ParticleID>> trackPID=trackPIDAssn.at(k_can_trkid);
        if(trackPID.size()!=0){
          std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
          float k_Brg_F_mu=-9999;float k_Brg_F_p=-9999;float k_Brg_F_pi=-9999;float k_Brg_F_k=-9999;float k_Brg_F_mip=-9999;
          float k_Brg_B_mu=-9999;float k_Brg_B_p=-9999;float k_Brg_B_pi=-9999;float k_Brg_B_k=-9999;float k_Brg_B_mip=-9999;
          for(size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
            anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
            //int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneID);
            //int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
            if(AlgScore.fAlgName == "Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                if(AlgScore.fAssumedPdg==13) k_chi_mu=AlgScore.fValue;
                if(AlgScore.fAssumedPdg==2212) k_chi_p=AlgScore.fValue;
                if(AlgScore.fAssumedPdg==211) k_chi_pi=AlgScore.fValue;
                if(AlgScore.fAssumedPdg==321) k_chi_k=AlgScore.fValue;
              }
            } 
            if(AlgScore.fAlgName == "BraggPeakLLH" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                if(anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
                  if(AlgScore.fAssumedPdg==13) k_Brg_F_mu=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==2212) k_Brg_F_p=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==211) k_Brg_F_pi=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==321) k_Brg_F_k=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==0) k_Brg_F_mip=AlgScore.fValue;
                }
                else if(anab::kTrackDir(AlgScore.fTrackDir) == anab::kBackward){
                  if(AlgScore.fAssumedPdg==13) k_Brg_B_mu=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==2212) k_Brg_B_p=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==211) k_Brg_B_pi=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==321) k_Brg_B_k=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==0) k_Brg_B_mip=AlgScore.fValue;
                }
              }
            }
            if(AlgScore.fAlgName == "PIDA_mean" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                k_pida_mean=AlgScore.fValue;
              }
            }
            if(AlgScore.fAlgName == "PIDA_median" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                k_pida_med=AlgScore.fValue;
              }
            }
            if(AlgScore.fAlgName == "PIDA_kde" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                k_kde=AlgScore.fValue;
              }
            }
            double dQdxcalibval = 198.;
            if(AlgScore.fAlgName == "TruncatedMean"){
              if(anab::kVariableType(AlgScore.fVariableType) == anab::kdEdxtruncmean) k_trm_dedx=AlgScore.fValue;
              if(anab::kVariableType(AlgScore.fVariableType) == anab::kdQdxtruncmean) k_trm_dqdx=AlgScore.fValue*dQdxcalibval;
            }
          } // loop over algo score vector...
          k_p_max=TMath::Max(k_Brg_F_p,k_Brg_B_p);
          k_k_max=TMath::Max(k_Brg_F_k,k_Brg_B_k);
          k_pi_max=TMath::Max(k_Brg_F_pi,k_Brg_B_pi);
          k_mu_max=TMath::Max(k_Brg_F_mu,k_Brg_B_mu);
          k_mip_max=TMath::Max(k_Brg_F_mip,k_Brg_B_mip);
          k_L1_ratio=float(k_mip_max)/(k_p_max+k_k_max);
          k_LL1=TMath::Log(k_L1_ratio);
          k_L2_ratio=float(k_mip_max+k_mu_max+k_pi_max)/(k_mip_max+k_mu_max+k_pi_max+k_p_max+k_k_max);
          k_LL2=TMath::Log(k_L2_ratio);
          k_Lp_ratio=float(k_mip_max+k_mu_max)/(k_mip_max+k_mu_max+k_p_max);
          k_LLp=TMath::Log(k_Lp_ratio);
          k_Lk_ratio=float(k_mip_max+k_mu_max)/(k_mip_max+k_mu_max+k_k_max);
          k_LLk=TMath::Log(k_Lk_ratio);
        } // pid vector non-empty
      } // associations are there

    } // valid Kaon ID

    if(mu_can_trkid!=-9999){

      cout << "Find true particle for muon track " << mu_can_trkid << endl;
      simb::MCParticle const* matched_mcparticle = NULL;
      std::unordered_map<int,double> trkide;
      double maxe=-1, tote=0;
      std::vector<simb::MCParticle const*> particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(mu_can_trkid);
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle,evt,fHitTruthAssns);

      for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
        particle_vec.clear(); match_vec.clear();
        particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
        for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
          trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
          tote += match_vec[i_p]->energy;
          if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
            maxe = trkide[ particle_vec[i_p]->TrackId() ];
            matched_mcparticle = particle_vec[i_p];
          }
        }
      }

      std::string dec("Decay");
      if(matched_mcparticle){
        int mu_geant_id=matched_mcparticle->TrackId();
        const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",mu_geant_id);
        mu_origin=int(mc_truth->Origin());
        mu_pdg=matched_mcparticle->PdgCode();
        //cout << "Found truth particle " << mu_pdg << endl;
        if(matched_mcparticle->Process()==dec) mu_isDec=1;
        if(matched_mcparticle->Mother()==k_geant_ID) mu_k_is_Mother=1;

        std::cout << "******************** Kaon track true matched ID : " << k_geant_ID << std::endl;
        std::cout << "******************** Muon track true mother matched ID : " << matched_mcparticle->Mother() << std::endl;

        if(mu_origin==1 && mu_pdg==-13 && mu_isDec==1 && mu_k_is_Mother==1) mu_ness=1;

        if(mu_origin==1 && (mu_pdg==-13 || mu_pdg==211) && mu_isDec==1 && mu_k_is_Mother!=1){
          for(auto const& pPart : ptList){
            int mom_trkid=matched_mcparticle->Mother();
            if(pPart->TrackId()==mom_trkid){ 
              mu_mom_process=pPart->Process();
              std::cout << "********************** Process name of the New Mother : " << mu_mom_process << std::endl;
              std::cout << "********************** ID of the New Mother : " << pPart->Mother() << std::endl;
              if(pPart->Process()=="kaon+Inelastic" && pPart->Mother()==k_geant_ID) mu_mom_k_inelas=1;
              break;
            }
          }
        }

      }

      cout << "Check muon track position" << endl;
      art::Ptr<recob::Track> ptrack(trackListHandle,mu_can_trkid);
      const recob::Track& track = *ptrack;
      TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
      TVector3 end(track.End().X(),track.End().Y(),track.End().Z());
      TVector3 dir_start(track.VertexDirection().X(),track.VertexDirection().Y(),track.VertexDirection().Z());	
      mu_plen=track.Length();
      mu_phi=dir_start.Phi();
      mu_theta=dir_start.Theta();
      if (isInsideVolume("5cmTPC",pos.X(),pos.Y(),pos.Z()) &&
          isInsideVolume("5cmTPC",end.X(),end.Y(),end.Z())) {
        mu_in_5_TPC=1;
      }
      if (isInsideVolume("CCInclusiveTPC",pos.X(),pos.Y(),pos.Z()) &&
          isInsideVolume("CCInclusiveTPC",end.X(),end.Y(),end.Z())) {
        mu_in_CC_TPC=1; 
      }

      cout << "Check muon track calorimetry" << endl;
      float KE_Mu_p2=0;
      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(mu_can_trkid);
      for(unsigned int ical=0; ical<calos.size(); ++ical){
        if(!calos[ical]) continue;
        if(!calos[ical]->PlaneID().isValid) continue;
        int planenum = calos[ical]->PlaneID().Plane;
        if(planenum<0||planenum>2) continue; 
        if(planenum==fplane){
          const size_t NHits=calos[ical] -> dEdx().size();
          mu_hit=int(NHits);
          mu_range=calos[ical]->Range();
          float bin_size=0;
          if(mu_range<10) bin_size=float(mu_range)/2;
          if(mu_range>=10) bin_size=5.0;
          std::vector<float> st_vec,en_vec;
          float rr_x0=-1;float rr_y0=-1;float rr_z0=-1;
          for(size_t iHit = 0; iHit < NHits; ++iHit){
            const auto& TrkPos=(calos[ical] -> XYZ())[iHit];
            if(iHit==0){
              rr_x0=TrkPos.X();rr_y0=TrkPos.Y();rr_z0=TrkPos.Z();
            }
            if((calos[ical]->ResidualRange())[iHit]<=bin_size) en_vec.push_back((calos[ical] -> dEdx())[iHit]);
            if((mu_range-(calos[ical]->ResidualRange())[iHit])<=bin_size) st_vec.push_back((calos[ical] -> dEdx())[iHit]);
            if(iHit>=1){
              if((calos[ical]->ResidualRange())[iHit]>=0 && (calos[ical]->ResidualRange())[iHit-1]>=0)
                KE_Mu_p2=KE_Mu_p2+TMath::Abs(((calos[ical]->ResidualRange())[iHit]-(calos[ical]->ResidualRange())[iHit-1])*(calos[ical] -> dEdx())[iHit]);
            }
            if(iHit==0){
              if((calos[ical]->ResidualRange())[iHit]>=0 && (calos[ical]->ResidualRange())[iHit+1]>=0)
                KE_Mu_p2=KE_Mu_p2+TMath::Abs(((calos[ical]->ResidualRange())[iHit]-(calos[ical]->ResidualRange())[iHit+1])*(calos[ical] -> dEdx())[iHit]);
            }
          } 

          if(st_vec.size() && en_vec.size()){
            float en_median=TMath::Median(en_vec.size(),&en_vec[0]);
            float st_median=TMath::Median(st_vec.size(),&st_vec[0]);
            float st_hit0=TMath::Sqrt((pos.X()-rr_x0)*(pos.X()-rr_x0) + (pos.Y()-rr_y0)*(pos.Y()-rr_y0) + (pos.Z()-rr_z0)*(pos.Z()-rr_z0));
            float en_hit0=TMath::Sqrt((end.X()-rr_x0)*(end.X()-rr_x0) + (end.Y()-rr_y0)*(end.Y()-rr_y0) + (end.Z()-rr_z0)*(end.Z()-rr_z0));

            if(st_hit0<en_hit0){
              mu_start_dedx=en_median;
              mu_end_dedx=st_median;
            }

            if(en_hit0<st_hit0){
              mu_end_dedx=en_median;
              mu_start_dedx=st_median;
            }

            if(en_hit0==st_hit0){
              mu_end_dedx=en_median;
              mu_start_dedx=st_median;
            }

            if(en_median>st_median){ 
              mu_large_dedx=en_median;
              mu_small_dedx=st_median;
            }
            if(en_median<st_median){ 
              mu_large_dedx=st_median;
              mu_small_dedx=en_median;
            }
            if(en_median==st_median){ 
              mu_large_dedx=en_median;
              mu_small_dedx=st_median;
            }
          }
        } // grab plane 2
      } // loop over calorimetry objects

      mu_KE=KE_Mu_p2;

      cout << "Check muon track PID" << endl;
      if(trackPIDAssn.isValid()){
        std::vector<art::Ptr<anab::ParticleID>> trackPID=trackPIDAssn.at(mu_can_trkid);
        if(trackPID.size()!=0){
          std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
          float mu_Brg_F_mu=-9999;float mu_Brg_F_p=-9999;float mu_Brg_F_pi=-9999;float mu_Brg_F_k=-9999;float mu_Brg_F_mip=-9999;
          float mu_Brg_B_mu=-9999;float mu_Brg_B_p=-9999;float mu_Brg_B_pi=-9999;float mu_Brg_B_k=-9999;float mu_Brg_B_mip=-9999;
          for(size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
            anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
            //int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneID);
            //int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
            if(AlgScore.fAlgName == "Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                if(AlgScore.fAssumedPdg==13) mu_chi_mu=AlgScore.fValue;
                if(AlgScore.fAssumedPdg==2212) mu_chi_p=AlgScore.fValue;
                if(AlgScore.fAssumedPdg==211) mu_chi_pi=AlgScore.fValue;
                if(AlgScore.fAssumedPdg==321) mu_chi_k=AlgScore.fValue;
              }
            } 
            if(AlgScore.fAlgName == "BraggPeakLLH" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                if(anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
                  if(AlgScore.fAssumedPdg==13) mu_Brg_F_mu=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==2212) mu_Brg_F_p=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==211) mu_Brg_F_pi=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==321) mu_Brg_F_k=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==0) mu_Brg_F_mip=AlgScore.fValue;
                }
                else if(anab::kTrackDir(AlgScore.fTrackDir) == anab::kBackward){
                  if(AlgScore.fAssumedPdg==13) mu_Brg_B_mu=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==2212) mu_Brg_B_p=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==211) mu_Brg_B_pi=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==321) mu_Brg_B_k=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==0) mu_Brg_B_mip=AlgScore.fValue;
                }
              }
            }
            if(AlgScore.fAlgName == "PIDA_mean" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                mu_pida_mean=AlgScore.fValue;
              }
            }
            if(AlgScore.fAlgName == "PIDA_median" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                mu_pida_med=AlgScore.fValue;
              }
            }
            if(AlgScore.fAlgName == "PIDA_kde" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                mu_kde=AlgScore.fValue;
              }
            }
            double dQdxcalibval = 198.;
            if(AlgScore.fAlgName == "TruncatedMean"){
              if(anab::kVariableType(AlgScore.fVariableType) == anab::kdEdxtruncmean) mu_trm_dedx=AlgScore.fValue;
              if(anab::kVariableType(AlgScore.fVariableType) == anab::kdQdxtruncmean) mu_trm_dqdx=AlgScore.fValue*dQdxcalibval;
            }
          } // loop over algo score vector...
          mu_p_max=TMath::Max(mu_Brg_F_p,mu_Brg_B_p);
          mu_k_max=TMath::Max(mu_Brg_F_k,mu_Brg_B_k);
          mu_pi_max=TMath::Max(mu_Brg_F_pi,mu_Brg_B_pi);
          mu_mu_max=TMath::Max(mu_Brg_F_mu,mu_Brg_B_mu);
          mu_mip_max=TMath::Max(mu_Brg_F_mip,mu_Brg_B_mip);
          mu_L1_ratio=float(mu_mip_max)/(mu_p_max+mu_k_max);
          mu_LL1=TMath::Log(mu_L1_ratio);
          mu_L2_ratio=float(mu_mip_max+mu_mu_max+mu_pi_max)/(mu_mip_max+mu_mu_max+mu_pi_max+mu_p_max+mu_k_max);
          mu_LL2=TMath::Log(mu_L1_ratio);
          mu_Lp_ratio=float(mu_mip_max+mu_mu_max)/(mu_mip_max+mu_mu_max+mu_p_max);
          mu_LLp=TMath::Log(mu_Lp_ratio);
          mu_Lk_ratio=float(mu_mip_max+mu_mu_max)/(mu_mip_max+mu_mu_max+mu_k_max);
          mu_LLk=TMath::Log(mu_Lk_ratio); 
        } // pid vector non-empty
      } // associations are there
    } // valid Muon ID
    //} // kaon and muon ID are valid

    //////////////////////////////////////////////////////////// Looking into CC  Muon

    cout << "Find true particle for CC muon track" << endl;
    for(int k=0; k<NTracks; k++){

      art::Ptr<recob::Track> ptrack(trackListHandle,k);
      //const recob::Track& track = *ptrack;

      //cout << "track " << k << " track id " << track.ID() << " track key " << ptrack.key() << endl;

      simb::MCParticle const* matched_mcparticle = NULL;
      std::unordered_map<int,double> trkide;

      double maxe=-1, tote=0;
      std::vector<simb::MCParticle const*> particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle,evt,fHitTruthAssns);

      for (size_t i_h=0; i_h<hits_from_track.size(); i_h++) {
        particle_vec.clear(); match_vec.clear();
        particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
        for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
          trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
          tote += match_vec[i_p]->energy;
          if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
            maxe = trkide[ particle_vec[i_p]->TrackId() ];
            matched_mcparticle = particle_vec[i_p];
          }
        }
      }

      std::string pri("primary");

      if(matched_mcparticle){
        int mu_geant_id=matched_mcparticle->TrackId();
        const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",mu_geant_id);
        int cc_mu_origin=int(mc_truth->Origin());
        cc_mu_pdg=matched_mcparticle->PdgCode();
        cout << "Found truth particle " << cc_mu_pdg << endl;
        if(matched_mcparticle->Process()==pri && TMath::Abs(cc_mu_pdg)==13 && cc_mu_origin==1){ 
          pri_Mu_is=1;
          cc_mu_trkid=ptrack.key();
          break;
        }
      }

    }

    ////////////////////////////////////////////// End of truth Muon inforamtion  //////////////////////////////////////////////////////////////////////////////

    if(cc_mu_trkid!=-9999){

      cout << "Check CC muon track position" << endl;
      art::Ptr<recob::Track> ptrack(trackListHandle,cc_mu_trkid);
      const recob::Track& track = *ptrack;
      cc_mu_tlen=track.Length();
      TVector3 dir_start(track.VertexDirection().X(),track.VertexDirection().Y(),track.VertexDirection().Z());	
      cc_mu_phi=dir_start.Phi();
      cc_mu_theta=dir_start.Theta();
      TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
      TVector3 end(track.End().X(),track.End().Y(),track.End().Z());
      float st_vtx=TMath::Sqrt((reco_ccmu_vtx_x-pos.X())*(reco_ccmu_vtx_x-pos.X()) +
                               (reco_ccmu_vtx_y-pos.Y())*(reco_ccmu_vtx_y-pos.Y()) +
                               (reco_ccmu_vtx_z-pos.Z())*(reco_ccmu_vtx_z-pos.Z()));
      cc_dis_vtx=st_vtx;

      cout << "Check CC muon track calorimetry" << endl;
      float KE_Mu_p2=0;
      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(cc_mu_trkid);
      for(unsigned int ical=0; ical<calos.size(); ++ical){
        if(!calos[ical]) continue;
        if(!calos[ical]->PlaneID().isValid) continue;
        int planenum = calos[ical]->PlaneID().Plane;
        if(planenum<0||planenum>2) continue; 
        if(planenum==fplane){
          const size_t NHits = calos[ical] -> dEdx().size();
          cc_mu_range=calos[ical]->Range();
          cc_mu_hit=int(NHits);
          float bin_size=0;
          if(cc_mu_range<10) bin_size=float(cc_mu_range)/2;
          if(cc_mu_range>=10) bin_size=5.0;
          std::vector<float> st_vec,en_vec;
          for(size_t iHit = 0; iHit < NHits; ++iHit){
            if((calos[ical]->ResidualRange())[iHit]<=bin_size) en_vec.push_back((calos[ical] -> dEdx())[iHit]);
            if((cc_mu_range-(calos[ical]->ResidualRange())[iHit])<=bin_size) st_vec.push_back((calos[ical] -> dEdx())[iHit]);
            if(iHit>=1){
              if((calos[ical]->ResidualRange())[iHit]>=0 && (calos[ical]->ResidualRange())[iHit-1]>=0)
                KE_Mu_p2=KE_Mu_p2+TMath::Abs(((calos[ical]->ResidualRange())[iHit]-(calos[ical]->ResidualRange())[iHit-1])*(calos[ical] -> dEdx())[iHit]);
            }
            if(iHit==0){
              if((calos[ical]->ResidualRange())[iHit]>=0 && (calos[ical]->ResidualRange())[iHit+1]>=0)
                KE_Mu_p2=KE_Mu_p2+TMath::Abs(((calos[ical]->ResidualRange())[iHit]-(calos[ical]->ResidualRange())[iHit+1])*(calos[ical] -> dEdx())[iHit]);
            }
          }

          if(st_vec.size() && en_vec.size()){
            float en_median=TMath::Median(en_vec.size(),&en_vec[0]);
            float st_median=TMath::Median(st_vec.size(),&st_vec[0]);
            if(en_median>st_median){ 
              cc_mu_large_dedx=en_median;
              cc_mu_small_dedx=st_median;
            }
            if(en_median<st_median){ 
              cc_mu_large_dedx=st_median;
              cc_mu_small_dedx=en_median;
            }
            if(en_median==st_median){ 
              cc_mu_large_dedx=en_median;
              cc_mu_small_dedx=st_median;
            }
          }
        } // grab plane 2
      } // loop over calorimetry objects

      cc_mu_KE=KE_Mu_p2;

      cout << "Check CC muon track PID" << endl;
      if(trackPIDAssn.isValid()){
        std::vector<art::Ptr<anab::ParticleID>> trackPID=trackPIDAssn.at(cc_mu_trkid);
        if(trackPID.size()!=0){
          std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
          float cc_mu_Brg_F_mu=-9999;float cc_mu_Brg_F_p=-9999;float cc_mu_Brg_F_pi=-9999;float cc_mu_Brg_F_k=-9999;float cc_mu_Brg_F_mip=-9999;
          float cc_mu_Brg_B_mu=-9999;float cc_mu_Brg_B_p=-9999;float cc_mu_Brg_B_pi=-9999;float cc_mu_Brg_B_k=-9999;float cc_mu_Brg_B_mip=-9999;
          for(size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
            anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
            //int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneID);
            //int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
            if(AlgScore.fAlgName == "Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                if(AlgScore.fAssumedPdg==13) cc_mu_chi_mu=AlgScore.fValue;
                if(AlgScore.fAssumedPdg==2212) cc_mu_chi_p=AlgScore.fValue;
                if(AlgScore.fAssumedPdg==211) cc_mu_chi_pi=AlgScore.fValue;
                if(AlgScore.fAssumedPdg==321) cc_mu_chi_k=AlgScore.fValue;
              }
            } 
            if(AlgScore.fAlgName == "BraggPeakLLH" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                if(anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
                  if(AlgScore.fAssumedPdg==13) cc_mu_Brg_F_mu=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==2212) cc_mu_Brg_F_p=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==211) cc_mu_Brg_F_pi=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==321) cc_mu_Brg_F_k=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==0) cc_mu_Brg_F_mip=AlgScore.fValue;
                }
                else if(anab::kTrackDir(AlgScore.fTrackDir) == anab::kBackward){
                  if(AlgScore.fAssumedPdg==13) cc_mu_Brg_B_mu=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==2212) cc_mu_Brg_B_p=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==211) cc_mu_Brg_B_pi=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==321) cc_mu_Brg_B_k=AlgScore.fValue;
                  if(AlgScore.fAssumedPdg==0) cc_mu_Brg_B_mip=AlgScore.fValue;
                }
              }
            }
            if(AlgScore.fAlgName == "PIDA_mean" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                cc_mu_pida_mean=AlgScore.fValue;
              }
            }
            if(AlgScore.fAlgName == "PIDA_median" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                cc_mu_pida_med=AlgScore.fValue;
              }
            }
            if(AlgScore.fAlgName == "PIDA_kde" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                cc_mu_kde=AlgScore.fValue;
              }
            }
            double dQdxcalibval = 198.;
            if(AlgScore.fAlgName == "TruncatedMean"){
              if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==fplane) {
                if(anab::kVariableType(AlgScore.fVariableType) == anab::kdEdxtruncmean) cc_mu_trm_dedx=AlgScore.fValue;
                if(anab::kVariableType(AlgScore.fVariableType) == anab::kdQdxtruncmean) cc_mu_trm_dqdx=AlgScore.fValue*dQdxcalibval;
              }
            }
          } // loop over algo score vector...
          cc_mu_p_max=TMath::Max(cc_mu_Brg_F_p,cc_mu_Brg_B_p);
          cc_mu_k_max=TMath::Max(cc_mu_Brg_F_k,cc_mu_Brg_B_k);
          cc_mu_pi_max=TMath::Max(cc_mu_Brg_F_pi,cc_mu_Brg_B_pi);
          cc_mu_mu_max=TMath::Max(cc_mu_Brg_F_mu,cc_mu_Brg_B_mu);
          cc_mu_mip_max=TMath::Max(cc_mu_Brg_F_mip,cc_mu_Brg_B_mip);
          cc_mu_L1_ratio=float(cc_mu_mip_max)/(cc_mu_p_max+cc_mu_k_max);
          cc_mu_LL1=TMath::Log(cc_mu_L1_ratio);
          cc_mu_L2_ratio=float(cc_mu_mip_max+cc_mu_mu_max+cc_mu_pi_max)/(cc_mu_mip_max+cc_mu_mu_max+cc_mu_pi_max+cc_mu_p_max+cc_mu_k_max);
          cc_mu_LL2=TMath::Log(cc_mu_L1_ratio);
          cc_mu_Lp_ratio=float(cc_mu_mip_max+cc_mu_mu_max)/(cc_mu_mip_max+cc_mu_mu_max+cc_mu_p_max);
          cc_mu_LLp=TMath::Log(cc_mu_Lp_ratio);
          cc_mu_Lk_ratio=float(cc_mu_mip_max+cc_mu_mu_max)/(cc_mu_mip_max+cc_mu_mu_max+cc_mu_k_max);
          cc_mu_LLk=TMath::Log(cc_mu_Lk_ratio); 
        } // pid vector non-empty
      } // associations are there
    } // cc mu trkid is valid

  }  // long track clost to vtx and multiplicity >2
  //} // more than 3 track reconstructed in the event

  //////////////////////////////////////////////////////////// Looking into the  truth level information //////////////////////////////////////////////////////////////////

  cout << "Looking into the truth level information" << endl;

  int gk_trk_ID=-9999;
  float k_px=-9999;float k_py=-9999; float k_pz=-9999;
  float kp=-9999;
  for(auto const& pPart : ptList){
    std::string pri("primary");
    bool isPrimary=0;
    isPrimary=pPart->Process()==pri;
    TLorentzVector mcstart, mcend;
    unsigned int pstarti, pendi;
    double plen=length(*pPart, mcstart, mcend, pstarti, pendi);
    if(isPrimary && TMath::Abs(pPart->PdgCode())==13){
      const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",pPart->TrackId());
      int origin=int(mc_truth->Origin());
      if(origin==1){ 
        Tr_pri_mu_is=1; 
        Tr_pri_mu_pdg=pPart->PdgCode();
      } 
    }
    if(isPrimary && pPart->PdgCode()==321){
      const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",pPart->TrackId());
      int origin=int(mc_truth->Origin());
      if(origin==1){ 
        if(pPart->EndE()*1000<520) Tr_pri_st_k_is=1; // initially this was 510
        gk_trk_ID=pPart->TrackId();
        std::cout << "********************** Truth level Kaon trk ID : " << gk_trk_ID <<  std::endl;
        std::cout << "****************************** End Energy of Primary Kaon : " << pPart->EndE()*1000 << std::endl;
        float ksx=pPart->Vx();float ksy=pPart->Vy();float ksz=pPart->Vz();
        float kex=pPart->EndX();float key=pPart->EndY();float kez=pPart->EndZ();
        Tr_k_plen=plen;
        Tr_k_endE=pPart->EndE();
        Tr_k_phi=pPart->Momentum().Phi();
        Tr_k_theta=pPart->Momentum().Theta();
        k_px=pPart->Px(); k_py=pPart->Py(); k_pz=pPart->Pz();
        kp=pPart->Momentum().Vect().Mag();
        if (isInsideVolume("TPC",ksx,ksy,ksz) &&
            isInsideVolume("TPC",kex,key,kez)) {
          Tr_k_inTPC=1;
        }
        if (isInsideVolume("5cmTPC",ksx,ksy,ksz) &&
            isInsideVolume("5cmTPC",kex,key,kez)) {
          Tr_k_in_5_TPC=1;
        }
        if (isInsideVolume("CCInclusiveTPC",ksx,ksy,ksz) &&
            isInsideVolume("CCInclusiveTPC",kex,key,kez)) {
          Tr_k_in_CC_TPC=1;
        }
      }  
    }
  }// loop over plist

  for (auto const& sPart : ptList) {
    TLorentzVector mcstart, mcend;
    unsigned int pstarti, pendi;
    double plen=length(*sPart, mcstart, mcend, pstarti, pendi);
    if(sPart->Process()=="Decay" && (sPart->PdgCode()==-13 || sPart->PdgCode()==211)){ 
      std::cout << "********** Truth info Found a decaying " << sPart->PdgCode() << std::endl;
      std::cout << "********** Mother trkID of the Found decay particle " << sPart->Mother() << std::endl;
    }
    if(sPart->Process()=="Decay" && (sPart->PdgCode()==-13 || sPart->PdgCode()==211) && sPart->Mother()==gk_trk_ID){
      std::cout << "****************************** We found typical kaon decaying event ********************************" << std::endl;
      const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",sPart->TrackId());
      int origin=int(mc_truth->Origin());
      if(origin==1){
        Tr_dec_mu_is=1;
        Tr_dec_mu_pi_pdg=sPart->PdgCode();
        float musx=sPart->Vx();float musy=sPart->Vy();float musz=sPart->Vz();
        float muex=sPart->EndX();float muey=sPart->EndY();float muez=sPart->EndZ();
        Tr_mu_plen=plen;
        Tr_mu_phi=sPart->Momentum().Phi();
        Tr_mu_theta=sPart->Momentum().Theta();
        Tr_kmu_open_ang=(TMath::ACos(float(sPart->Px()*k_px + sPart->Py()*k_py + sPart->Pz()*k_pz)/(kp*sPart->Momentum().Vect().Mag())));
        if (isInsideVolume("TPC",musx,musy,musz) &&
            isInsideVolume("TPC",muex,muey,muez)) {
          Tr_mu_inTPC=1;
        }
        if (isInsideVolume("5cmTPC",musx,musy,musz) &&
            isInsideVolume("5cmTPC",muex,muey,muez)) {
          Tr_mu_in_5_TPC=1;
        }
        if (isInsideVolume("CCInclusiveTPC",musx,musy,musz) &&
            isInsideVolume("CCInclusiveTPC",muex,muey,muez)) {
          Tr_mu_in_CC_TPC=1;
        }
        break;
      }
    }

    if(sPart->Process()=="Decay" && (sPart->PdgCode()==-13 || sPart->PdgCode()==211) && sPart->Mother()!=gk_trk_ID){
      for(auto const& pPart : ptList){
        int mom_trkid=sPart->Mother();
        if(pPart->TrackId()==mom_trkid){ 
          std::cout << "************** Process Name : " << pPart->Process() << std::endl;
          std::cout << "************** Mother ID of Mother : " << pPart->Mother() << std::endl;	  
          if(pPart->Process()=="kaon+Inelastic" && pPart->Mother()==gk_trk_ID){ // initiall (pPart->Process()=="kaon+Inel")

            ///////////////////////////////////////// checking whetehr this inelastic kaon is making it's own track ////////////////

            for(int i=0; i<NTracks; i++){
              art::Ptr<recob::Track> ptrack(trackListHandle,i);
              const recob::Track& track = *ptrack;

              simb::MCParticle const* matched_mcparticle = NULL;
              std::unordered_map<int,double> trkide;
              double maxe=-1, tote=0;
              std::vector<simb::MCParticle const*> particle_vec;
              std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
              //std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(track.ID());
              std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
              art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle,evt,fHitTruthAssns);

              for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
                particle_vec.clear(); match_vec.clear();
                particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
                for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                  trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
                  tote += match_vec[i_p]->energy;
                  if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
                    maxe = trkide[ particle_vec[i_p]->TrackId() ];
                    matched_mcparticle = particle_vec[i_p];
                  }
                }
              }

              if(matched_mcparticle){
                if(matched_mcparticle->TrackId()==pPart->TrackId()){
                  kinelas_has_traks=1;
                  //kinelas_reco_trkID=track.ID();
                  kinelas_reco_trkID=ptrack.key();
                  kinelas_tlen=track.Length();
                  break;
                }
              }

            }     

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            std::cout << "****************************** We found non-typical kaon decaying event ********************************" << std::endl;
            std::cout << "****************************** End Energy of Secondary Kaon : " << pPart->EndE()*1000 << std::endl;

            True_kinelas_KE=(pPart->E()-0.493);
            TLorentzVector mcstart, mcend;
            unsigned int pstarti, pendi;
            True_kinelas_tlen=length(*pPart, mcstart, mcend, pstarti, pendi);
            Tr_K_Inelas=1;    

            std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Testing Kaon Inelastic interaction &&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
            std::cout << "Inelastic Kaon has his own reco track : " << kinelas_has_traks << std::endl;
            std::cout << "Reco track ID of Inelastic Kaon and Primary Kaon : " << kinelas_reco_trkID << "  " << k_can_trkid << std::endl;
            std::cout << "True KE of Inelastic Kaon : " << True_kinelas_KE << std::endl;
            std::cout << "True tracklength of Inelastic Kaon : " << True_kinelas_tlen << std::endl;
            std::cout << "Kaon Inelastic swith is " << Tr_K_Inelas << std::endl;

            const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",sPart->TrackId());
            int origin=int(mc_truth->Origin());
            if(origin==1){
              Tr_dec_mu_is=1;
              Tr_dec_mu_pi_pdg=sPart->PdgCode();
              float musx=sPart->Vx();float musy=sPart->Vy();float musz=sPart->Vz();
              float muex=sPart->EndX();float muey=sPart->EndY();float muez=sPart->EndZ();
              Tr_mu_plen=plen;
              Tr_mu_phi=sPart->Momentum().Phi();
              Tr_mu_theta=sPart->Momentum().Theta();
              Tr_kmu_open_ang=(TMath::ACos(float(sPart->Px()*k_px + sPart->Py()*k_py + sPart->Pz()*k_pz)/(kp*sPart->Momentum().Vect().Mag())));
              if (isInsideVolume("TPC",musx,musy,musz) &&
                  isInsideVolume("TPC",muex,muey,muez)) {
                Tr_mu_inTPC=1;
              }
              if (isInsideVolume("5cmTPC",musx,musy,musz) &&
                  isInsideVolume("5cmTPC",muex,muey,muez)) {
                Tr_mu_in_5_TPC=1;
              }
              if (isInsideVolume("CCInclusiveTPC",musx,musy,musz) &&
                  isInsideVolume("CCInclusiveTPC",muex,muey,muez)) {
                Tr_mu_in_CC_TPC=1;
              }
              break;
            }
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////// Getting Calorimetric information of Kaon and Muon //////////////////////////////////////////////////


  if(k_can_trkid!=-9999 && mu_can_trkid!=-9999){
    cout << "Getting Calorimetric information of Kaon and Muon" << endl;
    std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(k_can_trkid);
    for(size_t ical = 0; ical<calos.size(); ++ical){
      if(!calos[ical]) continue;
      if(!calos[ical]->PlaneID().isValid) continue;
      int planenum = calos[ical]->PlaneID().Plane;
      if(planenum<0||planenum>2) continue;
      const size_t NHits = calos[ical] -> dEdx().size();
      for(size_t iHit = 0; iHit < NHits; ++iHit){
        k_dedx[planenum][iHit]=(calos[ical] -> dEdx())[iHit];
        k_rr[planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
      }
    }
    calos=fmcal.at(mu_can_trkid);
    for(size_t ical = 0; ical<calos.size(); ++ical){
      if(!calos[ical]) continue;
      if(!calos[ical]->PlaneID().isValid) continue;
      int planenum = calos[ical]->PlaneID().Plane;
      if(planenum<0||planenum>2) continue;
      const size_t NHits = calos[ical] -> dEdx().size();
      for(size_t iHit = 0; iHit < NHits; ++iHit){
        mu_dedx[planenum][iHit]=(calos[ical] -> dEdx())[iHit];
        mu_rr[planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
      }
    }
  }

  cout << "Filling tree" << endl;
  fEventTree->Fill();

} // end of analyze function
 
 /////////////////////////////////////////// Reset Function ///////////////////////////////
 
 void CCKaonFilter::reset(){

      run=-9999;
      subrun=-9999;
      event=-9999;

      true_nu_energy=-9999;
      true_nu_pdg=-9999;
      true_nu_mode=-9999;
      true_nu_ccnc=-9999;
      true_nu_vtx_x=-9999;
      true_nu_vtx_y=-9999;
      true_nu_vtx_z=-9999;
      true_nu_vtx_inTPC=false;
      true_nu_vtx_in5cmTPC=false;
      true_nu_vtx_inCCInclusiveTPC=false;
      true_nu_vtx_inOldCCInclusiveTPC=false;

      true_lepton_pdg=-9999;
      true_lepton_p=-9999;
      true_lepton_theta=-9999;
      true_lepton_phi=-9999;

      true_kaon_length=-9999;
      true_kaon_ke=-9999;
      true_kaon_theta=-9999;
      true_kaon_phi=-9999;
      true_kaon_lepton_angle=-9999;
      true_kaon_end_process=-9999;
      true_kaon_end_ke=-9999;
      true_kaon_end_x=-9999;
      true_kaon_end_y=-9999;
      true_kaon_end_z=-9999;
      true_kaon_end_inTPC=false;
      true_kaon_end_in5cmTPC=false;
      true_kaon_end_inCCInclusiveTPC=false;
      true_kaon_end_inOldCCInclusiveTPC=false;

      true_kaon_daughter_length=-9999;
      true_kaon_daughter_ke=-9999;
      true_kaon_daughter_angle=-9999;
      true_kaon_daughter_pdg=-9999;
      true_kaon_daughter_end_x=-9999;
      true_kaon_daughter_end_y=-9999;
      true_kaon_daughter_end_z=-9999;
      true_kaon_daughter_end_inTPC=false;
      true_kaon_daughter_end_in5cmTPC=false;
      true_kaon_daughter_end_inCCInclusiveTPC=false;
      true_kaon_daughter_end_inOldCCInclusiveTPC=false;

      reco_nu_cc_filter=false;

      reco_nu_vtx_x=-9999;
      reco_nu_vtx_y=-9999;
      reco_nu_vtx_z=-9999;
      reco_nu_vtx_inTPC=false;
      reco_nu_vtx_in5cmTPC=false;
      reco_nu_vtx_inCCInclusiveTPC=false;
      reco_nu_ndaughters=-9999;

      reco_ccmu_vtx_x=-9999;
      reco_ccmu_vtx_y=-9999;
      reco_ccmu_vtx_z=-9999;
      reco_ccmu_vtx_inTPC=false;
      reco_ccmu_vtx_in5cmTPC=false;
      reco_ccmu_vtx_inCCInclusiveTPC=false;

      reco_ntracks=0;
      for(int k=0; k<kMaxTracks; k++){
        reco_track_distance[k]=-9999;
        reco_track_nhits0[k]=-9999;
        reco_track_nhits1[k]=-9999;
        reco_track_nhits2[k]=-9999;
        reco_track_length[k]=-9999;
        reco_track_dir[k]=false;
        reco_track_chi2k[k]=-9999;
        reco_track_chi2p[k]=-9999;
        reco_track_chi2pi[k]=-9999;
        reco_track_chi2mu[k]=-9999;
        reco_track_vtx_inTPC[k]=false;
        reco_track_vtx_in5cmTPC[k]=false;
        reco_track_vtx_inCCInclusiveTPC[k]=false;
        reco_track_end_inTPC[k]=false;
        reco_track_end_in5cmTPC[k]=false;
        reco_track_end_inCCInclusiveTPC[k]=false;
        reco_track_pdg[k]=-9999;
        reco_track_origin[k]=-9999;
        reco_track_primary[k]=false;
        reco_track_ndaughters[k]=0;
        for(int m=0; m<kMaxTracks; m++){
          reco_track_daughter_distance[k][m]=-9999;
          reco_track_daughter_vtx_distance[k][m]=-9999;
          reco_track_daughter_nhits0[k][m]=-9999;
          reco_track_daughter_nhits1[k][m]=-9999;
          reco_track_daughter_nhits2[k][m]=-9999;
          reco_track_daughter_length[k][m]=-9999;
          reco_track_daughter_chi2k[k][m]=-9999;
          reco_track_daughter_chi2p[k][m]=-9999;
          reco_track_daughter_chi2pi[k][m]=-9999;
          reco_track_daughter_chi2mu[k][m]=-9999;
          reco_track_daughter_vtx_inTPC[k][m]=false;
          reco_track_daughter_vtx_in5cmTPC[k][m]=false;
          reco_track_daughter_vtx_inCCInclusiveTPC[k][m]=false;
          reco_track_daughter_end_inTPC[k][m]=false;
          reco_track_daughter_end_in5cmTPC[k][m]=false;
          reco_track_daughter_end_inCCInclusiveTPC[k][m]=false;
          reco_track_daughter_pdg[k][m]=-9999;
          reco_track_daughter_origin[k][m]=-9999;
          reco_track_daughter_primary[k][m]=false;
        }
      }
      reco_track_ndaughters2=-9999;

      k_can_trkid=-9999;
      mu_can_trkid=-9999;
      k_mu_can_dis=-9999;
      k_mu_open_angle=-9999;
      k_vtx_dis=-9999;
      //mtc_k_5cm=-9999;
      //mtc_k_10cm=-9999;
      //mtc_k_endE=-9999;
      //mtc_mu_5cm=-9999;
      //mtc_mu_10cm=-9999;	
      //mtc_mu_pid=-9999;
      k_geant_ID=-9999;
      k_origin=-9999;
      k_pdg=-9999;
      k_isPri=-9999;
      k_endE=-9999;
      k_ness=-9999;
      kaon_vtx_dis=-9999;
      k_plen=-9999;
      k_phi=-9999;
      k_theta=-9999;
      k_in_5_TPC=-9999;
      k_in_CC_TPC=-9999;
      k_hit=-9999;
      k_range=-9999;
      k_KE=-9999;
      k_large_dedx=-9999;
      k_small_dedx=-9999;
      k_chi_p=-9999;
      k_chi_k=-9999;
      k_chi_pi=-9999;
      k_chi_mu=-9999;
      k_p_max=-9999;
      k_k_max=-9999;
      k_pi_max=-9999;
      k_mu_max=-9999;
      k_mip_max=-9999;
      k_L1_ratio=-9999;
      k_LL1=-9999;
      k_L2_ratio=-9999;
      k_LL2=-9999;
      k_Lp_ratio=-9999;
      k_LLp=-9999;
      k_Lk_ratio=-9999;
      k_LLk=-9999;
      k_pida_mean=-9999;
      k_pida_med=-9999;
      k_kde=-9999;
      k_trm_dqdx=-9999;
      k_trm_dedx=-9999;
      mu_pdg=-9999;
      mu_isDec=-9999;
      mu_origin=-9999;
      mu_k_is_Mother=-9999;
      mu_mom_k_inelas=-9999;
      mu_ness=-9999;
      mu_plen=-9999;
      mu_phi=-9999;
      mu_theta=-9999;
      mu_in_5_TPC=-9999;
      mu_in_CC_TPC=-9999;
      mu_KE=-9999;
      mu_hit=-9999;
      mu_range=-9999;
      mu_large_dedx=-9999;
      mu_small_dedx=-9999;
      mu_chi_p=-9999;
      mu_chi_k=-9999;
      mu_chi_pi=-9999;
      mu_chi_mu=-9999;
      mu_p_max=-9999;
      mu_k_max=-9999;
      mu_pi_max=-9999;
      mu_mu_max=-9999;
      mu_mip_max=-9999;
      mu_L1_ratio=-9999;
      mu_LL1=-9999;
      mu_L2_ratio=-9999;
      mu_LL2=-9999;
      mu_Lp_ratio=-9999;
      mu_LLp=-9999;
      mu_Lk_ratio=-9999;
      mu_LLk=-9999;
      mu_pida_mean=-9999;
      mu_pida_med=-9999;
      mu_kde=-9999;
      mu_trm_dqdx=-9999;
      mu_trm_dedx=-9999;
      mu_mom_process="NA";
      cc_mu_trkid=-9999;
      cc_mu_tlen=-9999;
      cc_mu_phi=-9999;
      cc_mu_pdg=-9999;
      cc_mu_theta=-9999;
      cc_mu_range=-9999;
      cc_mu_KE=-9999;
      cc_mu_hit=-9999;
      cc_mu_large_dedx=-9999;
      cc_mu_small_dedx=-9999;
      cc_dis_vtx=-9999;
      cc_mu_chi_p=-9999;
      cc_mu_chi_k=-9999;
      cc_mu_chi_pi=-9999;
      cc_mu_chi_mu=-9999;
      cc_mu_p_max=-9999;
      cc_mu_k_max=-9999;
      cc_mu_pi_max=-9999;
      cc_mu_mu_max=-9999;
      cc_mu_mip_max=-9999;
      cc_mu_L1_ratio=-9999;
      cc_mu_LL1=-9999;
      cc_mu_L2_ratio=-9999;
      cc_mu_LL2=-9999;
      cc_mu_Lp_ratio=-9999;
      cc_mu_LLp=-9999;
      cc_mu_Lk_ratio=-9999;
      cc_mu_LLk=-9999;
      cc_mu_pida_mean=-9999;
      cc_mu_pida_med=-9999;
      cc_mu_kde=-9999;
      cc_mu_trm_dqdx=-9999;
      cc_mu_trm_dedx=-9999;
      longest_trkid=-9999;
      longest_trklen=-9999;
      pri_Mu_is=-9999;
      Tr_pri_mu_pdg=-9999;
      Tr_pri_mu_is=-9999;
      Tr_pri_st_k_is=-9999;
      Tr_K_Inelas=-9999;
      Tr_k_plen=-9999;
      Tr_k_endE=-9999;
      Tr_k_theta=-9999;
      Tr_k_phi=-9999;
      Tr_dec_mu_is=-9999;
      Tr_dec_mu_pi_pdg=-9999;
      Tr_mu_plen=-9999;
      Tr_mu_theta=-9999;
      Tr_mu_phi=-9999;
      Tr_k_inTPC=-9999;
      Tr_mu_inTPC=-9999;
      Tr_k_in_5_TPC=-9999;
      Tr_k_in_CC_TPC=-9999;
      Tr_mu_in_5_TPC=-9999;
      Tr_mu_in_CC_TPC=-9999;
      Tr_kmu_open_ang=-9999;
      vtx_5cm_mult=-9999;
      k_start_dedx=-9999;
      k_end_dedx=-9999;
      mu_start_dedx=-9999;
      mu_end_dedx=-9999;
      cut_1=-9999;
      cut_2=-9999;
      cut_3=-9999;
      cut_4=-9999;
      cut_5=-9999;
      cut_6=-9999;
      cut_7=-9999;
      cut_8=-9999;
      cut_9=-9999;
      cut_10=-9999;
      cut_11=-9999;
      cut_12=-9999;
      kinelas_has_traks=-9999;
      kinelas_reco_trkID=-9999;
      kinelas_tlen=-9999;
      True_kinelas_KE=-9999;
      True_kinelas_tlen=-9999;
      
      for(int j=0; j<3; j++){
	  for(int k=0; k<3000; k++){
	      k_dedx[j][k]=-9999;
              k_rr[j][k]=-9999;
	      mu_dedx[j][k]=-9999;
              mu_rr[j][k]=-9999;
          }
      }
 }
 
 /////////////////////////////////////////////////////////////////////////////////////////
 
bool CCKaonFilter::isInsideVolume(string volume, double x, double y, double z)
{
  if (volume=="TPC") {
    if (x>0 && x<256.35 && y>-116.5 && y<116.5 && z>0 && z<1036.8) {
      return true;
    }
  }
  else if (volume=="5cmTPC") {
    if (x>5 && x<251.35 && y>-111.5 && y<111.5 && z>5 && z<1031.8) {
      return true;
    }
  }
  else if (volume=="CCInclusiveTPC") {
    if (x>10 && x<246.35 && y>-106.5 && y<106.5 && z>10 && z<986.8) {
      return true;
    }
  }
  else if (volume=="OldCCInclusiveTPC") {
    if (x>12 && x<244.35 && y>-81.5 && y<81.5 && z>25 && z<951.8 && !(z>675 && z<775)) {
      return true;
    }
  }
  return false;
}

 double CCKaonFilter::length(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
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
 
 DEFINE_ART_MODULE(CCKaonFilter)
}


