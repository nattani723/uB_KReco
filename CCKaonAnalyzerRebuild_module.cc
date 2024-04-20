#include "CCKaonAnalyzerRebuild_module.h"

//#include "ubana/SinglePhotonAnalysis/SinglePhoton_module.h"
#include "headers/analyze_Slice.h"
#include "headers/analyze_Tracks.h"
#include "headers/analyze_Showers.h"
#include "headers/analyze_MCTruth.h"
#include "headers/analyze_OpFlashes.h"

//#include "headers/particle_split_basetool_14Aug.h"
#include "headers/particle_split_basetool.h"
#include "headers/track_production.h"


//#include "gallery/Event.h"
//#include "gallery/ValidHandle.h"
//#include "headers/analyze_Template.h"
//#include "headers/second_shower_search.h"
//#include "headers/analyze_EventWeight.h"

/*
//#include "ubana/SinglePhotonAnalysis/fiducial_volume.h"
//#include "ubana/SinglePhotonAnalysis/isolation.h"
*/

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/TriggerResults.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "PID/LLR_PID.h"
#include "PID/LLRPID_proton_muon_lookup.h"

#include "PID_K/LLR_PID_K.h"
#include "PID_K/LLRPID_kaon_proton_lookup.h"

//#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "TTree.h"
#include "TMath.h"

#include <array>
#include <vector>
#include <map>
//#include "LinkDef.h"


#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/make_tool.h"


#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"
#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
    

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "canvas/Persistency/Common/TriggerResults.h" 
#include "fhiclcpp/ParameterSetRegistry.h" 

#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>

#include "TTree.h"
#include "TTimeStamp.h"

//#include "ubana/SinglePhotonAnalysis/SinglePhoton_module.h"

//#ifdef __MAKECINT__
#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+; 
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+; 
#endif

//const int kMaxTracks=20;

using namespace std;
namespace Kaon_Analyzer{

  //namespace microboone{
//namespace single_photon{
	
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
/*
class CCKaonAnalyzer : public art::EDAnalyzer {
  public:

    explicit CCKaonAnalyzer(fhicl::ParameterSet const& pset);
    virtual ~CCKaonAnalyzer();

    void endSubRun(const art::SubRun &subrun);
    void beginJob();
    void analyze(const art::Event& evt);
    void reset();

    //void fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, const recob::Track trk, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart, int track_i=-1, int daughter_i=-1);
    void fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, int track_i=-1, int daughter_i=-1);
    //double ModBoxCorrection(const double dQdx, const float x, const float y, const float z);
    //float GetLocalEFieldMag(const float x, const float y, const float z);
    void fillPID(const std::vector<art::Ptr<anab::ParticleID>> &trackPID, double angle_y, int track_i=-1, int daughter_i=-1);
    void fillTrueMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_track,
                          art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                          int track_i=-1,
                          int daughter_i=-1);
  //void GetPFParticleIdMap(const Kaon_Analyzer::CCKaonAnalyzer::PFParticleHandle &pfParticleHandle, Kaon_Analyzer::CCKaonAnalyzer::PFParticleIdMap &pfParticleMap);
  
    double length(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);

    bool isInsideVolume(string volume, double x, double y, double z);
    bool isInsideVolume(string volume, const TVector3& v) {
      return isInsideVolume(volume, v.X(), v.Y(), v.Z());
    }

  //void 

  private:

    searchingfornues::LLRPID llr_pid_calculator;
    searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;

    searchingfornuesk::LLRPIDK llr_pid_calculator_k;
    searchingfornuesk::KaonProtonLookUpParameters protonmuon_parameters_k;

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

    Int_t   true_lepton_pdg;
    Float_t true_lepton_p;
    Float_t true_lepton_ke;
    Float_t true_lepton_theta;
    Float_t true_lepton_costheta;
    Float_t true_lepton_phi;

    Int_t   true_nkaons;
    Float_t true_kaon_length;
    Float_t true_kaon_p;
    Float_t true_kaon_ke;
    Float_t true_kaon_theta;
    Float_t true_kaon_costheta;
    Float_t true_kaon_phi;
    Float_t true_kaon_ccmuon_angle;
    Float_t true_kaon_ccmuon_cosangle;
    Int_t   true_kaon_end_process;
    Float_t true_kaon_end_ke;
    Float_t true_kaon_end_x;
    Float_t true_kaon_end_y;
    Float_t true_kaon_end_z;
    Bool_t  true_kaon_end_inTPC;
    Bool_t  true_kaon_end_in5cmTPC;
    Bool_t  true_kaon_end_inCCInclusiveTPC;

    Int_t   true_kaon_ndaughters;
    Int_t   true_kaon_ndaughters_decay;
    Int_t   true_kaon_ndaughters_inelastic;
    Int_t   true_kaon_ndecmup;
    Int_t   true_kaon_ndecpip;
    Int_t   true_kaon_ninekap;
    Int_t   true_kaon_ninepip;
    Int_t   true_kaon_ninepro;
    Float_t true_kaon_daughter_length;
    Float_t true_kaon_daughter_p;
    Float_t true_kaon_daughter_ke;
    Float_t true_kaon_daughter_theta;
    Float_t true_kaon_daughter_costheta;
    Float_t true_kaon_daughter_angle;
    Float_t true_kaon_daughter_cosangle;
    Int_t   true_kaon_daughter_pdg;
    Float_t true_kaon_daughter_end_x;
    Float_t true_kaon_daughter_end_y;
    Float_t true_kaon_daughter_end_z;
    Bool_t  true_kaon_daughter_end_inTPC;
    Bool_t  true_kaon_daughter_end_in5cmTPC;
    Bool_t  true_kaon_daughter_end_inCCInclusiveTPC;

    Bool_t  reco_nu_cc_filter;

    Float_t reco_nu_vtx_x;
    Float_t reco_nu_vtx_y;
    Float_t reco_nu_vtx_z;
    Bool_t  reco_nu_vtx_inTPC;
    Bool_t  reco_nu_vtx_in5cmTPC;
    Bool_t  reco_nu_vtx_inCCInclusiveTPC;
    Int_t   reco_nu_ndaughters;
    Int_t   reco_nu_cc_nmue;

    Float_t reco_ccmu_vtx_x;
    Float_t reco_ccmu_vtx_y;
    Float_t reco_ccmu_vtx_z;
    Bool_t  reco_ccmu_vtx_inTPC;
    Bool_t  reco_ccmu_vtx_in5cmTPC;
    Bool_t  reco_ccmu_vtx_inCCInclusiveTPC;
    Int_t   reco_ccmu_true_pdg;
    Int_t   reco_ccmu_true_origin;
    Bool_t  reco_ccmu_true_primary;
    Bool_t  reco_ccmu_true_end_inTPC;
    Bool_t  reco_ccmu_true_end_in5cmTPC;
    Bool_t  reco_ccmu_true_end_inCCInclusiveTPC;
    Float_t reco_ccmu_true_length;

    Int_t   reco_ntracks;
    Int_t   m_reco_track_sliceId[kMaxTracks];
    Int_t   m_reco_track_is_nuslice[kMaxTracks];
    Float_t reco_track_distance[kMaxTracks];
    Int_t   reco_track_nhits0[kMaxTracks];
    Int_t   reco_track_nhits1[kMaxTracks];
    Int_t   reco_track_nhits2[kMaxTracks];

    Float_t   reco_track_kin0[kMaxTracks];
    Float_t   reco_track_kin1[kMaxTracks];
    Float_t   reco_track_kin2[kMaxTracks];
    //vector<vector<Float_t>> reco_track_dEdx;
    //vector<vector<Float_t>> reco_track_ResRan;

    //Float_t reco_track_dEdx_pl0[kMaxTracks][2000];
    //Float_t reco_track_ResRan_pl0[kMaxTracks][2000];

    //Float_t reco_track_dEdx_pl1[kMaxTracks][2000];
    //Float_t reco_track_ResRan_pl1[kMaxTracks][2000];

    //Float_t reco_track_dEdx_pl2[kMaxTracks][2000];
    //Float_t reco_track_ResRan_pl2[kMaxTracks][2000];
    //array<vector<Float_t>,kMaxTracks> reco_track_ResRan_arr;
    //array<vector<Float_t>,kMaxTracks> reco_track_dEdx_arr;

    Float_t reco_track_length[kMaxTracks];
    Float_t reco_track_theta[kMaxTracks];
    Float_t reco_track_phi[kMaxTracks];
    Bool_t  reco_track_dir[kMaxTracks];

    Float_t reco_track_P_vtx[kMaxTracks];
    Float_t reco_track_P_str[kMaxTracks];
    Float_t reco_track_P_end[kMaxTracks];

    Float_t reco_track_chi2ka_pl0[kMaxTracks];
    Float_t reco_track_chi2pr_pl0[kMaxTracks];
    Float_t reco_track_chi2pi_pl0[kMaxTracks];
    Float_t reco_track_chi2mu_pl0[kMaxTracks];
    Float_t reco_track_chi2ka_pl1[kMaxTracks];
    Float_t reco_track_chi2pr_pl1[kMaxTracks];
    Float_t reco_track_chi2pi_pl1[kMaxTracks];
    Float_t reco_track_chi2mu_pl1[kMaxTracks];
    Float_t reco_track_chi2ka_pl2[kMaxTracks];
    Float_t reco_track_chi2pr_pl2[kMaxTracks];
    Float_t reco_track_chi2pi_pl2[kMaxTracks];
    Float_t reco_track_chi2mu_pl2[kMaxTracks];
    Float_t reco_track_chi2ka_3pl[kMaxTracks];
    Float_t reco_track_chi2pr_3pl[kMaxTracks];
    Float_t reco_track_chi2pi_3pl[kMaxTracks];
    Float_t reco_track_chi2mu_3pl[kMaxTracks];
    Float_t reco_track_likepr_3pl[kMaxTracks];


    Float_t reco_track_Bragg_fwd_ka_pl0[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pr_pl0[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pi_pl0[kMaxTracks];
    Float_t reco_track_Bragg_fwd_mu_pl0[kMaxTracks];
    Float_t reco_track_Bragg_fwd_ka_pl1[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pr_pl1[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pi_pl1[kMaxTracks];
    Float_t reco_track_Bragg_fwd_mu_pl1[kMaxTracks];
    Float_t reco_track_Bragg_fwd_ka_pl2[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pr_pl2[kMaxTracks];
    Float_t reco_track_Bragg_fwd_pi_pl2[kMaxTracks];
    Float_t reco_track_Bragg_fwd_mu_pl2[kMaxTracks];

    Float_t reco_track_MIP_pl0[kMaxTracks];
    Float_t reco_track_MIP_pl1[kMaxTracks];
    Float_t reco_track_MIP_pl2[kMaxTracks];


    Float_t reco_track_daughter_Bragg_fwd_ka_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pr_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pi_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_mu_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_ka_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pr_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pi_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_mu_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_ka_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pr_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_pi_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_Bragg_fwd_mu_pl2[kMaxTracks][kMaxTracks];

    Float_t reco_track_daughter_MIP_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_MIP_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_MIP_pl2[kMaxTracks][kMaxTracks];


    Float_t reco_track_llrpid_3pl[kMaxTracks];
    Float_t reco_track_total_llrpid_3pl[kMaxTracks];
    Float_t reco_track_llrpid_k_3pl[kMaxTracks];
    Bool_t  reco_track_vtx_inTPC[kMaxTracks];
    Bool_t  reco_track_vtx_in5cmTPC[kMaxTracks];
    Bool_t  reco_track_vtx_inCCInclusiveTPC[kMaxTracks];
    Bool_t  reco_track_end_inTPC[kMaxTracks];
    Bool_t  reco_track_end_in5cmTPC[kMaxTracks];
    Bool_t  reco_track_end_inCCInclusiveTPC[kMaxTracks];
    Int_t   reco_track_true_pdg[kMaxTracks];
    Int_t   reco_track_true_origin[kMaxTracks];
    Bool_t  reco_track_true_primary[kMaxTracks];
    Bool_t  reco_track_true_end_inTPC[kMaxTracks];
    Bool_t  reco_track_true_end_in5cmTPC[kMaxTracks];
    Bool_t  reco_track_true_end_inCCInclusiveTPC[kMaxTracks];
    Float_t reco_track_true_length[kMaxTracks];

    Int_t   reco_track_ndaughters[kMaxTracks];
    Float_t reco_track_daughter_distance[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_vtx_distance[kMaxTracks][kMaxTracks];
    Float_t reco_angle_track_daughter[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits0[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits1[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits2[kMaxTracks][kMaxTracks];


    //vector<vector<vector<Float_t>>> reco_track_daughter_dEdx;
    //vector<vector<vector<Float_t>>> reco_track_daughter_ResRan;

    //Float_t reco_track_daughter_dEdx_pl0[kMaxTracks][kMaxTracks][2000];
    //Float_t reco_track_daughter_ResRan_pl0[kMaxTracks][kMaxTracks][2000];

    //Float_t reco_track_daughter_dEdx_pl1[kMaxTracks][kMaxTracks][2000];
    //Float_t reco_track_daughter_ResRan_pl1[kMaxTracks][kMaxTracks][2000];

    //Float_t reco_track_daughter_dEdx_pl2[kMaxTracks][kMaxTracks][2000];
    //Float_t reco_track_daughter_ResRan_pl2[kMaxTracks][kMaxTracks][2000];



    Float_t reco_track_daughter_length[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_theta[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_phi[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_pl0[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_pl1[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_pl2[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2ka_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pr_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2pi_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_chi2mu_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_likepr_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_llrpid_3pl[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_llrpid_k_3pl[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_inTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_in5cmTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_vtx_inCCInclusiveTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_inTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_in5cmTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_end_inCCInclusiveTPC[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_true_pdg[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_true_origin[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_primary[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_end_inTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_end_in5cmTPC[kMaxTracks][kMaxTracks];
    Bool_t  reco_track_daughter_true_end_inCCInclusiveTPC[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_true_length[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_true_mother[kMaxTracks][kMaxTracks];

    Int_t   k_can_trkid; //track id of reco kaon
    Int_t   mu_can_trkid; //track id of reco muon
    Float_t k_mu_can_dis; //distance between kaon and muon end
    Float_t k_mu_open_angle; //angle between kaon and muon end
    Float_t k_vtx_dis; //distance between kaon and muon start

    //Int_t k_geant_ID; //kaon track matched true particle geant ID
    //Int_t k_origin; //kaon track matched true particle origin
    //Int_t k_pdg; //kaon track matched true particle pdg
    //Int_t k_isPri; //kaon track matched true particle is primary
    //Float_t k_endE; //kaon track matched true particle end energy
    //Int_t k_ness; //kaon track matched true particle is (pdg==321,primary,endE<510)

    //Float_t kaon_vtx_dis; //kaon track and vertex distance (same as k_vtx_dis?)
    //Float_t k_plen; //kaon track length
    //Float_t k_phi; //kaon track phi
    //Float_t k_theta; //kaon track theta
    //Int_t k_in_5_TPC; //kaon track is within 5 cm from the TPC edges
    //Int_t k_in_CC_TPC; //kaon track is within CC fiducial volume
    //Int_t k_hit; //number of kaon dEdx hits from calorimetry
    //Float_t k_range; //kaon range from calorimetry
    //Float_t k_KE; //kaon kinectic energy from calorimetry
    //Float_t k_large_dedx; //dedx median?
    //Float_t k_small_dedx; //dedx median?
    //Float_t k_chi_p; //kaon track proton chi2 score
    //Float_t k_chi_k; //kaon track kaon chi2 score
    //Float_t k_chi_pi; //kaon track pion chi2 score
    //Float_t k_chi_mu; //kaon track muon chi2 score
    //Float_t k_p_max; //kaon track proton bragg peak likelihood maximum forward or backward
    //Float_t k_k_max; //kaon track kaon bragg peak likelihood maximum forward or backward
    //Float_t k_pi_max; //kaon track pion bragg peak likelihood maximum forward or backward
    //Float_t k_mu_max; //kaon track muon bragg peak likelihood maximum forward or backward
    //Float_t k_mip_max; //kaon track mip bragg peak likelihood maximum forward or backward
    //Float_t k_L1_ratio; //kaon track mip / (proton+kaon)
    //Float_t k_LL1; //kaon track log( mip /(proton+kaon) )
    //Float_t k_L2_ratio; //kaon track mip+muon+pion / (mip+muon+pion+kaon)
    //Float_t k_LL2; //kaon track log( mip+muon+pion / (mip+muon+pion+kaon) )
    //Float_t k_Lp_ratio; //kaon track mip+muon / (mip+muon+proton)
    //Float_t k_LLp;
    //Float_t k_Lk_ratio; //kaon track mip+muon / (mip+muon+kaon)
    //Float_t k_LLk; //kaon track log ( (mip+muon) / (mip+muon+kaon) )
    //Float_t k_pida_mean; //PIDA
    //Float_t k_pida_med; //PIDA
    //Float_t k_kde; //PIDA
    //Float_t k_trm_dqdx; //truncated mean
    //Float_t k_trm_dedx; //truncated mean
    //Float_t k_dedx[3][3000]; //kaon track dedx per plane per hit
    //Float_t k_rr[3][3000]; //kaon track residual range per plane per hit
    //Float_t mu_dedx[3][3000]; //muon track dedx per plane per hit
    //Float_t mu_rr[3][3000]; //muon track residual range per plane per hit
    //
    //Int_t  mu_pdg; //muon track true matched particle pdg
    //Int_t  mu_isDec; //muon track true matched particle is true decay?
    //Int_t  mu_origin; //muon track true matched particle event origin (unknown=0,beam=1,cosmic=2,SN=3,single=4)
    //Int_t  mu_k_is_Mother; //muon track true matched particle is daugther of true kaon
    //Int_t  mu_mom_k_inelas; //muon track true matched particle from inelastic K interaction
    //Int_t  mu_ness; //muon track true matched particle is an interesting true muon (mu_origin==1 && mu_pdg==-13 && mu_isDec==1 && mu_k_is_Mother==1)
    //Float_t mu_plen; //muon track track length
    //Float_t mu_phi; //muon track phi
    //Float_t mu_theta; //muon track theta
    //Int_t   mu_in_5_TPC; //muon track is within 5 cm from TPC edges
    //Int_t   mu_in_CC_TPC; //muon track is within CC fiducial volume
    //Float_t mu_KE; //muon track reconstructed kinetic energy
    //Int_t   mu_hit; //muon track number of hits
    //Float_t mu_range; //muon range from calorimetry
    //Float_t mu_large_dedx;
    //Float_t mu_small_dedx;
    //Float_t mu_chi_p;
    //Float_t mu_chi_k;
    //Float_t mu_chi_pi;
    //Float_t mu_chi_mu;
    //Float_t mu_p_max;
    //Float_t mu_k_max;
    //Float_t mu_pi_max;
    //Float_t mu_mu_max;
    //Float_t mu_mip_max;
    //Float_t mu_L1_ratio;
    //Float_t mu_LL1;
    //Float_t mu_L2_ratio;
    //Float_t mu_LL2;
    //Float_t mu_Lp_ratio;
    //Float_t mu_LLp;
    //Float_t mu_Lk_ratio;
    //Float_t mu_LLk;
    //Float_t mu_pida_mean;
    //Float_t mu_pida_med;
    //Float_t mu_kde;
    //Float_t mu_trm_dqdx;
    //Float_t mu_trm_dedx;
    //std::string mu_mom_process; //muon track true matched particle creation process
    //
    //Int_t cc_mu_trkid;
    //Float_t cc_mu_tlen;
    //Float_t cc_mu_phi;
    //Float_t cc_mu_theta;
    //Float_t cc_mu_range;
    //Float_t cc_mu_KE;
    //Int_t   cc_mu_hit;
    //Float_t cc_mu_large_dedx;
    //Float_t cc_mu_small_dedx;
    //Float_t cc_dis_vtx;
    //Int_t cc_mu_pdg;
    //Float_t cc_mu_chi_p;
    //Float_t cc_mu_chi_k;
    //Float_t cc_mu_chi_pi;
    //Float_t cc_mu_chi_mu;
    //Float_t cc_mu_p_max;
    //Float_t cc_mu_k_max;
    //Float_t cc_mu_pi_max;
    //Float_t cc_mu_mu_max;
    //Float_t cc_mu_mip_max;
    //Float_t cc_mu_L1_ratio;
    //Float_t cc_mu_LL1;
    //Float_t cc_mu_L2_ratio;
    //Float_t cc_mu_LL2;
    //Float_t cc_mu_Lp_ratio;
    //Float_t cc_mu_LLp;
    //Float_t cc_mu_Lk_ratio;
    //Float_t cc_mu_LLk;
    //Float_t cc_mu_pida_mean;
    //Float_t cc_mu_pida_med;
    //Float_t cc_mu_kde;
    //Float_t cc_mu_trm_dqdx;
    //Float_t cc_mu_trm_dedx;
    //Int_t pri_Mu_is;
    //
    //Int_t Tr_pri_mu_pdg;
    //Int_t Tr_pri_mu_is;
    //Int_t Tr_pri_st_k_is;
    //Int_t Tr_K_Inelas;
    //Float_t Tr_k_plen;
    //Float_t Tr_k_endE;
    //Float_t Tr_k_theta;
    //Float_t Tr_k_phi;
    //Int_t Tr_dec_mu_is;
    //Int_t Tr_dec_mu_pi_pdg;
    //Float_t Tr_mu_plen;
    //Float_t Tr_mu_theta;
    //Float_t Tr_mu_phi;
    //Int_t Tr_k_inTPC;
    //Int_t Tr_mu_inTPC;
    //Int_t Tr_k_in_5_TPC;
    //Int_t Tr_k_in_CC_TPC;
    //Int_t Tr_mu_in_5_TPC;
    //Int_t Tr_mu_in_CC_TPC;
    //Float_t Tr_kmu_open_ang;
    //
    //Int_t longest_trkid;
    //Float_t longest_trklen;
    //Int_t vtx_5cm_mult;
    //Float_t k_start_dedx;
    //Float_t k_end_dedx;
    //Float_t mu_start_dedx;
    //Float_t mu_end_dedx;
    //
    //Int_t  kinelas_has_traks;
    //Int_t  kinelas_reco_trkID;
    //Float_t kinelas_tlen;
    //Float_t True_kinelas_KE;
    //Float_t True_kinelas_tlen;

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
    Int_t cut_13;

    //const trkf::TrackMomentumCalculator trkmom;

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
    std::string fPFParticleLabel;
  //    std::string m_pandoraLabel;
  //   std::string m_is_verbose;  

    //std::vector<Float_t> adEdx;
  std::vector<Float_t> dv0;
  std::vector<Float_t> rv0;
  std::vector<Float_t> dv1;
  std::vector<Float_t> rv1;
  std::vector<Float_t> dv2;
  std::vector<Float_t> rv2;



    bool isMC;

    //// Tree for the POT subrun info
    TTree *fSubrunTree;
    uint m_run, m_subrun;
    float m_pot;
    //float event_weight;

    std::map<std::string, std::vector<float>> event_weight;

    std::vector<std::string> evtwgt_funcname;          // the name of the functions used
    std::vector<std::vector<float>> evtwgt_weight;    // the weights (a vector for each function used)
    std::vector<int> evtwgt_nweight;                   // number of weights for each function
    //Int_t evtwgt_nfunc;                                // number of functions used

}; // class CCKaonAnalyzerRebuild
*/

//========================================================================
CCKaonAnalyzerRebuild::CCKaonAnalyzerRebuild(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","gaushit")),
  fLArG4ModuleLabel         (pset.get< std::string >("LArG4ModuleLabel","largeant")),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","generator")),  
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","pandora")),
  fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel","pandora")),
  //fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","pandoracaliSCE")),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","pandoraKalmanShowercali")),
  fPIDLabel                 (pset.get< std::string >("PIDLabel","pandorapid")),
  fHitTruthAssns            (pset.get< std::string >("HitTruthAssn","gaushitTruthMatch")), 
  fHitTrackAssns            (pset.get< std::string >("HitTrackAssn","pandora")), 
  fHitShowerAssns            (pset.get< std::string >("HitShowerAssn","pandora")), 
  m_pfp_producer            (pset.get< std::string >("pfp_producer","pandora")),
  fPFParticleLabel          (pset.get< std::string >("PFParticleLabel", "pandora")),
  fSpacePointproducer       (pset.get< std::string >("SpacePointproducer", "pandora")),
  //fSpacePointproducer  = p.get< art::InputTag >("SpacePointproducer");
//  m_pandoraLabel            (pset.get< std::string >("PandoraLabel")),
//  m_is_verbose              (pset.get<bool>("Verbose",false)),
  isMC                      (pset.get< bool >("IsMC",true))

  //reco_track_dEdx(nullptr),
  //reco_track_ResRan(nullptr)
{
  //  fm_piPFParticleLabel = pset.get<std::string>("PFParticleLabel");
  //  Kaon_Analyzer::CCKaonAnalyzerRebuild tmp;
  //tmp.m_pandoraLabel = pset.get<std::string>("PandoraLabel");
  //tmp.m_is_verbose = pset.get<bool>("Verbose",false);
  
  theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
  SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  geom = lar::providerFrom<geo::Geometry>();

        m_is_verbose = pset.get<bool>("Verbose",false);
        m_use_delaunay = pset.get<bool>("useDelaunay",false);
        m_is_data = pset.get<bool>("isData",false);
        m_is_overlayed = pset.get<bool>("isOverlayed",false);
        m_fill_trees = pset.get<bool>("FillTrees",true);
        m_run_pi0_filter = pset.get<bool>("RunPi0Filter",false);
        if(m_run_pi0_filter) m_is_data = true;// If running in filter mode, treat all as data

        m_pandoraLabel = pset.get<std::string>("PandoraLabel");
        m_trackLabel = pset.get<std::string>("TrackLabel");
	//        m_trackLabel_old = pset.get<std::string>("TrackLabelOld");
        m_sliceLabel = pset.get<std::string>("SliceLabel","pandora");
        m_showerLabel = pset.get<std::string>("ShowerLabel");
        m_caloLabel = pset.get<std::string>("CaloLabel");
        m_flashLabel = pset.get<std::string>("FlashLabel");

        m_hitfinderLabel = pset.get<std::string>("HitFinderModule", "gaushit");
        m_badChannelLabel = pset.get<std::string>("BadChannelLabel","badmasks");
        m_badChannelProducer = pset.get<std::string>("BadChannelProducer","simnfspl1");

        m_generatorLabel = pset.get<std::string>("GeneratorLabel","generator");
        m_mcTrackLabel = pset.get<std::string>("MCTrackLabel","mcreco");
        m_mcShowerLabel = pset.get<std::string>("MCShowerLabel","mcreco");
        m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");
        m_backtrackerLabel = pset.get<std::string>("BackTrackerModule","gaushitTruthMatch");
        m_hitMCParticleAssnsLabel = pset.get<std::string>("HitMCParticleAssnLabel","gaushitTruthMatch");

	m_CRTTzeroLabel = pset.get<std::string>("CRTTzeroLabel","crttzero");
	m_runCRT = pset.get<bool>("runCRT",false);
	m_CRTHitProducer = pset.get<std::string>("CRTHitProducer", "crthitcorr");

	m_gain_mc =pset.get<std::vector<double>>("gain_mc");
	m_wire_spacing = pset.get<double>("wire_spacing");
	m_width_dqdx_box = pset.get<double>("width_box");
	m_length_dqdx_box = pset.get<double>("length_box");
	m_truthmatching_signaldef = pset.get<std::string>("truthmatching_signaldef");
	m_pidLabel = pset.get<std::string>("ParticleIDLabel","particleid");
	m_shower3dLabel = pset.get<std::string>("Shower3DLabel","shrreco3d");

        m_run_all_pfps = pset.get<bool>("runAllPFPs",false);
        rangen = new TRandom3(22);
        bool_make_sss_plots = true;



  /*
  m_pandoraLabel = pset.get<std::string>("PandoraLabel");
  m_is_verbose = pset.get<bool>("Verbose",false);
  m_is_data = pset.get<bool>("isData",false); 
  m_use_delaunay = pset.get<bool>("useDelaunay",false);
  m_is_overlayed = pset.get<bool>("isOverlayed",false); 
  m_fill_trees = pset.get<bool>("FillTrees",true);   

  //m_mcTrackLabel = pset.get<std::string>("MCTrackLabel","mcreco");
  m_showerLabel = pset.get<std::string>("ShowerLabel");
  m_trackLabel = pset.get<std::string>("TrackLabel");
  //m_mcShowerLabel = pset.get<std::string>("MCShowerLabel","mcreco");
  m_shower3dLabel = pset.get<std::string>("Shower3DLabel","shrreco3d");  
  m_run_all_pfps = pset.get<bool>("runAllPFPs",false); 
  //m_truthmatching_signaldef = pset.get<std::string>("truthmatching_signaldef");

  m_hitfinderLabel = pset.get<std::string>("HitFinderModule", "gaushit");
  m_flashLabel = pset.get<std::string>("FlashLabel");
  m_shower3dLabel = pset.get<std::string>("Shower3DLabel","shrreco3d");
  m_generatorLabel = pset.get<std::string>("GeneratorLabel","generator");
  m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");
  m_hitMCParticleAssnsLabel = pset.get<std::string>("HitMCParticleAssnLabel","gaushitTruthMatch");
  m_truthmatching_signaldef = pset.get<std::string>("truthmatching_signaldef");
  */


  // set dedx pdf parameters
  llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
  llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
  llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

  llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
  llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
  llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

  llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
  llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
  llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);


  llr_pid_calculator_k.set_dedx_binning(0, protonmuon_parameters_k.dedx_edges_pl_0);
  llr_pid_calculator_k.set_par_binning(0, protonmuon_parameters_k.parameters_edges_pl_0);
  llr_pid_calculator_k.set_lookup_tables(0, protonmuon_parameters_k.dedx_pdf_pl_0);

  llr_pid_calculator_k.set_dedx_binning(1, protonmuon_parameters_k.dedx_edges_pl_1);
  llr_pid_calculator_k.set_par_binning(1, protonmuon_parameters_k.parameters_edges_pl_1);
  llr_pid_calculator_k.set_lookup_tables(1, protonmuon_parameters_k.dedx_pdf_pl_1);

  llr_pid_calculator_k.set_dedx_binning(2, protonmuon_parameters_k.dedx_edges_pl_2);
  llr_pid_calculator_k.set_par_binning(2, protonmuon_parameters_k.parameters_edges_pl_2);
  llr_pid_calculator_k.set_lookup_tables(2, protonmuon_parameters_k.dedx_pdf_pl_2);
}
 
//========================================================================
CCKaonAnalyzerRebuild::~CCKaonAnalyzerRebuild()
{
  //destructor
}

//========================================================================
void CCKaonAnalyzerRebuild::endSubRun(const art::SubRun &subrun)
{
  //if (!m_isData)
  //{
    art::Handle<sumdata::POTSummary> potSummaryHandle;
    m_pot = subrun.getByLabel("generator", potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    // -- std::cout << "[CCKaonAnalyzerRebuild::endSubRun] Storing POT info!" << std::endl;
  //}

  m_run = subrun.run();
  m_subrun = subrun.subRun();
  fSubrunTree->Fill();
}

//========================================================================
  void CCKaonAnalyzerRebuild::filter(const art::Event &evt)
{

  /*
  auto const TPC = (*geom).begin_TPC();
  auto ID = TPC.ID();
  m_Cryostat = ID.Cryostat;
  m_TPC = ID.TPC;
  */

  _time2cm = theDetector->SamplingRate() / 1000.0 * theDetector->DriftVelocity( theDetector->Efield(), theDetector->Temperature() );//found in ProtoShowerPandora_tool.cc  

  this->ClearVertex();


  //Collect the PFParticles from the event. This is the core!                                                                                                                                                                                                                                                                 
        // Collect all the hits. We will need these. Lets grab both the handle as well as a vector of art::Ptr as I like both. 
        art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(m_hitfinderLabel); 
        std::vector<art::Ptr<recob::Hit>> hitVector;
        art::fill_ptr_vector(hitVector,hitHandle);

        //Lets do "THE EXACT SAME STUFF" for Optical Flashes
        art::ValidHandle<std::vector<recob::OpFlash>> const & flashHandle  = evt.getValidHandle<std::vector<recob::OpFlash>>(m_flashLabel);
        std::vector<art::Ptr<recob::OpFlash>> flashVector;
        art::fill_ptr_vector(flashVector,flashHandle);

        //tracks
        art::ValidHandle<std::vector<recob::Track>> const & trackHandle  = evt.getValidHandle<std::vector<recob::Track>>(m_trackLabel);
        std::vector<art::Ptr<recob::Track>> trackVector;
        art::fill_ptr_vector(trackVector,trackHandle);

        //BadChannels
        art::Handle<std::vector<int> > badChannelHandle;
        std::vector<int> badChannelVector;
        if(evt.getByLabel(m_badChannelProducer, m_badChannelLabel, badChannelHandle)){
            badChannelVector            = *(badChannelHandle);
        }


        //Collect the PFParticles from the event. This is the core!

  art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);
  std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
  art::fill_ptr_vector(pfParticleVector,pfParticleHandle);

  //get the cluster handle for the dQ/dx calc
  art::ValidHandle<std::vector<recob::Cluster>> const & clusterHandle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pandoraLabel);
  std::vector< art::Ptr<recob::Cluster> > clusterVector;
  art::fill_ptr_vector(clusterVector,clusterHandle);


    //So a cross check
    /*                                                                                                                                                                                                                                                                                                                      
    if (!pfParticleHandle.isValid())
      {
	mf::LogDebug("SinglePhoton") << "  Failed to find the PFParticles.\n";
	if(m_run_pi0_filter)
	  return false;
	else
	  return true;
      }
    */

    //This is another pandora helper. I don't like PFParticle ID lookups but I guess lets keep for now;                                                                                                                                                                                                                                     
    // Produce a map of the PFParticle IDs for fast navigation through the hierarchy                                                                                                                                                                                                                                                        
    //    typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;
    //    CCKaonAnalyzerRebuild::
  
    
    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);
  
    //Slices                                                                                                                                                                                                                                                                                                                                
    art::ValidHandle<std::vector<recob::Slice>> const & sliceHandle  = evt.getValidHandle<std::vector<recob::Slice>>(m_pandoraLabel);
    std::vector<art::Ptr<recob::Slice>> sliceVector;
    art::fill_ptr_vector(sliceVector,sliceHandle);
    //And some associations                                                                                                                                                                                                                                                                                                                 
    art::FindManyP<recob::PFParticle> pfparticles_per_slice(sliceHandle, evt, m_pandoraLabel);
    art::FindManyP<recob::Hit> hits_per_slice(sliceHandle, evt, m_pandoraLabel);

    std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::PFParticle>> > sliceToPFParticlesMap;
    std::map<int, std::vector<art::Ptr<recob::PFParticle>> > sliceIDToPFParticlesMap;
    for(size_t i=0; i< sliceVector.size(); ++i){
      auto slice = sliceVector[i];
      sliceToPFParticlesMap[slice] =pfparticles_per_slice.at(slice.key());
      sliceIDToPFParticlesMap[slice->ID()] = pfparticles_per_slice.at(slice.key());
    }

    std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::Hit>> > sliceToHitsMap;
    std::map<int, std::vector<art::Ptr<recob::Hit>> > sliceIDToHitsMap;
    for(size_t i=0; i< sliceVector.size(); ++i){
      auto slice = sliceVector[i];
      sliceToHitsMap[slice] =hits_per_slice.at(slice.key());
      sliceIDToHitsMap[slice->ID()] = hits_per_slice.at(slice.key());
    }

        //And some verticies.        
        art::ValidHandle<std::vector<recob::Vertex>> const & vertexHandle = evt.getValidHandle<std::vector<recob::Vertex>>(m_pandoraLabel);
        std::vector<art::Ptr<recob::Vertex>> vertexVector;
        art::fill_ptr_vector(vertexVector,vertexHandle);
        art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
        std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > pfParticlesToVerticesMap;
        for(size_t i=0; i< pfParticleVector.size(); ++i){
            auto pfp = pfParticleVector[i];
            pfParticlesToVerticesMap[pfp] =vertices_per_pfparticle.at(pfp.key());
        }

        //------- 3D showers
	
        art::FindOneP<recob::Shower> showerreco3D_per_pfparticle(pfParticleHandle, evt, m_shower3dLabel);
        std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>> pfParticlesToShowerReco3DMap;
        for(size_t i=0; i< pfParticleVector.size(); ++i){
            auto pfp = pfParticleVector[i];
            if(!showerreco3D_per_pfparticle.at(pfp.key()).isNull()){
                pfParticlesToShowerReco3DMap[pfp] = showerreco3D_per_pfparticle.at(pfp.key());
            }

        }
	

        // Once we have actual verticies, lets concentrate on JUST the neutrino PFParticles for now:
        //--------------------------------
        // Produce two PFParticle vectors containing final-state particles:
        // 1. Particles identified as cosmic-rays - recontructed under cosmic-hypothesis
        // 2. Daughters of the neutrino PFParticle - reconstructed under the neutrino hypothesis
        std::vector< art::Ptr<recob::PFParticle> > crParticles;
        std::vector< art::Ptr<recob::PFParticle> > nuParticles;
        this->GetFinalStatePFParticleVectors(pfParticleMap, pfParticlesToVerticesMap, crParticles, nuParticles);


        if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Spacepoints"<<std::endl;
        //Look, here is a map that I just forced myself rather than build using helpers. Not that different is it. But for somereason I only use PFParticles.. huh,
        //Spacepoint associaitions
        art::FindManyP<recob::SpacePoint> spacePoints_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
        std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap;
        for(size_t i=0; i< nuParticles.size(); ++i){
            const art::Ptr<recob::PFParticle> pfp = nuParticles[i];
            pfParticleToSpacePointsMap[pfp] = spacePoints_per_pfparticle.at(pfp.key());
        }


    if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get PandoraMetadata"<<std::endl;
    //add the associaton between PFP and metadata, this is important to look at the slices and scores                                                         
              
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt, m_pandoraLabel);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > pfParticleToMetadataMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
      const art::Ptr<recob::PFParticle> pfp = pfParticleVector[i];
      pfParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
    }


        if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Clusters"<<std::endl;

        //Get a map between the PFP's and the clusters. Although Mark isn't a fan of clusters, they're imporant for the shower dQ/dx
        //Also need a map between clusters and hits
        art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
        art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, m_pandoraLabel);
        std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap;
        std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap;
        //fill map PFP to Clusters
        for(size_t i=0; i< nuParticles.size(); ++i){
            auto pfp = nuParticles[i];
            pfParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());
            // pfParticleToSpacePointsMap[pfp] = spacePoints_per_pfparticle.at(pfp.key());
            // pfParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
        }
        //fill map Cluster to Hits
        for(size_t i=0; i< clusterVector.size(); ++i){
            auto cluster = clusterVector[i];
            clusterToHitsMap[cluster] = hits_per_cluster.at(cluster.key());
        }

        if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Build hits to PFP Maps"<<std::endl;



        //taking out the Larpandora helper functions here because they don't match to non-neutrino slice hits for some reason

        //OK Here we build two IMPORTANT maps for the analysis, (a) given a PFParticle get a vector of hits..
        //and (b) given a single hit, get the PFParticle it is in (MARK: is it only one? always? RE-MARK: Yes)
        std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
        //        std::map<art::Ptr<recob::Hit>, art::Ptr<recob::PFParticle>>                hitToPFParticleMap;
        //Using a pandora helper here, but to be honest we should probably just build using normal associations so keep independant if pssoble
        // lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(evt, m_pandoraLabel, pfParticleToHitsMap, hitToPFParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters);

        //use pfp->cluster and cluster->hit to build pfp->hit map
        //for each PFP
        for(size_t i=0; i<nuParticles.size(); ++i){
            auto pfp = nuParticles[i];

            // std::cout<<"starting to match to hits for pfp "<<pfp->Self()<<std::endl;
            //get the associated clusters
            std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[pfp] ;

            //make empty vector to store hits
            std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

            // std::cout<<"-- there are "<<clusters_vec.size()<<" associated clusters"<<std::endl;

            //for each cluster, get the associated hits
            for (art::Ptr<recob::Cluster> cluster: clusters_vec){
                std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];

                //   std::cout<<"looking at cluster in pfp "<<pfp->Self()<<" with "<<hits_vec.size() <<" hits"<<std::endl;
                //insert hits into vector
                hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
            }

            //fill the map
            pfParticleToHitsMap[pfp] = hits_for_pfp;
            //std::cout<<"saving a total of "<<hits_for_pfp.size()<<" hits for pfp "<<pfp->Self()<<std::endl;

        }//for each pfp


	    
    //these are all filled in analyze slice                                                                                                                                                                                                                                                                                             
    std::vector< art::Ptr<recob::Track> > tracks;
    std::vector< art::Ptr<recob::Shower> > showers;
    std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle >> trackToNuPFParticleMap;
    std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap;


    if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Tracks and Showers"<<std::endl; 
    this->CollectTracksAndShowers(nuParticles, pfParticleMap,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);




        //**********************************************************************************************/
        //**********************************************************************************************/
        //---------------------------------- MC TRUTH Data Only---------------------------
        //**********************************************************************************************/
        //**********************************************************************************************/

        //Get the MCtruth handles and vectors
        std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;
        std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

        //Then build a map from MCparticles to Hits and vice versa
        std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  mcParticleToHitsMap;
        std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  hitToMCParticleMap;


        //Apparrently a MCParticle doesn't know its origin (thanks Andy!)
        //I would also like a map from MCparticle to MCtruth and then I will be done.  and Vice Versa
        //Note which map is which!       //First  is one-to-many.         //Second is one-to-one
        std::map< art::Ptr<simb::MCTruth>,    std::vector<art::Ptr<simb::MCParticle>>>  MCTruthToMCParticlesMap;
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>                  MCParticleToMCTruthMap;
        std::map<int, art::Ptr<simb::MCParticle> >                                     MCParticleToTrackIdMap;

        std::vector<art::Ptr<sim::MCTrack>> mcTrackVector;
        std::vector<art::Ptr<sim::MCShower>> mcShowerVector;

        std::vector<art::Ptr<simb::MCParticle>> matchedMCParticleVector;
        std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > trackToMCParticleMap;
        std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > showerToMCParticleMap;

        //Given a simb::MCParticle we would like a map to either a sim::MCTrack or sim::MCShower
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCTrack> > MCParticleToMCTrackMap;
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCShower> > MCParticleToMCShowerMap;


        //**********************************************************************************************/
        //**********************************************************************************************/
        //Some event based properties



    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind
    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> primaryPFPSliceIdVec; //stores a pair of only the primary PFP's in the event and the slice ind                                                                                                                                                                              
    std::map<int, double> sliceIdToNuScoreMap; //map between a slice Id and neutrino score                                               
    std::map<art::Ptr<recob::PFParticle>, bool> PFPToClearCosmicMap; //returns true for clear cosmic, false otherwise                    
    std::map<art::Ptr<recob::PFParticle>, int> PFPToSliceIdMap; //returns the slice id for all PFP's                                     
    std::map<art::Ptr<recob::PFParticle>,bool> PFPToNuSliceMap;
    std::map<art::Ptr<recob::PFParticle>,double> PFPToTrackScoreMap;
    std::map<int, int> sliceIdToNumPFPsMap;
    std::cout<<"SinglePhoton::analyze::AnalyzeSlice()\t||\t Starting"<<std::endl;


    this->AnalyzeSlices(pfParticleToMetadataMap, pfParticleMap,  primaryPFPSliceIdVec, sliceIdToNuScoreMap, PFPToClearCosmicMap, PFPToSliceIdMap, PFPToNuSliceMap, PFPToTrackScoreMap);
    //std::cout<<"There are "<< allPFPSliceIdVec.size()<<" pfp-slice id matches stored in the vector"<<std::endl;
    std::cout<<"SinglePhoton::analyze\t||\tthe number of PPF's with stored clear cosmic info is "<<PFPToClearCosmicMap.size()<<std::endl;
    std::cout<<"SinglePhoton::analyze\t||\tthe number of PFP's stored in the PFPToSliceIdMap is "<<PFPToSliceIdMap.size()<<std::endl;

    if (PFPToSliceIdMap.size() < 1){
      std::cout<<"ERROR, not storing PFP's in PFPToSliceIdMap"<<std::endl;
    }


    for (auto pair:sliceIDToPFParticlesMap){ 
      std::vector<art::Ptr<recob::PFParticle>> pfp_vec = pair.second;
      int slice_id = pair.first;
      //if (slice_vec[0]->Slice() != PFPToSliceIdMap[pfp] )                                                                                                                                                                                                                                                                           
      for(auto pfp: pfp_vec){
	if (slice_id != PFPToSliceIdMap[pfp] && PFPToSliceIdMap[pfp]>=0){
	  std::cout<<"sliceIDToPFParticlesMap[slice->ID()] for pfp "<<pfp->Self()<<" is slice "<< slice_id<< "but PFPToSliceIdMap[pfp] = "<<PFPToSliceIdMap[pfp]<<std::endl;
	}
      }

    }


    //if CRT info, get CRT hits  
    art::Handle<std::vector<crt::CRTHit>> crthit_h; //only filled when there are hits, otherwise empty 
    art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
    double evt_timeGPS_nsec = -999;
    if(m_runCRT){
      evt.getByLabel(m_DAQHeaderProducer, rawHandle_DAQHeader);
      evt.getByLabel(m_CRTHitProducer, crthit_h);
      raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
      art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
      evt_timeGPS_nsec = evtTimeGPS.timeLow();

      std::cout<<"SinglePhoton::analyze \t||\t Got CRT hits"<<std::endl;
    }

    this->AnalyzeFlashes(flashVector, crthit_h, evt_timeGPS_nsec);

    this->AnalyzeTracks(tracks, trackToNuPFParticleMap, pfParticleToSpacePointsMap,  MCParticleToTrackIdMap, sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap,  PFPToTrackScoreMap, PFPToNuSliceMap,pfParticleMap);

    this->AnalyzeShowers(showers,showerToNuPFParticleMap, pfParticleToHitsMap, pfParticleToClustersMap, clusterToHitsMap,sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap, PFPToNuSliceMap, PFPToTrackScoreMap,pfParticleMap,pfParticlesToShowerReco3DMap); 
    



        if(!m_is_data){

            art::ValidHandle<std::vector<simb::MCTruth>> const & mcTruthHandle= evt.getValidHandle<std::vector<simb::MCTruth>>(m_generatorLabel);
            art::fill_ptr_vector(mcTruthVector,mcTruthHandle);
            art::ValidHandle<std::vector<simb::MCParticle>> const & mcParticleHandle= evt.getValidHandle<std::vector<simb::MCParticle>>(m_geantModuleLabel);
            art::fill_ptr_vector(mcParticleVector,mcParticleHandle);

	    this->CollectMCParticles(evt, m_geantModuleLabel, MCTruthToMCParticlesMap, MCParticleToMCTruthMap, MCParticleToTrackIdMap);

            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcparticles_per_hit(hitHandle, evt, m_hitMCParticleAssnsLabel);

            //mcc9 march miniretreat fix
            std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP
            std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing

            m_test_matched_hits = 0;

            for(size_t j=0; j<hitVector.size();j++){
                const art::Ptr<recob::Hit> hit = hitVector[j];

                particle_vec.clear(); match_vec.clear(); //only store per hit

                mcparticles_per_hit.get(hit.key(), particle_vec, match_vec);

                if(particle_vec.size() > 0){
                    m_test_matched_hits++;
                }

            }

            this->BuildMCParticleHitMaps(evt, m_geantModuleLabel, hitVector,  mcParticleToHitsMap, hitToMCParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters,  MCParticleToTrackIdMap);

            //std::cout<<"SinglePhoton\t||\t Starting backtracker on recob::track"<<std::endl;
            //std::vector<double> trk_overlay_vec = recoMCmatching<art::Ptr<recob::Track>>( tracks, trackToMCParticleMap, trackToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, matchedMCParticleVector);

    std::cout<<"SinglePhoton\t||\t Starting backtracker on recob::shower"<<std::endl;
    this->showerRecoMCmatching(showers, showerToMCParticleMap, showerToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, matchedMCParticleVector, pfParticleMap,  MCParticleToTrackIdMap, sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap, PFPToNuSliceMap);



    std::cout<<"filling info in ncdelta slice tree"<<std::endl;
    this->AnalyzeRecoMCSlices( m_truthmatching_signaldef, MCParticleToTrackIdMap, showerToNuPFParticleMap , allPFPSliceIdVec, showerToMCParticleMap, trackToNuPFParticleMap, trackToMCParticleMap,  PFPToSliceIdMap);
	}

    /*
    this->AnalyzeSlices(pfParticleToMetadataMap, pfParticleMap,  primaryPFPSliceIdVec, sliceIdToNuScoreMap, PFPToClearCosmicMap, PFPToSliceIdMap, PFPToNuSliceMap, PFPToTrackScoreMap);
    //std::cout<<"There are "<< allPFPSliceIdVec.size()<<" pfp-slice id matches stored in the vector"<<std::endl;                                                                                                                                                                                                                       
    std::cout<<"SinglePhoton::analyze\t||\tthe number of PPF's with stored clear cosmic info is "<<PFPToClearCosmicMap.size()<<std::endl;
    std::cout<<"SinglePhoton::analyze\t||\tthe number of PFP's stored in the PFPToSliceIdMap is "<<PFPToSliceIdMap.size()<<std::endl;
    if (PFPToSliceIdMap.size() < 1){
      std::cout<<"ERROR, not storing PFP's in PFPToSliceIdMap"<<std::endl;
    }
    */

    for (auto pair:PFPToNuSliceMap){
      auto pfp = pair.first;
      auto is_nuslice = pair.second;
      if (is_nuslice){
	std::cout<<"pfp in nuslice "<<pfp->Self()<<std::endl;
      }

    }


    //---------------------- END OF LOOP, fill vertex ---------------------
    ncdelta_slice_tree->Fill();
    vertex_tree->Fill();
    eventweight_tree->Fill();


}

//========================================================================
void CCKaonAnalyzerRebuild::beginJob()
{

  //initialiseCanvas();

  art::ServiceHandle<art::TFileService> tfs;

  fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");

  //gInterpreter->GenerateDictionary("vector<vector<vector<Float_t>> >", "vector");
  fEventTree->Branch("event", &event, "event/I");
  fEventTree->Branch("run", &run, "run/I");
  fEventTree->Branch("subrun", &subrun, "surbrun/I");

  fEventTree->Branch("IsKaon", &IsKaon, "IsKaon/I");
  fEventTree->Branch("IsSingleKaon", &IsSingleKaon, "IsSingleKaon/I");
  fEventTree->Branch("IsAssociatedKaon", &IsAssociatedKaon, "IsAssociatedKaon/I");
  fEventTree->Branch("IsMuBR", &IsMuBR, "IsMuBR/I");
  fEventTree->Branch("IsPiBR", &IsPiBR, "IsPiBR/I");
  fEventTree->Branch("prip", &prip, "prip/I");
  fEventTree->Branch("prip_k_dau", &prip_k_dau, "prip_k_dau/I");


  fEventTree->Branch("process",&Process);
  fEventTree->Branch("endprocess",&EndProcess);

  fEventTree->Branch("event_weight" ,&event_weight);

  fEventTree->Branch("evtwgt_funcname" ,&evtwgt_funcname);
  fEventTree->Branch("evtwgt_weight" ,&evtwgt_weight);
  fEventTree->Branch("evtwgt_nweight" ,&evtwgt_nweight);

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

  fEventTree->Branch("true_lepton_pdg", &true_lepton_pdg, "true_lepton_pdg/I");
  fEventTree->Branch("true_lepton_p", &true_lepton_p, "true_lepton_p/F");
  fEventTree->Branch("true_lepton_ke", &true_lepton_ke, "true_lepton_ke/F");
  fEventTree->Branch("true_lepton_theta", &true_lepton_theta, "true_lepton_theta/F");
  fEventTree->Branch("true_lepton_costheta", &true_lepton_costheta, "true_lepton_costheta/F");
  fEventTree->Branch("true_lepton_phi", &true_lepton_phi, "true_lepton_phi/F");

  fEventTree->Branch("true_nkaons", &true_nkaons, "true_nkaons/I");

  fEventTree->Branch("true_kaon_length", &true_kaon_length, "true_kaon_length/F");
  fEventTree->Branch("true_kaon_p", &true_kaon_p, "true_kaon_p/F");
  fEventTree->Branch("true_kaon_ke", &true_kaon_ke, "true_kaon_ke/F");
  fEventTree->Branch("true_kaon_theta", &true_kaon_theta, "true_kaon_theta/F");
  fEventTree->Branch("true_kaon_costheta", &true_kaon_costheta, "true_kaon_costheta/F");
  fEventTree->Branch("true_kaon_phi", &true_kaon_phi, "true_kaon_phi/F");
  fEventTree->Branch("true_kaon_ccmuon_angle", &true_kaon_ccmuon_angle, "true_kaon_ccmuon_angle/F");
  fEventTree->Branch("true_kaon_ccmuon_cosangle", &true_kaon_ccmuon_cosangle, "true_kaon_ccmuon_cosangle/F");
  fEventTree->Branch("true_kaon_end_process", &true_kaon_end_process, "true_kaon_end_process/I");
  fEventTree->Branch("true_kaon_end_ke", &true_kaon_end_ke, "true_kaon_end_ke/F");
  fEventTree->Branch("true_kaon_start_x", &true_kaon_start_x, "true_kaon_start_x/F");
  fEventTree->Branch("true_kaon_start_y", &true_kaon_start_y, "true_kaon_start_y/F");
  fEventTree->Branch("true_kaon_start_z", &true_kaon_start_z, "true_kaon_start_z/F");
  fEventTree->Branch("true_kaon_end_x", &true_kaon_end_x, "true_kaon_end_x/F");
  fEventTree->Branch("true_kaon_end_y", &true_kaon_end_y, "true_kaon_end_y/F");
  fEventTree->Branch("true_kaon_end_z", &true_kaon_end_z, "true_kaon_end_z/F");
  fEventTree->Branch("true_kaon_end_inTPC", &true_kaon_end_inTPC, "true_kaon_end_inTPC/O");
  fEventTree->Branch("true_kaon_end_in5cmTPC", &true_kaon_end_in5cmTPC, "true_kaon_end_in5cmTPC/O");
  fEventTree->Branch("true_kaon_end_inCCInclusiveTPC", &true_kaon_end_inCCInclusiveTPC, "true_kaon_end_inCCInclusiveTPC/O");


  fEventTree->Branch("true_dau_muon_length", &true_dau_muon_length, "true_dau_muon_length/F");
  fEventTree->Branch("true_dau_muon_p", &true_dau_muon_p, "true_dau_muon_p/F");
  fEventTree->Branch("true_dau_muon_ke", &true_dau_muon_ke, "true_dau_muon_ke/F");
  fEventTree->Branch("true_dau_muon_theta", &true_dau_muon_theta, "true_dau_muon_theta/F");
  fEventTree->Branch("true_dau_muon_costheta", &true_dau_muon_costheta, "true_dau_muon_costheta/F");
  fEventTree->Branch("true_dau_muon_phi", &true_dau_muon_phi, "true_dau_muon_phi/F");
  fEventTree->Branch("true_dau_muon_ccmuon_angle", &true_dau_muon_ccmuon_angle, "true_dau_muon_ccmuon_angle/F");
  fEventTree->Branch("true_dau_muon_ccmuon_cosangle", &true_dau_muon_ccmuon_cosangle, "true_dau_muon_ccmuon_cosangle/F");
  fEventTree->Branch("true_dau_muon_end_process", &true_dau_muon_end_process, "true_dau_muon_end_process/I");
  fEventTree->Branch("true_dau_muon_end_ke", &true_dau_muon_end_ke, "true_dau_muon_end_ke/F");
  fEventTree->Branch("true_dau_muon_start_x", &true_dau_muon_start_x, "true_dau_muon_start_x/F");
  fEventTree->Branch("true_dau_muon_start_y", &true_dau_muon_start_y, "true_dau_muon_start_y/F");
  fEventTree->Branch("true_dau_muon_start_z", &true_dau_muon_start_z, "true_dau_muon_start_z/F");
  fEventTree->Branch("true_dau_muon_end_x", &true_dau_muon_end_x, "true_dau_muon_end_x/F");
  fEventTree->Branch("true_dau_muon_end_y", &true_dau_muon_end_y, "true_dau_muon_end_y/F");
  fEventTree->Branch("true_dau_muon_end_z", &true_dau_muon_end_z, "true_dau_muon_end_z/F");

  fEventTree->Branch("true_dau_pip_length", &true_dau_pip_length, "true_dau_pip_length/F");
  fEventTree->Branch("true_dau_pip_p", &true_dau_pip_p, "true_dau_pip_p/F");
  fEventTree->Branch("true_dau_pip_ke", &true_dau_pip_ke, "true_dau_pip_ke/F");
  fEventTree->Branch("true_dau_pip_theta", &true_dau_pip_theta, "true_dau_pip_theta/F");
  fEventTree->Branch("true_dau_pip_costheta", &true_dau_pip_costheta, "true_dau_pip_costheta/F");
  fEventTree->Branch("true_dau_pip_phi", &true_dau_pip_phi, "true_dau_pip_phi/F");
  fEventTree->Branch("true_dau_pip_ccmuon_angle", &true_dau_pip_ccmuon_angle, "true_dau_pip_ccmuon_angle/F");
  fEventTree->Branch("true_dau_pip_ccmuon_cosangle", &true_dau_pip_ccmuon_cosangle, "true_dau_pip_ccmuon_cosangle/F");
  fEventTree->Branch("true_dau_pip_end_process", &true_dau_pip_end_process, "true_dau_pip_end_process/I");
  fEventTree->Branch("true_dau_pip_end_ke", &true_dau_pip_end_ke, "true_dau_pip_end_ke/F");
  fEventTree->Branch("true_dau_pip_start_x", &true_dau_pip_start_x, "true_dau_pip_start_x/F");
  fEventTree->Branch("true_dau_pip_start_y", &true_dau_pip_start_y, "true_dau_pip_start_y/F");
  fEventTree->Branch("true_dau_pip_start_z", &true_dau_pip_start_z, "true_dau_pip_start_z/F");
  fEventTree->Branch("true_dau_pip_end_x", &true_dau_pip_end_x, "true_dau_pip_end_x/F");
  fEventTree->Branch("true_dau_pip_end_y", &true_dau_pip_end_y, "true_dau_pip_end_y/F");
  fEventTree->Branch("true_dau_pip_end_z", &true_dau_pip_end_z, "true_dau_pip_end_z/F");


  fEventTree->Branch("true_dau_pin_length", &true_dau_pin_length, "true_dau_pin_length/F");
  fEventTree->Branch("true_dau_pin_p", &true_dau_pin_p, "true_dau_pin_p/F");
  fEventTree->Branch("true_dau_pin_ke", &true_dau_pin_ke, "true_dau_pin_ke/F");
  fEventTree->Branch("true_dau_pin_theta", &true_dau_pin_theta, "true_dau_pin_theta/F");
  fEventTree->Branch("true_dau_pin_costheta", &true_dau_pin_costheta, "true_dau_pin_costheta/F");
  fEventTree->Branch("true_dau_pin_phi", &true_dau_pin_phi, "true_dau_pin_phi/F");
  fEventTree->Branch("true_dau_pin_ccmuon_angle", &true_dau_pin_ccmuon_angle, "true_dau_pin_ccmuon_angle/F");
  fEventTree->Branch("true_dau_pin_ccmuon_cosangle", &true_dau_pin_ccmuon_cosangle, "true_dau_pin_ccmuon_cosangle/F");
  fEventTree->Branch("true_dau_pin_end_process", &true_dau_pin_end_process, "true_dau_pin_end_process/I");
  fEventTree->Branch("true_dau_pin_end_ke", &true_dau_pin_end_ke, "true_dau_pin_end_ke/F");
  fEventTree->Branch("true_dau_pin_start_x", &true_dau_pin_start_x, "true_dau_pin_start_x/F");
  fEventTree->Branch("true_dau_pin_start_y", &true_dau_pin_start_y, "true_dau_pin_start_y/F");
  fEventTree->Branch("true_dau_pin_start_z", &true_dau_pin_start_z, "true_dau_pin_start_z/F");
  fEventTree->Branch("true_dau_pin_end_x", &true_dau_pin_end_x, "true_dau_pin_end_x/F");
  fEventTree->Branch("true_dau_pin_end_y", &true_dau_pin_end_y, "true_dau_pin_end_y/F");
  fEventTree->Branch("true_dau_pin_end_z", &true_dau_pin_end_z, "true_dau_pin_end_z/F");

  fEventTree->Branch("cheat_num_hits", &cheat_num_hits, "cheat_num_hits[20][10]/I");
  fEventTree->Branch("cheat_ini_hits", &cheat_ini_hits, "cheat_ini_hits[20][10]/I");
  fEventTree->Branch("cheat_closest_distance", &cheat_closest_distance, "cheat_closest_distance[20][10]/F");

  fEventTree->Branch("cheat_peak_pdg", &cheat_peak_pdg, "cheat_peak_pdg[20][10]/F");
  fEventTree->Branch("cheat_peak_theta", &cheat_peak_theta, "cheat_peak_theta[20][10]/F");
  fEventTree->Branch("cheat_pip_trkln", &cheat_pip_trkln, "cheat_pip_trkln/F");
  fEventTree->Branch("cheat_mup_trkln", &cheat_mup_trkln, "cheat_mup_trkln/F");

  fEventTree->Branch("cheat_peak_phi", &cheat_peak_phi, "cheat_peak_phi[20][10]/F");
  fEventTree->Branch("best_peak_theta", &best_peak_theta, "best_peak_theta[20][10]/F");
  fEventTree->Branch("best_peak_phi", &best_peak_phi, "best_peak_phi[20][10]/F");
  fEventTree->Branch("best_peak_trkln", &best_peak_trkln, "best_peak_trkln[20][10]/F");

  fEventTree->Branch("true_length", &true_length, "true_length[10]/F");
  fEventTree->Branch("true_p", &true_p, "true_p[10]/F");
  fEventTree->Branch("true_ke", &true_ke, "true_ke[10]/F");
  fEventTree->Branch("true_theta", &true_theta, "true_theta[10]/F");
  fEventTree->Branch("true_costheta", &true_costheta, "true_costheta[10]/F");
  fEventTree->Branch("true_phi", &true_phi, "true_phi[10]/F");
  fEventTree->Branch("true_ccmuon_angle", &true_ccmuon_angle, "true_ccmuon_angle[10]/F");
  fEventTree->Branch("true_ccmuon_cosangle", &true_ccmuon_cosangle, "true_ccmuon_cosangle[10]/F");
  fEventTree->Branch("true_end_process", &true_end_process, "true_end_process[10]/I");
  fEventTree->Branch("true_end_ke", &true_end_ke, "true_end_ke[10]/F");
  fEventTree->Branch("true_end_x", &true_end_x, "true_end_x[10]/F");
  fEventTree->Branch("true_end_y", &true_end_y, "true_end_y[10]/F");
  fEventTree->Branch("true_end_z", &true_end_z, "true_end_z[10]/F");
  fEventTree->Branch("true_end_inTPC", &true_end_inTPC, "true_end_inTPC[10]/O");
  fEventTree->Branch("true_end_in5cmTPC", &true_end_in5cmTPC, "true_end_in5cmTPC[10]/O");
  fEventTree->Branch("true_end_inCCInclusiveTPC", &true_end_inCCInclusiveTPC, "true_end_inCCInclusiveTPC[10]/O");

  fEventTree->Branch("true_kaon_ndaughters", &true_kaon_ndaughters, "true_kaon_ndaughters/I");
  fEventTree->Branch("true_kaon_ndaughters_decay", &true_kaon_ndaughters_decay, "true_kaon_ndaughters_decay/I");
  fEventTree->Branch("true_kaon_ndaughters_inelastic", &true_kaon_ndaughters_inelastic, "true_kaon_ndaughters_inelastic/I");
  fEventTree->Branch("true_kaon_ndecmup", &true_kaon_ndecmup, "true_kaon_ndecmup/I");
  fEventTree->Branch("true_kaon_ndecpip", &true_kaon_ndecpip, "true_kaon_ndecpip/I");
  fEventTree->Branch("true_kaon_ninekap", &true_kaon_ninekap, "true_kaon_ninekap/I");
  fEventTree->Branch("true_kaon_ninepip", &true_kaon_ninepip, "true_kaon_ninepip/I");
  fEventTree->Branch("true_kaon_ninepro", &true_kaon_ninepro, "true_kaon_ninepro/I");
  fEventTree->Branch("true_kaon_daughter_length", &true_kaon_daughter_length, "true_kaon_daughter_length/F");
  fEventTree->Branch("true_kaon_daughter_p", &true_kaon_daughter_p, "true_kaon_daughter_p/F");
  fEventTree->Branch("true_kaon_daughter_ke", &true_kaon_daughter_ke, "true_kaon_daughter_ke/F");
  fEventTree->Branch("true_kaon_daughter_theta", &true_kaon_daughter_theta, "true_kaon_daughter_theta/F");
  fEventTree->Branch("true_kaon_daughter_costheta", &true_kaon_daughter_costheta, "true_kaon_daughter_costheta/F");
  fEventTree->Branch("true_kaon_daughter_angle", &true_kaon_daughter_angle, "true_kaon_daughter_angle/F");
  fEventTree->Branch("true_kaon_daughter_cosangle", &true_kaon_daughter_cosangle, "true_kaon_daughter_cosangle/F");
  fEventTree->Branch("true_kaon_daughter_pdg", &true_kaon_daughter_pdg, "true_kaon_daughter_pdg/I");
  fEventTree->Branch("true_kaon_daughter_start_x", &true_kaon_daughter_start_x, "true_kaon_daughter_start_x/F");
  fEventTree->Branch("true_kaon_daughter_start_y", &true_kaon_daughter_start_y, "true_kaon_daughter_start_y/F");
  fEventTree->Branch("true_kaon_daughter_start_z", &true_kaon_daughter_start_z, "true_kaon_daughter_start_z/F");
  fEventTree->Branch("true_kaon_daughter_end_x", &true_kaon_daughter_end_x, "true_kaon_daughter_end_x/F");
  fEventTree->Branch("true_kaon_daughter_end_y", &true_kaon_daughter_end_y, "true_kaon_daughter_end_y/F");
  fEventTree->Branch("true_kaon_daughter_end_z", &true_kaon_daughter_end_z, "true_kaon_daughter_end_z/F");
  fEventTree->Branch("true_kaon_daughter_end_inTPC", &true_kaon_daughter_end_inTPC, "true_kaon_daughter_endInTPC/O");
  fEventTree->Branch("true_kaon_daughter_end_in5cmTPC", &true_kaon_daughter_end_in5cmTPC, "true_kaon_daughter_end_in5cmTPC/O");
  fEventTree->Branch("true_kaon_daughter_end_inCCInclusiveTPC", &true_kaon_daughter_end_inCCInclusiveTPC, "true_kaon_daughter_end_inCCInclusiveTPC/O");

  fEventTree->Branch("true_nhyperons", &true_nhyperons, "true_nhyperons/I");

  fEventTree->Branch("true_nkaons_anti", &true_nkaons_anti, "true_nkaons_anti/I");
  fEventTree->Branch("true_kaon_length_anti", &true_kaon_length_anti, "true_kaon_length_anti/F");
  fEventTree->Branch("true_kaon_p_anti", &true_kaon_p_anti, "true_kaon_p_anti/F");
  fEventTree->Branch("true_kaon_ke_anti", &true_kaon_ke_anti, "true_kaon_ke_anti/F");
  fEventTree->Branch("true_kaon_theta_anti", &true_kaon_theta_anti, "true_kaon_theta_anti/F");
  fEventTree->Branch("true_kaon_costheta_anti", &true_kaon_costheta_anti, "true_kaon_costheta_anti/F");
  fEventTree->Branch("true_kaon_phi_anti", &true_kaon_phi_anti, "true_kaon_phi_anti/F");
  fEventTree->Branch("true_kaon_ccmuon_angle_anti", &true_kaon_ccmuon_angle_anti, "true_kaon_ccmuon_angle_anti/F");
  fEventTree->Branch("true_kaon_ccmuon_cosangle_anti", &true_kaon_ccmuon_cosangle_anti, "true_kaon_ccmuon_cosangle_anti/F");
  fEventTree->Branch("true_kaon_end_process_anti", &true_kaon_end_process_anti, "true_kaon_end_process_anti/I");
  fEventTree->Branch("true_kaon_end_ke_anti", &true_kaon_end_ke_anti, "true_kaon_end_ke_anti/F");
  fEventTree->Branch("true_kaon_end_x_anti", &true_kaon_end_x_anti, "true_kaon_end_x_anti/F");
  fEventTree->Branch("true_kaon_end_y_anti", &true_kaon_end_y_anti, "true_kaon_end_y_anti/F");
  fEventTree->Branch("true_kaon_end_z_anti", &true_kaon_end_z_anti, "true_kaon_end_z_anti/F");
  fEventTree->Branch("true_kaon_end_inTPC_anti", &true_kaon_end_inTPC_anti, "true_kaon_end_inTPC_anti/O");
  fEventTree->Branch("true_kaon_end_in5cmTPC_anti", &true_kaon_end_in5cmTPC_anti, "true_kaon_end_in5cmTPC_anti/O");
  fEventTree->Branch("true_kaon_end_inCCInclusiveTPC_anti", &true_kaon_end_inCCInclusiveTPC_anti, "true_kaon_end_inCCInclusiveTPC_anti/O");

  fEventTree->Branch("true_kaon_ndaughters_anti", &true_kaon_ndaughters_anti, "true_kaon_ndaughters_anti/I");
  fEventTree->Branch("true_kaon_ndaughters_decay_anti", &true_kaon_ndaughters_decay_anti, "true_kaon_ndaughters_decay_anti/I");
  fEventTree->Branch("true_kaon_ndaughters_inelastic_anti", &true_kaon_ndaughters_inelastic_anti, "true_kaon_ndaughters_inelastic_anti/I");
  fEventTree->Branch("true_kaon_ndecmup_anti", &true_kaon_ndecmup_anti, "true_kaon_ndecmup_anti/I");
  fEventTree->Branch("true_kaon_ndecpip_anti", &true_kaon_ndecpip_anti, "true_kaon_ndecpip_anti/I");
  fEventTree->Branch("true_kaon_ninekap_anti", &true_kaon_ninekap_anti, "true_kaon_ninekap_anti/I");
  fEventTree->Branch("true_kaon_ninepip_anti", &true_kaon_ninepip_anti, "true_kaon_ninepip_anti/I");
  fEventTree->Branch("true_kaon_ninepro_anti", &true_kaon_ninepro_anti, "true_kaon_ninepro_anti/I");
  fEventTree->Branch("true_kaon_daughter_length_anti", &true_kaon_daughter_length_anti, "true_kaon_daughter_length_anti/F");
  fEventTree->Branch("true_kaon_daughter_p_anti", &true_kaon_daughter_p_anti, "true_kaon_daughter_p_anti/F");
  fEventTree->Branch("true_kaon_daughter_ke_anti", &true_kaon_daughter_ke_anti, "true_kaon_daughter_ke_anti/F");
  fEventTree->Branch("true_kaon_daughter_theta_anti", &true_kaon_daughter_theta_anti, "true_kaon_daughter_theta_anti/F");
  fEventTree->Branch("true_kaon_daughter_costheta_anti", &true_kaon_daughter_costheta_anti, "true_kaon_daughter_costheta_anti/F");
  fEventTree->Branch("true_kaon_daughter_angle_anti", &true_kaon_daughter_angle_anti, "true_kaon_daughter_angle_anti/F");
  fEventTree->Branch("true_kaon_daughter_cosangle_anti", &true_kaon_daughter_cosangle_anti, "true_kaon_daughter_cosangle_anti/F");
  fEventTree->Branch("true_kaon_daughter_pdg_anti", &true_kaon_daughter_pdg_anti, "true_kaon_daughter_pdg_anti/I");
  fEventTree->Branch("true_kaon_daughter_end_x_anti", &true_kaon_daughter_end_x_anti, "true_kaon_daughter_end_x_anti/F");
  fEventTree->Branch("true_kaon_daughter_end_y_anti", &true_kaon_daughter_end_y_anti, "true_kaon_daughter_end_y_anti/F");
  fEventTree->Branch("true_kaon_daughter_end_z_anti", &true_kaon_daughter_end_z_anti, "true_kaon_daughter_end_z_anti/F");
  fEventTree->Branch("true_kaon_daughter_end_inTPC_anti", &true_kaon_daughter_end_inTPC_anti, "true_kaon_daughter_endInTPC_anti/O");
  fEventTree->Branch("true_kaon_daughter_end_in5cmTPC_anti", &true_kaon_daughter_end_in5cmTPC_anti, "true_kaon_daughter_end_in5cmTPC_anti/O");
  fEventTree->Branch("true_kaon_daughter_end_inCCInclusiveTPC_anti", &true_kaon_daughter_end_inCCInclusiveTPC_anti, "true_kaon_daughter_end_inCCInclusiveTPC_anti/O");


  fEventTree->Branch("reco_nu_cc_filter", &reco_nu_cc_filter, "reco_nu_cc_filter/O");

  fEventTree->Branch("reco_nu_vtx_x", &reco_nu_vtx_x, "reco_nu_vtx_x/F");
  fEventTree->Branch("reco_nu_vtx_y", &reco_nu_vtx_y, "reco_nu_vtx_y/F");
  fEventTree->Branch("reco_nu_vtx_z", &reco_nu_vtx_z, "reco_nu_vtx_z/F");
  fEventTree->Branch("reco_nu_vtx_inTPC", &reco_nu_vtx_inTPC, "reco_nu_vtx_inTPC/O");
  fEventTree->Branch("reco_nu_vtx_in5cmTPC", &reco_nu_vtx_in5cmTPC, "reco_nu_vtx_in5cmTPC/O");
  fEventTree->Branch("reco_nu_vtx_inCCInclusiveTPC", &reco_nu_vtx_inCCInclusiveTPC, "reco_nu_vtx_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_nu_ndaughters", &reco_nu_ndaughters, "reco_nu_ndaughters/I"          );
  fEventTree->Branch("reco_nu_cc_nmue", &reco_nu_cc_nmue, "reco_nu_cc_nmue/I"          );

  fEventTree->Branch("reco_ccmu_vtx_x", &reco_ccmu_vtx_x, "reco_ccmu_vtx_x/F");
  fEventTree->Branch("reco_ccmu_vtx_y", &reco_ccmu_vtx_y, "reco_ccmu_vtx_y/F");
  fEventTree->Branch("reco_ccmu_vtx_z", &reco_ccmu_vtx_z, "reco_ccmu_vtx_z/F");
  fEventTree->Branch("reco_ccmu_vtx_inTPC",&reco_ccmu_vtx_inTPC, "reco_ccmu_vtx_inTPC/O");
  fEventTree->Branch("reco_ccmu_vtx_in5cmTPC", &reco_ccmu_vtx_in5cmTPC, "reco_ccmu_vtx_in5cmTPC/O");
  fEventTree->Branch("reco_ccmu_vtx_inCCInclusiveTPC", &reco_ccmu_vtx_inCCInclusiveTPC, "reco_ccmu_vtx_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_ccmu_true_pdg", &reco_ccmu_true_pdg, "reco_ccmu_true_pdg/I");
  fEventTree->Branch("reco_ccmu_true_origin", &reco_ccmu_true_origin, "reco_ccmu_true_origin/I");
  fEventTree->Branch("reco_ccmu_true_primary", &reco_ccmu_true_primary, "reco_ccmu_true_primary/O");
  fEventTree->Branch("reco_ccmu_true_end_inTPC", &reco_ccmu_true_end_inTPC, "reco_ccmu_true_end_inTPC/O");
  fEventTree->Branch("reco_ccmu_true_end_in5cmTPC", &reco_ccmu_true_end_in5cmTPC, "reco_ccmu_true_end_in5cmTPC/O");
  fEventTree->Branch("reco_ccmu_true_end_inCCInclusiveTPC", &reco_ccmu_true_end_inCCInclusiveTPC, "reco_ccmu_true_end_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_ccmu_true_length", &reco_ccmu_true_length, "reco_ccmu_true_length/F");

  fEventTree->Branch("reco_ntracks", &reco_ntracks, "reco_ntracks/I");
  fEventTree->Branch("reco_nshowers", &reco_ntracks, "reco_nshowers/I");
  fEventTree->Branch("nTracks", &nTracks, "nTracks/I");
  fEventTree->Branch("nShowers", &nShowers, "nShowers/I");
  //  fEventTree->Branch("m_reco_track_sliceId", &m_reco_track_sliceId, "m_reco_track_sliceId[20]/F");
  //fEventTree->Branch("m_reco_track_is_nuslice", &m_reco_track_is_nuslice, "m_reco_track_is_nuslice[20]/F");

  fEventTree->Branch("reco_track_start_x", &reco_track_start_x, "reco_track_start_x[20]/F");
  fEventTree->Branch("reco_track_start_y", &reco_track_start_y, "reco_track_start_y[20]/F");
  fEventTree->Branch("reco_track_start_z", &reco_track_start_z, "reco_track_start_z[20]/F");
  fEventTree->Branch("reco_track_end_x", &reco_track_end_x, "reco_track_end_x[20]/F");
  fEventTree->Branch("reco_track_end_y", &reco_track_end_y, "reco_track_end_y[20]/F");
  fEventTree->Branch("reco_track_end_z", &reco_track_end_z, "reco_track_end_z[20]/F");


  fEventTree->Branch("reco_track_distance", &reco_track_distance, "reco_track_distance[20]/F");
  fEventTree->Branch("reco_track_nhits0", &reco_track_nhits0, "reco_track_nhits0[20]/I");
  fEventTree->Branch("reco_track_nhits1", &reco_track_nhits1, "reco_track_nhits1[20]/I");
  fEventTree->Branch("reco_track_nhits2", &reco_track_nhits2, "reco_track_nhits2[20]/I");

  fEventTree->Branch("reco_track_kin0", &reco_track_kin0, "reco_track_kin0[20]/F");
  fEventTree->Branch("reco_track_kin1", &reco_track_kin1, "reco_track_kin1[20]/F");
  fEventTree->Branch("reco_track_kin2", &reco_track_kin2, "reco_track_kin2[20]/F");
  //fEventTree->Branch("reco_track_dEdx" ,&reco_track_dEdx);
  //fEventTree->Branch("reco_track_ResRan" ,&reco_track_ResRan);
  //

  fEventTree->Branch("reco_track_dEdx_pl0", &reco_track_dEdx_pl0,"reco_track_dEdx_pl0[20][2000]");
  fEventTree->Branch("reco_track_ResRan_pl0", &reco_track_ResRan_pl0,"reco_track_ResRan_pl0[20][2000]");

  fEventTree->Branch("reco_track_dEdx_pl1", &reco_track_dEdx_pl1,"reco_track_dEdx_pl1[20][2000]");
  fEventTree->Branch("reco_track_ResRan_pl1", &reco_track_ResRan_pl1,"reco_track_ResRan_pl1[20][2000]");

  fEventTree->Branch("reco_track_dEdx_pl2", &reco_track_dEdx_pl2,"reco_track_dEdx_pl2[20][2000]");
  fEventTree->Branch("reco_track_ResRan_pl2", &reco_track_ResRan_pl2,"reco_track_ResRan_pl2[20][2000]");

  fEventTree->Branch("reco_track_length", &reco_track_length, "reco_track_length[20]/F");
  fEventTree->Branch("reco_track_theta", &reco_track_theta, "reco_track_theta[20]/F");
  fEventTree->Branch("reco_track_phi", &reco_track_phi, "reco_track_phi[20]/F");
  fEventTree->Branch("reco_track_dir", &reco_track_dir, "reco_track_dir[20]/O");

  fEventTree->Branch("reco_track_P_vtx", &reco_track_P_vtx, "reco_track_P_vtx[20]/F");
  fEventTree->Branch("reco_track_P_str", &reco_track_P_str, "reco_track_P_str[20]/F");
  fEventTree->Branch("reco_track_P_end", &reco_track_P_end, "reco_track_P_end[20]/F");

  fEventTree->Branch("reco_track_chi2ka_pl0", &reco_track_chi2ka_pl0, "reco_track_chi2ka_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2pr_pl0", &reco_track_chi2pr_pl0, "reco_track_chi2pr_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2pi_pl0", &reco_track_chi2pi_pl0, "reco_track_chi2pi_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2mu_pl0", &reco_track_chi2mu_pl0, "reco_track_chi2mu_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2ka_pl1", &reco_track_chi2ka_pl1, "reco_track_chi2ka_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2pr_pl1", &reco_track_chi2pr_pl1, "reco_track_chi2pr_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2pi_pl1", &reco_track_chi2pi_pl1, "reco_track_chi2pi_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2mu_pl1", &reco_track_chi2mu_pl1, "reco_track_chi2mu_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2ka_pl2", &reco_track_chi2ka_pl2, "reco_track_chi2ka_pl2[20]/F");
  fEventTree->Branch("reco_track_chi2pr_pl2", &reco_track_chi2pr_pl2, "reco_track_chi2pr_pl2[20]/F");
  fEventTree->Branch("reco_track_chi2pi_pl2", &reco_track_chi2pi_pl2, "reco_track_chi2pi_pl2[20]/F");
  fEventTree->Branch("reco_track_chi2mu_pl2", &reco_track_chi2mu_pl2, "reco_track_chi2mu_pl2[20]/F");

  fEventTree->Branch("reco_track_Bragg_fwd_ka_pl0", &reco_track_Bragg_fwd_ka_pl0, "reco_track_Bragg_fwd_ka_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pr_pl0", &reco_track_Bragg_fwd_pr_pl0, "reco_track_Bragg_fwd_pr_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pi_pl0", &reco_track_Bragg_fwd_pi_pl0, "reco_track_Bragg_fwd_pi_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_mu_pl0", &reco_track_Bragg_fwd_mu_pl0, "reco_track_Bragg_fwd_mu_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_ka_pl1", &reco_track_Bragg_fwd_ka_pl1, "reco_track_Bragg_fwd_ka_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pr_pl1", &reco_track_Bragg_fwd_pr_pl1, "reco_track_Bragg_fwd_pr_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pi_pl1", &reco_track_Bragg_fwd_pi_pl1, "reco_track_Bragg_fwd_pi_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_mu_pl1", &reco_track_Bragg_fwd_mu_pl1, "reco_track_Bragg_fwd_mu_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_ka_pl2", &reco_track_Bragg_fwd_ka_pl2, "reco_track_Bragg_fwd_ka_pl2[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pr_pl2", &reco_track_Bragg_fwd_pr_pl2, "reco_track_Bragg_fwd_pr_pl2[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pi_pl2", &reco_track_Bragg_fwd_pi_pl2, "reco_track_Bragg_fwd_pi_pl2[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_mu_pl2", &reco_track_Bragg_fwd_mu_pl2, "reco_track_Bragg_fwd_mu_pl2[20]/F");


  fEventTree->Branch("reco_track_MIP_pl0", &reco_track_MIP_pl0, "reco_track_MIP_pl0[20]/F");
  fEventTree->Branch("reco_track_MIP_pl1", &reco_track_MIP_pl1, "reco_track_MIP_pl1[20]/F");
  fEventTree->Branch("reco_track_MIP_pl2", &reco_track_MIP_pl2, "reco_track_MIP_pl2[20]/F");


  fEventTree->Branch("reco_track_chi2ka_3pl", &reco_track_chi2ka_3pl, "reco_track_chi2ka_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2pr_3pl", &reco_track_chi2pr_3pl, "reco_track_chi2pr_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2pi_3pl", &reco_track_chi2pi_3pl, "reco_track_chi2pi_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2mu_3pl", &reco_track_chi2mu_3pl, "reco_track_chi2mu_3pl[20]/F");
  fEventTree->Branch("reco_track_likepr_3pl", &reco_track_likepr_3pl, "reco_track_likepr_3pl[20]/F");
  fEventTree->Branch("reco_track_llrpid_3pl", &reco_track_llrpid_3pl, "reco_track_llrpid_3pl[20]/F");
  fEventTree->Branch("reco_track_total_llrpid_3pl", &reco_track_total_llrpid_3pl, "reco_track_total_llrpid_3pl[20]/F");
  fEventTree->Branch("reco_track_llrpid_k_3pl", &reco_track_llrpid_k_3pl, "reco_track_llrpid_k_3pl[20]/F");
  fEventTree->Branch("reco_track_vtx_inTPC", &reco_track_vtx_inTPC, "reco_track_vtx_inTPC[20]/O");
  fEventTree->Branch("reco_track_vtx_in5cmTPC", &reco_track_vtx_in5cmTPC, "reco_track_vtx_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_vtx_inCCInclusiveTPC", &reco_track_vtx_inCCInclusiveTPC, "reco_track_vtx_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_end_inTPC", &reco_track_end_inTPC, "reco_track_end_inTPC[20]/O");
  fEventTree->Branch("reco_track_end_in5cmTPC", &reco_track_end_in5cmTPC, "reco_track_end_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_end_inCCInclusiveTPC", &reco_track_end_inCCInclusiveTPC, "reco_track_end_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_true_pdg", &reco_track_true_pdg, "reco_track_true_pdg[20]/I");
  fEventTree->Branch("reco_track_true_origin", &reco_track_true_origin, "reco_track_true_origin[20]/I");
  fEventTree->Branch("reco_track_true_primary", &reco_track_true_primary, "reco_track_true_primary[20]/O");
  fEventTree->Branch("reco_track_true_end_inTPC", &reco_track_true_end_inTPC, "reco_track_true_end_inTPC[20]/O");
  fEventTree->Branch("reco_track_true_end_in5cmTPC", &reco_track_true_end_in5cmTPC, "reco_track_true_end_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_true_end_inCCInclusiveTPC", &reco_track_true_end_inCCInclusiveTPC, "reco_track_true_end_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_true_length", &reco_track_true_length, "reco_track_true_length[20]/F");


  fEventTree->Branch("reco_track_match_e", &reco_track_match_e, "reco_track_match_e[20][10]/F");
  fEventTree->Branch("reco_track_match_hit", &reco_track_match_hit, "reco_track_match_hit[20][10]/F");
  fEventTree->Branch("reco_track_match_epdg", &reco_track_match_epdg, "reco_track_match_epdg[20][10]/I");
  fEventTree->Branch("reco_track_match_hitpdg", &reco_track_match_hitpdg, "reco_track_match_hitpdg[20][10]/I");

  fEventTree->Branch("reco_track_daughter_match_e", &reco_track_daughter_match_e, "reco_track_daughter_match_e[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_match_hit", &reco_track_daughter_match_hit, "reco_track_daughter_match_hit[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_match_epdg", &reco_track_daughter_match_epdg, "reco_track_daughter_match_epdg[20][20][10]/I");
  fEventTree->Branch("reco_track_daughter_match_hitpdg", &reco_track_daughter_match_hitpdg, "reco_track_daughter_match_hitpdg[20][20][10]/I");

  fEventTree->Branch("reco_track_daughter_shower_match_e", &reco_track_daughter_shower_match_e, "reco_track_daughter_shower_match_e[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_shower_match_hit", &reco_track_daughter_shower_match_hit, "reco_track_daughter_shower_match_hit[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_shower_match_epdg", &reco_track_daughter_shower_match_epdg, "reco_track_daughter_shower_match_epdg[20][20][10]/I");
  fEventTree->Branch("reco_track_daughter_shower_match_hitpdg", &reco_track_daughter_shower_match_hitpdg, "reco_track_daughter_shower_match_hitpdg[20][20][10]/I");


  fEventTree->Branch("reco_track_true_pdg_sh", &reco_track_true_pdg_sh, "reco_track_true_pdg_sh[20]/I");
  fEventTree->Branch("reco_track_true_origin_sh", &reco_track_true_origin_sh, "reco_track_true_origin_sh[20]/I");
  fEventTree->Branch("reco_track_true_primary_sh", &reco_track_true_primary_sh, "reco_track_true_primary_sh[20]/O");
  fEventTree->Branch("reco_track_true_end_inTPC_sh", &reco_track_true_end_inTPC_sh, "reco_track_true_end_inTPC_sh[20]/O");
  fEventTree->Branch("reco_track_true_end_in5cmTPC_sh", &reco_track_true_end_in5cmTPC_sh, "reco_track_true_end_in5cmTPC_sh[20]/O");
  fEventTree->Branch("reco_track_true_end_inCCInclusiveTPC_sh", &reco_track_true_end_inCCInclusiveTPC_sh, "reco_track_true_end_inCCInclusiveTPC[20]_sh/O");
  fEventTree->Branch("reco_track_true_length_sh", &reco_track_true_length_sh, "reco_track_true_length_sh[20]/F");

  fEventTree->Branch("reco_track_daughter_true_pdg_sh", &reco_track_daughter_true_pdg_sh, "reco_track_daughter_true_pdg_sh[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_origin_sh", &reco_track_daughter_true_origin_sh, "reco_track_daughter_true_origin_sh[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_primary_sh", &reco_track_daughter_true_primary_sh, "reco_track_daughter_true_primary[20][20]_sh/O");
  fEventTree->Branch("reco_track_daughter_true_end_inTPC_sh", &reco_track_daughter_true_end_inTPC_sh, "reco_track_daughter_true_end_inTPC[20][20]_sh/O");
  fEventTree->Branch("reco_track_daughter_true_end_in5cmTPC_sh", &reco_track_daughter_true_end_in5cmTPC_sh, "reco_track_daughter_true_end_in5cmTPC[20][20]_sh/O");
  fEventTree->Branch("reco_track_daughter_true_end_inCCInclusiveTPC_sh", &reco_track_daughter_true_end_inCCInclusiveTPC_sh, "reco_track_daughter_true_end_inCCInclusiveTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_length_sh", &reco_track_daughter_true_length_sh, "reco_track_daughter_true_length_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_true_mother_sh", &reco_track_daughter_true_mother_sh, "reco_track_daughter_true_mother_sh[20][20]/I");

  fEventTree->Branch("reco_track_daughter_start_x", &reco_track_daughter_start_x, "reco_track_daughter_start_x[20][20]/F"); 
  fEventTree->Branch("reco_track_daughter_start_y", &reco_track_daughter_start_y, "reco_track_daughter_start_y[20][20]/F"); 
  fEventTree->Branch("reco_track_daughter_start_z", &reco_track_daughter_start_z, "reco_track_daughter_start_z[20][20]/F"); 
  fEventTree->Branch("reco_track_daughter_end_x", &reco_track_daughter_end_x, "reco_track_daughter_end_x[20][20]/F");
  fEventTree->Branch("reco_track_daughter_end_y", &reco_track_daughter_end_y, "reco_track_daughter_end_y[20][20]/F");
  fEventTree->Branch("reco_track_daughter_end_z", &reco_track_daughter_end_z, "reco_track_daughter_end_z[20][20]/F");

  fEventTree->Branch("reco_shower_ndaughters", &reco_shower_ndaughters, "reco_shower_ndaughters[20]/I");
  fEventTree->Branch("reco_track_daughter_distance_sh", &reco_track_daughter_distance_sh, "reco_track_daughter_distance_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_distance_sh", &reco_track_daughter_vtx_distance_sh, "reco_track_daughter_vtx_distance_sh[20][20]/F");
  fEventTree->Branch("reco_angle_track_daughter_sh", &reco_angle_track_daughter_sh, "reco_angle_track_daughter_sh[20][20]/F");
  fEventTree->Branch("reco_angle_daughter_track_daughter_sh", &reco_angle_daughter_track_daughter_sh, "reco_angle_daughter_track_daughter_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_length_sh", &reco_track_daughter_length_sh, "reco_track_daughter_length_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_theta_sh", &reco_track_daughter_theta_sh, "reco_track_daughter_theta_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_phi_sh", &reco_track_daughter_phi_sh, "reco_track_daughter_phi_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_open_angle_sh", &reco_track_daughter_open_angle_sh, "reco_track_daughter_open_angle_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_dedx_pl0_sh", &reco_track_daughter_dedx_pl0_sh, "reco_track_daughter_dedx_pl0_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_dedx_pl1_sh", &reco_track_daughter_dedx_pl1_sh, "reco_track_daughter_dedx_pl1_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_dedx_pl2_sh", &reco_track_daughter_dedx_pl2_sh, "reco_track_daughter_dedx_pl2_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_energy_pl0_sh", &reco_track_daughter_energy_pl0_sh, "reco_track_daughter_energy_pl0_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_energy_pl1_sh", &reco_track_daughter_energy_pl1_sh, "reco_track_daughter_energy_pl1_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_energy_pl2_sh", &reco_track_daughter_energy_pl2_sh, "reco_track_daughter_energy_pl2_sh[20][20]/F");

  fEventTree->Branch("reco_track_daughter_vtx_inTPC_sh", &reco_track_daughter_vtx_inTPC_sh, "reco_track_daughter_vtx_inTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_in5cmTPC_sh", &reco_track_daughter_vtx_in5cmTPC_sh, "reco_track_daughter_vtx_in5cmTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_inCCInclusiveTPC_sh", &reco_track_daughter_vtx_inCCInclusiveTPC_sh, "reco_track_daughter_vtx_inCCInclusiveTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inTPC_sh", &reco_track_daughter_end_inTPC_sh, "reco_track_daughter_end_inTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_in5cmTPC_sh", &reco_track_daughter_end_in5cmTPC_sh, "reco_track_daughter_end_in5cmTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inCCInclusiveTPC_sh", &reco_track_daughter_end_inCCInclusiveTPC_sh, "reco_track_daughter_end_inCCInclusiveTPC_sh[20][20]/O");

  fEventTree->Branch("reco_track_ndaughters", &reco_track_ndaughters, "reco_track_ndaughters[20]/I");
  fEventTree->Branch("reco_track_daughter_distance", &reco_track_daughter_distance, "reco_track_daughter_distance[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_distance", &reco_track_daughter_vtx_distance, "reco_track_daughter_vtx_distance[20][20]/F");
  fEventTree->Branch("reco_angle_track_daughter", &reco_angle_track_daughter, "reco_angle_track_daughter[20][20]/F");
  fEventTree->Branch("reco_track_daughter_nhits0", &reco_track_daughter_nhits0, "reco_track_daughter_nhits0[20][20]/I");
  fEventTree->Branch("reco_track_daughter_nhits1", &reco_track_daughter_nhits1, "reco_track_daughter_nhits1[20][20]/I");
  fEventTree->Branch("reco_track_daughter_nhits2", &reco_track_daughter_nhits2, "reco_track_daughter_nhits2[20][20]/I");

  //fEventTree->Branch("reco_track_daughter_dEdx", &reco_track_daughter_dEdx);
  //fEventTree->Branch("reco_track_daughter_ResRan", &reco_track_daughter_ResRan);
  //
  fEventTree->Branch("reco_track_daughter_dEdx_pl0", &reco_track_daughter_dEdx_pl0, "reco_track_daughter_dEdx_pl0[20][20][2000]/F");
  fEventTree->Branch("reco_track_daughter_ResRan_pl0", &reco_track_daughter_ResRan_pl0, "reco_track_daughter_ResRan_pl0[20][20][2000]/F");

  fEventTree->Branch("reco_track_daughter_dEdx_pl1", &reco_track_daughter_dEdx_pl1, "reco_track_daughter_dEdx_pl1[20][20][2000]/F");
  fEventTree->Branch("reco_track_daughter_ResRan_pl1", &reco_track_daughter_ResRan_pl1, "reco_track_daughter_ResRan_pl1[20][20][2000]/F");

  fEventTree->Branch("reco_track_daughter_dEdx_pl2", &reco_track_daughter_dEdx_pl2, "reco_track_daughter_dEdx_pl2[20][20][2000]/F");
  fEventTree->Branch("reco_track_daughter_ResRan_pl2", &reco_track_daughter_ResRan_pl2, "reco_track_daughter_ResRan_pl2[20][20][2000]/F");

  fEventTree->Branch("reco_track_daughter_length", &reco_track_daughter_length, "reco_track_daughter_length[20][20]/F");
  fEventTree->Branch("reco_track_daughter_theta", &reco_track_daughter_theta, "reco_track_daughter_theta[20][20]/F");
  fEventTree->Branch("reco_track_daughter_phi", &reco_track_daughter_phi, "reco_track_daughter_phi[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl0", &reco_track_daughter_chi2ka_pl0, "reco_track_daughter_chi2ka_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl0", &reco_track_daughter_chi2pr_pl0, "reco_track_daughter_chi2pr_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl0", &reco_track_daughter_chi2pi_pl0, "reco_track_daughter_chi2pi_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl0", &reco_track_daughter_chi2mu_pl0, "reco_track_daughter_chi2mu_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl1", &reco_track_daughter_chi2ka_pl1, "reco_track_daughter_chi2ka_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl1", &reco_track_daughter_chi2pr_pl1, "reco_track_daughter_chi2pr_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl1", &reco_track_daughter_chi2pi_pl1, "reco_track_daughter_chi2pi_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl1", &reco_track_daughter_chi2mu_pl1, "reco_track_daughter_chi2mu_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl2", &reco_track_daughter_chi2ka_pl2, "reco_track_daughter_chi2ka_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl2", &reco_track_daughter_chi2pr_pl2, "reco_track_daughter_chi2pr_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl2", &reco_track_daughter_chi2pi_pl2, "reco_track_daughter_chi2pi_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl2", &reco_track_daughter_chi2mu_pl2, "reco_track_daughter_chi2mu_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_3pl", &reco_track_daughter_chi2ka_3pl, "reco_track_daughter_chi2ka_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_3pl", &reco_track_daughter_chi2pr_3pl, "reco_track_daughter_chi2pr_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_3pl", &reco_track_daughter_chi2pi_3pl, "reco_track_daughter_chi2pi_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_3pl", &reco_track_daughter_chi2mu_3pl, "reco_track_daughter_chi2mu_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_likepr_3pl", &reco_track_daughter_likepr_3pl, "reco_track_daughter_likepr_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_3pl", &reco_track_daughter_llrpid_3pl, "reco_track_daughter_llrpid_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_k_3pl", &reco_track_daughter_llrpid_k_3pl, "reco_track_daughter_llrpid_k_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_inTPC", &reco_track_daughter_vtx_inTPC, "reco_track_daughter_vtx_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_in5cmTPC", &reco_track_daughter_vtx_in5cmTPC, "reco_track_daughter_vtx_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_inCCInclusiveTPC", &reco_track_daughter_vtx_inCCInclusiveTPC, "reco_track_daughter_vtx_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inTPC", &reco_track_daughter_end_inTPC, "reco_track_daughter_end_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_in5cmTPC", &reco_track_daughter_end_in5cmTPC, "reco_track_daughter_end_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inCCInclusiveTPC", &reco_track_daughter_end_inCCInclusiveTPC, "reco_track_daughter_end_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_pdg", &reco_track_daughter_true_pdg, "reco_track_daughter_true_pdg[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_origin", &reco_track_daughter_true_origin, "reco_track_daughter_true_origin[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_primary", &reco_track_daughter_true_primary, "reco_track_daughter_true_primary[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_inTPC", &reco_track_daughter_true_end_inTPC, "reco_track_daughter_true_end_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_in5cmTPC", &reco_track_daughter_true_end_in5cmTPC, "reco_track_daughter_true_end_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_inCCInclusiveTPC", &reco_track_daughter_true_end_inCCInclusiveTPC, "reco_track_daughter_true_end_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_length", &reco_track_daughter_true_length, "reco_track_daughter_true_length[20][20]/F");
  fEventTree->Branch("reco_track_daughter_true_mother", &reco_track_daughter_true_mother, "reco_track_daughter_true_mother[20][20]/I");


  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl0", &reco_track_daughter_Bragg_fwd_ka_pl0, "reco_track_daughter_Bragg_fwd_ka_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl0", &reco_track_daughter_Bragg_fwd_pr_pl0, "reco_track_daughter_Bragg_fwd_pr_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl0", &reco_track_daughter_Bragg_fwd_pi_pl0, "reco_track_daughter_Bragg_fwd_pi_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl0", &reco_track_daughter_Bragg_fwd_mu_pl0, "reco_track_daughter_Bragg_fwd_mu_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl1", &reco_track_daughter_Bragg_fwd_ka_pl1, "reco_track_daughter_Bragg_fwd_ka_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl1", &reco_track_daughter_Bragg_fwd_pr_pl1, "reco_track_daughter_Bragg_fwd_pr_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl1", &reco_track_daughter_Bragg_fwd_pi_pl1, "reco_track_daughter_Bragg_fwd_pi_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl1", &reco_track_daughter_Bragg_fwd_mu_pl1, "reco_track_daughter_Bragg_fwd_mu_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl2", &reco_track_daughter_Bragg_fwd_ka_pl2, "reco_track_daughter_Bragg_fwd_ka_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl2", &reco_track_daughter_Bragg_fwd_pr_pl2, "reco_track_daughter_Bragg_fwd_pr_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl2", &reco_track_daughter_Bragg_fwd_pi_pl2, "reco_track_daughter_Bragg_fwd_pi_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl2", &reco_track_daughter_Bragg_fwd_mu_pl2, "reco_track_daughter_Bragg_fwd_mu_pl2[20][20]/F");


  fEventTree->Branch("reco_track_daughter_MIP_pl0", &reco_track_daughter_MIP_pl0, "reco_track_daughter_MIP_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_MIP_pl1", &reco_track_daughter_MIP_pl1, "reco_track_daughter_MIP_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_MIP_pl2", &reco_track_daughter_MIP_pl2, "reco_track_daughter_MIP_pl2[20][20]/F");



  fEventTree->Branch("k_can_trkid", &k_can_trkid,"k_can_trkid/I");
  fEventTree->Branch("mu_can_trkid", &mu_can_trkid,"mu_can_trkid/I");
  fEventTree->Branch("k_mu_can_dis", &k_mu_can_dis,"k_mu_can_dis/F");
  fEventTree->Branch("k_mu_open_angle", &k_mu_open_angle,"k_mu_open_angle/F");
  fEventTree->Branch("k_vtx_dis", &k_vtx_dis,"k_vtx_dis/F");
  //fEventTree->Branch("k_geant_ID", &k_geant_ID,"k_geant_ID/I");
  //fEventTree->Branch("k_origin", &k_origin,"k_origin/I");
  //fEventTree->Branch("k_pdg", &k_pdg,"k_pdg/I");
  //fEventTree->Branch("k_isPri", &k_isPri,"k_isPri/I");
  //fEventTree->Branch("k_endE", &k_endE,"k_endE/F");
  //fEventTree->Branch("k_ness", &k_ness,"k_ness/I");
  //fEventTree->Branch("kaon_vtx_dis", &kaon_vtx_dis,"kaon_vtx_dis/F");
  //fEventTree->Branch("k_plen", &k_plen,"k_plen/F");
  //fEventTree->Branch("k_phi", &k_phi,"k_phi/F");
  //fEventTree->Branch("k_theta", &k_theta,"k_theta/F");
  //fEventTree->Branch("k_in_5_TPC", &k_in_5_TPC,"k_in_5_TPC/I");
  //fEventTree->Branch("k_in_CC_TPC", &k_in_CC_TPC,"k_in_CC_TPC/I");
  //fEventTree->Branch("k_hit", &k_hit,"k_hit/I");
  //fEventTree->Branch("k_range", &k_range,"k_range/F");
  //fEventTree->Branch("k_KE", &k_KE,"k_KE/F");
  //fEventTree->Branch("k_large_dedx", &k_large_dedx,"k_large_dedx/F");
  //fEventTree->Branch("k_small_dedx", &k_small_dedx,"k_small_dedx/F");
  //fEventTree->Branch("k_chi_p", &k_chi_p,"k_chi_p/F");
  //fEventTree->Branch("k_chi_k", &k_chi_k,"k_chi_k/F");
  //fEventTree->Branch("k_chi_pi", &k_chi_pi,"k_chi_pi/F");
  //fEventTree->Branch("k_chi_mu", &k_chi_mu,"k_chi_mu/F");
  //fEventTree->Branch("k_p_max", &k_p_max,"k_p_max/F");
  //fEventTree->Branch("k_k_max", &k_k_max,"k_k_max/F");
  //fEventTree->Branch("k_pi_max", &k_pi_max,"k_pi_max/F");
  //fEventTree->Branch("k_mu_max", &k_mu_max,"k_mu_max/F");
  //fEventTree->Branch("k_mip_max", &k_mip_max,"k_mip_max/F");
  //fEventTree->Branch("k_L1_ratio", &k_L1_ratio,"k_L1_ratio/F");
  //fEventTree->Branch("k_LL1", &k_LL1,"k_LL1/F");
  //fEventTree->Branch("k_L2_ratio", &k_L2_ratio,"k_L2_ratio/F");
  //fEventTree->Branch("k_LL2", &k_LL2,"k_LL2/F");
  //fEventTree->Branch("k_Lp_ratio", &k_Lp_ratio,"k_Lp_ratio/F");
  //fEventTree->Branch("k_LLp", &k_LLp,"k_LLp/F");
  //fEventTree->Branch("k_Lk_ratio", &k_Lk_ratio,"k_Lk_ratio/F");
  //fEventTree->Branch("k_LLk", &k_LLk,"k_LLk/F");
  //fEventTree->Branch("k_pida_mean", &k_pida_mean,"k_pida_mean/F");
  //fEventTree->Branch("k_pida_med", &k_pida_med,"k_pida_med/F");
  //fEventTree->Branch("k_kde", &k_kde,"k_kde/F");
  //fEventTree->Branch("k_trm_dqdx", &k_trm_dqdx,"k_trm_dqdx/F");
  //fEventTree->Branch("k_trm_dedx", &k_trm_dedx,"k_trm_dedx/F");
  //fEventTree->Branch("k_dedx",k_dedx,"k_dedx[3][3000]/F");
  //fEventTree->Branch("k_rr",k_rr,"k_rr[3][3000]/F");
  //fEventTree->Branch("mu_pdg", &mu_pdg,"mu_pdg/I");
  //fEventTree->Branch("mu_isDec", &mu_isDec,"mu_isDec/I");
  //fEventTree->Branch("mu_origin", &mu_origin,"mu_origin/I");
  //fEventTree->Branch("mu_k_is_Mother", &mu_k_is_Mother,"mu_k_is_Mother/I");
  //fEventTree->Branch("mu_mom_k_inelas", &mu_mom_k_inelas,"mu_mom_k_inelas/I");
  //fEventTree->Branch("mu_ness", &mu_ness,"mu_ness/I");
  //fEventTree->Branch("mu_plen", &mu_plen,"mu_plen/F");
  //fEventTree->Branch("mu_phi", &mu_phi,"mu_phi/F");
  //fEventTree->Branch("mu_theta", &mu_theta,"mu_theta/F");
  //fEventTree->Branch("mu_in_5_TPC", &mu_in_5_TPC,"mu_in_5_TPC/I");
  //fEventTree->Branch("mu_in_CC_TPC", &mu_in_CC_TPC,"mu_in_CC_TPC/I");
  //fEventTree->Branch("mu_KE", &mu_KE,"mu_KE/F");
  //fEventTree->Branch("mu_hit", &mu_hit,"mu_hit/I");
  //fEventTree->Branch("mu_range", &mu_range,"mu_range/F");
  //fEventTree->Branch("mu_large_dedx", &mu_large_dedx,"mu_large_dedx/F");
  //fEventTree->Branch("mu_small_dedx", &mu_small_dedx,"mu_small_dedx/F");
  //fEventTree->Branch("mu_chi_p", &mu_chi_p,"mu_chi_p/F");
  //fEventTree->Branch("mu_chi_k", &mu_chi_k,"mu_chi_k/F");
  //fEventTree->Branch("mu_chi_pi", &mu_chi_pi,"mu_chi_pi/F");
  //fEventTree->Branch("mu_chi_mu", &mu_chi_mu,"mu_chi_mu/F");
  //fEventTree->Branch("mu_p_max", &mu_p_max,"mu_p_max/F");
  //fEventTree->Branch("mu_k_max", &mu_k_max,"mu_k_max/F");
  //fEventTree->Branch("mu_pi_max", &mu_pi_max,"mu_pi_max/F");
  //fEventTree->Branch("mu_mu_max", &mu_mu_max,"mu_mu_max/F");
  //fEventTree->Branch("mu_mip_max", &mu_mip_max,"mu_mip_max/F");
  //fEventTree->Branch("mu_L1_ratio", &mu_L1_ratio,"mu_L1_ratio/F");
  //fEventTree->Branch("mu_LL1", &mu_LL1,"mu_LL1/F");
  //fEventTree->Branch("mu_L2_ratio", &mu_L2_ratio,"mu_L2_ratio/F");
  //fEventTree->Branch("mu_LL2", &mu_LL2,"mu_LL2/F");
  //fEventTree->Branch("mu_Lp_ratio", &mu_Lp_ratio,"mu_Lp_ratio/F");
  //fEventTree->Branch("mu_LLp", &mu_LLp,"mu_LLp/F");
  //fEventTree->Branch("mu_Lk_ratio", &mu_Lk_ratio,"mu_Lk_ratio/F");
  //fEventTree->Branch("mu_LLk", &mu_LLk,"mu_LLk/F");
  //fEventTree->Branch("mu_pida_mean", &mu_pida_mean,"mu_pida_mean/F");
  //fEventTree->Branch("mu_pida_med", &mu_pida_med,"mu_pida_med/F");
  //fEventTree->Branch("mu_kde", &mu_kde,"mu_kde/F");
  //fEventTree->Branch("mu_trm_dqdx", &mu_trm_dqdx,"mu_trm_dqdx/F");
  //fEventTree->Branch("mu_trm_dedx", &mu_trm_dedx,"mu_trm_dedx/F");
  //fEventTree->Branch("mu_dedx",mu_dedx,"mu_dedx[3][3000]/F");
  //fEventTree->Branch("mu_rr",mu_rr,"mu_rr[3][3000]/F");
  //fEventTree->Branch("mu_mom_process",&mu_mom_process);
  //fEventTree->Branch("cc_mu_trkid", &cc_mu_trkid,"cc_mu_trkid/I");
  //fEventTree->Branch("cc_mu_tlen", &cc_mu_tlen,"cc_mu_tlen/F");
  //fEventTree->Branch("cc_mu_phi", &cc_mu_phi,"cc_mu_phi/F");
  //fEventTree->Branch("cc_mu_theta", &cc_mu_theta,"cc_mu_theta/F");
  //fEventTree->Branch("cc_mu_range", &cc_mu_range,"cc_mu_range/F");
  //fEventTree->Branch("cc_mu_KE", &cc_mu_KE,"cc_mu_KE/F");
  //fEventTree->Branch("cc_mu_hit", &cc_mu_hit,"cc_mu_hit/I");
  //fEventTree->Branch("cc_mu_large_dedx", &cc_mu_large_dedx,"cc_mu_large_dedx/F");
  //fEventTree->Branch("cc_mu_small_dedx", &cc_mu_small_dedx,"cc_mu_small_dedx/F");
  //fEventTree->Branch("cc_dis_vtx", &cc_dis_vtx,"cc_dis_vtx/F");
  //fEventTree->Branch("cc_mu_pdg", &cc_mu_pdg,"cc_mu_pdg/I");
  //fEventTree->Branch("cc_mu_chi_p", &cc_mu_chi_p,"cc_mu_chi_p/F");
  //fEventTree->Branch("cc_mu_chi_k", &cc_mu_chi_k,"cc_mu_chi_k/F");
  //fEventTree->Branch("cc_mu_chi_pi", &cc_mu_chi_pi,"cc_mu_chi_pi/F");
  //fEventTree->Branch("cc_mu_chi_mu", &cc_mu_chi_mu,"cc_mu_chi_mu/F");
  //fEventTree->Branch("cc_mu_p_max", &cc_mu_p_max,"cc_mu_p_max/F");
  //fEventTree->Branch("cc_mu_k_max", &cc_mu_k_max,"cc_mu_k_max/F");
  //fEventTree->Branch("cc_mu_pi_max", &cc_mu_pi_max,"cc_mu_pi_max/F");
  //fEventTree->Branch("cc_mu_mu_max", &cc_mu_mu_max,"cc_mu_mu_max/F");
  //fEventTree->Branch("cc_mu_mip_max", &cc_mu_mip_max,"cc_mu_mip_max/F");
  //fEventTree->Branch("cc_mu_L1_ratio", &cc_mu_L1_ratio,"cc_mu_L1_ratio/F");
  //fEventTree->Branch("cc_mu_LL1", &cc_mu_LL1,"cc_mu_LL1/F");
  //fEventTree->Branch("cc_mu_L2_ratio", &cc_mu_L2_ratio,"cc_mu_L2_ratio/F");
  //fEventTree->Branch("cc_mu_LL2", &cc_mu_LL2,"cc_mu_LL2/F");
  //fEventTree->Branch("cc_mu_Lp_ratio", &cc_mu_Lp_ratio,"cc_mu_Lp_ratio/F");
  //fEventTree->Branch("cc_mu_LLp", &cc_mu_LLp,"cc_mu_LLp/F");
  //fEventTree->Branch("cc_mu_Lk_ratio", &cc_mu_Lk_ratio,"cc_mu_Lk_ratio/F");
  //fEventTree->Branch("cc_mu_LLk", &cc_mu_LLk,"cc_mu_LLk/F");
  //fEventTree->Branch("cc_mu_pida_mean", &cc_mu_pida_mean,"cc_mu_pida_mean/F");
  //fEventTree->Branch("cc_mu_pida_med", &cc_mu_pida_med,"cc_mu_pida_med/F");
  //fEventTree->Branch("cc_mu_kde", &cc_mu_kde,"cc_mu_kde/F");
  //fEventTree->Branch("cc_mu_trm_dqdx", &cc_mu_trm_dqdx,"cc_mu_trm_dqdx/F");
  //fEventTree->Branch("cc_mu_trm_dedx", &cc_mu_trm_dedx,"cc_mu_trm_dedx/F");
  //fEventTree->Branch("longest_trkid", &longest_trkid,"longest_trkid/I");
  //fEventTree->Branch("longest_trklen", &longest_trklen,"longest_trklen/F");
  //fEventTree->Branch("Tr_pri_mu_pdg", &Tr_pri_mu_pdg,"Tr_pri_mu_pdg/I");
  //fEventTree->Branch("pri_Mu_is", &pri_Mu_is,"pri_Mu_is/I");
  //fEventTree->Branch("Tr_pri_st_k_is", &Tr_pri_st_k_is,"Tr_pri_st_k_is/I");
  //fEventTree->Branch("Tr_K_Inelas", &Tr_K_Inelas,"Tr_K_Inelas/I");
  //fEventTree->Branch("Tr_k_plen", &Tr_k_plen,"Tr_k_plen/F");
  //fEventTree->Branch("Tr_k_theta", &Tr_k_theta,"Tr_k_theta/F");
  //fEventTree->Branch("Tr_k_phi", &Tr_k_phi,"Tr_k_phi/F");
  //fEventTree->Branch("Tr_dec_mu_is", &Tr_dec_mu_is,"Tr_dec_mu_is/I");
  //fEventTree->Branch("Tr_dec_mu_pi_pdg", &Tr_dec_mu_pi_pdg,"Tr_dec_mu_pi_pdg/I");
  //fEventTree->Branch("Tr_mu_plen", &Tr_mu_plen,"Tr_mu_plen/F");
  //fEventTree->Branch("Tr_k_endE", &Tr_k_endE,"Tr_k_endE/F");
  //fEventTree->Branch("Tr_mu_theta", &Tr_mu_theta,"Tr_mu_theta/F");
  //fEventTree->Branch("Tr_mu_phi", &Tr_mu_phi,"Tr_mu_phi/F");
  //fEventTree->Branch("Tr_k_inTPC", &Tr_k_inTPC,"Tr_k_inTPC/I");
  //fEventTree->Branch("Tr_mu_inTPC", &Tr_mu_inTPC,"Tr_mu_inTPC/I");
  //fEventTree->Branch("Tr_k_in_5_TPC", &Tr_k_in_5_TPC,"Tr_k_in_5_TPC/I");
  //fEventTree->Branch("Tr_k_in_CC_TPC", &Tr_k_in_CC_TPC,"Tr_k_in_CC_TPC/I");
  //fEventTree->Branch("Tr_mu_in_5_TPC", &Tr_mu_in_5_TPC,"Tr_mu_in_5_TPC/I");
  //fEventTree->Branch("Tr_mu_in_CC_TPC", &Tr_mu_in_CC_TPC,"Tr_mu_in_CC_TPC/I");
  //fEventTree->Branch("Tr_kmu_open_ang", &Tr_kmu_open_ang,"Tr_kmu_open_ang/F");
  //fEventTree->Branch("vtx_5cm_mult", &vtx_5cm_mult,"vtx_5cm_mult/I");
  //fEventTree->Branch("k_start_dedx", &k_start_dedx,"k_start_dedx/F");
  //fEventTree->Branch("k_end_dedx", &k_end_dedx,"k_end_dedx/F");
  //fEventTree->Branch("mu_start_dedx", &mu_start_dedx,"mu_start_dedx/F");
  //fEventTree->Branch("mu_end_dedx", &mu_end_dedx,"mu_end_dedx/F");
  fEventTree->Branch("cut_1", &cut_1,"cut_1/I");
  fEventTree->Branch("cut_2", &cut_2,"cut_2/I");
  fEventTree->Branch("cut_3", &cut_3,"cut_3/I");
  fEventTree->Branch("cut_4", &cut_4,"cut_4/I");
  fEventTree->Branch("cut_5", &cut_5,"cut_5/I");
  fEventTree->Branch("cut_6", &cut_6,"cut_6/I");
  fEventTree->Branch("cut_7", &cut_7,"cut_7/I");
  fEventTree->Branch("cut_8", &cut_8,"cut_8/I");
  fEventTree->Branch("cut_9", &cut_9,"cut_9/I");
  fEventTree->Branch("cut_10", &cut_10,"cut_10/I");
  fEventTree->Branch("cut_11", &cut_11,"cut_11/I");
  fEventTree->Branch("cut_12", &cut_12,"cut_12/I");
  fEventTree->Branch("cut_13", &cut_13,"cut_13/I");
  fEventTree->Branch("PFP_have_nuslice",& PFP_have_nuslice);

  //fEventTree->Branch("kinelas_has_traks", &kinelas_has_traks,"kinelas_has_traks/I");
  //fEventTree->Branch("kinelas_reco_trkID", &kinelas_reco_trkID,"kinelas_reco_trkID/I");
  //fEventTree->Branch("kinelas_tlen", &kinelas_tlen,"kinelas_tlen/F");
  //fEventTree->Branch("True_kinelas_KE", &True_kinelas_KE,"True_kinelas_KE/F");
  //fEventTree->Branch("True_kinelas_tlen", &True_kinelas_tlen,"True_kinelas_tlen/F");

  fSubrunTree = tfs->make<TTree>("subruns", "SubRun Tree");
  fSubrunTree->Branch("run", &m_run, "run/i");
  fSubrunTree->Branch("subRun", &m_subrun, "subRun/i");
  //if (!m_isData)
  fSubrunTree->Branch("pot", &m_pot, "pot/F");
    
  vertex_tree = tfs->make<TTree>("vertex_tree", "vertex_tree");
  pot_tree = tfs->make<TTree>("pot_tree", "pot_tree");
  eventweight_tree = tfs->make<TTree>("eventweight_tree", "eventweight_tree");
  ncdelta_slice_tree = tfs->make<TTree>("ncdelta_slice_tree", "ncdelta_slice_tree");

  this->CreateTrackBranches();
  this->CreateShowerBranches();
  this->CreateMCTruthBranches();
  this->CreateSliceBranches();

}

void CCKaonAnalyzerRebuild::analyze( const art::Event& evt){
  reset();  
  run=evt.run();
  subrun=evt.subRun();
  event=evt.id().event(); 
  cout << "Run " << run;
  cout << " Subrun " << subrun;
  cout << " Event " << event << endl;

  art::Handle< std::vector<simb::MCTruth> > mctruths;
  std::vector< art::Ptr<simb::MCParticle> > ptList;


  //-----------scan--------------
  
    n_prip=0;
    n_prip_k_dau=0;
    Int_t prim_k_id = -999;
    prip = {-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999};
    prip_k_dau = {-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999};
    IsKaon = false;
    IsSingleKaon = false;
    IsAssociatedKaon = false;
    IsMuBR = false;
    IsPiBR = false;
    flag_prim_k = true;
    

    // get MCTruth     
    evt.getByLabel("generator", mctruths);
    if (mctruths->size()!=1) {
      //std::cout << "Number of MCTruths objects in event " << mctruths->size() << std::endl;
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
    true_nu_vtx_inTPC = isInsideVolume("TPC",mctruth.GetNeutrino().Nu().Position().Vect());
    true_nu_vtx_in5cmTPC = isInsideVolume("5cmTPC",mctruth.GetNeutrino().Nu().Position().Vect());
    true_nu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",mctruth.GetNeutrino().Nu().Position().Vect());
    

    /*
    //sim::ParticleList::const_iterator itPart = ptList.begin(),
    for(size_t iPart = 0; iPart < ptList.size(); ++iPart){
      const simb::MCParticle* pPart = (itPart++)->second;
      if (!pPart) {
	throw art::Exception(art::errors::LogicError)
	  << "GEANT particle #" << iPart << " returned a null pointer";
      }
      TrackIDtoIndex.emplace(TrackID, iPart);
      gpdg.push_back(pPart->PdgCode());
      gmother.push_back(pPart->Mother());
      
    }
    */

    // print true information
    for (auto const& pPart : ptList) {

      //if (pPart->Process()!="Primary" && pPart->E()-pPart->Mass()<0.001) continue;
      
      cout << "trackId " << pPart->TrackId();
      cout << ", mother " << pPart->Mother();
      cout << ", pdg " << pPart->PdgCode();
      cout << ", process " << pPart->Process() << " - " << pPart->EndProcess();
      cout << ", ke " << pPart->E()-pPart->Mass() << " - " << pPart->EndE()-pPart->Mass();
      cout << ", P " << pPart->Momentum().Vect().Mag();
      cout << ", Eng " << pPart->E() << ", EndE " << pPart->EndE();
      cout << ", X " << pPart->Vx() << " - " << pPart->EndX();
      cout << ", Y " << pPart->Vy() << " - " << pPart->EndY();
      cout << ", Z " << pPart->Vz() << " - " << pPart->EndZ();
      cout << ", EndPointx " << pPart->EndPosition()[0];
      cout << ", EndPointy " << pPart->EndPosition()[1];
      cout << ", EndPointz " << pPart->EndPosition()[2];
      cout << ", delta Z " << pPart->EndZ()-pPart->Vz();
      cout << ", ndaughters " << pPart->NumberDaughters();
      cout << endl;
      

      /*
      cout << "trackId " << pPart->TrackId();
      cout << ", mother " << pPart->Mother();
      cout << ", pdg " << pPart->PdgCode();
      cout << ", process " << pPart->Process() << " - " << pPart->EndProcess();
      cout << ", X " << pPart->Vz() << " - " << pPart->EndZ();
      cout << ", EndPointx " << pPart->EndPosition()[2];
      cout << ", Y " << pPart->Vx() << " - " << pPart->EndX();
      cout << ", EndPointy " << pPart->EndPosition()[0];
      cout << endl;
      */
      
      //cout << n_prip << endl;
      if(n_prip<20 && pPart->Mother()==0){//primary interaction

	prip.at(n_prip) = pPart->PdgCode();

	if( (flag_prim_k ==true)  && (prip.at(n_prip) == 321)){//kaon
	  prim_k_id = pPart->TrackId();
	  flag_prim_k = false;
	}
	++n_prip;

      }

      else continue;
      //++n_prip;  

    }
    n_prip=0;



    for (auto const& pPart : ptList) {  

      if(n_prip_k_dau<20 && pPart->Mother()==prim_k_id){
	prip_k_dau.at(n_prip_k_dau) = pPart->PdgCode();
	++n_prip_k_dau;
      }

      else continue;
      //++n_prip_k_dau;      

    }
    n_prip_k_dau=0;
    flag_prim_k = true;


    /*
    if(IsKaon == false){

      //if(std::find(prip.begin(),prip.end(),321)!=prip.end() && std::find(prip.begin(),prip.end(),-321)==prip.end() && std::find(prip.begin(),prip.end(),13)!=prip.end()){
    if(std::find(prip.begin(),prip.end(),321)!=prip.end() && std::find(prip.begin(),prip.end(),13)!=prip.end()){
      IsKaon = true;
      cout << "IS K+" << endl;

      if( std::find(prip.begin(),prip.end(),3212)!=prip.end() || std::find(prip.begin(),prip.end(),3222)!=prip.end() || std::find(prip.begin(),prip.end(),3112)!=prip.end() || std::find(prip.begin(),prip.end(),3122)!=prip.end()){
	IsAssociatedKaon = true;
	
	if(true_nu_pdg==14 && true_nu_ccnc==0 && true_nu_vtx_inCCInclusiveTPC==1){
	  
	  cout << "IS ASSOCIATED KAON" << endl;
	  cout << "Run " << run;
	  cout << " Subrun " << subrun;
	  cout << " Event " << event << endl;
	  cout << "Primary K+ ID: " << prim_k_id << endl;
	  
	}
      }
      else{
	IsSingleKaon = true;
	if(true_nu_pdg==14 && true_nu_ccnc==0 && true_nu_vtx_inCCInclusiveTPC==1){
	//if(true_nu_pdg==14 && true_nu_ccnc==0 && true_kaon_ke>=0 && true_nu_vtx_inCCInclusiveTPC==1){
	  cout << "IS SINGLE KAON!!!!" << endl;
	  cout << "Run " << run;
	  cout << " Subrun " << subrun;
	  cout << " Event " << event << endl;
	  cout << "Primary K+ ID: " << prim_k_id << endl;
	  
	  }
      }

      if( std::find(prip_k_dau.begin(),prip_k_dau.end(),-13)!=prip_k_dau.end()){
	IsMuBR = true;
	cout << "IS KAON -> MU" << endl;
      }
      else if( std::find(prip_k_dau.begin(),prip_k_dau.end(),211)!=prip_k_dau.end() ){
	IsPiBR = true;
	cout << "IS KAON -> PI" << endl;
      }

    }

    }
*/

  //-----------------------------


  if (isMC) {

    // get MC event weights
    art::InputTag eventweight_tag("eventweight");
    art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
    
    if (evt.getByLabel(eventweight_tag, eventweights_handle)) {
      std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
      art::fill_ptr_vector(eventweights, eventweights_handle);
      std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;
      //std::cout << "Weight: " <<std::endl;
      for (auto const& x : evtwgt_map) {
        //std::cout << "Weight - Size fun:\t" << x.first << " " << x.second.size() << std::endl;
        std::vector<float> wg;
      //}
      const std::vector<double> &weights = evtwgt_map.at(x.first);
      //event_weight = weights.front();
      //std::cout << "Event Weight:  " << event_weight << std::endl;
      for (unsigned int i=0;i<weights.size();i++) {
        //cout << i << "\t" << weights[i] << endl;
        wg.push_back(weights[i]);
      }

      event_weight[x.first] = wg;
      evtwgt_funcname.push_back(x.first);
      evtwgt_weight.push_back(wg);
      evtwgt_nweight.push_back(x.second.size());

      }
    }
    else {
      //std::cout << "***** Failed obtaining eventweight ******" << std::endl;
    }

    // get MCTruth
    evt.getByLabel("generator", mctruths);
    if (mctruths->size()!=1) {
      //std::cout << "Number of MCTruths objects in event " << mctruths->size() << std::endl;
      return;
      }
      simb::MCTruth mctruth = mctruths->at(0);

    // get MCParticles
    /*
    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
    if (evt.getByLabel(fLArG4ModuleLabel, mcParticleHandle)){
      art::fill_ptr_vector(ptList, mcParticleHandle); 
    }
    */

    // true neutrino information
    true_nu_energy = mctruth.GetNeutrino().Nu().E();
    true_nu_pdg = mctruth.GetNeutrino().Nu().PdgCode();
    true_nu_ccnc = mctruth.GetNeutrino().CCNC();//0=CC,1=NC
    true_nu_mode = mctruth.GetNeutrino().Mode();//0=QE,1=RES,2=DIS
    true_nu_vtx_x = mctruth.GetNeutrino().Nu().Vx();
    true_nu_vtx_y = mctruth.GetNeutrino().Nu().Vy();
    true_nu_vtx_z = mctruth.GetNeutrino().Nu().Vz();
    true_nu_vtx_inTPC = isInsideVolume("TPC",mctruth.GetNeutrino().Nu().Position().Vect());
    true_nu_vtx_in5cmTPC = isInsideVolume("5cmTPC",mctruth.GetNeutrino().Nu().Position().Vect());
    true_nu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",mctruth.GetNeutrino().Nu().Position().Vect());

    true_lepton_pdg = mctruth.GetNeutrino().Lepton().PdgCode();
    true_lepton_p = mctruth.GetNeutrino().Lepton().P();
    true_lepton_ke = mctruth.GetNeutrino().Lepton().E()-mctruth.GetNeutrino().Lepton().Mass();
    true_lepton_theta = mctruth.GetNeutrino().Lepton().Momentum().Theta();
    true_lepton_costheta = mctruth.GetNeutrino().Lepton().Momentum().CosTheta();
    true_lepton_phi = mctruth.GetNeutrino().Lepton().Momentum().Phi();
    TVector3 true_lepton_pvector = mctruth.GetNeutrino().Lepton().Momentum().Vect();

    // print true information
    /*
    cout << mctruths->at(0) << endl;
    for (auto const& pPart : ptList) {
      if (pPart->Process()!="Primary" && pPart->E()-pPart->Mass()<0.001) continue;
      cout << "trackId " << pPart->TrackId();
      cout << ", mother " << pPart->Mother();
      cout << ", pdg " << pPart->PdgCode();
      cout << ", process " << pPart->Process() << " - " << pPart->EndProcess();
      cout << ", ke " << pPart->E()-pPart->Mass() << " - " << pPart->EndE()-pPart->Mass();
      cout << ", X " << pPart->Vx() << " - " << pPart->EndX();
      cout << ", Z " << pPart->Vz() << " - " << pPart->EndZ();
      cout << ", delta Z " << pPart->EndZ()-pPart->Vz();
      cout << ", ndoughters " << pPart->NumberDaughters();
      cout << endl;
    }
    */

    // find the the highest momentum k+
    double true_kaon_maxp = -1;
    int true_kaon_maxp_id = -1;
    true_nkaons = 0;
    true_nprimary = 0;
    true_nhyperons = 0;
    //double true_hyperon_maxp = -1;
    //int true_hyperon_maxp_id = -1;

    int nparts = 0; 
    for (auto const& pPart : ptList) {
      if(pPart->Process()=="primary" && int(mctruths->at(0).Origin())==1){
	TLorentzVector mcstart, mcend;
	unsigned int pstarti, pendi;
	true_length[nparts] = length(*pPart, mcstart, mcend, pstarti, pendi);
	true_p[nparts]  = pPart->P();
	true_ke[nparts] = pPart->E()-pPart->Mass();
	true_theta[nparts] = pPart->Momentum().Theta();
	true_costheta[nparts] = pPart->Momentum().CosTheta();
	true_phi[nparts] = pPart->Momentum().Phi();
	true_ccmuon_angle[nparts] = pPart->Momentum().Angle(true_lepton_pvector);
	true_ccmuon_cosangle[nparts] = TMath::Cos(pPart->Momentum().Angle(true_lepton_pvector));
	//true_pvector = pPart->Momentum().Vect();
	true_end_ke[nparts] = pPart->EndE()-pPart->Mass();
	true_end_x[nparts] = pPart->EndX();
	true_end_y[nparts] = pPart->EndY();
	true_end_z[nparts] = pPart->EndZ();
	true_end_inTPC[nparts] = isInsideVolume("TPC",pPart->EndPosition().Vect());
	true_end_in5cmTPC[nparts] = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
	true_end_inCCInclusiveTPC[nparts] = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
      }

      if (pPart->Process()=="primary" && pPart->PdgCode()==321) {
	Process = pPart->Process();
	EndProcess = pPart->EndProcess();
        if(int(mctruths->at(0).Origin())==1) { // make sure origin is neutrino beam
          if (true_kaon_maxp < pPart->P()) {
            true_kaon_maxp = pPart->P();
            true_kaon_maxp_id = pPart->TrackId();
            true_nkaons++;
          }
        }  
      }

      if (pPart->Process()=="primary") {
	//Process = pPart->Process();
	//EndProcess = pPart->EndProcess();
        if(int(mctruths->at(0).Origin())==1) { // make sure origin is neutrino beam
            true_nprimary++;
        }  
      }

      if (pPart->Process()=="primary" && (pPart->PdgCode()==3212 || pPart->PdgCode()==3122 || pPart->PdgCode()==3222 || pPart->PdgCode()==3112)) {
	//Process = pPart->Process();
	//EndProcess = pPart->EndProcess();
        if(int(mctruths->at(0).Origin())==1) { // make sure origin is neutrino beam
          //if (true_hyperon_maxp < pPart->P()) {
            //true_hyperon_maxp = pPart->P();
            //true_hyperon_maxp_id = pPart->TrackId();
            true_nhyperons++;
	    //}
        }  
      }
      nparts++;
    }

    // true signal interaction: nu CC with at least one k+ // selection with true info
    if (true_nu_pdg==14 && true_nu_ccnc==0 && true_nkaons>0) {

      vector<int> decay_daughters;
      vector<int> decay_muplus_id;
      vector<int> decay_piplus_id;
      vector<int> inelastic_daughters;
      vector<int> inelastic_kaplus_id;
      TVector3 true_kaon_pvector;
      TVector3 true_dau_muon_pvector;
      TVector3 true_dau_pip_pvector;
      TVector3 true_dau_pin_pvector;

      cout << "true_kaon_maxp_id: " << true_kaon_maxp_id << endl;

      for (auto const& pPart : ptList) {

        // primary kaon+ (w/ highest momentum)
        if (pPart->TrackId()==true_kaon_maxp_id) {
          TLorentzVector mcstart, mcend;
          unsigned int pstarti, pendi;
          true_kaon_length = length(*pPart, mcstart, mcend, pstarti, pendi);
          true_kaon_p  = pPart->P();
          true_kaon_ke = pPart->E()-pPart->Mass();
          true_kaon_theta = pPart->Momentum().Theta();
          true_kaon_costheta = pPart->Momentum().CosTheta();
          true_kaon_phi = pPart->Momentum().Phi();
          true_kaon_ccmuon_angle = pPart->Momentum().Angle(true_lepton_pvector);
          true_kaon_ccmuon_cosangle = TMath::Cos(pPart->Momentum().Angle(true_lepton_pvector));
          true_kaon_pvector = pPart->Momentum().Vect();
          true_kaon_end_ke = pPart->EndE()-pPart->Mass();

	  true_kaon_start_x = pPart->Vx();
	  true_kaon_start_y = pPart->Vy();
	  true_kaon_start_z = pPart->Vz();
          true_kaon_end_x = pPart->EndX();
          true_kaon_end_y = pPart->EndY();
          true_kaon_end_z = pPart->EndZ();
          true_kaon_end_inTPC = isInsideVolume("TPC",pPart->EndPosition().Vect());
          true_kaon_end_in5cmTPC = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
          true_kaon_end_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
        }

        // find kaon+ daughters
        else if (pPart->Mother()==true_kaon_maxp_id) {

	  //cout << "pPart->TrackId(): " << pPart->TrackId() << ", it's PDG: " << pPart->PdgCode() << endl;

          // decay daughters
          if (pPart->Process()=="Decay") {
            decay_daughters.push_back(pPart->PdgCode());
	    //cout << "decay_daughters: " << pPart->PdgCode() << endl;
            if (pPart->PdgCode()==-13) {
              decay_muplus_id.push_back(pPart->TrackId());

	      TLorentzVector mcstart, mcend;
	      unsigned int pstarti, pendi;
	      true_dau_muon_length = length(*pPart, mcstart, mcend, pstarti, pendi);
	      true_dau_muon_p  = pPart->P();
	      true_dau_muon_ke = pPart->E()-pPart->Mass();
	      true_dau_muon_theta = pPart->Momentum().Theta();
	      true_dau_muon_costheta = pPart->Momentum().CosTheta();
	      true_dau_muon_phi = pPart->Momentum().Phi();
	      true_dau_muon_pvector = pPart->Momentum().Vect();
	      true_dau_muon_end_ke = pPart->EndE()-pPart->Mass();
	      
	      true_dau_muon_start_x = pPart->Vx();
	      true_dau_muon_start_y = pPart->Vy();
	      true_dau_muon_start_z = pPart->Vz();
	      true_dau_muon_end_x = pPart->EndX();
	      true_dau_muon_end_y = pPart->EndY();
	      true_dau_muon_end_z = pPart->EndZ();

            }
            else if (pPart->PdgCode()==211) {
              decay_piplus_id.push_back(pPart->TrackId());

	      TLorentzVector mcstart, mcend;
	      unsigned int pstarti, pendi;
	      true_dau_pip_length = length(*pPart, mcstart, mcend, pstarti, pendi);
	      true_dau_pip_p  = pPart->P();
	      true_dau_pip_ke = pPart->E()-pPart->Mass();
	      true_dau_pip_theta = pPart->Momentum().Theta();
	      true_dau_pip_costheta = pPart->Momentum().CosTheta();
	      true_dau_pip_phi = pPart->Momentum().Phi();
	      true_dau_pip_pvector = pPart->Momentum().Vect();
	      true_dau_pip_end_ke = pPart->EndE()-pPart->Mass();
	      
	      true_dau_pip_start_x = pPart->Vx();
	      true_dau_pip_start_y = pPart->Vy();
	      true_dau_pip_start_z = pPart->Vz();
	      true_dau_pip_end_x = pPart->EndX();
	      true_dau_pip_end_y = pPart->EndY();
	      true_dau_pip_end_z = pPart->EndZ();

            }
	    else if (pPart->PdgCode()==111) {

	      TLorentzVector mcstart, mcend;
	      unsigned int pstarti, pendi;
	      true_dau_pin_length = length(*pPart, mcstart, mcend, pstarti, pendi);
	      true_dau_pin_p  = pPart->P();
	      true_dau_pin_ke = pPart->E()-pPart->Mass();
	      true_dau_pin_theta = pPart->Momentum().Theta();
	      true_dau_pin_costheta = pPart->Momentum().CosTheta();
	      true_dau_pin_phi = pPart->Momentum().Phi();
	      true_dau_pin_pvector = pPart->Momentum().Vect();
	      true_dau_pin_end_ke = pPart->EndE()-pPart->Mass();
	      
	      true_dau_pin_start_x = pPart->Vx();
	      true_dau_pin_start_y = pPart->Vy();
	      true_dau_pin_start_z = pPart->Vz();
	      true_dau_pin_end_x = pPart->EndX();
	      true_dau_pin_end_y = pPart->EndY();
	      true_dau_pin_end_z = pPart->EndZ();

    cout << "true_track_length: " << true_dau_pin_length << endl;

	    }
            //cout << "Decay " << pPart->PdgCode() << endl;
          }

          // inelastic interaction daughters
          if (pPart->Process()=="kaon+Inelastic") {
            inelastic_daughters.push_back(pPart->PdgCode());
            if (pPart->PdgCode()==321) {
              inelastic_kaplus_id.push_back(pPart->TrackId());
            }
            //cout << "Inelastic " << pPart->PdgCode() << endl;
          }

        }

      }//MC particles loop

      true_kaon_ndaughters = decay_daughters.size() + inelastic_daughters.size();
      true_kaon_ndaughters_decay = decay_daughters.size();
      true_kaon_ndaughters_inelastic = inelastic_daughters.size();

      int n_decay_muplus     = decay_muplus_id.size();
      int n_decay_piplus     = decay_piplus_id.size();
      int n_decay_numu       = count(decay_daughters.begin(),decay_daughters.end(),14);
      int n_decay_pi0        = count(decay_daughters.begin(),decay_daughters.end(),111);
      int n_inelastic_kaplus = inelastic_kaplus_id.size();
      int n_inelastic_piplus = count(inelastic_daughters.begin(),inelastic_daughters.end(),211);
      int n_inelastic_proton = count(inelastic_daughters.begin(),inelastic_daughters.end(),2212);

      true_kaon_ndecmup = n_decay_muplus;
      true_kaon_ndecpip = n_decay_piplus;
      true_kaon_ninekap = n_inelastic_kaplus;
      true_kaon_ninepip = n_inelastic_piplus;
      true_kaon_ninepro = n_inelastic_proton;
    
      // define kaon end process
      int daughter_id = -1;
      if (decay_daughters.size()) {
        //kaon+ -> muon+ numu
        if (decay_daughters.size()==2 && n_decay_muplus==1 && n_decay_numu==1) {
          true_kaon_end_process = 0;
          daughter_id = decay_muplus_id.front();
        }
        //kaon+ -> pion+ + pi0
        else if (decay_daughters.size()==2 && n_decay_piplus==1 && n_decay_pi0==1) {
          true_kaon_end_process = 1;
          daughter_id = decay_piplus_id.front();
        }
        //kaon+ -> muon+ nopion+
        else if (n_decay_muplus==1 && n_decay_piplus==0) {
          true_kaon_end_process = 2;
          daughter_id = decay_muplus_id.front();
        }
        //kaon+ -> pion+ and nomuon+
        else if (n_decay_muplus==0 && n_decay_piplus==1) {
          true_kaon_end_process = 3;
          daughter_id = decay_piplus_id.front();
        }
        //other decays
        else {
          true_kaon_end_process = 4;
        }
      }

      if (inelastic_daughters.size()) {
        //inelasitic with kaon+
        if (n_inelastic_kaplus==1) {
          true_kaon_end_process = 5;
          daughter_id = inelastic_kaplus_id.front();
        }
        //other inelastic
        else {
          true_kaon_end_process = 6;
        }
      }

      //this shouldn't happen but it does! rarely
      if (decay_daughters.size() && inelastic_daughters.size()) {
        // -- cout << "Primary kaon+ decay and inelastic at the same time!!!" << endl;
        true_kaon_end_process = 7;
      }

      // -- cout << "True kaon end process " << true_kaon_end_process << endl;

      //check kaon daughter if it exists
      //only if there is one pion+, one muon+ or one kaon+
      if (daughter_id!=-1) {
        for (auto const& pPart : ptList) {
      
          // kaon+ daughter
          if (pPart->TrackId()==daughter_id) {
            TLorentzVector mcstart, mcend;
            unsigned int pstarti, pendi;
            true_kaon_daughter_length = length(*pPart, mcstart, mcend, pstarti, pendi);
            true_kaon_daughter_p = pPart->P();
            true_kaon_daughter_ke = pPart->E()-pPart->Mass();
            true_kaon_daughter_theta = pPart->Momentum().Theta();
            true_kaon_daughter_costheta = pPart->Momentum().CosTheta();
            true_kaon_daughter_angle = pPart->Momentum().Angle(true_kaon_pvector);
            true_kaon_daughter_cosangle = TMath::Cos(pPart->Momentum().Angle(true_kaon_pvector));
            true_kaon_daughter_pdg = pPart->PdgCode();

	    true_kaon_daughter_start_x = pPart->Vx();
	    true_kaon_daughter_start_y = pPart->Vy();
	    true_kaon_daughter_start_z = pPart->Vz();
            true_kaon_daughter_end_x = pPart->EndX();
            true_kaon_daughter_end_y = pPart->EndY();
            true_kaon_daughter_end_z = pPart->EndZ();
            true_kaon_daughter_end_inTPC = isInsideVolume("TPC",pPart->EndPosition().Vect());
            true_kaon_daughter_end_in5cmTPC = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
            true_kaon_daughter_end_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
            break;
          }
      
        }//MC particles loop
      }//kaon daughter exists

    }//is true signal

  
    // find the the highest momentum k-
    double true_kaon_maxp_anti = -1;
    int true_kaon_maxp_id_anti = -1;
    true_nkaons_anti = 0;

    for (auto const& pPart : ptList) {
      if (pPart->Process()=="primary" && pPart->PdgCode()==-321) {
        if(int(mctruths->at(0).Origin())==1) { // make sure origin is neutrino beam
          if (true_kaon_maxp_anti < pPart->P()) {
            true_kaon_maxp_anti = pPart->P();
            true_kaon_maxp_id_anti = pPart->TrackId();
            true_nkaons_anti++;
          }
        }  
      }
    }
  
    // true signal interaction: nu CC with at least one k- // selection with true info
    if (true_nu_pdg==14 && true_nu_ccnc==0 && true_nkaons_anti>0) {

      vector<int> decay_daughters_anti;
      vector<int> decay_muplus_id_anti;
      vector<int> decay_piplus_id_anti;
      vector<int> inelastic_daughters_anti;
      vector<int> inelastic_kaplus_id_anti;
      TVector3 true_kaon_pvector_anti;

      for (auto const& pPart : ptList) {

        // primary kaon- (w/ highest momentum)
        if (pPart->TrackId()==true_kaon_maxp_id_anti) {
          TLorentzVector mcstart_anti, mcend_anti;
          unsigned int pstarti_anti, pendi_anti;
          true_kaon_length_anti = length(*pPart, mcstart_anti, mcend_anti, pstarti_anti, pendi_anti);
          true_kaon_p_anti  = pPart->P();
          true_kaon_ke_anti = pPart->E()-pPart->Mass();
          true_kaon_theta_anti = pPart->Momentum().Theta();
          true_kaon_costheta_anti = pPart->Momentum().CosTheta();
          true_kaon_phi_anti = pPart->Momentum().Phi();
          true_kaon_ccmuon_angle_anti = pPart->Momentum().Angle(true_lepton_pvector);
          true_kaon_ccmuon_cosangle_anti = TMath::Cos(pPart->Momentum().Angle(true_lepton_pvector));
          true_kaon_pvector_anti = pPart->Momentum().Vect();
          true_kaon_end_ke_anti = pPart->EndE()-pPart->Mass();
          true_kaon_end_x_anti = pPart->EndX();
          true_kaon_end_y_anti = pPart->EndY();
          true_kaon_end_z_anti = pPart->EndZ();
          true_kaon_end_inTPC_anti = isInsideVolume("TPC",pPart->EndPosition().Vect());
          true_kaon_end_in5cmTPC_anti = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
          true_kaon_end_inCCInclusiveTPC_anti = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
        }

        // find kaon- daughters
        else if (pPart->Mother()==true_kaon_maxp_id_anti) {

          // decay daughters
          if (pPart->Process()=="Decay") {
            decay_daughters_anti.push_back(pPart->PdgCode());
            if (pPart->PdgCode()==13) {
              decay_muplus_id_anti.push_back(pPart->TrackId());
            }
            else if (pPart->PdgCode()==-211) {
              decay_piplus_id_anti.push_back(pPart->TrackId());
            }
            //cout << "Decay " << pPart->PdgCode() << endl;
          }

          // inelastic interaction daughters
          if (pPart->Process()=="kaon+Inelastic") {
            inelastic_daughters_anti.push_back(pPart->PdgCode());
            if (pPart->PdgCode()==-321) {
              inelastic_kaplus_id_anti.push_back(pPart->TrackId());
            }
            //cout << "Inelastic " << pPart->PdgCode() << endl;
          }

        }

      }//MC particles loop
    
      true_kaon_ndaughters_anti = decay_daughters_anti.size() + inelastic_daughters_anti.size();
      true_kaon_ndaughters_decay_anti = decay_daughters_anti.size();
      true_kaon_ndaughters_inelastic_anti = inelastic_daughters_anti.size();

      int n_decay_muplus_anti     = decay_muplus_id_anti.size();
      int n_decay_piplus_anti     = decay_piplus_id_anti.size();
      int n_decay_numu_anti       = count(decay_daughters_anti.begin(),decay_daughters_anti.end(),-14);
      int n_decay_pi0_anti        = count(decay_daughters_anti.begin(),decay_daughters_anti.end(),111);
      int n_inelastic_kaplus_anti = inelastic_kaplus_id_anti.size();
      int n_inelastic_piplus_anti = count(inelastic_daughters_anti.begin(),inelastic_daughters_anti.end(),-211);
      int n_inelastic_proton_anti = count(inelastic_daughters_anti.begin(),inelastic_daughters_anti.end(),-2212);

      true_kaon_ndecmup_anti = n_decay_muplus_anti;
      true_kaon_ndecpip_anti = n_decay_piplus_anti;
      true_kaon_ninekap_anti = n_inelastic_kaplus_anti;
      true_kaon_ninepip_anti = n_inelastic_piplus_anti;
      true_kaon_ninepro_anti = n_inelastic_proton_anti;

      // define kaon end process
      int daughter_id_anti = -1;
      if (decay_daughters_anti.size()) {
        //kaon+ -> muon+ numu
        if (decay_daughters_anti.size()==2 && n_decay_muplus_anti==1 && n_decay_numu_anti==1) {
          true_kaon_end_process_anti = 0;
          daughter_id_anti = decay_muplus_id_anti.front();
        }
        //kaon+ -> pion+ + pi0
        else if (decay_daughters_anti.size()==2 && n_decay_piplus_anti==1 && n_decay_pi0_anti==1) {
          true_kaon_end_process_anti = 1;
          daughter_id_anti = decay_piplus_id_anti.front();
        }
        //kaon+ -> muon+ nopion+
        else if (n_decay_muplus_anti==1 && n_decay_piplus_anti==0) {
          true_kaon_end_process_anti = 2;
          daughter_id_anti = decay_muplus_id_anti.front();
        }
        //kaon+ -> pion+ and nomuon+
        else if (n_decay_muplus_anti==0 && n_decay_piplus_anti==1) {
          true_kaon_end_process_anti = 3;
          daughter_id_anti = decay_piplus_id_anti.front();
        }
        //other decays
        else {
          true_kaon_end_process_anti = 4;
        }
      }

      if (inelastic_daughters_anti.size()) {
        //inelasitic with kaon-
        if (n_inelastic_kaplus_anti==1) {
          true_kaon_end_process_anti = 5;
          daughter_id_anti = inelastic_kaplus_id_anti.front();
        }
        //other inelastic
        else {
          true_kaon_end_process_anti = 6;
        }
      }

      //this shouldn't happen but it does! rarely
      if (decay_daughters_anti.size() && inelastic_daughters_anti.size()) {
        // -- cout << "Primary kaon+ decay and inelastic at the same time!!!" << endl;
        true_kaon_end_process_anti = 7;
      }

      // -- cout << "True kaon end process " << true_kaon_end_process << endl;

      //check kaon daughter if it exists
      //only if there is one pion+, one muon+ or one kaon+
      if (daughter_id_anti!=-1) {
        for (auto const& pPart : ptList) {
      
          // kaon- daughter
          if (pPart->TrackId()==daughter_id_anti) {
            TLorentzVector mcstart, mcend;
            unsigned int pstarti, pendi;
            true_kaon_daughter_length_anti = length(*pPart, mcstart, mcend, pstarti, pendi);
            true_kaon_daughter_p_anti = pPart->P();
            true_kaon_daughter_ke_anti = pPart->E()-pPart->Mass();
            true_kaon_daughter_theta_anti = pPart->Momentum().Theta();
            true_kaon_daughter_costheta_anti = pPart->Momentum().CosTheta();
            true_kaon_daughter_angle_anti = pPart->Momentum().Angle(true_kaon_pvector_anti);
            true_kaon_daughter_cosangle_anti = TMath::Cos(pPart->Momentum().Angle(true_kaon_pvector_anti));
            true_kaon_daughter_pdg_anti = pPart->PdgCode();
            true_kaon_daughter_end_x_anti = pPart->EndX();
            true_kaon_daughter_end_y_anti = pPart->EndY();
            true_kaon_daughter_end_z_anti = pPart->EndZ();
            true_kaon_daughter_end_inTPC_anti = isInsideVolume("TPC",pPart->EndPosition().Vect());
            true_kaon_daughter_end_in5cmTPC_anti = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
            true_kaon_daughter_end_inCCInclusiveTPC_anti = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
            break;
          }
      
        }//MC particles loop
      }//kaon daughter exists

    }//is true signal

    this->filter(evt);
  
    # if 0
    //Collect the PFParticles from the event. This is the core!                                                                                                                                                                                                                                                                             
    art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);
    std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
    art::fill_ptr_vector(pfParticleVector,pfParticleHandle);
    //So a cross check
    /*                                                                                                                                                                                                                                                                                                                      
    if (!pfParticleHandle.isValid())
      {
	mf::LogDebug("SinglePhoton") << "  Failed to find the PFParticles.\n";
	if(m_run_pi0_filter)
	  return false;
	else
	  return true;
      }
    */

    //This is another pandora helper. I don't like PFParticle ID lookups but I guess lets keep for now;                                                                                                                                                                                                                                     
    // Produce a map of the PFParticle IDs for fast navigation through the hierarchy                                                                                                                                                                                                                                                        
    //    typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;
    //    CCKaonAnalyzerRebuild::
  
    
    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);
  
    //Slices                                                                                                                                                                                                                                                                                                                                
    art::ValidHandle<std::vector<recob::Slice>> const & sliceHandle  = evt.getValidHandle<std::vector<recob::Slice>>(m_pandoraLabel);
    std::vector<art::Ptr<recob::Slice>> sliceVector;
    art::fill_ptr_vector(sliceVector,sliceHandle);
    //And some associations                                                                                                                                                                                                                                                                                                                 
    art::FindManyP<recob::PFParticle> pfparticles_per_slice(sliceHandle, evt, m_pandoraLabel);
    art::FindManyP<recob::Hit> hits_per_slice(sliceHandle, evt, m_pandoraLabel);

    std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::PFParticle>> > sliceToPFParticlesMap;
    std::map<int, std::vector<art::Ptr<recob::PFParticle>> > sliceIDToPFParticlesMap;
    for(size_t i=0; i< sliceVector.size(); ++i){
      auto slice = sliceVector[i];
      sliceToPFParticlesMap[slice] =pfparticles_per_slice.at(slice.key());
      sliceIDToPFParticlesMap[slice->ID()] = pfparticles_per_slice.at(slice.key());
    }

    std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::Hit>> > sliceToHitsMap;
    std::map<int, std::vector<art::Ptr<recob::Hit>> > sliceIDToHitsMap;
    for(size_t i=0; i< sliceVector.size(); ++i){
      auto slice = sliceVector[i];
      sliceToHitsMap[slice] =hits_per_slice.at(slice.key());
      sliceIDToHitsMap[slice->ID()] = hits_per_slice.at(slice.key());
    }

    if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get PandoraMetadata"<<std::endl;
    //add the associaton between PFP and metadata, this is important to look at the slices and scores                                                                                                                                                                                                                                       
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt, m_pandoraLabel);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > pfParticleToMetadataMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
      const art::Ptr<recob::PFParticle> pfp = pfParticleVector[i];
      pfParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
    }
    
    //these are all filled in analyze slice                                                                                                                                                                                                                                                                                             
    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind                                                                                                                                                                                               
    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> primaryPFPSliceIdVec; //stores a pair of only the primary PFP's in the event and the slice ind                                                                                                                                                                              
    std::map<int, double> sliceIdToNuScoreMap; //map between a slice Id and neutrino score                                                                                                                                                                                                                                              
    std::map<art::Ptr<recob::PFParticle>, bool> PFPToClearCosmicMap; //returns true for clear cosmic, false otherwise                                                                                                                                                                                                                   
    std::map<art::Ptr<recob::PFParticle>, int> PFPToSliceIdMap; //returns the slice id for all PFP's                                                                                                                                                                                                                                    
    std::map<art::Ptr<recob::PFParticle>,bool> PFPToNuSliceMap;
    std::map<art::Ptr<recob::PFParticle>,double> PFPToTrackScoreMap;
    std::map<int, int> sliceIdToNumPFPsMap;
    //std::cout<<"SinglePhoton::analyze::AnalyzeSlice()\t||\t Starting"<<std::endl;

    
    this->AnalyzeSlices(pfParticleToMetadataMap, pfParticleMap,  primaryPFPSliceIdVec, sliceIdToNuScoreMap, PFPToClearCosmicMap, PFPToSliceIdMap, PFPToNuSliceMap, PFPToTrackScoreMap);
    //std::cout<<"There are "<< allPFPSliceIdVec.size()<<" pfp-slice id matches stored in the vector"<<std::endl;                                                                                                                                                                                                                       
    //std::cout<<"SinglePhoton::analyze\t||\tthe number of PPF's with stored clear cosmic info is "<<PFPToClearCosmicMap.size()<<std::endl;
    //std::cout<<"SinglePhoton::analyze\t||\tthe number of PFP's stored in the PFPToSliceIdMap is "<<PFPToSliceIdMap.size()<<std::endl;
    if (PFPToSliceIdMap.size() < 1){
      std::cout<<"ERROR, not storing PFP's in PFPToSliceIdMap"<<std::endl;
    }

    for (auto pair:PFPToNuSliceMap){
      auto pfp = pair.first;
      auto is_nuslice = pair.second;
      if (is_nuslice){
	std::cout<<"pfp in nuslice "<<pfp->Self()<<std::endl;
      }

    }

    for (auto pair:sliceIDToPFParticlesMap){
      std::vector<art::Ptr<recob::PFParticle>> pfp_vec = pair.second;
      int slice_id = pair.first;
      //if (slice_vec[0]->Slice() != PFPToSliceIdMap[pfp] )                                                                                                                                                                                                                                                                           
      for(auto pfp: pfp_vec){
	if (slice_id != PFPToSliceIdMap[pfp] && PFPToSliceIdMap[pfp]>=0){
	  //std::cout<<"sliceIDToPFParticlesMap[slice->ID()] for pfp "<<pfp->Self()<<" is slice "<< slice_id<< "but PFPToSliceIdMap[pfp] = "<<PFPToSliceIdMap[pfp]<<std::endl;
	}
      }

    }    
    #endif

  }//isMC

  // Check if event passed the NuCC inclusive filter
  reco_nu_cc_filter = false;
  string process(isMC ? "OverlayFiltersPostStage2" : "DataFiltersPostStage2");
  art::InputTag trigResInputTag("TriggerResults","",process.data()); // the last is the name of process where the filters were run
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
  /*
  art::Handle<std::vector< recob::PFParticle > > pfphandle;
  vector< art::Ptr< recob::PFParticle > > pfpvector;
  vector< art::Ptr< recob::PFParticle > > pfpvectorPandora;
  if(evt.getByLabel(fPFParticleLabel,pfphandle)) art::fill_ptr_vector(pfpvector,pfphandle);
  if(!pfpvector.size()) return;

  art::FindManyP< larpandoraobj::PFParticleMetadata > pfpmeta(pfphandle,evt,"pandora");
  art::FindManyP<recob::Slice> nuslice(pfpvector, evt, "pandora");

  size_t pfpkey = 999999; 

  for(const art::Ptr< recob::PFParticle > &pfp : pfpvector){

    const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfpmeta.at(pfp.key()));

    if (pfParticleMetadataList.empty()) continue;

    for (const art::Ptr<larpandoraobj::PFParticleMetadata> &pfpmd : pfParticleMetadataList ) {
      auto pfParticlePropertiesMap = pfpmd->GetPropertiesMap();
      if (pfParticlePropertiesMap.empty()) continue;
      for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	if(!(it->first == "TrackScore")) continue; 
	cout << "PFParticle properties: " << it->first << " Key: " << pfpmd.key() << " Value: " << it->second << endl;
	pfpkey = pfp.key();
	pfpvectorPandora.push_back(pfp);
      }
    }
    auto slices = nuslice.at(pfp.key());
    cout << "Slices: " << slices.size() << endl;
  }
  */

  art::Handle< vector< recob::PFParticle > > pfphandle;
  vector< art::Ptr< recob::PFParticle > > pfpvector;
  vector< art::Ptr< recob::PFParticle > > pfpvectorPandora;

  if(evt.getByLabel(fPFParticleLabel,pfphandle)) art::fill_ptr_vector(pfpvector,pfphandle);
  if(!pfpvector.size()) return;

  art::FindManyP< larpandoraobj::PFParticleMetadata > pfpmeta(pfphandle,evt,"pandora");
  art::FindManyP<recob::Slice> nuslice(pfpvector, evt, "pandora");

  size_t pfpkey = 999999;

  for(const art::Ptr< recob::PFParticle > &pfp : pfpvector){

    const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfpmeta.at(pfp.key()));

    if (pfParticleMetadataList.empty()) continue;


    for (const art::Ptr<larpandoraobj::PFParticleMetadata> &pfpmd : pfParticleMetadataList ) {
      auto pfParticlePropertiesMap = pfpmd->GetPropertiesMap();
      //for (const larpandoraobj::PFParticleMetadata pfpmd : pfParticleMetadataList ) {
      //auto pfParticlePropertiesMap = pfpmd.GetPropertiesMap();
      if (pfParticlePropertiesMap.empty()) continue; //cout << pfParticlePropertiesMap.size() << endl;
      for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	if(!(it->first == "TrackScore")) continue;
	//cout << "PFParticle properties: " << it->first << " Key: " << pfpmd.key() << " Value: " << it->second << endl;
	pfpkey = pfp.key();
	pfpvectorPandora.push_back(pfp);
      }
    }// for all metadata items in the particle metadata
    if(pfp.key() == pfpkey){
      /*
      cout << "PFP ID: " <<  pfp->Self() <<endl;
      cout << "PFP IsPrimary: " << pfp->IsPrimary() <<endl;
      cout << "PFP parent: " << pfp->Parent() <<endl;
      cout << "PFP PDG: " << pfp->PdgCode() <<endl;
      cout << "PFP key: " << pfp.key() << "=====================\n" <<endl;
      */
    }

    auto slices = nuslice.at(pfp.key());
    //cout << "Slices: " << slices.size() << endl;
    //" ID: " << slice->ID() endl;

    //for(size_t i=0; i< slices.size(); ++i){
    //auto slice = slices[i];
    //}

  }// for entries in list



  if (!reco_nu_cc_filter) {
    // -- std::cout << "Event didn't pass NuCC filter" << std::endl;
    fEventTree->Fill();
    return;
  }

  // -- std::cout << "Event passed NuCC filter" << std::endl;
  cut_1=1;

  // Collect all recontructed particles
  art::Handle<std::vector<recob::PFParticle>> pfparticles;
  evt.getByLabel(m_pfp_producer, pfparticles);

  //cout << "pfpvector.size(): " << pfpvector.size() << endl;
  //cout << "pfparticles->size(): " << pfparticles->size() << endl;
  //cout << "pfphandle->size(): " << pfphandle->size() << endl;

  if (pfparticles->size()==0) {
    std::cout << "No PFParticles found" << std::endl;
    fEventTree->Fill();
    return;
  }
  // -- std::cout << "Number of PFParticles " << pfparticles->size() << std::endl;
  cut_2=1;	  

  // Get PFParticle associations
  art::FindManyP<anab::T0> pfp_muon_assn(pfparticles, evt, "NuCCproducer");
  if(!pfp_muon_assn.isValid()){
    // -- cout << "PFParticle-T0 associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  cut_3=1;	  

  art::FindManyP<recob::Track> pfparticleTrackAssn(pfparticles, evt, "pandora");
  if(!pfparticleTrackAssn.isValid()){
    // -- cout << "PFParticle-Track associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  cut_4=1;	  

  art::FindManyP<recob::Vertex> pfparticleVertexAssn(pfparticles, evt, "pandora");
  if(!pfparticleVertexAssn.isValid()){
    // -- cout << "PFParticle-Vertex associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  cut_5=1;	  

  // Find recontructed neutrino (there should be one)
  lar_pandora::PFParticleVector pfneutrinos(0);
  for (unsigned int i=0; i<pfparticles->size(); ++i) {

    art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);

    if (pfparticle->IsPrimary() && pfparticle->PdgCode()==14) {
      pfneutrinos.push_back(pfparticle);
    }

  }

  if (pfneutrinos.size() != 1) {
    // --cout << "Number of neutrinos is not one" << endl;
    cout << "Number of neutrinos " << pfneutrinos.size() << endl;
    fEventTree->Fill();
    return;
  }
  cut_6=1;	  

  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
  //--cout << "Found one neutrino";
  //--cout << " ID " << pfnu->Self();
  //--cout << " PDG " << pfnu->PdgCode();
  //--cout << " Number of daughters " << pfnu->Daughters().size() << endl;

  reco_nu_vtx_x = pfparticleVertexAssn.at(pfnu.key()).front()->position().X();
  reco_nu_vtx_y = pfparticleVertexAssn.at(pfnu.key()).front()->position().Y();
  reco_nu_vtx_z = pfparticleVertexAssn.at(pfnu.key()).front()->position().Z();
  //--cout << "Neutrino vertex (x,y,z) = " << reco_nu_vtx_x << ", " << reco_nu_vtx_y << ", " << reco_nu_vtx_z << endl;
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
      //--cout << "Neutrino daughter";
      //--cout << " ID: " << pfparticle->Self();
     //-- cout << " PDG: " << pfparticle->PdgCode();
     //-- cout << " Track key: " << track.key() << endl;
     //-- cout << " Z : " << track->Start().Z() << " - " << track->End().Z() << endl;
     //-- cout << " length : " << track->Length() << endl;

      reco_nu_daughters_id.push_back(track.key());

      // CC muon has a T0 associated
      if (pfp_muon_assn.at(i).size()==1) {
        pfmuons.push_back(pfparticle);
      }

    }

  }

  reco_nu_ndaughters = reco_nu_daughters_id.size();
  //--std::cout << "Number of neutrino daughters with one associated track " << reco_nu_ndaughters << std::endl;
  reco_nu_cc_nmue = pfmuons.size();

  if (pfmuons.size()!=1) {
    //--std::cout << "Number of CC inclusive muons is not 1: " << pfmuons.size() << std::endl;
    fEventTree->Fill();
    return;
  }
  cut_7=1;	  

  art::Ptr<recob::PFParticle> pfmuon = pfmuons.front();
  art::Ptr<recob::Track> trkmuon = pfparticleTrackAssn.at(pfmuon.key()).front();
  //--cout << "Found CC muon";
  //--cout << " ID " << pfmuon->Self();
  //--cout << " PDG " << pfmuon->PdgCode();
  //--cout << " Parent " << pfmuon->Parent();
  //--cout << " Track key " << trkmuon.key() << endl;

  reco_ccmu_vtx_x = pfparticleVertexAssn.at(pfmuon.key()).front()->position().X();
  reco_ccmu_vtx_y = pfparticleVertexAssn.at(pfmuon.key()).front()->position().Y();
  reco_ccmu_vtx_z = pfparticleVertexAssn.at(pfmuon.key()).front()->position().Z();
  //--cout << "CC muon start (x,y,z) = " << reco_ccmu_vtx_x << ", " << reco_ccmu_vtx_y << ", " << reco_ccmu_vtx_z << endl;
  reco_ccmu_vtx_inTPC = isInsideVolume("TPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);
  reco_ccmu_vtx_in5cmTPC = isInsideVolume("5cmTPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);
  reco_ccmu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);

  // get hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitlist, hitListHandle);
  }
  //--cout << "Found " << hitlist.size() << " hits" << endl;

  // get tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) {
    art::fill_ptr_vector(tracklist, trackListHandle);
  }
  //--cout << "Found " << tracklist.size() << " tracks" << endl;

  // get showers
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if(evt.getByLabel(fShowerModuleLabel,showerListHandle)) {
    art::fill_ptr_vector(showerlist, showerListHandle);
  }

  // get track associations
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  if(!fmcal.isValid()){
    //--cout << "Track-Calorimetry associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  cut_8=1;	  

  art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, evt, fHitTrackAssns);
  art::FindManyP<recob::Hit> hits_from_showers(showerListHandle, evt, fHitShowerAssns);
  if(!hits_from_tracks.isValid()){
    //--cout << "Track-Hit associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  cut_9=1;	  

  art::FindManyP<anab::ParticleID> trackPIDAssn(trackListHandle, evt, fPIDLabel);
  if(!trackPIDAssn.isValid()){
    //--cout << "Track PID associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  cut_10=1;


  // find track multiplicity around vertex
  int NTracks=tracklist.size();
  int ntracks = 0;

  int NShowers=showerlist.size();
  //if (showerListHandle.isValid()) NShowers=showerlist.size();
  //int nshowers = 0;

  //cout << "NTracks: " << NTracks << ", NShowers: " << NShowers << endl;

  //selection cuts
  std::vector<int> kaon_can_trkID;
  std::vector<int> muon_can_trkID;

  /*
  std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> trackToPFParticleMap;
  const std::map<art::Ptr<recob::PFParticle>, int> pfParticleToSliceIDMap;
  const std::map<int, std::vector<art::Ptr<recob::Hit>>> sliceIDToHitsMap;
  const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
      
  for(size_t t =0; t< tracklist.size(); t++){
    art::Ptr<recob::Track> track = trackslist[t];
    art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap[track];
    int sliceid = pfParticleToSliceIDMap.at(pfp);
    auto slicehits = sliceIDToHitsMap.at(sliceid);
    auto trackhits = pfParticleToHitsMap.at(pfp);

    std::cout<<"SinglePhoton::SSS\t||\ttrack "<<t<<" is in slice "<<sliceid<<" which has "<<slicehits.size()<<" hits. This track has  "<<trackhits.size()<<" of them. "<<std::endl;
    total_track_hits+=trackhits.size();
    if(nu_slice_id !=  sliceid && nu_slice_id != -999){
      std::cout<<"ERROR!! In Second Shower Search, the neutrino slice ID changed? this: "<<sliceid<<", last: "<<nu_slice_id<<std::endl;
      exit(EXIT_FAILURE);
    }
    nu_slice_id = sliceid;


    for(auto &h: trackhits){
      associated_hits.push_back(h);
    }

  }
  */

  // loop over tracks again to look for kaon track //->using reconinfo not true?
  //--cout << "Looking for kaon tracks from " << NTracks << " tracks available" << endl;
  //
std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>> assocMCPartt;
    if (isMC) {
      assocMCPartt = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> (new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, evt, fHitTruthAssns));
    }
  // loop over nu daughters
  trkf::TrackMomentumCalculator trkmom;
  //trkmom.SetMinLength(50);

  //initialiseCanvas();

  cout << "reco_nu_ndaughters: " << reco_nu_ndaughters << endl;
  for (int i=0; i<reco_nu_ndaughters; i++) {
  //for (int i=0; i<NTracks; i++) {

    std::map<int,int> mrgidpdg;
    std::map<int,int> mrgidpdg_3D;
    std::map<int,TVector3> mrgidmom;
    vector<TH1D*> h_angular_distribution_pfparticle_cheat;
    vector<TH2D*> h_angular_distribution_pfparticle_cheat_3D;
    vector<TH1D*> h_angular_distribution_pfparticle;
    vector<TH2D*> h_angular_distribution_pfparticle_3D;
    vector<TH2D*> h_angular_distribution_pfparticle_surface;
    vector<TH3D*> h_angular_distribution_pfparticle_sphere;
    vector<bool> v_trk_flg;
    vector<bool> v_trk_flg_3D;
    vector<bool> v_trk_flg_surface;
    vector<bool> v_trk_flg_peak;
    vector<int> v_pdg;
    vector<int> v_pdg_3D;
    vector<int> v_pdg_peak;
    //vector<TVector2> view_peak_vector_cheat;
    std::map<double, TVector2, std::greater<>> view_peak_map;
    std::map<int, std::map<double, TVector2, std::greater<>>> view_peak_map_cheat;
    //vector<TVector2> view_peak_vector;
    vector<TVector2> best_peak_bins;
    std::map<int, vector<TVector2>> best_peak_bins_cheat;

    std::map<int, std::map<int, double>> angular_distribution_mrgid_map;
    std::map<int, std::map<int, std::map<int, double>>> angular_distribution_mrgid_map_3D;
    std::map<int, double> angular_distribution_map_track;
    std::map<int, double> angular_distribution_map_shower;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_track;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_shower;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_pfparticle;

    std::vector<art::Ptr<recob::Hit>> unavailable_hit_list;
    std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list;
    std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list_cheat;
    //std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list_cheat_pip;
    //std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list_cheat_mup;

    std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_vector;
    std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_cheat_vector;
    std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_cheat_pip_vector;
    std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_cheat_mup_vector;
    std::vector<art::Ptr<recob::Hit>> true_pi0_hit_list;

    std::vector<art::Ptr<recob::Hit>> hits_from_reco;
    std::vector<art::Ptr<recob::Hit>> hits_from_true_pip;
    std::vector<art::Ptr<recob::Hit>> hits_from_true_mup;

    //art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);

    int currentMergedId = 1;
    int currentMergedId_3D = 1;


    art::Ptr<recob::Track> ptrack(trackListHandle,reco_nu_daughters_id[i]);
    //art::Ptr<recob::Track> ptrack(trackListHandle,i);
    const recob::Track& track = *ptrack;

    /*
    std::map<art::Ptr<recob::PFParticle>,bool> PFPToNuSliceMap;
    std::map<art::Ptr<recob::PFParticle>, int> PFPToSliceIdMap;
    std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> trackToNuPFParticleMap;

    const art::Ptr<recob::PFParticle> pfps = trackToNuPFParticleMap[ptrack];
    m_reco_track_is_nuslice[ntracks] = PFPToNuSliceMap[pfps];
    m_reco_track_sliceId[ntracks] = PFPToSliceIdMap[pfps];
    cout << PFPToNuSliceMap[pfps] << " " << PFPToSliceIdMap[pfps] << endl;
    */

    // skip cc muon track
    if (ptrack.key()==trkmuon.key()) cout << "HEY THIS IS CC MUON" << endl;
    if (ptrack.key()==trkmuon.key()) continue;

    cout << track.Vertex().X() << endl;
    // check track start and end
    TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());

    TVector3 end(track.End().X(),track.End().Y(),track.End().Z());

    reco_track_start_x[ntracks] = track.Vertex().X();
    reco_track_start_y[ntracks] = track.Vertex().Y();
    reco_track_start_z[ntracks] = track.Vertex().Z();

    reco_track_end_x[ntracks] = track.End().X();
    reco_track_end_y[ntracks] = track.End().Y();
    reco_track_end_z[ntracks] = track.End().Z();

    double st_vtx=TMath::Sqrt((reco_nu_vtx_x-pos.X())*(reco_nu_vtx_x-pos.X()) +
                              (reco_nu_vtx_y-pos.Y())*(reco_nu_vtx_y-pos.Y()) +
                              (reco_nu_vtx_z-pos.Z())*(reco_nu_vtx_z-pos.Z()));
    reco_track_distance[ntracks] = st_vtx;//distance between vtx and track start position

    reco_track_vtx_inTPC[ntracks] = isInsideVolume("TPC",pos);
    reco_track_vtx_in5cmTPC[ntracks] = isInsideVolume("5cmTPC",pos);
    reco_track_vtx_inCCInclusiveTPC[ntracks] = isInsideVolume("CCInclusiveTPC",pos);

    reco_track_end_inTPC[ntracks] = isInsideVolume("TPC",end);
    reco_track_end_in5cmTPC[ntracks] = isInsideVolume("5cmTPC",end);
    reco_track_end_inCCInclusiveTPC[ntracks] = isInsideVolume("CCInclusiveTPC",end);

    // track length and angles
    double trklen=track.Length();
    reco_track_length[ntracks] = trklen;
    reco_track_theta[ntracks] = track.Theta();
    reco_track_phi[ntracks] = track.Phi();

    reco_track_P_vtx[ntracks] = track.VertexMomentum();
    reco_track_P_str[ntracks] = track.StartMomentum();
    reco_track_P_end[ntracks] = track.EndMomentum();

    //TVector3 endmom = track.EndMomentumVector<TVector3>(); 
    //TVector3 strmom = track.StartMomentumVector<TVector3>(); 
    //TVector3 vtxmom = track.VertexMomentumVector<TVector3>(); 

    //auto trkf::TrackMomentumCalculator trkm;
    //cout << "mom calc: " << trkmom.GetTrackMomentum(trklen, 321) <<endl;
    //
    //double llhmom =  trkmom.GetMomentumMultiScatterLLHD(ptrack);

    // check kaon track start and end distance from vertex
    double end_dis=TMath::Sqrt((reco_nu_vtx_x-end.X())*(reco_nu_vtx_x-end.X()) +
                               (reco_nu_vtx_y-end.Y())*(reco_nu_vtx_y-end.Y()) +
                               (reco_nu_vtx_z-end.Z())*(reco_nu_vtx_z-end.Z()));
    reco_track_dir[ntracks] = (st_vtx<end_dis);

    //cout << "Track id " << ptrack.key();
    //cout << " X " << pos.X() << " - " << end.X();
    //cout << " Z " << pos.Z() << " - " << end.Z();
    //cout << " length " << trklen << endl;

    //fillCalorimetry(fmcal.at(ptrack.key()),track,assocMCPartt,ntracks);
    fillCalorimetry(fmcal.at(ptrack.key()),ntracks);

    /*
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
    */
  
    // check PID
    if (trackPIDAssn.isValid()){
      double delta_z = end.Z()-pos.Z();
      double delta_y = end.Y()-pos.Y();
      double angle_y = TMath::ATan2(delta_z, delta_y);
      fillPID(trackPIDAssn.at(ptrack.key()), angle_y, ntracks);
    }

    /*
      std::vector<art::Ptr<anab::ParticleID>> trackPID=trackPIDAssn.at(ptrack.key());
      if (trackPID.size()>0){
        double chi2ka[3] = {0,0,0};
        double chi2pr[3] = {0,0,0};
        double chi2pi[3] = {0,0,0};
        double chi2mu[3] = {0,0,0};
        std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
        for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
          anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
          if (AlgScore.fAlgName == "Chi2") {
            if (anab::kVariableType(AlgScore.fVariableType) == anab::kGOF) {
              for (int pl=0; pl<3; pl++) {
                if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==pl) {
                  if (AlgScore.fAssumedPdg==321)  chi2ka[pl]=AlgScore.fValue;
                  if (AlgScore.fAssumedPdg==2212) chi2pr[pl]=AlgScore.fValue;
                  if (AlgScore.fAssumedPdg==211)  chi2pi[pl]=AlgScore.fValue;
                  if (AlgScore.fAssumedPdg==13)   chi2mu[pl]=AlgScore.fValue;
                }
              }
            }
          }
          if (AlgScore.fAlgName == "ThreePlaneProtonPID") {
            if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood) {
              if (AlgScore.fPlaneMask==UBPID::uB_SinglePlaneGetBitset(2)) {
                if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward) {
                  if (AlgScore.fAssumedPdg==2212) reco_track_3pidpr[ntracks]=AlgScore.fValue;
                }
              }
            }
          } 
        }
        if (chi2ka[2]>0) reco_track_chi2ka[ntracks] = chi2ka[2];
        if (chi2pr[2]>0) reco_track_chi2pr[ntracks] = chi2pr[2];
        if (chi2pi[2]>0) reco_track_chi2pi[ntracks] = chi2pi[2];
        if (chi2mu[2]>0) reco_track_chi2mu[ntracks] = chi2mu[2];
        double delta_z = end.Z()-pos.Z();
        double delta_y = end.Y()-pos.Y();
        double theta0 = TMath::ATan2(delta_z,delta_y) - TMath::Pi()/3;
        double theta1 = TMath::ATan2(delta_z,delta_y) + TMath::Pi()/3.;
        double theta2 = TMath::ATan2(delta_z,delta_y);
        double wpl0 = TMath::Power(TMath::Sin(theta0),2) >= 0.05 ? 1 : 0;
        double wpl1 = TMath::Power(TMath::Sin(theta1),2) >= 0.05 ? 1 : 0;
        double wpl2 = TMath::Power(TMath::Sin(theta2),2) >= 0.05 ? 1 : 0;
        double chi2ka_3pl = (wpl0*chi2ka[0] + wpl1*chi2ka[1] + wpl2*chi2ka[2])/(wpl0 + wpl1 + wpl2);
        double chi2pr_3pl = (wpl0*chi2pr[0] + wpl1*chi2pr[1] + wpl2*chi2pr[2])/(wpl0 + wpl1 + wpl2);
        double chi2pi_3pl = (wpl0*chi2pi[0] + wpl1*chi2pi[1] + wpl2*chi2pi[2])/(wpl0 + wpl1 + wpl2);
        double chi2mu_3pl = (wpl0*chi2mu[0] + wpl1*chi2mu[1] + wpl2*chi2mu[2])/(wpl0 + wpl1 + wpl2);
        if (chi2ka_3pl>0) reco_track_chi2ka_3pl[ntracks] = chi2ka_3pl;
        if (chi2pr_3pl>0) reco_track_chi2pr_3pl[ntracks] = chi2pr_3pl;
        if (chi2pi_3pl>0) reco_track_chi2pi_3pl[ntracks] = chi2pi_3pl;
        if (chi2mu_3pl>0) reco_track_chi2mu_3pl[ntracks] = chi2mu_3pl;
      }
    }
    */

    // find true matched particle
      //std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>> assocMCPart;
    if (isMC) {
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);
      //assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> (new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, evt, fHitTruthAssns));
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
      fillTrueMatching(hits_from_track, particles_per_hit, ntracks);
      fillTrackMatching(hits_from_track, particles_per_hit, ntracks);
      /*
      simb::MCParticle const* matched_mcparticle = NULL;
      std::unordered_map<int,double> trkide;
      double maxe=-1, tote=0;
      std::vector<simb::MCParticle const*> particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
    
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
        reco_track_true_pdg[ntracks] = matched_mcparticle->PdgCode();
        const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
        reco_track_true_origin[ntracks]=int(mc_truth->Origin());
        reco_track_true_primary[ntracks]=matched_mcparticle->Process()=="primary";
        double x = matched_mcparticle->EndX();
        double y = matched_mcparticle->EndY();
        double z = matched_mcparticle->EndZ();
        reco_track_true_end_inTPC[ntracks] = isInsideVolume("TPC",x,y,z);
        reco_track_true_end_in5cmTPC[ntracks] = isInsideVolume("5cmTPC",x,y,z);
        reco_track_true_end_inCCInclusiveTPC[ntracks] = isInsideVolume("CCInclusiveTPC",x,y,z);
        TLorentzVector mcstart, mcend;
        unsigned int pstarti, pendi;
        reco_track_true_length[ntracks] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
      }
      */

    }//isMC
  
    // check if there is a track at the end
    int ndaughters = 0;
    int ndaughters_sh = 0;

    double dgvecX_tr=0;
    double dgvecY_tr=0;
    double dgvecZ_tr=0; 
  
    //cout << "NTracks: " << NTracks << endl;
    for (int j=0; j<NTracks; j++) {
    
      art::Ptr<recob::Track> ptrack2(trackListHandle,j);
      const recob::Track& track2 = *ptrack2;
    
      // skip all primary tracks
      
      bool skip = false;
      for (int k=0; k<reco_nu_ndaughters; k++) {
        if (int(ptrack2.key())==reco_nu_daughters_id[k]) {
          skip=true;
          break;
        }
      }
      if (skip) continue;
      

      TVector3 pos2(track2.Vertex().X(),track2.Vertex().Y(),track2.Vertex().Z());

      double track2_distance=TMath::Sqrt((end.X()-pos2.X())*(end.X()-pos2.X()) +
                                         (end.Y()-pos2.Y())*(end.Y()-pos2.Y()) +
                                         (end.Z()-pos2.Z())*(end.Z()-pos2.Z()));

      //cout << "distance to vertex: " << track2_distance << endl;

      // check distance to vertex
      //if (track2_distance<10) { //7cm, 
      if (track2_distance<40) { //7cm, 
      //if (track2_distance<30) { //7cm, 
      //if (track2_distance<100) { //7cm,
	cout << "this is daughter track" << endl;
        reco_track_daughter_distance[ntracks][ndaughters] = track2_distance;

        TVector3 end2(track2.End().X(),track2.End().Y(),track2.End().Z());

	reco_track_daughter_start_x[ntracks][ndaughters] = track2.Vertex().X();
	reco_track_daughter_start_y[ntracks][ndaughters] = track2.Vertex().Y();
	reco_track_daughter_start_z[ntracks][ndaughters] = track2.Vertex().Z();
	reco_track_daughter_end_x[ntracks][ndaughters] = track2.End().X();
	reco_track_daughter_end_y[ntracks][ndaughters] = track2.End().Y();
	reco_track_daughter_end_z[ntracks][ndaughters] = track2.End().Z();

        reco_track_daughter_vtx_inTPC[ntracks][ndaughters] = isInsideVolume("TPC",pos2);
        reco_track_daughter_vtx_in5cmTPC[ntracks][ndaughters] = isInsideVolume("5cmTPC",pos2);
        reco_track_daughter_vtx_inCCInclusiveTPC[ntracks][ndaughters] = isInsideVolume("CCInclusiveTPC",pos2);
      
        reco_track_daughter_end_inTPC[ntracks][ndaughters] = isInsideVolume("TPC",end2);
        reco_track_daughter_end_in5cmTPC[ntracks][ndaughters] = isInsideVolume("5cmTPC",end2);
        reco_track_daughter_end_inCCInclusiveTPC[ntracks][ndaughters] = isInsideVolume("CCInclusiveTPC",end2);

        // track length and angles
        reco_track_daughter_length[ntracks][ndaughters] = track2.Length();
        reco_track_daughter_theta[ntracks][ndaughters] = track2.Theta();
        reco_track_daughter_phi[ntracks][ndaughters] = track2.Phi();

        double start_dis=TMath::Sqrt((reco_nu_vtx_x-pos2.X())*(reco_nu_vtx_x-pos2.X()) +
                                     (reco_nu_vtx_y-pos2.Y())*(reco_nu_vtx_y-pos2.Y()) +
                                     (reco_nu_vtx_z-pos2.Z())*(reco_nu_vtx_z-pos2.Z()));    

        reco_track_daughter_vtx_distance[ntracks][ndaughters] = start_dis;

        double trvecX = end.X() - pos.X();
        double trvecY = end.Y() - pos.Y();
        double trvecZ = end.Z() - pos.Z();

        double dgvecX = end2.X() - pos2.X();
        double dgvecY = end2.Y() - pos2.Y();
        double dgvecZ = end2.Z() - pos2.Z();

        dgvecX_tr = dgvecX;
        dgvecY_tr = dgvecY;
        dgvecZ_tr = dgvecZ;

        double nortr=TMath::Sqrt( trvecX*trvecX + trvecY*trvecY + trvecZ*trvecZ  );
        double nordg=TMath::Sqrt( dgvecX*dgvecX + dgvecY*dgvecY + dgvecZ*dgvecZ  );

        double cosTh = (trvecX*dgvecX + trvecY*dgvecY + trvecZ*dgvecZ)/(nortr*nordg);

        double angTrDg = acos(cosTh);

        reco_angle_track_daughter[ntracks][ndaughters] = angTrDg;
        cout << " ---------------- Angle btw track and daughter: " << angTrDg << endl;

        //cout << "Daughter id " << ptrack2.key();
        //cout << " X " << pos2.X() << " - " << end2.X();
        //cout << " Z " << pos2.Z() << " - " << end2.Z();
        //cout << " length " << track2.Length() << endl;

        //fillCalorimetry(fmcal.at(ptrack2.key()),track2,assocMCPartt,ntracks,ndaughters);
        fillCalorimetry(fmcal.at(ptrack2.key()),ntracks,ndaughters);
/*
        cout << "Tracks: " << ntracks << " Daughters: " << ndaughters <<endl;

        cout << "++++++++++ Size of tracks: " << reco_track_daughter_ResRan.size() <<endl;

        for(unsigned int m=0; m<reco_track_daughter_ResRan.size(); m++){
         
            cout << "+++++++++++++ Size of daughters: " << reco_track_daughter_ResRan[m].size() << endl;
            for(unsigned int n=0; n<reco_track_daughter_ResRan[m].size(); n++){
            
               cout << "\n================ Size of hits: " <<reco_track_daughter_ResRan[m][n].size() << endl;
            
            }
        
        }

        
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
        */

        // check PID
        if (trackPIDAssn.isValid()) {
          double delta_z = end2.Z()-pos2.Z();
          double delta_y = end2.Y()-pos2.Y();
          double angle_y = TMath::ATan2(delta_z, delta_y);
          fillPID(trackPIDAssn.at(ptrack2.key()), angle_y, ntracks, ndaughters);
        }
      
        // find true matched particle
        if (isMC) {
          art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);
          std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack2.key());
	  art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
	  hits_from_reco.insert(hits_from_reco.end(), hits_from_track.begin(), hits_from_track.end());

          fillTrueMatching(hits_from_track, particles_per_hit, ntracks, ndaughters);
  	  fillTrackMatching(hits_from_track, particles_per_hit, ntracks, ndaughters);



	  std::vector<simb::MCParticle const*> particle_vec;
	  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
	  for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
	    particle_vec.clear();
	    match_vec.clear();
	    particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
	    
	    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
	      if(particle_vec[i_p]->PdgCode() == 211) hits_from_true_pip.push_back(hits_from_track[i_h]);
	      if(particle_vec[i_p]->PdgCode() == -13) hits_from_true_mup.push_back(hits_from_track[i_h]);
	    }
	  }


	  //if(true_kaon_end_process == 0 || true_kaon_end_process == 1){
	  //if(reco_track_true_pdg[ntracks] == 321){
	      fillAngularDistributionMap(hits_from_track, end, spacepoint_per_hit, angular_distribution_map_track);
	      fillAngularDistributionMap3D(hits_from_track, end, spacepoint_per_hit, angular_distribution_map_3D_track);
	      //smoothAngularDistributionMap(angular_distribution_map_track); 
	      fillAngularDistributionMapCheat(hits_from_track, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map, mrgidpdg, currentMergedId);
	      fillAngularDistributionMapCheat3D(hits_from_track, true_pi0_hit_list, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map_3D, mrgidpdg_3D, mrgidmom, currentMergedId_3D);

	      //currentMergedId = fillAngularDistributionMapCheat(hits_from_track, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map, mrgidpdg, currentMergedId);
	      //  }
	      //}
	 
        }//isMC
      
        ndaughters++;

      } // remove condition on gap between K+ end and daughter start
    
      reco_track_ndaughters[ntracks] = ndaughters;
          
      if(prim_k_id>0 && IsKaon == true){
	if(reco_track_true_pdg[ntracks] != 321){
	  flg_reco_k = true;
	  if(reco_track_daughter_true_pdg[ntracks][0] == -13) flg_reco_mu = true; 
	  if(reco_track_daughter_true_pdg[ntracks][0] == 211) flg_reco_pi = true; 
	}
      }
	
    }//NTrack loop
  
    
    for (int j_s=0; j_s<NShowers; j_s++) {
      art::Ptr<recob::Shower> pshower2(showerListHandle,j_s);
      const recob::Shower& shower2 = *pshower2;
      
      // skip all primary showers                                                                                                                           
      bool skip2 = false;
      for (int k_sh=0; k_sh<reco_nu_ndaughters; k_sh++) {
	if (int(pshower2.key())==reco_nu_daughters_id[k_sh]) {
	  skip2=true;
	  break;
	}
      }
      if (skip2) continue;
      TVector3 pos2_sh(shower2.ShowerStart().X(),shower2.ShowerStart().Y(),shower2.ShowerStart().Z());

      double track2_distance_sh=TMath::Sqrt((end.X()-pos2_sh.X())*(end.X()-pos2_sh.X()) +
                                            (end.Y()-pos2_sh.Y())*(end.Y()-pos2_sh.Y()) +
                                            (end.Z()-pos2_sh.Z())*(end.Z()-pos2_sh.Z()));

      //cout << "distance to shower vertex: " << track2_distance_sh << endl;

      // check distance to vertex
      //if (ndaughters_sh<=50 && track2_distance_sh<100) { //50cm	-> 15m cm
      if (ndaughters_sh<=20 && track2_distance_sh<40) { //50cm -> 15m cm
      //if (ndaughters_sh<=20 && track2_distance_sh<30) { //50cm -> 15m cm

        reco_track_daughter_distance_sh[ntracks][ndaughters_sh] = track2_distance_sh;

        //TVector3 end2_sh(shower2.Direction().X(),shower2.Direction().Y(),shower2.Direction().Z());
        TVector3 end2_sh(shower2.ShowerStart().X() + shower2.Length() * shower2.Direction().X(),
                         shower2.ShowerStart().Y() + shower2.Length() * shower2.Direction().Y(),
                         shower2.ShowerStart().Z() + shower2.Length() * shower2.Direction().Z());

	//ls.OpenAngle(); shr->dEdx()[0])  shr->Energy()[1])
	//cout << "SHOWER ENERGY: " << shower2.Energy()[2] << endl;
	/*
	reco_track_daughter_open_angle_sh[ntracks][ndaughters_sh] = shower2.OpenAngle();
	reco_track_daughter_dedx_pl0_sh[ntracks][ndaughters_sh] = shower2.dEdx()[0];
	reco_track_daughter_dedx_pl1_sh[ntracks][ndaughters_sh] = shower2.dEdx()[1];
	reco_track_daughter_dedx_pl2_sh[ntracks][ndaughters_sh] = shower2.dEdx()[2];
	reco_track_daughter_energy_pl0_sh[ntracks][ndaughters_sh] = shower2.Energy()[0];
	reco_track_daughter_energy_pl1_sh[ntracks][ndaughters_sh] = shower2.Energy()[1];
	reco_track_daughter_energy_pl2_sh[ntracks][ndaughters_sh] = shower2.Energy()[2];
	*/

        reco_track_daughter_vtx_inTPC_sh[ntracks][ndaughters_sh] = isInsideVolume("TPC",pos2_sh);
        reco_track_daughter_vtx_in5cmTPC_sh[ntracks][ndaughters_sh] = isInsideVolume("5cmTPC",pos2_sh);
        reco_track_daughter_vtx_inCCInclusiveTPC_sh[ntracks][ndaughters_sh] = isInsideVolume("CCInclusiveTPC",pos2_sh);

        reco_track_daughter_end_inTPC_sh[ntracks][ndaughters] = isInsideVolume("TPC",end2_sh);
        reco_track_daughter_end_in5cmTPC_sh[ntracks][ndaughters] = isInsideVolume("5cmTPC",end2_sh);
        reco_track_daughter_end_inCCInclusiveTPC_sh[ntracks][ndaughters] = isInsideVolume("CCInclusiveTPC",end2_sh);

        // track length and angles
        reco_track_daughter_length_sh[ntracks][ndaughters_sh] = shower2.Length();
        reco_track_daughter_theta_sh[ntracks][ndaughters_sh] = shower2.Direction().Theta();
        reco_track_daughter_phi_sh[ntracks][ndaughters_sh] = shower2.Direction().Phi();
	//reco_track_daughter_true_pdg_sh[ntracks][ndaughters_sh] = shower2.PdgCode();

        double start_dis_sh=TMath::Sqrt((reco_nu_vtx_x-pos2_sh.X())*(reco_nu_vtx_x-pos2_sh.X()) +
                                        (reco_nu_vtx_y-pos2_sh.Y())*(reco_nu_vtx_y-pos2_sh.Y()) +
                                        (reco_nu_vtx_z-pos2_sh.Z())*(reco_nu_vtx_z-pos2_sh.Z()));    

        reco_track_daughter_vtx_distance_sh[ntracks][ndaughters_sh] = start_dis_sh;

	double trvecX = end.X() - pos.X();
        double trvecY = end.Y() - pos.Y();
        double trvecZ = end.Z() - pos.Z();

        double dgvecX_sh = end2_sh.X() - pos2_sh.X();
        double dgvecY_sh = end2_sh.Y() - pos2_sh.Y();
        double dgvecZ_sh = end2_sh.Z() - pos2_sh.Z();

        double nortr=TMath::Sqrt( trvecX*trvecX + trvecY*trvecY + trvecZ*trvecZ  );
        double nordg_tr=TMath::Sqrt( dgvecX_tr*dgvecX_tr + dgvecY_tr*dgvecY_tr + dgvecZ_tr*dgvecZ_tr  );
        double nordg_sh=TMath::Sqrt( dgvecX_sh*dgvecX_sh + dgvecY_sh*dgvecY_sh + dgvecZ_sh*dgvecZ_sh  );

        double cosTh_sh = (trvecX*dgvecX_sh + trvecY*dgvecY_sh + trvecZ*dgvecZ_sh)/(nortr*nordg_sh);
        double cosTh_tr_sh = (dgvecX_tr*dgvecX_sh + dgvecX_tr*dgvecY_sh + dgvecX_tr*dgvecZ_sh)/(nordg_tr*nordg_sh);

        double angTrDg_sh = acos(cosTh_sh);
        double angTrDg_tr_sh = acos(cosTh_tr_sh);

        reco_angle_track_daughter_sh[ntracks][ndaughters_sh] = angTrDg_sh;
        if(ndaughters==1) reco_angle_daughter_track_daughter_sh[ntracks][ndaughters_sh] = angTrDg_tr_sh;
        cout << " ---------------- Angle btw track and shower: " << angTrDg_sh << endl;

        // find true matched particle
        if (isMC) {
	  //cout << "find true matched particle of showers" << endl;
          art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);
          std::vector<art::Ptr<recob::Hit>> hits_from_shower = hits_from_showers.at(pshower2.key()); // can it be hits_from_tracks?
	  art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
	  hits_from_reco.insert(hits_from_reco.end(), hits_from_shower.begin(), hits_from_shower.end()); 

          fillTrueMatching_sh(hits_from_shower, particles_per_hit, ntracks, ndaughters_sh);
          fillShowerMatching(hits_from_shower, particles_per_hit, ntracks, ndaughters_sh);


	  std::vector<simb::MCParticle const*> particle_vec;
	  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
	  for(size_t i_h=0; i_h<hits_from_shower.size(); i_h++){
	    particle_vec.clear();
	    match_vec.clear();
	    particles_per_hit.get(hits_from_shower[i_h].key(),particle_vec,match_vec);
	    
	    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
	      if(particle_vec[i_p]->PdgCode() == 211) hits_from_true_pip.push_back(hits_from_shower[i_h]);
	      if(particle_vec[i_p]->PdgCode() == -13) hits_from_true_mup.push_back(hits_from_shower[i_h]);
	    }
	  }


	  // if(true_kaon_end_process == 0 || true_kaon_end_process == 1){
	  //if(reco_track_true_pdg[ntracks] == 321){
	      
	      fillAngularDistributionMap(hits_from_shower, end, spacepoint_per_hit, angular_distribution_map_shower);
	      fillAngularDistributionMap3D(hits_from_shower, end, spacepoint_per_hit, angular_distribution_map_3D_shower);
	      //smoothAngularDistributionMap(angular_distribution_map_shower);
	      fillAngularDistributionMapCheat(hits_from_shower, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map, mrgidpdg, currentMergedId);
	      fillAngularDistributionMapCheat3D(hits_from_shower, true_pi0_hit_list, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map_3D, mrgidpdg_3D, mrgidmom, currentMergedId_3D);
	      //currentMergedId = fillAngularDistributionMapCheat(hits_from_shower, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map, mrgidpdg, currentMergedId);

	      //   }
	      // }

          //fillAngularDistributionMapCheat(hits_from_shower, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map, mrgidpdg, currentMergedId);
	  //fillAngularDistributionMap(hits_from_shower, end, spacepoint_per_hit, angular_distribution_map_shower);
	  //cout <<"angular_distribution_map_shower.size() inside: " << angular_distribution_map_shower.size() << endl;
	  //fillHistAngularDistributionMap(angular_distribution_map,  h_angular_dirtribution_pfparticle);
	  //cout << "after fillTrue the pdg is: " << reco_track_daughter_true_pdg_sh[ntracks][ndaughters_sh] << endl;
        }

	ndaughters_sh++;

      }
    
      reco_shower_ndaughters[ntracks] = ndaughters_sh;
    }
  
    smoothAngularDistributionMapCheat3D(angular_distribution_mrgid_map_3D);
    smoothAngularDistributionMap3D(angular_distribution_map_3D_track);
    smoothAngularDistributionMap3D(angular_distribution_map_3D_shower);

    accumulateAngularDistributionMap3D(angular_distribution_map_3D_track, angular_distribution_map_3D_shower, angular_distribution_map_3D_pfparticle);

    obtainPeakVectorCheat3D(angular_distribution_mrgid_map_3D, mrgidpdg_3D, view_peak_map_cheat, v_pdg_peak);
    obtainPeakVector3D(angular_distribution_map_3D_track, v_trk_flg_peak, view_peak_map, true);
    obtainPeakVector3D(angular_distribution_map_3D_shower, v_trk_flg_peak, view_peak_map, false);

    findBestAngularPeakCheat3D(view_peak_map_cheat, best_peak_bins_cheat);
    findBestAngularPeak3D(angular_distribution_map_3D_pfparticle, view_peak_map, best_peak_bins);

    art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);

    cout << "reco_track_true_pdg is " << reco_track_true_pdg[ntracks] << ", and process is" << true_kaon_end_process << endl;
    cout << "findshowerspine3D for RECO" << endl;
    findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_vector, end, view_peak_map, best_peak_bins);

    //if(best_peak_bins.size()!=0) findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_vector, end, view_peak_map, best_peak_bins);
    //if(best_peak_bins.size()!=0) findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list, end, view_peak_map, best_peak_bins);


    //hits_from_true_pip


    if(best_peak_bins_cheat.size()!=0){
      int i=0;
      for(auto& entry : best_peak_bins_cheat){
	//cout << "PEAK OF PDG IS " << entry.first << endl;
	if(reco_track_true_pdg[ntracks] == 321){
	  if(entry.first == 211){
	    
	    cout << "findshowerspine3D for PI+" << endl;
	    cout << "true_track_length of pion: " << true_dau_pip_length << endl;
	    findShowerSpine3D(hits_from_true_pip, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_pip_vector, end, view_peak_map, entry.second);
	    //findShowerSpine3D(hits_from_true_pip, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_pip, end, view_peak_map, entry.second);
	    //findShowerSpine3D(hits_from_true_pip_vector, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_pip, end, view_peak_map, entry.second);
	    
	    //trackRebuid(shower_spine_hit_list_cheat_pip, spacepoint_per_hit, track); 

	    //findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_pip, end, view_peak_map, entry.second);
	  }
	  if(entry.first == -13){
	    cout << "findshowerspine3D for MU+" << endl;
	    cout << "true_track_length of muon: " << true_dau_muon_length << endl;
	    findShowerSpine3D(hits_from_true_mup, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_mup_vector, end, view_peak_map, entry.second);
	    //findShowerSpine3D(hits_from_true_mup, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_mup, end, view_peak_map, entry.second);
	    //findShowerSpine3D(hits_from_true_mup_vector, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_mup, end, view_peak_map, entry.second);
	    //findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_mup, end, view_peak_map, entry.second);
	  }
	}
	//findShowerSpine3D(true_pi0_hit_list, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat, end, view_peak_map, entry.second);
	//findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat, end, view_peak_map, entry.second);

	//findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat, end, view_peak_map, entry.second, true_dau_pip_length);


	findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_cheat_vector, end, view_peak_map, entry.second, cheat_num_hits[ntracks][i], cheat_ini_hits[ntracks][i], cheat_closest_distance[ntracks][i]);

	i++;
      }
    }


    //cout << "shower_spine_hit_list.size() " << shower_spine_hit_list.size() << endl;
 
    //for(auto& shower_spine_hit_list : shower_spine_hit_list_vector){
      //if(shower_spine_hit_list.size()!=0) trackRebuid(shower_spine_hit_list, spacepoint_per_hit, track);
    //}

    //if(shower_spine_hit_list.size()!=0) obtainLongitudinalDecomposition(shower_spine_hit_list, spacepoint_per_hit);   
    //if(shower_spine_hit_list.size()!=0) obtainLength(shower_spine_hit_list, spacepoint_per_hit, end);
    //if(shower_spine_hit_list.size()!=0) trackRebuid(shower_spine_hit_list, spacepoint_per_hit, track); 

    /*
    if(shower_spine_hit_list_cheat.size()!=0) obtainLongitudinalDecomposition(shower_spine_hit_list_cheat, spacepoint_per_hit);
    if(shower_spine_hit_list_cheat.size()!=0) obtainLength(shower_spine_hit_list_cheat, spacepoint_per_hit, end);
    if(shower_spine_hit_list_cheat.size()!=0) trackRebuid(shower_spine_hit_list_cheat, spacepoint_per_hit, track); 
    */

    /*
    if(true_pi0_hit_list.size()!=0) obtainLongitudinalDecomposition(true_pi0_hit_list, spacepoint_per_hit);
    if(true_pi0_hit_list.size()!=0) obtainLength(true_pi0_hit_list, spacepoint_per_hit, end);
    if(true_pi0_hit_list.size()!=0) trackRebuid(true_pi0_hit_list, spacepoint_per_hit, track);
    */

    //cout << "true_track_length of muon: " << true_dau_muon_length << endl;
    //cout << "true_track_length of pip: " << true_dau_pip_length << endl;


    if(best_peak_bins_cheat.size()!=0){
      int i=0;
      for(auto const& entry : best_peak_bins_cheat){
	for(auto const& entry2 : entry.second){
	//if(entry.)
	//for(auto const& entry2 : entry.second){
	  //entry.second.begin()->second
	  cheat_peak_pdg[ntracks][i] = entry.first;
	  //cheat_peak_height[ntracks][i] = 
	  cheat_peak_theta[ntracks][i] = entry2.X() * thetaBinSize;
	  cheat_peak_phi[ntracks][i] = entry2.Y() * phiBinSize;
	  //cheat_peak_theta[ntracks][i] = entry.second.begin()->second.X() * thetaBinSize;
	  //cheat_peak_phi[ntracks][i] = entry.second.begin()->second.Y() * phiBinSize;
	  //}
	i++;
	}
      }
    }

    
    cout << "hits_from_true_pip.size(): " << hits_from_true_pip.size() << endl;
    cout << "hits_from_true_mup.size(): " << hits_from_true_mup.size() << endl;

    //if(hits_from_true_pip.size()<20) cheat_pip_trkln = -9999;
    //if(shower_spine_hit_list_cheat_pip.size()==0) cheat_pip_trkln = -9999;
    if(shower_spine_hit_list_cheat_pip_vector.size()==0 || shower_spine_hit_list_cheat_pip_vector[0].size()==0) cheat_pip_trkln = -9999;
    else{
      cheat_pip_trkln = trackRebuid(hits_from_true_pip, spacepoint_per_hit, track);
      cout << "full cheat pip trkln " << cheat_pip_trkln << endl;
     cheat_pip_trkln = trackRebuid(shower_spine_hit_list_cheat_pip_vector[0], spacepoint_per_hit, track);
      cout << "half cheat pip trkln " << cheat_pip_trkln << endl;
    }

    //if(hits_from_true_mup.size()<20) cheat_mup_trkln = -9999;
    //if(shower_spine_hit_list_cheat_mup.size()==0) cheat_mup_trkln = -9999;
    if(shower_spine_hit_list_cheat_mup_vector.size()==0 || shower_spine_hit_list_cheat_mup_vector[0].size()==0) cheat_mup_trkln = -9999;
    else{
      cheat_mup_trkln = trackRebuid(hits_from_true_mup, spacepoint_per_hit, track);
      cout << "full cheat mup trkln " << cheat_mup_trkln << endl;
      cheat_mup_trkln = trackRebuid(shower_spine_hit_list_cheat_mup_vector[0], spacepoint_per_hit, track);
      //cheat_mup_trkln = trackRebuid(shower_spine_hit_list_cheat_mup, spacepoint_per_hit, track);
      cout << "half cheat mup trkln " << cheat_mup_trkln << endl; 
    }

    art::Handle< std::vector<recob::Track> > trackListHandle_new;
    std::vector<art::Ptr<recob::Track> > tracklist_new;
    if(evt.getByLabel(fTrackModuleLabel,trackListHandle_new)){ 
      art::fill_ptr_vector(tracklist_new, trackListHandle_new);
    }
    
    //const std::vector<recob::Track> * tracklistvec = trackListHandle_new.product();
    //const_cast<std::vector<recob::Track>*>(tracklistvec) = trackListHandle_new.product();


    //std::vector<recob::Track> * tracklistvec = const_cast<std::vector<recob::Track>*>(trackListHandle_new.product());

    /*
    std::type_info const& typeInfoOfWrapper{typeid(art::Wrapper< std::vector<recob::Track> >)};
    gallery::Event gal_evt = evt;
    auto res = evt.getByLabel(typeInfoOfWrapper, fTrackModuleLabel);
    */

    double trkln;
    //if(true_kaon_end_process==0) cout << "best_peak_bins.size() of PI+ decay is " << best_peak_bins.size() << endl;
    if(best_peak_bins.size()!=0){
	int i=0;
      for(auto const& entry : best_peak_bins){
	best_peak_theta[ntracks][i] = entry.X() * thetaBinSize;
	best_peak_phi[ntracks][i] = entry.Y() * phiBinSize;

	
	if(shower_spine_hit_list_vector.size()!=0){
	  //if(shower_spine_hit_list_vector[i].size()==0 || shower_spine_hit_list_vector[i].size() > 3000 || shower_spine_hit_list_vector[i].size()<10) continue;
	  if( shower_spine_hit_list_vector[i].size()==0 ) continue;
	  //cout << "shower_spine_hit_list_vector[i]: " << shower_spine_hit_list_vector[i].size() << endl;
	  cout << "CCKAON PROCESS IS " << true_kaon_end_process << endl;
	  trkln = trackRebuid(shower_spine_hit_list_vector[i], spacepoint_per_hit, track);

	  //recob::Track reco_track = trackRebuid2(shower_spine_hit_list_vector[i], spacepoint_per_hit, track);
	  //cout << "BEFORE push_back tracklistvec.size(): " << tracklistvec->size() << endl;   
	  //tracklistvec->push_back(reco_track);
	  //cout << "AFTER push_back tracklistvec.size(): " << tracklistvec->size() << endl;  

	  cout << "trackRebuid: " << trkln << endl;
	  best_peak_trkln[ntracks][i] = trkln;
	  //best_peak_trkln[ntracks][i] = trackRebuid(shower_spine_hit_list_vector[i], spacepoint_per_hit, track);
	}
	
	i++;
      }
    }
    if(true_kaon_end_process==0) cout << "true_track_length of muon: " << true_dau_muon_length << endl;
    if(true_kaon_end_process==1) cout << "true_track_length of pip: " << true_dau_pip_length << endl;


    //trackListHandle_new = art::Handle< std::vector<recob::Track> >{tracklistvec, fTrackModuleLabel};

    //art::Provenance prov = &const_cast<art::Provenance*>(trackListHandle_new.provenance());

    //art::Provenance * prov = const_cast<art::Provenance*>(trackListHandle_new.provenance());
    //art::ValidHandle< std::vector<recob::Track> > trackListHandle_valid(tracklistvec, *prov );

    //art::ValidHandle< std::vector<recob::Track> > trackListHandle_valid(tracklistvec, const_cast<art::Provenance*>(trackListHandle_new.provenance()) );
    //art::ValidHandle< std::vector<recob::Track> > trackListHandle_valid = evt.getValidHandle< std::vector<recob::Track> >(fTrackModuleLabel);

    //trackListHandle_new = evt.DataViewImpl::getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
    //getHandle<PROD>(tag);


    //art::FindManyP<anab::Calorimetry> fmcal_new(trackListHandle_valid, evt, fCalorimetryModuleLabel);
    //cout << "tracklistvec.size()" << tracklistvec->size() << endl;
    //if(!fmcal_new.isValid()) cout << "fmcal_new is INVALID" << endl;

    /*
    int ndaughters_new=0;
    for (int j=0; j<NTracks; j++) {
    
      art::Ptr<recob::Track> ptrack2(trackListHandle_valid,j);
      const recob::Track& track2 = *ptrack2;
   
      TVector3 pos2(track2.Vertex().X(),track2.Vertex().Y(),track2.Vertex().Z());

      double track2_distance=TMath::Sqrt((end.X()-pos2.X())*(end.X()-pos2.X()) +
                                         (end.Y()-pos2.Y())*(end.Y()-pos2.Y()) +
                                         (end.Z()-pos2.Z())*(end.Z()-pos2.Z()));

      cout << "ndaughters_new: " << ndaughters_new << ", track2_distance: " << track2_distance << endl;
      fillCalorimetry(fmcal_new.at(ptrack2.key()),ntracks,ndaughters_new);
      //if (track2_distance<300) fillCalorimetry(fmcal_new.at(ptrack2.key()),ntracks,ndaughters_new);
      ndaughters_new++;
    }
*/

    /*
    fillHistAngularDistributionMapCheat(angular_distribution_mrgid_map, mrgidpdg, h_angular_distribution_pfparticle_cheat, v_pdg);
    fillHistAngularDistributionMapCheat3D(angular_distribution_mrgid_map_3D, mrgidpdg_3D, h_angular_distribution_pfparticle_cheat_3D, v_pdg_3D);
    fillHistAngularDistributionMap(angular_distribution_map_track,  h_angular_distribution_pfparticle, v_trk_flg, true);
    fillHistAngularDistributionMap(angular_distribution_map_shower,  h_angular_distribution_pfparticle, v_trk_flg, false);

    fillHistAngularDistributionMap3D(angular_distribution_map_3D_track,  h_angular_distribution_pfparticle_3D, v_trk_flg_3D, true); 
    fillHistAngularDistributionMap3D(angular_distribution_map_3D_shower,  h_angular_distribution_pfparticle_3D, v_trk_flg_3D, false);
    */

    //fillHistAngularDistributionMapSurface(angular_distribution_map_3D_track,  h_angular_distribution_pfparticle_surface, v_trk_flg_surface, true);
    //fillHistAngularDistributionMapSurface(angular_distribution_map_3D_shower,  h_angular_distribution_pfparticle_surface, v_trk_flg_surface, true); 
    //fillHistAngularDistributionMapSphere(angular_distribution_map_3D_track,  h_angular_distribution_pfparticle_sphere, v_trk_flg_3D, true); 
    //fillHistAngularDistributionMapSphere(angular_distribution_map_3D_shower,  h_angular_distribution_pfparticle_sphere, v_trk_flg_3D, false); 


    /*
    drawHistAngularDistributionMapCheat(h_angular_distribution_pfparticle_cheat, v_pdg, "cheat_angle_distribution.root", c);
    drawHistAngularDistributionMap(h_angular_distribution_pfparticle, v_trk_flg, "reco_angle_distribution.root", c);
    */
    /*
    fillHistAngularDistributionMapCheat3D(angular_distribution_mrgid_map_3D, mrgidpdg_3D, h_angular_distribution_pfparticle_cheat_3D, v_pdg_3D);
    fillHistAngularDistributionMap3D(angular_distribution_map_3D_track,  h_angular_distribution_pfparticle_3D, v_trk_flg_3D, true); 
    fillHistAngularDistributionMap3D(angular_distribution_map_3D_shower,  h_angular_distribution_pfparticle_3D, v_trk_flg_3D, false);
    drawHistAngularDistributionMapCheat3D(h_angular_distribution_pfparticle_cheat_3D, v_pdg_3D, "cheat_angle_distribution_3D.root", c);
    drawHistAngularDistributionMap3D(h_angular_distribution_pfparticle_3D, v_trk_flg_3D, "reco_angle_distribution_3D.root", c);
    */
    /*
    drawHistAngularDistributionMapSurface(h_angular_distribution_pfparticle_surface, v_trk_flg_surface, "reco_angle_distribution_surface.root", c);
    drawHistAngularDistributionMapSphere(h_angular_distribution_pfparticle_sphere, v_trk_flg_3D, "reco_angle_distribution_sphere.root", c);
    */

    mrgidpdg.clear();
    angular_distribution_mrgid_map.clear();
    angular_distribution_map_track.clear();
    angular_distribution_map_shower.clear();
    h_angular_distribution_pfparticle_cheat.clear();
    h_angular_distribution_pfparticle.clear();
    h_angular_distribution_pfparticle_3D.clear();
    v_trk_flg.clear();
    v_trk_flg_3D.clear();
    v_pdg.clear();
    v_pdg_3D.clear();
   
    ntracks++;
  
    if (reco_track_end_inTPC && reco_track_vtx_inTPC) {

      int hits_p2 = reco_track_nhits2[ntracks];
      if(hits_p2>=5){
        if(trklen>=5){ // initiall this was 5 cm
          if (st_vtx<=7 && (st_vtx<end_dis)) { // 2.5 // 5
            //save track ID for kaon candidate
            kaon_can_trkID.push_back(ptrack.key());
            cut_11=1;
            //--cout << "kaon track candidate " << ptrack.key() << endl;
          } // vtx distance to the K track < 2.5
        } // trklen > 7 cm
      } // more than 5 hits in the collection plane
    }// track is contained

    //} // vertex distance
  } // loop over K  reco trks	

  reco_ntracks = ntracks;

  if(flg_reco_k == true){
    cout << "Event has reconstrcted K+ track" << endl;
    //if(true_kaon_end_process == 0){
    //cout << "K+ -> Mu+ decay" << endl;
      if(flg_reco_mu == true) cout << "Event has reconstrcted mu+ track" << endl;
      else cout << "Event DOES NOT have reconstrcted mu+ track" << endl;
      //}
    //if(true_kaon_end_process == 1){
      //cout << "K+ -> Pi+ decay" << endl;
      if(flg_reco_pi == true) cout << "Event has reconstrcted pi+ track" << endl; 
      else cout << "Event DOES NOT have reconstrcted pi+ track" << endl;
      //}
  }
  else cout << "Event DO NOT have reconstrcted K+ track" << endl;

  flg_reco_k  = false;
  flg_reco_mu  = false;
  flg_reco_pi  = false;

  /*
  for (int i=0; i<reco_nu_ndaughters; i++) {
    for (int j=0; j<NShowers; j++) {
      art::Ptr<recob::Shower> pshower2(showerListHandle,j);
      //const recob::Track& shower2 = *pshower2;
      
      // skip all primary showers                                                                                                                                                                                                                                                                                                   
      bool skip2 = false;
      for (int k=0; k<reco_nu_ndaughters; k++) {
	if (int(pshower2.key())==reco_nu_daughters_id[k]) {
	  skip2=true;
	  break;
	}
      }
      if (skip2) continue;
      ndaughters_sh++;
      reco_shower_ndaughters[ntracks] = ndaughters_sh;
    }
    nshowers++;
  }

  reco_nshowers = nshowers;
  */

  //--cout << "Number of kaon candidate tracks " << reco_nu_ndaughters << endl;

  // find true matched particle for ccmuon
  if (isMC) {

    simb::MCParticle const* matched_mcparticle = NULL;
    std::unordered_map<int,double> trkide;
    double maxe=-1, tote=0;
    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
    std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(trkmuon.key());
    art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle,evt,fHitTruthAssns);


    //spacepoint
    art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer); 
    //std::vector<recob::SpacePoint const*> spacepoint_vec;
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;

    for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
      particle_vec.clear(); match_vec.clear();
      spacepoint_vec.clear();

      particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
      //spacepoint_per_hit.get(hits_from_track[i_h].key(), spacepoint_vec);

      //const std::vector< art::Ptr<recob::SpacePoint> > spacepoint_vec = spacepoint_per_hit.at(hits_from_track[i_h].key());
      spacepoint_vec = spacepoint_per_hit.at(hits_from_track[i_h].key());

      //for(size_t i_s=0; i_s<spacepoint_vec.size(); ++i_s){
	//cout << "3D XYZ " << spacepoint_vec[i_s]->XYZ()[0] << " " << spacepoint_vec[i_s]->XYZ()[1] << " " << spacepoint_vec[i_s]->XYZ()[2] << endl;
      //}  

      for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
        trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
        //trkide[ particle_vec[i_p]->TrackId() ] ++;
        tote += match_vec[i_p]->energy;
        if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
          maxe = trkide[ particle_vec[i_p]->TrackId() ];
          matched_mcparticle = particle_vec[i_p];
        }
      }
    }

    if(matched_mcparticle){
      reco_ccmu_true_pdg = matched_mcparticle->PdgCode();
      const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
      reco_ccmu_true_origin=int(mc_truth->Origin());
      reco_ccmu_true_primary=matched_mcparticle->Process()=="primary";
      double x = matched_mcparticle->EndX();
      double y = matched_mcparticle->EndY();
      double z = matched_mcparticle->EndZ();
      reco_ccmu_true_end_inTPC = isInsideVolume("TPC",x,y,z);
      reco_ccmu_true_end_in5cmTPC = isInsideVolume("5cmTPC",x,y,z);
      reco_ccmu_true_end_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",x,y,z);
      TLorentzVector mcstart, mcend;
      unsigned int pstarti, pendi;
      reco_ccmu_true_length = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
    }
  }//isMC

  if (reco_ntracks>0) {

    // loop over tracks again to look for muon track
    //--cout << "Looking for muon tracks from " << NTracks << " tracks available" << endl;
    for(int j=0; j<NTracks; j++){

      art::Ptr<recob::Track> ptrack(trackListHandle,j);
      const recob::Track& track = *ptrack;

      std::map<art::Ptr<recob::PFParticle>,bool> PFPToNuSliceMap;
      std::map<art::Ptr<recob::PFParticle>, int> PFPToSliceIdMap;
      std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> trackToNuPFParticleMap;

	//Need Analyze track
	/*
      const art::Ptr<recob::PFParticle> pfps = trackToNuPFParticleMap[ptrack];
      m_reco_track_is_nuslice[j] = PFPToNuSliceMap[pfps];
      m_reco_track_sliceId[j] = PFPToSliceIdMap[pfps];
      cout << PFPToNuSliceMap[pfps] << " " << PFPToSliceIdMap[pfps] << endl;
	cout << "ggggggggg" << endl;
	*/

      // skip cc muon track
      if (ptrack.key()==trkmuon.key()) continue;

      TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
      TVector3 end(track.End().X(),track.End().Y(),track.End().Z());

      // check track is contained?
      if (isInsideVolume("TPC",pos) &&
          isInsideVolume("TPC",end)) {

        // check track start is far away from vertex
        float start_dis=TMath::Sqrt((reco_nu_vtx_x-pos.X())*(reco_nu_vtx_x-pos.X()) +
                                    (reco_nu_vtx_y-pos.Y())*(reco_nu_vtx_y-pos.Y()) +
                                    (reco_nu_vtx_z-pos.Z())*(reco_nu_vtx_z-pos.Z()));    
        if(start_dis>5){   // 5

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

            // check track length
            float trklen=track.Length();
            if(trklen>5){ // trklen>40 && trklen<60 // trklen>30 && trklen<70 // trklen>20 && trklen<70

              //save track ID for muon candidate
              muon_can_trkID.push_back(ptrack.key());
              cut_12=1;
              //--cout << "muon track candidate " << ptrack.key() << endl;

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

    //--cout << "Number of candidate muon tracks found: " << close_dis.size() << endl;

    //found muon track 5cm away from kaon track end
    if (mu_mindis_can.size() && k_mindis_can.size()) {

      cut_13=1;

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

      k_vtx_dis=TMath::Sqrt((pos_k.X()-reco_nu_vtx_x)*(pos_k.X()-reco_nu_vtx_x) +
                            (pos_k.Y()-reco_nu_vtx_y)*(pos_k.Y()-reco_nu_vtx_y) +
                            (pos_k.Z()-reco_nu_vtx_z)*(pos_k.Z()-reco_nu_vtx_z));

      //--cout << "Found muon track closest to kaon track end" << endl;
      //--cout << "kaon track key " << k_can_trkid << endl;
      //--cout << "muon track key " << mu_can_trkid << endl;
      //--cout << "distance between kaon and vertex " << k_vtx_dis << " cm" << endl;
      //--cout << "distance between muon and kaon " << k_mu_can_dis << " cm" << endl;

      ////////////////////////////////////////////////////////////////////// Doing End point matching for reco and true tracks /////////////////////////////////////
      ////////////////////////////////////////////////////////////////// End of end point matching //////////////////////////////
    } // minimum distance vectors are filled
    //} // kaon and muon ID vectors are filled

    //////////////////////////////////////////////////////// RECO-TRUTH MATCHING /////////////////////////////////////////////////////

/*
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
      kaon_vtx_dis=TMath::Sqrt((reco_nu_vtx_x-pos.X())*(reco_nu_vtx_x-pos.X()) +
                               (reco_nu_vtx_y-pos.Y())*(reco_nu_vtx_y-pos.Y()) +
                               (reco_nu_vtx_z-pos.Z())*(reco_nu_vtx_z-pos.Z()));

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

        //std::cout << "******************** Kaon track true matched ID : " << k_geant_ID << std::endl;
        //std::cout << "******************** Muon track true mother matched ID : " << matched_mcparticle->Mother() << std::endl;

        if(mu_origin==1 && mu_pdg==-13 && mu_isDec==1 && mu_k_is_Mother==1) mu_ness=1;

        if(mu_origin==1 && (mu_pdg==-13 || mu_pdg==211) && mu_isDec==1 && mu_k_is_Mother!=1){
          for(auto const& pPart : ptList){
            int mom_trkid=matched_mcparticle->Mother();
            if(pPart->TrackId()==mom_trkid){ 
              mu_mom_process=pPart->Process();
              //std::cout << "********************** Process name of the New Mother : " << mu_mom_process << std::endl;
              //std::cout << "********************** ID of the New Mother : " << pPart->Mother() << std::endl;
              if(pPart->Process()=="kaon+Inelastic" && pPart->Mother()==k_geant_ID) mu_mom_k_inelas=1;
              break;
            }
          }
        }

      }

      //cout << "Check muon track position" << endl;
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

      //cout << "Check muon track calorimetry" << endl;
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

      //cout << "Check muon track PID" << endl;
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
    //for(int k=0; k<NTracks; k++){
    {

      //art::Ptr<recob::Track> ptrack(trackListHandle,trkccmu);
      //const recob::Track& track = *ptrack;

      //cout << "track " << k << " track id " << track.ID() << " track key " << ptrack.key() << endl;

      simb::MCParticle const* matched_mcparticle = NULL;
      std::unordered_map<int,double> trkide;

      double maxe=-1, tote=0;
      std::vector<simb::MCParticle const*> particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
      //std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(trkmuon.key());
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
        cc_mu_trkid=trkmuon.key();
        if(matched_mcparticle->Process()==pri && TMath::Abs(cc_mu_pdg)==13 && cc_mu_origin==1){ 
          pri_Mu_is=1;
          //break;
        }
      }

    }

    ////////////////////////////////////////////// End of truth Muon inforamtion  //////////////////////////////////////////////////////////////////////////////

    if(cc_mu_trkid!=-9999){

      //cout << "Check CC muon track position" << endl;
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

      //cout << "Check CC muon track calorimetry" << endl;
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

  //cout << "Looking into the truth level information" << endl;

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
        //std::cout << "********************** Truth level Kaon trk ID : " << gk_trk_ID <<  std::endl;
        //std::cout << "****************************** End Energy of Primary Kaon : " << pPart->EndE()*1000 << std::endl;
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
      //std::cout << "********** Truth info Found a decaying " << sPart->PdgCode() << std::endl;
      //std::cout << "********** Mother trkID of the Found decay particle " << sPart->Mother() << std::endl;
    }
    if(sPart->Process()=="Decay" && (sPart->PdgCode()==-13 || sPart->PdgCode()==211) && sPart->Mother()==gk_trk_ID){
      //std::cout << "****************************** We found typical kaon decaying event ********************************" << std::endl;
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
          //std::cout << "************** Process Name : " << pPart->Process() << std::endl;
          //std::cout << "************** Mother ID of Mother : " << pPart->Mother() << std::endl;	  
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

            //std::cout << "****************************** We found non-typical kaon decaying event ********************************" << std::endl;
            //std::cout << "****************************** End Energy of Secondary Kaon : " << pPart->EndE()*1000 << std::endl;

            True_kinelas_KE=(pPart->E()-0.493);
            TLorentzVector mcstart, mcend;
            unsigned int pstarti, pendi;
            True_kinelas_tlen=length(*pPart, mcstart, mcend, pstarti, pendi);
            Tr_K_Inelas=1;    

            //std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Testing Kaon Inelastic interaction &&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
            //std::cout << "Inelastic Kaon has his own reco track : " << kinelas_has_traks << std::endl;
            //std::cout << "Reco track ID of Inelastic Kaon and Primary Kaon : " << kinelas_reco_trkID << "  " << k_can_trkid << std::endl;
            //std::cout << "True KE of Inelastic Kaon : " << True_kinelas_KE << std::endl;
            //std::cout << "True tracklength of Inelastic Kaon : " << True_kinelas_tlen << std::endl;
            //std::cout << "Kaon Inelastic swith is " << Tr_K_Inelas << std::endl;

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
    //cout << "Getting Calorimetric information of Kaon and Muon" << endl;
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
*/
  }

  //endCanvas();
  cout << "Filling tree ----------------" << endl;


  //cout << "before filling reco_track_daughter_true_pdg_sh[track_i][daughter_i]: " << reco_track_daughter_true_pdg_sh[0][0] << endl;

  fEventTree->Fill();

  //cout << "after filling reco_track_daughter_true_pdg_sh[track_i][daughter_i]: " << reco_track_daughter_true_pdg_sh[0][0] << endl;


} // end of analyze function
 
 /////////////////////////////////////////// Reset Function ///////////////////////////////
 
 void CCKaonAnalyzerRebuild::reset(){

   run=-9999;
   subrun=-9999;
   event=-9999;
   
   event_weight.clear();
   evtwgt_funcname.clear();
   evtwgt_weight.clear();
   evtwgt_nweight.clear();
   
   
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
   
   IsKaon=false;
   IsSingleKaon=false;
   IsAssociatedKaon=false;
   IsMuBR=false;
   IsPiBR=false;
   flag_prim_k=false;
   
   true_lepton_pdg=-9999;
   true_lepton_p=-9999;
   true_lepton_ke=-9999;
   true_lepton_theta=-9999;
   true_lepton_costheta=-9999;
   true_lepton_phi=-9999;
   
   true_nkaons=-9999;
   true_kaon_length=-9999;
   true_kaon_p=-9999;
   true_kaon_ke=-9999;
   true_kaon_theta=-9999;
   true_kaon_costheta=-9999;
   true_kaon_phi=-9999;
   true_kaon_ccmuon_angle=-9999;
   true_kaon_ccmuon_cosangle=-9999;
   true_kaon_end_process=-9999;
   true_kaon_end_ke=-9999;
   true_kaon_start_x=-9999;
   true_kaon_start_y=-9999;
   true_kaon_start_z=-9999;
   true_kaon_end_x=-9999;
   true_kaon_end_y=-9999;
   true_kaon_end_z=-9999;
   true_kaon_end_inTPC=false;
   true_kaon_end_in5cmTPC=false;
   true_kaon_end_inCCInclusiveTPC=false;
   
   true_dau_muon_length=-9999;
   true_dau_muon_p=-9999;
   true_dau_muon_ke=-9999;
   true_dau_muon_theta=-9999;
   true_dau_muon_costheta=-9999;
   true_dau_muon_phi=-9999;
   true_dau_muon_ccmuon_angle=-9999;
   true_dau_muon_ccmuon_cosangle=-9999;
   true_dau_muon_end_process=-9999;
   true_dau_muon_end_ke=-9999;
   true_dau_muon_start_x=-9999;
   true_dau_muon_start_y=-9999;
   true_dau_muon_start_z=-9999;
   true_dau_muon_end_x=-9999;
   true_dau_muon_end_y=-9999;
   true_dau_muon_end_z=-9999;
   
   true_dau_pip_length=-9999;
   true_dau_pip_p=-9999;
   true_dau_pip_ke=-9999;
   true_dau_pip_theta=-9999;
   true_dau_pip_costheta=-9999;
   true_dau_pip_phi=-9999;
   true_dau_pip_ccmuon_angle=-9999;
   true_dau_pip_ccmuon_cosangle=-9999;
   true_dau_pip_end_process=-9999;
   true_dau_pip_end_ke=-9999;
   true_dau_pip_start_x=-9999;
   true_dau_pip_start_y=-9999;
   true_dau_pip_start_z=-9999;
   true_dau_pip_end_x=-9999;
   true_dau_pip_end_y=-9999;
   true_dau_pip_end_z=-9999;
   
   
   true_dau_pin_length=-9999;
   true_dau_pin_p=-9999;
   true_dau_pin_ke=-9999;
   true_dau_pin_theta=-9999;
   true_dau_pin_costheta=-9999;
   true_dau_pin_phi=-9999;
   true_dau_pin_ccmuon_angle=-9999;
   true_dau_pin_ccmuon_cosangle=-9999;
   true_dau_pin_end_process=-9999;
   true_dau_pin_end_ke=-9999;
   true_dau_pin_start_x=-9999;
   true_dau_pin_start_y=-9999;
   true_dau_pin_start_z=-9999;
   true_dau_pin_end_x=-9999;
   true_dau_pin_end_y=-9999;
   true_dau_pin_end_z=-9999;


   for(int k=0; k<kMaxParticles; k++){ 
     true_length[k]=-9999;
     true_p[k]=-9999;
     true_ke[k]=-9999;
     true_theta[k]=-9999;
     true_costheta[k]=-9999;
     true_phi[k]=-9999;
     true_ccmuon_angle[k]=-9999;
     true_ccmuon_cosangle[k]=-9999;
     true_end_process[k]=-9999;
     true_end_ke[k]=-9999;
     true_end_x[k]=-9999;
     true_end_y[k]=-9999;
     true_end_z[k]=-9999;
     true_end_inTPC[k]=false;
     true_end_in5cmTPC[k]=false;
     true_end_inCCInclusiveTPC[k]=false;
   }

   true_kaon_ndaughters=-9999;
   true_kaon_ndaughters_decay=-9999;
   true_kaon_ndaughters_inelastic=-9999;
   true_kaon_ndecmup=-9999;
   true_kaon_ndecpip=-9999;
   true_kaon_ninekap=-9999;
   true_kaon_ninepip=-9999;
   true_kaon_ninepro=-9999;
   true_kaon_daughter_length=-9999;
   true_kaon_daughter_p=-9999;
   true_kaon_daughter_ke=-9999;
   true_kaon_daughter_theta=-9999;
   true_kaon_daughter_costheta=-9999;
   true_kaon_daughter_angle=-9999;
   true_kaon_daughter_cosangle=-9999;
   true_kaon_daughter_pdg=-9999;
   true_kaon_daughter_start_x=-9999;
   true_kaon_daughter_start_y=-9999;
   true_kaon_daughter_start_z=-9999;
   true_kaon_daughter_end_x=-9999;
   true_kaon_daughter_end_y=-9999;
   true_kaon_daughter_end_z=-9999;
   true_kaon_daughter_end_inTPC=false;
   true_kaon_daughter_end_in5cmTPC=false;
   true_kaon_daughter_end_inCCInclusiveTPC=false;
   
   true_nhyperons=-9999;
   
   true_nkaons_anti=-9999;
   true_kaon_length_anti=-9999;
   true_kaon_p_anti=-9999;
   true_kaon_ke_anti=-9999;
   true_kaon_theta_anti=-9999;
   true_kaon_costheta_anti=-9999;
   true_kaon_phi_anti=-9999;
   true_kaon_ccmuon_angle_anti=-9999;
   true_kaon_ccmuon_cosangle_anti=-9999;
   true_kaon_end_process_anti=-9999;
   true_kaon_end_ke_anti=-9999;
   true_kaon_end_x_anti=-9999;
   true_kaon_end_y_anti=-9999;
   true_kaon_end_z_anti=-9999;
   true_kaon_end_inTPC_anti=false;
   true_kaon_end_in5cmTPC_anti=false;
   true_kaon_end_inCCInclusiveTPC_anti=false;
   
   true_kaon_ndaughters_anti=-9999;
   true_kaon_ndaughters_decay_anti=-9999;
   true_kaon_ndaughters_inelastic_anti=-9999;
   true_kaon_ndecmup_anti=-9999;
   true_kaon_ndecpip_anti=-9999;
   true_kaon_ninekap_anti=-9999;
   true_kaon_ninepip_anti=-9999;
   true_kaon_ninepro_anti=-9999;
   true_kaon_daughter_length_anti=-9999;
   true_kaon_daughter_p_anti=-9999;
   true_kaon_daughter_ke_anti=-9999;
   true_kaon_daughter_theta_anti=-9999;
   true_kaon_daughter_costheta_anti=-9999;
   true_kaon_daughter_angle_anti=-9999;
   true_kaon_daughter_cosangle_anti=-9999;
   true_kaon_daughter_pdg_anti=-9999;
   true_kaon_daughter_end_x_anti=-9999;
   true_kaon_daughter_end_y_anti=-9999;
   true_kaon_daughter_end_z_anti=-9999;
   true_kaon_daughter_end_inTPC_anti=false;
   true_kaon_daughter_end_in5cmTPC_anti=false;
   true_kaon_daughter_end_inCCInclusiveTPC_anti=false;
   
   cheat_pip_trkln = -9999;
   cheat_mup_trkln = -9999;

   reco_nu_cc_filter=false;
   
   reco_nu_vtx_x=-9999;
   reco_nu_vtx_y=-9999;
   reco_nu_vtx_z=-9999;
   reco_nu_vtx_inTPC=false;
   reco_nu_vtx_in5cmTPC=false;
   reco_nu_vtx_inCCInclusiveTPC=false;
   reco_nu_ndaughters=-9999;
   reco_nu_cc_nmue=-9999;
   
   reco_ccmu_vtx_x=-9999;
   reco_ccmu_vtx_y=-9999;
   reco_ccmu_vtx_z=-9999;
   reco_ccmu_vtx_inTPC=false;
   reco_ccmu_vtx_in5cmTPC=false;
   reco_ccmu_vtx_inCCInclusiveTPC=false;
   reco_ccmu_true_pdg=-9999;
   reco_ccmu_true_origin=-9999;
   reco_ccmu_true_primary=false;
   reco_ccmu_true_end_inTPC=false;
   reco_ccmu_true_end_in5cmTPC=false;
   reco_ccmu_true_end_inCCInclusiveTPC=false;
   reco_ccmu_true_length=-9999;
   
   //reco_track_dEdx.clear();
   //reco_track_ResRan.clear();
   
   //reco_track_daughter_dEdx.clear();
   //reco_track_daughter_ResRan.clear();
   //vector<Float_t> v;
   
   reco_ntracks=0;
   reco_nshowers=0;
   
   for(int k=0; k<kMaxTracks; k++){


     for(int l=0; l<kMaxParticles; l++){ 
       cheat_num_hits[k][l] = -9999;
       cheat_ini_hits[k][l] = -9999;
       cheat_closest_distance[k][l] = -9999;       

       cheat_peak_pdg[k][l] = -9999;
       cheat_peak_theta[k][l] = -9999;
       cheat_peak_phi[k][l] = -9999;
       best_peak_theta[k][l] = -9999;
       best_peak_phi[k][l] = -9999;
       best_peak_trkln[k][l] = -9999;
     }  

     //m_reco_track_sliceId[k]=-9;
     //m_reco_track_is_nuslice[k]=-9;
     reco_track_start_x[k]=-9;
     reco_track_start_y[k]=-9;
     reco_track_start_z[k]=-9;
     reco_track_end_x[k]=-9;
     reco_track_end_y[k]=-9;
     reco_track_end_z[k]=-9;
     
     reco_track_distance[k]=-9;
     reco_track_nhits0[k]=-9;
     reco_track_nhits1[k]=-9;
     reco_track_nhits2[k]=-9;
     reco_track_kin0[k]=-9;
     reco_track_kin1[k]=-9;
     reco_track_kin2[k]=-9;
     
     //reco_track_dEdx_arr.at(k).clear();
     //reco_track_ResRan_arr.at(k).clear();
     
     reco_track_length[k]=-9;
     reco_track_theta[k]=-9;
     reco_track_phi[k]=-9;
     
     reco_track_P_vtx[k]=-9;
     reco_track_P_str[k]=-9;
     reco_track_P_end[k]=-9;
     
     reco_track_dir[k]=false;
     reco_track_chi2ka_pl0[k]=-9;
     reco_track_chi2pr_pl0[k]=-9;
     reco_track_chi2pi_pl0[k]=-9;
     reco_track_chi2mu_pl0[k]=-9;
     reco_track_chi2ka_pl1[k]=-9;
     reco_track_chi2pr_pl1[k]=-9;
     reco_track_chi2pi_pl1[k]=-9;
     reco_track_chi2mu_pl1[k]=-9;
     reco_track_chi2ka_pl2[k]=-9;
     reco_track_chi2pr_pl2[k]=-9;
     reco_track_chi2pi_pl2[k]=-9;
     reco_track_chi2mu_pl2[k]=-9;
     
     
     reco_track_Bragg_fwd_ka_pl0[k] = -9; 
     reco_track_Bragg_fwd_pr_pl0[k] = -9; 
     reco_track_Bragg_fwd_pi_pl0[k] = -9; 
     reco_track_Bragg_fwd_mu_pl0[k] = -9; 
     reco_track_Bragg_fwd_ka_pl1[k] = -9;  
     reco_track_Bragg_fwd_pr_pl1[k] = -9;  
     reco_track_Bragg_fwd_pi_pl1[k] = -9; 
     reco_track_Bragg_fwd_mu_pl1[k] = -9; 
     reco_track_Bragg_fwd_ka_pl2[k] = -9; 
     reco_track_Bragg_fwd_pr_pl2[k] = -9; 
     reco_track_Bragg_fwd_pi_pl2[k] = -9; 
     reco_track_Bragg_fwd_mu_pl2[k] = -9; 
     
     reco_track_MIP_pl0[k] = -9; 
     reco_track_MIP_pl1[k] = -9; 
     reco_track_MIP_pl2[k] = -9; 
     
     reco_track_chi2ka_3pl[k]=-99;
     reco_track_chi2pr_3pl[k]=-99;
     reco_track_chi2pi_3pl[k]=-99;
     reco_track_chi2mu_3pl[k]=-99;
     reco_track_likepr_3pl[k]=-99;
     reco_track_llrpid_3pl[k] = -9;
     reco_track_total_llrpid_3pl[k] = -9;
     reco_track_llrpid_k_3pl[k] = -9;
     reco_track_vtx_inTPC[k]=false;
     reco_track_vtx_in5cmTPC[k]=false;
     reco_track_vtx_inCCInclusiveTPC[k]=false;
     reco_track_end_inTPC[k]=false;
     reco_track_end_in5cmTPC[k]=false;
     reco_track_end_inCCInclusiveTPC[k]=false;
     reco_track_true_pdg[k]=-999;
     reco_track_true_origin[k]=-99;
     reco_track_true_primary[k]=false;
     reco_track_true_end_inTPC[k]=false;
     reco_track_true_end_in5cmTPC[k]=false;
     reco_track_true_end_inCCInclusiveTPC[k]=false;
     reco_track_true_length[k]=-9;
     reco_track_ndaughters[k]=0;
     reco_shower_ndaughters[k]=0;
     
     for(int l=0; l<kMaxMerge; l++){ 
       reco_track_match_e[k][l] = -999;
       reco_track_match_hit[k][l] = -999;
       reco_track_match_epdg[k][l] = -999;
       reco_track_match_hitpdg[k][l] = -999;
     }
     
     
       for(int d=0; d<2000; d++){
       reco_track_dEdx_pl0[k][d]=-9;
       reco_track_ResRan_pl0[k][d]=-9;
       reco_track_dEdx_pl1[k][d]=-9;
       reco_track_ResRan_pl1[k][d]=-9;
       reco_track_dEdx_pl2[k][d]=-9;
       reco_track_ResRan_pl2[k][d]=-9;
       }
     
      
     for(int m=0; m<kMaxShowers; m++){ 
       
       reco_track_daughter_distance_sh[k][m]=-9;
       reco_track_daughter_vtx_distance_sh[k][m]=-9;
       reco_angle_track_daughter_sh[k][m]=-9;
       reco_angle_daughter_track_daughter_sh[k][m]=-9;
       
       reco_track_daughter_length_sh[k][m]=-99;
       reco_track_daughter_theta_sh[k][m]=-99;
       reco_track_daughter_phi_sh[k][m]=-99;
       reco_track_daughter_open_angle_sh[k][m]=-99;
       reco_track_daughter_dedx_pl0_sh[k][m]=-99;
       reco_track_daughter_dedx_pl1_sh[k][m]=-99;
       reco_track_daughter_dedx_pl2_sh[k][m]=-99;
       reco_track_daughter_energy_pl0_sh[k][m]=-99;
       reco_track_daughter_energy_pl1_sh[k][m]=-99;
       reco_track_daughter_energy_pl2_sh[k][m]=-99;
       
       reco_track_true_pdg_sh[k]=-999;
       reco_track_true_origin_sh[k]=-99;
       reco_track_true_primary_sh[k]=false;
       reco_track_true_end_inTPC_sh[k]=false;
       reco_track_true_end_in5cmTPC_sh[k]=false;
       reco_track_true_end_inCCInclusiveTPC_sh[k]=false;
       reco_track_true_length_sh[k]=-9;
       
       reco_track_daughter_true_pdg_sh[k][m]=-999;
       reco_track_daughter_true_origin_sh[k][m]=-99;
       reco_track_daughter_true_primary_sh[k][m]=false;
       reco_track_daughter_true_end_inTPC_sh[k][m]=false;
       reco_track_daughter_true_end_in5cmTPC_sh[k][m]=false;
       reco_track_daughter_true_end_inCCInclusiveTPC_sh[k][m]=false;
       reco_track_daughter_true_length_sh[k][m]=-9999;
       reco_track_daughter_true_mother_sh[k][m]=-9999;
       
       reco_track_daughter_vtx_inTPC_sh[k][m]=false;
       reco_track_daughter_vtx_in5cmTPC_sh[k][m]=false;
       reco_track_daughter_vtx_inCCInclusiveTPC_sh[k][m]=false;
       reco_track_daughter_end_inTPC_sh[k][m]=false;
       reco_track_daughter_end_in5cmTPC_sh[k][m]=false;
       reco_track_daughter_end_inCCInclusiveTPC_sh[k][m]=false;
       
       for(int l=0; l<kMaxMerge; l++){ 
	 reco_track_daughter_shower_match_e[k][m][l] = -999;
	 reco_track_daughter_shower_match_hit[k][m][l] = -999;
	 reco_track_daughter_shower_match_epdg[k][m][l] = -999;
	 reco_track_daughter_shower_match_hitpdg[k][m][l] = -999;
       }
     }
     
     for(int m=0; m<kMaxTracks; m++){
       reco_track_daughter_start_x[k][m]=-9;
       reco_track_daughter_start_y[k][m]=-9;
       reco_track_daughter_start_z[k][m]=-9;
       reco_track_daughter_end_x[k][m]=-9;
       reco_track_daughter_end_y[k][m]=-9;
       reco_track_daughter_end_z[k][m]=-9;
       

       reco_track_daughter_distance[k][m]=-9;
       reco_track_daughter_vtx_distance[k][m]=-9;
       reco_angle_track_daughter[k][m]=-9;
       reco_track_daughter_nhits0[k][m]=-9;
       reco_track_daughter_nhits1[k][m]=-9;
       reco_track_daughter_nhits2[k][m]=-9;
       
          for(int d=0; d<2000; d++){
          reco_track_daughter_dEdx_pl0[k][m][d]=-9;
          reco_track_daughter_ResRan_pl0[k][m][d]=-9;
          reco_track_daughter_dEdx_pl1[k][m][d]=-9;
          reco_track_daughter_ResRan_pl1[k][m][d]=-9;
          reco_track_daughter_dEdx_pl2[k][m][d]=-9;
          reco_track_daughter_ResRan_pl2[k][m][d]=-9;
	  }
       
       
       for(int l=0; l<kMaxMerge; l++){ 
	 reco_track_daughter_match_e[k][m][l] = -999;
	 reco_track_daughter_match_hit[k][m][l] = -999;
	 reco_track_daughter_match_epdg[k][m][l] = -999;
	 reco_track_daughter_match_hitpdg[k][m][l] = -999;
       }

       reco_track_daughter_length[k][m]=-99;
       reco_track_daughter_theta[k][m]=-99;
       reco_track_daughter_phi[k][m]=-99;
       reco_track_daughter_chi2ka_pl0[k][m]=-9;
       reco_track_daughter_chi2pr_pl0[k][m]=-9;
       reco_track_daughter_chi2pi_pl0[k][m]=-9;
       reco_track_daughter_chi2mu_pl0[k][m]=-9;
       reco_track_daughter_chi2ka_pl1[k][m]=-9;
       reco_track_daughter_chi2pr_pl1[k][m]=-9;
       reco_track_daughter_chi2pi_pl1[k][m]=-9;
       reco_track_daughter_chi2mu_pl1[k][m]=-9;
       reco_track_daughter_chi2ka_pl2[k][m]=-9;
       reco_track_daughter_chi2pr_pl2[k][m]=-9;
       reco_track_daughter_chi2pi_pl2[k][m]=-9;
       reco_track_daughter_chi2mu_pl2[k][m]=-9;
       reco_track_daughter_chi2ka_3pl[k][m]=-99;
       reco_track_daughter_chi2pr_3pl[k][m]=-99;
       reco_track_daughter_chi2pi_3pl[k][m]=-99;
       reco_track_daughter_chi2mu_3pl[k][m]=-99;
       reco_track_daughter_likepr_3pl[k][m]=-99;
       reco_track_daughter_llrpid_3pl[k][m]=-9;
       reco_track_daughter_llrpid_k_3pl[k][m]=-9;
       reco_track_daughter_vtx_inTPC[k][m]=false;
       reco_track_daughter_vtx_in5cmTPC[k][m]=false;
       reco_track_daughter_vtx_inCCInclusiveTPC[k][m]=false;
       reco_track_daughter_end_inTPC[k][m]=false;
       reco_track_daughter_end_in5cmTPC[k][m]=false;
       reco_track_daughter_end_inCCInclusiveTPC[k][m]=false;
       reco_track_daughter_true_pdg[k][m]=-999;
       reco_track_daughter_true_origin[k][m]=-99;
       reco_track_daughter_true_primary[k][m]=false;
       reco_track_daughter_true_end_inTPC[k][m]=false;
       reco_track_daughter_true_end_in5cmTPC[k][m]=false;
       reco_track_daughter_true_end_inCCInclusiveTPC[k][m]=false;
       reco_track_daughter_true_length[k][m]=-9999;
       reco_track_daughter_true_mother[k][m]=-9999;
       
       reco_track_daughter_Bragg_fwd_ka_pl0[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pr_pl0[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pi_pl0[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_mu_pl0[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_ka_pl1[k][m] = -999;  
       reco_track_daughter_Bragg_fwd_pr_pl1[k][m] = -999;  
       reco_track_daughter_Bragg_fwd_pi_pl1[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_mu_pl1[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_ka_pl2[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pr_pl2[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pi_pl2[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_mu_pl2[k][m] = -999; 
       
       reco_track_daughter_MIP_pl0[k][m] = -999; 
       reco_track_daughter_MIP_pl1[k][m] = -999; 
       reco_track_daughter_MIP_pl2[k][m] = -999; 
       
       
     }
   }
   
   k_can_trkid=-999;
   mu_can_trkid=-999;
   k_mu_can_dis=-999;
   k_mu_open_angle=-999;
   k_vtx_dis=-999;
   //k_geant_ID=-9999;
   //k_origin=-9999;
   //k_pdg=-9999;
   //k_isPri=-9999;
   //k_endE=-9999;
   //k_ness=-9999;
   //kaon_vtx_dis=-9999;
   //k_plen=-9999;
   //k_phi=-9999;
   //k_theta=-9999;
   //k_in_5_TPC=-9999;
   //k_in_CC_TPC=-9999;
   //k_hit=-9999;
   //k_range=-9999;
   //k_KE=-9999;
   //k_large_dedx=-9999;
   //k_small_dedx=-9999;
   //k_chi_p=-9999;
   //k_chi_k=-9999;
   //k_chi_pi=-9999;
   //k_chi_mu=-9999;
   //k_p_max=-9999;
   //k_k_max=-9999;
   //k_pi_max=-9999;
   //k_mu_max=-9999;
   //k_mip_max=-9999;
   //k_L1_ratio=-9999;
   //k_LL1=-9999;
   //k_L2_ratio=-9999;
   //k_LL2=-9999;
   //k_Lp_ratio=-9999;
   //k_LLp=-9999;
   //k_Lk_ratio=-9999;
   //k_LLk=-9999;
   //k_pida_mean=-9999;
   //k_pida_med=-9999;
   //k_kde=-9999;
   //k_trm_dqdx=-9999;
   //k_trm_dedx=-9999;
   //mu_pdg=-9999;
   //mu_isDec=-9999;
   //mu_origin=-9999;
   //mu_k_is_Mother=-9999;
   //mu_mom_k_inelas=-9999;
   //mu_ness=-9999;
   //mu_plen=-9999;
   //mu_phi=-9999;
   //mu_theta=-9999;
   //mu_in_5_TPC=-9999;
   //mu_in_CC_TPC=-9999;
   //mu_KE=-9999;
   //mu_hit=-9999;
   //mu_range=-9999;
   //mu_large_dedx=-9999;
   //mu_small_dedx=-9999;
   //mu_chi_p=-9999;
   //mu_chi_k=-9999;
   //mu_chi_pi=-9999;
   //mu_chi_mu=-9999;
   //mu_p_max=-9999;
   //mu_k_max=-9999;
   //mu_pi_max=-9999;
   //mu_mu_max=-9999;
   //mu_mip_max=-9999;
   //mu_L1_ratio=-9999;
   //mu_LL1=-9999;
   //mu_L2_ratio=-9999;
   //mu_LL2=-9999;
   //mu_Lp_ratio=-9999;
   //mu_LLp=-9999;
   //mu_Lk_ratio=-9999;
   //mu_LLk=-9999;
   //mu_pida_mean=-9999;
   //mu_pida_med=-9999;
   //mu_kde=-9999;
   //mu_trm_dqdx=-9999;
   //mu_trm_dedx=-9999;
   //mu_mom_process="NA";
   //cc_mu_trkid=-9999;
   //cc_mu_tlen=-9999;
   //cc_mu_phi=-9999;
   //cc_mu_pdg=-9999;
   //cc_mu_theta=-9999;
   //cc_mu_range=-9999;
   //cc_mu_KE=-9999;
   //cc_mu_hit=-9999;
   //cc_mu_large_dedx=-9999;
   //cc_mu_small_dedx=-9999;
   //cc_dis_vtx=-9999;
   //cc_mu_chi_p=-9999;
   //cc_mu_chi_k=-9999;
   //cc_mu_chi_pi=-9999;
   //cc_mu_chi_mu=-9999;
   //cc_mu_p_max=-9999;
   //cc_mu_k_max=-9999;
   //cc_mu_pi_max=-9999;
   //cc_mu_mu_max=-9999;
   //cc_mu_mip_max=-9999;
   //cc_mu_L1_ratio=-9999;
   //cc_mu_LL1=-9999;
   //cc_mu_L2_ratio=-9999;
   //cc_mu_LL2=-9999;
   //cc_mu_Lp_ratio=-9999;
   //cc_mu_LLp=-9999;
   //cc_mu_Lk_ratio=-9999;
   //cc_mu_LLk=-9999;
   //cc_mu_pida_mean=-9999;
   //cc_mu_pida_med=-9999;
   //cc_mu_kde=-9999;
   //cc_mu_trm_dqdx=-9999;
      //cc_mu_trm_dedx=-9999;
      //longest_trkid=-9999;
      //longest_trklen=-9999;
      //pri_Mu_is=-9999;
      //Tr_pri_mu_pdg=-9999;
      //Tr_pri_mu_is=-9999;
      //Tr_pri_st_k_is=-9999;
      //Tr_K_Inelas=-9999;
      //Tr_k_plen=-9999;
      //Tr_k_endE=-9999;
      //Tr_k_theta=-9999;
      //Tr_k_phi=-9999;
      //Tr_dec_mu_is=-9999;
      //Tr_dec_mu_pi_pdg=-9999;
      //Tr_mu_plen=-9999;
      //Tr_mu_theta=-9999;
      //Tr_mu_phi=-9999;
      //Tr_k_inTPC=-9999;
      //Tr_mu_inTPC=-9999;
      //Tr_k_in_5_TPC=-9999;
      //Tr_k_in_CC_TPC=-9999;
      //Tr_mu_in_5_TPC=-9999;
      //Tr_mu_in_CC_TPC=-9999;
      //Tr_kmu_open_ang=-9999;
      //vtx_5cm_mult=-9999;
      //k_start_dedx=-9999;
      //k_end_dedx=-9999;
      //mu_start_dedx=-9999;
      //mu_end_dedx=-9999;
      cut_1=-9;
      cut_2=-9;
      cut_3=-9;
      cut_4=-9;
      cut_5=-9;
      cut_6=-9;
      cut_7=-9;
      cut_8=-9;
      cut_9=-9;
      cut_10=-9;
      cut_11=-9;
      cut_12=-9;
      cut_13=-9;
      //kinelas_has_traks=-9999;
      //kinelas_reco_trkID=-9999;
      //kinelas_tlen=-9999;
      //True_kinelas_KE=-9999;
      //True_kinelas_tlen=-9999;
      
      //for(int j=0; j<3; j++){
      //  for(int k=0; k<3000; k++){
      //    k_dedx[j][k]=-9999;
      //    k_rr[j][k]=-9999;
      //    mu_dedx[j][k]=-9999;
      //    mu_rr[j][k]=-9999;
      //  }
      //}
      //
      //delete [] reco_track_daughter_dEdx_a;
      
      rv0.clear();
      dv0.clear();

      rv1.clear();
      dv1.clear();

      rv2.clear();
      dv2.clear();

 }
 
 /////////////////////////////////////////////////////////////////////////////////////////
 
bool CCKaonAnalyzerRebuild::isInsideVolume(string volume, double x, double y, double z)
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
  return false;
}

double CCKaonAnalyzerRebuild::length(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
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

//void CCKaonAnalyzerRebuild::fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, const recob::Track trk, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart, int track_i, int daughter_i)
//{
void CCKaonAnalyzerRebuild::fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, int track_i, int daughter_i)
{

  int hits_p0=0;
  int hits_p1=0;
  int hits_p2=0;
  float kin_p0=0;
  float kin_p1=0;
  float kin_p2=0;
    //vector<vector<Float_t>> dvv;
  //vector<vector<Float_t>> rvv;
  //double fADCtoE[]={232,249,243.7};



  double llr_pid_total = 0;
  double llr_pid_total_k = 0;

  //for(unsigned int ical=0; ical<calos.size(); ++ical){
  for (auto const &calo : calos) {

    //if(!calo) continue;
    //if(!calo->PlaneID().isValid) continue;
    int planenum = calo->PlaneID().Plane;
    //if(planenum<0||planenum>2) continue; 
    if(planenum==0) {
            hits_p0=calo->dEdx().size();
            kin_p0=calo->KineticEnergy();
    }

    if(planenum==1) {
            hits_p1=calo->dEdx().size();
            kin_p1=calo->KineticEnergy();
    }

    if(planenum==2) {
            hits_p2=calo->dEdx().size();
            kin_p2=calo->KineticEnergy();
    }

    //auto calo = calos[ical];
    auto const &plane = calo->PlaneID().Plane;
    auto const &dedx_values = calo->dEdx();
    //auto const &dqdx_values = calo->dQdx();

    if(planenum==0){
      for(auto& dedx : calo->dEdx()){
	cout << "calo->dEdx(): " << dedx << endl;
      }
    }


    if(planenum==0) dv0 = calo->dEdx();
    if(planenum==0) rv0 = calo->ResidualRange();

    if(planenum==1) dv1 = calo->dEdx();
    if(planenum==1) rv1 = calo->ResidualRange();

    if(planenum==2) dv2 = calo->dEdx();
    if(planenum==2) rv2 = calo->ResidualRange();

    auto const &rr = calo->ResidualRange();
    auto const &pitch = calo->TrkPitchVec();
    //auto const &xyz_v = calo->XYZ();

    std::vector<std::vector<float>> par_values;
    par_values.push_back(rr);
    par_values.push_back(pitch);
/*
    std::vector<float> dqdx_values_corrected;//, dedx_values_corrected;

    if (!isMC){
        dqdx_values_corrected = dqdx_values;
    }
    else{
        dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(calo, trk, assocMCPart, true, 0.1, false);
    }
    //dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(calo, value, assocMCPart, true, 0.1, false);

    for (size_t i = 0; i < dqdx_values_corrected.size(); i++){
    float aux_dedx;
    aux_dedx = ModBoxCorrection(dqdx_values_corrected[i]*fADCtoE[plane], xyz_v[i].X(), xyz_v[i].Y(), xyz_v[i].Z());
    dedx_values_corrected.push_back(aux_dedx);
    calo_energy += aux_dedx * pitch[i];
    }

    if(planenum==2){ 
            for(auto& n:dqdx_values){cout << "dQdx before: " << n <<endl;}
            for(auto& n:dqdx_values_corrected){cout << "dQdx after: " << n <<endl;}
    }
*/
    if (calo->ResidualRange().size() == 0) continue;

    //llr_pid_total += llr_pid_calculator.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);
    //llr_pid_total_k += llr_pid_calculator_k.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);
    llr_pid_total += llr_pid_calculator.LLR_many_hits_one_plane(dedx_values, par_values, plane);
    llr_pid_total_k += llr_pid_calculator_k.LLR_many_hits_one_plane(dedx_values, par_values, plane);

  }

  double llr_pid_score = atan(llr_pid_total / 100.) * 2 / 3.14159266;
  double llr_pid_score_k = atan(llr_pid_total_k / 100.) * 2 / 3.14159266;



  if (daughter_i<0) {
    reco_track_nhits0[track_i] = hits_p0;
    reco_track_nhits1[track_i] = hits_p1;
    reco_track_nhits2[track_i] = hits_p2;
    reco_track_kin0[track_i] = kin_p0;
    reco_track_kin1[track_i] = kin_p1;
    reco_track_kin2[track_i] = kin_p2;
    cout << "Kin p0: " << kin_p0 << " Kin p1: " << kin_p1 << " Kin p2: " << kin_p2 <<endl;
    reco_track_llrpid_3pl[track_i] = llr_pid_score;
    reco_track_total_llrpid_3pl[track_i] = llr_pid_total;
    reco_track_llrpid_k_3pl[track_i] = llr_pid_score_k;

    //reco_track_dEdx.push_back(dv);
    //unsigned int ti = track_i;
    ////reco_track_ResRan.resize(kMaxTracks, std::vector<Float_t>(hits_p2));
    ////reco_track_ResRan[track_i] = rv;
    //reco_track_ResRan.push_back(rv);
    
    //if(dv.size()!=0){
    ////reco_track_dEdx.resize(kMaxTracks, std::vector<Float_t>(hits_p2));
    ////reco_track_dEdx[track_i] = dv;
    //reco_track_ResRan[track_i] = rv;
    //} else{    
    //reco_track_dEdx_arr = dvv;
    //cout << "dEdx has nothing: " << reco_track_dEdx_arr.at(track_i).at(0);
    //}
    
    for(unsigned int d=0; d<dv0.size(); d++){       
       reco_track_dEdx_pl0[track_i][d] = dv0[d];
       reco_track_ResRan_pl0[track_i][d] = rv0[d];
    }
    for(unsigned int d=0; d<dv1.size(); d++){       
       reco_track_dEdx_pl1[track_i][d] = dv1[d];
       reco_track_ResRan_pl1[track_i][d] = rv1[d];
    }
    for(unsigned int d=0; d<dv2.size(); d++){       
       reco_track_dEdx_pl2[track_i][d] = dv2[d];
       reco_track_ResRan_pl2[track_i][d] = rv2[d];
    }


  }
  else {
    reco_track_daughter_nhits0[track_i][daughter_i] = hits_p0;
    reco_track_daughter_nhits1[track_i][daughter_i] = hits_p1;
    reco_track_daughter_nhits2[track_i][daughter_i] = hits_p2;
    reco_track_daughter_llrpid_3pl[track_i][daughter_i] = llr_pid_score;
    reco_track_daughter_llrpid_k_3pl[track_i][daughter_i] = llr_pid_score_k;
  
    cout << "llr_pid_score is " << reco_track_daughter_llrpid_3pl[track_i][daughter_i] << endl;

    //reco_track_daughter_dEdx = vector<vector<vector<Float_t>>>(kMaxTracks,vector<vector<Float_t>>(kMaxTracks,vector<Float_t>(hits_p2,-9999)));

    //dvv.resize(kMaxTracks,vector<Float_t>(hits_p2));
    //reco_track_daughter_dEdx.resize(kMaxTracks, dvv);
    //reco_track_daughter_dEdx[track_i][daughter_i] = dv;

    //rvv.resize(kMaxTracks,vector<Float_t>(hits_p2));
    //reco_track_daughter_ResRan.resize(kMaxTracks, rvv);
    //reco_track_daughter_ResRan[track_i][daughter_i] = rv;

    for(unsigned int d=0; d<dv0.size(); d++){       
       reco_track_daughter_dEdx_pl0[track_i][daughter_i][d] = dv0[d];
       reco_track_daughter_ResRan_pl0[track_i][daughter_i][d] = rv0[d];
    }
    for(unsigned int d=0; d<dv1.size(); d++){       
       reco_track_daughter_dEdx_pl1[track_i][daughter_i][d] = dv1[d];
       reco_track_daughter_ResRan_pl1[track_i][daughter_i][d] = rv1[d];
    }
    for(unsigned int d=0; d<dv2.size(); d++){       
       reco_track_daughter_dEdx_pl2[track_i][daughter_i][d] = dv2[d];
       reco_track_daughter_ResRan_pl2[track_i][daughter_i][d] = rv2[d];
    }


  }

}

/*
float GetLocalEFieldMag(const float x, const float y, const float z){

   const detinfo::DetectorProperties* detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
   auto const *sce = lar::providerFrom<spacecharge::SpaceChargeService>();

   double E_field_nominal = detprop->Efield();        // Electric Field in the drift region in KV/cm

   //correct Efield for SCE
   geo::Vector_t E_field_offsets = {0.,0.,0.};
   E_field_offsets = sce->GetCalEfieldOffsets(geo::Point_t{x,y, z});
   TVector3 E_field_vector = {E_field_nominal*(1 + E_field_offsets.X()), E_field_nominal*E_field_offsets.Y(), E_field_nominal*E_field_offsets.Z()};
   float E_field = E_field_vector.Mag();
   return E_field;
}


double CCKaonAnalyzerRebuild::ModBoxCorrection(const double dQdx, const float x, const float y, const float z) {
              
double rho = 1.383;//detprop->Density();            // LAr density in g/cm^3
double Wion = 23.6/1e6;//util::kGeVToElectrons;  // 23.6 eV = 1e, Wion in MeV/e
                          
auto E_field = GetLocalEFieldMag(x,y,z); // kV / cm
                                
double fModBoxA = 0.930;
double fModBoxB = 0.212;
                                            
double Beta = fModBoxB / (rho * E_field);
double Alpha = fModBoxA;
double dEdx = (exp(Beta * Wion * dQdx ) - Alpha) / Beta;
                                                              
return dEdx;
                                                                    
}

*/
void CCKaonAnalyzerRebuild::fillPID(const std::vector<art::Ptr<anab::ParticleID>> &trackPID, double angle_y, int track_i, int daughter_i)
{

  // check PID
  if (trackPID.size()==0) return;

  double chi2ka[3] = {-9.999,-9.999,-9.999};
  double chi2pr[3] = {-9.999,-9.999,-9.999};
  double chi2pi[3] = {-9.999,-9.999,-9.999};
  double chi2mu[3] = {-9.999,-9.999,-9.999};
  Float_t Bragg_fwd_ka[3] = {-1,-1,-1};
  Float_t Bragg_fwd_pr[3] = {-1,-1,-1};
  Float_t Bragg_fwd_pi[3] = {-1,-1,-1};
  Float_t Bragg_fwd_mu[3] = {-1,-1,-1};
  Float_t No_Bragg[3] = {-1,-1,-1};
  double likepr_3pl = -1;

  std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
  for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){


    anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
    if (AlgScore.fAlgName == "Chi2") {
      if (anab::kVariableType(AlgScore.fVariableType) == anab::kGOF) {
        for (int pl=0; pl<3; pl++) {
          if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==pl) {
            if (AlgScore.fAssumedPdg==321)  chi2ka[pl]=AlgScore.fValue;
            if (AlgScore.fAssumedPdg==2212) chi2pr[pl]=AlgScore.fValue;
            if (AlgScore.fAssumedPdg==211)  chi2pi[pl]=AlgScore.fValue;
            if (AlgScore.fAssumedPdg==13)   chi2mu[pl]=AlgScore.fValue;
          }
        }
      }
    }
    if (AlgScore.fAlgName == "ThreePlaneProtonPID") {
      if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood) {
        if (AlgScore.fPlaneMask==UBPID::uB_SinglePlaneGetBitset(2)) {
          if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward) {
            if (AlgScore.fAssumedPdg==2212) {
              likepr_3pl = AlgScore.fValue;
            }
          }
        }
      }
    }


    if (AlgScore.fAlgName == "BraggPeakLLH") {
      if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood) {
          //if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward) {
           
              
            for(int pl=0; pl<3; pl++){

            if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==pl) {
            if (TMath::Abs(AlgScore.fAssumedPdg)==2212) Bragg_fwd_pr[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==321) Bragg_fwd_ka[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==13) Bragg_fwd_mu[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==211) Bragg_fwd_pi[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==0) No_Bragg[pl] = AlgScore.fValue;
            }

            }
          //}
      }
    }

  }

  //double delta_z = end.Z()-pos.Z();
  //double delta_y = end.Y()-pos.Y();
  //double theta0 = TMath::ATan2(delta_z,delta_y) - TMath::Pi()/3;
  //double theta1 = TMath::ATan2(delta_z,delta_y) + TMath::Pi()/3.;
  //double theta2 = TMath::ATan2(delta_z,delta_y);
  double theta0 = angle_y - TMath::Pi()/3;
  double theta1 = angle_y + TMath::Pi()/3.;
  double theta2 = angle_y;
  double wpl0 = TMath::Power(TMath::Sin(theta0),2) >= 0.05 ? 1 : 0;
  double wpl1 = TMath::Power(TMath::Sin(theta1),2) >= 0.05 ? 1 : 0;
  double wpl2 = TMath::Power(TMath::Sin(theta2),2) >= 0.05 ? 1 : 0;
  double chi2ka_3pl = (wpl0*chi2ka[0] + wpl1*chi2ka[1] + wpl2*chi2ka[2])/(wpl0 + wpl1 + wpl2);
  double chi2pr_3pl = (wpl0*chi2pr[0] + wpl1*chi2pr[1] + wpl2*chi2pr[2])/(wpl0 + wpl1 + wpl2);
  double chi2pi_3pl = (wpl0*chi2pi[0] + wpl1*chi2pi[1] + wpl2*chi2pi[2])/(wpl0 + wpl1 + wpl2);
  double chi2mu_3pl = (wpl0*chi2mu[0] + wpl1*chi2mu[1] + wpl2*chi2mu[2])/(wpl0 + wpl1 + wpl2);

  if (daughter_i<0) {
    reco_track_chi2ka_pl0[track_i] = chi2ka[0];
    reco_track_chi2pr_pl0[track_i] = chi2pr[0];
    reco_track_chi2pi_pl0[track_i] = chi2pi[0];
    reco_track_chi2mu_pl0[track_i] = chi2mu[0];
    reco_track_chi2ka_pl1[track_i] = chi2ka[1];
    reco_track_chi2pr_pl1[track_i] = chi2pr[1];
    reco_track_chi2pi_pl1[track_i] = chi2pi[1];
    reco_track_chi2mu_pl1[track_i] = chi2mu[1];
    reco_track_chi2ka_pl2[track_i] = chi2ka[2];
    reco_track_chi2pr_pl2[track_i] = chi2pr[2];
    reco_track_chi2pi_pl2[track_i] = chi2pi[2];
    reco_track_chi2mu_pl2[track_i] = chi2mu[2];


    reco_track_Bragg_fwd_ka_pl0[track_i] = Bragg_fwd_ka[0];
    reco_track_Bragg_fwd_pr_pl0[track_i] = Bragg_fwd_pr[0];
    reco_track_Bragg_fwd_pi_pl0[track_i] = Bragg_fwd_pi[0];
    reco_track_Bragg_fwd_mu_pl0[track_i] = Bragg_fwd_mu[0];
    reco_track_Bragg_fwd_ka_pl1[track_i] = Bragg_fwd_ka[1];
    reco_track_Bragg_fwd_pr_pl1[track_i] = Bragg_fwd_pr[1];
    reco_track_Bragg_fwd_pi_pl1[track_i] = Bragg_fwd_pi[1];
    reco_track_Bragg_fwd_mu_pl1[track_i] = Bragg_fwd_mu[1];
    reco_track_Bragg_fwd_ka_pl2[track_i] = Bragg_fwd_ka[2];
    reco_track_Bragg_fwd_pr_pl2[track_i] = Bragg_fwd_pr[2];
    reco_track_Bragg_fwd_pi_pl2[track_i] = Bragg_fwd_pi[2];
    reco_track_Bragg_fwd_mu_pl2[track_i] = Bragg_fwd_mu[2];

    std::cout << "Track: " << track_i << "\t\tBraggPeakLLH: " << reco_track_Bragg_fwd_pr_pl2[track_i] <<std::endl;

    reco_track_MIP_pl0[track_i] = No_Bragg[0];
    reco_track_MIP_pl1[track_i] = No_Bragg[1];
    reco_track_MIP_pl2[track_i] = No_Bragg[2];

    reco_track_chi2ka_3pl[track_i] = chi2ka_3pl;
    reco_track_chi2pr_3pl[track_i] = chi2pr_3pl;
    reco_track_chi2pi_3pl[track_i] = chi2pi_3pl;
    reco_track_chi2mu_3pl[track_i] = chi2mu_3pl;
    reco_track_likepr_3pl[track_i] = likepr_3pl;
  }
  else {
    reco_track_daughter_chi2ka_pl0[track_i][daughter_i] = chi2ka[0];
    reco_track_daughter_chi2pr_pl0[track_i][daughter_i] = chi2pr[0];
    reco_track_daughter_chi2pi_pl0[track_i][daughter_i] = chi2pi[0];
    reco_track_daughter_chi2mu_pl0[track_i][daughter_i] = chi2mu[0];
    reco_track_daughter_chi2ka_pl1[track_i][daughter_i] = chi2ka[1];
    reco_track_daughter_chi2pr_pl1[track_i][daughter_i] = chi2pr[1];
    reco_track_daughter_chi2pi_pl1[track_i][daughter_i] = chi2pi[1];
    reco_track_daughter_chi2mu_pl1[track_i][daughter_i] = chi2mu[1];
    reco_track_daughter_chi2ka_pl2[track_i][daughter_i] = chi2ka[2];
    reco_track_daughter_chi2pr_pl2[track_i][daughter_i] = chi2pr[2];
    reco_track_daughter_chi2pi_pl2[track_i][daughter_i] = chi2pi[2];
    reco_track_daughter_chi2mu_pl2[track_i][daughter_i] = chi2mu[2];
    reco_track_daughter_chi2ka_3pl[track_i][daughter_i] = chi2ka_3pl;
    reco_track_daughter_chi2pr_3pl[track_i][daughter_i] = chi2pr_3pl;
    reco_track_daughter_chi2pi_3pl[track_i][daughter_i] = chi2pi_3pl;
    reco_track_daughter_chi2mu_3pl[track_i][daughter_i] = chi2mu_3pl;
    reco_track_daughter_likepr_3pl[track_i][daughter_i] = likepr_3pl;

    reco_track_daughter_Bragg_fwd_ka_pl0[track_i][daughter_i] = Bragg_fwd_ka[0];
    reco_track_daughter_Bragg_fwd_pr_pl0[track_i][daughter_i] = Bragg_fwd_pr[0];
    reco_track_daughter_Bragg_fwd_pi_pl0[track_i][daughter_i] = Bragg_fwd_pi[0];
    reco_track_daughter_Bragg_fwd_mu_pl0[track_i][daughter_i] = Bragg_fwd_mu[0];
    reco_track_daughter_Bragg_fwd_ka_pl1[track_i][daughter_i] = Bragg_fwd_ka[1];
    reco_track_daughter_Bragg_fwd_pr_pl1[track_i][daughter_i] = Bragg_fwd_pr[1];
    reco_track_daughter_Bragg_fwd_pi_pl1[track_i][daughter_i] = Bragg_fwd_pi[1];
    reco_track_daughter_Bragg_fwd_mu_pl1[track_i][daughter_i] = Bragg_fwd_mu[1];
    reco_track_daughter_Bragg_fwd_ka_pl2[track_i][daughter_i] = Bragg_fwd_ka[2];
    reco_track_daughter_Bragg_fwd_pr_pl2[track_i][daughter_i] = Bragg_fwd_pr[2];
    reco_track_daughter_Bragg_fwd_pi_pl2[track_i][daughter_i] = Bragg_fwd_pi[2];
    reco_track_daughter_Bragg_fwd_mu_pl2[track_i][daughter_i] = Bragg_fwd_mu[2];

    reco_track_daughter_MIP_pl0[track_i][daughter_i] = No_Bragg[0];
    reco_track_daughter_MIP_pl1[track_i][daughter_i] = No_Bragg[1];
    reco_track_daughter_MIP_pl2[track_i][daughter_i] = No_Bragg[2];



  }

/*
            Bragg_fwd_pr[0] = -9999;
            Bragg_fwd_ka[0] = -9999;
            Bragg_fwd_pi[0] = -9999;
            Bragg_fwd_mu[0] = -9999;
            No_Bragg[0] = -9999;

            Bragg_fwd_pr[1] = -9999;
            Bragg_fwd_ka[1] = -9999;
            Bragg_fwd_pi[1] = -9999;
            Bragg_fwd_mu[1] = -9999;
            No_Bragg[1] = -9999;

            Bragg_fwd_pr[2] = -9999;
            Bragg_fwd_ka[2] = -9999;
            Bragg_fwd_pi[2] = -9999;
            Bragg_fwd_mu[2] = -9999;
            No_Bragg[2] = -9999;

*/


}

void CCKaonAnalyzerRebuild::fillTrueMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_track,
                                      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                      int track_i,
                                      int daughter_i)
{

  //std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
  //art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);

  //cout << "hits_from_track.size(): " << hits_from_track.size() << endl;

  simb::MCParticle const* matched_mcparticle = NULL;
  std::unordered_map<int,double> trkide;
  double maxe=-1, tote=0;
  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
  for(size_t i_h=0; i_h<hits_from_track.size(); i_h++) {
    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hits_from_track[i_h].key(), particle_vec, match_vec);

    //cout << "track's particle_vec.size(): " << particle_vec.size() << endl;

    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
      //trkide[ particle_vec[i_p]->TrackId() ] ++;
      tote += match_vec[i_p]->energy;
      if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
        maxe = trkide[ particle_vec[i_p]->TrackId() ];
        matched_mcparticle = particle_vec[i_p];
      }
    }
  }
 

  if(matched_mcparticle){
    //const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
    TLorentzVector mcstart, mcend;
    unsigned int pstarti, pendi;

    if (daughter_i<0) {
      //cout << "this case daughter_i is: " << daughter_i << endl;
      reco_track_true_pdg[track_i] = matched_mcparticle->PdgCode();
      reco_track_true_origin[track_i] = 1;//int(mc_truth->Origin());
      reco_track_true_primary[track_i] = matched_mcparticle->Process()=="primary";
      reco_track_true_end_inTPC[track_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
      reco_track_true_end_in5cmTPC[track_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_true_end_inCCInclusiveTPC[track_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_true_length[track_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
    }
    else {
      reco_track_daughter_true_pdg[track_i][daughter_i] = matched_mcparticle->PdgCode();
      reco_track_daughter_true_origin[track_i][daughter_i] = 1;//int(mc_truth->Origin());
      reco_track_daughter_true_primary[track_i][daughter_i] = matched_mcparticle->Process()=="primary";
      reco_track_daughter_true_end_inTPC[track_i][daughter_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_end_in5cmTPC[track_i][daughter_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_end_inCCInclusiveTPC[track_i][daughter_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_length[track_i][daughter_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
      //for (auto const& pPart : ptList) {
      //  if (pPart->TrackId()==matched_mcparticle2->Mother()) {
      //    reco_track_daughter_true_mother[ntracks][ndaughters] = pPart->PdgCode();
      //    break;
      //  }
      //}
    }
  }else{
      cout << "no matched_mcparticle found in fillTrue" << endl;
  }

}


void CCKaonAnalyzerRebuild::fillTrackMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_track,
                                      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                      int track_i,
                                      int daughter_i)
{
  std::map<int,double> trkide;
  std::map<int,double> trkidhit;
  std::map<int,int> trkidpdg;    
  std::map<int,int> trkidmrgid;
  std::map<int,int> trkidmother;

  map<int, int>::iterator it_trkidpdg;
  map<int, int>::iterator it_trkidmrgid;
  map<int, int>::iterator it_trkidmrgid2;
  map<int, double>::iterator it_trkide;
  map<int, double>::iterator it_trkidhit;

  std::map<double, int, std::greater<int>> epdg;
  std::map<double, int, std::greater<int>> hitpdg;

  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
  for(size_t i_h=0; i_h<hits_from_track.size(); i_h++) {
    particle_vec.clear(); match_vec.clear();

    particles_per_hit.get(hits_from_track[i_h].key(), particle_vec, match_vec);

    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
      trkidhit[ particle_vec[i_p]->TrackId() ]++;
      trkidpdg[ particle_vec[i_p]->TrackId() ] = particle_vec[i_p]->PdgCode();
      trkidmrgid[ particle_vec[i_p]->TrackId() ] = -9;
      trkidmother[ particle_vec[i_p]->TrackId() ] = particle_vec[i_p]->Mother();
    }
  }


  //get merged id for two generations

  int currentMergedId = 1; 

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){
        
    if(trkidmrgid[it_trkidpdg->first] != -9) continue;  //not filled yet
    trkidmrgid[it_trkidpdg->first] = currentMergedId; //fill current id
    int currentMotherTrackId = trkidmother[it_trkidpdg->first];
    
    while (currentMotherTrackId > 0) {//have mother

      if(trkidpdg.find(currentMotherTrackId)==trkidpdg.end()) break;
      if(trkidpdg[currentMotherTrackId]!=it_trkidpdg->second) break;
      
	trkidmrgid[currentMotherTrackId] = currentMergedId;
	currentMotherTrackId = trkidmother[currentMotherTrackId];
      }
    ++currentMergedId; 
  }                                                                                             
 
  
  for (it_trkidmrgid = trkidmrgid.begin(); it_trkidmrgid != trkidmrgid.end(); it_trkidmrgid++){
    for (it_trkidmrgid2 = trkidmrgid.begin(); it_trkidmrgid2 != trkidmrgid.end(); it_trkidmrgid2++){
      
      if(it_trkidmrgid->first == it_trkidmrgid2->first) continue; //skip if they have the same trackID
      
      if(it_trkidmrgid->second == it_trkidmrgid2->second){ //check if they have same MergedId 
	//trkide.find( it_trkidmrgid->first );
	trkide[ it_trkidmrgid->first ] += trkide[ it_trkidmrgid2->first ];
	trkidhit[ it_trkidmrgid->first ] += trkidhit[ it_trkidmrgid2->first ];
	
	trkidmrgid[ it_trkidmrgid2->first ] = -1; //removed if it's been merged 
	trkide[ it_trkidmrgid2->first ] = -1;
	trkidhit[ it_trkidmrgid2->first ] = -1;
      }      
    }	
  }
  

  it_trkide = trkide.begin();
  it_trkidhit = trkidhit.begin();
  it_trkidmrgid = trkidmrgid.begin();

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){

    if(it_trkidmrgid->second ==-1) continue;
    epdg[ it_trkide->second ] = it_trkidpdg->second;
    hitpdg[ it_trkidhit->second ] = it_trkidpdg->second;

    it_trkide++;//? how to initialise
    it_trkidhit++;
    it_trkidmrgid++;
  }


  Int_t merge_i = 0;  
 
  /* 
cout << "TRKID AND MRGID" << endl;

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){
    if(trkidpdg.find( trkidmother[ it_trkidpdg->first] ) == trkidpdg.end()) continue;
    cout << "trkid: " << it_trkidpdg->first << ", mrgid: " << trkidmrgid[ it_trkidpdg->first] << ", PDG: " << it_trkidpdg->second << endl;
    cout << "trkid: " << it_trkidpdg->first << ", mother: " << trkidmother[ it_trkidpdg->first] << ", PDG: " << trkidpdg[ trkidmother[ it_trkidpdg->first] ] << endl;
  }

  cout << "ENERGY AND ITS PDG" << endl;
  */

  for (auto const& x : epdg){

    if(daughter_i<0){
      reco_track_match_e[track_i][merge_i] = x.first;
      //cout << "after fill reco_track_match_e[track_i][merge_i]: " << reco_track_match_e[track_i][merge_i] << endl;
      reco_track_match_epdg[track_i][merge_i] = x.second;

    }else{
      reco_track_daughter_match_e[track_i][daughter_i][merge_i] = x.first;
      reco_track_daughter_match_epdg[track_i][daughter_i][merge_i] = x.second;
    }
    merge_i++;   
  }
  
  merge_i=0;
  //cout << "HITS AND ITS PDG" << endl;
  for (auto const& y : hitpdg){
    if(daughter_i<0){ 
      reco_track_match_hit[track_i][merge_i] = y.first;     
      reco_track_match_hitpdg[track_i][merge_i] = y.second;
    }else{
      reco_track_daughter_match_e[track_i][daughter_i][merge_i] = y.first;
      reco_track_daughter_match_epdg[track_i][daughter_i][merge_i] = y.second;
    }
    merge_i++;
  }
  
}




void CCKaonAnalyzerRebuild::fillTrueMatching_sh(std::vector<art::Ptr<recob::Hit>>& hits_from_shower,
                                         art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                         int track_i,
                                         int daughter_i)
{

  simb::MCParticle const* matched_mcparticle = NULL;
  std::unordered_map<int,double> shwride;
  double maxe=-1, tote=0;
  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
  for(size_t i_h=0; i_h<hits_from_shower.size(); i_h++) {
    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hits_from_shower[i_h].key(), particle_vec, match_vec);

    //cout << "shower's particle_vec.size(): " << particle_vec.size() << endl;

    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      shwride[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
      tote += match_vec[i_p]->energy;

      //cout << "particle_vec[i_p]->TrackId(): " << particle_vec[i_p]->TrackId() << endl;
      //cout << "match_vec[i_p]->energy: " << match_vec[i_p]->energy << endl;
      //cout << "shwride[ particle_vec[i_p]->TrackId() ]: " << shwride[ particle_vec[i_p]->TrackId() ] << endl;
      if( shwride[ particle_vec[i_p]->TrackId() ] > maxe ){
        maxe = shwride[ particle_vec[i_p]->TrackId() ];
        matched_mcparticle = particle_vec[i_p];
      }
    }
  }
  //cout << "MAXE: " << maxe << endl;


  if(matched_mcparticle){
    //const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
    TLorentzVector mcstart, mcend;
    unsigned int pstarti, pendi;

    //cout << "found matched_mcparticle" << endl;
    //cout << "and its pdg: " << matched_mcparticle->PdgCode() << endl;

    if(track_i>=0 && daughter_i>=0){

      //cout << "filling reco_track_daughter_true_pdg_sh[track_i][daughter_i] with track_i: " << track_i << ", daughter_i: " << daughter_i << ", pdg" << matched_mcparticle->PdgCode() << endl;
      //cout << "before filling reco_track_daughter_true_pdg_sh[track_i][daughter_i]: " << reco_track_daughter_true_pdg_sh[track_i][daughter_i] << endl;
      //cout << "after filling reco_track_daughter_true_pdg_sh[track_i][daughter_i]: " << reco_track_daughter_true_pdg_sh[track_i][daughter_i] << endl;

      reco_track_daughter_true_pdg_sh[track_i][daughter_i] = matched_mcparticle->PdgCode();

      reco_track_daughter_true_origin_sh[track_i][daughter_i] = 1;//int(mc_truth->Origin());
      reco_track_daughter_true_primary_sh[track_i][daughter_i] = matched_mcparticle->Process()=="primary";
      reco_track_daughter_true_end_inTPC_sh[track_i][daughter_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_end_in5cmTPC_sh[track_i][daughter_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_end_inCCInclusiveTPC_sh[track_i][daughter_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_length_sh[track_i][daughter_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
      //}
    }
  }
}


void CCKaonAnalyzerRebuild::fillShowerMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_shower,
                                      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                      int track_i,
                                      int daughter_i)
{
  std::map<int,double> trkide;
  std::map<int,double> trkidhit;
  std::map<int,int> trkidpdg;
  std::map<int,int> trkidmrgid;
  std::map<int,int> trkidmother;

  map<int, int>::iterator it_trkidpdg;
  map<int, int>::iterator it_trkidmrgid;
  map<int, int>::iterator it_trkidmrgid2;
  map<int, double>::iterator it_trkide;
  map<int, double>::iterator it_trkidhit;

  std::map<double, int, std::greater<int>> epdg;
  std::map<double, int, std::greater<int>> hitpdg;

  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
  for(size_t i_h=0; i_h<hits_from_shower.size(); i_h++) {
    particle_vec.clear(); match_vec.clear();

    particles_per_hit.get(hits_from_shower[i_h].key(), particle_vec, match_vec);

    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
      trkidhit[ particle_vec[i_p]->TrackId() ]++;
      trkidpdg[ particle_vec[i_p]->TrackId() ] = particle_vec[i_p]->PdgCode();
      trkidmrgid[ particle_vec[i_p]->TrackId() ] = -9;
      trkidmother[ particle_vec[i_p]->TrackId() ] = particle_vec[i_p]->Mother();
    }
  }
 
  //get merged id for two generations

  int currentMergedId = 1; 

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){
        
    if(trkidmrgid[it_trkidpdg->first] != -9) continue;  //not filled yet
    trkidmrgid[it_trkidpdg->first] = currentMergedId; //fill current id
    int currentMotherTrackId = trkidmother[it_trkidpdg->first];
    
    while (currentMotherTrackId > 0) {//have mother

      if(trkidpdg.find(currentMotherTrackId)==trkidpdg.end()) break;
      if(trkidpdg[currentMotherTrackId]!=it_trkidpdg->second) break;
      
	trkidmrgid[currentMotherTrackId] = currentMergedId;
	currentMotherTrackId = trkidmother[currentMotherTrackId];
      }
    ++currentMergedId; 
  }                                                                                             
 
  
  for (it_trkidmrgid = trkidmrgid.begin(); it_trkidmrgid != trkidmrgid.end(); it_trkidmrgid++){
    for (it_trkidmrgid2 = trkidmrgid.begin(); it_trkidmrgid2 != trkidmrgid.end(); it_trkidmrgid2++){
      
      if(it_trkidmrgid->first == it_trkidmrgid2->first) continue; //skip if they have the same trackID
      
      if(it_trkidmrgid->second == it_trkidmrgid2->second){ //check if they have same MergedId 
	trkide[ it_trkidmrgid->first ] += trkide[ it_trkidmrgid2->first ];
	trkidhit[ it_trkidmrgid->first ] += trkidhit[ it_trkidmrgid2->first ];
	
	trkidmrgid[ it_trkidmrgid2->first ] = -1; //removed if it's been merged 
	trkide[ it_trkidmrgid2->first ] = -1;
	trkidhit[ it_trkidmrgid2->first ] = -1;
      }      
    }	
  }
  

  it_trkide = trkide.begin();
  it_trkidhit = trkidhit.begin();
  it_trkidmrgid = trkidmrgid.begin();

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){

    if(it_trkidmrgid->second ==-1) continue;
    epdg[ it_trkide->second ] = it_trkidpdg->second;
    hitpdg[ it_trkidhit->second ] = it_trkidpdg->second;

    it_trkide++;
    it_trkidhit++;
    it_trkidmrgid++;
  }

  Int_t merge_i = 0; 

  /*
  cout << "TRKID AND MRGID" << endl;

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){
    cout << "shwid: " << it_trkidpdg->first << ", mrgid: " << trkidmrgid[ it_trkidpdg->first] << ", PDG: " << it_trkidpdg->second << endl;
    if(trkidpdg.find( trkidmother[ it_trkidpdg->first] ) == trkidpdg.end()) continue;
    cout << "shwid: " << it_trkidpdg->first << ", mother: " << trkidmother[ it_trkidpdg->first] << ", PDG: " << trkidpdg[ trkidmother[ it_trkidpdg->first] ] << endl;
  }

  cout << "ENERGY AND ITS PDG" << endl;
  */

  for (auto const& x : epdg){

    if(track_i>=0 && daughter_i>=0){
      //cout << "track_i: " << track_i << ", daughter_i: " << daughter_i << ", merge_i: " << merge_i << endl;
      //cout << "before fill reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i]: " << reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i] << endl;
      reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i] = x.first;
      //cout << "after fill reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i]: " << reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i] << endl;

      reco_track_daughter_shower_match_epdg[track_i][daughter_i][merge_i] = x.second;
    }
    merge_i++;   
  }
  
  merge_i=0;

  //cout << "HITS AND ITS PDG" << endl;

  for (auto const& y : hitpdg){
    if(track_i>=0 && daughter_i>=0){ 
      reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i] = y.first;
      reco_track_daughter_shower_match_epdg[track_i][daughter_i][merge_i] = y.second;
    }
    merge_i++;
  }
  
}




///////////////////////////////////////////////////////////////////////////////////////////

//========================================================================

    void CCKaonAnalyzerRebuild::ClearVertex(){

        //------------ Event related Variables -------------
        m_event_number = -99;
        m_subrun_number = -99;
        m_run_number = -99;
        m_test_matched_hits = 0;

        m_genie_spline_weight = 1.0;

        //------------ Vertex related Variables -------------
        m_reco_vertex_size = 0;
        m_vertex_pos_x=-99999;
        m_vertex_pos_y=-99999;
        m_vertex_pos_z=-99999;
        m_vertex_pos_tick=-9999;
        m_vertex_pos_wire_p0=-9999;
        m_vertex_pos_wire_p1=-9999;
        m_vertex_pos_wire_p2=-9999;

        m_reco_vertex_to_nearest_dead_wire_plane0=-99999;
        m_reco_vertex_to_nearest_dead_wire_plane1=-99999;
        m_reco_vertex_to_nearest_dead_wire_plane2=-99999;

        m_reco_slice_objects = 0;


        //this->ClearIsolation();

        //this->ClearSecondShowers();
        //------------- Flash related Variables ------------------
        this->ClearFlashes();

        //------------- Track Related Variables -----------------
        this->ClearTracks();

        //------------- Track Related Variables -----------------
        this->ClearShowers();
        this->ClearMCTruths();

        //------------- EventWeight Related Variables -----------------
        //this->ClearEventWeightBranches();



        //MetaData Related Varibles
        this->ClearSlices();


    }


  void CCKaonAnalyzerRebuild::GetVertex(const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, const art::Ptr<recob::PFParticle> & particle ){

        if(m_is_verbose) std::cout<<"SinglePhoton::Getvertex()\t||\t Starting to analyze recob::Vertex\n";
        int n_vert =0;

        //std::cout<<"There are "<<pfParticlesToVerticesMap.count(particle)<<" verticies associated with this particle"<<std::endl;

        lar_pandora::PFParticlesToVertices::const_iterator vIter = pfParticlesToVerticesMap.find(particle);
        if (pfParticlesToVerticesMap.end() != vIter)
        {
            const lar_pandora::VertexVector &vertexVector = vIter->second;
            if (!vertexVector.empty())
            {
                if (vertexVector.size() !=1)
                    std::cout << " Warning: Found particle with more than one associated vertex " << "\n";

                const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
                double xyz[3] = {0.0, 0.0, 0.0} ;
                vertex->XYZ(xyz);

                n_vert++;
                //std::cout<<"Vertex!"<<"\t "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<"\n";

                m_vertex_pos_x = xyz[0];
                m_vertex_pos_y = xyz[1];
                m_vertex_pos_z = xyz[2];

                if(!m_run_pi0_filter){
                    m_reco_vertex_to_nearest_dead_wire_plane0 = distanceToNearestDeadWire(0, m_vertex_pos_y, m_vertex_pos_z,geom,bad_channel_list_fixed_mcc9);
                    m_reco_vertex_to_nearest_dead_wire_plane1 = distanceToNearestDeadWire(1, m_vertex_pos_y, m_vertex_pos_z,geom,bad_channel_list_fixed_mcc9);
                    m_reco_vertex_to_nearest_dead_wire_plane2 = distanceToNearestDeadWire(2, m_vertex_pos_y, m_vertex_pos_z,geom,bad_channel_list_fixed_mcc9);
                }

            }else{
                std::cout << " Error: vertexVector associated with this particle is empty " << "\n";
                std::cerr << " Error: vertexVector associated with this particle is empty " << "\n";
                //exit(0);

            }
        }

        if(m_is_verbose) std::cout<<"SinglePhoton::Getvertex()\t||\t Finished. Found "<<n_vert<<" vertices.\n";
    }


  void CCKaonAnalyzerRebuild::GetPFParticleIdMap(const Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleHandle &pfParticleHandle, Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleIdMap &pfParticleMap)
//  void GetPFParticleIdMap(const Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleHandle &pfParticleHandle, Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleIdMap &pfParticleMap)
  {
    //std::cout<<"Filling pfParticleMap with from the handle with total number "<<pfParticleHandle->size()<<std::endl;                                                                                                                                                                                                                      
    for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
      {
	const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
	// std::cout<<"Adding PFP to pfParticleMap with pfp id  "<<pParticle->Self()<<std::endl;                                                                                                                                                                                                                                            
	if (!pfParticleMap.insert(Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second)
	  {
	    throw cet::exception("SinglePhoton") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
	  }
      }
  }


  void CCKaonAnalyzerRebuild::GetFinalStatePFParticleVectors(const Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleIdMap &pfParticleMap, const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleVector &crParticles, Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleVector &nuParticles )
    {

        int found = 0;
        int primaries = 0;
        int full = 0;
        for (Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)
        {
            const art::Ptr<recob::PFParticle> pParticle(it->second);

            full++;
            // Only look for primary particles
            if (!pParticle->IsPrimary()) continue;

            // Check if this particle is identified as the neutrino
            const int pdg(pParticle->PdgCode());
            const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);


            primaries++;
            // If it is, lets get the vertex position
            if(isNeutrino){
                found++;
                this->GetVertex(pfParticlesToVerticesMap, pParticle );

            }

            // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
            if (!isNeutrino)
            {
                crParticles.push_back(pParticle);
                continue;
            }

            // ATTN. We are filling nuParticles under the assumption that there is only one reconstructed neutrino identified per event.
            //       If this is not the case please handle accordingly
            if (!nuParticles.empty())
            {
                throw cet::exception("SinglePhoton") << "  This event contains multiple reconstructed neutrinos!";
            }

            // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
            for (const size_t daughterId : pParticle->Daughters())
            {
                if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                    throw cet::exception("SinglePhoton") << "  Invalid PFParticle collection!";

                nuParticles.push_back(pfParticleMap.at(daughterId));
            }
        }
        std::cout<<"SinglePhoton::GetFinalStatePFParticleVectors()\t||\t Found "<<primaries<<" primary PFParticles (out of "<<full<<") of which: "<<found<<" were neutrinos."<<std::endl;
        m_reco_vertex_size = found;

    }

 
  void CCKaonAnalyzerRebuild::CollectTracksAndShowers(const Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleVector &particles,const Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleIdMap pfParticleMap, const Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleHandle &pfParticleHandle, const art::Event &evt, Kaon_Analyzer::CCKaonAnalyzerRebuild::TrackVector &tracks, Kaon_Analyzer::CCKaonAnalyzerRebuild::ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap)
  {
    
    
    // Get the associations between PFParticles and tracks/showers from the event
    art::FindManyP< recob::Track     > pfPartToTrackAssoc(pfParticleHandle, evt, m_trackLabel);
    art::FindManyP< recob::Shower    > pfPartToShowerAssoc(pfParticleHandle, evt, m_showerLabel);
    
    //if running over the neutrino slice only 
    if (m_run_all_pfps == false){ 
      for (const art::Ptr<recob::PFParticle> &pParticle : particles) {
	const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
	const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));
	
	CCKaonAnalyzerRebuild::FillTracksAndShowers(associatedTracks, associatedShowers, pParticle,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);
      }
    } else{ //if running over all slices
      std::cout<<"SinglePhoton\t||\tThe total number of PFP's in the map is "<<pfParticleMap.size()<<std::endl;
      //            std::cout<<"The total number of PFP's in the vector is "<< particles.size()<<std::endl;
      for (auto pair : pfParticleMap){
	const art::Ptr<recob::PFParticle> &pParticle = pair.second;
	
	const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
	const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));
	
	CCKaonAnalyzerRebuild::FillTracksAndShowers(associatedTracks, associatedShowers, pParticle,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);
	
      }
    }
    
    
  }
  
 
  void CCKaonAnalyzerRebuild::FillTracksAndShowers( const std::vector< art::Ptr<recob::Track> > & associatedTracks, const std::vector< art::Ptr<recob::Shower> > & associatedShowers, const art::Ptr<recob::PFParticle> &pParticle , const Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleHandle &pfParticleHandle, const art::Event &evt, Kaon_Analyzer::CCKaonAnalyzerRebuild::TrackVector &tracks, Kaon_Analyzer::CCKaonAnalyzerRebuild::ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap)
  {

    const unsigned int nTracks(associatedTracks.size());
    const unsigned int nShowers(associatedShowers.size());
    //const unsigned int 
    //nTracks(associatedTracks.size());
    //const unsigned int 
    //nShowers(associatedShowers.size());


    // Check if the PFParticle has no associated tracks or showers                                                                                                                                                                                                                                                                                       
    if (nTracks == 0 && nShowers == 0)
      {
	mf::LogDebug("SinglePhoton") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << "\n";
	return;
      }

    // Check if there is an associated track                                                                                                                                                                                                                                                                                                             
    if (nTracks == 1 && nShowers == 0)
      {
	tracks.push_back(associatedTracks.front());
	trackToNuPFParticleMap[tracks.back()]= pParticle;
	return;
      }

    // Check if there is an associated shower                                                                                                                                                                                                                                                                                                            
    if (nTracks == 0 && nShowers == 1)
      {
	showers.push_back(associatedShowers.front());
	showerToNuPFParticleMap[showers.back()] = pParticle;
	return;
      }

    throw cet::exception("CCKaonAnalyzerRebuild") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();

  }


    double CCKaonAnalyzerRebuild::triangle_area(double a1, double a2, double b1, double b2, double c1, double c2){
        return fabs((a1*(b2-c2)+b1*(c2-a2)+c1*(a2-b2))/2.0);
    }

    int CCKaonAnalyzerRebuild::quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area){

        std::vector<double> z(n,0.0);

        TGraph2D *g = new TGraph2D(n,X,Y,&z[0]);
        TGraphDelaunay delan(g);
        delan.SetMarginBinsContent(0);
        delan.ComputeZ(0,0);
        delan.FindAllTriangles();
        (*num_triangles)=delan.GetNdt();

        //Grab the locations of all the trianges. These will be intergers referencing to position in X,Y arrays
        Int_t *MT = delan.GetMTried();
        Int_t *NT = delan.GetNTried();
        Int_t *PT = delan.GetPTried();

        (*area)=0.0;
        for(int i = 0; i<delan.GetNdt(); i++){
            (*area)+=triangle_area(X[MT[i]-1],Y[MT[i]-1],X[NT[i]-1],Y[NT[i]-1],X[PT[i]-1],Y[PT[i]-1]);
        }

        delete g;
        return 0;
    }

    int CCKaonAnalyzerRebuild::delaunay_hit_wrapper(const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<int> & num_hits, std::vector<int>& num_triangles, std::vector<double> & area){

        int n = hits.size();
        std::vector<double> C0,T0;
        std::vector<double> C1,T1;
        std::vector<double> C2,T2;
        size_t n_0=0;
        size_t n_1=0;
        size_t n_2=0;

        for(int i=0;i<n; i++){
            const art::Ptr<recob::Hit> hit = hits[i];
            switch(hit->View()){
                case 0:
                    C0.push_back((double)hit->Channel());         
                    T0.push_back(hit->PeakTime());         
                    n_0++;
                    break;
                case 1:
                    C1.push_back((double)hit->Channel());         
                    T1.push_back(hit->PeakTime());         
                    n_1++;
                    break;
                case 2:
                    C2.push_back((double)hit->Channel());         
                    T2.push_back(hit->PeakTime());         
                    n_2++;
                    break;
                default:
                    break;
            }
        }
        if(m_use_delaunay){
            if(n_0>0) this->quick_delaunay_fit(n_0, &C0[0]  , &T0[0]  , &num_triangles[0],&area[0]);
            if(n_1>0) this->quick_delaunay_fit(n_1, &C1[0]  , &T1[0]  , &num_triangles[1],&area[1]);
            if(n_2>0) this->quick_delaunay_fit(n_2, &C2[0]  , &T2[0]  , &num_triangles[2],&area[2]);
        }
        num_hits[0] = n_0;
        num_hits[1] = n_1;
        num_hits[2] = n_2;

        //std::cout<<"Plane 0: "<<n_0<<" hits with "<<num_triangles[0]<<" triangles of area: "<< area[0]<<std::endl;
        //std::cout<<"Plane 1: "<<n_1<<" hits with "<<num_triangles[1]<<" triangles of area: "<< area[1]<<std::endl;
        //std::cout<<"Plane 2: "<<n_2<<" hits with "<<num_triangles[2]<<" triangles of area: "<< area[2]<<std::endl;

        return 0;
    }

    int CCKaonAnalyzerRebuild::spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected, std::vector<double> & input){
        corrected.resize(3);

        double kx = input[0];
        double ky = input[1];
        double kz = input[2];

        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks->TPCG4Time2Tick(mcparticle->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();

        double xtimeoffset = theDetector->ConvertTicksToX(g4Ticks,0,0,0);

        //        double xOffset = -scecorr.X() +xtimeoffset+0.6;
        double yOffset = scecorr.Y();
        double zOffset = scecorr.Z();

        corrected[0]=kx - scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
        corrected[1]=ky+yOffset;
        corrected[2]=kz+zOffset;

        //std::cout<<"CCKaonAnalyzerRebuild\t||\tTRIGGER_OFF: "<<kx<<" "<<xOffset<<" "<<theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)<<" "<<scecorr.X()<<std::endl;
        //std::cout<<"CCKaonAnalyzerRebuild\t||\tTRIGGER_OFF: "<<xOffset<<" "<<yOffset<<" "<<zOffset<<std::endl;
        //std::cout<<"CCKaonAnalyzerRebuild\t||\tTRIGGER_OFF: mcp->T(): "<<mcparticle->T()<<" TPCG4Time2Tick(): "<<detClocks->TPCG4Time2Tick(mcparticle->T())<<". "<<theDetector->GetXTicksOffset(0,0,0)<<" "<<theDetector->TriggerOffset()<<std::endl;
        return 0;
    }





    int CCKaonAnalyzerRebuild::spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected){
        corrected.resize(3);

        double kx = mcparticle->Vx();
        double ky = mcparticle->Vy();
        double kz = mcparticle->Vz();

        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks->TPCG4Time2Tick(mcparticle->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();

        double xtimeoffset = theDetector->ConvertTicksToX(g4Ticks,0,0,0);

        //double xOffset = -scecorr.X() +xtimeoffset+0.6;
        double yOffset = scecorr.Y();
        double zOffset = scecorr.Z();

        corrected[0]=kx - scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
        corrected[1]=ky+yOffset;
        corrected[2]=kz+zOffset;

        //std::cout<<"CCKaonAnalyzerRebuild\t||\tTRIGGER_OFF: "<<kx<<" "<<xOffset<<" "<<theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)<<" "<<scecorr.X()<<std::endl;
        //std::cout<<"CCKaonAnalyzerRebuild\t||\tTRIGGER_OFF: "<<xOffset<<" "<<yOffset<<" "<<zOffset<<std::endl;
        //std::cout<<"CCKaonAnalyzerRebuild\t||\tTRIGGER_OFF: mcp->T(): "<<mcparticle->T()<<" TPCG4Time2Tick(): "<<detClocks->TPCG4Time2Tick(mcparticle->T())<<". "<<theDetector->GetXTicksOffset(0,0,0)<<" "<<theDetector->TriggerOffset()<<std::endl;
        return 0;
    }










    int CCKaonAnalyzerRebuild::spacecharge_correction(const simb::MCParticle & mcparticle, std::vector<double> & corrected){
        corrected.resize(3);
        //Space Charge Effect! functionize this soon.
        double kx = mcparticle.Vx();
        double ky = mcparticle.Vy();
        double kz = mcparticle.Vz();

        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks->TPCG4Time2Tick(mcparticle.T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();

        double xtimeoffset = theDetector->ConvertTicksToX(g4Ticks,0,0,0);

        corrected[0]=kx - scecorr.X() +xtimeoffset+0.6;
        corrected[1]=ky + scecorr.Y();
        corrected[2]=kz + scecorr.Z();
        return 0;
    }

    void CCKaonAnalyzerRebuild::CollectMCParticles(const art::Event &evt, const std::string &label, std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> &truthToParticles,        std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>              &particlesToTruth, std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap)
    {

        //    if (evt.isRealData())
        //      throw cet::exception("LArPandora") << " PandoraCollector::CollectMCParticles --- Trying to access MC truth from real data ";

        art::Handle< std::vector< simb::MCParticle>  > theParticles;
        evt.getByLabel(label, theParticles);

        if (!theParticles.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles " << std::endl;
        }

        art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, evt, label);

        for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i)
        {
            const art::Ptr<simb::MCParticle> particle(theParticles, i);
            const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));
            truthToParticles[truth].push_back(particle);
            particlesToTruth[particle] = truth;
            MCParticleToTrackIdMap[particle->TrackId()] = particle;
        }

        std::cout<<"CCKaonAnalyzerRebuild::CollectMCParticles() \t||\t the number of MCParticles in the event is "<<theParticles->size()<<std::endl;
    }

    void CCKaonAnalyzerRebuild::CollectSimChannels(const art::Event &evt, const std::string &label,  std::vector< art::Ptr<sim::SimChannel> >  &simChannelVector)
    {
        //    if (evt.isRealData())
        //      throw cet::exception("LArPandora") << " PandoraCollector::CollectSimChannels --- Trying to access MC truth from real data ";

        art::Handle< std::vector<sim::SimChannel> > theSimChannels;
        evt.getByLabel(label, theSimChannels);

        if (!theSimChannels.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find sim channels... " << std::endl;
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theSimChannels->size() << " SimChannels " << std::endl;
        }

        for (unsigned int i = 0; i < theSimChannels->size(); ++i)
        {
            const art::Ptr<sim::SimChannel> channel(theSimChannels, i);
            simChannelVector.push_back(channel);
        }
    }


    void CCKaonAnalyzerRebuild::BuildMCParticleHitMaps(const art::Event &evt, const std::string &label, const std::vector<art::Ptr<recob::Hit>> &hitVector,   std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  &particlesToHits,         std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  &hitsToParticles, const lar_pandora::LArPandoraHelper::DaughterMode daughterMode, std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap)
    {
        std::vector< art::Ptr<sim::SimChannel> >   simChannelVector;
        std::map< art::Ptr<simb::MCTruth>,     std::vector<art::Ptr<simb::MCParticle>>  >    truthToParticles;
        std::map< art::Ptr<simb::MCParticle>,  art::Ptr<simb::MCTruth> > particlesToTruth;
        std::map< art::Ptr<recob::Hit>,    std::vector< sim::TrackIDE >    >               hitsToTrackIDEs;

        this->CollectSimChannels(evt, label, simChannelVector);
        this->CollectMCParticles(evt, label, truthToParticles, particlesToTruth, MCParticleToTrackIdMap);
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(hitVector, simChannelVector, hitsToTrackIDEs);
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(hitsToTrackIDEs, truthToParticles, particlesToHits, hitsToParticles, daughterMode);


    }



  // DEFINE_ART_MODULE(CCKaonAnalyzerRebuild)
}
//}
