#ifndef CCKAON_PRODUCTION
#define CCKAON_PRODUCTION

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/View.h"
#include "art/Utilities/make_tool.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h" 

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h" 
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "/uboone/app/users/taniuchi/51_pandora/srcs/larpandoracontent/larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h" 

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Deprecated/BezierTrack.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larcoreobj/SummaryData/POTSummary.h"


#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"


#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"

//#include "CCKaonAnalyzer_module.h"
/*
#include "headers/analyze_Slice.h"
#include "headers/analyze_Tracks.h"
#include "headers/analyze_Showers.h"
#include "headers/analyze_MCTruth.h"
#include "headers/analyze_OpFlashes.h"
#include "headers/particle_split_basetool_producer.h"
#include "headers/track_production_producer.h"
*/

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "PID/LLR_PID.h"
#include "PID/LLRPID_proton_muon_lookup.h"
#include "PID_K/LLR_PID_K.h"
#include "PID_K/LLRPID_kaon_proton_lookup.h"

#include "Pandora/PdgTable.h" 

#include "TCanvas.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TTimeStamp.h"

#include <iostream>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TVectorT.h>
#include <THStack.h>
#include <TPDF.h>
#include <TLegend.h>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <array>
#include <map>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <string>
#include <numeric>
#include <algorithm>
#include <functional>
#include <typeinfo>
#include <memory>
#include <chrono>
#include <sys/stat.h>


//#ifdef __MAKECINT__
#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+; 
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+; 
#endif

const int kMaxTracks=20;
//const int kMaxTracks=50;

using namespace std;
namespace Kaon_Analyzer{


  class CCKaonProducer : public art::EDProducer {

  public:
    explicit CCKaonProducer(fhicl::ParameterSet const& pset);
    virtual ~CCKaonProducer();

    void beginJob();
    void produce(art::Event& evt) override;
    //void produce(const art::Event& evt) override;


     void fillAngularDistributionMapCheat(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
				          TVector3 Kend_candidate,
				          art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                          art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
                                          std::map<int, std::map<int, double>>& angular_distribution_mrgid_map,
					  std::map<int,int> &mrgidpdg,
					  int& currentMergedId);

     void fillAngularDistributionMapCheat3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
					    std::vector<art::Ptr<recob::Hit>>& true_pi0_hit_list,
					    TVector3 Kend_candidate,
					    art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
					    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
					    std::map<int, std::map<int, std::map<int, double>>>& angular_distribution_mrgid_map_3D,
					    std::map<int,int> &mrgidpdg,
					    std::map<int,TVector3> &mrgidmom,
					    int& currentMergedId);

     void fillAngularDistributionMap(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
				     TVector3 Kend_candidate,
                                     art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
                                     std::map<int, double>& angular_distribution_map);


     void fillAngularDistributionMap3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
				       TVector3 Kend_candidate,
				       art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
				       std::map<int, std::map<int, double>> &angular_distribution_map_3D);

     int fillHistAngularDistributionMapCheat( std::map<int, std::map<int, double>>& angular_distribution_mrgid_map,
					      std::map<int,int>& mrgidpdg,
					      std::vector<TH1D*>& h_angular_distribution_pfparticle_cheat,
					      vector<int>& v_pdg);


     int fillHistAngularDistributionMapCheat3D( std::map<int, std::map<int, std::map<int, double>>>& angular_distribution_mrgid_map_3D,
						std::map<int,int>& mrgidpdg,
						std::vector<TH2D*>& h_angular_distribution_pfparticle_cheat_3D,
						vector<int>& v_pdg);
     
     int fillHistAngularDistributionMap( std::map<int, double>& angular_distribution_map,
					 vector<TH1D*>& h_angular_distribution_pfparticle,
					 vector<bool>& v_trk_flg,
					 bool trk_flg);

     int fillHistAngularDistributionMap3D( std::map<int, std::map<int, double>>& angular_distribution_map_3D,
					   vector<TH2D*>& h_angular_distribution_pfparticle_3D,
					   vector<bool>& v_trk_flg,
					   bool trk_flg);

     int fillHistAngularDistributionMapSurface( std::map<int, std::map<int, double>>& angular_distribution_map_3D,
					       vector<TH2D*>& h_angular_distribution_pfparticle_surface,
					       vector<bool>& v_trk_flg,
					       bool trk_flg);

     int fillHistAngularDistributionMapSphere( std::map<int, std::map<int, double>>& angular_distribution_map_3D,
					       vector<TH3D*>& h_angular_distribution_pfparticle_sphere,
					       vector<bool>& v_trk_flg,
					       bool trk_flg);

     void smoothAngularDistributionMap(std::map<int, double> &angular_distribution_map);

     void smoothAngularDistributionMapCheat3D(std::map<int, std::map<int, std::map<int, double>>> &angular_distribution_mrgid_map_3D);

     void smoothAngularDistributionMap3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D);

     void accumulateAngularDistributionMap3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D_1,
					     std::map<int, std::map<int, double>> &angular_distribution_map_3D_2,
					     std::map<int, std::map<int, double>> &angular_distribution_map_3D);

     void obtainPeakVectorCheat3D(std::map<int, std::map<int, std::map<int, double>>> &angular_distribution_mrgid_map_3D, 
				  std::map<int,int>& mrgidpdg,
				  std::map<int, std::map<double, TVector2, std::greater<>>>& view_peak_map_cheat,
				  vector<int>& v_pdg_peak);

     void obtainPeakVector3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D,
			     vector<bool>& v_trk_flg,
			     std::map<double, TVector2, std::greater<>>& view_peak_map,
			     bool trk_flg);

     void findBestAngularPeakCheat3D( std::map<int, std::map<double, TVector2, std::greater<>>>& view_peak_map_cheat,
				     std::map<int, vector<TVector2>> &best_peak_bins_cheat);


     void findBestAngularPeak3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D,
				std::map<double, TVector2, std::greater<>>& view_peak_vector,
				vector<TVector2> &best_peak_bins);

     void findShowerSpine3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
			    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			    std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			    std::vector<std::vector<art::Ptr<recob::Hit>>>& shower_spine_hit_list_vector,
			    TVector3 Kend_candidate,
			    std::map<double, TVector2, std::greater<>>& view_peak_map,
			    vector<TVector2> &best_peak_bins);

     void findShowerSpine3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
			    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			    std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			    std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
			    TVector3 Kend_candidate,
			    std::map<double, TVector2, std::greater<>>& view_peak_map,
			    vector<TVector2> &best_peak_bins);

     void findShowerSpine3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
			    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			    std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			    std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
			    TVector3 Kend_candidate,
			    std::map<double, TVector2, std::greater<>>& view_peak_map,
			    vector<TVector2> &best_peak_bins,
			    double true_length);  
        
     bool collectSubsectionHits(const lar_content::ThreeDSlidingFitResult &extrapolated_fit,
				const TVector3 &extrapolated_start_position,
				const TVector3 &extrapolated_end_position,
				const TVector3 &extrapolated_direction,
				const bool is_end_downstream,
				std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
				art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
				vector<TVector3> &running_fit_position_vec,
				pandora::CartesianPointVector &pandora_running_fit_position_vec,
				std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
				std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list);
     

     bool collectSubsectionHits(const lar_content::ThreeDSlidingFitResult &extrapolated_fit,
				const TVector3 &extrapolated_start_position,
				const TVector3 &extrapolated_end_position,
				const TVector3 &extrapolated_direction,
				const bool is_end_downstream,
				std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
				art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
				vector<TVector3> &running_fit_position_vec,
				pandora::CartesianPointVector &pandora_running_fit_position_vec,
				std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
				std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
				TVector3 Kend_candidate,
				double true_length);
     
     void collectConnectedHits(std::vector<art::Ptr<recob::Hit>>& collected_hits,
			       art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			       const TVector3 &extrapolated_start_position,
			       const TVector3 &extrapolated_direction,
			       vector<TVector3> &running_fit_position_vec,
			       pandora::CartesianPointVector &pandora_running_fit_position_vec,
			       std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list);

     TVector3 getClosestPointToLine3D(const TVector3 &extrapolated_start_position,
				      const TVector3 &extrapolated_direction,
				      art::Ptr<recob::Hit>& collected_hit,
				      const TVector3 &hit_position);

     bool isCloseToLine(const TVector3 &hit_position,
			const TVector3 &line_start,
			const TVector3 &line_direction,
			const double distance_to_line);

     double getClosestDistance(art::Ptr<recob::Hit>& collected_hit,
			       art::Ptr<recob::Hit>& shower_spine_hit,
			       art::FindManyP<recob::SpacePoint>& spacepoint_per_hit);

     double getClosestDistance(art::Ptr<recob::Hit>& collected_hit,
			       std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
			       art::FindManyP<recob::SpacePoint>& spacepoint_per_hit);

     double getClosestDistance(const TVector3 &position,
			       vector<TVector3> &test_positions);

     void obtainLength(std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
		       art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
		       TVector3 Kend_candidate);

     void obtainLongitudinalDecomposition(std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
					  art::FindManyP<recob::SpacePoint>& spacepoint_per_hit);
     
     int drawHistAngularDistributionMapCheat( std::vector<TH1D*> &h_angular_distribution_pfparticle_cheat,
					      std::vector<int> &v_pdg,
					      TString outfile_name,
					      TCanvas*& c);

     int drawHistAngularDistributionMapCheat3D( std::vector<TH2D*> &h_angular_distribution_pfparticle_cheat_3D,
						std::vector<int> &v_pdg,
						TString outfile_name,
						TCanvas*& c);
     
     int drawHistAngularDistributionMap( std::vector<TH1D*>& h_angular_distribution_pfparticle,
					 std::vector<bool>& v_trk_flg,
					 TString outfile_name,
					 TCanvas*& c);

     int drawHistAngularDistributionMap3D( std::vector<TH2D*>& h_angular_distribution_pfparticle_3D,
					   std::vector<bool>& v_trk_flg,
					   TString outfile_name,
					   TCanvas*& c);

     int drawHistAngularDistributionMapSurface( std::vector<TH2D*> &h_angular_distribution_pfparticle_surface,
						std::vector<bool> &v_trk_flg,
						TString outfile_name,
						TCanvas* &c);

     int drawHistAngularDistributionMapSphere( std::vector<TH3D*> &h_angular_distribution_pfparticle_sphere,
					       std::vector<bool> &v_trk_flg,
					       TString outfile_name,
					       TCanvas* &c);

     void fillTrueMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_track, 
                           art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit, 
                           int track_i=-1, 
                           int daughter_i=-1); 
     void fillTrueMatching_sh(std::vector<art::Ptr<recob::Hit>>& hits_from_shower, 
                           art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit, 
                           int track_i=-1, 
                           int daughter_i=-1);
     void fillTrackMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_track, 
                           art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit, 
                           int track_i=-1, 
                           int daughter_i=-1); 

     void fillShowerMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_shower, 
                           art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit, 
                           int track_i=-1, 
                           int daughter_i=-1); 



  private:
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fShowerModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fPIDLabel;
    std::string fHitTruthAssns;
    std::string fHitTrackAssns;
    std::string fHitShowerAssns;
    std::string m_pfp_producer;
    std::string fPFParticleLabel;
    std::string fSpacePointproducer;
    std::string fPandoraLabel;
    bool isMC;

    detinfo::DetectorProperties const * theDetector ;// = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detinfo::DetectorClocks    const *  detClocks   ;//= lar::providerFrom<detinfo::DetectorClocksService>();  
    spacecharge::SpaceCharge const * SCE;
    geo::GeometryCore const * geom;

    Bool_t  reco_nu_cc_filter;

    Float_t reco_nu_vtx_x;
    Float_t reco_nu_vtx_y;
    Float_t reco_nu_vtx_z;
    Int_t   reco_nu_ndaughters;
    Int_t   reco_nu_cc_nmue;

    Float_t reco_ccmu_vtx_x;
    Float_t reco_ccmu_vtx_y;
    Float_t reco_ccmu_vtx_z;

    Int_t   reco_track_true_pdg[kMaxTracks];
 
  };

  /*
  CCKaonProducer::CCKaonProducer(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    produces< std::vector<recob::Track> >();
  }
  */

  
  DEFINE_ART_MODULE(CCKaonProducer)     
  
}

#endif

