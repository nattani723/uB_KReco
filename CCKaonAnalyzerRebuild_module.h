#ifndef CCKAON_ANALYSIS
#define CCKAON_ANALYSIS

#include "TH1.h"
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
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <string.h>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 

#include "larcoreobj/SummaryData/POTSummary.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include  "nusimdata/SimulationBase/GTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "/uboone/app/users/taniuchi/51_pandora/srcs/larpandoracontent/larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

//#include "larpandora/LArPandoraEventBuilding/LArPandoraTrackCreation_module.cc"

#include "larcore/Geometry/Geometry.h"

#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Helper function for PID stuff
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TRandom3.h"

#include "Pandora/PdgTable.h"
#include <chrono>

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <sys/stat.h>

#include "bad_channel_matching.h"
#include "sssVeto_BDT.class.h"
#include "DBSCAN.h"
//------------------------------------------------------------------------------------------------------------------------------------------


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

const int kMaxTracks=20;
const int kMaxShowers=20;
//const int kMaxTracks=50;
//const int kMaxShowers=50;
const int kMaxParticles=10;
const int kMaxMerge=10;

namespace Kaon_Analyzer
{
  //namespace microboone{

  using namespace std;
  
   class CCKaonAnalyzerRebuild : public art::EDAnalyzer { 
   public: 

     explicit CCKaonAnalyzerRebuild(fhicl::ParameterSet const& pset); 
     virtual ~CCKaonAnalyzerRebuild(); 

     void endSubRun(const art::SubRun &subrun); 
     void beginJob(); 
     void filter(const art::Event &evt); 
     void analyze(const art::Event& evt); 
     void reset(); 

     //void fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, const recob::Track trk, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart, int track_i=-1, int daughter_i=-1); 
     void fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, int track_i=-1, int daughter_i=-1); 
     //double ModBoxCorrection(const double dQdx, const float x, const float y, const float z); 
     //float GetLocalEFieldMag(const float x, const float y, const float z); 
     void fillPID(const std::vector<art::Ptr<anab::ParticleID>> &trackPID, double angle_y, int track_i=-1, int daughter_i=-1); 


     /*
     int fillAngularDistributionMapCheat(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
				          TVector3 Kend_candidate,
				          art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                          art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
                                          std::map<int, std::map<int, double>>& angular_distribution_mrgid_map,
					  std::map<int,int> &mrgidpdg,
					  int currentMergedId);
     */

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
			    std::vector<std::vector<art::Ptr<recob::Hit>>>& shower_spine_hit_list_vector,
			    TVector3 Kend_candidate,
			    std::map<double, TVector2, std::greater<>>& view_peak_map,
			    vector<TVector2> &best_peak_bins,
			    int& ini_hits,
			    int& num_hits,
			    float& closest_distance);

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


     /*     
     void trackRebuid(std::vector<art::Ptr<recob::Hit>>& hit_list,
		      art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
		      recob::Track& track);

     recob::Track buildTrack(int id,
			     lar_content::LArTrackStateVector& trackStateVector);
     */

   //void GetPFParticleIdMap(const Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleHandle &pfParticleHandle, Kaon_Analyzer::CCKaonAnalyzerRebuild::PFParticleIdMap &pfParticleMap); 
  
     double length(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi); 

     bool isInsideVolume(string volume, double x, double y, double z); 
     bool isInsideVolume(string volume, const TVector3& v) { 
       return isInsideVolume(volume, v.X(), v.Y(), v.Z()); 
      } 


    typedef art::ValidHandle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
    typedef std::vector< art::Ptr<recob::Track> > TrackVector;
    typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;
    typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;

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

     Int_t n_prip;
     std::vector<int> prip{vector<int>(20,-999)}; 
     Int_t n_prip_k_dau;
     std::vector<int> prip_k_dau{vector<int>(20,-999)};
     Bool_t IsKaon;
     Bool_t IsSingleKaon;
     Bool_t IsAssociatedKaon;
     Bool_t IsMuBR; 
     Bool_t IsPiBR; 
     Bool_t flag_prim_k;

     Bool_t flg_reco_k  = false;
     Bool_t flg_reco_mu  = false;
     Bool_t flg_reco_pi  = false;

     std::string Process;
     std::string EndProcess;

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
    Int_t   true_nprimary;
    Int_t   true_nhyperons;

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

    Float_t true_dau_muon_length;
    Float_t true_dau_muon_p;
    Float_t true_dau_muon_ke;
    Float_t true_dau_muon_theta;
    Float_t true_dau_muon_costheta;
    Float_t true_dau_muon_phi;
    Float_t true_dau_muon_ccmuon_angle;
    Float_t true_dau_muon_ccmuon_cosangle;
    Int_t   true_dau_muon_end_process;
    Float_t true_dau_muon_end_ke;
    Float_t true_dau_muon_start_x;
    Float_t true_dau_muon_start_y;
    Float_t true_dau_muon_start_z;
    Float_t true_dau_muon_end_x;
    Float_t true_dau_muon_end_y;
    Float_t true_dau_muon_end_z;
    Bool_t  true_dau_muon_end_inTPC;
    Bool_t  true_dau_muon_end_in5cmTPC;
    Bool_t  true_dau_muon_end_inCCInclusiveTPC;

    Float_t true_dau_pip_length;
    Float_t true_dau_pip_p;
    Float_t true_dau_pip_ke;
    Float_t true_dau_pip_theta;
    Float_t true_dau_pip_costheta;
    Float_t true_dau_pip_phi;
    Float_t true_dau_pip_ccmuon_angle;
    Float_t true_dau_pip_ccmuon_cosangle;
    Int_t   true_dau_pip_end_process;
    Float_t true_dau_pip_end_ke;
    Float_t true_dau_pip_start_x;
    Float_t true_dau_pip_start_y;
    Float_t true_dau_pip_start_z;
    Float_t true_dau_pip_end_x;
    Float_t true_dau_pip_end_y;
    Float_t true_dau_pip_end_z;
    Bool_t  true_dau_pip_end_inTPC;
    Bool_t  true_dau_pip_end_in5cmTPC;
    Bool_t  true_dau_pip_end_inCCInclusiveTPC;

    Float_t true_dau_pin_length;
    Float_t true_dau_pin_p;
    Float_t true_dau_pin_ke;
    Float_t true_dau_pin_theta;
    Float_t true_dau_pin_costheta;
    Float_t true_dau_pin_phi;
    Float_t true_dau_pin_ccmuon_angle;
    Float_t true_dau_pin_ccmuon_cosangle;
    Int_t   true_dau_pin_end_process;
    Float_t true_dau_pin_end_ke;
    Float_t true_dau_pin_start_x;
    Float_t true_dau_pin_start_y;
    Float_t true_dau_pin_start_z;
    Float_t true_dau_pin_end_x;
    Float_t true_dau_pin_end_y;
    Float_t true_dau_pin_end_z;
    Bool_t  true_dau_pin_end_inTPC;
    Bool_t  true_dau_pin_end_in5cmTPC;
    Bool_t  true_dau_pin_end_inCCInclusiveTPC;

    //should be [kMaxParticles] -> [kMaxTracks][kMaxParticles]
    Int_t cheat_num_hits[kMaxTracks][kMaxParticles];
    Int_t cheat_ini_hits[kMaxTracks][kMaxParticles];
    Float_t cheat_closest_distance[kMaxTracks][kMaxParticles];
    Float_t cheat_peak_pdg[kMaxTracks][kMaxParticles];
    Float_t cheat_peak_theta[kMaxTracks][kMaxParticles];
    Float_t cheat_peak_phi[kMaxTracks][kMaxParticles];
    Float_t best_peak_theta[kMaxTracks][kMaxParticles];
    Float_t best_peak_phi[kMaxTracks][kMaxParticles];
    Float_t best_peak_trkln[kMaxTracks][kMaxParticles];
    Float_t cheat_pip_trkln;
    Float_t cheat_mup_trkln;

    Float_t true_length[kMaxParticles];
    Float_t true_p[kMaxParticles];
    Float_t true_ke[kMaxParticles];
    Float_t true_theta[kMaxParticles];
    Float_t true_costheta[kMaxParticles];
    Float_t true_phi[kMaxParticles];
    Float_t true_ccmuon_angle[kMaxParticles];
    Float_t true_ccmuon_cosangle[kMaxParticles];
    Int_t   true_end_process[kMaxParticles];
    Float_t true_end_ke[kMaxParticles];
    Float_t true_end_x[kMaxParticles];
    Float_t true_end_y[kMaxParticles];
    Float_t true_end_z[kMaxParticles];
    Bool_t  true_end_inTPC[kMaxParticles];
    Bool_t  true_end_in5cmTPC[kMaxParticles];
    Bool_t  true_end_inCCInclusiveTPC[kMaxParticles];

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


    Int_t   true_nkaons_anti;
    Float_t true_kaon_length_anti;
    Float_t true_kaon_p_anti;
    Float_t true_kaon_ke_anti;
    Float_t true_kaon_theta_anti;
    Float_t true_kaon_costheta_anti;
    Float_t true_kaon_phi_anti;
    Float_t true_kaon_ccmuon_angle_anti;
    Float_t true_kaon_ccmuon_cosangle_anti;
    Int_t   true_kaon_end_process_anti;
    Float_t true_kaon_end_ke_anti;
    Float_t true_kaon_end_x_anti;
    Float_t true_kaon_end_y_anti;
    Float_t true_kaon_end_z_anti;
    Bool_t  true_kaon_end_inTPC_anti;
    Bool_t  true_kaon_end_in5cmTPC_anti;
    Bool_t  true_kaon_end_inCCInclusiveTPC_anti;

    Int_t   true_kaon_ndaughters_anti;
    Int_t   true_kaon_ndaughters_decay_anti;
    Int_t   true_kaon_ndaughters_inelastic_anti;
    Int_t   true_kaon_ndecmup_anti;
    Int_t   true_kaon_ndecpip_anti;
    Int_t   true_kaon_ninekap_anti;
    Int_t   true_kaon_ninepip_anti;
    Int_t   true_kaon_ninepro_anti;
    Float_t true_kaon_daughter_length_anti;
    Float_t true_kaon_daughter_p_anti;
    Float_t true_kaon_daughter_ke_anti;
    Float_t true_kaon_daughter_theta_anti;
    Float_t true_kaon_daughter_costheta_anti;
    Float_t true_kaon_daughter_angle_anti;
    Float_t true_kaon_daughter_cosangle_anti;
    Int_t   true_kaon_daughter_pdg_anti;
    Float_t true_kaon_daughter_end_x_anti;
    Float_t true_kaon_daughter_end_y_anti;
    Float_t true_kaon_daughter_end_z_anti;
    Bool_t  true_kaon_daughter_end_inTPC_anti;
    Bool_t  true_kaon_daughter_end_in5cmTPC_anti;
    Bool_t  true_kaon_daughter_end_inCCInclusiveTPC_anti;

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
     Int_t   reco_nshowers;
     Int_t   nTracks;
     Int_t   nShowers; 
    //    Int_t   m_reco_track_sliceId[kMaxTracks];
    //    Int_t   m_reco_track_is_nuslice[kMaxTracks];

     Float_t true_kaon_start_x;
     Float_t true_kaon_start_y;
     Float_t true_kaon_start_z;

     Float_t true_kaon_daughter_start_x;
     Float_t true_kaon_daughter_start_y;
     Float_t true_kaon_daughter_start_z;
     
     Float_t reco_track_start_x[kMaxTracks];//
     Float_t reco_track_start_y[kMaxTracks];
    Float_t reco_track_start_z[kMaxTracks];
    Float_t reco_track_end_x[kMaxTracks];
    Float_t reco_track_end_y[kMaxTracks];
    Float_t reco_track_end_z[kMaxTracks];


    Float_t reco_track_daughter_start_x[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_start_y[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_start_z[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_end_x[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_end_y[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_end_z[kMaxTracks][kMaxTracks];


    Float_t reco_track_distance[kMaxTracks];
    Int_t   reco_track_nhits0[kMaxTracks];
    Int_t   reco_track_nhits1[kMaxTracks];
    Int_t   reco_track_nhits2[kMaxTracks];

    Float_t   reco_track_kin0[kMaxTracks];
    Float_t   reco_track_kin1[kMaxTracks];
    Float_t   reco_track_kin2[kMaxTracks];
    //vector<vector<Float_t>> reco_track_dEdx;
    //vector<vector<Float_t>> reco_track_ResRan;

    Float_t reco_track_dEdx_pl0[kMaxTracks][2000];
    Float_t reco_track_ResRan_pl0[kMaxTracks][2000];

    Float_t reco_track_dEdx_pl1[kMaxTracks][2000];
    Float_t reco_track_ResRan_pl1[kMaxTracks][2000];

    Float_t reco_track_dEdx_pl2[kMaxTracks][2000];
    Float_t reco_track_ResRan_pl2[kMaxTracks][2000];
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

    Float_t reco_track_match_e[kMaxTracks][kMaxMerge];
    Float_t reco_track_match_hit[kMaxTracks][kMaxMerge];
    Int_t reco_track_match_epdg[kMaxTracks][kMaxMerge];
    Int_t reco_track_match_hitpdg[kMaxTracks][kMaxMerge];

    Float_t reco_track_daughter_match_e[kMaxTracks][kMaxTracks][kMaxMerge];
    Float_t reco_track_daughter_match_hit[kMaxTracks][kMaxTracks][kMaxMerge];
    Int_t reco_track_daughter_match_epdg[kMaxTracks][kMaxTracks][kMaxMerge];
    Int_t reco_track_daughter_match_hitpdg[kMaxTracks][kMaxTracks][kMaxMerge];

    Float_t reco_track_daughter_shower_match_e[kMaxTracks][kMaxShowers][kMaxMerge];
    Float_t reco_track_daughter_shower_match_hit[kMaxTracks][kMaxShowers][kMaxMerge];
    Int_t reco_track_daughter_shower_match_epdg[kMaxTracks][kMaxShowers][kMaxMerge];
    Int_t reco_track_daughter_shower_match_hitpdg[kMaxTracks][kMaxShowers][kMaxMerge];


    Int_t   reco_track_ndaughters[kMaxTracks];
    Float_t reco_track_daughter_distance[kMaxTracks][kMaxTracks];
    Float_t reco_track_daughter_vtx_distance[kMaxTracks][kMaxTracks];
    Float_t reco_angle_track_daughter[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits0[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits1[kMaxTracks][kMaxTracks];
    Int_t   reco_track_daughter_nhits2[kMaxTracks][kMaxTracks];


     Int_t   reco_shower_ndaughters[kMaxTracks]; 
     Float_t reco_track_daughter_distance_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_vtx_distance_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_angle_track_daughter_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_angle_daughter_track_daughter_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_length_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_theta_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_phi_sh[kMaxTracks][kMaxShowers]; 
     Int_t reco_track_daughter_true_pdg_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_open_angle_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_dedx_pl0_sh[kMaxTracks][kMaxShowers];
     Float_t reco_track_daughter_dedx_pl1_sh[kMaxTracks][kMaxShowers];
     Float_t reco_track_daughter_dedx_pl2_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_energy_pl0_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_energy_pl1_sh[kMaxTracks][kMaxShowers]; 
     Float_t reco_track_daughter_energy_pl2_sh[kMaxTracks][kMaxShowers]; 
 
     Bool_t  reco_track_daughter_vtx_inTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_vtx_in5cmTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_vtx_inCCInclusiveTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_end_inTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_end_in5cmTPC_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_end_inCCInclusiveTPC_sh[kMaxTracks][kMaxShowers]; 

     Int_t   reco_track_true_pdg_sh[kMaxShowers]; 
     Int_t   reco_track_true_origin_sh[kMaxShowers];
     Bool_t  reco_track_true_primary_sh[kMaxShowers]; 
     Bool_t  reco_track_true_end_inTPC_sh[kMaxShowers];
     Bool_t  reco_track_true_end_in5cmTPC_sh[kMaxShowers];
     Bool_t  reco_track_true_end_inCCInclusiveTPC_sh[kMaxShowers]; 
     Float_t reco_track_true_length_sh[kMaxShowers]; 

     Int_t   reco_track_daughter_true_origin_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_true_primary_sh[kMaxTracks][kMaxShowers]; 
     Bool_t  reco_track_daughter_true_end_inTPC_sh[kMaxTracks][kMaxShowers];
     Bool_t  reco_track_daughter_true_end_in5cmTPC_sh[kMaxTracks][kMaxShowers];
     Bool_t  reco_track_daughter_true_end_inCCInclusiveTPC_sh[kMaxTracks][kMaxShowers];
     Float_t reco_track_daughter_true_length_sh[kMaxTracks][kMaxShowers];
     Int_t   reco_track_daughter_true_mother_sh[kMaxTracks][kMaxShowers];

    //vector<vector<vector<Float_t>>> reco_track_daughter_dEdx;
    //vector<vector<vector<Float_t>>> reco_track_daughter_ResRan;

    Float_t reco_track_daughter_dEdx_pl0[kMaxTracks][kMaxTracks][2000];
    Float_t reco_track_daughter_ResRan_pl0[kMaxTracks][kMaxTracks][2000];

    Float_t reco_track_daughter_dEdx_pl1[kMaxTracks][kMaxTracks][2000];
    Float_t reco_track_daughter_ResRan_pl1[kMaxTracks][kMaxTracks][2000];

    Float_t reco_track_daughter_dEdx_pl2[kMaxTracks][kMaxTracks][2000];
    Float_t reco_track_daughter_ResRan_pl2[kMaxTracks][kMaxTracks][2000];


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
    std::string fShowerModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fPIDLabel;
    std::string fHitTruthAssns;
    std::string fHitTrackAssns;
    std::string fHitShowerAssns;
    std::string m_pfp_producer;
    std::string fPFParticleLabel;
    std::string fSpacePointproducer;
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

    //TString output_pdf;
    TCanvas * c = new TCanvas("c", "c", 800, 600); 
    //bool isFirstpdf = true;
    //bool isFirstpdf_cheat = true;
    //int i_pdf=0;
    //int i_pdf_cheat=0;
    //const int num_pdf=3;
    //const int num_pdf_cheat=3;

    //std::map<int, int> angular_distribution_map_track;
    //std::map<int, int> angular_distribution_map_shower;
    //std::map<int, std::map<int, int>> angular_distribution_mrgid_map;
    //vector<TH1D*> h_angular_distribution_pfparticle;
    //vector<bool> v_trk_flg;

    std::map<std::string, std::vector<float>> event_weight;

    std::vector<std::string> evtwgt_funcname;          // the name of the functions used
    std::vector<std::vector<float>> evtwgt_weight;    // the weights (a vector for each function used)
    std::vector<int> evtwgt_nweight;                   // number of weights for each function
    //Int_t evtwgt_nfunc;                                // number of functions used


            void ClearVertex();
            /**
             *  @brief  Produce a mapping from PFParticle ID to the art ptr to the PFParticle itself for fast navigation
             *
             *  @param  pfParticleHandle the handle for the PFParticle collection
             *  @param  pfParticleMap the mapping from ID to PFParticle
             */
	    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

            /**
             * @brief Print out scores in PFParticleMetadata
             *
             * @param evt the art event to analyze
             * @param pfParticleHandle the handle for the PFParticle collection
             */
            void PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const;

            /**
             *  @brief  Produce a mapping from PFParticle ID to the art ptr to the PFParticle itself for fast navigation
             *
             *  @param  pfParticleMap the mapping from ID to PFParticle
             *  @param  crParticles a vector to hold the top-level PFParticles reconstructed under the cosmic hypothesis
             *  @param  nuParticles a vector to hold the final-states of the reconstruced neutrino
             */
            void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap,const lar_pandora::PFParticlesToVertices &particlesToVertices, PFParticleVector &crParticles, PFParticleVector &nuParticles);

            /**
             *  @brief  Collect associated tracks and showers to particles in an input particle vector
             *
             *  @param  particles a vector holding PFParticles from which to find the associated tracks and showers
             *  @param  pfParticleHandle the handle for the PFParticle collection
             *  @param  evt the art event to analyze
             *  @param  tracks a vector to hold the associated tracks
             *  @param  showers a vector to hold the associated showers
             */
            void CollectTracksAndShowers(const PFParticleVector &particles,const PFParticleIdMap pfParticleMap,  const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap);

            void FillTracksAndShowers( const std::vector< art::Ptr<recob::Track> > & associatedTracks, const std::vector< art::Ptr<recob::Shower> > & associatedShowers, const art::Ptr<recob::PFParticle> &pParticle , const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap);

            void GetVertex(const lar_pandora::PFParticlesToVertices & particlestoVertices, const art::Ptr<recob::PFParticle> & particle );

            void CollectCalo(const art::Event &evt,const art::Ptr<recob::Shower> &shower);


            /*
             *@brief Calculated the shower energy by looping over all the hits and summing the charge
             *@param hits -  an art pointer of all the hits in a shower
             *
             *
             * */
            double CalcEShower(const std::vector<art::Ptr<recob::Hit>> &hits);

            double CalcEShowerPlane(const std::vector<art::Ptr<recob::Hit>>& hits, int plane);

            int getNHitsPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane);


            /**
             *@brief Takes a hit and multiplies the charge by the gain
             *@param thishitptr art pointer to a hit
             *@param plane the plane the hit is on
             **/
            double GetQHit(art::Ptr<recob::Hit> thishitptr, int plane);

            /**
             * @brief Calculate the E value in MeV for a given hit
             * @param thishit - an individual hit 
             * 
             *
             * */
            double QtoEConversionHit(art::Ptr<recob::Hit> thishitptr, int plane);

            /**
             * @brief Calculate the E value in MeV from a given Q value
             * @param q - the charge value
             * 
             * */
            double QtoEConversion(double q);


            /**
             *@brief Takes a vector of dQ/dx values and converts to dE/dx
             *@param dqdx - vector of dqdx points
             *
             * */
            std::vector<double> CalcdEdxFromdQdx(std::vector<double> dqdx);

            /**
             *
             *@brief For a single shower, calculates the dQdx for each hit in the clusters in the shower for a single plane
             *@param shower - a Pandora shower
             *@param clusters - all of the clusters in the shower
             *@param clusterToHitMap - a map between each cluster and all of the hits in the cluster
             *@param plane - a single plane
             * * */

            std::vector<double> CalcdQdxShower(
                    const art::Ptr<recob::Shower>& shower,
                    const std::vector<art::Ptr<recob::Cluster>> & clusters, 
                    std::map<art::Ptr<recob::Cluster>,    std::vector<art::Ptr<recob::Hit>> > &  clusterToHitMap ,int plane);
            /**
             *@brief Gets the pitch between the 3D reconstructed shower direction and the wires for a given plane (the dx in dQdx)
             *@param shower_dir - the 3D shower direction
             *@param plane - a single plane
             * */
            double getPitch(TVector3 shower_dir, int plane);
            TVector3 getWireVec(int plane);
            double getCoswrtWires(TVector3 shower_dir, TVector3 wire_dir);
            double getAnglewrtWires(TVector3 shower_dir, int plane);

            double getAmalgamateddEdx(double angle_wrt_plane0, double angle_wrt_plane1, double angle_wrt_plane2, double median_plane0, double median_plane1, double median_plane2, int plane0_nhits, int plane1_nhits, int plane2_nhits);
            int getAmalgamateddEdxNHits(double amalgamateddEdx, double median_plane0, double median_plane1, double median_plane2, int plane0_nhits, int plane1_nhits, int plane2_nhits);
            double degToRad(double deg);
            double radToDeg(double rad);
            /**
             *@brief Calculates the four corners of a box of given length and width around a cluster given the start point and axis direction
             *@param cluster_start - the start position of a cluster in CM
             *@param cluster_axis - calculated from the cluster end minus the cluster start
             *@param width - typically ~1cm
             *@param length - typically a few cm
             *
             * */
            std::vector<std::vector<double>> buildRectangle(std::vector<double> cluster_start, std::vector<double> cluster_axis, double width, double length);

            /**
             *@brief For a 2d point on a plane in cm and a rectangle, returns true if the point is inside of the rectangle
             *@param thishit_pos - 2d location of a hit in cm
             *@param rectangle - vector of the positions of the four corners of the rectangle
             *
             * */
            bool insideBox(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle);

            /**
             *
             *@brief For a 2d point on a plane in cm and a rectangle, returns true if ponint is inside or on the boundary
             *uses triangle area check
             *
             * */
            bool isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle);

            double areaTriangle(double x1, double y1, double x2, double y2, double x3, double y3);

            /***
             *@brief returns the value at the median position in a vector of doubles, returns nan for vector of size <= 0
             *@param thisvector - vector of doubles
             *
             * */
            double getMedian(std::vector<double> thisvector);


            //----------------  Templatees ----------------------------
            void AnalyzeTemplates();
            void ClearTemplates();
            void ResizeTemplates(size_t);
            void CreateTemplateBranches();


            //---------------- SecondShower----
            void ClearSecondShowers();
            void ResizeSecondShowers(size_t size);

            void CreateSecondShowerBranches();

            void SecondShowerSearch(
                    const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap,
                    const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,
                    const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,  
                    const std::map<art::Ptr<recob::PFParticle>, int> & pfParticleToSliceIDMap, const std::map<int, std::vector<art::Ptr<recob::Hit>>>& sliceIDToHitsMap,
                    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
                    std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
                    std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
                    std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap);

            std::vector<double>SecondShowerMatching(std::vector<art::Ptr<recob::Hit>>& hitz,
                    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
                    std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
                    std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
                    std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap);



	    //            sss_score ScoreCluster(int,int,std::vector<art::Ptr<recob::Hit>>&,double,double, const art::Ptr<recob::Shower>&);
            TGraph* GetNearestNpts(int,int,std::vector<art::Ptr<recob::Hit>>&,double,double,int);
            int CompareToShowers(int,int,std::vector<art::Ptr<recob::Hit>>&,double,double,
                    const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & showertopfparticlemap,      const   std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfparticletohitsmap,                    double eps);
            //---------------- Isolation ----------------- 

            void ClearIsolation();
            void CreateIsolationBranches();
            void IsolationStudy(
                    const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap,
                    const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,
                    const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,  
                    const std::map<art::Ptr<recob::PFParticle>, int> & pfParticleToSliceIDMap, const std::map<int, std::vector<art::Ptr<recob::Hit>>>& sliceIDToHitsMap);



            //----------------  Flashes ----------------------------
                void AnalyzeFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes,  art::Handle<std::vector<crt::CRTHit>> crthit_h, double evt_timeGPS_nsec);
                     
          //  void AnalyzeFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes,  art::Handle<std::vector<crt::CRTHit>> crthit_h);
            void ClearFlashes();
            void ResizeFlashes(size_t);
            void CreateFlashBranches();


             //----------------  Tracks ---------------------------- 
             void AnalyzeTracks(const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & tracktopfparticlemap, std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>>> & pfparticletospacepointmap , std::map<int, art::Ptr<simb::MCParticle> > &  MCParticleToTrackIdMap, std::map<int, double> &sliceIdToNuScoreMap, 
                     std::map<art::Ptr<recob::PFParticle>,bool> &PFPToClearCosmicMap, 
                     std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap, 
                     std::map<art::Ptr<recob::PFParticle>,double> &PFPToTrackScoreMap, 
                     std::map<art::Ptr<recob::PFParticle>,bool> &PFPToNuSliceMap, 
                     PFParticleIdMap &pfParticleMap 
                     ); 

             void ClearTracks(); 
             void ResizeTracks(size_t); 
             void CreateTrackBranches(); 
             void AnalyzeTrackCalo(const std::vector<art::Ptr<recob::Track>> &tracks, std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<anab::Calorimetry>>> &trackToCaloMap); 
            void RecoMCTracks(const std::vector<art::Ptr<recob::Track>>& tracks,  std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle>> & trackToPFParticleMap, std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > & trackToMCParticleMap,  std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap,std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector, std::map< int, art::Ptr<simb::MCParticle> > &      MCParticleToTrackIdMap, 
                    std::map<int, double>& sliceIdToNuScoreMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
                    std::vector<double>& vec);

            void CollectPID(std::vector<art::Ptr<recob::Track>> & tracks,std::map< art::Ptr<recob::Track>, art::Ptr<anab::ParticleID>> & trackToPIDMap);
            TGraph proton_length2energy_tgraph;

            //----------------  Showers ----------------------------

            void AnalyzeShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> & showerToPFParticleMap, std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap,std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > & pfParticleToClusterMap, std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > & clusterToHitMap,
                    std::map<int, double> &sliceIdToNuScoreMap,
                    std::map<art::Ptr<recob::PFParticle>,bool> &PFPToClearCosmicMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap, 
                    std::map<art::Ptr<recob::PFParticle>,bool> &PFPToNuSliceMap,
                    std::map<art::Ptr<recob::PFParticle>,double> &PFPToTrackScoreMap,
                    PFParticleIdMap &pfParticleMap,
                    std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>> &PFPtoShowerReco3DMap
                    );
            void ClearShowers();
            void ResizeShowers(size_t);
            void CreateShowerBranches();
            void AnalyzeKalmanShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> & showerToPFParticleMap,
                    std::map<art::Ptr<recob::PFParticle>,art::Ptr<recob::Track>> & pfptotrkmap,
                    std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<anab::Calorimetry>>> & trktocalomap,
                    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap
                    );

            void RecoMCShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> & showerToPFParticleMap, std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,  std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap,
                    std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector);

            std::vector<double> showerRecoMCmatching(std::vector<art::Ptr<recob::Shower>>& objectVector,
                    std::map<art::Ptr<recob::Shower>,art::Ptr<simb::MCParticle>>& objectToMCParticleMap,
                    std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>>& objectToPFParticleMap,
                    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> >& pfParticleToHitsMap,
                    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
                    std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
                    std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
                    std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap,
                    std::map<int, double> & sliceIdToNuScoreMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToNuSliceMap);

  


            //---------------- MCTruths ----------------------------

            void AnalyzeMCTruths(std::vector<art::Ptr<simb::MCTruth>> & mcTruthVector,  std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector );
            void ClearMCTruths();
            void ResizeMCTruths(size_t);
            void CreateMCTruthBranches();

            std::map<int,std::string> is_delta_map;

            //---------------- EventWeight ----------------------------

            void AnalyzeEventWeight(art::Event const & e );
            void ClearEventWeightBranches();
            void CreateEventWeightBranches();

            //These three are shameless steals from LArPandorHelper But overlays dont work so this is a direct clone. We will filter out later.




            void CollectSimChannels(const art::Event &evt, const std::string &label,  std::vector< art::Ptr<sim::SimChannel> >  &simChannelVector);
            void CollectMCParticles(const art::Event &evt, const std::string &label, std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> &truthToParticles,        std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>              &particlesToTruth, std::map< int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIdMap);
            void BuildMCParticleHitMaps(const art::Event &evt, const std::string &label, const std::vector<art::Ptr<recob::Hit>> &hitVector,   std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  &particlesToHits,         std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  &hitsToParticles, const lar_pandora::LArPandoraHelper::DaughterMode daughterMode, std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap);


            //-------------- Slices/Pandora Metadata ---------------//
            void  ClearSlices();
            void  ResizeSlices(size_t size); 
            //void  ResizeMatchedSlices(size_t size_shower ,size_t size_track); 
            void CreateSliceBranches();
            //void CreateMatchedSliceBranches();
	    
            void AnalyzeSlices(std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > & pfParticleToMetadataMap,
                    PFParticleIdMap &pfParticleMap,
                    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & primaryPFPSliceIdVec,   
                    std::map<int, double> & sliceIdToNuScoreMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,   
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToNuSliceMap,
                    std::map<art::Ptr<recob::PFParticle>,double>& PFPToTrackScoreMap);
	    

            // void GetPFPsPerSlice( std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
            //        std::map<int, int>& sliceIdToNumPFPsMap );
            std::vector<int>  GetPFPsPerSlice( std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap );
            //void  GetPFPsPerSlice( std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap , std::vector<int> &sliceIdToNumPFPsvec);

            int GetShowerSlice(art::Ptr<recob::Shower>& this_shower, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>>& showerToPFParticleMap, std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec);

            int GetTrackSlice(art::Ptr<recob::Track>& this_track, std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>& trackToPFParticleMap, std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec);
            //can also look at things like shower energy, conversion length, etc.

            std::vector<int> GetNumShowersPerSlice(std::map< art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>>& showerToPFParticleMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap );
            //        void GetNumShowersPerSlice(std::map< art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>>& showerToPFParticleMap,
            //                 std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
            //                std::map<int, int>& sliceIdToNumShowersMap );


            std::vector<int> GetNumTracksPerSlice(std::map< art::Ptr<recob::Track>,art::Ptr<recob::PFParticle>>& trackToPFParticleMap,   
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap);

            void AnalyzeRecoMCSlices(std::string signal_def, std::map<int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIDMap,
                    std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle> > & showerToPFParticleMap, 
                    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec, 
                    std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,
                    std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle> > & trackToNuPFParticleMap,
                    std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > &trackToMCParticleMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap);


            void FindSignalSlice(std::string signal_def, std::map<int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIDMap,std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle> > & showerToPFParticleMap,  std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec, std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap, std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle> > & trackToNuPFParticleMap, std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > &trackToMCParticleMap);

            int  m_reco_slice_num; //total number of slices in the event
            std::vector<double> m_reco_slice_nuscore; //vector of the neutrino score for each slice in an event
            int m_reco_slice_shower_num_matched_signal; //the number of sim showers matched an MCP in the signal def
            int m_reco_slice_track_num_matched_signal; //the number of sim showers matched an MCP in the signal def
            std::vector<int> m_reco_slice_shower_matched_sliceId; //the slice id for each matched shower
            std::vector<int> m_reco_slice_track_matched_sliceId; //the slice id for each matched track

            std::vector<int> m_reco_slice_num_pfps; //the total number of PFP's per slice
            std::vector<int> m_reco_slice_num_showers; //the subset of PFP's that are showers
            std::vector<int> m_reco_slice_num_tracks; //the subset of PFP's that are tracks

            std::vector<double> m_reco_slice_shower_matched_energy; //the energy for each matched shower
            std::vector<double> m_reco_slice_track_matched_energy; //the energy for each matched track
            std::vector<double> m_reco_slice_shower_matched_conversion; //the conversion distance for each matched shower
            std::vector<double> m_reco_slice_shower_matched_overlay_frac; //fraction of overlay hits for each matched shower
            //std::map<art::Ptr<recob::PFParticle>, double > & pfParticleToNuScoreMap;//is filled during analyze slices

            std::vector<double> m_matched_signal_shower_overlay_fraction;
            //std::vector<double> m_matched_signal_shower_conversion_length;
            std::vector<double> m_matched_signal_shower_true_E;
            std::vector<double> m_matched_signal_shower_nuscore;
            std::vector<int> m_matched_signal_shower_sliceId;
            std::vector<bool> m_matched_signal_shower_is_clearcosmic;
            int m_matched_signal_shower_num = 0;
            std::vector<bool> m_matched_signal_shower_is_nuslice;
            std::vector<int> m_matched_signal_shower_tracks_in_slice;
            std::vector<int> m_matched_signal_shower_showers_in_slice;


            std::vector<double> m_matched_signal_track_true_E;
            std::vector<double> m_matched_signal_track_nuscore;
            std::vector<int> m_matched_signal_track_sliceId;
            std::vector<bool> m_matched_signal_track_is_clearcosmic;
            //  std::vector<bool> m_matched_signal_track_is_nuslice;
            std::vector<bool> m_matched_signal_track_is_nuslice;
            std::vector<int> m_matched_signal_track_tracks_in_slice;
            std::vector<int> m_matched_signal_track_showers_in_slice;

            int m_matched_signal_track_num = 0;   

            //int m_matched_signal_total_num_slices;

            bool m_reco_1g1p_is_same_slice;
            bool m_reco_1g1p_is_multiple_slices;
            bool m_reco_1g1p_is_nuslice;
            bool m_reco_1g0p_is_nuslice;
            double m_reco_1g1p_nuscore;
            double  m_reco_1g0p_nuscore;
            bool m_is_matched_1g1p;
            bool m_is_matched_1g0p;
            bool m_no_matched_showers;
            bool m_multiple_matched_showers;
            bool m_multiple_matched_tracks;



            //------------------ Delaunay triangle tools -----------//

            double triangle_area(double a1, double a2, double b1, double b2, double c1, double c2);
            int quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area);
            int delaunay_hit_wrapper(const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<int> & num_hits, std::vector<int>& num_triangles, std::vector<double> & area);


            int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected);
            int spacecharge_correction(const simb::MCParticle & mcparticle, std::vector<double> & corrected);

            int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected, std::vector<double> & input);

            //databased http://dbdata0vm.fnal.gov:8186/uboonecon_prod/app/data?f=channelstatus_data&t=357812824
            std::vector<std::pair<int,int>> bad_channel_list_fixed_mcc9;
            std::map<int,bool> bad_channel_map_fixed_mcc9;

            TRandom3 *rangen;
            std::string m_shower3dLabel;
            std::string m_showerKalmanLabel;
            std::string m_showerKalmanCaloLabel;
	    std::string m_pandoraLabel;         ///< The label for the pandora producer
            std::string m_trackLabel;           ///< The label for the track producer from PFParticles
	    //            std::string m_trackLabel_old;
            std::string m_sliceLabel;          
            std::string m_showerLabel;          ///< The label for the shower producer from PFParticles
            std::string m_caloLabel;            ///< The label for calorimetry associations producer
            std::string m_flashLabel;
            std::string m_geantModuleLabel;
            std::string m_backtrackerLabel;
            std::string m_hitfinderLabel;
            std::string m_hitMCParticleAssnsLabel;
            std::string m_potLabel;
            std::string m_generatorLabel;
            std::string m_badChannelLabel;
            std::string m_badChannelProducer;
            std::string m_mcTrackLabel;
            std::string m_mcShowerLabel;
            std::string m_pidLabel;            ///< For PID stuff
            std::string m_CRTTzeroLabel;
            std::string m_CRTHitProducer;
            bool m_use_PID_algorithms;
            bool m_use_delaunay;
            bool m_is_verbose;
	    bool m_print_out_event;
            bool m_is_data;
            bool m_is_overlayed;
            bool m_run_all_pfps;
            bool m_has_CRT;
            bool m_fill_trees;
            bool m_run_pi0_filter;

            bool m_runCRT;
            double m_DTOffset;
            double  m_Resolution;
            std::string  m_DAQHeaderProducer;
            std::ofstream out_stream;

            double m_exiting_photon_energy_threshold ;
            double m_exiting_proton_energy_threshold ;

            double m_track_calo_min_dEdx;
            double m_track_calo_max_dEdx;
            double m_track_calo_min_dEdx_hits;
            double m_track_calo_trunc_fraction;

            detinfo::DetectorProperties const * theDetector ;// = lar::providerFrom<detinfo::DetectorPropertiesService>();
            detinfo::DetectorClocks    const *  detClocks   ;//= lar::providerFrom<detinfo::DetectorClocksService>();
            spacecharge::SpaceCharge const * SCE;
            geo::GeometryCore const * geom;
            double m_work_function;
            double m_recombination_factor;
            //double m_gain;
            std::vector<double> m_gain_mc; 
            std::vector<double> m_gain_data; 
            double m_wire_spacing;

            int m_Cryostat;
            int m_TPC;

            double m_width_dqdx_box;
            double m_length_dqdx_box;

            TTree* pot_tree;
            TTree* vertex_tree;
            TTree* eventweight_tree;
            TTree* ncdelta_slice_tree;

            //------------ POT related variables --------------
            int m_number_of_events;
            double m_pot_count;
            int m_number_of_vertices;

            //------------ Event Related Variables -------------
            int m_run_number;
            int m_subrun_number;
            int m_event_number;

            int m_test_matched_hits;
            int m_reco_slice_objects;
            //------- Second shower related variables ----
            int m_sss_num_unassociated_hits;
            int m_sss_num_associated_hits;

            int m_sss_num_candidates;

            ReadBDT * sssVetov1;

            std::vector<int> m_sss_candidate_num_hits;
            std::vector<int> m_sss_candidate_num_wires;
            std::vector<int>  m_sss_candidate_num_ticks;
            std::vector<int>  m_sss_candidate_plane;
            std::vector<double> m_sss_candidate_PCA;
            std::vector<double> m_sss_candidate_impact_parameter;
            std::vector<double> m_sss_candidate_fit_slope;
            std::vector<double> m_sss_candidate_veto_score;
            std::vector<double> m_sss_candidate_fit_constant;
            std::vector<double>  m_sss_candidate_mean_tick;
            std::vector<double> m_sss_candidate_max_tick;
            std::vector<double> m_sss_candidate_min_tick;
            std::vector<double> m_sss_candidate_min_wire;
            std::vector<double> m_sss_candidate_max_wire;
            std::vector<double> m_sss_candidate_mean_wire;
            std::vector<double> m_sss_candidate_min_dist;
            std::vector<double> m_sss_candidate_energy;
            std::vector<double> m_sss_candidate_angle_to_shower;
            std::vector<double> m_sss_candidate_closest_neighbour;
            std::vector<int>   m_sss_candidate_matched;
            std::vector<int>       m_sss_candidate_pdg;
            std::vector<int>        m_sss_candidate_parent_pdg;
            std::vector<int>    m_sss_candidate_trackid;
            std::vector<double>       m_sss_candidate_overlay_fraction;

            bool bool_make_sss_plots;

            //------------ Vertex Related variables -------------
            int m_reco_vertex_size;
            double m_vertex_pos_x;
            double m_vertex_pos_y;
            double m_vertex_pos_z;
            double m_vertex_pos_tick;
            double m_vertex_pos_wire_p0;
            double m_vertex_pos_wire_p2;
            double m_vertex_pos_wire_p1;

            int m_reco_asso_showers;

            double m_reco_vertex_to_nearest_dead_wire_plane0;
            double m_reco_vertex_to_nearest_dead_wire_plane1;
            double m_reco_vertex_to_nearest_dead_wire_plane2;

            //added eventweight
            //-------------- EventWeight related variables -------------
            static const int k_max_mc_particles=100;

            int m_run_number_eventweight;
            int m_subrun_number_eventweight;
            int m_event_number_eventweight;

            double m_mcflux_nu_pos_x;
            double m_mcflux_nu_pos_y;
            double m_mcflux_nu_pos_z;
            double m_mcflux_nu_mom_x;
            double m_mcflux_nu_mom_y;
            double m_mcflux_nu_mom_z;
            double m_mcflux_nu_mom_E;
            int m_mcflux_ntype;
            int m_mcflux_ptype;
            double m_mcflux_nimpwt;
            double m_mcflux_dk2gen;
            double m_mcflux_nenergyn;
            double m_mcflux_tpx;
            double m_mcflux_tpy;
            double m_mcflux_tpz;
            double m_mcflux_vx;
            double m_mcflux_vy;
            double m_mcflux_vz;
            int m_mcflux_tptype;
            int m_mctruth_nparticles;
            int m_mctruth_particles_track_Id[k_max_mc_particles];
            int m_mctruth_particles_pdg_code[k_max_mc_particles];
            int m_mctruth_particles_mother[k_max_mc_particles];
            int m_mctruth_particles_status_code[k_max_mc_particles];
            int m_mctruth_particles_num_daughters[k_max_mc_particles]; //other similar variables
            int m_mctruth_particles_daughters[100][100];
            double m_mctruth_particles_Gvx[k_max_mc_particles];
            double m_mctruth_particles_Gvy[k_max_mc_particles];
            double m_mctruth_particles_Gvz[k_max_mc_particles];
            double m_mctruth_particles_Gvt[k_max_mc_particles];
            double m_mctruth_particles_px0[k_max_mc_particles];
            double m_mctruth_particles_py0[k_max_mc_particles];
            double m_mctruth_particles_pz0[k_max_mc_particles];
            double m_mctruth_particles_e0[k_max_mc_particles];
            int m_mctruth_particles_rescatter[k_max_mc_particles];
            double m_mctruth_particles_polx[k_max_mc_particles];
            double m_mctruth_particles_poly[k_max_mc_particles];
            double m_mctruth_particles_polz[k_max_mc_particles];
            int m_mctruth_neutrino_ccnc;
            int m_mctruth_neutrino_mode;
            int m_mctruth_neutrino_interaction_type;
            int m_mctruth_neutrino_target;
            int m_mctruth_neutrino_nucleon;
            int m_mctruth_neutrino_quark;
            double m_mctruth_neutrino_w;
            double m_mctruth_neutrino_x;
            double m_mctruth_neutrino_y;
            double m_mctruth_neutrino_qsqr;
            bool m_gtruth_is_sea_quark;
            int m_gtruth_tgt_pdg;
            double m_gtruth_weight;
            double m_gtruth_probability;
            double m_gtruth_xsec;
            double m_gtruth_diff_xsec;
            int m_gtruth_gphase_space;
            double m_gtruth_vertex_x;
            double m_gtruth_vertex_y;
            double m_gtruth_vertex_z;
            double m_gtruth_vertex_T;
            int m_gtruth_gscatter;
            int m_gtruth_gint;
            int m_gtruth_res_num;
            int m_gtruth_num_piplus;
            int m_gtruth_num_pi0;
            int m_gtruth_num_piminus;
            int m_gtruth_num_proton;
            int m_gtruth_num_neutron;
            bool m_gtruth_is_charm;
            double m_gtruth_gx;
            double m_gtruth_gy;
            double m_gtruth_gt;
            double m_gtruth_gw;
            double m_gtruth_gQ2;
            double m_gtruth_gq2;
            int m_gtruth_probe_pdg;
            double m_gtruth_probe_p4_x;
            double m_gtruth_probe_p4_y;
            double m_gtruth_probe_p4_z;
            double m_gtruth_probe_p4_E;
            double m_gtruth_hit_nuc_p4_x;
            double m_gtruth_hit_nuc_p4_y;
            double m_gtruth_hit_nuc_p4_z;
            double m_gtruth_hit_nuc_p4_E;
            double m_gtruth_fs_had_syst_p4_x;
            double m_gtruth_fs_had_syst_p4_y;
            double m_gtruth_fs_had_syst_p4_z;
            double m_gtruth_fs_had_syst_p4_E;

            //-------------- Flash related variables -------------
            int m_reco_num_templates;
            std::vector<double> m_reco_template;


            //-------------- Flash related variables -------------
            std::vector<double> m_reco_flash_total_pe;
            std::vector<double> m_reco_flash_time;
            std::vector<double> m_reco_flash_time_width;
            std::vector<double> m_reco_flash_abs_time;
            std::vector<int>    m_reco_flash_frame;
            std::vector<double> m_reco_flash_ycenter;
            std::vector<double> m_reco_flash_ywidth;
            std::vector<double> m_reco_flash_zcenter;
            std::vector<double> m_reco_flash_zwidth;
            std::vector<double> m_reco_flash_total_pe_in_beamgate;
            std::vector<double> m_reco_flash_time_in_beamgate;
            std::vector<double> m_reco_flash_ycenter_in_beamgate;
            std::vector<double> m_reco_flash_zcenter_in_beamgate;

            int m_reco_num_flashes;
            int m_reco_num_flashes_in_beamgate;

            double m_beamgate_flash_start;
            double m_beamgate_flash_end;


            //----------- CRT related variables -----------------
            //fields storing information about the CRT hit closest to the flash
            double m_CRT_min_hit_time;
            double m_CRT_min_hit_PE;
            double m_CRT_min_hit_x;
            double m_CRT_min_hit_y;
            double m_CRT_min_hit_z;

            //Fields storing information about all CRT hits in event
            std::vector<double> m_CRT_hits_time;
            std::vector<double> m_CRT_hits_PE;
            std::vector<double> m_CRT_hits_x;
            std::vector<double> m_CRT_hits_y;
            std::vector<double> m_CRT_hits_z;

            double m_CRT_dt; //time between flash and nearest CRT hit

            //------------ Track related Variables -------------
            int m_reco_asso_tracks;
            std::vector<int> m_reco_track_num_daughters;
            std::vector<double> m_reco_track_daughter_trackscore;
            std::vector<double> m_reco_track_length;
            std::vector<double> m_reco_track_dirx;
            std::vector<double> m_reco_track_diry;
            std::vector<double> m_reco_track_dirz;
            std::vector<double> m_reco_track_startx;
            std::vector<double> m_reco_track_starty;
            std::vector<double> m_reco_track_startz;
            std::vector<double> m_reco_track_endx;
            std::vector<double> m_reco_track_endy;
            std::vector<double> m_reco_track_endz;
            std::vector<double>   m_reco_track_theta_yz;
            std::vector<double>   m_reco_track_phi_yx;


            std::vector<int> m_reco_track_num_trajpoints;
            std::vector<int> m_reco_track_num_spacepoints;
            std::vector<double> m_reco_track_proton_kinetic_energy;
            std::vector<size_t>  m_reco_track_ordered_energy_index;
            std::vector<size_t>  m_reco_track_ordered_displacement_index;
            std::vector<double> m_reco_track_spacepoint_principal0;
            std::vector<double> m_reco_track_spacepoint_principal1;
            std::vector<double> m_reco_track_spacepoint_principal2;

            std::vector<double> m_reco_track_spacepoint_chi;
            std::vector<double> m_reco_track_spacepoint_max_dist;

            std::vector<int> m_reco_track_best_calo_plane;

            std::vector<double> m_reco_track_mean_dEdx_best_plane;
            std::vector<double> m_reco_track_mean_dEdx_start_half_best_plane;
            std::vector<double> m_reco_track_mean_dEdx_end_half_best_plane;
            std::vector<int> m_reco_track_good_calo_best_plane;
            std::vector<std::vector<double>> m_reco_track_trunc_dEdx_best_plane;
            std::vector<double> m_reco_track_mean_trunc_dEdx_best_plane;
            std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_best_plane;
            std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_best_plane;
            std::vector<double> m_reco_track_trunc_PIDA_best_plane;
            std::vector<std::vector<double>> m_reco_track_resrange_best_plane;
            std::vector<std::vector<double>> m_reco_track_dEdx_best_plane;



            std::vector<double> m_reco_track_mean_dEdx_p0;
            std::vector<double> m_reco_track_mean_dEdx_start_half_p0;
            std::vector<double> m_reco_track_mean_dEdx_end_half_p0;
            std::vector<int> m_reco_track_good_calo_p0;
            std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p0;
            std::vector<double> m_reco_track_mean_trunc_dEdx_p0;
            std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p0;
            std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p0;
            std::vector<double> m_reco_track_trunc_PIDA_p0;
            std::vector<std::vector<double>> m_reco_track_resrange_p0;
            std::vector<std::vector<double>> m_reco_track_dEdx_p0;

            std::vector<double> m_reco_track_mean_dEdx_p1;
            std::vector<double> m_reco_track_mean_dEdx_start_half_p1;
            std::vector<double> m_reco_track_mean_dEdx_end_half_p1;
            std::vector<int> m_reco_track_good_calo_p1;
            std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p1;
            std::vector<double> m_reco_track_mean_trunc_dEdx_p1;
            std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p1;
            std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p1;
            std::vector<double> m_reco_track_trunc_PIDA_p1;
            std::vector<std::vector<double>> m_reco_track_resrange_p1;
            std::vector<std::vector<double>> m_reco_track_dEdx_p1;

            std::vector<double> m_reco_track_mean_dEdx_p2;
            std::vector<double> m_reco_track_mean_dEdx_start_half_p2;
            std::vector<double> m_reco_track_mean_dEdx_end_half_p2;
            std::vector<int> m_reco_track_good_calo_p2;
            std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p2;
            std::vector<double> m_reco_track_mean_trunc_dEdx_p2;
            std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p2;
            std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p2;
            std::vector<double> m_reco_track_trunc_PIDA_p2;
            std::vector<std::vector<double>> m_reco_track_resrange_p2;
            std::vector<std::vector<double>> m_reco_track_dEdx_p2;


            std::vector<int> m_reco_track_num_calo_hits_p0;
            std::vector<int> m_reco_track_num_calo_hits_p1;
            std::vector<int> m_reco_track_num_calo_hits_p2;


            std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane0;
            std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane1;
            std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane2;

	    std::vector<int> m_reco_track_sliceId; //the slice id for the slice continaing the reco track
	    std::vector<double> m_reco_track_nuscore; //the neutrino score of the slice containing the reco track
            std::vector<bool> m_reco_track_isclearcosmic;//true if reco track is in a clear cosmic slice
            std::vector<double> m_reco_track_trackscore;
            std::vector<bool> m_reco_track_is_nuslice;


            std::vector<int> m_sim_track_matched;
            std::vector<double> m_sim_track_overlay_fraction;
            std::vector<double> m_sim_track_energy;
            std::vector<double> m_sim_track_mass;
            std::vector<double> m_sim_track_kinetic_energy;
            std::vector<int> m_sim_track_pdg;
            std::vector<int> m_sim_track_parent_pdg;
            std::vector<int> m_sim_track_origin;
            std::vector<std::string> m_sim_track_process;
            std::vector<double> m_sim_track_startx;
            std::vector<double> m_sim_track_starty;
            std::vector<double> m_sim_track_startz;
            std::vector<int> m_sim_track_trackID;

            std::vector<int> m_sim_track_sliceId; //the slice id for the slice continaing the sim track
            std::vector<double> m_sim_track_nuscore; //the neutrino score of the slice containing the sim track
            std::vector<bool> m_sim_track_isclearcosmic;//true if sim track is in a clear cosmic slice

            /*-------------------------------------------------------------------------------------*/
            std::vector<double> m_isolation_min_dist_trk_shr;	
            std::vector<double> m_isolation_min_dist_trk_unassoc;

            std::vector<double> m_isolation_num_shr_hits_win_1cm_trk;	    
            std::vector<double> m_isolation_num_shr_hits_win_2cm_trk;	    
            std::vector<double> m_isolation_num_shr_hits_win_5cm_trk;	    
            std::vector<double> m_isolation_num_shr_hits_win_10cm_trk;	    

            std::vector<double> m_isolation_num_unassoc_hits_win_1cm_trk;	    
            std::vector<double> m_isolation_num_unassoc_hits_win_2cm_trk;	    
            std::vector<double> m_isolation_num_unassoc_hits_win_5cm_trk;	    
            std::vector<double> m_isolation_num_unassoc_hits_win_10cm_trk;	    

            std::vector<double> m_isolation_nearest_shr_hit_to_trk_wire;	
            std::vector<double> m_isolation_nearest_shr_hit_to_trk_time;	

            std::vector<double> m_isolation_nearest_unassoc_hit_to_trk_wire;	
            std::vector<double> m_isolation_nearest_unassoc_hit_to_trk_time;

            /*-------------------------------------------------------------------------------------*/

            //------------ Shower related Variables  -------------

            std::vector<int> m_reco_shower_num_daughters;
            std::vector<double> m_reco_shower_daughter_trackscore;

            std::vector<int>   m_reco_shower3d_exists;
            std::vector<double>   m_reco_shower3d_startx;
            std::vector<double>   m_reco_shower3d_starty;
            std::vector<double>   m_reco_shower3d_startz;
            std::vector<double>   m_reco_shower3d_dirx;
            std::vector<double>   m_reco_shower3d_diry;
            std::vector<double>   m_reco_shower3d_dirz;
            std::vector<double>   m_reco_shower3d_theta_yz;
            std::vector<double>   m_reco_shower3d_phi_yx;

            std::vector<double> m_reco_shower3d_openingangle;
            std::vector<double> m_reco_shower3d_length;
            std::vector<double> m_reco_shower3d_conversion_distance;

            std::vector<double>   m_reco_shower3d_impact_parameter;
            std::vector<double>    m_reco_shower3d_implied_dirx;
            std::vector<double>     m_reco_shower3d_implied_diry;
            std::vector<double>     m_reco_shower3d_implied_dirz;

            std::vector<double> m_reco_shower3d_energy_plane0;
            std::vector<double> m_reco_shower3d_energy_plane1;
            std::vector<double> m_reco_shower3d_energy_plane2;

            std::vector<double> m_reco_shower3d_dEdx_plane0;
            std::vector<double> m_reco_shower3d_dEdx_plane1;
            std::vector<double> m_reco_shower3d_dEdx_plane2;



            std::vector<double>   m_reco_shower_startx;
            std::vector<double>   m_reco_shower_starty;
            std::vector<double>   m_reco_shower_startz;
            std::vector<double>   m_reco_shower_dirx;
            std::vector<double>   m_reco_shower_diry;
            std::vector<double>   m_reco_shower_dirz;
            std::vector<double>   m_reco_shower_theta_yz;
            std::vector<double>   m_reco_shower_phi_yx;

            std::vector<double> m_reco_shower_openingangle;
            std::vector<double> m_reco_shower_length;
            std::vector<double> m_reco_shower_conversion_distance;

            std::vector<double>   m_reco_shower_impact_parameter;
            std::vector<double>    m_reco_shower_implied_dirx;
            std::vector<double>     m_reco_shower_implied_diry;
            std::vector<double>     m_reco_shower_implied_dirz;

            std::vector<int> m_reco_shower_delaunay_num_triangles_plane0;
            std::vector<int> m_reco_shower_delaunay_num_triangles_plane1;
            std::vector<int> m_reco_shower_delaunay_num_triangles_plane2;

            std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane0;
            std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane1;
            std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane2;


            //shower flash matching

            std::vector<double> m_reco_shower_flash_shortest_distz;
            std::vector<double> m_reco_shower_flash_shortest_disty;
            std::vector<double> m_reco_shower_flash_shortest_distyz;

            std::vector<int> m_reco_shower_flash_shortest_index_z;
            std::vector<int> m_reco_shower_flash_shortest_index_y;
            std::vector<int> m_reco_shower_flash_shortest_index_yz;

            //end flash matching
            std::vector<int> m_reco_shower_num_hits_plane0;
            std::vector<int> m_reco_shower_num_hits_plane1;
            std::vector<int> m_reco_shower_num_hits_plane2;


            std::vector<double> m_reco_shower_delaunay_area_plane0;
            std::vector<double> m_reco_shower_delaunay_area_plane1;
            std::vector<double> m_reco_shower_delaunay_area_plane2;

            std::vector<int> m_reco_shower_sliceId; //the slice id for the slice continaing the reco shower
            std::vector<double> m_reco_shower_nuscore; //the neutrino score of the slice containing the reco shower
            std::vector<bool> m_reco_shower_isclearcosmic;//true if reco shower is in a clear cosmic slice
            std::vector<bool> m_reco_shower_is_nuslice;//true if reco shower is in a clear cosmic slice
            std::vector<double> m_reco_shower_trackscore;

            std::vector<double> m_reco_shower_kalman_exists;
            std::vector<double>   m_reco_shower_kalman_median_dEdx_plane0;
            std::vector<double>     m_reco_shower_kalman_median_dEdx_plane1;
            std::vector<double>   m_reco_shower_kalman_median_dEdx_plane2;
            std::vector<double>   m_reco_shower_kalman_median_dEdx_allplane;
            std::vector<double>      m_reco_shower_kalman_mean_dEdx_plane0;
            std::vector<double>    m_reco_shower_kalman_mean_dEdx_plane1;
            std::vector<double>    m_reco_shower_kalman_mean_dEdx_plane2;




            std::vector<int> m_sim_shower_matched;
            std::vector<double> m_sim_shower_energy;
            std::vector<double> m_sim_shower_kinetic_energy;
            std::vector<double> m_sim_shower_mass;
            std::vector<int> m_sim_shower_pdg;
            std::vector<int> m_sim_shower_trackID;
            std::vector<int> m_sim_shower_parent_pdg;
            std::vector<int> m_sim_shower_parent_trackID;
            std::vector<int> m_sim_shower_origin;
            std::vector<std::string> m_sim_shower_process;
            std::vector<std::string> m_sim_shower_end_process;
            std::vector<double> m_sim_shower_start_x;
            std::vector<double> m_sim_shower_start_y;
            std::vector<double> m_sim_shower_start_z;
            std::vector<double> m_sim_shower_vertex_x;
            std::vector<double> m_sim_shower_vertex_y;
            std::vector<double> m_sim_shower_vertex_z;

            std::vector<double> m_sim_shower_px;
            std::vector<double> m_sim_shower_py;
            std::vector<double> m_sim_shower_pz;


            std::vector<int> m_sim_shower_is_true_shower;
            std::vector<int> m_sim_shower_best_matched_plane;
            std::vector<double> m_sim_shower_matched_energy_fraction_plane0;
            std::vector<double> m_sim_shower_matched_energy_fraction_plane1;
            std::vector<double> m_sim_shower_matched_energy_fraction_plane2;
            std::vector<double> m_sim_shower_overlay_fraction;

            std::vector<int> m_sim_shower_sliceId; //the slice id for the slice continaing the sim shower matched to reco
            std::vector<double> m_sim_shower_nuscore; //the neutrino score of the slice containing the sim shower matched to reco
            std::vector<bool> m_sim_shower_isclearcosmic;//true if sim shower matched to reco is in a clear cosmic slice
            std::vector<bool> m_sim_shower_is_nuslice;//true if sim shower matched to reco is in a clear cosmic slice



            //------------ MCTruth related Variables  -------------
            int m_mctruth_num;
            int m_mctruth_origin;
            double m_mctruth_nu_E;
            double m_mctruth_nu_vertex_x;
            double m_mctruth_nu_vertex_y;
            double m_mctruth_nu_vertex_z;
            double m_mctruth_lepton_E;
            int m_mctruth_nu_pdg;
            int m_mctruth_lepton_pdg;
            int m_mctruth_mode ;
            int m_mctruth_interaction_type ;
            int m_mctruth_ccnc;
            double m_mctruth_qsqr;

            int m_mctruth_num_daughter_particles;
            std::vector<int> m_mctruth_daughters_pdg;
            std::vector<double> m_mctruth_daughters_E;

            int     m_mctruth_num_exiting_photons ;
            int      m_mctruth_num_exiting_protons ;
            int    m_mctruth_num_exiting_pi0 ;
            int   m_mctruth_num_exiting_pipm ;
            int   m_mctruth_num_exiting_neutrons; 
            int   m_mctruth_num_exiting_delta0; 
            int   m_mctruth_num_exiting_deltapm; 
            int   m_mctruth_num_exiting_deltapp; 

            double m_mctruth_leading_exiting_proton_energy;

            int m_mctruth_is_delta_radiative;
            int m_mctruth_delta_radiative_1g1p_or_1g1n;
            double m_mctruth_delta_photon_energy;
            double m_mctruth_delta_proton_energy;
            double m_mctruth_delta_neutron_energy;
            std::vector<int> m_mctruth_exiting_delta0_num_daughters;

            std::vector<int> m_mctruth_exiting_photon_trackID;
            std::vector<int> m_mctruth_exiting_photon_mother_trackID;
            std::vector<int> m_mctruth_exiting_photon_from_delta_decay;
            std::vector<double> m_mctruth_exiting_photon_energy;
            std::vector<int> m_mctruth_exiting_proton_trackID;
            std::vector<int> m_mctruth_exiting_proton_mother_trackID;
            std::vector<int> m_mctruth_exiting_proton_from_delta_decay;
            std::vector<double> m_mctruth_exiting_proton_energy;

            int  m_mctruth_num_reconstructable_protons;

            bool  m_mctruth_is_reconstructable_1g1p;
            bool  m_mctruth_is_reconstructable_1g0p;


            std::vector<double>        m_mctruth_exiting_pi0_E;
            std::vector<double>        m_mctruth_exiting_pi0_px;
            std::vector<double>        m_mctruth_exiting_pi0_py;
            std::vector<double>        m_mctruth_exiting_pi0_pz;

            double m_mctruth_pi0_leading_photon_energy;
            double m_mctruth_pi0_subleading_photon_energy;

            std::string  m_truthmatching_signaldef;

            //the calo calculated quantities 
            std::vector<double> m_reco_shower_energy_max; //for each hit in a shower, converts Q->E, and sums
            std::vector<double> m_reco_shower_energy_plane0;
            std::vector<double> m_reco_shower_energy_plane1;
            std::vector<double> m_reco_shower_energy_plane2;

            std::vector<double> m_reco_shower_plane0;
            std::vector<double> m_reco_shower_plane1;
            std::vector<double> m_reco_shower_plane2;

            std::vector<double> m_reco_shower_plane0_nhits;
            std::vector<double> m_reco_shower_plane1_nhits;
            std::vector<double> m_reco_shower_plane2_nhits;


            std::vector<size_t>  m_reco_shower_ordered_energy_index;
            std::vector<std::vector<double>> m_reco_shower_dQdx_plane0; //for each shower, looks at the hits for all clusters in the plane, stores the dQ/dx for each hit 
            std::vector<std::vector<double>> m_reco_shower_dQdx_plane1;
            std::vector<std::vector<double>> m_reco_shower_dQdx_plane2;
            std::vector<std::vector<double>> m_reco_shower_dEdx_plane0; //dE/dx from the calculated dQ/dx for each hit in shower on plane 	
            std::vector<std::vector<double>> m_reco_shower_dEdx_plane1;
            std::vector<std::vector<double>> m_reco_shower_dEdx_plane2;

            std::vector<double> m_reco_shower_dEdx_plane0_mean;
            std::vector<double> m_reco_shower_dEdx_plane1_mean;
            std::vector<double> m_reco_shower_dEdx_plane2_mean;
            std::vector<double> m_reco_shower_dEdx_plane0_max;
            std::vector<double> m_reco_shower_dEdx_plane1_max;
            std::vector<double> m_reco_shower_dEdx_plane2_max;
            std::vector<double> m_reco_shower_dEdx_plane0_min;
            std::vector<double> m_reco_shower_dEdx_plane1_min;
            std::vector<double> m_reco_shower_dEdx_plane2_min;
            std::vector<double> m_reco_shower_dEdx_plane0_median;
            std::vector<double> m_reco_shower_dEdx_plane1_median;
            std::vector<double> m_reco_shower_dEdx_plane2_median;

            std::vector<double>  m_reco_shower_angle_wrt_wires_plane0;
            std::vector<double>  m_reco_shower_angle_wrt_wires_plane1;
            std::vector<double>  m_reco_shower_angle_wrt_wires_plane2;

            std::vector<double>  m_reco_shower_dEdx_amalgamated;
            std::vector<int>  m_reco_shower_dEdx_amalgamated_nhits;

            std::vector<double> m_reco_shower_dQdx_plane0_median;
            std::vector<double> m_reco_shower_dQdx_plane1_median;
            std::vector<double> m_reco_shower_dQdx_plane2_median;

            std::vector<double> m_reco_shower_dEdx_plane0_nhits;
            std::vector<double> m_reco_shower_dEdx_plane1_nhits;
            std::vector<double> m_reco_shower_dEdx_plane2_nhits;

            double _time2cm;//value modeled from David's shower code

            // PID-related variables
            std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane0;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane1;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane2;
            std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane0;
            std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane1;
            std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane2;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane0;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane1;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane2;
            std::vector<double> m_reco_track_pid_pida_plane0;
            std::vector<double> m_reco_track_pid_pida_plane1;
            std::vector<double> m_reco_track_pid_pida_plane2;
            std::vector<double> m_reco_track_pid_chi2_mu_plane0;
            std::vector<double> m_reco_track_pid_chi2_mu_plane1;
            std::vector<double> m_reco_track_pid_chi2_mu_plane2;
            std::vector<double> m_reco_track_pid_chi2_p_plane0;
            std::vector<double> m_reco_track_pid_chi2_p_plane1;
            std::vector<double> m_reco_track_pid_chi2_p_plane2;
            std::vector<double> m_reco_track_pid_three_plane_proton_pid;


            double m_genie_spline_weight;
            bool Pi0PreselectionFilter();

	    // newly added
	    bool PFP_have_nuslice;

  };


    template <typename T>
        std::vector<size_t> sort_indexes(const std::vector<T> &v) {

            std::vector<size_t> idx(v.size());
            std::iota(idx.begin(), idx.end(), 0);

            // sort indexes based on comparing values in v
            std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

            return idx;
        }

    double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo ){
        double wire = geo.WireCoordinate(Y, Z, plane, fTPC, fCryostat);
        return wire;
    }


    double calcTime(double X,int plane,int fTPC,int fCryostat, detinfo::DetectorProperties const& detprop){
        double time = detprop.ConvertXToTicks(X, plane, fTPC,fCryostat);
        return time;
    }

    struct sss_score{
        int plane;
        int cluster_label;
        double point_score;
        int n_hits;

        double mean_wire;
        double max_wire;
        double min_wire;
        double mean_tick;
        double max_tick;
        double min_tick;

        double close_tick;
        double close_wire;
        double angle;//w.r.t shower primary

        double impact_parameter;

        double max_dist_tick;
        double mean_dist_tick;
        double min_dist_tick;
        double max_dist_wire;
        double mean_dist_wire;
        double min_dist_wire;

        double mean_dist;
        double max_dist;
        double min_dist;

        double pca_0;
        double pca_1;
        double pca_theta;

        int n_wires;
        int n_ticks;

        bool pass;

        sss_score(int ip, int cl): plane(ip), cluster_label(cl){};
    };


    double dist_line_point( std::vector<double>&X1, std::vector<double>& X2, std::vector<double>& point){
      double x1 =X1.at(0);
      double y1 =X1.at(1);
      double z1 =X1.at(2);

      double x2 =X2.at(0);
      double y2 =X2.at(1);
      double z2 =X2.at(2);

      double x0 =point.at(0);
      double y0 =point.at(1);
      double z0 =point.at(2);

      double x10 = x1-x0;
      double y10 = y1-y0;
      double z10 = z1-z0;

      double x21 = x2-x1;
      double y21 = y2-y1;
      double z21 = z2-z1;

      double t = -(x10*x21+y10*y21+z10*z21)/fabs(x21*x21+y21*y21+z21*z21 );

      double d2 = pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2)+2*t*((x2-x1)*(x1-x0)+(y2-y1)*(y1-y0)+(z2-z1)*(z1-z0))+t*t*( pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));


      return sqrt(d2);

    }

    double distanceToNearestDeadWire(int plane, double Ypoint, double Zpoint,
				     const geo::GeometryCore * geom,
				     std::vector<std::pair<int,int>> & bad_channel_list_fixed_mcc9 ){

      double min_dist = 999999;

      for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
	int channel = bad_channel_list_fixed_mcc9[i].first;
	int is_ok = bad_channel_list_fixed_mcc9[i].second;
	if(is_ok>1)continue;

	auto wireids = geom->ChannelToWire(channel);
	auto result = geom->WireEndPoints(wireids[0]);

	//std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<<result.start().X()<<" "<<result.start().Y()<<" "<<result.start().Z()<<" "<<result.end().X()<<" "<<result.end().Y()<<" "<<result.end().Z()<<std::endl;                                                                                                               
	//                    std::cout<<wireids[0].Plane<<" "<<result.start().X()<<std::endl;                                                                                                                                                                                                                                             
	if(plane != (int)wireids[0].Plane) continue;

	std::vector<double> start = {0.0,result.start().Y(),result.start().Z()};
	std::vector<double> end = {0.0,result.end().Y(),result.end().Z()};
	std::vector<double> point = {0.0,Ypoint,Zpoint};
	double dist = dist_line_point(start,end,point);
	min_dist = std::min(dist,min_dist);
      }

      return min_dist;

    }



    class cluster {

        public:

            cluster(int ID, int plane, std::vector<std::vector<double>> &pts, std::vector<art::Ptr<recob::Hit>> &hits) :f_ID(ID), f_plane(plane), f_pts(pts), f_hits(hits) {

                f_npts = f_pts.size();
                if(pts.size() != hits.size()){
                    std::cerr<<"seaviewer::cluster, input hits and pts not alligned"<<std::endl;
                }
                std::vector<double> wires(f_npts);
                std::vector<double> ticks(f_npts);
                for(int p =0; p< f_npts; ++p){
                    wires[p]=f_pts[p][0];
                    ticks[p]=f_pts[p][1];
                }
                TGraph af_graph(f_npts,&wires[0],&ticks[0]);
                f_graph = af_graph;

            };

            int getID() {return f_ID;}
            int getN() {return f_npts;}
            int getPlane(){ return f_plane;}
            TGraph * getGraph(){ return &f_graph;}
            std::vector<art::Ptr<recob::Hit>>  getHits(){return f_hits;}
            int setSSScore(sss_score & scorein){ f_SSScore = &scorein; return 0;}
            sss_score * getSSScore(){return f_SSScore;}
        private:

            int f_ID;
            int f_npts;
            int f_plane;
            std::vector<std::vector<double>> f_pts;
            std::vector<art::Ptr<recob::Hit>> f_hits;
            TGraph f_graph;
            sss_score *f_SSScore;
    };



    DEFINE_ART_MODULE(CCKaonAnalyzerRebuild)

      } // namespace lar_pandora

#endif
//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows
