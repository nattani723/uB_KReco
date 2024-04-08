#ifndef PARTICLE_PRODUCER
#define PARTICLE_PRODUCER 

#include <bits/c++config.h>
#undef _GLIBCXX14_CONSTEXPR
#define _GLIBCXX14_CONSTEXPR

#include "/uboone/app/users/taniuchi/51_pandora/srcs/ubana/ubana/Filters/CCKaonFilter/CCKaonProducer_module.h"
#include "TMath.h"
#include <algorithm>
#include <functional>

#include "Api/PandoraApi.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArContent.h"

#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "/uboone/app/users/taniuchi/51_pandora/srcs/larpandoracontent/larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"


using namespace pandora;
namespace Kaon_Analyzer
{

  //void Kaon_Analyzer::CCKaonAnalyzer::ClearMCTruths(){
  //}

  //const double region_of_interest = 35;
  const double region_of_interest = 20;
  //const double theta0XZBinSize = 0.02;
  const double theta0XZBinSize = 0.005;
  //const int smoothing_window = 15;
  const int smoothing_window = 3;

  /*
  const double thetaBinSize = 0.03;
  const double costhetaBinSize = 0.03;
  const double phiBinSize = 0.03;
  */
  const double thetaBinSize = 0.06;
  const double costhetaBinSize = 0.06;
  const double phiBinSize = 0.06;

  const int num_bin = (int)(6.28/theta0XZBinSize);
  //const int num_bin_theta = (int)(6.28/thetaBinSize);
  const int num_bin_theta = (int)(3.14/thetaBinSize);
  const int num_bin_phi = (int)(6.28/phiBinSize);
  //const int num_bin_phi = (int)(6.28/phiBinSize);


  const double growing_fit_initial_length_munu = 10.;//20
  const double initial_fit_distance_to_line_munu = 1;//2

  //const double growing_fit_initial_length_pip = 60.;
  //const double initial_fit_distance_to_line_pip = 5;

  const double growing_fit_initial_length_pip = 10.;//15
  const double initial_fit_distance_to_line_pip =1;//1.5

  //const double growing_fit_initial_length = 55.;
  //const double initial_fit_distance_to_line = 10;

  //const double growing_fit_initial_length = 5.;
  //const double initial_fit_distance_to_line = 0.5;

  //const double growing_fit_initial_length = 1500.;
  //const double initial_fit_distance_to_line = 3000;

  const int min_initial_hits_found = 5;//10
  const int max_fitting_hits_munu = 15;
  const int max_fitting_hits_pip = 15;
  const int max_fitting_hits_micro = 15;//15

  // const int min_initial_hits_found = 20;
  //const int max_fitting_hits = 25;
  const int macro_sliding_fit_window =10;//5
  const int micro_sliding_fit_window = 5;
  const double growing_fit_segment_length = 5.;

    
  const double distance_to_line_munu = 0.75;
  const double hit_connection_distance_munu = 0.75;
  const double distance_to_line_pip = 0.75;
  const double hit_connection_distance_pip = 0.75;
  
  
  /*
  const double distance_to_line_munu = 1;
  const double hit_connection_distance_munu = 1.5;
  const double distance_to_line_pip = 1;
  const double hit_connection_distance_pip = 1.5;
  */

  //const double distance_to_line = 3;
  //const double hit_connection_distance = 1.5;

  const double distance_to_line_spineall = 1;//2
  const double hit_connection_distance_spineall = 5;//2 
  const double distance_to_line_micro = 1;
  const double hit_connection_distance_micro = 5;

  double distance_to_line = 1;
  double hit_connection_distance = 1;
  double max_fitting_hits = 15;

  /*
  const double growing_fit_initial_length = 50.;
  const double initial_fit_distance_to_line = 5.;

  //const double growing_fit_initial_length = 1500.;
  //const double initial_fit_distance_to_line = 3000;

  const int min_initial_hits_found = 10;
  const int max_fitting_hits = 20;
  const int macro_sliding_fit_window = 1000;//5
  const int micro_sliding_fit_window = 20;
  const double growing_fit_segment_length = 5.;
  const double distance_to_line = 0.5;
  const double hit_connection_distance = 0.5;
  */

  //------------------------------------------------------------------------------------------------------------------------------------------

  void CCKaonProducer::fillAngularDistributionMap3D( std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle, 
						     TVector3 Kend_candidate,
						     art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
						     std::map<int, std::map<int, double>> &angular_distribution_map_3D){
    

    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;

    const TVector3 Xaxis(1.,0.,0.); 

    for(size_t i_h=0; i_h<hits_from_pfparticle.size(); i_h++){
      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle[i_h].key());

      if(spacepoint_vec.size()!=1) continue;
      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      const TVector3 distance_vector = hit_position - Kend_candidate;
      
      if(distance_vector.Mag() > region_of_interest) continue;
      
      double theta = distance_vector.Theta();
      double sintheta = TMath::Sin(theta);
      double phi = distance_vector.Phi();

      int theta_factor = (int)(std::floor(theta / thetaBinSize));
      int phi_factor = (int)(std::floor(phi / phiBinSize));

      //cout << "theta: " <<  distance_vector.Theta() << ", factor: " << theta_factor << endl;
      //cout << "phi: " <<  distance_vector.Phi() << ", factor: " << phi_factor << endl;

      if(angular_distribution_map_3D.find(theta_factor) == angular_distribution_map_3D.end()){
	if(angular_distribution_map_3D[theta_factor].find(phi_factor) == angular_distribution_map_3D[theta_factor].end()){
	  angular_distribution_map_3D[theta_factor][phi_factor] = sintheta;
	}
      }
      else angular_distribution_map_3D[theta_factor][phi_factor] += sintheta;
      
    }

  }


  //-----------------------------------------------------------------------------------------------------------------------------------------------------


  void CCKaonProducer::smoothAngularDistributionMap3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D){
   

    const int loop_min = (-1) * (smoothing_window - 1) / 2;
    const int loop_max = (smoothing_window - 1) / 2;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_temp(angular_distribution_map_3D);
    angular_distribution_map_3D.clear();

    for (auto const& angular_distribution_theta : angular_distribution_map_3D_temp){
      const int current_bin_theta(angular_distribution_theta.first);
      
      for(auto const& angular_distribution_phi : angular_distribution_theta.second){
	const int current_bin_phi(angular_distribution_phi.first);
	int bin_count = 0;
	double total = 0;
	
	for (int bin_offset_theta = loop_min; bin_offset_theta <= loop_max; ++bin_offset_theta){
	  const int contributing_bin_theta = current_bin_theta + bin_offset_theta;
	  
	  for (int bin_offset_phi = loop_min; bin_offset_phi <= loop_max; ++bin_offset_phi){
	    ++bin_count;
	    const int contributing_bin_phi = current_bin_phi + bin_offset_phi;
	    
	    if(angular_distribution_map_3D_temp[ contributing_bin_theta ].find( contributing_bin_phi ) == angular_distribution_map_3D_temp[ contributing_bin_theta ].end()) total += 0;
	    else total += angular_distribution_map_3D_temp[ contributing_bin_theta ].at( contributing_bin_phi );
	    
	  }
	}
	angular_distribution_map_3D[ current_bin_theta ][ current_bin_phi ] = total/(double)bin_count;
      }
    }
    
  }


  //------------------------------------------------------------------------------------------------------------------------------------------



  void CCKaonProducer::accumulateAngularDistributionMap3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D_1,
							  std::map<int, std::map<int, double>> &angular_distribution_map_3D_2,
							  std::map<int, std::map<int, double>> &angular_distribution_map_3D){
    
    //std::map<int, std::map<int, double>> angular_distribution_map_3D = std::accumulate( angular_distribution_map_3D_1.begin(), angular_distribution_map_3D_1.end(), 

    angular_distribution_map_3D = angular_distribution_map_3D_1; //copy

    for( auto const& angular_distribution_map_3D_theta_2 : angular_distribution_map_3D_2){
      for( auto const& angular_distribution_map_3D_phi_2 : angular_distribution_map_3D_theta_2.second){
	if( angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first].find( angular_distribution_map_3D_phi_2.first ) == angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first].end() ) angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] = angular_distribution_map_3D_phi_2.second;
	else{
	  angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] += angular_distribution_map_3D_phi_2.second;
	  cout << "ACCUMULATE "  << angular_distribution_map_3D_1[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] << " " << angular_distribution_map_3D_2[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] << " " << angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] << endl;

	}
      }
    }
  }



  //------------------------------------------------------------------------------------------------------------------------------------------



  void CCKaonProducer::obtainPeakVector3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D,
					  vector<bool>& v_trk_flg_peak,
					  std::map<double, TVector2, std::greater<>>& view_peak_map,
					  bool trk_flg){
    //vector<TVector2>& view_peak_vector,

    
    for (auto const& angular_distribution_theta : angular_distribution_map_3D){
      
      for(auto const& angular_distribution_phi : angular_distribution_theta.second){
	  
	const int bin_theta = angular_distribution_theta.first;
	const int bin_phi = angular_distribution_phi.first;
	const double bin_weight = angular_distribution_phi.second ;
	
	int preceding_bin_theta = bin_theta - 1;
	int preceding_bin_phi = bin_phi - 1;
	int following_bin_theta = bin_theta + 1;
	int following_bin_phi = bin_phi + 1;
	bool FoundBin = false;
	bool MaxWeight = true;
	  
	for (int scan_bin_theta = preceding_bin_theta; scan_bin_theta <= following_bin_theta; ++scan_bin_theta){
	  for (int scan_bin_phi = preceding_bin_phi; scan_bin_phi <= following_bin_phi; ++scan_bin_phi){

	    if(scan_bin_theta == bin_theta && scan_bin_phi == bin_phi) continue;
	    
	    if(angular_distribution_map_3D[ scan_bin_theta ].find( scan_bin_phi ) == angular_distribution_map_3D[ scan_bin_theta ].end() ||
	       (std::fabs(angular_distribution_map_3D[ scan_bin_theta ].at( scan_bin_phi ) - angular_distribution_map_3D[ bin_theta ].at( bin_phi )) > std::numeric_limits<double>::epsilon())){
	      FoundBin = true;
	    }else{
	      FoundBin = false;
	      break;
	    }
	  }
	  if(FoundBin == false) break;
	}
	
	
	for (int scan_bin_theta = preceding_bin_theta; scan_bin_theta <= following_bin_theta; ++scan_bin_theta){
	  for (int scan_bin_phi = preceding_bin_phi; scan_bin_phi <= following_bin_phi; ++scan_bin_phi){
	    double scan_bin_weight;
	    if( angular_distribution_map_3D[ scan_bin_theta ].find( scan_bin_phi ) == angular_distribution_map_3D[ scan_bin_theta ].end() ){
	      scan_bin_weight = 0;
	    }else{
	      scan_bin_weight = angular_distribution_map_3D[scan_bin_theta].at(scan_bin_phi);
	    }
	    if(bin_weight < scan_bin_weight) MaxWeight = false;
	  }
	}
	
	if(MaxWeight == false) continue;
	

	TVector2 peak_bins;
	peak_bins.Set((double)bin_theta, (double)bin_phi);
	//std::map<int, int> peak_bins;
	//peak_bins[bin_theta] = bin_phi;
	//view_peak_vector.push_back(peak_bins);
	view_peak_map[bin_weight] = peak_bins;
	v_trk_flg_peak.push_back(trk_flg);
	
      }
    }
  }
  

  //------------------------------------------------------------------------------------------------------------------------------------------

  
  void CCKaonProducer::findBestAngularPeak3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D,
					     std::map<double, TVector2, std::greater<>>& view_peak_map,
					     vector<TVector2> &best_peak_bins){

    //double open_angle;
    if(view_peak_map.size()==0) return;

    const TVector2 best_peak_bin = view_peak_map.begin()->second;
    if(view_peak_map.begin()->first < 1) return;
    best_peak_bins.push_back(best_peak_bin);
    cout << "view_peak_map[0] is " << view_peak_map.begin()->first << " " << view_peak_map.begin()->second.X() << endl;
    cout << "best_peak_bin is " << best_peak_bin.X()*thetaBinSize << " " << best_peak_bin.Y()*thetaBinSize << endl;

    TVector3 best_peak_dir;
    Double_t theta_best = best_peak_bin.X()*thetaBinSize;
    Double_t phi_best = best_peak_bin.Y()*phiBinSize;

    best_peak_dir.SetMagThetaPhi(1., theta_best, phi_best);
    //best_peak_dir.SetTheta(theta_best);
    //best_peak_dir.SetPhi(phi_best);

    for (auto const& entry : view_peak_map){
      
    TVector2 theta_phi_bin = entry.second;
    TVector3 view_peak_dir;
    Double_t theta_view = theta_phi_bin.X()*thetaBinSize;
    Double_t phi_view = theta_phi_bin.Y()*phiBinSize;

    view_peak_dir.SetMagThetaPhi(1., theta_view, phi_view);
    //view_peak_dir.SetMag(1.);
    //view_peak_dir.SetTheta(theta_view);
    //view_peak_dir.SetPhi(phi_view);
      //view_peak_dir.SetXYZ(TMath::Sin(theta_phi_bin.X())*TMath::Cos(theta_phi_bin.Y()), TMath::Sin(theta_phi_bin.X())*TMath::Sin(theta_phi_bin.Y()), TMath::Cos(theta_phi_bin.X()));

      //for(unsigned int i_peak=0; i_peak<best_peak_bins.size(); i_peak++){
      if(!(theta_phi_bin.X() == best_peak_bin.X() && theta_phi_bin.Y() == best_peak_bin.Y())){
	//TVector2 best_peak_bin = best_peak_bins[i_peak];
	//TVector3 best_peak_dir;
	//best_peak_bin.SetXYZ( TMath::Sin(best_peak_bin.X())*TMath::Cos(best_peak_bin.Y()), TMath::Sin(best_peak_bin.X())*TMath::Sin(best_peak_bin.Y()), TMath::Cos(best_peak_bin.X()) );

	double open_angle = view_peak_dir.Angle(best_peak_dir);
	//if(open_angle > TMath::Pi()*1/2 && entry.first>0.4) best_peak_bins.push_back(theta_phi_bin);
	if(open_angle > TMath::Pi()*2/3 && entry.first>0.4) best_peak_bins.push_back(theta_phi_bin);
	//best_peak_bins.push_back(theta_phi_bin);

	//if(open_angle > TMath::Pi()*0.5) best_peak_bins.push_back(theta_phi_bin);
      }

      if(best_peak_bins.size() >= 3) break; //just look second highest peak with large opening angle

    }

    for(auto const& entry : best_peak_bins){
      cout << "best peaks are: " << entry.X()*thetaBinSize << " " << entry.Y()*phiBinSize << endl;
    }
    /*
    if(best_peak_bins.size()==2){
      TVector3 best1;
      TVector3 best2;
      best1.SetMagThetaPhi(1, best_peak_bins[0].X()*thetaBinSize, best_peak_bins[0].Y()*thetaBinSize);
      best2.SetMagThetaPhi(1, best_peak_bins[1].X()*thetaBinSize, best_peak_bins[1].Y()*thetaBinSize);
      cout << "TMath::Pi()*0.5 is " << TMath::Pi()*0.5 << " and open_angle is " << best1.Angle(best2) << endl;
    }
    */

  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  
  void CCKaonProducer::findShowerSpine3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
					 art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
					 std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
					 std::vector<std::vector<art::Ptr<recob::Hit>>>& shower_spine_hit_list_vector,
					 TVector3 Kend_candidate,
					 std::map<double, TVector2, std::greater<>>& view_peak_map,
					 vector<TVector2> &best_peak_bins){

    art::ServiceHandle<geo::Geometry> theGeometry;
    const unsigned int nWirePlanes(theGeometry->MaxPlanes());

    if (nWirePlanes > 3)
        throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- More than three wire planes present ";

    if ((0 == theGeometry->Ncryostats()) || (0 == theGeometry->NTPC(0)))
        throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- unable to access first tpc in first cryostat ";

    std::unordered_set<geo::_plane_proj> planeSet;
    for (unsigned int iPlane = 0; iPlane < nWirePlanes; ++iPlane)
        (void) planeSet.insert(theGeometry->TPC(0, 0).Plane(iPlane).View());

    if ((nWirePlanes != planeSet.size()) || !planeSet.count(geo::kU) || !planeSet.count(geo::kV) || (planeSet.count(geo::kW) && planeSet.count(geo::kY)))
        throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- expect to find u and v views; if there is one further view, it must be w or y ";

    const bool useYPlane((nWirePlanes > 2) && planeSet.count(geo::kY));  

    const float wirePitchU(theGeometry->WirePitch(geo::kU));
    const float wirePitchV(theGeometry->WirePitch(geo::kV));
    const float wirePitchW((nWirePlanes < 3) ? 0.5f * (wirePitchU + wirePitchV) : (useYPlane) ? theGeometry->WirePitch(geo::kY) : theGeometry->WirePitch(geo::kW));

    const float sliding_fit_pitch = wirePitchW;

    if(best_peak_bins.size()==0 || hits_from_pfparticle.size()==0 || spacepoint_per_hit.size()==0) return;

    for(auto& entry : best_peak_bins){

    double highest_l = 0.;
    vector<TVector3> running_fit_position_vec;
    pandora::CartesianPointVector pandora_running_fit_position_vec;
    TVector3 initial_direction;
    std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list;

    double theta = entry.X()*thetaBinSize;
    double phi = entry.Y()*phiBinSize;
    //cout << "theta: " << theta << ", phi: " << phi << endl;
    initial_direction.SetMagThetaPhi(1, theta, phi);
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;

    cout << "initial_direction: " << initial_direction.X() << " " << initial_direction.Y() << " " << initial_direction.Z() << endl;

    cout << "hits_from_pfparticle.size(): " << hits_from_pfparticle.size() << endl;
    for(size_t i_h=0; i_h<hits_from_pfparticle.size(); i_h++){

      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle[i_h].key());
      if(spacepoint_vec.size()!=1) continue;

      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
      const TVector3 distance_vector = hit_position - Kend_candidate;
      const double l = initial_direction.Dot(distance_vector);
      const double t = initial_direction.Cross(distance_vector).Mag();//Mag2??

      cout << "candidate of all hits" << endl;
      cout << "l: " << l << ", t: " << t << endl;
    
      if( (best_peak_bins.size()==1 &&  (l<growing_fit_initial_length_munu) && (l>0.) && (t<initial_fit_distance_to_line_munu)) ||
	  (best_peak_bins.size()>1 &&  (l<growing_fit_initial_length_pip) && (l>0.) && (t<initial_fit_distance_to_line_pip)) ){

      //if( (l<growing_fit_initial_length) && (l>0.) && (t<initial_fit_distance_to_line) ){
	if(l>highest_l) highest_l = l;

	if(std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), hits_from_pfparticle[i_h]) == unavailable_hit_list.end()) shower_spine_hit_list.push_back(hits_from_pfparticle[i_h]);

	running_fit_position_vec.push_back(hit_position);
	pandora_running_fit_position_vec.push_back(pandora_hit_position);
      }
    }

    cout << "running_fit_position_vec.size(): " << running_fit_position_vec.size() << endl;
    cout << "pandora_running_fit_position_vec.size(): " << pandora_running_fit_position_vec.size() << endl;

    // Require significant number of initial hits
    //initial_running_fit_position_vec_size = running_fit_position_vec.size();
    if (running_fit_position_vec.size() < min_initial_hits_found){
      cout << "Requiring significant number of initial hits" << endl;
      shower_spine_hit_list.clear();
      shower_spine_hit_list_vector.push_back(shower_spine_hit_list);
      continue;
      //return;
    }

    // Perform a running fit to follow pathway
    bool is_end_downstream = false;
    if(initial_direction.Z() > 0.) is_end_downstream = true;


    //const float sliding_fit_pitch = lar_content::LArGeometryHelper::GetWireZPitch(this->GetPandora());
    //const float sliding_fit_pitch = 0.30;

    TVector3 extrapolated_direction;
    //extrapolated_direction.SetMagThetaPhi(1, best_peak_bins[0].X() * thetaBinSize, best_peak_bins[0].Y() *phiBinSize);
    extrapolated_direction.SetMagThetaPhi(1, theta, phi);
    TVector3 extrapolated_start_position = Kend_candidate;
    TVector3 extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * highest_l);

    const pandora::CartesianVector pandora_vertex_position(extrapolated_start_position.X(), extrapolated_start_position.Y(), extrapolated_start_position.Z());
    lar_content::LArConnectionPathwayHelper::SortByDistanceToPoint vtx_distance(pandora_vertex_position);
    std::sort(pandora_running_fit_position_vec.begin(), pandora_running_fit_position_vec.end(), vtx_distance);
    TVector3 closest_hit;
    closest_hit.SetXYZ(pandora_running_fit_position_vec[0].GetX(), pandora_running_fit_position_vec[0].GetY(), pandora_running_fit_position_vec[0].GetZ());
    cout << "closest distance from Kend_candidate is " << (Kend_candidate - closest_hit).Mag() << endl;
        
    if( (Kend_candidate - closest_hit).Mag()>8 ){
      shower_spine_hit_list.clear();
      shower_spine_hit_list_vector.push_back(shower_spine_hit_list);
      continue;
    }


    unsigned int count = 0;
    bool hits_collected = true;

    while (hits_collected){

      ++count;

      //try{
	//const int excess_hits_in_fit = running_fit_position_vec.size() - max_fitting_hits;
	const int excess_hits_in_fit = pandora_running_fit_position_vec.size() - max_fitting_hits;
	cout << "excess_hits_in_fit: " << excess_hits_in_fit << endl;

	  const pandora::CartesianVector pandora_extrapolated_end_position(extrapolated_end_position.X(), extrapolated_end_position.Y(), extrapolated_end_position.Z());
	  lar_content::LArConnectionPathwayHelper::SortByDistanceToPoint smpl(pandora_extrapolated_end_position);
	  std::sort(pandora_running_fit_position_vec.begin(), pandora_running_fit_position_vec.end(), smpl);
	
	if(excess_hits_in_fit>0){

	  // Remove furthest away hits
	  cout << "before loop pandora_running_fit_position_vec.size(): " << pandora_running_fit_position_vec.size() << endl; 
	  for (int i = 0; i < excess_hits_in_fit; ++i){
	    //running_fit_position_vec.erase(std::prev(running_fit_position_vec.end()));
	    pandora_running_fit_position_vec.erase(std::prev(pandora_running_fit_position_vec.end()));
	  }

	  running_fit_position_vec.clear();
	  TVector3 hit_position_tmp;
	  for(auto pandora_running_fit_position : pandora_running_fit_position_vec){
	    hit_position_tmp.SetXYZ(pandora_running_fit_position.GetX(), pandora_running_fit_position.GetY(), pandora_running_fit_position.GetZ());
	    running_fit_position_vec.push_back(hit_position_tmp);
	  }

	}
	//cout << "after loop pandora_running_fit_position_vec.size(): " << pandora_running_fit_position_vec.size() << endl; 

	distance_to_line = distance_to_line_munu;
	hit_connection_distance = hit_connection_distance_munu;

	const lar_content::ThreeDSlidingFitResult extrapolated_fit(&pandora_running_fit_position_vec, macro_sliding_fit_window, sliding_fit_pitch);
	if(count == 1){
	  extrapolated_start_position = extrapolated_end_position;
	}else{
	  if(is_end_downstream){
	    extrapolated_start_position.SetXYZ(extrapolated_fit.GetGlobalMaxLayerPosition().GetX(), extrapolated_fit.GetGlobalMaxLayerPosition().GetY(), extrapolated_fit.GetGlobalMaxLayerPosition().GetZ());
	  }else{
	    extrapolated_start_position.SetXYZ(extrapolated_fit.GetGlobalMinLayerPosition().GetX(), extrapolated_fit.GetGlobalMinLayerPosition().GetY(), extrapolated_fit.GetGlobalMinLayerPosition().GetZ());
	  }
	}

	if(is_end_downstream){
	  extrapolated_direction.SetXYZ(extrapolated_fit.GetGlobalMaxLayerDirection().GetX(), extrapolated_fit.GetGlobalMaxLayerDirection().GetY(), extrapolated_fit.GetGlobalMaxLayerDirection().GetZ());
	  if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
	    extrapolated_direction *= -1.;
	    is_end_downstream = false;
	  }

	}else{
	  extrapolated_direction.SetXYZ(extrapolated_fit.GetGlobalMinLayerDirection().GetX()*(-1.), extrapolated_fit.GetGlobalMinLayerDirection().GetY()*(-1.), extrapolated_fit.GetGlobalMinLayerDirection().GetZ()*(-1.));
	  if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
	    extrapolated_direction *= -1.;
	    is_end_downstream = true;
	  }
	}

	extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * growing_fit_segment_length);
	cout << "extrapolated_direction for count " << count << " is " << extrapolated_direction.X() << " " << extrapolated_direction.Y() << " "  << extrapolated_direction.Z() << endl;


	//cout << "before collectSubsectionHits" << endl;
	hits_collected = this->collectSubsectionHits(extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, hits_from_pfparticle, spacepoint_per_hit, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, shower_spine_hit_list);
	//cout << "after collectSubsectionHits" << endl;


	///// as a final effort, fit with all spine hit list
	if(!hits_collected){

	distance_to_line = distance_to_line_spineall;
	hit_connection_distance = hit_connection_distance_spineall;

	  pandora::CartesianPointVector pandora_running_fit_position_spineall_vec;

	  for(size_t i_h=0; i_h<shower_spine_hit_list.size(); i_h++){
	    spacepoint_vec.clear();
	    spacepoint_vec = spacepoint_per_hit.at(shower_spine_hit_list[i_h].key());
	    if(spacepoint_vec.size()!=1) continue;

	    const TVector3 hit_position = spacepoint_vec[0]->XYZ();
	    const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
	    pandora_running_fit_position_spineall_vec.push_back(pandora_hit_position);
	  }

	  //const lar_content::ThreeDSlidingFitResult micro_extrapolated_fit(&pandora_running_fit_position_spineall_vec, micro_sliding_fit_window, sliding_fit_pitch);
	  const lar_content::ThreeDSlidingFitResult spineall_extrapolated_fit(&pandora_running_fit_position_spineall_vec, 20, sliding_fit_pitch);
	  
	  if(count == 1)
	    {
	      extrapolated_start_position = extrapolated_start_position;
	    }else{
	    if(is_end_downstream){
	      extrapolated_start_position.SetXYZ(spineall_extrapolated_fit.GetGlobalMaxLayerPosition().GetX(), spineall_extrapolated_fit.GetGlobalMaxLayerPosition().GetY(), spineall_extrapolated_fit.GetGlobalMaxLayerPosition().GetZ());
	    }else{
	      extrapolated_start_position.SetXYZ(spineall_extrapolated_fit.GetGlobalMinLayerPosition().GetX(), spineall_extrapolated_fit.GetGlobalMinLayerPosition().GetY(), spineall_extrapolated_fit.GetGlobalMinLayerPosition().GetZ());	
	    }
	  }
	  
	    if(is_end_downstream){
	      extrapolated_direction.SetXYZ(spineall_extrapolated_fit.GetGlobalMaxLayerDirection().GetX(), spineall_extrapolated_fit.GetGlobalMaxLayerDirection().GetY(), spineall_extrapolated_fit.GetGlobalMaxLayerDirection().GetZ());
	      if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
		extrapolated_direction *= -1.;
		is_end_downstream = false;
	      }
	    }else{
	      extrapolated_direction.SetXYZ(spineall_extrapolated_fit.GetGlobalMinLayerDirection().GetX()*(-1.), spineall_extrapolated_fit.GetGlobalMinLayerDirection().GetY()*(-1.), spineall_extrapolated_fit.GetGlobalMinLayerDirection().GetZ()*(-1.));
	      if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
		extrapolated_direction *= -1.;
		is_end_downstream = true;
	      }
	    }

	    extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * growing_fit_segment_length);
   cout << "extrapolated_direction for count " << count << " is " << extrapolated_direction.X() << " " << extrapolated_direction.Y() << " "  << extrapolated_direction.Z() << endl;
	    
	    hits_collected = this->collectSubsectionHits(spineall_extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, hits_from_pfparticle, spacepoint_per_hit, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, shower_spine_hit_list);

	}


	// As a final effort, reduce the sliding fit window
	if (!hits_collected)
	  {

	    const int excess_hits_in_fit_micro = pandora_running_fit_position_vec.size() - max_fitting_hits_micro;
	    cout << "excess_hits_in_fit_micro: " << excess_hits_in_fit_micro << endl;
	    
	    pandora::CartesianPointVector pandora_running_fit_position_micro_vec;
	    pandora_running_fit_position_micro_vec = pandora_running_fit_position_vec;

	if(excess_hits_in_fit_micro>0){
	  // Remove furthest away hits
	  for (int i = 0; i < excess_hits_in_fit_micro; ++i){
	    pandora_running_fit_position_micro_vec.erase(std::prev(pandora_running_fit_position_micro_vec.end()));
	  }
	}

	distance_to_line = distance_to_line_micro;
	hit_connection_distance = hit_connection_distance_micro;

	const lar_content::ThreeDSlidingFitResult micro_extrapolated_fit(&pandora_running_fit_position_micro_vec, micro_sliding_fit_window, sliding_fit_pitch);
	// const lar_content::ThreeDSlidingFitResult micro_extrapolated_fit(&pandora_running_fit_position_vec, 5., sliding_fit_pitch);

	    if(count == 1)
	      {
		extrapolated_start_position = extrapolated_start_position;
	    }else{
	      if(is_end_downstream){
		extrapolated_start_position.SetXYZ(micro_extrapolated_fit.GetGlobalMaxLayerPosition().GetX(), micro_extrapolated_fit.GetGlobalMaxLayerPosition().GetY(), micro_extrapolated_fit.GetGlobalMaxLayerPosition().GetZ());
	      }else{
		extrapolated_start_position.SetXYZ(micro_extrapolated_fit.GetGlobalMinLayerPosition().GetX(), micro_extrapolated_fit.GetGlobalMinLayerPosition().GetY(), micro_extrapolated_fit.GetGlobalMinLayerPosition().GetZ());	
	      }
	    }

	    if(is_end_downstream){
	      extrapolated_direction.SetXYZ(micro_extrapolated_fit.GetGlobalMaxLayerDirection().GetX(), micro_extrapolated_fit.GetGlobalMaxLayerDirection().GetY(), micro_extrapolated_fit.GetGlobalMaxLayerDirection().GetZ());
	      if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
		extrapolated_direction *= -1.;
		is_end_downstream = false;
	      }

	    }else{
	      extrapolated_direction.SetXYZ(micro_extrapolated_fit.GetGlobalMinLayerDirection().GetX()*(-1.), micro_extrapolated_fit.GetGlobalMinLayerDirection().GetY()*(-1.), micro_extrapolated_fit.GetGlobalMinLayerDirection().GetZ()*(-1.));
	      if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
		extrapolated_direction *= -1.;
		is_end_downstream = true;
	      }

	    }

	    extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * growing_fit_segment_length);
	    
	    hits_collected = this->collectSubsectionHits(micro_extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, hits_from_pfparticle, spacepoint_per_hit, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, shower_spine_hit_list);
	}

	//}
	/*
      catch (const StatusCodeException &){
	cout << "StatusCodeException was called in particle_split_basetool_producer" << endl;
	return;
      }
	*/
    }

    if(shower_spine_hit_list.size()==0) cout << "shower_spine_hit_list.size() is 0" << endl;
    cout << "the number of hits in shower spine is " << shower_spine_hit_list.size() << endl;
    if(shower_spine_hit_list.size()<50){
      cout << "Requiring significant number of final hits" << endl;
      shower_spine_hit_list.clear();
      shower_spine_hit_list_vector.push_back(shower_spine_hit_list);
      continue;
    }


    shower_spine_hit_list_vector.push_back(shower_spine_hit_list);

    }
  }




  //----------------------------------------------------------------------------------------------------------------


  bool CCKaonProducer::collectSubsectionHits(const lar_content::ThreeDSlidingFitResult &extrapolated_fit,
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
					     double true_length)
  {
    
    float extrapolated_start_l = 0.;
    float extrapolated_start_t = 0.;
    float extrapolated_start_t2 = 0.;

    float extrapolated_end_l = 0.;
    float extrapolated_end_t = 0.;
    float extrapolated_end_t2 = 0.;
    
    float hit_l = 0.;
    float hit_t = 0.;
    float hit_t2 = 0.;


    pandora::CartesianVector pandora_extrapolated_start_position(extrapolated_start_position.X(), extrapolated_start_position.Y(), extrapolated_start_position.Z());
    pandora::CartesianVector pandora_extrapolated_end_position(extrapolated_end_position.X(), extrapolated_end_position.Y(), extrapolated_end_position.Z());

    extrapolated_fit.GetLocalPosition(pandora_extrapolated_start_position, extrapolated_start_l, extrapolated_start_t, extrapolated_start_t2);
    extrapolated_fit.GetLocalPosition(pandora_extrapolated_end_position, extrapolated_end_l, extrapolated_end_t, extrapolated_end_t2);

    cout << "extrapolated_start_l: " << extrapolated_start_l << endl;;
    cout << "extrapolated_end_l: " << extrapolated_end_l << endl;;

    std::vector<art::Ptr<recob::Hit>> collected_hits;
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;
    pandora::CartesianPointVector pandora_hit_positions;

    for(size_t i_h=0; i_h<hits_from_pfparticle.size(); i_h++){

      spacepoint_vec.clear();

      if(std::find(shower_spine_hit_list.begin(), shower_spine_hit_list.end(), hits_from_pfparticle[i_h]) != shower_spine_hit_list.end()) continue;
      if(std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), hits_from_pfparticle[i_h]) != unavailable_hit_list.end()) continue;

      //particles_per_hit.get(hits_from_pfparticle[i_h].key(),particle_vec_cheat,match_vec_cheat);
      spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle[i_h].key());
      if(spacepoint_vec.size()!=1) continue;

      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
      pandora_hit_positions.push_back(pandora_hit_position);

      extrapolated_fit.GetLocalPosition(pandora_hit_position, hit_l, hit_t, hit_t2);

      //cout << "is_end_downstream: " << is_end_downstream << ", hit_l: " << hit_l << endl;

      // Assess whether hit is within section boundaries
      if(is_end_downstream && ((hit_l < extrapolated_start_l) || (hit_l > extrapolated_end_l))) continue;
      if(!is_end_downstream && ((hit_l > extrapolated_start_l) || (hit_l < extrapolated_end_l))) continue;


      // Assess whether hit is close to connecting line - taking account hit width if necessary
      cout << "Assess whether hit is close to connecting line" << endl;
      if (this->isCloseToLine(hit_position, extrapolated_start_position, extrapolated_direction, distance_to_line)){
	collected_hits.push_back(hits_from_pfparticle[i_h]);
      }
      else{
	const TVector3 closest_point_in_hit(getClosestPointToLine3D(extrapolated_start_position, extrapolated_direction, hits_from_pfparticle[i_h], hit_position));

	if (this->isCloseToLine(closest_point_in_hit, extrapolated_start_position, extrapolated_direction, distance_to_line)) collected_hits.push_back(hits_from_pfparticle[i_h]);
      }

    }

      TVector3 hit_position;
      TVector3 distance;

    for(unsigned int i = 0; i < pandora_hit_positions.size(); i++){

      hit_position.SetXYZ(pandora_hit_positions[i].GetX(), pandora_hit_positions[i].GetY(), pandora_hit_positions[i].GetZ());
      distance = Kend_candidate - hit_position;

      if(distance.Mag() > true_length){

	//cout << "distance from K end is exceeding true pi+ track length" << endl;
	extrapolated_fit.GetLocalPosition(pandora_hit_positions[i], hit_l, hit_t, hit_t2);
	//cout << "is_end_downstream: " << is_end_downstream << ", hit_l: " << hit_l << endl;

	if( ( is_end_downstream && ((hit_l < extrapolated_start_l) || (hit_l > extrapolated_end_l)) ) ||
	    ( !is_end_downstream && ((hit_l > extrapolated_start_l) || (hit_l < extrapolated_end_l)) )
	    ){
	  //cout << "this hit would be excluded by continue" << endl;
	}else cout << "INCLUDED" << endl;

      }else cout << "NOT EXCEEDING" << endl;

    }
    
    const int n_initial_hits(shower_spine_hit_list.size());
    this->collectConnectedHits(collected_hits, spacepoint_per_hit, extrapolated_start_position, extrapolated_direction, running_fit_position_vec, pandora_running_fit_position_vec, shower_spine_hit_list);
    const int n_final_hits(shower_spine_hit_list.size());

    cout << "n_initial_hits: " << n_initial_hits << ", n_final_hits: " << n_final_hits << endl;

    return (n_final_hits != n_initial_hits);
  }



  //------------------------------------------------------------------------------------------------------------------------------------------

  

  bool CCKaonProducer::collectSubsectionHits(const lar_content::ThreeDSlidingFitResult &extrapolated_fit,
					     const TVector3 &extrapolated_start_position,
					     const TVector3 &extrapolated_end_position,
					     const TVector3 &extrapolated_direction,
					     const bool is_end_downstream,
					     std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
					     art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
					     vector<TVector3> &running_fit_position_vec,
					     pandora::CartesianPointVector &pandora_running_fit_position_vec,
					     std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
					     std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list)
  {
    
    float extrapolated_start_l = 0.;
    float extrapolated_start_t = 0.;
    float extrapolated_start_t2 = 0.;

    pandora::CartesianVector pandora_extrapolated_start_position(extrapolated_start_position.X(), extrapolated_start_position.Y(), extrapolated_start_position.Z());
    pandora::CartesianVector pandora_extrapolated_end_position(extrapolated_end_position.X(), extrapolated_end_position.Y(), extrapolated_end_position.Z());

    extrapolated_fit.GetLocalPosition(pandora_extrapolated_start_position, extrapolated_start_l, extrapolated_start_t, extrapolated_start_t2);

    float extrapolated_end_l = 0.;
    float extrapolated_end_t = 0.;
    float extrapolated_end_t2 = 0.;
    extrapolated_fit.GetLocalPosition(pandora_extrapolated_end_position, extrapolated_end_l, extrapolated_end_t, extrapolated_end_t2);

    cout << "extrapolated_start_l: " << extrapolated_start_l << endl;;
    cout << "extrapolated_end_l: " << extrapolated_end_l << endl;;

    std::vector<art::Ptr<recob::Hit>> collected_hits;
    //std::vector<simb::MCParticle const*> particle_vec_cheat;
    //std::vector<anab::BackTrackerHitMatchingData const*> match_vec_cheat;
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;

    pandora::CartesianPointVector pandora_hits_from_pfparticle_position_vec;
    for(size_t i_h=0; i_h<hits_from_pfparticle.size(); i_h++){
      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle[i_h].key());
      if(spacepoint_vec.size()==0) continue;
      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      pandora::CartesianVector pandora_hit_from_pfparticle_position(hit_position.X(), hit_position.Y(), hit_position.Z());
      pandora_hits_from_pfparticle_position_vec.push_back(pandora_hit_from_pfparticle_position);
    }
    lar_content::LArConnectionPathwayHelper::SortByDistanceToPoint smpl(pandora_extrapolated_start_position);
    std::sort(pandora_hits_from_pfparticle_position_vec.begin(), pandora_hits_from_pfparticle_position_vec.end(), smpl);

    //for(auto ptr : pandora_hits_from_pfparticle_position_vec)
    std::vector<art::Ptr<recob::Hit>> hits_from_pfparticle_sort;
    for(size_t i_p=0; i_p<pandora_hits_from_pfparticle_position_vec.size(); i_p++){
      TVector3 sort_position;
      sort_position.SetXYZ(pandora_hits_from_pfparticle_position_vec[i_p].GetX(), pandora_hits_from_pfparticle_position_vec[i_p].GetY(), pandora_hits_from_pfparticle_position_vec[i_p].GetZ());

      for(size_t i_h=0; i_h<hits_from_pfparticle.size(); i_h++){
	spacepoint_vec.clear();
	spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle[i_h].key());
	const TVector3 hit_position = spacepoint_vec[0]->XYZ();
	if(sort_position == hit_position){
	  hits_from_pfparticle_sort.push_back(hits_from_pfparticle[i_h]);
	  break;
	}
      }
    }
    
    
    for(size_t i_h=0; i_h<hits_from_pfparticle_sort.size(); i_h++){

      spacepoint_vec.clear();

      if(std::find(shower_spine_hit_list.begin(), shower_spine_hit_list.end(), hits_from_pfparticle_sort[i_h]) != shower_spine_hit_list.end()) continue;
      if(std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), hits_from_pfparticle_sort[i_h]) != unavailable_hit_list.end()) continue;

      //particles_per_hit.get(hits_from_pfparticle[i_h].key(),particle_vec_cheat,match_vec_cheat);
      spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle_sort[i_h].key());
      if(spacepoint_vec.size()!=1) continue;
      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());

      float hit_l = 0.;
      float hit_t = 0.;
      float hit_t2 = 0.;
      extrapolated_fit.GetLocalPosition(pandora_hit_position, hit_l, hit_t, hit_t2);

      //cout << "is_end_downstream: " << is_end_downstream << ", hit_l: " << hit_l << endl;

      // Assess whether hit is within section boundaries

      if(is_end_downstream && ((hit_l < extrapolated_start_l) || (hit_l > extrapolated_end_l))) continue;
      if(!is_end_downstream && ((hit_l > extrapolated_start_l) || (hit_l < extrapolated_end_l))) continue;


      // Assess whether hit is close to connecting line - taking account hit width if necessary
      cout << "Assess whether hit is close to connecting line" << this->isCloseToLine(hit_position, extrapolated_start_position, extrapolated_direction, distance_to_line)  << endl;

      this->isCloseToLine(hit_position, extrapolated_start_position, extrapolated_direction, distance_to_line);

      if (this->isCloseToLine(hit_position, extrapolated_start_position, extrapolated_direction, distance_to_line)){
	cout << "called collected_hits.push_back" << endl;
	collected_hits.push_back(hits_from_pfparticle_sort[i_h]);
      }
      /*
      else{
	const TVector3 closest_point_in_hit(getClosestPointToLine3D(extrapolated_start_position, extrapolated_direction, hits_from_pfparticle_sort[i_h], hit_position));

	if (this->isCloseToLine(closest_point_in_hit, extrapolated_start_position, extrapolated_direction, distance_to_line)) collected_hits.push_back(hits_from_pfparticle_sort[i_h]);
      }
      */

    }

    /*
    for(size_t i_h=0; i_h<hits_from_pfparticle.size(); i_h++){

      spacepoint_vec.clear();

      if(std::find(shower_spine_hit_list.begin(), shower_spine_hit_list.end(), hits_from_pfparticle[i_h]) != shower_spine_hit_list.end()) continue;
      if(std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), hits_from_pfparticle[i_h]) != unavailable_hit_list.end()) continue;

      //particles_per_hit.get(hits_from_pfparticle[i_h].key(),particle_vec_cheat,match_vec_cheat);
      spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle[i_h].key());
      if(spacepoint_vec.size()!=1) continue;
      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());

      float hit_l = 0.;
      float hit_t = 0.;
      float hit_t2 = 0.;
      extrapolated_fit.GetLocalPosition(pandora_hit_position, hit_l, hit_t, hit_t2);

      //cout << "is_end_downstream: " << is_end_downstream << ", hit_l: " << hit_l << endl;

      // Assess whether hit is within section boundaries
      if(is_end_downstream && ((hit_l < extrapolated_start_l) || (hit_l > extrapolated_end_l))) continue;
      if(!is_end_downstream && ((hit_l > extrapolated_start_l) || (hit_l < extrapolated_end_l))) continue;


      // Assess whether hit is close to connecting line - taking account hit width if necessary
      cout << "Assess whether hit is close to connecting line" << endl;
      if (this->isCloseToLine(hit_position, extrapolated_start_position, extrapolated_direction, distance_to_line)){
	cout << "called collected_hits.push_back" << endl;
	collected_hits.push_back(hits_from_pfparticle[i_h]);
      }
      else{
	
	const TVector3 closest_point_in_hit(getClosestPointToLine3D(extrapolated_start_position, extrapolated_direction, hits_from_pfparticle[i_h], hit_position));

	if (this->isCloseToLine(closest_point_in_hit, extrapolated_start_position, extrapolated_direction, distance_to_line)) collected_hits.push_back(hits_from_pfparticle[i_h]);
      }

    }
    */

    const int n_initial_hits(shower_spine_hit_list.size());
    this->collectConnectedHits(collected_hits, spacepoint_per_hit, extrapolated_start_position, extrapolated_direction, running_fit_position_vec, pandora_running_fit_position_vec, shower_spine_hit_list);
    const int n_final_hits(shower_spine_hit_list.size());

    return (n_final_hits != n_initial_hits);
  }

  

  //------------------------------------------------------------------------------------------------------------------------------------------


  void CCKaonProducer::collectConnectedHits(std::vector<art::Ptr<recob::Hit>>& collected_hits,
					    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
					    const TVector3 &extrapolated_start_position,
					    const TVector3 &extrapolated_direction,
					    vector<TVector3> &running_fit_position_vec,
					    pandora::CartesianPointVector &pandora_running_fit_position_vec,
					    std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list)
  {
    // Now add connected hits - taking into account hit width if necessary

    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;
    bool found = true;

    while(found)
      {
	found = false;

	vector<TVector3> shower_spine_hit_list_position;
	std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_2;

	for(size_t i_h=0; i_h<shower_spine_hit_list.size(); i_h++){
	  spacepoint_vec_2 = spacepoint_per_hit.at(shower_spine_hit_list[i_h].key());
	  const TVector3 hit2_position = spacepoint_vec_2[0]->XYZ();
	  shower_spine_hit_list_position.push_back(hit2_position);
	}


	for(size_t i_h=0; i_h<collected_hits.size(); i_h++){


	  if(std::find(shower_spine_hit_list.begin(), shower_spine_hit_list.end(), collected_hits[i_h]) != shower_spine_hit_list.end()) continue;

	  spacepoint_vec.clear();
	  spacepoint_vec = spacepoint_per_hit.at(collected_hits[i_h].key());
	  TVector3 hit_position = spacepoint_vec[0]->XYZ();
	  pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z()); 

	  if(this->getClosestDistance(hit_position, running_fit_position_vec) > hit_connection_distance){


	    //pandora::CartesianVector pandora_extrapolated_start_position(extrapolated_start_position.X(), extrapolated_start_position.Y(), extrapolated_start_position.Z());
	    //pandora::CartesianVector pandora_extrapolated_direction(extrapolated_direction.X(), extrapolated_direction.Y(), extrapolated_direction.Z());
	    //pandora::CaloHit pandora_collected_hit(collected_hits[i_h]);

	    //hit_position = getClosestPointToLine3D(extrapolated_start_position, extrapolated_direction, collected_hits[i_h], hit_position);

	    //if(getClosestDistance(collected_hits[i_h], shower_spine_hit_list, spacepoint_per_hit) > hit_connection_distance) continue;
	    //if(this->getClosestDistance(hit_position, shower_spine_hit_list_position) > hit_connection_distance) continue;
	    continue;
	  }

	  found = true;
	  running_fit_position_vec.push_back(hit_position);
	  pandora_running_fit_position_vec.push_back(pandora_hit_position);
	  shower_spine_hit_list.push_back(collected_hits[i_h]);
	}
      }

  }


  //------------------------------------------------------------------------------------------------------------------------------------------     

  TVector3 CCKaonProducer::getClosestPointToLine3D(const TVector3 &extrapolated_start_position,
						   const TVector3 &extrapolated_direction,
						   art::Ptr<recob::Hit>& collected_hit,
						   const TVector3 &hit_position)
  {

    if (std::fabs(extrapolated_direction.Z()) < std::numeric_limits<float>::epsilon()) return hit_position;

    double x_on_line = extrapolated_start_position.X();
    double y_on_line = extrapolated_start_position.Y();

    if (std::fabs(extrapolated_direction.X()) > std::numeric_limits<float>::epsilon()){
      const double gradient_x = extrapolated_direction.Z() / extrapolated_direction.X();
      x_on_line += ((hit_position.Z() - extrapolated_start_position.Z()) / gradient_x);
    }

    if (std::fabs(extrapolated_direction.Y()) > std::numeric_limits<float>::epsilon()){
      const double gradient_y = extrapolated_direction.Z() / extrapolated_direction.Y();
      y_on_line += ((hit_position.Z() - extrapolated_start_position.Z()) / gradient_y);
    }

    double hit_startT = collected_hit->PeakTimeMinusRMS();;
    double hit_endT = collected_hit->PeakTimePlusRMS();;
    const double hit_width = hit_startT - hit_endT;

    const double hit_low_x_edge = hit_position.X() - hit_width * 0.5;
    const double hit_high_x_edge = hit_position.X() + hit_width * 0.5;
    double closest_point_x;

    if(x_on_line < hit_low_x_edge) closest_point_x = hit_low_x_edge;
    else if(x_on_line < hit_high_x_edge) closest_point_x = hit_high_x_edge;
    else closest_point_x = x_on_line;

    const double hit_low_y_edge = hit_position.Y() - hit_width * 0.5;
    const double hit_high_y_edge = hit_position.Y() + hit_width * 0.5;
    double closest_point_y;

    if(y_on_line < hit_low_y_edge) closest_point_y = hit_low_y_edge;
    else if(y_on_line < hit_high_y_edge) closest_point_y = hit_high_y_edge;
    else closest_point_y = y_on_line;

    return TVector3(closest_point_x, closest_point_y, hit_position.Z());

  }


  //------------------------------------------------------------------------------------------------------------------------------------------


  bool CCKaonProducer::isCloseToLine(const TVector3 &hit_position,
				     const TVector3 &line_start,
				     const TVector3 &line_direction,
				     const double distance_to_line)
  {
    const double transverse_distance_from_line = line_direction.Cross(hit_position - line_start).Mag();
    if(transverse_distance_from_line > distance_to_line) return false;
    cout << "transverse_distance_from_line is " << transverse_distance_from_line << ", and distance_to_line is " << distance_to_line << endl;
    return true;
  }



  //------------------------------------------------------------------------------------------------------------------------------------------  

  double CCKaonProducer::getClosestDistance(art::Ptr<recob::Hit>& collected_hit,
					    art::Ptr<recob::Hit>& shower_spine_hit,
					    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit)
  {
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_1;
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_2;

    spacepoint_vec_1 = spacepoint_per_hit.at(collected_hit.key());
    spacepoint_vec_2 = spacepoint_per_hit.at(shower_spine_hit.key());

    const TVector3 hit1_position = spacepoint_vec_1[0]->XYZ();
    const TVector3 hit2_position = spacepoint_vec_2[0]->XYZ();

    double hit1_startT = collected_hit->PeakTimeMinusRMS();;
    double hit1_endT = collected_hit->PeakTimePlusRMS();;
    const double hit1_width = hit1_startT - hit1_endT;
    //cout << "hit1_width: " <<  hit1_width << endl;

    const double hit1_low_x_edge = hit1_position.X() - hit1_width * 0.5;
    const double hit1_high_x_edge = hit1_position.X() + hit1_width * 0.5;
    const double hit1_low_y_edge = hit1_position.Y() - hit1_width * 0.5;
    const double hit1_high_y_edge = hit1_position.Y() + hit1_width * 0.5;

    double hit2_startT = shower_spine_hit->PeakTimeMinusRMS();;
    double hit2_endT = shower_spine_hit->PeakTimePlusRMS();;
    const double hit2_width = hit2_startT - hit2_endT;

    const double hit2_low_x_edge = hit2_position.X() - hit2_width * 0.5;
    const double hit2_high_x_edge = hit2_position.X() + hit2_width * 0.5;
    const double hit2_low_y_edge = hit2_position.Y() - hit2_width * 0.5;
    const double hit2_high_y_edge = hit2_position.Y() + hit2_width * 0.5;

    const double mod_delta_z = std::fabs(hit1_position.Z() - hit2_position.Z());
    double delta_x, delta_y;

    if(!((hit1_high_x_edge < hit2_low_x_edge) || (hit2_high_x_edge < hit1_low_x_edge) || (hit1_high_y_edge < hit2_low_y_edge) || (hit2_high_y_edge < hit1_low_y_edge))) return mod_delta_z;

    if(hit1_low_x_edge < hit2_low_x_edge) delta_x = hit2_low_x_edge - hit1_high_x_edge;
    else delta_x = hit1_low_x_edge - hit2_high_x_edge;

    if(hit1_low_y_edge < hit2_low_y_edge) delta_y = hit2_low_y_edge - hit1_high_y_edge;
    else delta_y = hit1_low_y_edge - hit2_high_y_edge;

    //cout << "connectedness: " << std::sqrt((delta_x*delta_x) + (delta_y*delta_y) + (mod_delta_z*mod_delta_z)) << endl;
    return std::sqrt((delta_x*delta_x) + (delta_y*delta_y) + (mod_delta_z*mod_delta_z));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------  

  double CCKaonProducer::getClosestDistance(art::Ptr<recob::Hit>& collected_hit,
					    std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
					    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit)
  {
    double closest_distance(std::numeric_limits<double>::max());

    for(size_t i_h=0; i_h<shower_spine_hit_list.size(); i_h++){
      const double separation = getClosestDistance(collected_hit, shower_spine_hit_list[i_h], spacepoint_per_hit);
      if(separation < closest_distance) closest_distance = separation;
    }
    return closest_distance;
  }
    
  //------------------------------------------------------------------------------------------------------------------------------------------


  double CCKaonProducer::getClosestDistance(const TVector3 &position,
					    vector<TVector3> &test_positions)
  {
    double closest_distance_sqaured(std::numeric_limits<double>::max());

    for(const TVector3 test_position : test_positions){
      const double separation_squared = (test_position - position).Mag2();
      if (separation_squared < closest_distance_sqaured) closest_distance_sqaured = separation_squared;
    }

    return std::sqrt(closest_distance_sqaured);
  }


  //------------------------------------------------------------------------------------------------------------------------------------------

  void CCKaonProducer::obtainLength(std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
				    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
				    TVector3 Kend_candidate)
  {

    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;
    pandora::CartesianPointVector pandora_hit_positions;
    //vector<TVector3> hit_positions;

    for(size_t i_h=0; i_h<shower_spine_hit_list.size(); i_h++){

      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(shower_spine_hit_list[i_h].key());
      if(spacepoint_vec.size()!=1) continue;

      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
      pandora_hit_positions.push_back(pandora_hit_position);
      //hit_positions.push_back(hit_position);
    }

    pandora::CartesianVector pandora_start_position(Kend_candidate.X(), Kend_candidate.Y(), Kend_candidate.Z());

    lar_content::LArConnectionPathwayHelper::SortByDistanceToPoint smpl(pandora_start_position);
    std::sort(pandora_hit_positions.begin(), pandora_hit_positions.end(), smpl);

    double length = 0;
    TVector3 hit_position;
    TVector3 hit_position_prev;
    TVector3 distance;
    TVector3 displacemnet;
    
    for(unsigned int i = 0; i < pandora_hit_positions.size(); i++){

      hit_position.SetXYZ(pandora_hit_positions[i].GetX(), pandora_hit_positions[i].GetY(), pandora_hit_positions[i].GetZ());
      distance = Kend_candidate - hit_position;

      cout << "distance from K end is " << distance.Mag() << pandora_hit_positions[i].GetX() << " " << pandora_hit_positions[i].GetY() << " " << pandora_hit_positions[i].GetZ() << endl;

      if(i>0){
	hit_position_prev.SetXYZ(pandora_hit_positions[i-1].GetX(), pandora_hit_positions[i-1].GetY(), pandora_hit_positions[i-1].GetZ());
	displacemnet = hit_position - hit_position_prev;
	length += displacemnet.Mag();
      }
    }

    cout << "rebuilt length is " << length << endl;

  }
  

  //------------------------------------------------------------------------------------------------------------------------------------------





  //------------------------------------------------------------------------------------------------------------------------------------------

  void CCKaonProducer::obtainLongitudinalDecomposition(std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
						       art::FindManyP<recob::SpacePoint>& spacepoint_per_hit)
  {

    //cout << "GetWireZPitch "  << lar_content::LArGeometryHelper::GetWireZPitch(this->GetPandora()) << endl;

    const float sliding_fit_pitch = 0.30;
    //const float sliding_fit_pitch = 0.030;
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;
    pandora::CartesianPointVector pandora_hit_positions;

    cout << "inside obtainLongitudinalDecomposition shower_spine_hit_list.size() is" << shower_spine_hit_list.size() << endl;

    for(size_t i_h=0; i_h<shower_spine_hit_list.size(); i_h++){

      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(shower_spine_hit_list[i_h].key());
      if(spacepoint_vec.size()!=1) continue;

      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
      pandora_hit_positions.push_back(pandora_hit_position);

    }

    cout << "pandora_hit_positions.size() is " << pandora_hit_positions.size() << endl;

    const lar_content::ThreeDSlidingFitResult threeD_sliding_fit(&pandora_hit_positions, micro_sliding_fit_window, sliding_fit_pitch);
    std::vector<lar_content::LayerFitResultMap> layer_fit_result_map(threeD_sliding_fit.GetLayerFitResultMap());
    ////std::map<std::map<int,int> pandora::CaloHitList> layer_to_hit_map;
    ////std::map<int, std::map<int, pandora::CaloHitList>> layer_to_hit_map;
    std::map< std::vector<int>,  art::Ptr<recob::Hit>> layer_to_hit_map;
    
    for(size_t i_h=0; i_h<shower_spine_hit_list.size(); i_h++){

      float hit_l = 0; float hit_t1 = 0; float hit_t2 = 0;
      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(shower_spine_hit_list[i_h].key()); 
      if(spacepoint_vec.size()!=1) continue;
      const TVector3 hit_position = spacepoint_vec[0]->XYZ();
      const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z()); 

      //cout << "hit_position is " << hit_position.X() << " " << hit_position.Y() << " " << hit_position.Z() << endl;
      //cout << "pandora_hit_position is " << pandora_hit_position.GetX() << " " << pandora_hit_position.GetY() << " " << pandora_hit_position.GetZ() << endl;
      //int layer_first; int layer_second;

      threeD_sliding_fit.GetLocalPosition(pandora_hit_position, hit_l, hit_t1, hit_t2);
      std::vector<int> layers = {threeD_sliding_fit.GetLayers(hit_l)[0], threeD_sliding_fit.GetLayers(hit_l)[1]};
      layer_to_hit_map[layers] = shower_spine_hit_list[i_h];
    }

    double running_distance = 0;
    std::map< std::vector<int>, art::Ptr<recob::Hit>>::iterator iter;
    
    for(auto iter = layer_to_hit_map.begin(); iter != layer_to_hit_map.end(); iter++){
    
    std::vector<int> layer = iter->first;
    //cout << "layer = iter->first loop: " << layer[0] << " " << layer[1] << endl;
    //double layer_l_first = layer_fit_result_map[0].at(layer[0]).GetL();
    //double layer_l_second = layer_fit_result_map[1].at(layer[1]).GetL();
    
    std::vector<int> higher_layer;
    std::vector<int> middle_layer;
    std::vector<int> lower_layer;
    
    //if(std::next(iter) == layer_to_hit_map.end()) higher_layer = layer;
    //else ;
    
    std::next(iter) == layer_to_hit_map.end() ? higher_layer = layer : higher_layer = std::next(iter)->first;
    middle_layer = layer;
    iter == layer_to_hit_map.begin() ? lower_layer = layer : lower_layer = std::prev(iter)->first;
    
    pandora::CartesianVector lower_layer_position(0., 0., 0.), middle_layer_position(0., 0., 0.), higher_layer_position(0., 0., 0.);
    
    threeD_sliding_fit.GetGlobalPosition(layer_fit_result_map[0].at(lower_layer[0]).GetL(), layer_fit_result_map[0].at(lower_layer[0]).GetFitT(), layer_fit_result_map[1].at(lower_layer[1]).GetFitT(), lower_layer_position);
    threeD_sliding_fit.GetGlobalPosition(layer_fit_result_map[0].at(middle_layer[0]).GetL(), layer_fit_result_map[0].at(middle_layer[0]).GetFitT(), layer_fit_result_map[1].at(middle_layer[1]).GetFitT(), middle_layer_position);
    threeD_sliding_fit.GetGlobalPosition(layer_fit_result_map[0].at(higher_layer[0]).GetL(), layer_fit_result_map[0].at(higher_layer[0]).GetFitT(), layer_fit_result_map[1].at(higher_layer[1]).GetFitT(), higher_layer_position);

    /*    
    cout << "layer_fit_result_map[0].at(lower_layer[0]).GetL(): " << layer_fit_result_map[0].at(lower_layer[0]).GetL() << endl;
    cout << "layer_fit_result_map[0].at(lower_layer[0]).GetT(): " << layer_fit_result_map[0].at(lower_layer[0]).GetFitT() << endl;
    cout << "layer_fit_result_map[1].at(lower_layer[1]).GetT(): " << layer_fit_result_map[1].at(lower_layer[1]).GetFitT() << endl;
    cout << " " << endl;
    cout << "layer_fit_result_map[0].at(middle_layer[0]).GetL(): " << layer_fit_result_map[0].at(middle_layer[0]).GetL() << endl;
    cout << "layer_fit_result_map[0].at(middle_layer[0]).GetT(): " << layer_fit_result_map[0].at(middle_layer[0]).GetFitT() << endl;
    cout << "layer_fit_result_map[1].at(middle_layer[1]).GetT(): " << layer_fit_result_map[1].at(middle_layer[1]).GetFitT() << endl;
    */

    double layer_length;
    std::next(iter) == layer_to_hit_map.end() ? layer_length = 0. : iter == layer_to_hit_map.begin() ? layer_length = 0. : layer_length = (middle_layer_position - lower_layer_position).GetMagnitude();
    
    running_distance += layer_length;
    
    }

    cout << "running_distance: " << running_distance << endl;
    
  }


  //------------------------------------------------------------------------------------------------------------------------------------------



  int CCKaonProducer::fillHistAngularDistributionMapCheat( std::map<int, std::map<int, double>> &angular_distribution_mrgid_map,
							   std::map<int,int> &mrgidpdg,
							   std::vector<TH1D*> &h_angular_distribution_pfparticle_cheat,
							   vector<int>& v_pdg
							   ){

    //TH1D * h = new TH1D("", "", num_bin, -3.14, 3.14);

    map<int, int>::iterator it_mrgidpdg;
    it_mrgidpdg = mrgidpdg.begin();
 
    cout << "CHEATING angular_distribution_mrgid_map.size(): " << angular_distribution_mrgid_map.size() << endl;

    if(angular_distribution_mrgid_map.size()==0) return 0;

    for (auto const& angular_distribution_map : angular_distribution_mrgid_map) {//for each mrgid
      
      TH1D * h = new TH1D("", "", num_bin, -3.14, 3.14); 
      for (auto const& x : angular_distribution_map.second) {
	h->Fill(theta0XZBinSize * x.first, x.second);
      }
      h_angular_distribution_pfparticle_cheat.push_back(h);
      v_pdg.push_back(mrgidpdg[it_mrgidpdg->first]);
      cout << "mrgidpdg[it_mrgidpdg->first]: " << mrgidpdg[it_mrgidpdg->first] << endl;
      cout << "angular_distribution_map.second hist entries: " << h->GetEntries() << endl;
      //h->Reset();
      it_mrgidpdg++;
    }

    cout << "h_angular_distribution_pfparticle_cheat.size(): " << h_angular_distribution_pfparticle_cheat.size() <<endl;
    cout << "v_pdg.size(): " << v_pdg.size() << endl;
    for (int i: v_pdg){
      cout << "content of v_pdg is " << i << endl;
    }

    return 0;
  }


  //------------------------------------------------------------------------------------------------------------------------------------------



  int CCKaonProducer::fillHistAngularDistributionMapCheat3D( std::map<int, std::map<int, std::map<int, double>>> &angular_distribution_mrgid_map_3D,
							     std::map<int,int> &mrgidpdg,
							     std::vector<TH2D*> &h_angular_distribution_pfparticle_cheat_3D,
							     vector<int>& v_pdg
							     ){

    map<int, int>::iterator it_mrgidpdg;
    it_mrgidpdg = mrgidpdg.begin();
 
    if(angular_distribution_mrgid_map_3D.size()==0) return 0;

    for (auto const& angular_distribution_map_3D : angular_distribution_mrgid_map_3D) {//for each mrgid
      
    TH2D * h = new TH2D("", "", num_bin_theta, 0, 3.14, num_bin_phi, -3.14, 3.14); 
    for (auto const& x : angular_distribution_map_3D.second) {
      for (auto const& y : x.second) {
	h->Fill(thetaBinSize * x.first, phiBinSize * y.first, y.second);
      }
    }
    h_angular_distribution_pfparticle_cheat_3D.push_back(h);
    v_pdg.push_back(mrgidpdg[it_mrgidpdg->first]);
    it_mrgidpdg++;
    }

    return 0;
  }



  //------------------------------------------------------------------------------------------------------------------------------------------

  
  int CCKaonProducer::fillHistAngularDistributionMap( std::map<int, double> &angular_distribution_map,
						      std::vector<TH1D*> &h_angular_distribution_pfparticle,
						      vector<bool>& v_trk_flg,
						      bool trk_flg){

    //int num_bin = (int)(6.28/theta0XZBinSize);
    //TH1D * h = new TH1D("", "", num_bin, -3.14, 3.14);

    cout << "RECO angular_distribution_map.size(): " << angular_distribution_map.size() << endl;

    if(angular_distribution_map.size()==0) return 0;
    
    TH1D * h = new TH1D("", "", num_bin, -3.14, 3.14); 
    for (auto const& x : angular_distribution_map) {
      h->Fill(theta0XZBinSize * x.first, x.second);
    }
    h_angular_distribution_pfparticle.push_back(h);
    v_trk_flg.push_back(trk_flg);
    //h->Reset();
    
    return 0;
  }


  //------------------------------------------------------------------------------------------------------------------------------------------


  int CCKaonProducer::fillHistAngularDistributionMap3D( std::map<int, std::map<int, double>> &angular_distribution_map_3D,
							std::vector<TH2D*> &h_angular_distribution_pfparticle_3D,
							vector<bool>& v_trk_flg,
							bool trk_flg){

    //int num_bin_theta = (int)(6.28/thetaBinSize);
    //int num_bin_phi = (int)(6.28/phiBinSize);

    cout << "RECO angular_distribution_map_3D.size(): " << angular_distribution_map_3D.size() << endl;

    if(angular_distribution_map_3D.size()==0) return 0;
    
    TH2D * h = new TH2D("", "", num_bin_theta, 0, 3.14, num_bin_phi, -3.14, 3.14); 
    for (auto const& x : angular_distribution_map_3D) {
      for (auto const& y : x.second) {
	h->Fill(thetaBinSize * x.first, phiBinSize * y.first, y.second);
      }
    }
    h_angular_distribution_pfparticle_3D.push_back(h);
    v_trk_flg.push_back(trk_flg);
    //h->Reset();
    
    return 0;
  }



  //------------------------------------------------------------------------------------------------------------------------------------------


  int CCKaonProducer::fillHistAngularDistributionMapSurface( std::map<int, std::map<int, double>> &angular_distribution_map_3D,
							     std::vector<TH2D*> &h_angular_distribution_pfparticle_surface,
							     vector<bool>& v_trk_flg,
							     bool trk_flg){

    //int num_bin_theta = (int)(6.28/thetaBinSize);
    //int num_bin_phi = (int)(6.28/phiBinSize);
    //int num_bin_theta = 30;
    //int num_bin_phi = 30;

    cout << "RECO angular_distribution_map_3D.size(): " << angular_distribution_map_3D.size() << endl;

    if(angular_distribution_map_3D.size()==0) return 0;
    
    TH2D * h = new TH2D("", "", 30, 0, 3.14, 60, -3.14, 3.14); 
    for (auto const& x : angular_distribution_map_3D) {
      for (auto const& y : x.second) {
	h->Fill(thetaBinSize * x.first, phiBinSize * y.first, y.second);
      }
    }
    h_angular_distribution_pfparticle_surface.push_back(h);
    v_trk_flg.push_back(trk_flg);
    //h->Reset();
    
    return 0;
  }



  //------------------------------------------------------------------------------------------------------------------------------------------



  int CCKaonProducer::fillHistAngularDistributionMapSphere( std::map<int, std::map<int, double>> &angular_distribution_map_3D,
							    std::vector<TH3D*> &h_angular_distribution_pfparticle_sphere,
							    vector<bool>& v_trk_flg,
							    bool trk_flg){

    if(angular_distribution_map_3D.size()==0) return 0;
    
    //TH3D * h = new TH3D("", "", num_bin_theta, 0, 1.5, num_bin_theta, -3.14, 3.14, num_bin_phi, -3.14, 3.14); 
    TH3D * h = new TH3D("", "", 50, -1., 1., 50, -1., 1., 50, -1., 1.); 
    for (auto const& x : angular_distribution_map_3D) {
      for (auto const& y : x.second) {
	double theta = thetaBinSize * x.first;
	double phi = phiBinSize * y.first;
	double weight = y.second;
	h->Fill(TMath::Sin(theta)*TMath::Cos(phi), TMath::Sin(theta)*TMath::Sin(phi), TMath::Cos(theta), weight);
      }
    }
    h_angular_distribution_pfparticle_sphere.push_back(h);
    v_trk_flg.push_back(trk_flg);
    //h->Reset();
    
    return 0;
  }


  //------------------------------------------------------------------------------------------------------------------------------------------



  int CCKaonProducer::drawHistAngularDistributionMapCheat( std::vector<TH1D*> &h_angular_distribution_pfparticle_cheat,
							   std::vector<int> &v_pdg,
							   TString outfile_name,
							   TCanvas* &c){

    TFile outfile(outfile_name, "update");
    THStack * hs = new THStack("hs", "");  
    TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");

    std::vector<TH1D*> h_pdg(6);
    //int num_bin = (int)(6.28/theta0XZBinSize);
    for(int i_pdg=0; i_pdg<6; i_pdg++){
      h_pdg[i_pdg] = new TH1D("", "", num_bin, -3.14, 3.14);
    }

    c->SetFillStyle(4000);

    //cout << "h_angular_distribution_pfparticle_cheat.size(): " << h_angular_distribution_pfparticle_cheat.size() << endl;

    if(h_angular_distribution_pfparticle_cheat.size()==0){
      cout << "h_angular_distribution_pfparticle_cheat.size()==0" << endl;
      return 0;
    }

    for(unsigned int i = 0; i < h_angular_distribution_pfparticle_cheat.size(); i++){
      h_angular_distribution_pfparticle_cheat[i]->SetFillStyle(4000);

      if(v_pdg[i]==321){
	h_pdg[0]->SetLineColor(kBlue);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "K+");
	h_pdg[0]->Add(h_angular_distribution_pfparticle_cheat[i]);
      }
      else if(v_pdg[i]==-13){
	h_pdg[1]->SetLineColor(kCyan);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "mu+");
	h_pdg[1]->Add(h_angular_distribution_pfparticle_cheat[i]);
      }
      else if(v_pdg[i]==211){
	h_pdg[2]->SetLineColor(kMagenta);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "pi+");
	h_pdg[2]->Add(h_angular_distribution_pfparticle_cheat[i]);
      }
      else if(v_pdg[i]==999){
	//else if(v_pdg[i]==22 || v_pdg[i]==11 || v_pdg[i]==-11){
	h_pdg[3]->SetLineColor(kGreen+2);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "gamma/e");
	h_pdg[3]->Add(h_angular_distribution_pfparticle_cheat[i]);
      }
      else if(v_pdg[i]==2212 || v_pdg[i]==2112){
	h_pdg[4]->SetLineColor(kRed);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "nucleon");
	h_pdg[4]->Add(h_angular_distribution_pfparticle_cheat[i]);
      }
      else{
	h_pdg[5]->SetLineColor(kBlack);
	h_pdg[5]->Add(h_angular_distribution_pfparticle_cheat[i]); 
      }
      //cout << "adding stack histos: " << h_angular_distribution_pfparticle_cheat[i]->GetEntries() << endl;
    }

    leg->AddEntry(h_pdg[0], "K+");
    leg->AddEntry(h_pdg[1], "mu+");
    leg->AddEntry(h_pdg[2], "pi+");
    leg->AddEntry(h_pdg[3], "gamma/e");
    leg->AddEntry(h_pdg[4], "nucleon");
    leg->AddEntry(h_pdg[5], "others");

    for(int i_pdg=0; i_pdg<6; i_pdg++){
      hs->Add(h_pdg[i_pdg]);
    }


    hs->Draw("hist");
    //hs->Draw("hist nostack");
    leg->Draw();
    c->Write();

    cout << "cheat c written" << endl;

    //c->Close();

    delete hs;
    delete leg;

    return 0;
  }



  //------------------------------------------------------------------------------------------------------------------------------------------



  int CCKaonProducer::drawHistAngularDistributionMapCheat3D( std::vector<TH2D*> &h_angular_distribution_pfparticle_cheat_3D,
							     std::vector<int> &v_pdg,
							     TString outfile_name,
							     TCanvas* &c){
    
    TFile outfile(outfile_name, "update");
    THStack * hs = new THStack("hs", "");  
    TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");

    std::vector<TH2D*> h_pdg(6);
    //int num_bin_theta = (int)(6.28/thetaBinSize);
    //int num_bin_phi = (int)(6.28/phiBinSize);
    for(int i_pdg=0; i_pdg<6; i_pdg++){
      h_pdg[i_pdg] = new TH2D("", "", num_bin_theta, 0, 3.14, num_bin_phi, -3.14, 3.14);
    }

    c->SetFillStyle(1001);

    //cout << "h_angular_distribution_pfparticle_cheat.size(): " << h_angular_distribution_pfparticle_cheat.size() << endl;

    if(h_angular_distribution_pfparticle_cheat_3D.size()==0){
      cout << "h_angular_distribution_pfparticle_cheat.size()==0" << endl;
      return 0;
    }

    for(unsigned int i = 0; i < h_angular_distribution_pfparticle_cheat_3D.size(); i++){
      h_angular_distribution_pfparticle_cheat_3D[i]->SetFillStyle(1001);

      if(v_pdg[i]==321){
	h_pdg[0]->SetFillColor(kBlue);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "K+");
	h_pdg[0]->Add(h_angular_distribution_pfparticle_cheat_3D[i]);
      }
      else if(v_pdg[i]==-13){
	h_pdg[1]->SetFillColor(kCyan);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "mu+");
	h_pdg[1]->Add(h_angular_distribution_pfparticle_cheat_3D[i]);
      }
      else if(v_pdg[i]==211){
	h_pdg[2]->SetFillColor(kMagenta);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "pi+");
	h_pdg[2]->Add(h_angular_distribution_pfparticle_cheat_3D[i]);
      }
      else if(v_pdg[i]==999){
	//else if(v_pdg[i]==22 || v_pdg[i]==11 || v_pdg[i]==-11){
	h_pdg[3]->SetFillColor(kGreen+2);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "gamma/e");
	h_pdg[3]->Add(h_angular_distribution_pfparticle_cheat_3D[i]);
      }
      else if(v_pdg[i]==2212 || v_pdg[i]==2112){
	h_pdg[4]->SetFillColor(kRed);
	//leg->AddEntry(h_angular_distribution_pfparticle_cheat[i], "nucleon");
	h_pdg[4]->Add(h_angular_distribution_pfparticle_cheat_3D[i]);
      }
      else{
	h_pdg[5]->SetFillColor(kBlack);
	h_pdg[5]->Add(h_angular_distribution_pfparticle_cheat_3D[i]); 
      }
      //cout << "adding stack histos: " << h_angular_distribution_pfparticle_cheat[i]->GetEntries() << endl;
    }

    leg->AddEntry(h_pdg[0], "K+");
    leg->AddEntry(h_pdg[1], "mu+");
    leg->AddEntry(h_pdg[2], "pi+");
    leg->AddEntry(h_pdg[3], "gamma/e");
    leg->AddEntry(h_pdg[4], "nucleon");
    leg->AddEntry(h_pdg[5], "others");

    for(int i_pdg=0; i_pdg<6; i_pdg++){
      hs->Add(h_pdg[i_pdg]);
    }


    hs->Draw("hist""lego3 0");
    //hs->Draw("hist nostack");
    leg->Draw();
    c->Write();

    delete hs;
    delete leg;

    return 0;
  }


  //------------------------------------------------------------------------------------------------------------------------------------------


  int CCKaonProducer::drawHistAngularDistributionMap( std::vector<TH1D*> &h_angular_distribution_pfparticle,
						      std::vector<bool> &v_trk_flg,
						      TString outfile_name,
						      TCanvas* &c){

    TFile outfile(outfile_name, "update");
    THStack * hs = new THStack("hs", "");
    TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");

    //output_pdf.Form("reco_angle_distribution.pdf");
    c->SetFillStyle(4000);
    //gStyle->SetOptStat(0);

    //cout << "h_angular_distribution_pfparticle.size(): " << h_angular_distribution_pfparticle.size() << endl;

    if(h_angular_distribution_pfparticle.size()==0){
      cout << "h_angular_distribution_pfparticle.size()==0" << endl;
      return 0;
    }

    for(unsigned int i = 0; i < h_angular_distribution_pfparticle.size(); i++){
      h_angular_distribution_pfparticle[i]->SetFillStyle(4000);
      if(v_trk_flg[i]==true){
	h_angular_distribution_pfparticle[i]->SetLineColor(kRed);
	leg->AddEntry(h_angular_distribution_pfparticle[i], "reco track");
      }else{
	h_angular_distribution_pfparticle[i]->SetLineColor(kBlue);
	leg->AddEntry(h_angular_distribution_pfparticle[i], "reco shower");
      }
      cout << "adding stack histos: " << h_angular_distribution_pfparticle[i]->GetEntries() << endl;
      hs->Add(h_angular_distribution_pfparticle[i]);
    }

    hs->Draw("hist");
    //hs->Draw("hist nostack");
    leg->Draw();
    c->Write();

    cout << "reco c written" << endl;
    //c->Close();
    delete hs;
    delete leg;

    return 0;
  }



  //------------------------------------------------------------------------------------------------------------------------------------------


  int CCKaonProducer::drawHistAngularDistributionMap3D( std::vector<TH2D*> &h_angular_distribution_pfparticle_3D,
							std::vector<bool> &v_trk_flg,
							TString outfile_name,
							TCanvas* &c){

    TFile outfile(outfile_name, "update");
    THStack * hs = new THStack("hs", "");
    TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");

    //output_pdf.Form("reco_angle_distribution.pdf");
    c->SetFillStyle(4000);
    //gStyle->SetOptStat(0);

    //cout << "h_angular_distribution_pfparticle.size(): " << h_angular_distribution_pfparticle.size() << endl;

    if(h_angular_distribution_pfparticle_3D.size()==0){
      cout << "h_angular_distribution_pfparticle_3D.size()==0" << endl;
      return 0;
    }

    for(unsigned int i = 0; i < h_angular_distribution_pfparticle_3D.size(); i++){
      h_angular_distribution_pfparticle_3D[i]->SetFillStyle(1001);
      if(v_trk_flg[i]==true){
	h_angular_distribution_pfparticle_3D[i]->SetFillColor(kRed);
	leg->AddEntry(h_angular_distribution_pfparticle_3D[i], "reco track");
      }else{
	h_angular_distribution_pfparticle_3D[i]->SetFillColor(kBlue);
	leg->AddEntry(h_angular_distribution_pfparticle_3D[i], "reco shower");
      }
      cout << "adding stack histos: " << h_angular_distribution_pfparticle_3D[i]->GetEntries() << endl;
      hs->Add(h_angular_distribution_pfparticle_3D[i]);
    }

    hs->Draw("hist""lego3 0");
    //hs->Draw("hist nostack");
    leg->Draw();
    c->Write();

    delete hs;
    delete leg;

    return 0;
  }


  //------------------------------------------------------------------------------------------------------------------------------------------


  int CCKaonProducer::drawHistAngularDistributionMapSurface( std::vector<TH2D*> &h_angular_distribution_pfparticle_surface,
							     std::vector<bool> &v_trk_flg,
							     TString outfile_name,
							     TCanvas* &c){
    
    TFile outfile(outfile_name, "update");
    THStack * hs = new THStack("hs", "");
    TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");

    //output_pdf.Form("reco_angle_distribution.pdf");
    c->SetFillStyle(4000);
    //gStyle->SetOptStat(0);

    //cout << "h_angular_distribution_pfparticle.size(): " << h_angular_distribution_pfparticle.size() << endl;

    if(h_angular_distribution_pfparticle_surface.size()==0){
      cout << "h_angular_distribution_pfparticle_surface.size()==0" << endl;
      return 0;
    }

    for(unsigned int i = 0; i < h_angular_distribution_pfparticle_surface.size(); i++){
      h_angular_distribution_pfparticle_surface[i]->SetFillStyle(1001);
      if(v_trk_flg[i]==true){
	h_angular_distribution_pfparticle_surface[i]->SetFillColor(kRed);
	leg->AddEntry(h_angular_distribution_pfparticle_surface[i], "reco track");
      }else{
	h_angular_distribution_pfparticle_surface[i]->SetFillColor(kBlue);
	leg->AddEntry(h_angular_distribution_pfparticle_surface[i], "reco shower");
      }
      hs->Add(h_angular_distribution_pfparticle_surface[i]);
    }

    hs->Draw("hist""LEGO1 SPH");
    //hs->Draw("hist nostack");
    leg->Draw();
    c->Write();

    delete hs;
    delete leg;

    return 0;
  }



  //------------------------------------------------------------------------------------------------------------------------------------------



  int CCKaonProducer::drawHistAngularDistributionMapSphere( std::vector<TH3D*> &h_angular_distribution_pfparticle_sphere,
							    std::vector<bool> &v_trk_flg,
							    TString outfile_name,
							    TCanvas* &c){
    TFile outfile(outfile_name, "update");
    TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");

    std::vector<TH3D*> h_pdg(2);

    //int num_bin_theta = (int)(6.28/thetaBinSize);
    //int num_bin_phi = (int)(6.28/phiBinSize);
    for(int i_pdg=0; i_pdg<2; i_pdg++){
      h_pdg[i_pdg] = new TH3D("", "", 50, -1., 1., 50, -1., 1., 50, -1., 1.);
    }

    c->SetFillStyle(4000);

    if(h_angular_distribution_pfparticle_sphere.size()==0){
      cout << "h_angular_distribution_pfparticle_sphere.size()==0" << endl;
      return 0;
    }

    for(unsigned int i = 0; i < h_angular_distribution_pfparticle_sphere.size(); i++){
      h_angular_distribution_pfparticle_sphere[i]->SetFillStyle(1001);
      cout << "v_trk_flg[i] is " << v_trk_flg[i] << endl;
      if(v_trk_flg[i]==true){
	h_angular_distribution_pfparticle_sphere[i]->SetFillColor(kRed);
	h_pdg[0]->Add(h_angular_distribution_pfparticle_sphere[i]);
	//leg->AddEntry(h_angular_distribution_pfparticle_sphere[i], "reco track");
      }else{
	h_angular_distribution_pfparticle_sphere[i]->SetFillColor(kBlue);
	h_pdg[1]->Add(h_angular_distribution_pfparticle_sphere[i]);
	//leg->AddEntry(h_angular_distribution_pfparticle_sphere[i], "reco shower");
      }
      //hs->Add(h_angular_distribution_pfparticle_sphere[i]);
    }

    leg->AddEntry(h_pdg[0], "reco track");    
    leg->AddEntry(h_pdg[1], "reco shower");

    //hs->Draw("");
    //hs->Draw("hist nostack");
    h_pdg[0]->SetFillColor(kRed);
    h_pdg[1]->SetFillColor(kBlue);
    h_pdg[0]->Draw("BOX3");
    h_pdg[1]->Draw("BOX3""same");
    leg->Draw();
    c->Write();

    //delete hs;
    delete leg;

    return 0;
  }


}

#endif
