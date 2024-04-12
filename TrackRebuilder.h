#include "ParticleDirectionFinder.h"
#include "TrackHitCollector.h"

#include "art/Framework/Core/EDProducer.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h" 

#include <iostream>

#ifndef TRACK_REBUILDER
#define TRACK_REBUILDER 1

using namespace pandora;

namespace kaon_reconstruction 
{
  class TrackRebuilder
  {
  public:
    TrackRebuilder();
    pandora::STATUSCODE Run(HitList& track_hit_list, const recob::Track& primary_track);

  private:
    recob::Track track_rebuild(HitList& track_hit_list, const recob::Track& primary_track) const;
    recob::Track build_track(int track_id, lar_content::LArTrackStateVector& track_state_vector) const;

  };


}//namespace kaon_reconstruction

#endif
