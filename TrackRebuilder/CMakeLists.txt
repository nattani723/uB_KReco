link_libraries( -L$ENV{BOOST_LIB} -lboost_system -lboost_filesystem)


art_make( 
	  BASENAME_ONLY
	  LIBRARY_NAME
		CCKaonHeaders
          LIB_LIBRARIES
		ubraw_RawData
		ubobj_Trigger
            	larcorealg_Geometry
            	larcore_Geometry_Geometry_service
		ubana_BasicTool_GeoAlgo
            	lardata_Utilities
            	larevt_Filters
            	lardataobj_RawData
		lardataobj_RecoBase
            	lardataobj_MCBase
		larcoreobj_SummaryData
		larsim_MCCheater_BackTrackerService_service
		larsim_MCCheater_ParticleInventoryService_service
		larpandora_LArPandoraInterface
            	${PANDORASDK}
            	${PANDORAMONITORING}
            	LArPandoraContent 
		nusimdata_SimulationBase
		${UB_DATA_TYPES}
		${Boost_SERIALIZATION_LIBRARY}
		nutools_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                art_Persistency_Common canvas
                art_Persistency_Provenance canvas
                art_Utilities canvas
                ${MF_MESSAGELOGGER}
                        
                ${FHICLCPP}
		cetlib cetlib_except
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
                ${ROOT_GEOM} 
          MODULE_LIBRARIES
			ubraw_RawData
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			ubana_BasicTool_GeoAlgo
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
			lardataobj_RecoBase
            		lardataobj_MCBase
            		larcoreobj_SummaryData	
			larsim_MCCheater_BackTrackerService_service
			larpandora_LArPandoraInterface
            		SignalShapingServiceMicroBooNE_service
			nusimdata_SimulationBase
			ubana_LEEPhotonAnalysis
			${UB_DATA_TYPES}
			${Boost_SERIALIZATION_LIBRARY}
                        nutools_RandomUtils_NuRandomService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
					cetlib cetlib_except
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
)
                      
install_headers()
install_fhicl()
install_source()
