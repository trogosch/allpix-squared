/**
* @file
* @brief Definition of [CaloOutputWriter] module
* @copyright Copyright (c) 2017-2020 CERN and the Allpix Squared authors.
* This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
* In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
* Intergovernmental Organization or submit itself to any jurisdiction.
*
* Contains minimal dummy module to use as a start for the development of your own module
*
* Refer to the User's Manual for more details.
*/

#include <memory>

#include <map>
#include <string>
#include <vector>

#include <TH1I.h>
#include <TH2I.h>
#include <TH3I.h>
#include <TTree.h>

#include "core/config/Configuration.hpp"
#include "core/geometry/GeometryManager.hpp"
#include "core/messenger/Messenger.hpp"
#include "core/module/Module.hpp"

#include "Cluster.hpp"
#include "objects/PixelHit.hpp"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;
#endif

namespace allpix {
  /**
  * @ingroup Modules
  * @brief Module to do function
  *
  * More detailed explanation of module
  */
  class CaloOutputWriterModule : public Module {
  public:
    /**
    * @brief Constructor for this unique module
    * @param config Configuration object for this module as retrieved from the steering file
    * @param messenger Pointer to the messenger object to allow binding to messages on the bus
    * @param geo_manager Pointer to the geometry manager, containing the detectors
    */
    CaloOutputWriterModule(Configuration& config, Messenger* messenger, GeometryManager* geo_manager);
    /**
    * @brief Destructor deletes the internal objects used to build the ROOT Tree
    */
    // ~CaloOutputWriterModule() override;

    /**
    * @brief [Initialise this module]
    */
    void init() override;

    /**
    * @brief [Run the function of this module]
    */
    void run(unsigned int) override;

    /**
    * @brief Write the TTree and histograms to the output file
    */
    void finalize() override;

  private:
    // General module members
    GeometryManager* geo_manager_;
    Messenger* messenger_;
    std::vector<std::shared_ptr<PixelHitMessage>> pixel_messages_;
    std::vector<std::shared_ptr<MCParticleMessage>> mcparticle_messages_;
    std::vector<std::shared_ptr<MCTrackMessage>> mctrack_messages_;

    /**
     * @brief Perform a sparse clustering on the PixelHits
     */
    std::vector<Cluster> doClustering(int lay);


    std::vector<const MCParticle*> getPrimaryParticles(int lay) const;

    void addMCParticlesToList(int lay);

    void clearVectors();

    int getLayerNumber(std::string name);
    int getLane(int layer);
    int getLaneFORONECHIP(int layer, int idx);


    //Histogramms
    TH1D* h1PixelHits;
    TH1D* h1PixelHits2;
    TH1D* h1Clusters;
    TH1D* h1Clusters2;
    TH1D* h1ParticlesContributing2PixelHit;
    static const int Nhistlayer = 48;
    TH1D* h1ParticlesContributing2PixelHitLayer[Nhistlayer];
    TH1D* h1PixelHitsPerLayer;

    TH1D* h1MCParticles;
    TH1D* h1MCElectronsPositrons;
    TH1D* h1MCPhotons;
    TH1D* h1MCElectrons;
    TH1D* h1MCPositrons;
    TH1D* h1MCParticlesAngle2zAxis;
    TH1D* h1MCPrimariessAngle2zAxisLayer0;
    TH1D* h1MCParticleInLayerID;
    TH1D* h1MCParticlesPerLayer;

    TH3D* h3MCParticlesAngleEnergyLayer;
    TH3D* h3MCElectronAngleEnergyLayer;
    TH3D* h3MCPositronAngleEnergyLayer;
    TH3D* h3MCPhotonAngleEnergyLayer;

    TH1D* h1MCPrimaries;
    TH1D* h1MCPrimariesAngle2zAxis;
    TH1D* h1MCelectronAngle2zAxis;
    TH1D* h1MCpositronAngle2zAxis;
    TH1D* h1MCphotonAngle2zAxis;
    static const int NAnglebins = 3;
    TH1D* h1MCPrimariesEgyInAngleBin[NAnglebins];
    TH1D* h1MCelectronEgyInAngleBin[NAnglebins];
    TH1D* h1MCpositronEgyInAngleBin[NAnglebins];
    TH1D* h1MCphotonEgyInAngleBin[NAnglebins];



    TH1D* h1MCPrimariesID;
    TH1D* h1MCPrimariesPerLayer;
    TH1D* h1electronEgy;
    TH1D* h1positronEgy;
    TH1D* h1photonEgy;
    TH1D* h1neutronEgy;
    TH1D* h1electronEgyLow;
    TH1D* h1positronEgyLow;
    TH1D* h1photonEgyLow;
    TH1D* h1neutronEgyLow;
    TH1D* h1electronInLayerEgyLow;
    TH1D* h1positronInLayerEgyLow;
    TH1D* h1photonInLayerEgyLow;
    TH1D* h1neutronInLayerEgyLow;

    TH1D* h1MCParticleTrackLength;
    TH1D* h1MCParticleTrackLengthZ;
    TH1D* h1ElectronTrackLength;
    TH1D* h1ElectronTrackLengthZ;


    // CaloOutput tree like in testbeam data called Frames with vectors
    TTree* Frames;
    std::vector<Int_t> lane;
    std::vector<Int_t> column;
    std::vector<Int_t> row;
    int eventNumber;
    int nHits;
    std::vector<int> nMCParticles;
    std::vector<int> particleID;
    std::vector<double> particleE;
    std::vector<std::string> process;
    std::vector<std::vector<int>> particleIDs;
    std::vector<std::vector<ROOT::Math::XYZPoint>> LocalStartPoint;
    std::vector<std::vector<ROOT::Math::XYZPoint>> LocalEndPoint;
    std::vector<const MCParticle*> MCPartList;
    std::vector<const MCParticle*> MCPartElectronPositronList;
    std::vector<const MCParticle*> MCPartElectronList;
    std::vector<const MCParticle*> MCPartPositronList;
    std::vector<const MCParticle*> MCPartPhotonList;
    std::vector<const MCParticle*> MCPrimariesList;
  };
} // namespace allpix
