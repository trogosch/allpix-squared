/**
* @file
* @brief Implementation of [CaloOutputWriter] module
* @copyright Copyright (c) 2017-2020 CERN and the Allpix Squared authors.
* This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
* In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
* Intergovernmental Organization or submit itself to any jurisdiction.
*/

#include "CaloOutputWriterModule.hpp"


#include <algorithm>
#include <memory>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include "TFile.h"


#include <TBranchElement.h>
#include <TClass.h>

#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include "TFile.h"
#include "TTree.h"
#include "TFrame.h"
#include "TBenchmark.h"
#include "TSystem.h"

#include "core/config/ConfigReader.hpp"
#include "core/utils/file.h"
#include "core/utils/log.h"
#include "core/utils/type.h"
#include "core/messenger/Messenger.hpp"
#include "core/geometry/HybridPixelDetectorModel.hpp"

#include "objects/PixelHit.hpp"
#include "Cluster.hpp"
#include "objects/Object.hpp"
#include "objects/objects.h"

#include "tools/ROOT.h"

// #ifdef __MAKECINT__
//  #pragma link C++ class vector<int> +;
//  #pragma link C++ class vector<vector<int>> +;
//  #pragma link C++ class vector<vector<ROOT::Math::XYZPoint>> +;
//  #endif

using namespace allpix;

CaloOutputWriterModule::CaloOutputWriterModule(Configuration& config, Messenger* messenger, GeometryManager* geo_manager)
: Module(config), geo_manager_(geo_manager), messenger_(messenger) {

  // ... Implement ... (Typically bounds the required messages and optionally sets configuration defaults)
  // Input required by this module
  messenger_->bindMulti(this, &CaloOutputWriterModule::pixel_messages_);
  messenger_->bindMulti(this, &CaloOutputWriterModule::mcparticle_messages_, MsgFlags::REQUIRED);
  // messenger_->bindMulti(this, &CaloOutputWriterModule::mctrack_messages_, MsgFlags::REQUIRED);

}

// ROOTObjectWriterModule::~CaloOutputWriterModule() {
//     // Delete all object pointers
//     for(auto& index_data : write_list_) {
//         delete index_data.second;
//     }
// }


void CaloOutputWriterModule::init() {
  // gROOT->ProcessLine("#include <vector>");
  using namespace ROOT::Math;

  // gStyle->SetPalette(kBird);
  gStyle->SetPalette(kRainBow);
  TColor::InvertPalette();

  // Loop over detectors and do something
  std::vector<std::shared_ptr<Detector>> detectors = geo_manager_->getDetectors();
  for(auto& detector : detectors) {
    // Get the detector name
    std::string detectorName = detector->getName();
    LOG(DEBUG) << "Detector with name " << detectorName;
    // std::cout << "Detector with name " << detectorName << std::endl;
  }

  //Root TTree for output similar to testbeam data
  Frames = new TTree("Frames", "Frames");
  Frames ->Branch("eventNumber", &eventNumber);
  Frames ->Branch("nHits", &nHits);
  Frames ->Branch("lane", &lane);
  Frames ->Branch("column", &column);
  Frames ->Branch("row", &row);
  // Frames ->Branch("nMCParticles", &nMCParticles);
  // Frames ->Branch("particleID", &particleID);
  // Frames ->Branch("particleE", &particleE);
  // Frames ->Branch("process", &process);

  // Frames ->Branch("particleIDs", &particleIDs);
  // Frames ->Branch("LocalStartPoint", &LocalStartPoint);
  // Frames ->Branch("LocalEndPoint", &LocalEndPoint);

  //Histograms Pixel Hits
  h1PixelHits = new TH1D("h1PixelHits", "h1PixelHits;pixel hits; N",3000, 0, 3000);
  h1PixelHits->Sumw2();
  h1PixelHits2 = new TH1D("h1PixelHits2", "h1PixelHits2;pixel hits; N",5000, 25000, 30000);
  h1PixelHits2->Sumw2();
  h1Clusters = new TH1D("h1Clusters", "h1Clusters;clusters; N",1000, 0, 1000);
  h1Clusters->Sumw2();
  h1Clusters2 = new TH1D("h1Clusters2", "h1Clusters2;clusters; N",30000, 0, 30000);
  h1Clusters2->Sumw2();
  h1ParticlesContributing2PixelHit = new TH1D("h1ParticlesContributing2PixelHit", "h1ParticlesContributing2PixelHit;MCParticles per pixel hit; N",15, 0.5, 15.5);
  h1ParticlesContributing2PixelHit->Sumw2();
  // for(int ith1=0;ith1<Nhistlayer;ith1++){
  //   h1ParticlesContributing2PixelHitLayer[ith1] = new TH1D(Form("h1ParticlesContributing2PixelHitLayer[%i]",ith1), Form("h1ParticlesContributing2PixelHitLayer[%i];MCParticles per pixel hit; N",ith1),15, 0.5, 15.5);
  //   h1ParticlesContributing2PixelHitLayer[ith1]->Sumw2();
  // }
  h1PixelHitsPerLayer = new TH1D("h1PixelHitsPerLayer", "h1PixelHitsPerLayer; layer; pixel hits per layer", Nhistlayer, -0.5, 23.5);
  h1PixelHitsPerLayer->Sumw2();

  //Histograms MC Particles
  h1MCParticles = new TH1D("h1MCParticles", "h1MCParticles; MC particles; N",200, 0, 6000);
  h1MCElectronsPositrons = new TH1D("h1MCElectronsPositrons", "h1MCElectronsPositrons; MC electron/positrons; N",200, 0, 3000);
  h1MCPhotons = new TH1D("h1MCPhotons", "h1MCPhotons; MC photons; N",200, 0, 3000);
  h1MCElectrons = new TH1D("h1MCElectrons", "h1MCElectrons; MC electrons; N",200, 0, 3000);
  h1MCPositrons = new TH1D("h1MCPositrons", "h1MCPositrons; MC positrons; N",200, 0, 3000);
  h1MCParticlesAngle2zAxis = new TH1D("h1MCParticlesAngle2zAxis", "h1MCParticlesAngle2zAxis; Theta ; N",1000, 0, TMath::Pi());
  h1MCPrimariessAngle2zAxisLayer0 = new TH1D("h1MCPrimariessAngle2zAxisLayer0", "h1MCPrimariessAngle2zAxisLayer0; Theta ; N",1000, 0, TMath::Pi());
  h1MCParticleInLayerID = new TH1D("h1MCParticleInLayerID", "h1MCParticleInLayerID; ID ; N",200, -30, 30);
  h1MCParticlesPerLayer = new TH1D("h1MCParticlesPerLayer", "h1MCParticlesPerLayer; layer; MC Particles per layer", Nhistlayer, -0.5, 23.5);
  h1MCParticles->Sumw2();
  h1MCElectronsPositrons->Sumw2();
  h1MCPhotons->Sumw2();
  h1MCElectrons->Sumw2();
  h1MCPositrons->Sumw2();
  h1MCParticlesAngle2zAxis->Sumw2();
  h1MCPrimariessAngle2zAxisLayer0->Sumw2();
  h1MCParticleInLayerID->Sumw2();
  h1MCParticlesPerLayer->Sumw2();

  h1electronEgy = new TH1D("h1electronEgy", "h1electronEgyDistrb; E (MeV); dN/dE (MeV^(-1))", 500, 0, 1000);
  h1positronEgy = new TH1D("h1positronEgy", "h1positronEgyDistrb; E (MeV); dN/dE (MeV^(-1))", 500, 0, 1000);
  h1photonEgy = new TH1D("h1photonEgy", "h1photonEgyDistrb; E (MeV); dN/dE (MeV^(-1))", 1000, 0, 1000);
  h1neutronEgy = new TH1D("h1neutronEgy", "h1neutronEgyDistrb; E (MeV); dN/dE (MeV^(-1))", 500, 0, 1000);
  h1electronEgy->Sumw2();
  h1positronEgy->Sumw2();
  h1photonEgy->Sumw2();
  h1neutronEgy->Sumw2();

  h1electronEgyLow = new TH1D("h1electronEgyLow", "h1electronEgyLow; E (MeV); dN/dE (MeV^(-1))", 1000, 0, 20);
  h1positronEgyLow = new TH1D("h1positronEgyLow", "h1positronEgyLow; E (MeV); dN/dE (MeV^(-1))", 1000, 0, 20);
  h1photonEgyLow = new TH1D("h1photonEgyLow", "h1photonEgyLow; E (MeV); dN/dE (MeV^(-1))", 1000, 0, 20);
  h1neutronEgyLow = new TH1D("h1neutronEgyLow", "h1neutronEgyLow; E (MeV); dN/dE (MeV^(-1))", 1000, 0, 20);
  h1electronEgyLow->Sumw2();
  h1positronEgyLow->Sumw2();
  h1photonEgyLow->Sumw2();
  h1neutronEgyLow->Sumw2();

  h1electronInLayerEgyLow = new TH1D("h1electronInLayerEgyLow", "h1electronInLayerEgyLow; E (MeV); dN/dE (MeV^(-1))", 1000, 0, 20);
  h1positronInLayerEgyLow = new TH1D("h1positronInLayerEgyLow", "h1positronInLayerEgyLow; E (MeV); dN/dE (MeV^(-1))", 1000, 0, 20);
  h1photonInLayerEgyLow = new TH1D("h1photonInLayerEgyLow", "h1photonInLayerEgyLow; E (MeV); dN/dE (MeV^(-1))", 1000, 0, 20);
  h1neutronInLayerEgyLow = new TH1D("h1neutronInLayerEgyLow", "h1neutronInLayerEgyLow; E (MeV); dN/dE (MeV^(-1))", 1000, 0, 20);
  h1electronInLayerEgyLow->Sumw2();
  h1positronInLayerEgyLow->Sumw2();
  h1photonInLayerEgyLow->Sumw2();
  h1neutronInLayerEgyLow->Sumw2();


  //Histograms MC Primary Particles
  h1MCPrimaries = new TH1D("h1MCPrimaries", "h1MCPrimaries;primaries;dN/dprimaries",200, 0, 6000);
  h1MCPrimariesID = new TH1D("h1MCPrimariesID", "h1MCPrimariesID;ID;dN/dID",8000, -4000, 4000);
  h1MCPrimariesAngle2zAxis = new TH1D("h1MCPrimariesAngle2zAxis", "h1MCPrimariesAngle2zAxis; Theta ; N",500, 0, TMath::Pi());
  h1MCelectronAngle2zAxis = new TH1D("h1MCelectronAngle2zAxis", "h1MCelectronAngle2zAxis; #theta ; dN/d#theta",500, 0, TMath::Pi());
  h1MCpositronAngle2zAxis = new TH1D("h1MCpositronAngle2zAxis", "h1MCpositronAngle2zAxis; #theta ; dN/d#theta",500, 0, TMath::Pi());
  h1MCphotonAngle2zAxis = new TH1D("h1MCphotonAngle2zAxis", "h1MCphotonAngle2zAxis; #theta ; dN/d#theta",500, 0, TMath::Pi());
  // for(int angleBin = 0; angleBin<NAnglebins; angleBin++){
  //   h1MCPrimariesEgyInAngleBin[angleBin] = new TH1D(Form("h1MCPrimariesEgyInAngleBin[%i]",angleBin), Form("h1MCPrimariesEgyInAngleBin[%i]; E (MeV); dN/dE (MeV^(-1))", angleBin),1000, 0, 25);
  //   h1MCPrimariesEgyInAngleBin[angleBin]->Sumw2();
  //   h1MCelectronEgyInAngleBin[angleBin] = new TH1D(Form("h1MCelectronEgyInAngleBin[%i]",angleBin), Form("h1MCelectronEgyInAngleBin[%i]; E (MeV); dN/dE (MeV^(-1))", angleBin),1000, 0, 25);
  //   h1MCelectronEgyInAngleBin[angleBin]->Sumw2();
  //   h1MCpositronEgyInAngleBin[angleBin] = new TH1D(Form("h1MCpositronEgyInAngleBin[%i]",angleBin), Form("h1MCpositronEgyInAngleBin[%i]; E (MeV); dN/dE (MeV^(-1))", angleBin),1000, 0, 25);
  //   h1MCpositronEgyInAngleBin[angleBin]->Sumw2();
  //   h1MCphotonEgyInAngleBin[angleBin] = new TH1D(Form("h1MCphotonEgyInAngleBin[%i]",angleBin), Form("h1MCphotonEgyInAngleBin[%i]; E (MeV); dN/dE (MeV^(-1))", angleBin),1000, 0, 25);
  //   h1MCphotonEgyInAngleBin[angleBin]->Sumw2();
  // }
  h1MCPrimariesPerLayer = new TH1D("h1MCPrimariesPerLayer", "h1MCPrimariesPerLayer;layer;primaries per layer", Nhistlayer, -0.5, 23.5);
  //
  //
  // h3MCParticlesAngleEnergyLayer = new TH3D("h3MCParticlesAngleEnergyLayer", "h3MCParticlesAngleEnergyLayer; #theta ; energy (MeV/c^{2}) ; layer", 500, 0, TMath::Pi(), 500, 0, 25, Nhistlayer, -0.5, 47.5);
  // h3MCParticlesAngleEnergyLayer->Sumw2();
  // h3MCElectronAngleEnergyLayer = new TH3D("h3MCElectronAngleEnergyLayer", "h3MCElectronAngleEnergyLayer; #theta ; energy (MeV/c^{2}) ; layer", 500, 0, TMath::Pi(), 500, 0, 25, Nhistlayer, -0.5, 47.5);
  // h3MCElectronAngleEnergyLayer->Sumw2();
  // h3MCPositronAngleEnergyLayer = new TH3D("h3MCPositronAngleEnergyLayer", "h3MCPositronAngleEnergyLayer; #theta ; energy (MeV/c^{2}) ; layer", 500, 0, TMath::Pi(), 500, 0, 25, Nhistlayer, -0.5, 47.5);
  // h3MCPositronAngleEnergyLayer->Sumw2();
  // h3MCPhotonAngleEnergyLayer = new TH3D("h3MCPhotonAngleEnergyLayer", "h3MCPhotonAngleEnergyLayer; #theta ; energy (MeV/c^{2}) ; layer", 500, 0, TMath::Pi(), 500, 0, 25, Nhistlayer, -0.5, 47.5);
  // h3MCPhotonAngleEnergyLayer->Sumw2();






  h1MCPrimaries->Sumw2();
  h1MCPrimariesID->Sumw2();
  h1MCPrimariesAngle2zAxis->Sumw2();
  h1MCelectronAngle2zAxis->Sumw2();
  h1MCpositronAngle2zAxis->Sumw2();
  h1MCphotonAngle2zAxis->Sumw2();
  h1MCPrimariesPerLayer->Sumw2();
  //Histograms MC Primary Particle ID
  // h1MCPrimariesPerLayer = new TH1D("h1MCPrimariesPerLayer", "h1MCPrimariesPerLayer;layer;primaries per layer", Nhistlayer, -0.5, 23.5);
  // h1MCPrimariesPerLayer->Sumw2();

  h1MCParticleTrackLength = new TH1D("h1MCParticleTrackLength", "h1MCParticleTrackLength; length (mm); N", 400, 0, 200);
  h1MCParticleTrackLengthZ = new TH1D("h1MCParticleTrackLengthZ", "h1MCParticleTrackLengthZ; length (mm); N", 800, -200, 200);
  h1ElectronTrackLength = new TH1D("h1ElectronTrackLength", "h1ElectronTrackLength; length (mm); N", 400, 0, 200);
  h1ElectronTrackLengthZ = new TH1D("h1ElectronTrackLengthZ", "h1ElectronTrackLengthZ; length (mm); N", 800, -200, 200);
  h1MCParticleTrackLength->Sumw2();
  h1MCParticleTrackLengthZ->Sumw2();
  h1ElectronTrackLength->Sumw2();
  h1ElectronTrackLengthZ->Sumw2();

}

void CaloOutputWriterModule::run(unsigned int event_num) {
  gROOT->ProcessLine("#include <vector>");
  // ... Implement ... (Typically uses the configuration to execute function and outputs an message)
  // Loop through all received messages and print some information

  int ilayer = 0;
  int HitsPerEvent = 0;
  double ClustersPerEvent = 0.0;

  bool basic = 0;

  if(basic==1){
    for(auto& message : pixel_messages_) { //loop over the messages per layer
      std::string detectorName = message->getDetector()->getName();
      LOG(DEBUG) << "Picked up " << message->getData().size() << " objects from detector " << detectorName;
      // std::cout  << "\nPicked up " << message->getData().size() << " objects from detector " << detectorName << std::endl;

      // std::cout << "\nilayer befoer call" << getLayerNumber(detectorName) << std::endl;
      // ilayer = getLayerNumber(detectorName);
      // std::cout << "\nilayer after call" << getLayerNumber(detectorName) << std::endl;

      addMCParticlesToList(ilayer);

      // Loop through pixel hits in a detector layer
      int HitsPerLayer = 0;
      LOG(TRACE) << "Loop over pixel Hits in Layer\t" << getLayerNumber(detectorName) << " (internal layer " << ilayer << ")" << std::endl;
      // std::cout << "\nLoop over pixel Hits in Layer\t" << getLayerNumber(detectorName) << " (internal layer " << ilayer << ")" << std::endl;
      for(auto& pixel_hit : message->getData()) {
        auto pixel_idx = pixel_hit.getPixel().getIndex();

        // Add pixel
        // if(pixel_idx.x()==512 || pixel_idx.x()==513 || pixel_idx.x()==514 || pixel_idx.x()==515) continue;
        HitsPerLayer++;
        HitsPerEvent++;

        // lane.push_back( getLaneFORONECHIP(ilayer, static_cast<int>(pixel_idx.x())) );
        // if(static_cast<int>(pixel_idx.x())<512) row.push_back( 511-static_cast<int>(pixel_idx.x()) );
        // if(static_cast<int>(pixel_idx.x())>512) row.push_back( static_cast<int>(pixel_idx.x())-516 );
        lane.push_back( getLane(getLayerNumber(detectorName)) );
        row.push_back( static_cast<int>(pixel_idx.x()));
        column.push_back( static_cast<int>(pixel_idx.y()) );
        nMCParticles.push_back( static_cast<int>(pixel_hit.getMCParticles().size()) );

      } //END: for(auto& pixel_hit : message->getData())

      h1PixelHitsPerLayer->AddBinContent(getLayerNumber(detectorName)+1, HitsPerLayer);

      // Perform a clustering
      std::vector<Cluster> clusters = doClustering(ilayer);
      // std::cout << "Number of Clusters in current layer\t" << static_cast<double>(clusters.size()) << std::endl;
      ClustersPerEvent += static_cast<double>(clusters.size());

      ilayer++;
    }
  }



  if(basic==0){
    for(auto& message : pixel_messages_) { //loop over the messages per layer
      std::string detectorName = message->getDetector()->getName();
      LOG(DEBUG) << "Picked up " << message->getData().size() << " objects from detector " << detectorName;
      // std::cout  << "\nPicked up " << message->getData().size() << " objects from detector " << detectorName << std::endl;

      // std::cout << "\nilayer befoer call" << getLayerNumber(detectorName) << std::endl;
      // ilayer = getLayerNumber(detectorName);
      // std::cout << "\nilayer after call" << getLayerNumber(detectorName) << std::endl;

      addMCParticlesToList(ilayer);


      std::vector<int> curr_particleIDs;
      std::vector<ROOT::Math::XYZPoint> curr_LocalStartPoint;
      std::vector<ROOT::Math::XYZPoint> curr_LocalEndPoint;
      // Loop through pixel hits in a detector layer
      int HitsPerLayer = 0;
      LOG(TRACE) << "Loop over pixel Hits in Layer\t" << getLayerNumber(detectorName) << " (internal layer " << ilayer << ")" << std::endl;
      // std::cout << "\nLoop over pixel Hits in Layer\t" << getLayerNumber(detectorName) << " (internal layer " << ilayer << ")" << std::endl;
      for(auto& pixel_hit : message->getData()) {
        auto pixel_idx = pixel_hit.getPixel().getIndex();

        // Add pixel
        // if(pixel_idx.x()==512 || pixel_idx.x()==513 || pixel_idx.x()==514 || pixel_idx.x()==515) continue;
        HitsPerLayer++;
        HitsPerEvent++;

        //Loop through list of particles associated to a pixel hit
        //If particle in list has a parent in the list, it is not counted as particle contributing to this pixel hit
        std::vector<bool> ParticleHasParentInPixelHit;
        for(int ipart=0; ipart<static_cast<int>(pixel_hit.getMCParticles().size()); ipart++){

          // Cout welche Teilchen unter welchem Prozess entstehen
          // std::cout << "Teilchen " << pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getTrack()->getParticleID() << "\t from \t" << pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getTrack()->getCreationProcessName() << std::endl;

          if(pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent() == nullptr){
            ParticleHasParentInPixelHit.push_back(0);
            continue;
          }
          bool ipart_hasParent=0;

          for(int jpart=0; jpart<static_cast<int>(pixel_hit.getMCParticles().size()); jpart++){
            if(ipart==jpart) continue;

            if(pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getTrack()->getParticleID() == pixel_hit.getMCParticles().at(static_cast<unsigned long>(jpart))->getTrack()->getParticleID()
            && pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getTrack()->getKineticEnergyInitial() == pixel_hit.getMCParticles().at(static_cast<unsigned long>(jpart))->getTrack()->getKineticEnergyInitial()){
              ipart_hasParent = 1;
              break;
            }
            if(pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getParent() != nullptr) {
              if(pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getParent()->getTrack()->getParticleID() == pixel_hit.getMCParticles().at(static_cast<unsigned long>(jpart))->getTrack()->getParticleID()
              && pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getParent()->getTrack()->getKineticEnergyInitial() == pixel_hit.getMCParticles().at(static_cast<unsigned long>(jpart))->getTrack()->getKineticEnergyInitial()){
                ipart_hasParent = 1;
                break;
              }
            }

          }
          ParticleHasParentInPixelHit.push_back(ipart_hasParent);
        }

        int nMcPartsPerPixelHit = 0;
        for(unsigned long i=0; i<ParticleHasParentInPixelHit.size(); i++){
          if(ParticleHasParentInPixelHit.at(i)==0) nMcPartsPerPixelHit++;
        }

        // parent could be outside the pixel looked at, if particles have the same parent they are counted as one particle
        int remainingSecParts=0;
        int sameParentSecParticles=0;
        if(nMcPartsPerPixelHit>1){
          for(int ipart=0; ipart<static_cast<int>(pixel_hit.getMCParticles().size()); ipart++){
            if(ParticleHasParentInPixelHit.at(static_cast<unsigned long>(ipart)) != 0) continue;
            if(pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent() == nullptr) continue;
            remainingSecParts++;
            // sameParentSecParticles++;

            for(int jpart=0; jpart<static_cast<int>(pixel_hit.getMCParticles().size()); jpart++){
              if(ipart==jpart) continue;
              if(ParticleHasParentInPixelHit.at(static_cast<unsigned long>(jpart)) != 0) continue;
              if(pixel_hit.getMCParticles().at(static_cast<unsigned long>(jpart))->getParent() == nullptr) continue;
              if(pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getTrack()->getParticleID() == pixel_hit.getMCParticles().at(static_cast<unsigned long>(jpart))->getParent()->getTrack()->getParticleID()
              && pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getTrack()->getKineticEnergyInitial() == pixel_hit.getMCParticles().at(static_cast<unsigned long>(jpart))->getParent()->getTrack()->getKineticEnergyInitial()){
                sameParentSecParticles++;
                break;
              }
              // if(pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getParent() != nullptr) {
              //   if(pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getParent()->getTrack()->getParticleID() == pixel_hit.getMCParticles().at(static_cast<unsigned long>(jpart))->getParent()->getParent()->getTrack()->getParticleID()
              //   && pixel_hit.getMCParticles().at(static_cast<unsigned long>(ipart))->getParent()->getParent()->getTrack()->getKineticEnergyInitial() == pixel_hit.getMCParticles().at(static_cast<unsigned long>(jpart))->getParent()->getParent()->getTrack()->getKineticEnergyInitial()){
              //     remainingSecParts--;
              //     break;
              //   }
              // }

            }
          }
        }


        if(nMcPartsPerPixelHit>1 && sameParentSecParticles!=0){
          nMcPartsPerPixelHit = nMcPartsPerPixelHit-sameParentSecParticles+1;
        }


        ParticleHasParentInPixelHit.clear();
        ParticleHasParentInPixelHit.shrink_to_fit();


        // h1ParticlesContributing2PixelHitLayer[ilayer]->Fill(nMcPartsPerPixelHit);
        h1ParticlesContributing2PixelHit->Fill(nMcPartsPerPixelHit);

        // if(nMcPartsPerPixelHit>3){
        //   std::cout <<"pixel hit\t" <<  pixel_idx.x() << "\t" << pixel_idx.y() << "\twith N MCParticle hits: " << pixel_hit.getMCParticles().size()<< " in Layer " << ilayer << std::endl;
        //   std::cout <<" counted hits with no parent in pixel: " << nMcPartsPerPixelHit << std::endl;
        //   for(auto& mcparticle : pixel_hit.getMCParticles()){
        //     if(mcparticle->getParent() == nullptr){
        //       std::cout << "Primary Particle:\t";
        //       std::cout << mcparticle->getTrack()->getParticleID() << "\t";
        //       std::cout << mcparticle->getTrack()->getCreationProcessName() << "\t";
        //       std::cout << mcparticle->getTrack()->getKineticEnergyInitial() << "\n";
        //     }
        //     if(mcparticle->getParent() != nullptr){
        //       std::cout << "Secondary Particle:\t";
        //       std::cout << mcparticle->getTrack()->getParticleID() << "\t";
        //       std::cout << mcparticle->getTrack()->getCreationProcessName() << "\t";
        //       std::cout << mcparticle->getTrack()->getKineticEnergyInitial() << "\t";
        //       std::cout << "with parent: " <<  mcparticle->getParent()->getTrack()->getParticleID() << "\t";
        //       std::cout << mcparticle->getParent()->getTrack()->getCreationProcessName() << "\t";
        //       std::cout << mcparticle->getParent()->getTrack()->getKineticEnergyInitial()<< "\t";
        //       if(mcparticle->getParent()->getParent() == nullptr) std::cout << std::endl;
        //       if(mcparticle->getParent()->getParent() != nullptr) {
        //         std::cout << "with parent " << mcparticle->getParent()->getParent()->getTrack()->getParticleID() << "\t";
        //         std::cout << mcparticle->getParent()->getParent()->getTrack()->getCreationProcessName() << "\t";
        //         std::cout << mcparticle->getParent()->getParent()->getTrack()->getKineticEnergyInitial() << "\n";
        //       }
        //     }
        //   }
        //   std::cout << std::endl;
        // }







        if(pixel_hit.getMCParticles().size()==1){
          particleID.push_back(pixel_hit.getMCParticles().at(0)->getParticleID());
          particleE.push_back(pixel_hit.getMCParticles().at(0)->getTrack()->getKineticEnergyInitial());
          process.push_back(pixel_hit.getMCParticles().at(0)->getTrack()->getCreationProcessName());
        }
        if(pixel_hit.getMCParticles().size()!=1){
          particleID.push_back(0);
          particleE.push_back(0);
          process.push_back("");
        }

        // lane.push_back( getLaneFORONECHIP(ilayer, static_cast<int>(pixel_idx.x())) );
        // if(static_cast<int>(pixel_idx.x())<512) row.push_back( 511-static_cast<int>(pixel_idx.x()) );
        // if(static_cast<int>(pixel_idx.x())>512) row.push_back( static_cast<int>(pixel_idx.x())-516 );
        lane.push_back( getLane(getLayerNumber(detectorName)) );
        row.push_back( static_cast<int>(pixel_idx.x()));
        column.push_back( static_cast<int>(pixel_idx.y()) );
        nMCParticles.push_back( static_cast<int>(pixel_hit.getMCParticles().size()) );

        for(auto& ipart : pixel_hit.getMCParticles()){
          curr_particleIDs.push_back( static_cast<int>(ipart->getParticleID()) );
          curr_LocalStartPoint.push_back(ipart->getLocalEndPoint());
          curr_LocalEndPoint.push_back(ipart->getLocalStartPoint());
        }

        particleIDs.push_back(curr_particleIDs);
        LocalStartPoint.push_back(curr_LocalStartPoint);
        LocalEndPoint.push_back(curr_LocalEndPoint);

        curr_particleIDs.clear();
        curr_particleIDs.shrink_to_fit();
        curr_LocalStartPoint.clear();
        curr_LocalStartPoint.shrink_to_fit();
        curr_LocalEndPoint.clear();
        curr_LocalEndPoint.shrink_to_fit();
      } //END: for(auto& pixel_hit : message->getData())

      h1PixelHitsPerLayer->AddBinContent(getLayerNumber(detectorName)+1, HitsPerLayer);

      LOG(TRACE) << "Filled histograms with primaries";



      // Perform a clustering
      std::vector<Cluster> clusters = doClustering(ilayer);
      // std::cout << "Number of Clusters in current layer\t" << static_cast<double>(clusters.size()) << std::endl;
      ClustersPerEvent += static_cast<double>(clusters.size());





      // std::cout << ilayer << std::endl;
      // std::cout << getLayerNumber(detectorName) << std::endl;
      // if(getLayerNumber(detectorName)==0 || getLayerNumber(detectorName)==1){
      // std::cout << mctrack_messages_.size() << std::endl;
      // for(const auto& mc_particle : mctrack_messages_.at(static_cast<unsigned long>(0))->getData()) {
      //
      //   // std::cout << "MCParticle " << mc_particle.getParticleID() << " produced in volume:\t " << mc_particle.getTrack()->getOriginatingVolumeName() << std::endl;
      //   std::cout << "Primary Particle:\t";
      //   std::cout << mc_particle.getParticleID() << "\t";
      //   std::cout << mc_particle.getCreationProcessName() << "\t";
      //   std::cout << mc_particle.getAngleToZAxis() << "\t";
      //   std::cout << mc_particle.getKineticEnergyInitial() << "\n";
      //   // Check for possible parents:
      //   const auto* parent = mc_particle.getParent();
      //   if(parent != nullptr) {
      //     continue;
      //   }
      //
      //
      // }



      // Retrieve all MC particles in this detector which are primary particles (not produced within the sensor):
      auto primary_particles = getPrimaryParticles(ilayer);
      // std::cout << "\nFound Hits" <<HitsPerLayer << " and primaries " << primary_particles.size() << " in layer" << getLayerNumber(detectorName) << std::endl;
      //
      // for(const auto& particle : primary_particles){
      //   std::cout << "Primary Particle:\t";
      //   std::cout << particle->getTrack()->getParticleID() << "\t";
      //   std::cout << particle->getTrack()->getCreationProcessName() << "\t";
      //   std::cout << particle->getTrack()->getAngleToZAxis() << "\t";
      //   std::cout << particle->getTrack()->getKineticEnergyInitial() << "\n";
      // }


      // if(mcparticle_messages_.at(static_cast<unsigned long>(ilayer)) != nullptr){
      //
      //   for(const auto& mc_particle : mcparticle_messages_.at(static_cast<unsigned long>(ilayer))->getData()) {
      //     std::cout << "Particle:\t";
      //     std::cout << mc_particle.getTrack()->getParticleID() << "\t";
      //     std::cout << mc_particle.getTrack()->getCreationProcessName() << "\t";
      //     std::cout << mc_particle.getTrack()->getAngleToZAxis() << "\t";
      //     std::cout << mc_particle.getTrack()->getKineticEnergyInitial() << "\n";
      //
      //   }
      // }




      // Evaluate the clusters
      // for(const auto& clus : clusters) {
      //
      //   auto clusterPos = clus.getPosition();
      //
      //   auto cluster_particles = clus.getMCParticles();
      //
      //
      //   // if( clusterPos.x() > 186 ||  (clusterPos.x() < 186 && clusterPos.y() < 341)  || (clusterPos.x() < 186 && clusterPos.y() > 683)  ) {
      //   if( clusterPos.x() > 190 ||  (clusterPos.x() < 190 && clusterPos.y() < 335)  || (clusterPos.x() < 190 && clusterPos.y() > 690)  ) {
      //     // 5mm/26.88um = 182
      //     // 5mm/29.24um = 171
      //
      //     std::cout << "\nCluster position (x,y):\t" << clusterPos.x() << "\t" <<  clusterPos.y() << "/t with " << cluster_particles.size() << " MC belonging to Cluster:" << std::endl;
      //     for(const auto& particle : cluster_particles) {
      //
      //       // std::cout << particle->getParticleID() << "\tfrom\t" << particle->getTrack()->getCreationProcessName() << "\tat angle\t" << particle->getTrack()->getAngleToZAxis() << std::endl;
      //       // std::cout << "MC Particle (x,y,z):\t" << particlePos.x() << "\t" <<  particlePos.y() << "\t" <<  particlePos.z()<< std::endl;
      //
      //
      //
      //       if(particle->getParent() == nullptr){
      //             std::cout << "Primary Particle:\t";
      //             std::cout << particle->getTrack()->getParticleID() << "\t";
      //             std::cout << particle->getTrack()->getCreationProcessName() << "\t";
      //             std::cout << particle->getTrack()->getAngleToZAxis() << "\t";
      //             std::cout << particle->getTrack()->getKineticEnergyInitial() << "\n";
      //           }
      //           if(particle->getParent() != nullptr){
      //             std::cout << "Secondary Particle:\t";
      //             std::cout << particle->getTrack()->getParticleID() << "\t";
      //             std::cout << particle->getTrack()->getCreationProcessName() << "\t";
      //             std::cout << particle->getTrack()->getAngleToZAxis() << "\t";
      //             std::cout << particle->getTrack()->getKineticEnergyInitial() << "\t";
      //             std::cout << "with parent: " <<  particle->getParent()->getTrack()->getParticleID() << "\t";
      //             std::cout << particle->getParent()->getTrack()->getCreationProcessName() << "\t";
      //             std::cout << particle->getParent()->getTrack()->getAngleToZAxis() << "\t";
      //             std::cout << particle->getParent()->getTrack()->getKineticEnergyInitial()<< "\t";
      //             if(particle->getParent()->getParent() == nullptr) std::cout << std::endl;
      //             if(particle->getParent()->getParent() != nullptr) {
      //               std::cout << "with parent " << particle->getParent()->getParent()->getTrack()->getParticleID() << "\t";
      //               std::cout << particle->getParent()->getParent()->getTrack()->getCreationProcessName() << "\t";
      //               std::cout << particle->getParent()->getParent()->getTrack()->getAngleToZAxis() << "\t";
      //               std::cout << particle->getParent()->getParent()->getTrack()->getKineticEnergyInitial() << "\n";
      //             }
      //           }
      //         }
      //         std::cout << std::endl;
      //
      //   }
      //
      //
      //   // std::cout << "This cluster is connected to " << cluster_particles.size() << " MC particles" << std::endl;
      //
      //   // Find all particles connected to this cluster which are also primaries:
      //   std::vector<const MCParticle*> intersection;
      //   std::set_intersection(primary_particles.begin(),
      //   primary_particles.end(),
      //   cluster_particles.begin(),
      //   cluster_particles.end(),
      //   std::back_inserter(intersection));
      //
      //   // std::cout << "Matching primaries: " << intersection.size() << std::endl;
      //   for(const auto& particle : intersection) {
      //
      //     // std::cout << "MC Particle " << particle->getParticleID() << "  gehoert zu Cluster" << std::endl;
      //
      //     auto pitch = message->getDetector()->getModel()->getPixelSize();
      //
      //     // auto particlePos = particle->getLocalReferencePoint() + track_smearing(track_resolution_);
      //     auto particlePos = particle->getLocalReferencePoint();
      //     // std::cout << "MCParticle at " << Units::display(particlePos, {"mm", "um"}) << std::endl;
      //
      //     // auto inPixelPos = ROOT::Math::XYVector(std::fmod(particlePos.x() + pitch.x() / 2, pitch.x()),
      //     // std::fmod(particlePos.y() + pitch.y() / 2, pitch.y()));
      //     // std::cout << "MCParticle in pixel at " << Units::display(inPixelPos, {"mm", "um"}) << std::endl;
      //
      //     h1MCPrimariessAngle2zAxisLayer0->Fill(particle->getTrack()->getAngleToZAxis());
      //
      //     // if(particle->getTrack()->getAngleToZAxis() > 1.570796327){
      //     // if( clusterPos.x() > 186 ||  (clusterPos.x() < 186 && clusterPos.y() < 341)  || (clusterPos.x() < 186 && clusterPos.y() > 683)  ) {
      //     if( clusterPos.x() > 190 ||  (clusterPos.x() < 190 && clusterPos.y() < 335)  || (clusterPos.x() < 190 && clusterPos.y() > 690)  ) {
      //
      //       //   // 5mm/26.88um = 182
      //       //   // 5mm/29.24um = 171
      //     std::cout << "primary Particle of Cluster:" << std::endl;
      //       std::cout << particle->getParticleID() << "\tfrom\t" << particle->getTrack()->getCreationProcessName() << "\tat angle\t" << particle->getTrack()->getAngleToZAxis() << std::endl;
      //       // std::cout << "Belongs to pixel with position (x,y):\t" << clusterPos.x() << "\t" <<  clusterPos.y() << std::endl;
      //       //   std::cout << "MC Particle (x,y,z):\t" << particlePos.x() << "\t" <<  particlePos.y() << "\t" <<  particlePos.z()<< std::endl;
      //       //
      //       // }
      //
      //
      //     }
      //   }
      //
      // }
      // }





      ilayer++;
    }
  }




  h1PixelHits->Fill(HitsPerEvent);
  h1PixelHits2->Fill(HitsPerEvent);
  h1Clusters->Fill(ClustersPerEvent);
  h1Clusters2->Fill(ClustersPerEvent);
  h1MCParticles->Fill(MCPartList.size());
  h1MCPrimaries->Fill(MCPrimariesList.size());
  h1MCElectronsPositrons->Fill(MCPartElectronPositronList.size());
  h1MCElectrons->Fill(MCPartElectronList.size());
  h1MCPositrons->Fill(MCPartPositronList.size());
  h1MCPhotons->Fill(MCPartPhotonList.size());

  LOG(DEBUG) << "Found " << HitsPerEvent << " HitsPerEvent in this event" << std::endl;
  LOG(DEBUG) << "Found " << ClustersPerEvent << " HitsPerEvent in this event" << std::endl;
  // std::cout << "Found " << ClustersPerEvent << " HitsPerEvent in this event" << std::endl;
  LOG(DEBUG) << "Found " << MCPartList.size()  << " MC particles by own loop in this event" << std::endl;
  LOG(DEBUG) << "Found " << MCPrimariesList.size()  << " MC primaries by own loop in this event" << std::endl;
  LOG(DEBUG) << "Found " << MCPartElectronPositronList.size()  << " MC electrons/positrons by own loop in this event" << std::endl;
  LOG(DEBUG) << "Found " << MCPartElectronList.size()  << " MC electrons by own loop in this event" << std::endl;
  LOG(DEBUG) << "Found " << MCPartPositronList.size()  << " MC positrons by own loop in this event" << std::endl;
  LOG(DEBUG) << "Found " << MCPartPhotonList.size()  << " MC photons by own loop in this event" << std::endl;


  if(basic==0){
  //loop over all MC primary particles in this event
  for(auto& particle : MCPrimariesList) {
    // std::cout << particle->getParticleID() << std::endl;
    // h3MCParticlesAngleEnergyLayer->Fill(particle->getTrack()->getAngleToZAxis() , particle->getTrack()->getKineticEnergyInitial() , ilayer-1);
    h1MCPrimariesAngle2zAxis->Fill(particle->getTrack()->getAngleToZAxis());
    // if(particle->getTrack()->getAngleToZAxis() < (TMath::Pi()/2.0-(7.0*TMath::TwoPi()/360.0)) ) h1MCPrimariesEgyInAngleBin[0]->Fill(particle->getTrack()->getKineticEnergyInitial());
    // if(particle->getTrack()->getAngleToZAxis() < (TMath::Pi()/2.0+(7.0*TMath::TwoPi()/360.0)) && particle->getTrack()->getAngleToZAxis() > (TMath::Pi()/2.0-(7.0*TMath::TwoPi()/360.0)) )h1MCPrimariesEgyInAngleBin[1]->Fill(particle->getTrack()->getKineticEnergyInitial());
    // if(particle->getTrack()->getAngleToZAxis() > (TMath::Pi()/2.0+(7.0*TMath::TwoPi()/360.0)) )h1MCPrimariesEgyInAngleBin[2]->Fill(particle->getTrack()->getKineticEnergyInitial());
    h1MCPrimariesID->Fill(particle->getParticleID());
  }
  for(auto& particle : MCPartElectronList) {
    // std::cout << particle->getParticleID() << std::endl;
    // h3MCElectronAngleEnergyLayer->Fill(particle->getTrack()->getAngleToZAxis() , particle->getTrack()->getKineticEnergyInitial() , ilayer-1);
    h1MCelectronAngle2zAxis->Fill(particle->getTrack()->getAngleToZAxis());
    // if(particle->getTrack()->getAngleToZAxis() < (TMath::Pi()/2.0-(7.0*TMath::TwoPi()/360.0)) ) h1MCelectronEgyInAngleBin[0]->Fill(particle->getTrack()->getKineticEnergyInitial());
    // if(particle->getTrack()->getAngleToZAxis() < (TMath::Pi()/2.0+(7.0*TMath::TwoPi()/360.0)) && particle->getTrack()->getAngleToZAxis() > (TMath::Pi()/2.0-(7.0*TMath::TwoPi()/360.0)) )h1MCelectronEgyInAngleBin[1]->Fill(particle->getTrack()->getKineticEnergyInitial());
    // if(particle->getTrack()->getAngleToZAxis() > (TMath::Pi()/2.0+(7.0*TMath::TwoPi()/360.0)) )h1MCelectronEgyInAngleBin[2]->Fill(particle->getTrack()->getKineticEnergyInitial());
  }
  for(auto& particle : MCPartPositronList) {
    // std::cout << particle->getParticleID() << std::endl;
    // h3MCPositronAngleEnergyLayer->Fill(particle->getTrack()->getAngleToZAxis() , particle->getTrack()->getKineticEnergyInitial() , ilayer-1);
    h1MCpositronAngle2zAxis->Fill(particle->getTrack()->getAngleToZAxis());
    // if(particle->getTrack()->getAngleToZAxis() < (TMath::Pi()/2.0-(7.0*TMath::TwoPi()/360.0)) ) h1MCpositronEgyInAngleBin[0]->Fill(particle->getTrack()->getKineticEnergyInitial());
    // if(particle->getTrack()->getAngleToZAxis() < (TMath::Pi()/2.0+(7.0*TMath::TwoPi()/360.0)) && particle->getTrack()->getAngleToZAxis() > (TMath::Pi()/2.0-(7.0*TMath::TwoPi()/360.0)) )h1MCpositronEgyInAngleBin[1]->Fill(particle->getTrack()->getKineticEnergyInitial());
    // if(particle->getTrack()->getAngleToZAxis() > (TMath::Pi()/2.0+(7.0*TMath::TwoPi()/360.0)) )h1MCpositronEgyInAngleBin[2]->Fill(particle->getTrack()->getKineticEnergyInitial());
  }
  for(auto& particle : MCPartPhotonList) {
    // std::cout << particle->getParticleID() << std::endl;
    // h3MCPhotonAngleEnergyLayer->Fill(particle->getTrack()->getAngleToZAxis() , particle->getTrack()->getKineticEnergyInitial() , ilayer-1);
    h1MCphotonAngle2zAxis->Fill(particle->getTrack()->getAngleToZAxis());
    // if(particle->getTrack()->getAngleToZAxis() < (TMath::Pi()/2.0-(7.0*TMath::TwoPi()/360.0)) ) h1MCphotonEgyInAngleBin[0]->Fill(particle->getTrack()->getKineticEnergyInitial());
    // if(particle->getTrack()->getAngleToZAxis() < (TMath::Pi()/2.0+(7.0*TMath::TwoPi()/360.0)) && particle->getTrack()->getAngleToZAxis() > (TMath::Pi()/2.0-(7.0*TMath::TwoPi()/360.0)) )h1MCphotonEgyInAngleBin[1]->Fill(particle->getTrack()->getKineticEnergyInitial());
    // if(particle->getTrack()->getAngleToZAxis() > (TMath::Pi()/2.0+(7.0*TMath::TwoPi()/360.0)) )h1MCphotonEgyInAngleBin[2]->Fill(particle->getTrack()->getKineticEnergyInitial());
  }

  //loop over all MC particles in this event
  for(auto& particle : MCPartList) {
    h1MCParticlesAngle2zAxis->Fill(particle->getTrack()->getAngleToZAxis());

    h1MCParticleTrackLength->Fill(particle->getTrack()->getTrackLength());
    h1MCParticleTrackLengthZ->Fill(particle->getTrack()->getTrackLengthZ());

    if(particle->getParticleID()==11){//electrons
      h1electronEgy->Fill(particle->getTrack()->getKineticEnergyInitial());
      h1electronEgyLow->Fill(particle->getTrack()->getKineticEnergyInitial());
      h1ElectronTrackLength->Fill(particle->getTrack()->getTrackLength());
      h1ElectronTrackLengthZ->Fill(particle->getTrack()->getTrackLengthZ());
    }
    if(particle->getParticleID()==-11){//positrons
      h1positronEgyLow->Fill(particle->getTrack()->getKineticEnergyInitial());
      h1positronEgy->Fill(particle->getTrack()->getKineticEnergyInitial());
      h1ElectronTrackLength->Fill(particle->getTrack()->getTrackLength());
      h1ElectronTrackLengthZ->Fill(particle->getTrack()->getTrackLengthZ());
    }
    if(particle->getParticleID()==22) {//photons
      h1photonEgy->Fill(particle->getTrack()->getKineticEnergyInitial());
      h1photonEgyLow->Fill(particle->getTrack()->getKineticEnergyInitial());
    }
    if(particle->getParticleID()==2112){//neutrons
      h1neutronEgy->Fill(particle->getTrack()->getKineticEnergyInitial());
      h1neutronEgyLow->Fill(particle->getTrack()->getKineticEnergyInitial());
    }
    if(particle->getTrack()->getAngleToZAxis() < TMath::ACos(0)+0.007 && particle->getTrack()->getAngleToZAxis() > TMath::ACos(0)-0.007){
      h1MCParticleInLayerID->Fill(particle->getTrack()->getParticleID());
      if(particle->getParticleID()==11)   h1electronInLayerEgyLow->Fill(particle->getTrack()->getKineticEnergyInitial());
      if(particle->getParticleID()==-11)  h1positronInLayerEgyLow->Fill(particle->getTrack()->getKineticEnergyInitial());
      if(particle->getParticleID()==22)   h1photonInLayerEgyLow->Fill(particle->getTrack()->getKineticEnergyInitial());
      if(particle->getParticleID()==2112) h1neutronInLayerEgyLow->Fill(particle->getTrack()->getKineticEnergyInitial());
    }
  }
}


  nHits=HitsPerEvent;
  eventNumber=static_cast<int>(event_num);
  Frames -> Fill();

  clearVectors();

  LOG(TRACE) << "End of current Event";
} //END: void CaloOutputWriterModule::run(unsigned int)



void CaloOutputWriterModule::clearVectors(){
  lane.clear();
  lane.shrink_to_fit();
  column.clear();
  column.shrink_to_fit();
  row.clear();
  row.shrink_to_fit();
  nMCParticles.clear();
  nMCParticles.shrink_to_fit();
  particleID.clear();
  particleID.shrink_to_fit();
  particleE.clear();
  particleE.shrink_to_fit();
  process.clear();
  process.shrink_to_fit();
  MCPartList.clear();
  MCPartList.shrink_to_fit();
  MCPrimariesList.clear();
  MCPrimariesList.shrink_to_fit();
  MCPartElectronPositronList.clear();
  MCPartElectronPositronList.shrink_to_fit();
  MCPartElectronList.clear();
  MCPartElectronList.shrink_to_fit();
  MCPartPositronList.clear();
  MCPartPositronList.shrink_to_fit();
  MCPartPhotonList.clear();
  MCPartPhotonList.shrink_to_fit();

  particleIDs.clear();
  particleIDs.shrink_to_fit();
  LocalStartPoint.clear();
  LocalStartPoint.shrink_to_fit();
  LocalEndPoint.clear();
  LocalEndPoint.shrink_to_fit();
}


int CaloOutputWriterModule::getLaneFORONECHIP(int layer, int idx){

  if(layer==0 && idx<512) return 35;
  if(layer==0 && idx>512) return 32;

  if(layer==1 && idx<512) return 56;
  if(layer==1 && idx>512) return 59;

  if(layer==2 && idx<512) return 33;
  if(layer==2 && idx>512) return 34;

  if(layer==3 && idx<512) return 58;
  if(layer==3 && idx>512) return 57;

  if(layer==4 && idx<512) return 37;
  if(layer==4 && idx>512) return 36;

  if(layer==5 && idx<512) return 60;
  if(layer==5 && idx>512) return 61;

  if(layer==6 && idx<512) return 55;
  if(layer==6 && idx>512) return 38;

  if(layer==7 && idx<512) return 62;
  if(layer==7 && idx>512) return 79;

  if(layer==8 && idx<512) return 53;
  if(layer==8 && idx>512) return 54;

  if(layer==9 && idx<512) return 78;
  if(layer==9 && idx>512) return 77;

  if(layer==10 && idx<512) return 51;
  if(layer==10 && idx>512) return 52;

  if(layer==11 && idx<512) return 76;
  if(layer==11 && idx>512) return 75;

  if(layer==12 && idx<512) return 49;
  if(layer==12 && idx>512) return 50;

  if(layer==13 && idx<512) return 74;
  if(layer==13 && idx>512) return 73;

  if(layer==14 && idx<512) return 47;
  if(layer==14 && idx>512) return 48;

  if(layer==15 && idx<512) return 72;
  if(layer==15 && idx>512) return 71;

  if(layer==16 && idx<512) return 45;
  if(layer==16 && idx>512) return 46;

  if(layer==17 && idx<512) return 70;
  if(layer==17 && idx>512) return 69;

  if(layer==18 && idx<512) return 43;
  if(layer==18 && idx>512) return 44;

  if(layer==19 && idx<512) return 68;
  if(layer==19 && idx>512) return 67;

  if(layer==20 && idx<512) return 41;
  if(layer==20 && idx>512) return 42;

  if(layer==21 && idx<512) return 66;
  if(layer==21 && idx>512) return 65;

  if(layer==22 && idx<512) return 39;
  if(layer==22 && idx>512) return 40;

  if(layer==23 && idx<512) return 64;
  if(layer==23 && idx>512) return 63;

  return 999;
}



int CaloOutputWriterModule::getLane(int layer){

  if(layer==0  ) return 32;
  if(layer==1  ) return 35;
  if(layer==2  ) return 59;
  if(layer==3  ) return 56;
  if(layer==4  ) return 34;
  if(layer==5  ) return 33;
  if(layer==6  ) return 57;
  if(layer==7  ) return 58;
  if(layer==8  ) return 36;
  if(layer==9  ) return 37;
  if(layer==10 ) return 61;
  if(layer==11 ) return 60;
  if(layer==12 ) return 38;
  if(layer==13 ) return 55;
  if(layer==14 ) return 79;
  if(layer==15 ) return 62;
  if(layer==16 ) return 54;
  if(layer==17 ) return 53;
  if(layer==18 ) return 77;
  if(layer==19 ) return 78;
  if(layer==20 ) return 52;
  if(layer==21 ) return 51;
  if(layer==22 ) return 75;
  if(layer==23 ) return 76;
  if(layer==24 ) return 50;
  if(layer==25 ) return 49;
  if(layer==26 ) return 73;
  if(layer==27 ) return 74;
  if(layer==28 ) return 48;
  if(layer==29 ) return 47;
  if(layer==30 ) return 71;
  if(layer==31 ) return 72;
  if(layer==32 ) return 46;
  if(layer==33 ) return 45;
  if(layer==34 ) return 69;
  if(layer==35 ) return 70;
  if(layer==36 ) return 44;
  if(layer==37 ) return 43;
  if(layer==38 ) return 67;
  if(layer==39 ) return 68;
  if(layer==40 ) return 42;
  if(layer==41 ) return 41;
  if(layer==42 ) return 65;
  if(layer==43 ) return 66;
  if(layer==44 ) return 40;
  if(layer==45 ) return 39;
  if(layer==46 ) return 63;
  if(layer==47 ) return 64;

  return 999;
}




int CaloOutputWriterModule::getLayerNumber(std::string name){
  if(!name.compare("layer0a")) return 0;
  if(!name.compare("layer0b")) return 1;
  if(!name.compare("layer1a")) return 2;
  if(!name.compare("layer1b")) return 3;
  if(!name.compare("layer2a")) return 4;
  if(!name.compare("layer2b")) return 5;
  if(!name.compare("layer3a")) return 6;
  if(!name.compare("layer3b")) return 7;
  if(!name.compare("layer4a")) return 8;
  if(!name.compare("layer4b")) return 9;
  if(!name.compare("layer5a")) return 10;
  if(!name.compare("layer5b")) return 11;
  if(!name.compare("layer6a")) return 12;
  if(!name.compare("layer6b")) return 13;
  if(!name.compare("layer7a")) return 14;
  if(!name.compare("layer7b")) return 15;
  if(!name.compare("layer8a")) return 16;
  if(!name.compare("layer8b")) return 17;
  if(!name.compare("layer9a")) return 18;
  if(!name.compare("layer9b")) return 19;
  if(!name.compare("layer10a")) return 20;
  if(!name.compare("layer10b")) return 21;
  if(!name.compare("layer11a")) return 22;
  if(!name.compare("layer11b")) return 23;
  if(!name.compare("layer12a")) return 24;
  if(!name.compare("layer12b")) return 25;
  if(!name.compare("layer13a")) return 26;
  if(!name.compare("layer13b")) return 27;
  if(!name.compare("layer14a")) return 28;
  if(!name.compare("layer14b")) return 29;
  if(!name.compare("layer15a")) return 30;
  if(!name.compare("layer15b")) return 31;
  if(!name.compare("layer16a")) return 32;
  if(!name.compare("layer16b")) return 33;
  if(!name.compare("layer17a")) return 34;
  if(!name.compare("layer17b")) return 35;
  if(!name.compare("layer18a")) return 36;
  if(!name.compare("layer18b")) return 37;
  if(!name.compare("layer19a")) return 38;
  if(!name.compare("layer19b")) return 39;
  if(!name.compare("layer20a")) return 40;
  if(!name.compare("layer20b")) return 41;
  if(!name.compare("layer21a")) return 42;
  if(!name.compare("layer21b")) return 43;
  if(!name.compare("layer22a")) return 44;
  if(!name.compare("layer22b")) return 45;
  if(!name.compare("layer23a")) return 46;
  if(!name.compare("layer23b")) return 47;

  return 99999;
}














































void CaloOutputWriterModule::finalize() {

  gStyle->SetPalette(kRainBow);
  TColor::InvertPalette();
  // TFile *calo_file = new TFile("/Users/trogo/Documents/simulation_studies/allpix_simulation_output/CaloTrees.root", "RECREATE");
  // calo_file->cd();
  // Frames->Print();
  // Frames->Show(0);
  Frames->Write();
  h1PixelHits->Write();
  h1PixelHits2->Write();
  h1Clusters->Write();
  h1Clusters2->Write();
  h1ParticlesContributing2PixelHit->Write();
  // for(int ith1=0;ith1<Nhistlayer;ith1++){
  //   h1ParticlesContributing2PixelHitLayer[ith1]->Write();
  // }
  h1PixelHitsPerLayer->Write();

  h1MCParticles->Write();
  h1MCElectronsPositrons->Write();
  h1MCPhotons->Write();
  h1MCElectrons->Write();
  h1MCPositrons->Write();
  h1MCParticlesAngle2zAxis->Write();
  h1MCPrimariessAngle2zAxisLayer0->Write();
  h1MCParticleInLayerID->Write();
  h1MCParticlesPerLayer->Write();
  h1MCPrimaries->Write();
  h1MCPrimariesAngle2zAxis->Write();
  h1MCelectronAngle2zAxis->Write();
  h1MCpositronAngle2zAxis->Write();
  h1MCphotonAngle2zAxis->Write();
  // for(int angleBin = 0; angleBin<NAnglebins; angleBin++){
  //   h1MCPrimariesEgyInAngleBin[angleBin]->Write();
  //   h1MCelectronEgyInAngleBin[angleBin]->Write();
  //   h1MCpositronEgyInAngleBin[angleBin]->Write();
  //   h1MCphotonEgyInAngleBin[angleBin]->Write();
  // }
  // h3MCParticlesAngleEnergyLayer->Write();
  // h3MCElectronAngleEnergyLayer->Write();
  // h3MCPositronAngleEnergyLayer->Write();
  // h3MCPhotonAngleEnergyLayer->Write();


  h1MCPrimariesID->Write();
  h1MCPrimariesPerLayer->Write();
  h1electronEgy->Write();
  h1positronEgy->Write();
  h1photonEgy->Write();
  h1neutronEgy->Write();
  h1electronEgyLow->Write();
  h1positronEgyLow->Write();
  h1photonEgyLow->Write();
  h1neutronEgyLow->Write();
  h1electronInLayerEgyLow->Write();
  h1positronInLayerEgyLow->Write();
  h1photonInLayerEgyLow->Write();
  h1neutronInLayerEgyLow->Write();
  // calo_file->Close();

  h1MCParticleTrackLength->Write();
  h1MCParticleTrackLengthZ->Write();
  h1ElectronTrackLength->Write();
  h1ElectronTrackLengthZ->Write();

} //END: void CaloOutputWriterModule::finalize()






std::vector<const MCParticle*> CaloOutputWriterModule::getPrimaryParticles(int lay) const {
  std::vector<const MCParticle*> primaries;

  LOG(TRACE) << "in getPrimaryParticles(" << lay << ")";
  LOG(TRACE) << "mcparticle_messages_.at(static_cast<unsigned long>(ilayer)).size() = " << mcparticle_messages_.at(static_cast<unsigned long>(lay))->getData().size();

  // Loop over all MCParticles available
  if(mcparticle_messages_.at(static_cast<unsigned long>(lay)) != nullptr){
    LOG(TRACE) << "mcparticle_messages_.at(static_cast<unsigned long>(lay)) is no nullptr";
    LOG(TRACE) << "mcparticle_messages_.at(static_cast<unsigned long>(lay)).size() = " << mcparticle_messages_.at(static_cast<unsigned long>(lay))->getData().size();

    for(const auto& mc_particle : mcparticle_messages_.at(static_cast<unsigned long>(lay))->getData()) {

      // std::cout << "MCParticle " << mc_particle.getParticleID() << " produced in volume:\t " << mc_particle.getTrack()->getOriginatingVolumeName() << std::endl;

      // Check for possible parents:
      const auto* parent = mc_particle.getParent();
      if(parent != nullptr) {
        LOG(TRACE) << "MCParticle " << mc_particle.getParticleID();
        continue;
      }

      // This particle has no parent particles in the regarded sensor, return it.
      LOG(TRACE) << "MCParticle " << mc_particle.getParticleID() << " (primary)";
      primaries.push_back(&mc_particle);
    }
    LOG(TRACE) << "out mc_particle primary loop";

  }
  LOG(TRACE) << "before return primaries";

  return primaries;
}





void CaloOutputWriterModule::addMCParticlesToList(int lay){

  // Loop over all MCParticles available
  if(mcparticle_messages_.at(static_cast<unsigned long>(lay)) != nullptr){
    for(auto& mc_particle : mcparticle_messages_.at(static_cast<unsigned long>(lay))->getData()) {


      // if(MCPartList.empty() == 1) MCPartList.push_back(&mc_particle);

      bool addpart = 1;
      for(auto& particle : MCPartList) {
        // std::cout <<  particle->getTrack()->getKineticEnergyInitial() << "\t" << mc_particle.getTrack()->getKineticEnergyInitial() << std::endl;
        // std::cout <<  particle->getTrack()->getParticleID() << "\t" << mc_particle.getTrack()->getParticleID() << std::endl;

        if( particle->getTrack()->getKineticEnergyInitial() == mc_particle.getTrack()->getKineticEnergyInitial()
        && particle->getTrack()->getParticleID() == mc_particle.getTrack()->getParticleID()
        && particle->getTrack()->getStartPoint() == mc_particle.getTrack()->getStartPoint()
        && particle->getTrack()->getEndPoint() == mc_particle.getTrack()->getEndPoint()
      ){
        addpart = 0;
        break;
      }
    }

    // std::cout <<  " endl loop over particles in List, " << "particle is in list?\t" << addpart << std::endl;
    if(addpart == 1) {
      MCPartList.push_back(&mc_particle);
      if(mc_particle.getTrack()->getParticleID() == 11 || mc_particle.getTrack()->getParticleID() == -11) MCPartElectronPositronList.push_back(&mc_particle);
      if(mc_particle.getTrack()->getParticleID() == 11) MCPartElectronList.push_back(&mc_particle);
      if(mc_particle.getTrack()->getParticleID() == -11) MCPartPositronList.push_back(&mc_particle);
      if(mc_particle.getTrack()->getParticleID() == 22 ) MCPartPhotonList.push_back(&mc_particle);
      if(mc_particle.getParent() == nullptr) MCPrimariesList.push_back(&mc_particle);
    }
  }
}
}



std::vector<Cluster> CaloOutputWriterModule::doClustering(int lay) {
  std::vector<Cluster> clusters;
  std::map<const PixelHit*, bool> usedPixel;

  if(pixel_messages_.at(lay) == nullptr) {
    return clusters;
  }

  auto pixel_it = pixel_messages_.at(lay)->getData().begin();
  for(; pixel_it != pixel_messages_.at(lay)->getData().end(); pixel_it++) {
    const PixelHit* pixel_hit = &(*pixel_it);

    // Check if the pixel has been used:
    if(usedPixel[pixel_hit]) {
      continue;
    }

    // Create new cluster
    Cluster cluster(pixel_hit);
    usedPixel[pixel_hit] = true;
    LOG(TRACE) << "Creating new cluster with seed: " << pixel_hit->getPixel().getIndex();

    auto touching = [&](const PixelHit* pixel) {
      auto pxi1 = pixel->getIndex();
      for(const auto& cluster_pixel : cluster.getPixelHits()) {

        auto distance = [](unsigned int lhs, unsigned int rhs) { return (lhs > rhs ? lhs - rhs : rhs - lhs); };

        auto pxi2 = cluster_pixel->getIndex();
        if(distance(pxi1.x(), pxi2.x()) <= 1 && distance(pxi1.y(), pxi2.y()) <= 1) {
          return true;
        }
      }
      return false;
    };

    // Keep adding pixels to the cluster:
    for(auto other_pixel = pixel_it + 1; other_pixel != pixel_messages_.at(lay)->getData().end(); other_pixel++) {
      const PixelHit* neighbor = &(*other_pixel);

      // Check if neighbor has been used or if it touches the current cluster:
      if(usedPixel[neighbor] || !touching(neighbor)) {
        continue;
      }

      cluster.addPixelHit(neighbor);
      LOG(TRACE) << "Adding pixel: " << neighbor->getPixel().getIndex();
      usedPixel[neighbor] = true;
      other_pixel = pixel_it;
    }
    clusters.push_back(cluster);
  }
  return clusters;
}
