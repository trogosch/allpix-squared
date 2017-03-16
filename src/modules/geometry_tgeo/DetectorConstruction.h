/// \file DetectorConstruction.h
/// \brief Definition of the DetectorConstruction class.
///
/// Builds the detector geometry according to user defined parameters.
///
/// \date     Feb. 13 2017
/// \version  0.9
/// \author N. Gauvin; CERN

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <vector>
#include <map>

// ROOT
#include "TGeoManager.h"


// AllPix2
#include "GeoDsc.h"

/*** Names of detector parts
 * These are extremely important and should be placed in a visible way,
 * as they will be used to retrieve the objects from the gGeoManager. 
 ***/
const TString WrapperName = "Wrapper";
const TString PCBName = "PCB";
const TString WaferName = "Wafer"; //Box in AllPix1
const TString CoverName = "Coverlayer";
const TString SliceName = "Slice";
const TString PixelName = "Pixel";
const TString ChipName = "Chip";
const TString BumpName = "Bump";
const TString GuardRingsName = "GuardRings";


class DetectorConstruction
{

 public:
  DetectorConstruction();
  ~DetectorConstruction();

  void run();
  void Construct();
  void BuildPixelDevices();
  void BuildMaterialsAndMedia();
  void BuildAppliances();
  void BuildTestStructure();
  void ReadDetectorDescriptions(); /// Debugging only !
  
  // Global variables
  TGeoMedium* m_fillingWorldMaterial; /// Medium to fill the World.
  vector<GeoDsc*> m_geoMap; /// the detector descriptions
  
  // User defined parameters
  /*
    Medium to fill the World. Available media :
    - Air
    - Vacuum
  */
  TString m_userDefinedWorldMaterial; 
  TString m_userDefinedGeoOutputFile;
  TString m_buildAppliancesFlag;
  int m_Appliances_type;
  TString m_buildTestStructureFlag;
  map<int,TVector3> m_vectorWrapperEnhancement;
  map<int,TGeoTranslation> m_posVector; // position of medipix(es), key is detector Id
  map<int,TGeoRotation> m_rotVector; // map<int, G4RotationMatrix *>
  map<int,TGeoTranslation> m_posVectorAppliances; // 
};

#endif /*DetectorConstruction_h*/
