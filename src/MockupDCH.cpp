// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author W. Armstrong
 *
 * \ingroup trackers
 *
 * @{
 */
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <array>
#include "DD4hepDetectorHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

/** Barrel Tracker with space frame.
 *
 * - Optional "support" tag within the detector element.
 *
 * The shapes are created using createShape which can be one of many basic geomtries.
 * See the examples Check_shape_*.xml in
 * [dd4hep's examples/ClientTests/compact](https://github.com/AIDASoft/DD4hep/tree/master/examples/ClientTests/compact)
 * directory.
 *
 *
 * - Optional "frame" tag within the module element.
 *
 * \ingroup trackers
 *
 * \code
 * \endcode
 *
 *
 * @author Whitney Armstrong
 */

const int NUMWIRETYPE = 6, /*S=0,*/ F2=1, F1=2, F3=3, F5=4, F4=5;
string tagWire[NUMWIRETYPE]={"s","f2","f1","f3","f5","f4"};

double alpha=0*dd4hep::deg;

Volume GetVesselAssembly(dd4hep::xml::Dimension dimensions,Detector& description,dd4hep::xml::Handle_t vesselParam, bool SHOWVESSEL);
VolPlane GetSensitiveSurface(Volume vol,PlacedVolume pv,SensitiveDetector sens,double halft);
void CheckRepeatedVolume(map<string, Volume> volumes,string m_nam);
Volume GetComponentVol(Detector& description,double r, xml_comp_t param,Material gas,bool SHOWSENSOR);
double Pitch_z0(double r_z0, int nwires);
double Stereoangle_z0(double r_z0,double Lhalf);
int StereoSign(int iLayer);
double GetAngleRotX(int ilayer,double rLayer, double layerLength);
double GetWireLength(double rLayer,double layerLength,double rWire,int nwire);
double GetBuffer(double rLayer,double rWire ,double lwire,double angle, double r0);
double fwire_phi_offset(double r_z0, double rWire,double safety_phi_interspace);
Volume GetWire(int ilayer,int iwire, double r_z0, double layerLength,int nwire, double rWire, Material mat, VisAttr vis);
Assembly GetMeshWire(int ilayer,int iwire, double rLayer, double layerLength, int nwire,double rWire, Material mat, VisAttr vis);
Transform3D GetWireTransform(int ilayer,int iwire,double r_z0,double rWire,double layerLength);

//-----------------------------------------------------------------------------------//
static Ref_t create_MockupDCH(Detector& description, xml_h e, SensitiveDetector sens)
{
  int i;
  
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det    = e;
  Material                     gas      = description.material(x_det.child("gas").attr<std::string>(_Unicode(material)));
  int                          det_id   = x_det.id();
  string                       det_name = x_det.nameStr();
  DetElement                   sdet(det_name, det_id);

  map<string, Volume>                volumes;
  map<string, Placements>            sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;

  PlacedVolume pv;

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto &params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params, "boundary_material");
  }

  //-----------------------
  // Read parameters
  // from the xml file
  //----------------------- 
  dd4hep::xml::Dimension dimensions(x_det.dimensions());
  // Tube topVolumeShape(dimensions.rmin(), dimensions.rmax(), dimensions.length() * 0.5);
  // Volume assembly(det_name,topVolumeShape,air);

  auto ConstructionFlag = x_det.child("constructionFlag");
  int  MAXLAYER         = ConstructionFlag.attr<int>(_Unicode(maxLayer));
  bool SHOWVESSEL       = ConstructionFlag.attr<bool>(_Unicode(showVessel));
  bool SHOWLAYER        = ConstructionFlag.attr<bool>(_Unicode(showLayer));
  bool SHOWSENSOR       = ConstructionFlag.attr<bool>(_Unicode(showSensor));

  double layerLength=description.constantAsDouble("Layer_length");
  double layer_rmin =description.constantAsDouble("Layer_rmin");
  double layerThickness = description.constantAsDouble("Layer_thickness");
  double layerBuffer=description.constantAsDouble("Layer_buffer");

  //double buffer     =description.constantAsDouble("buffer");

  xml_comp_t x_layer = x_det.child("layer");
  cout<<"Total number of layer = "<<x_layer.repeat()<<endl;
  //-----------------------
  // THE detector assembly
  //-----------------------
  Assembly assembly(det_name);

  sens.setType("tracker");

  //----------------------- 
  // Vessel volume
  //-----------------------
  auto vesselParam       = x_det.child("vessel");
  Volume vesselAssembly=GetVesselAssembly(dimensions,description,vesselParam,SHOWVESSEL);
  assembly.placeVolume(vesselAssembly,Position(0,0,0));

  //-----------------------
  // wire volumes
  //----------------------- 
  auto wireParam       = x_det.child("wires");

  double wireRadius[NUMWIRETYPE]  = {wireParam.attr<double>(_Unicode(SWire_thickness))/2.,
				     wireParam.attr<double>(_Unicode(FCentralWire_thickness))/2.0,
				     wireParam.attr<double>(_Unicode(FSideWire_thickness))/2.0,
				     wireParam.attr<double>(_Unicode(FSideWire_thickness))/2.0,
				     wireParam.attr<double>(_Unicode(FSideWire_thickness))/2.0,
				     wireParam.attr<double>(_Unicode(FSideWire_thickness))/2.0};

  Material wireMat[NUMWIRETYPE]   = {description.material(wireParam.attr<std::string>(_Unicode(SWire_material))),
				     description.material(wireParam.attr<std::string>(_Unicode(FCentralWire_material))),
				     description.material(wireParam.attr<std::string>(_Unicode(FSideWire_material))),
				     description.material(wireParam.attr<std::string>(_Unicode(FSideWire_material))),
				     description.material(wireParam.attr<std::string>(_Unicode(FSideWire_material))),
				     description.material(wireParam.attr<std::string>(_Unicode(FSideWire_material)))};
  
  VisAttr wireVis[NUMWIRETYPE]     = {description.visAttributes(wireParam.attr<std::string>(_Unicode(SWire_vis))),
				      description.visAttributes(wireParam.attr<std::string>(_Unicode(FWire_vis))),
				      description.visAttributes(wireParam.attr<std::string>(_Unicode(FWire_vis))),
				      description.visAttributes(wireParam.attr<std::string>(_Unicode(FWire_vis))),
				      description.visAttributes(wireParam.attr<std::string>(_Unicode(FWire_vis))),
				      description.visAttributes(wireParam.attr<std::string>(_Unicode(FWire_vis)))};
 
  bool BUILDWIRE[NUMWIRETYPE]      = {wireParam.attr<bool>(_Unicode(buildSenseWires)),
				      wireParam.attr<bool>(_Unicode(buildFieldWires)),
				      wireParam.attr<bool>(_Unicode(buildFieldWires)),
				      wireParam.attr<bool>(_Unicode(buildFieldWires)),
				      wireParam.attr<bool>(_Unicode(buildFieldWires)),
				      wireParam.attr<bool>(_Unicode(buildFieldWires))};
  
  bool SHOWWIRE[NUMWIRETYPE]       = {ConstructionFlag.attr<bool>(_Unicode(showSwire)),
				      ConstructionFlag.attr<bool>(_Unicode(showFwire)),
				      ConstructionFlag.attr<bool>(_Unicode(showFwire)),
				      ConstructionFlag.attr<bool>(_Unicode(showFwire)),
				      ConstructionFlag.attr<bool>(_Unicode(showFwire)),
				      ConstructionFlag.attr<bool>(_Unicode(showFwire))};

  alpha                            = wireParam.attr<double>(_Unicode(alpha));

  int nphi0         = 120;
  int nphiIncrement = 10;
  int nSupperLayer  = 7;

  vector<Volume> wire_vol[NUMWIRETYPE];
  vector<double> vec_rLayer;

  int ilayer=0;
  double rLayer=0.;
  for (ilayer=0;ilayer<x_layer.repeat();ilayer++) {
    if (ilayer>=MAXLAYER) break;

    if (ilayer==0) rLayer= (2*layer_rmin+(ilayer+1)*layerThickness)/2.;
    else rLayer= rLayer + layerThickness;
    vec_rLayer.push_back(rLayer);

    // 1 cell = 1 swire + 4 fwires
    int nphi=nphi0+nphiIncrement*((int) ilayer/nSupperLayer);  //ncells
    int nwire=nphi;
    
    for (i =0;i<NUMWIRETYPE;i++) {
      if (!SHOWWIRE[i]) wireVis[i]=description.invisible();
      
      double r=rLayer;
      if (i==F1 || i==F4) r=r-wireRadius[i]*2; 
      if (i==F3 || i==F5) r=r+wireRadius[i]*2;
      
      Volume meshWire=GetMeshWire(ilayer,i,r,layerLength*0.99,nwire,wireRadius[i],wireMat[i],wireVis[i]);   //*0.99 to avoid overlap
      wire_vol[i].push_back(meshWire);
    }
  }

  //----------------------- 
  // build sensitive module
  //----------------------- 
  auto moduleParam =x_det.child("module");
  double moduleThickness=moduleParam.attr<double>(_Unicode(thickness));

  for (ilayer=0;ilayer<x_layer.repeat();ilayer++) {
    string m_nam = moduleParam.attr<string>(_Unicode(name)) + _toString(ilayer+1,"%d");
    CheckRepeatedVolume(volumes,m_nam);
    
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;    

    Volume c_vol=GetComponentVol(description,vec_rLayer.at(ilayer)+layerThickness/2.-moduleThickness/2., moduleParam,gas,SHOWSENSOR);
    pv = m_vol.placeVolume(c_vol, Position(0, 0, 0));
    
    if (moduleParam.attr<bool>(_Unicode(sensitive))) {
      VolPlane surf=GetSensitiveSurface(c_vol,pv,sens,moduleThickness/2.);
      volplane_surfaces[m_nam].push_back(surf);
      sensitives[m_nam].push_back(pv);
    }
  }
  
  //-----------------------
  // build the layers
  //-----------------------
  for (ilayer=0;ilayer<x_layer.repeat();ilayer++) {
    if (ilayer>=MAXLAYER) break;

    int        lay_id   = ilayer+1; //x_layer.id();
    string     m_nam    = x_layer.moduleStr() + _toString(ilayer+1,"%d");
    string     lay_nam  = det_name + _toString(ilayer+1, "_layer%d");

    Tube       lay_tub(layer_rmin+ilayer*layerThickness+layerBuffer,layer_rmin+(ilayer+1)*layerThickness,layerLength);
    Volume     lay_vol(lay_nam, lay_tub, gas); // Create the layer envelope volume.
    Position   lay_pos(0, 0, 0);

    if (SHOWLAYER) lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));
    else lay_vol.setVisAttributes(description.invisible());

    DetElement  lay_elt(sdet, lay_nam, lay_id);

    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    auto &layerParams = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams, "layer_material");
    }

    //-----------------------
    // place module
    //-----------------------
    int    module   = 1;

    pv = lay_vol.placeVolume(volumes[m_nam], Position(0,0,0));
    pv.addPhysVolID("module", module); 
    DetElement mod_elt(lay_elt, "module1", module); 
    mod_elt.setPlacement(pv); 

    DetElement comp_de(mod_elt, std::string("de_") + sensitives[m_nam][0].volume().name(), module); 
    comp_de.setPlacement(sensitives[m_nam][0]); 

    auto &comp_de_params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
    comp_de_params.set<string>("axis_definitions", "XYZ");
    volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nam][0]); 

    //-----------------------
    //place the wires
    //-----------------------
    for (i=0;i<NUMWIRETYPE;i++) { 
      if (BUILDWIRE[i]) pv = lay_vol.placeVolume(wire_vol[i].at(ilayer),Position(0,0,0)); 
    }

    // Create the PhysicalVolume for the layer.
    pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
    lay_elt.setPlacement(pv);
  }

  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly,Position(0,0,description.constantAsDouble("DCH_z0")));
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);
  return sdet;
}
//-----------------------------------------------------------------------------------//
Volume GetComponentVol(Detector& description,double r, xml_comp_t param,Material gas,bool SHOWSENSOR)
{
  Tube c_tub(r-param.thickness()/2.,r+param.thickness()/2.,param.length());
  Volume c_vol("component1", c_tub, gas);

  c_vol.setRegion(description, param.regionStr());
  c_vol.setLimitSet(description, param.limitsStr());

  if (SHOWSENSOR) c_vol.setVisAttributes(description, param.visStr());
  else c_vol.setVisAttributes(description.invisible());

  return c_vol;
}
//-----------------------------------------------------------------------------------//
void CheckRepeatedVolume(map<string, Volume> volumes,string m_nam)
{
  if (volumes.find(m_nam) != volumes.end()) {                                                                                                                                                                         
    printout(ERROR, "MockDCH:: ",                                                                                                                                                                                     
	     string((string("Module with named ") + m_nam + string(" already exists."))).c_str());                                                                                                                    
    throw runtime_error("Logics error in building modules.");                                                                                                                                                         
  } 
}
//-----------------------------------------------------------------------------------//
VolPlane GetSensitiveSurface(Volume vol,PlacedVolume pv,SensitiveDetector sens,double halft)
//VolPlane GetSensitiveSurface(int i,Volume vol,PlacedVolume pv,SensitiveDetector sens,double t_inner, double t_outer)
{
  //pv.addPhysVolID("sensor", i);
  pv.addPhysVolID("sensor", 1);
  vol.setSensitiveDetector(sens);
  
  // -------- create a measurement plane for the tracking surface attched to the sensitive volume ----- //
  Vector3D u(-1., 0., 0.);
  Vector3D v(0., -1., 0.);
  Vector3D n(0., 0., 1.);

  SurfaceType type(SurfaceType::Sensitive);
  
  //VolPlane surf(vol, type, t_inner, t_outer, u, v, n);
  VolPlane surf(vol, type, halft, halft, u, v, n);

  return surf;
}
//-----------------------------------------------------------------------------------//
double Pitch_z0(double r_z0, int nwires) 
{
  return TMath::TwoPi()*r_z0/nwires;
}
//-----------------------------------------------------------------------------------//
double Stereoangle_z0(double r_z0,double Lhalf) 
{
  //double alpha=15*dd4hep::deg;
  return atan( r_z0/Lhalf*tan(alpha/2/dd4hep::rad));

}
//-----------------------------------------------------------------------------------// 
int StereoSign(int iLayer) 
{
  if (iLayer%2==1) return 1;
  else return -1;
}
//-----------------------------------------------------------------------------------// 
double fwire_phi_offset(double r_z0, double rWire,double safety_phi_interspace)
{
  return atan(rWire/r_z0)*dd4hep::rad + safety_phi_interspace;
}
//-----------------------------------------------------------------------------------// 
double GetWireLength(double rLayer,double layerLength,double rWire,int nwire)
{
  double safety_z_interspace=1*dd4hep::nm;

  auto pitch_z0 = Pitch_z0(rLayer,nwire);
  double l=2*layerLength/cos(atan(pitch_z0/(2*layerLength)))/cos(Stereoangle_z0(rLayer,layerLength)/dd4hep::rad) ;
  
  double wlength=0.5*l
    - rWire*cos(Stereoangle_z0(rLayer,layerLength))
    - safety_z_interspace;

  return wlength; 
}
//-----------------------------------------------------------------------------------// 
double GetBuffer(double rLayer,double rWire ,double lwire,double angle, double r0)
{
  return sqrt(pow(rLayer+rWire,2)+pow(lwire*sin(angle),2))-r0;
}
//-----------------------------------------------------------------------------------//  
Volume GetWire(int ilayer,int iwire, double r_z0, double layerLength, int nwire,double rWire, Material mat, VisAttr vis)
{
  double wlength=GetWireLength(r_z0,layerLength,rWire,nwire);

  Tube wireTub(0,rWire,wlength);
  Volume wireVol(Form("l%d_one_%s",ilayer,tagWire[iwire].c_str()),wireTub,mat);  
  wireVol.setVisAttributes(vis);

  return wireVol;
}
//-----------------------------------------------------------------------------------//
Assembly GetMeshWire(int ilayer,int iwire, double rLayer, double layerLength, int nwire,double rWire, Material mat, VisAttr vis)
{
  Assembly meshWire(Form("l%d_meshWire_%s",ilayer,tagWire[iwire].c_str()));

  double wlength=GetWireLength(rLayer,layerLength,rWire,nwire);

  Tube wireTub(0,rWire,wlength);
  Volume wireVol(Form("l%d_one_%s",ilayer,tagWire[iwire].c_str()),wireTub,mat);
  wireVol.setVisAttributes(vis);

  double phi_step = (TMath::TwoPi()/nwire/2.)*dd4hep::rad;
  Transform3D wireTr=GetWireTransform(ilayer,iwire,rLayer,rWire,layerLength);
  
  for (int iphi=0;iphi<nwire;iphi++) {
    double phi_angle = 2*phi_step * iphi;
    if (iwire==F2 || iwire==F1 || iwire==F3) phi_angle = phi_angle+phi_step;

    Transform3D cellTr { RotationZ(phi_angle) };
    
    meshWire.placeVolume(wireVol, cellTr * wireTr);
  }
  
  return meshWire;
}
//-----------------------------------------------------------------------------------//
double GetAngleRotX(int ilayer,double rLayer, double layerLength)  //stereo angle
{
  return (-1.)*StereoSign(ilayer)*Stereoangle_z0(rLayer,layerLength);
}
//-----------------------------------------------------------------------------------//
Transform3D GetWireTransform(int ilayer,int iwire,double r_z0,double rWire,double layerLength)
{
  double rLayer=r_z0;
  if (iwire==F1 || iwire==F4) rLayer=r_z0-rWire*2.;
  if (iwire==F3 || iwire==F5) rLayer=r_z0+rWire*2.;
  
  RotationX stereoTr(GetAngleRotX(ilayer,rLayer,layerLength));
  Transform3D wireTr(stereoTr * Translation3D(rLayer,0.,0.));

  return wireTr;
}

//-----------------------------------------------------------------------------------// 
Volume GetVesselAssembly(dd4hep::xml::Dimension dimensions,Detector& description,dd4hep::xml::Handle_t vesselParam, bool SHOWVESSEL)
{
  double tShell=vesselParam.attr<double>(_Unicode(tShell));
  double tFill=vesselParam.attr<double>(_Unicode(tFill));
  double tEndcap=2*tShell+tFill;

  //barrel
  Tube outbarrel_tub(dimensions.rmax()-tShell,dimensions.rmax(),(dimensions.length()/2.-tEndcap));
  Tube inbarrel_tub(dimensions.rmin(),dimensions.rmin()+tShell,(dimensions.length()/2.-tEndcap));

  //endcaps
  Tube endcap_tub(dimensions.rmin(),dimensions.rmax(),tEndcap);

  Tube filling_tub(dimensions.rmin()+tShell,dimensions.rmax()-tShell,tFill);
  Volume filling_vol("filling",filling_tub,description.material(vesselParam.attr<std::string>(_Unicode(fillMat))));
  if (!SHOWVESSEL) filling_vol.setVisAttributes(description.invisible());

  UnionSolid tmp1(outbarrel_tub,inbarrel_tub,Position(0,0,0));
  UnionSolid tmp2(tmp1,endcap_tub,Position(0, 0, -dimensions.length()/2.));//+tEndcap/2.));
  UnionSolid tmp3(tmp2,endcap_tub,Position(0, 0, dimensions.length()/2.));//-tEndcap/2.));

  Volume vessel("vessel",tmp3,description.material(vesselParam.attr<std::string>(_Unicode(shellMat))));
  if (SHOWVESSEL) vessel.setVisAttributes(description.visAttributes(vesselParam.attr<std::string>(_Unicode(vis))));
  else vessel.setVisAttributes(description.invisible());
  
  vessel.placeVolume(filling_vol,Position(0, 0, dimensions.length()/2.));
  vessel.placeVolume(filling_vol,Position(0, 0, -dimensions.length()/2.));

  return vessel;
}
//-----------------------------------------------------------------------------------//
//@}
// clang-format off
DECLARE_DETELEMENT(D2EIC_mockupDCH, create_MockupDCH)
