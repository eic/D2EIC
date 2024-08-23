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

const int NUMWIRETYPE = 6, S=0, F2=1, F1=2, F3=3, F5=4, F4=5;
string tagWire[NUMWIRETYPE]={"s","f2","f1","f3","f5","f4"};

Assembly GetVesselAssembly(dd4hep::xml::Dimension dimensions,Detector& description,dd4hep::xml::Handle_t vesselParam, bool SHOWVESSEL);
VolPlane GetSensitiveSurface(int i,Volume vol,PlacedVolume pv,SensitiveDetector sens,double t_inner, double t_outer);
Volume GetComponentVol(int i,Detector& description,Material gas,xml_comp_t x_comp,bool SHOWSENSOR);
double Pitch_z0(double r_z0, int nwires);
double Stereoangle_z0(double r_z0,double Lhalf);
int StereoSign(int iLayer);
double WireLength(double r_z0,double Lhalf,int nwire);
double fwire_phi_offset(double r_z0, double rWire,double safety_phi_interspace);
Volume GetWire(int ilayer,int iwire, double r_z0, double layerLength,int nwire, double rWire, Material mat, VisAttr vis);
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

  // constructin flag
  auto ConstructionFlag = x_det.child("constructionFlag");
  int  MAXLAYER         = ConstructionFlag.attr<int>(_Unicode(maxLayer));
  bool SHOWVESSEL       = ConstructionFlag.attr<bool>(_Unicode(showVessel));
  bool SHOWLAYER        = ConstructionFlag.attr<bool>(_Unicode(showLayer));
  bool SHOWSENSOR       = ConstructionFlag.attr<bool>(_Unicode(showSensor));

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto &params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params, "boundary_material");
  }

  dd4hep::xml::Dimension dimensions(x_det.dimensions());
  //cout<<"rmin="<<dimensions.rmin()/dd4hep::cm<<" (cm), rmax="<<dimensions.rmax()/dd4hep::cm<<" (cm), full length="<<dimensions.length()/dd4hep::cm<<" (cm)"<<endl;
  // Tube topVolumeShape(dimensions.rmin(), dimensions.rmax(), dimensions.length() * 0.5);
  // Volume assembly(det_name,topVolumeShape,air);
  Assembly assembly(det_name);

  sens.setType("tracker");

  //----------------------- 
  // Vessel volume
  //-----------------------
  auto vesselParam       = x_det.child("vessel");
  Assembly vesselAssembly=GetVesselAssembly(dimensions,description,vesselParam,SHOWVESSEL);
  assembly.placeVolume(vesselAssembly,Position(0,0,0));

  //-----------------------
  // loop over the modules
  //-----------------------
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "BarrelTrackerWithFrame",
               string((string("Module with named ") + m_nam + string(" already exists."))).c_str());
      throw runtime_error("Logics error in building modules.");
    }

    //-----------------------
    // Compute module total 
    // thickness from components
    //-----------------------
    xml_coll_t ci(x_mod, _U(module_component));
    
    double total_thickness = 0;
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) { total_thickness += xml_comp_t(ci).thickness(); }

    //-----------------------
    // the module assembly volume
    //-----------------------
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    
    double thickness_so_far = 0.0;
    double thickness_sum    = -total_thickness / 2.0;

    int    ncomponents      = 0;

    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t   x_comp = mci;
      const double zoff = thickness_sum + x_comp.thickness() / 2.0;
      
      Volume c_vol=GetComponentVol(ncomponents,description,gas,x_comp,SHOWSENSOR);

      pv = m_vol.placeVolume(c_vol, Position(0, 0, zoff));
      
      if (x_comp.isSensitive()) {
	double inner_thickness = thickness_so_far + x_comp.thickness() / 2.0;
        double outer_thickness = total_thickness - thickness_so_far - x_comp.thickness() / 2.0;

	VolPlane surf=GetSensitiveSurface(ncomponents,c_vol,pv,sens,inner_thickness,outer_thickness);
	volplane_surfaces[m_nam].push_back(surf);
	
	sensitives[m_nam].push_back(pv);
      }

      thickness_sum += x_comp.thickness();
      thickness_so_far += x_comp.thickness();
    }
  }

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

  int nphi0         = 60;
  int nphiIncrement = 10;
  int nSupperLayer  = 7;
  
  double layerLength=description.constantAsDouble("Layer_length");

  vector<Volume> wire_vol[NUMWIRETYPE];
  vector<double> vec_rLayer;

  int ilayer=0;
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    double rLayer= x_layout.rc()-40*dd4hep::um/2.-wireRadius[S]*2-wireRadius[F4]*3.;
    vec_rLayer.push_back(rLayer);

    // 1 cell = 1 swire + 4 fwires
    int nphi=nphi0+nphiIncrement*((int) ilayer/nSupperLayer);  //ncells
    int nwire=nphi;

    for (i =0;i<NUMWIRETYPE;i++) {
      if (!SHOWWIRE[i]) wireVis[i]=description.invisible();
      Volume singleWire=GetWire(ilayer,i,rLayer,layerLength,nwire,wireRadius[i],wireMat[i],wireVis[i]);
      wire_vol[i].push_back(singleWire);
    }
    
    ilayer++;
  }
  
  //
  //-----------------------
  // build the layers
  //-----------------------
  ilayer=0;
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    if (ilayer>=MAXLAYER) break;

    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    int        lay_id   = x_layer.id();
    string     m_nam    = x_layer.moduleStr();
    string     lay_nam  = det_name + _toString(x_layer.id(), "_layer%d");
    Tube       lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length());
    Volume     lay_vol(lay_nam, lay_tub, gas); // Create the layer envelope volume.
    Position   lay_pos(0, 0, getAttrOrDefault(x_barrel, _U(z0), 0.));
    if (SHOWLAYER) lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));
    else lay_vol.setVisAttributes(description.invisible());

    double phi0     = x_layout.phi0();     // Starting phi of first module.
    double phi_tilt = x_layout.phi_tilt(); // Phi tilt of a module.
    double rc       = x_layout.rc();       // Radius of the module center.
    int    nphi     = x_layout.nphi();     // Number of modules in phi.
    double rphi_dr  = x_layout.dr();       // The delta radius of every other module.
    double phi_incr = (M_PI * 2) / nphi;   // Phi increment for one module.
    double phic     = phi0;                // Phi of the module center.
    double z0       = z_layout.z0();       // Z position of first module in phi.
    double nz       = z_layout.nz();       // Number of modules to place in z.
    double z_dr     = z_layout.dr();       // Radial displacement parameter, of every other module.

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
    // Z increment for module placement along Z axis.
    // Adjust for z0 at center of module rather than
    // the end of cylindrical envelope.
    double z_incr = nz > 1 ? (2.0 * z0) / (nz - 1) : 0.0;
    // Starting z for module placement along Z axis.
    double module_z = -z0;
    int    module   = 1;

    // Loop over the number of modules in phi.
    for (int ii = 0; ii < nphi; ii++) {
      //place the sensor
      double dx = z_dr * std::cos(phic + phi_tilt); // Delta x of module position.
      double dy = z_dr * std::sin(phic + phi_tilt); // Delta y of module position.
      double x  = rc * std::cos(phic);              // Basic x module position.
      double y  = rc * std::sin(phic);              // Basic y module position.

      // Loop over the number of modules in z.
      for (int j = 0; j < nz; j++) {
	Transform3D tr(RotationZYX(0, ((M_PI / 2) - phic), -M_PI / 2), Position(x, y, module_z));
	pv = lay_vol.placeVolume(volumes[m_nam], tr);
	pv.addPhysVolID("module", module);

	DetElement mod_elt(lay_elt, Form("module%d",module), module);
	mod_elt.setPlacement(pv);

	for (size_t ic = 0; ic < sensitives[m_nam].size(); ++ic) {
	  DetElement   comp_de(mod_elt, std::string("de_") + sensitives[m_nam][ic].volume().name(), module);
	  comp_de.setPlacement(sensitives[m_nam][ic]);

	  auto &comp_de_params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
	  comp_de_params.set<string>("axis_definitions", "XYZ");
	  volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nam][ic]);
	}

	/// Increase counters etc.
	module++;
	// Adjust the x and y coordinates of the module.
	x += dx;
	y += dy;
	// Flip sign of x and y adjustments.
	dx *= -1;
	dy *= -1;
	// Add z increment to get next z placement pos.
	module_z += z_incr;
      }
      phic += phi_incr; // Increment the phi placement of module.
      rc += rphi_dr;    // Increment the center radius according to dr parameter.
      rphi_dr *= -1;    // Flip sign of dr parameter.
      module_z = -z0;   // Reset the Z placement parameter for module.
    }
    
    //-----------------------
    //place the wires
    //-----------------------
    int    nWirePhi = (nphi0+nphiIncrement*((int) ilayer/nSupperLayer)); 
    double phi_step = (TMath::TwoPi()/nWirePhi/2.)*dd4hep::rad;
    
    for (i=0;i<NUMWIRETYPE;i++) {
      Transform3D wireTr=GetWireTransform(ilayer,i,vec_rLayer.at(ilayer),wireRadius[i],layerLength);

      for (int iphi=0;iphi<nWirePhi;iphi++) {
	double phi_angle = 2*phi_step * iphi;
	if (i==F2 || i==F1 || i==F3) phi_angle = phi_angle+phi_step;
	
	Transform3D cellTr { RotationZ(phi_angle) };
	
	if (BUILDWIRE[i]) pv = lay_vol.placeVolume(wire_vol[i].at(ilayer),cellTr * wireTr);
      }
    }

    // Create the PhysicalVolume for the layer.
    pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
    lay_elt.setPlacement(pv);
    
    ilayer++;
  }

  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly,Position(0,0,description.constantAsDouble("DCH_z0")));
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);
  return sdet;
}
//-----------------------------------------------------------------------------------//
Volume GetComponentVol(int i,Detector& description,Material gas,xml_comp_t x_comp,bool SHOWSENSOR)
{
  Box    c_box(x_comp.width() / 2, x_comp.length(), x_comp.thickness() / 2);
  Volume c_vol(Form("component%d",i), c_box, gas);

  c_vol.setRegion(description, x_comp.regionStr());
  c_vol.setLimitSet(description, x_comp.limitsStr());

  if (SHOWSENSOR) c_vol.setVisAttributes(description, x_comp.visStr());
  else c_vol.setVisAttributes(description.invisible());

  return c_vol;
}
//-----------------------------------------------------------------------------------//
VolPlane GetSensitiveSurface(int i,Volume vol,PlacedVolume pv,SensitiveDetector sens,double t_inner, double t_outer)
{
  pv.addPhysVolID("sensor", i);
  vol.setSensitiveDetector(sens);
  
  // -------- create a measurement plane for the tracking surface attched to the sensitive volume ----- //
  Vector3D u(-1., 0., 0.);
  Vector3D v(0., -1., 0.);
  Vector3D n(0., 0., 1.);

  SurfaceType type(SurfaceType::Sensitive);
  
  VolPlane surf(vol, type, t_inner, t_outer, u, v, n);

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
  double twist_angle=15*dd4hep::deg;
  return atan( r_z0/Lhalf*tan(twist_angle/2/dd4hep::rad));
}
//-----------------------------------------------------------------------------------// 
double WireLength( double r_z0,double Lhalf,int nwire) 
{
  //auto pitch_z0 = database.at(nlayer).Pitch_z0(r_z0);
  auto pitch_z0 = Pitch_z0(r_z0,nwire); 
  return  2*Lhalf/cos(atan(pitch_z0/(2*Lhalf)))/cos(Stereoangle_z0(r_z0,Lhalf)/dd4hep::rad) ;
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
Volume GetWire(int ilayer,int iwire, double r_z0, double layerLength, int nwire,double rWire, Material mat, VisAttr vis)
{
  double safety_z_interspace=1*dd4hep::nm;

  double rLayer=r_z0;
  if (iwire==F1 || iwire==F4) rLayer=r_z0-rWire*2;
  if (iwire==F3 || iwire==F5) rLayer=r_z0+rWire*2;

  double wlength=0.5*WireLength(rLayer,layerLength,nwire)
    - rWire*cos(Stereoangle_z0(rLayer,layerLength))
    - safety_z_interspace;
  
  Tube wireTub(0,rWire,wlength);
  Volume wireVol(Form("l%d_one_%s",ilayer,tagWire[iwire].c_str()),wireTub,mat);  
  wireVol.setVisAttributes(vis);

  return wireVol;
}
//-----------------------------------------------------------------------------------//
Transform3D GetWireTransform(int ilayer,int iwire,double r_z0,double rWire,double layerLength)
{
  double rLayer=r_z0;
  if (iwire==F1 || iwire==F4) rLayer=r_z0-rWire*2.;
  if (iwire==F3 || iwire==F5) rLayer=r_z0+rWire*2.;
  
  //cout<<"ilayer="<<ilayer<<", layerLength="<<layerLength<<", rlayer="<<rLayer/dd4hep::cm<<" (cm)"<<endl;
  RotationX stereoTr((-1.)*StereoSign(ilayer)*Stereoangle_z0(rLayer,layerLength));
  Transform3D wireTr(stereoTr * Translation3D(rLayer,0.,0.));

  return wireTr;
}
//-----------------------------------------------------------------------------------// 
Assembly GetVesselAssembly(dd4hep::xml::Dimension dimensions,Detector& description,dd4hep::xml::Handle_t vesselParam, bool SHOWVESSEL)
{
  double tShell=vesselParam.attr<double>(_Unicode(tShell));
  double tFill=vesselParam.attr<double>(_Unicode(tFill));
  double tEndcap=2*tShell+tFill;

  //barrel
  Tube outbarrel_tub(dimensions.rmax()-tShell,dimensions.rmax(),(dimensions.length()/2.-tEndcap));
  Volume outbarrel_vol("out_vesselBarrel",outbarrel_tub,description.material(vesselParam.attr<std::string>(_Unicode(shellMat))));
  if (SHOWVESSEL) outbarrel_vol.setVisAttributes(description.visAttributes(vesselParam.attr<std::string>(_Unicode(vis))));
  else outbarrel_vol.setVisAttributes(description.invisible());

  Tube inbarrel_tub(dimensions.rmin(),dimensions.rmin()+tShell,(dimensions.length()/2.-tEndcap));
  Volume inbarrel_vol("in_vesselBarrel",inbarrel_tub,description.material(vesselParam.attr<std::string>(_Unicode(shellMat))));
  if (SHOWVESSEL) inbarrel_vol.setVisAttributes(description.visAttributes(vesselParam.attr<std::string>(_Unicode(vis))));
  else inbarrel_vol.setVisAttributes(description.invisible());

  //endcaps
  Tube endcap_tub(dimensions.rmin(),dimensions.rmax(),tEndcap);
  Volume endcap_vol("vesselEndcap",endcap_tub,description.material(vesselParam.attr<std::string>(_Unicode(shellMat))));
  if (SHOWVESSEL) endcap_vol.setVisAttributes(description.visAttributes(vesselParam.attr<std::string>(_Unicode(vis))));
  else endcap_vol.setVisAttributes(description.invisible());

  Tube filling_tub(dimensions.rmin()+tShell,dimensions.rmax()-tShell,tFill);
  Volume filling_vol("filling",filling_tub,description.material(vesselParam.attr<std::string>(_Unicode(fillMat))));
  if (!SHOWVESSEL) filling_vol.setVisAttributes(description.invisible());
  endcap_vol.placeVolume(filling_vol,Position(0, 0, 0));

  Assembly vessel("vessel");
  vessel.placeVolume(outbarrel_vol,Position(0, 0, 0));
  vessel.placeVolume(inbarrel_vol,Position(0, 0, 0));
  vessel.placeVolume(endcap_vol,Position(0, 0, -dimensions.length()/2.+tEndcap/2.));
  vessel.placeVolume(endcap_vol,Position(0, 0, dimensions.length()/2.-tEndcap/2.));

  return vessel;
}
//-----------------------------------------------------------------------------------//
//@}
// clang-format off
DECLARE_DETELEMENT(D2EIC_mockupDCH, create_MockupDCH)
