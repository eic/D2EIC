//----------------------------------
//         DCH detector v2
//----------------------------------

/*!
 *  \brief     Detector constructor of DCH v2
 *  \details   This code creates full geometry of IDEA DCH subdetector
 *  \author    Alvaro Tolosa-Delgado alvaro.tolosa.delgado@cern.ch
 *  \author    Brieuc Francois       brieuc.francois@cern.ch
 *  \version   2
 *  \date      2024
 *  \pre       DD4hep compiled with Geant4+Qt
 */


#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/Detector.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <array>
#include "DD4hepDetectorHelper.h"
#include "./DriftChamber_info.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

namespace DCH_v2 {

  using DCH_length_t = dd4hep::rec::DCH_info_struct::DCH_length_t;
  using DCH_angle_t  = dd4hep::rec::DCH_info_struct::DCH_angle_t;
  using DCH_layer    = dd4hep::rec::DCH_info_struct::DCH_layer;

  /// Function to build DCH
  static dd4hep::Ref_t create_DCH(dd4hep::Detector &description, dd4hep::xml::Handle_t e, dd4hep::SensitiveDetector sens)
  {
    typedef vector<dd4hep::PlacedVolume> Placements;

    xml_det_t    x_det    = e;
    int          det_id   = x_det.id();
    string       det_name = x_det.nameStr();
    DetElement   sdet(det_name, det_id);
    PlacedVolume pv;
    

    map<string, Volume>                        volumes;
    map<string, Placements>                    sensitives;
    map<string, vector<dd4hep::rec::VolPlane>> volplane_surfaces;
    map<string, array<double, 2>>              module_thicknesses;

    //----------------------------------
    // Set detector type flag
    //----------------------------------
    dd4hep::xml::setDetectorTypeFlag(x_det, sdet);

    //----------------------------------
    // Add the volume boundary material 
    // if configured
    //----------------------------------
    auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);
    xml_comp_t x_boundary_material= x_det.child(_Unicode(boundary_material));
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,"boundary_material");

    /*for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
      xml_comp_t x_boundary_material = bmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,"boundary_material");
      }*/
    
    xml_comp_t x_layer_material = x_det.child(_Unicode(layer_material));
    /*for (xml_coll_t lmat(x_det, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      //DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
      //                                                "layer_material");
      cout<<"lmat="<<lmat<<", x_layer_material: "<<x_layer_material.binning()<<endl;
      }*/

    Assembly assembly(det_name);
    sens.setType("tracker");

    //----------------------------------
    // initialize empty DCH_info object
    //----------------------------------
    // data extension mechanism requires it to be a raw pointer
    dd4hep::rec::DCH_info * DCH_i = new dd4hep::rec::DCH_info();

    // fill DCH_i with information from the XML file
    {
      // DCH outer geometry dimensions
      DCH_i->Set_rin  ( description.constantAsDouble("DCH_gas_inner_cyl_R") );
      DCH_i->Set_rout ( description.constantAsDouble("DCH_gas_outer_cyl_R") );
      DCH_i->Set_lhalf( description.constantAsDouble("DCH_gas_Lhalf")       );

      // guard wires position, fix position
      DCH_i->Set_guard_rin_at_z0  ( description.constantAsDouble("DCH_guard_inner_r_at_z0" ) );
      DCH_i->Set_guard_rout_at_zL2( description.constantAsDouble("DCH_guard_outer_r_at_zL2") );

      DCH_angle_t dch_alpha = description.constantAsDouble("DCH_alpha");
      DCH_i->Set_twist_angle( 2*dch_alpha );

      DCH_i->Set_nsuperlayers        ( description.constantAsLong("DCH_nsuperlayers")         );
      DCH_i->Set_nlayersPerSuperlayer( description.constantAsLong("DCH_nlayersPerSuperlayer") );

      DCH_i->Set_ncell0              ( description.constantAsLong("DCH_ncell")            );
      DCH_i->Set_ncell_increment     ( description.constantAsLong("DCH_ncell_increment")  );
      DCH_i->Set_ncell_per_sector    ( description.constantAsLong("DCH_ncell_per_sector") );

      DCH_i->Set_first_width         ( description.constantAsDouble("DCH_first_width")   );
      DCH_i->Set_first_sense_r       ( description.constantAsDouble("DCH_first_sense_r") );

      bool buildLayers = x_det.attr<bool>(_Unicode(buildLayers));
      if (buildLayers) {
	DCH_i->BuildLayerDatabase();
	// safety check just in case something went wrong...
	if (DCH_i->IsDatabaseEmpty()) throw std::runtime_error("Empty database");
      }

      bool printExcelTable = x_det.attr<bool>(_Unicode(printExcelTable));
      if(printExcelTable) DCH_i->Show_DCH_info_database(std::cout);
    }

    // gas
    bool debugGeometry = x_det.hasChild(_Unicode(debugGeometry));
    auto gasElem   = x_det.child("gas");
    auto gasvolMat = description.material(gasElem.attr<std::string>(_Unicode(material)));
    auto gasvolVis = description.visAttributes(gasElem.attr<std::string>(_Unicode(vis)));

    // vessel
    auto vesselElem                 = x_det.child("vessel");
    auto vesselSkinVis              = description.visAttributes(vesselElem.attr<std::string>(_Unicode(visSkin)));
    auto vesselBulkVis              = description.visAttributes(vesselElem.attr<std::string>(_Unicode(visBulk)));
    auto vessel_mainMaterial        = description.material(vesselElem.attr<std::string>(_Unicode(mainMaterial)));
    auto vessel_fillmaterial_outerR = description.material(vesselElem.attr<std::string>(_Unicode(fillmaterial_outerR)));
    auto vessel_fillmaterial_endcap = description.material(vesselElem.attr<std::string>(_Unicode(fillmaterial_endcap)));
    DCH_length_t vessel_fillmaterial_fraction_outerR = vesselElem.attr<double>(_Unicode(fillmaterial_fraction_outerR)) ;
    DCH_length_t vessel_fillmaterial_fraction_endcap = vesselElem.attr<double>(_Unicode(fillmaterial_fraction_endcap)) ;

    if( 0 > vessel_fillmaterial_fraction_outerR || 1 < vessel_fillmaterial_fraction_outerR )
      throw std::runtime_error("vessel_fillmaterial_fraction_outerR must be between 0 and 1");
    if( 0 > vessel_fillmaterial_fraction_endcap || 1 < vessel_fillmaterial_fraction_endcap )
      throw std::runtime_error("vessel_fillmaterial_fraction_z must be between 0 and 1");

    // more wire stuff
    auto wiresElem = x_det.child("wires");
    //auto wiresVis  = description.visAttributes(wiresElem.attr<std::string>(_Unicode(vis)));
    bool buildSenseWires  = wiresElem.attr<bool>(_Unicode(buildSenseWires));
    bool buildFieldWires  = wiresElem.attr<bool>(_Unicode(buildFieldWires));

    DCH_length_t dch_SWire_thickness        = wiresElem.attr<double>(_Unicode(SWire_thickness)) ;
    DCH_length_t dch_FSideWire_thickness    = wiresElem.attr<double>(_Unicode(FSideWire_thickness)) ;
    DCH_length_t dch_FCentralWire_thickness = wiresElem.attr<double>(_Unicode(FCentralWire_thickness)) ;

    auto dch_SWire_material        = description.material(wiresElem.attr<std::string>(_Unicode(SWire_material)) ) ;
    auto dch_FSideWire_material    = description.material(wiresElem.attr<std::string>(_Unicode(FSideWire_material)) ) ;
    auto dch_FCentralWire_material = description.material(wiresElem.attr<std::string>(_Unicode(FCentralWire_material)) ) ;

    /* Geometry tree:
     * Gas (tube) -> Layer_1 (hyp) -> cell_1 (twisted tube)
     *                             -> cell_... (twisted tube)
     *            -> Layer_... (hyp) -> cell_1 (twisted tube)
     *                               -> cell_... (twisted tube)
     *            -> Inner radius vessel wall
     *            -> Outer radius vessel wall -> fill made of foam
     *            -> Endcap disk  vessel wall -> fill made of kapton
     *
     * Layers represent a segmentation in radius
     * Sectors represent a segmentation in phi
     * Each cell corresponds to a Detector Element
     * Vessel wall has to be defined as 3 volumes to account for independent thickness and materials
     */

    DCH_length_t safety_r_interspace   = 1    * dd4hep::nm;
    DCH_length_t safety_z_interspace   = 1    * dd4hep::nm;
    DCH_length_t safety_phi_interspace = 1e-6 * dd4hep::rad;

    DCH_length_t vessel_thickness_innerR  = description.constantAsDouble("DCH_vessel_thickness_innerR");
    DCH_length_t vessel_thickness_outerR  = description.constantAsDouble("DCH_vessel_thickness_outerR");
    DCH_length_t vessel_endcapdisk_zmin   = description.constantAsDouble("DCH_vessel_disk_zmin");
    DCH_length_t vessel_endcapdisk_zmax   = description.constantAsDouble("DCH_vessel_disk_zmax");

    DCH_length_t z0= description.constantAsDouble("DCH_z0");
    // if( 0 > vessel_thickness_z )
    // throw std::runtime_error("vessel_thickness_z must be positive");
    if( 0 > vessel_thickness_innerR )
      throw std::runtime_error("vessel_thickness_innerR must be positive");
    if( 0 > vessel_thickness_outerR )
      throw std::runtime_error("vessel_thickness_outerR must be positive");

    //----------------------------------
    // build gas volume
    //----------------------------------
    dd4hep::Tube gas_s   ( DCH_i->rin   - vessel_thickness_innerR,
                           DCH_i->rout  + vessel_thickness_outerR,
                           vessel_endcapdisk_zmax );
    dd4hep::Volume gas_v ( det_name+"_gas", gas_s, gasvolMat );

    gas_v.setVisAttributes( gasvolVis );
    gas_v.setRegion  ( description, x_det.regionStr() );
    gas_v.setLimitSet( description, x_det.limitsStr() );
    gas_v.setSensitiveDetector(sens);
    pv = assembly.placeVolume(gas_v);

    dd4hep::DetElement gas_DE(sdet,"gas",0);
    gas_DE.setPlacement(pv);
    
    //----------------------------------
    // build vessel
    //----------------------------------
    
    DCH_length_t vessel_innerR_start = DCH_i->rin - vessel_thickness_innerR + safety_r_interspace;
    DCH_length_t vessel_innerR_end   = DCH_i->rin   ;
    DCH_length_t vessel_outerR_start = DCH_i->rout;
    DCH_length_t vessel_outerR_end   = DCH_i->rout  + vessel_thickness_outerR - safety_r_interspace;
    DCH_length_t vessel_R_zhalf      = vessel_endcapdisk_zmax - safety_z_interspace;

    // vessel: inner R wall 
    dd4hep::Tube   vessel_innerR_s  ( vessel_innerR_start, vessel_innerR_end, vessel_R_zhalf );
    dd4hep::Volume vessel_innerR_v  ( det_name+"_vessel_innerR", vessel_innerR_s,  vessel_mainMaterial );
    vessel_innerR_v.setVisAttributes( vesselSkinVis );
    gas_v.placeVolume( vessel_innerR_v);

    // vessel: outer R wall
    dd4hep::Tube   vessel_outerR_s  ( vessel_outerR_start, vessel_outerR_end, vessel_R_zhalf );
    dd4hep::Volume vessel_outerR_v  ( det_name+"_vessel_outerR", vessel_outerR_s,  vessel_mainMaterial );
    vessel_outerR_v.setVisAttributes( vesselSkinVis );

    // if thickness fraction of bulk material is defined, build the bulk material
    if(0 < vessel_fillmaterial_fraction_outerR) {
      double       f                      = vessel_fillmaterial_fraction_outerR;
      DCH_length_t fillmaterial_thickness = f * (vessel_outerR_end - vessel_outerR_start);
      DCH_length_t rstart                 = vessel_outerR_start + 0.5*(1-f)*fillmaterial_thickness;
      DCH_length_t rend                   = vessel_outerR_end   - 0.5*(1-f)*fillmaterial_thickness;
      dd4hep::Tube   vessel_fillmat_outerR_s ( rstart, rend, vessel_R_zhalf - safety_z_interspace );
      dd4hep::Volume vessel_fillmat_outerR_v ( det_name+"_vessel_fillmat_outerR", 
					       vessel_fillmat_outerR_s,
					       vessel_fillmaterial_outerR );
      vessel_fillmat_outerR_v.setVisAttributes( vesselBulkVis );
      vessel_outerR_v.placeVolume(vessel_fillmat_outerR_v);
    }

    gas_v.placeVolume(vessel_outerR_v);

    // vessel: endcap walls
    DCH_length_t vessel_endcap_thickness = vessel_endcapdisk_zmax - vessel_endcapdisk_zmin - 2*safety_z_interspace;
    DCH_length_t vessel_endcap_zpos      = 0.5*(vessel_endcapdisk_zmax + vessel_endcapdisk_zmin);
    DCH_length_t vessel_endcap_rstart    = vessel_innerR_end   + safety_r_interspace;
    DCH_length_t vessel_endcap_rend      = vessel_outerR_start - safety_r_interspace;
    dd4hep::Tube vessel_endcap_s    ( vessel_endcap_rstart, vessel_endcap_rend, 0.5*vessel_endcap_thickness );
    dd4hep::Volume vessel_endcap_v  ( det_name+"_vessel_endcap", vessel_endcap_s, vessel_mainMaterial );
    vessel_endcap_v.setVisAttributes( vesselSkinVis );

    // if thickness fraction of bulk material is defined, build the bulk material
    if(0 < vessel_fillmaterial_fraction_endcap) {
      double f = vessel_fillmaterial_fraction_endcap;
      DCH_length_t fillmaterial_thickness = f * vessel_endcap_thickness;
      dd4hep::Tube vessel_fillmat_endcap_s    ( vessel_endcap_rstart + safety_r_interspace,
						vessel_endcap_rend   - safety_r_interspace  ,
						0.5*fillmaterial_thickness );
      dd4hep::Volume vessel_fillmat_endcap_v  ( det_name+"_vessel_fillmat_endcap",
						vessel_fillmat_endcap_s,
						vessel_fillmaterial_endcap );
      vessel_fillmat_endcap_v.setVisAttributes( vesselBulkVis );
      vessel_endcap_v.placeVolume(vessel_fillmat_endcap_v);
    }
    // place endcap wall at +/- z
    gas_v.placeVolume(vessel_endcap_v, dd4hep::Position(0,0, vessel_endcap_zpos));
    gas_v.placeVolume(vessel_endcap_v, dd4hep::Position(0,0,-vessel_endcap_zpos));
    
    
    //---------------------------------- 
    // build gas layers 
    //---------------------------------- 
    int cnt=0;
    for(const auto& [ilayer, l]  : DCH_i->database ) {
      if (ilayer>10) continue;
      //----------------------------------
      // INITIALIZATION OF THE LAYER
      //----------------------------------

      // Hyperboloid parameters:
      /// inner radius at z=0
      DCH_length_t rin   = l.radius_fdw_z0+safety_r_interspace;
      /// inner stereoangle, calculated from rin(z=0)
      //DCH_angle_t  stin  = DCH_i->stereoangle_z0(rin);
      /// outer radius at z=0
      DCH_length_t rout  = l.radius_fuw_z0-safety_r_interspace;
      /// outer stereoangle, calculated from rout(z=0)
      //DCH_angle_t  stout = DCH_i->stereoangle_z0(rout);
      /// half-length
      DCH_length_t dz    = DCH_i->Lhalf + safety_z_interspace;
      
      DCH_length_t maxGap=50*dd4hep::mm;
      
      //DCH_length_t myrout=(rin+dz*tan(stin)/2.);
      //if (rin+dz*tan(stin)/2. > DCH_i->rout) break;
      if (rout+maxGap > DCH_i->rout) break;
      /*std::cout<<"cnt="<<cnt
	       <<", ncell="<< l.nwires/2
	       <<", rout="<<rout/dd4hep::mm <<"(mm) , max r="
	       <<DCH_i->rout/dd4hep::mm <<"(mm)"<<std::endl;
      */

      // !!!!!!!!!!!!!!!!!!!!!!!!!!! 
      //if (ilayer!=112 && ilayer!=111) continue;
      // !!!!!!!!!!!!!!!!!!!!!!!!!!! 

      //dd4hep::Hyperboloid layer_s(rin, stin, rout, stout, dz);
      dd4hep::Tube layer_s(rin, rout, dz);

      string layer_name = "DCH_layer"+std::to_string(ilayer);  //need to match layer_pattern in xml
      dd4hep::Volume layer_v ( layer_name , layer_s, gasvolMat );
      layer_v.setVisAttributes( description.visAttributes( Form("dch_layer_vis%d", ilayer%22) ) );
      layer_v.setSensitiveDetector(sens);
      
      auto layer_pv = gas_v.placeVolume(layer_v);
      layer_pv.addPhysVolID("layer",ilayer);
      
      // ilayer is a counter that runs from 1 to 112 (nsuperlayers * nlayersPerSuperlayer)
      // it seems more convenient to store the layer number within the superlayer         
      // ilayerWithinSuperlayer runs from 0 to 7 (nlayersPerSuperlayer-1)  
      //int ilayerWithinSuperlayer = (ilayer-1) % DCH_i->nlayersPerSuperlayer;
      //layer_pv.addPhysVolID("layer", ilayerWithinSuperlayer  );

      // add superlayer bitfield
      //int nsuperlayer_minus_1 = DCH_i->Get_nsuperlayer_minus_1(ilayer);
      //layer_pv.addPhysVolID("superlayer", nsuperlayer_minus_1 ); 


      //----------------------------------
      // test ACTS stuff
      //----------------------------------
      sensitives[layer_name].push_back(layer_pv);
      module_thicknesses[layer_name] = {rin, rout};

      dd4hep::DetElement layer_DE(gas_DE,layer_name, ilayer); 
      cout<<"layer_DE: "<<layer_DE.parent().parent().name()<<"/"<<layer_DE.parent().name()<<"/"<<layer_DE.name()<<endl;
      layer_DE.setPlacement(layer_pv);

      auto& layer_DE_params =	DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(layer_DE);
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layer_DE_params, "layer_material");
      /*layer_DE_params.set<string>("axis_definitions", "XYZ");
      layer_DE_params.set<double>("envelope_r_min", rin);
      layer_DE_params.set<double>("envelope_r_max", rout);
      layer_DE_params.set<double>("envelope_z_min", -dz);
      layer_DE_params.set<double>("envelope_z_max", dz);*/
      cout<<"layer param: "<<layer_DE_params<<endl; 

      //----------------------------------
      // create a measurement plane for 
      // the tracking surface attched to 
      // the sensitive volume
      //----------------------------------
      Vector3D u(-1., 0., 0.);
      Vector3D v(0., -1., 0.);
      Vector3D n(0., 0., 1.);

      //----------------------------------
      // add surface
      //----------------------------------
      SurfaceType type(rec::SurfaceType::Sensitive);
      VolPlane surf(layer_v, type, module_thicknesses[layer_name][0], module_thicknesses[layer_name][1], u, v, n);
      volplane_surfaces[layer_name].push_back(surf);
      volSurfaceList(layer_DE)->push_back(volplane_surfaces[layer_name][0]);
      if (ilayer==0) cout<<"volplane_surfaces["<<layer_name<<"]= "<<volplane_surfaces[layer_name][0]<<endl;

      //---------------------------------- 
      // SEGMENTATION OF THE LAYER
      // INTO CELLS (TWISTED TUBES)
      // !!!! NOT USED !!!!
      //---------------------------------- 
      // ncells in this layer = 2x number of wires
      int ncells = l.nwires/2;
      DCH_angle_t phi_step = (TMath::TwoPi()/ncells)*dd4hep::rad;

      // unitary cell (Twisted tube) is repeated for each layer l.nwires/2 times
      // Twisted tube parameters
      DCH_angle_t cell_twistangle    = l.StereoSign() * DCH_i->twist_angle;
      DCH_length_t cell_rin_z0       = l.radius_fdw_z0 + 2*safety_r_interspace;
      DCH_length_t cell_rout_z0      = l.radius_fuw_z0 - 2*safety_r_interspace;
      DCH_length_t cell_rin_zLhalf   = DCH_i->Radius_zLhalf(cell_rin_z0);
      DCH_length_t cell_rout_zLhalf  = DCH_i->Radius_zLhalf(cell_rout_z0);
      DCH_length_t cell_dz           = DCH_i->Lhalf;
      DCH_angle_t cell_phi_width     = phi_step - safety_phi_interspace;
      dd4hep::TwistedTube cell_s( cell_twistangle, cell_rin_zLhalf, cell_rout_zLhalf, cell_dz, 1, cell_phi_width);

      // initialize cell volume
      std::string cell_name = det_name+"_layer"+std::to_string(ilayer)+"_cell";
      
      //---------------------------------- 
      // Single sense wire
      //---------------------------------- 
      // average radius to position sense wire
      DCH_length_t cell_rave_z0 = 0.5*(cell_rin_z0+cell_rout_z0);
      DCH_length_t cell_swire_radius = dch_SWire_thickness/2;
      DCH_length_t swlength = 0.5*DCH_i->WireLength(ilayer,cell_rave_z0)
	- cell_swire_radius*cos(DCH_i->stereoangle_z0(cell_rave_z0))
	- safety_z_interspace;
      
      dd4hep::Tube swire_s(0., dch_SWire_thickness, swlength);
      dd4hep::Volume swire_v(cell_name+"_swire", swire_s, dch_SWire_material);
      swire_v.setVisAttributes(description.visAttributes( Form("dch_vis_cell_%d",(int) ilayer%2 )));

      // Change sign of stereo angle to place properly the wire inside the twisted tube      
      dd4hep::RotationX stereoTr( (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(cell_rave_z0) );
      //dd4hep::Transform3D swireTr ( stereoTr * dd4hep::Translation3D(cell_rave_z0,0.,0.) );
      
      //----------------------------------
      // Single field wire 
      //----------------------------------
      
      // POSITIONING OF FIELD WIRES
      //
      //  The following sketch represents the crossection of a DCH cell, where
      //      O symbol = Field wires, the number in parenthesis is used as ID 
      //      X symbol = sense wire                                           
      //
      //   ^ radius
      //   O(1)---O(4)---O(6)    radius_z0 = l.radius_fuw_z0 == (++l).radius_fdw_z0
      //   O(2)   X      O(7)    radius_z0 = average(l.radius_fuw_z0, l.radius_fdw_z0)
      //   O(3)---O(5)---O(8)    radius_z0 = l.radius_fdw_z0 == (--l).radius_fuw_z0   
      //
      //   --> phi axis
      //
      //  In the previous sketch, the wires are shared among several cells.
      //  Since we are using an actual shape to contain each cell,
      //  it is not feasible.
      //
      //   O(1)---O(4)---    radius_z0 = l.radius_fuw_z0 - wire_thickness/2
      //   O(2)   X          radius_z0 = average(l.radius_fuw_z0, l.radius_fdw_z0)                                                                                                   
      //   O(3)---O(5)---    radius_z0 = l.radius_fdw_z0 + wire_thickness/2 
      //
      //  encapsulate the calculation of the phi offset into a function
      //  since it will be different for each field wire
      //  it includes the safety phi distance
      
      auto fwire_phi_offset = [&](DCH_length_t radial_distance, DCH_length_t wire_radius)->DCH_angle_t {
	return atan(wire_radius/radial_distance)*dd4hep::rad + safety_phi_interspace;
      };
      
      //----------------------------------  
      // fwire 2
      //----------------------------------  
      DCH_length_t fwire_radius = dch_FCentralWire_thickness/2;
      DCH_length_t fwire_r_z0   = cell_rave_z0;
      DCH_angle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(fwire_r_z0);
      DCH_angle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
      DCH_length_t fwire_length = 0.5*DCH_i->WireLength(ilayer, fwire_r_z0)
	- fwire_radius*cos(DCH_i->stereoangle_z0(fwire_r_z0))
	- safety_z_interspace;

      dd4hep::Tube fwire2_s(0., fwire_radius, fwire_length);
      dd4hep::Volume fwire2_v(cell_name+"_f2wire", fwire2_s, dch_FCentralWire_material );
      fwire2_v.setVisAttributes(description.visAttributes( Form("dch_vis_cell_%d",(int) ilayer%2 )));

      // Change sign of stereo angle to place properly the wire inside the twisted tube                                                                                                     
      dd4hep::RotationX fwire2StereoTr( fwire_stereo );
      dd4hep::RotationZ fwire2PhoTr( fwire_phi );
      dd4hep::Transform3D fwire2Tr ( fwire2PhoTr * fwire2StereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );

      //----------------------------------
      // fwire 1
      //---------------------------------- 
      fwire_radius = dch_FSideWire_thickness/2;
      // decrease radial distance, move it closer to the sense wire                                                                                                                         
      fwire_r_z0   = cell_rout_z0 - fwire_radius;
      fwire_stereo =  (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(fwire_r_z0);
      fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
      fwire_length = 0.5*DCH_i->WireLength(ilayer, fwire_r_z0)
	- fwire_radius*cos(DCH_i->stereoangle_z0(fwire_r_z0))
	- safety_z_interspace;

      dd4hep::Tube fwire1_s(0., fwire_radius, fwire_length);
      dd4hep::Volume fwire1_v(cell_name+"_f1wire", fwire1_s, dch_FSideWire_material );
      fwire1_v.setVisAttributes(description.visAttributes( Form("dch_vis_cell_%d",(int) ilayer%2 )));

      // Change sign of stereo angle to place properly the wire inside the twisted tube                                                                                                     
      dd4hep::RotationX fwire1StereoTr( fwire_stereo );
      dd4hep::RotationZ fwire1PhoTr( fwire_phi );
      dd4hep::Transform3D fwire1Tr ( fwire1PhoTr * fwire1StereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );

      //----------------------------------
      // fwire 3
      //----------------------------------
      //same radius as fwire 1
      fwire_r_z0   = cell_rin_z0 + fwire_radius;
      fwire_stereo =  (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(fwire_r_z0);
      fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
      fwire_length = 0.5*DCH_i->WireLength(ilayer, fwire_r_z0)
	- fwire_radius*cos(DCH_i->stereoangle_z0(fwire_r_z0))
	- safety_z_interspace;

      dd4hep::Tube fwire3_s(0., fwire_radius, fwire_length);
      dd4hep::Volume fwire3_v(cell_name+"_f3wire", fwire3_s, dch_FSideWire_material );
      fwire3_v.setVisAttributes(description.visAttributes( Form("dch_vis_cell_%d",(int) ilayer%2 )));

      // Change sign of stereo angle to place properly the wire inside the twisted tube                                                                                                     
      dd4hep::RotationX fwire3StereoTr( fwire_stereo );
      dd4hep::RotationZ fwire3PhoTr( fwire_phi );
      dd4hep::Transform3D fwire3Tr ( fwire3PhoTr * fwire3StereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );

      //----------------------------------
      // fwire 5
      //----------------------------------
      //same parameter as fwire 3, just no rotation z.
      dd4hep::Tube fwire5_s(0., fwire_radius, fwire_length);
      dd4hep::Volume fwire5_v(cell_name+"_f5wire", fwire5_s, dch_FSideWire_material );
      fwire5_v.setVisAttributes(description.visAttributes( Form("dch_vis_cell_%d", (int) ilayer%2 )));

      // Change sign of stereo angle to place properly the wire inside the twisted tube                                                                                                     
      dd4hep::RotationX fwire5StereoTr( fwire_stereo );
      dd4hep::Transform3D fwire5Tr ( fwire5StereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );

      //----------------------------------
      // fwire 4
      //----------------------------------
      //same radio as fwire 1, 3, and 5
      fwire_r_z0   = cell_rout_z0 - fwire_radius;
      fwire_stereo =  (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(fwire_r_z0);
      fwire_length = 0.5*DCH_i->WireLength(ilayer, fwire_r_z0)
	- fwire_radius*cos(DCH_i->stereoangle_z0(fwire_r_z0))
	- safety_z_interspace;

      dd4hep::Tube fwire4_s(0., fwire_radius, fwire_length);
      dd4hep::Volume fwire4_v(cell_name+"_f4wire", fwire4_s, dch_FSideWire_material );
      fwire4_v.setVisAttributes(description.visAttributes( Form("dch_vis_cell_%d",(int) ilayer%2 )));

      // Change sign of stereo angle to place properly the wire inside the twisted tube                                                                                                     
      dd4hep::RotationX fwire4StereoTr( fwire_stereo );
      dd4hep::Transform3D fwire4Tr ( fwire4StereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );
      

      //----------------------------------  
      // put wires in layer
      //----------------------------------    
      int maxphi = ncells;
      if(debugGeometry) maxphi=3;
      //maxphi=20;
      
      for(int nphi = 0; nphi < maxphi; ++nphi) {
	DCH_angle_t cell_phi_angle = phi_step * nphi + 0.25*cell_phi_width*(ilayer%2);
	dd4hep::Transform3D cellTr { dd4hep::RotationZ(cell_phi_angle) };
	
	if(buildSenseWires) {
	  auto swire_pv = layer_v.placeVolume(swire_v, cellTr* stereoTr * dd4hep::Translation3D(cell_rave_z0,0.,0.)); 
	  swire_pv.addPhysVolID("nphi", nphi);                   
	  swire_pv.addPhysVolID("stereosign", l.StereoSign() );  
	}
	  
	if(buildFieldWires) {
	  // f wire 2
	  auto fware2_pv=layer_v.placeVolume(fwire2_v,cellTr*fwire2Tr);
	  fware2_pv.addPhysVolID("nphi", nphi);
	  fware2_pv.addPhysVolID("stereosign", l.StereoSign() );

	  // f wire 1
	  auto fware1_pv=layer_v.placeVolume(fwire1_v,cellTr*fwire1Tr);
	  fware1_pv.addPhysVolID("nphi", nphi);
	  fware1_pv.addPhysVolID("stereosign", l.StereoSign() );
	  
	  // f wire 3
	  auto fware3_pv=layer_v.placeVolume(fwire3_v,cellTr*fwire3Tr);
	  fware3_pv.addPhysVolID("nphi", nphi);
	  fware3_pv.addPhysVolID("stereosign", l.StereoSign() );

	  // f wire 5
	  auto fware5_pv=layer_v.placeVolume(fwire5_v,cellTr*fwire5Tr);
	  fware5_pv.addPhysVolID("nphi", nphi);
	  fware5_pv.addPhysVolID("stereosign", l.StereoSign() );
	  
	  // f wire 4
	  auto fware4_pv=layer_v.placeVolume(fwire4_v,cellTr*fwire4Tr);
	  fware4_pv.addPhysVolID("nphi", nphi);
	  fware4_pv.addPhysVolID("stereosign", l.StereoSign() );
	}// end building field wires
      }// end building wires
      cnt++;
    }// end building layers
    std::cout<<"DCH: number of layers= "<<cnt<<std::endl;

    // Place our mother volume in the world
    sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
    assembly.setVisAttributes(description.invisible());
    pv = description.pickMotherVolume(sdet).placeVolume(assembly,dd4hep::Position(0,0,z0));
    pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
    sdet.setPlacement(pv);
    return sdet;
  }
  
}; // end DCH_v2 namespace


DECLARE_DETELEMENT(DriftChamber, DCH_v2::create_DCH)
