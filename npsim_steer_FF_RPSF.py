#!/usr/bin/env python
"""
DD4hep simulation with some argument parsing
Based on M. Frank and F. Gaede runSim.py
   @author  A.Sailer
   @version 0.1
Modified with standard EIC EPIC requirements.
"""
################################################################
# Obtained steering file from Alex Jentsch
# Modified to use for IP-8 Far-forward detectors Jihee Kim
################################################################
from __future__ import absolute_import, unicode_literals
import logging
import sys
import math

from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV, radian

SIM = DD4hepSimulation()
#############################################################
#   Particle Gun simulations use these options
#   --> comment out if using MC input
#############################################################
SIM.enableG4GPS = False
SIM.enableG4Gun = False
SIM.enableGun = True
SIM.crossingAngleBoost = 0.035*radian
SIM.inputFiles = []
SIM.macroFile = ""
SIM.runType = "run"
SIM.skipNEvents = 0
SIM.vertexOffset = [0.0, 0.0, 0.0, 0.0]
SIM.vertexSigma = [0.0, 0.0, 0.0, 0.0]
################################################################################
## Configuration for the magnetic field (stepper)
################################################################################
SIM.field.delta_chord = 1e-03
SIM.field.delta_intersection = 1e-03
SIM.field.delta_one_step = .5e-01*mm
SIM.field.eps_max = 1e-03
SIM.field.eps_min = 1e-04
SIM.field.largest_step = 100.0*mm
SIM.field.min_chord_step = 1.e-2*mm
SIM.field.stepper = "HelixSimpleRunge"
################################################################################
## Configuration for the DDG4 ParticleGun
################################################################################
SIM.gun.distribution = 'uniform'
SIM.gun.energy = 275.0*GeV
SIM.gun.isotrop = True
SIM.gun.particle = "proton"
SIM.gun.position = (0.*mm,0.*mm,0.*mm)
SIM.gun.thetaMin = 0.000*radian
SIM.gun.thetaMax = 0.005*radian
################################################################################
## Configuration for the output levels of DDG4 components
################################################################################
SIM.output.inputStage = 3
SIM.output.kernel = 3
SIM.output.part = 3
SIM.output.random = 6
################################################################################
## Configuration for the Particle Handler/ MCTruth treatment
################################################################################
SIM.part.enableDetailedHitsAndParticleInfo = False
SIM.part.keepAllParticles = True
SIM.part.minDistToParentVertex = 2.2e-14
SIM.part.printEndTracking = False
SIM.part.printStartTracking = False
SIM.part.saveProcesses = ['Decay']
################################################################################
## Configuration for the PhysicsList
################################################################################
SIM.physics.decays = True
SIM.physics.list = "FTFP_BERT"
################################################################################
## Properties for the random number generator
################################################################################
SIM.random.enableEventSeed = False
SIM.random.file = None
SIM.random.luxury = 1
SIM.random.replace_gRandom = True
SIM.random.seed = None
SIM.random.type = None 
