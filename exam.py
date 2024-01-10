#!/usr/bin/env python

import sys
from geant4_pybind import *
import math 

class ExamDetectorConstruction(G4VUserDetectorConstruction):

  def __init__(self):
     super().__init__()
     self.fScoringVolume = None
 
  def Construct(self):
     nist = G4NistManager.Instance()

     envelop_x = 30*cm
     envelop_y = 30*cm
     envelop_z = 30*cm

     envelop_mat = nist.FindOrBuildMaterial("G4_AIR")

     box_x = 1.2*envelop_x
     box_y = 1.2*envelop_y
     box_z = 1.2*envelop_z

     zTrans = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0.1*envelop_x, 0, 0.05*envelop_z))

     mat_hip1 = nist.FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP")
      
     mat_hip2 = nist.FindOrBuildMaterial("G4_Ti")

     checkOverlaps = True

     sWorld = G4Box("World", 0.5*envelop_x, 0.5*envelop_y,  
                    0.5*envelop_z)
 
     lWorld = G4LogicalVolume(sWorld, envelop_mat, "World")
 
     pWorld = G4PVPlacement(None, G4ThreeVector(), 
                            lWorld, "World", None, False,
                            0, checkOverlaps)

     sBox = G4Box("Box", 0.5*box_x, 0.5*box_y,
                  0.5*box_z)
     lBox = G4LogicalVolume(sBox, envelop_mat, "Box")
     G4PVPlacement(None, G4ThreeVector(),
                           lBox, "Box", lWorld, False,
                           0, checkOverlaps)

     sHip1 = G4Tubs("Hip", 0, 0.3*envelop_x, 0.5*envelop_y, 2*math.pi, 2*math.pi)
     sHip2 = G4Tubs("Prothesis", 0, 0.05*envelop_x, 0.5*envelop_y, 2*math.pi, 2*math.pi)
     sCut = G4SubtractionSolid("Cut", sHip2, sHip1, zTrans)

     lHip1 = G4LogicalVolume(sHip1, mat_hip1, "Hip")
     lHip2 = G4LogicalVolume(sHip2, mat_hip2, "Prothesis")

     G4PVPlacement(None, G4ThreeVector(), lHip1,
                   "Hip", lBox, True, 0, checkOverlaps)

     G4PVPlacement(None, G4ThreeVector(0.1*envelop_x, 0.05*envelop_y, 0), lHip2,
                   "Prothesis", lHip1, True, 0, checkOverlaps)

     self.fScoringVolume = lBox

     return pWorld


class ExamActionInitialization(G4VUserActionInitialization):
  """
  Initialization of user code.
  """
  def BuildForMaster(self):
    self.SetUserAction(ExamRunAction())

  def Build(self):
    self.SetUserAction(ExamPrimaryGeneratorAction())

    runAction = ExamRunAction()
    self.SetUserAction(runAction)

    eventAction = ExamEventAction(runAction)
    self.SetUserAction(eventAction)

class ExamPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
    def __init__(self):
        super().__init__()
        self.fEnvelopeBox = None
        self.fParticleGun = G4ParticleGun(1)

        particleTable = G4ParticleTable.GetParticleTable()
        particle = particleTable.FindParticle("gamma")
        self.fParticleGun.SetParticleDefinition(particle)
        self.fParticleGun.SetParticleMomentumDirection(G4ThreeVector(0, 0, 1))
        self.fParticleGun.SetParticleEnergy(1*MeV)

    def GeneratePrimaries(self, anEvent):
        envSizeX = 0
        envSizeY = 0
        envSizeZ = 0
    
        if self.fEnvelopeBox == None:
            envLV = G4LogicalVolumeStore.GetInstance().GetVolume("Box")

            if envLV != None:
                self.fEnvelopeBox = envLV.GetSolid()
          
            if self.fEnvelopeBox != None:
                envSizeX = self.fEnvelopeBox.GetXHalfLength()*2
                envSizeY = self.fEnvelopeBox.GetYHalfLength()*2
                envSizeZ = self.fEnvelopeBox.GetZHalfLength()*2
            else:
                msg = "Envelope volume of box shape not found.\n"
                msg += "Perhaps you have changed geometry.\n"
                msg += "The gun will be place at the center."
                G4Exception("ExamPrimaryGeneratorAction::GeneratePrimaries()",
                            "MyCode0002", G4ExceptionSeverity.JustWarning, msg)

            #size = 0.3
            x0 = envSizeX * (G4UniformRand() - 0.5)
            y0 = 0
            z0 = -0.5     
    
            self.fParticleGun.SetParticlePosition(G4ThreeVector(x0, y0, z0))
            self.fParticleGun.GeneratePrimaryVertex(anEvent)

class ExamRunAction(G4UserRunAction):
    def __init__(self):
     super().__init__()
  
     milligray = 1.e-3*gray
     microgray = 1.e-6*gray
     nanogray = 1.e-9*gray
     picogray = 1.e-12*gray
  
     G4UnitDefinition("milligray", "milliGy", "Dose", milligray)
     G4UnitDefinition("microgray", "microGy", "Dose", microgray)
     G4UnitDefinition("nanogray", "nanoGy", "Dose", nanogray)
     G4UnitDefinition("picogray", "picoGy", "Dose", picogray)
 
     self.edep = G4Accumulable(0)
     self.edep2 = G4Accumulable(0)
 
     accumulableManager = G4AccumulableManager.Instance()
     accumulableManager.RegisterAccumulable(self.edep)
     accumulableManager.RegisterAccumulable(self.edep2)
 
    def BeginOfRunAction(self, aRun):
     G4RunManager.GetRunManager().SetRandomNumberStore(False)
 
     accumulableManager = G4AccumulableManager.Instance()
     accumulableManager.Reset()
 
    def EndOfRunAction(self, aRun):
     nofEvents = aRun.GetNumberOfEvent()
     if nofEvents == 0:
       return
 
     # Merge accumulables
     accumulableManager = G4AccumulableManager.Instance()
     accumulableManager.Merge()
 
     edep = self.edep.GetValue()
     edep2 = self.edep2.GetValue()
 
     # Compute dose = total energy deposit in a run and its variance
     rms = edep2 - edep*edep/nofEvents
     if rms > 0:
       rms = math.sqrt(rms)
     else:
       rms = 0
 
     detectorConstruction = G4RunManager.GetRunManager().GetUserDetectorConstruction()
     mass = detectorConstruction.fScoringVolume.GetMass()
     dose = edep/mass
     rmsDose = rms/mass
 
     generatorAction = G4RunManager.GetRunManager().GetUserPrimaryGeneratorAction()
     runCondition = ""
     if generatorAction != None and isinstance(generatorAction, ExamPrimaryGeneratorAction):
       particleGun = generatorAction.fParticleGun
       runCondition += particleGun.GetParticleDefinition().GetParticleName() + "(s)"
       runCondition += " of "
       particleEnergy = particleGun.GetParticleEnergy()
       runCondition += "{:.5g}".format(G4BestUnit(particleEnergy, "Energy"))
 
       if self.IsMaster():
         print("--------------------End of Global Run-----------------------")
       else:
         print("--------------------End of Local Run------------------------")
 
     print(" The run consists of", nofEvents, runCondition)
     print(" Cumulated dose per run, in scoring volume: ", end="")
     print("{:.5f} rms = {:.5f}".format(G4BestUnit(dose, "Dose"), G4BestUnit(rmsDose, "Dose")))
     print("------------------------------------------------------------")
     print("")
 
    def AddEdep(self, edep):
     self.edep += edep
     self.edep2 += edep*edep


class ExamEventAction(G4UserEventAction):
  def __init__(self, runAction):
    super().__init__()
    self.fRunAction = runAction

  def BeginOfEventAction(self, anEvent):
    self.fEdep = 0
 
  def EndOfEventAction(self, anEvent):
    self.fRunAction.AddEdep(self.fEdep)

  def AddEdep(self, edep):
    self.fEdep += edep

class ExamSteppingAction(G4UserSteppingAction):
  def __init__(self, eventAction):
    super().__init__()
    self.fEventAction = eventAction
    self.fScoringVolume = None
 
  def UserSteppingAction(self, aStep):
    if self.fScoringVolume == None:
      detectorConstruction = G4RunManager.GetRunManager().GetUserDetectorConstruction()
      self.fScoringVolume = detectorConstruction.fScoringVolume
 
    volume = aStep.GetPreStepPoint().GetTouchable().GetVolume().GetLogicalVolume()
 
    # check if we are in scoring volume
    if volume != self.fScoringVolume:
      return
 
    # collect energy deposited in this step
    edepStep = aStep.GetTotalEnergyDeposit()
    self.fEventAction.AddEdep(edepStep)




ui = None
if len(sys.argv) == 1:
  ui = G4UIExecutive(len(sys.argv), sys.argv)

# Optionally: choose a different Random engine...
# G4Random.setTheEngine(MTwistEngine())

runManager = G4RunManagerFactory.CreateRunManager(G4RunManagerType.Serial)

runManager.SetUserInitialization(ExamDetectorConstruction())

# Physics list
physicsList = QBBC()
physicsList.SetVerboseLevel(1)

runManager.SetUserInitialization(physicsList)

# User action initialization
runManager.SetUserInitialization(ExamActionInitialization())

visManager = G4VisExecutive()
# G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
# visManager = G4VisExecutive("Quiet")
visManager.Initialize()

# Get the User Interface manager
UImanager = G4UImanager.GetUIpointer()

# # Process macro or start UI session
if ui == None:
  # batch mode
  command = "/control/execute "
  fileName = sys.argv[1]
  UImanager.ApplyCommand(command + fileName)
else:
  # interactive mode
  UImanager.ApplyCommand("/control/execute init_vis.mac")
  ui.SessionStart()
