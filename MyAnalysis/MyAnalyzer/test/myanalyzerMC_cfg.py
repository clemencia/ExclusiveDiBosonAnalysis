import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    )

#process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound'),
#    )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    "file:AAWW_AOD_14_1_wiI.root",
    )
)


# particle flow isolation
#
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

#
# rho value for isolation
#
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond

process.GlobalTag.globaltag = cms.string('START53_V7F::All')
process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")



#### type 1 PFMET corrections
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

#### type0 PFMET corrections needed?? makes no difference in the MET on the ntuple!!!
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False) #why false? isn't it by default?
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')
    )

process.load("MyAnalysis.MyAnalyzer.myanalyzerMC_cfi")

process.p = cms.Path(process.kt6PFJetsForIsolation * process.pfiso *
                     ## for pileup met correction ##
                     process.type0PFMEtCorrection *
                     ## for JES met correction    ##
                     process.producePFMETCorrections *
                     process.demo)
