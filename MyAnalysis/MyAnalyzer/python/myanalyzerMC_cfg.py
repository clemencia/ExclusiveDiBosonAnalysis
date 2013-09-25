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
    "/store/user/cmora/AAWW_a0wp2aCwm1_GENSIM/AAWW_a0wp2aCwm1_RECO/dd8eebe4de974f72d97bd8c652e9f8d2/AAWW_AOD_14_1_wiI.root",
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


##PFMET corrections
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")




process.GlobalTag.globaltag = 'START53_V7F::All'

process.load("MyAnalysis.MyAnalyzer.myanalyzerMC_cfi")


process.p = cms.Path(process.kt6PFJetsForIsolation * process.pfiso * process.producePFMETCorrections * process.demo)
