import FWCore.ParameterSet.Config as cms

TFileService = cms.Service("TFileService", fileName = cms.string("myanalysis_ntuple.root") )

isMC=True
demo = cms.EDAnalyzer('MyAnalyzer',
                      ismc=cms.bool(isMC),
                      isoValInputTags=cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                                    cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                                    cms.InputTag('elPFIsoValueNeutral03PFIdPFIso'))
)

