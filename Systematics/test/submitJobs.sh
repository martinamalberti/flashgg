#to save on eos
#mkdir /eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/94X/treesDumper_20180720/Sig
#mkdir /eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/94X/treesDumper_20180720/Data

fggRunJobs.py --load all_sig_jobs_2017.json -H -D -P -n 10 -d jobs_treesDumper_20180720_Sig/ -x cmsRun workspaceStd.py doFiducial=False tthTagsOnly=True dumpWorkspace=False dumpTrees=True doSystematics=False maxEvents=-1 --no-use-tarball -q cmscaf1nd --stage-to=/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/94X/treesDumper_20180720/Sig/ --stage-cmd="eos cp"

#fggRunJobs.py --load data_jobs_2017.json -H -D -P -n 200 -d jobs_treesDumper_20180720_Data/ -x cmsRun workspaceStd.py doFiducial=False tthTagsOnly=True dumpWorkspace=False dumpTrees=True doSystematics=False lumiMask=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt maxEvents=-1 --no-use-tarball -q 8nh --stage-to=/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/ntuples/tth/94X/treesDumper_20180720/Data/ --stage-cmd="eos cp"


