fggRunJobs.py --load jobs_SingleElectron.json -H -D -P -n 10 -d jobsDiphotonDumper_data -x cmsRun diphotonsDumperZee_cfg.py maxEvents=-1 --no-use-tarball -q 8nh
fggRunJobs.py --load jobs_DY.json -H -D -P -n 50 -d jobsDiphotonDumper_mc -x cmsRun diphotonsDumperZee_cfg.py maxEvents=-1 --no-use-tarball -q 8nh


fggRunJobs.py --load jobsZee_DoubleEG.json -H -D -P -n 10 -d jobsZeeDumper_DoubleEG -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh
fggRunJobs.py --load jobsZee_SingleElectron.json -H -D -P -n 100 -d jobsZeeDumper_SingleElectron_v2 -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh
fggRunJobs.py --load jobsZee_MC.json -H -D -P -n 100 -d jobsZeeDumper_MC_v2 -x cmsRun zeeDumper_cfg.py maxEvents=-1 --no-use-tarball -q 8nh


# to save output on eos
mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/zee_validation/GoldenJsonJuly20/jobsZeeDumper_Data_SingleElectron
fggRunJobs.py --load jobsZee_SingleElectron.json -H -D -P -n 100 -d jobsZeeDumper_Data_SingleElectron_GoldenJsonJuly20/ -x cmsRun zeeDumper_withSmearSyst_cfg.py maxEvents=-1 --no-use-tarball -q cmscaf1nd --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/zee_validation/GoldenJsonJuly20/jobsZeeDumper_Data_SingleElectron


mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/zee_validation/GoldenJsonJuly20/jobsZeeDumper_MC_DYToEE_withSmearSyst
fggRunJobs.py --load jobsZee_MC.json -H -D -P -n 200 -d jobsZeeDumper_MC_DYToEE_GoldenJsonJuly20_withSmearSyst -x cmsRun zeeDumper_withSmearSyst_cfg.py maxEvents=-1 --no-use-tarball -q cmscaf1nd --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/zee_validation/GoldenJsonJuly20/jobsZeeDumper_MC_DYToEE_withSmearSyst


mkdir ~/eos/cms/store/group/phys_higgs/cmshgg/malberti/flashgg/zee_validation/GoldenJsonJuly20/jobsZeeDumper_MC_DYToEE_withScaleSyst
fggRunJobs.py --load jobsZee_MC.json -H -D -P -n 200 -d jobsZeeDumper_MC_DYToEE_GoldenJsonJuly20_withScaleSyst -x cmsRun zeeDumper_withScaleSyst_cfg.py maxEvents=-1 --no-use-tarball -q cmscaf1nd --stage-to=/store/group/phys_higgs/cmshgg/malberti/flashgg/zee_validation/GoldenJsonJuly20/jobsZeeDumper_MC_DYToEE_withScaleSyst
