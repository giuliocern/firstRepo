g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA

file_names=(
["DYJetstoLL_amc_0J"]=DYJetstoLL_amc_0J
["DYJetstoLL_amc_1J"]=DYJetstoLL_amc_1J
["DYJetstoLL_amc_2J"]=DYJetstoLL_amc_2J
["SingleMuon"]=SingleMuon_reminiaod   
["ST_tW_top"]=ST_tW_top
["ST_tW_antitop"]=ST_tW_antitop
["ST_s-channel"]=ST_s
["ST_t-channel_top_4f_inclusiveDecays"]=ST_t_top
["ST_t-channel_antitop_4f_inclusiveDecays"]=ST_t_antitop
["TT"]=TT
["WW"]=WW
["WZ"]=WZ
["ZZ"]=ZZ
["WJetsToLNu"]=WJetsToLnu_madgraph
["VBF_HToMuMu"]=VBF_HToMuMu  
["GluGlu_HToMuMu"]=GluGlu_HToMuMu  
)

path=/afs/cern.ch/user/g/gimandor/private/eos/cms/store/group/phys_higgs/vbfHbb/VBFllskim/
postfix='test'
v='v25'
ROOT=.root
region=mu
ending=v25_reskim
applyJESWeight=0
applyQCDWeight=0
JESWeightNom=nom
QCDWeightNom=nom
outputdir=.

for key in ${!file_names[@]}; do
	data=0
	if [ $key == SingleMuon ]
	then
	 		data=1
			ending=v25
	fi
	f=$path${file_names[${key}]}_$ending.root
	./plot $f ${key} $region $data 0 nom 0 nom $v $postfix .
done

