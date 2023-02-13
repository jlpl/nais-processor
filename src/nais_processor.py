import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from datetime import date, datetime, timedelta
import matplotlib.dates as dts
import pandas as pd
import os
import locale
import warnings
import yaml
import re
import sys
from dateutil.parser import parse
from tinydb import TinyDB, Query
from tinydb.operations import add
import time
import json
import aerosol_functions as af

__pdoc__ = {
    'tubeloss': False,
    'average_mob': False,
    'average_dp': False,
    'find_diagnostic_names': False,
    'process_data': False,
    'get_environmental_data': False,
    'choose_particle_polarity': False,
    'bring_to_sealevel': False,
    'correct_inlet_losses': False,
    'wagner_ion_mode_correction': False,
    'clean_elem_noise': False,
    'clean_corona_ions': False,
    'add_flags': False,
}


# The final geometric mean diameters of diameter and mobility bins
dp_ion = np.array([7.86360416e-10, 9.08232168e-10, 1.04902018e-09, 1.21167006e-09,
1.39958930e-09, 1.61672083e-09, 1.86762862e-09, 2.15759741e-09,
2.49274932e-09, 2.88018000e-09, 3.32811839e-09, 3.84611427e-09,
4.44525917e-09, 5.13844742e-09, 5.94068566e-09, 6.86946146e-09,
7.94518431e-09, 9.19171623e-09, 1.06370142e-08, 1.23139134e-08,
1.42610904e-08, 1.65242568e-08, 1.91576555e-08, 2.22259544e-08,
2.58066722e-08, 2.99933244e-08, 3.48995548e-08, 4.06646353e-08])*1e9

dp_par = np.array([7.498942093324539870e-01,8.659643233600640144e-01,
9.999999999999980016e-01,1.154781984689456031e+00,1.333521432163321974e+00,
1.539926526059490097e+00,1.778279410038920094e+00,2.053525026457140079e+00,
2.371373705661659947e+00,2.738419634264360081e+00,3.162277660168379967e+00,
3.651741272548380213e+00,4.216965034285819591e+00,4.869675251658620141e+00,
5.623413251903479626e+00,6.493816315762099833e+00,7.498942093324560076e+00,
8.659643233600640144e+00,1.000000000000000000e+01,1.154781984689457985e+01,
1.333521432163323972e+01,1.539926526059490008e+01,1.778279410038922137e+01,
2.053525026457139901e+01,2.371373705661660125e+01,2.738419634264360170e+01,
3.162277660168379967e+01,3.651741272548380124e+01,4.216965034285819769e+01])

mob_ion = np.array([3.162277660168379937e-04,2.371373705661659990e-04,
1.778279410038920258e-04,1.333521432163320159e-04,1.000000000000000048e-04,
7.498942093324559917e-05,5.623413251903490022e-05,4.216965034285820205e-05,
3.162277660168380208e-05,2.371373705661660125e-05,1.778279410038919852e-05,
1.333521432163319990e-05,1.000000000000000082e-05,7.498942093324561442e-06,
5.623413251903490361e-06,4.216965034285830030e-06,3.162277660168380038e-06,
2.371373705661659871e-06,1.778279410038920148e-06,1.333521432163330027e-06,
1.000000000000000167e-06,7.498942093324570124e-07,5.623413251903499890e-07,
4.216965034285829924e-07,3.162277660168379721e-07,2.371373705661660136e-07,
1.778279410038920042e-07,1.333521432163329868e-07])*1e4

mob_ion_geomeans=np.array([2.73841963e-04, 2.05352503e-04, 1.53992653e-04, 1.15478198e-04,
8.65964323e-05, 6.49381632e-05, 4.86967525e-05, 3.65174127e-05,
2.73841963e-05, 2.05352503e-05, 1.53992653e-05, 1.15478198e-05,
8.65964323e-06, 6.49381632e-06, 4.86967525e-06, 3.65174127e-06,
2.73841963e-06, 2.05352503e-06, 1.53992653e-06, 1.15478198e-06,
8.65964323e-07, 6.49381632e-07, 4.86967525e-07, 3.65174127e-07,
2.73841963e-07, 2.05352503e-07, 1.53992653e-07])*1e4

dp_par_geomeans=np.array([0.80584219,  0.93057204,  1.07460783,  1.24093776,  1.43301257,
1.6548171 ,  1.91095297,  2.20673407,  2.54829675,  2.94272718,
3.39820833,  3.92418976,  4.53158364,  5.23299115,  6.0429639 ,
6.97830585,  8.05842188,  9.30572041, 10.74607828, 12.40937761,
14.3301257 , 16.548171  , 19.10952975, 22.06734069, 25.48296748,
29.42727176, 33.98208329, 39.24189758])

dlogmob_ion=np.array([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
0.125])

dlogdp_ion = np.array([0.06257524, 0.0625811 , 0.06259375, 0.06260838, 0.06262533,
0.06264495, 0.06266769, 0.06269404, 0.06272461, 0.06276008,
0.06280128, 0.06284916, 0.06290487, 0.06296974, 0.06304539,
0.0631337 , 0.06323696, 0.06335788, 0.06349974, 0.0636665 ,
0.06386292, 0.06409481, 0.06436924, 0.06469482, 0.06508209,
0.06554394, 0.06609614, 0.06639699])

dlogdp_par=np.array([0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
0.0625, 0.0625, 0.0625, 0.0625, 0.0625])

filename_formats = [
["%Y-%m-%d.ions.nds","%Y-%m-%d.particles.nds","%Y-%m-%d.log"],
["%Y%m%d-block-ions.spectra","%Y%m%d-block-particles.spectra","%Y%m%d-block.records"],
["%Y%m%d-block-ions.spectra","%Y%m%d-block-particles.spectra","%Y%m%d-block.diagnostics"]]

possible_sampleflow_names1 = [
"sampleflow",
"Flowaer"]

possible_sampleflow_names2 = [
"pos_sampleflow.mean",
"neg_sampleflow.mean",
"pos_sampleflow",
"neg_sampleflow"
]

possible_temperature_names = [
"temperature.mean",
"temperature",
"temp"]

possible_pressure_names = [
"baro.mean",
"baro"]

# electrometer size ranges for different inverters (for the purpose of cleaning out electrometer noise):
ions_pos_v14_lrnd={"0": [7.16444775804687e-10, 1.0700473216535486e-09], "1": [8.766005865635541e-10, 1.2912139078236106e-09], "2":
 [1.0233784015731513e-09, 1.494607599390042e-09], "3": [1.167004143869059e-09, 1.6953050539978397e-09], "4":
 [1.3171140158277396e-09, 1.9129953633709412e-09], "5": [1.5010400712091295e-09, 2.196448726880819e-09], "6":
 [1.7374102549917467e-09, 2.5397495919423145e-09], "7": [1.9987455846433743e-09, 2.909827835540643e-09], "8":
 [2.308391221999045e-09, 3.399628318414748e-09], "9": [2.743654460662328e-09, 4.064377206025429e-09], "10":
 [3.235105579106799e-09, 4.708980203333236e-09], "11": [3.657786198230896e-09, 5.320149444473465e-09], "12":
 [4.40743242351629e-09, 7.235758783698336e-09], "13": [6.341141170615947e-09, 1.0173608443214825e-08], "14":
 [8.61139257420043e-09, 1.2899482426374689e-08], "15": [1.0474248763637253e-08, 1.556121847426194e-08], "16":
 [1.2937549036927316e-08, 2.07105715562843e-08], "17": [1.778277482687919e-08, 2.8143993934982054e-08], "18":
 [2.3703396062208645e-08, 3.638374245002531e-08], "19": [2.9466514534877575e-08, 4.3834357529742775e-08], "20":
 [3.3648711881521194e-08, 4.601489588650497e-08]}

particles_pos_v14_lrnd_elm25_chv={"0": [7.084021223198885e-10, 1.1041550022447838e-09], "1": [8.824680153980799e-10, 1.3381985992962885e-09], "2":
 [1.027268571461551e-09, 1.5535807431315568e-09], "3": [1.1710814400196452e-09, 1.7778453394121642e-09], "4":
 [1.3364125967253046e-09, 2.0521125317491067e-09], "5": [1.5165756393379476e-09, 2.364104261278049e-09], "6":
 [1.7469782917199996e-09, 2.7691235104649e-09], "7": [2.0030670784526445e-09, 3.1804078983335615e-09], "8":
 [2.3555174234652947e-09, 3.889745378517969e-09], "9": [2.8998890919509277e-09, 4.85643554134041e-09], "10":
 [3.4535426802312536e-09, 5.802269084443331e-09], "11": [3.952828282143858e-09, 6.706155173416103e-09], "12":
 [4.680496931685137e-09, 8.345157694215357e-09], "13": [5.908440036792274e-09, 1.0914143399425853e-08], "14":
 [7.25270144756779e-09, 1.3424351306635164e-08], "15": [8.481012704476841e-09, 1.5854436158750896e-08], "16":
 [9.92828580185976e-09, 1.8875874031548135e-08], "17": [1.2480795668183397e-08, 2.5489807672756865e-08], "18":
 [1.618264753193867e-08, 3.382694961210274e-08], "19": [1.9259686684447952e-08, 4.1933601744775866e-08], "20":
 [2.22037089383357e-08, 4.865271942977691e-08], "21": [2.4978007194778478e-08, 5.377912883862792e-08], "22":
 [2.7446273522212072e-08, 5.731891713329477e-08], "23": [2.966486504673139e-08, 6.006253448103395e-08], "24": 
 [3.1538173585926114e-08, 6.213974119411533e-08]}

particles_neg_v14_lrnd={"0": [7.071203017729596e-10, 1.1005013232430319e-09], "1": [8.797742415579472e-10, 1.3341534055269684e-09], "2":
 [1.0249360413125393e-09, 1.5500231709924995e-09], "3": [1.1663524905165027e-09, 1.76746761608056e-09], "4":
 [1.3247911530621861e-09, 2.0313528986656205e-09], "5": [1.4993872332769961e-09, 2.3322201548214762e-09], "6":
 [1.7211922995795243e-09, 2.7237415257553574e-09], "7": [1.9709779821521803e-09, 3.127660894213944e-09], "8":
 [2.2719521313715783e-09, 3.683906828365737e-09], "9": [2.6910604870370505e-09, 4.432425549741674e-09], "10":
 [3.1273717566860052e-09, 5.176267101656894e-09], "11": [3.5214138488359515e-09, 5.88784346872649e-09], "12":
 [4.344519392306214e-09, 8.075601998567308e-09], "13": [6.055647568060633e-09, 1.15707771780943e-08], "14":
 [7.891354399664643e-09, 1.4822317894774355e-08], "15": [9.488534907698969e-09, 1.7987926934224552e-08], "16": 
 [1.1972289469619363e-08, 2.4352144139698158e-08], "17": [1.6205665582181165e-08, 3.529177218334523e-08], "18": 
 [2.1726991804436577e-08, 4.890078120743212e-08], "19": [2.747153199034406e-08, 5.80468336800726e-08], "20": 
 [3.181625187151363e-08, 6.282803902888346e-08]}

ions_pos_v141_lrnd_elm25_chv={"0": [7.174408260504354e-10, 1.073624396936229e-09], "1": [8.795427161035404e-10, 1.2949773919916534e-09], "2":
 [1.0258582665799353e-09, 1.4979919549108165e-09], "3": [1.1719327677994699e-09, 1.70372171508208e-09], "4":
 [1.32799146465657e-09, 1.932594430907855e-09], "5": [1.520240256748169e-09, 2.2261429400706174e-09], "6":
 [1.7628304747501982e-09, 2.581861661002808e-09], "7": [2.0323715276272585e-09, 2.9598100485518586e-09], "8":
 [2.3816068293602794e-09, 3.5773872919384414e-09], "9": [2.962208652767496e-09, 4.432897017241381e-09], "10":
 [3.5762132418655334e-09, 5.239831433099603e-09], "11": [4.133555700629352e-09, 6.039237962304223e-09], "12": 
 [4.882640365333869e-09, 7.457408037033863e-09], "13": [6.2783199783114545e-09, 9.575956296394927e-09], "14": 
 [7.904204442058769e-09, 1.1694806122760362e-08], "15": [9.352755529575812e-09, 1.3732975020325302e-08], "16": 
 [1.0981897607602298e-08, 1.6276182350230556e-08], "17": [1.3500089105006028e-08, 2.1603165023592403e-08], "18": 
 [1.8022095523973226e-08, 2.708852306358357e-08], "19": [2.156970936611481e-08, 3.1921441726778743e-08], "20": 
 [2.4519917336200788e-08, 3.609220688998671e-08], "21": [2.7168015322179712e-08, 3.993591879861599e-08], "22": 
 [2.9639336081649807e-08, 4.339405159440388e-08], "23": [3.19559928473482e-08, 4.538489884003738e-08], "24":
 [3.362603005585794e-08, 4.572616914940904e-08]}

ions_pos_v141_hrnd_elm25_chv={"0": [7.177816842294488e-10, 1.072221766346416e-09], "1": [8.784324649240066e-10, 1.2952331371673022e-09], "2":
 [1.0282272593935292e-09, 1.496448130128779e-09], "3": [1.1709728034437686e-09, 1.70109575761716e-09], "4":
 [1.3300366642732846e-09, 1.9341044305153346e-09], "5": [1.521422514332051e-09, 2.220832669677381e-09], "6":
 [1.7632530836302399e-09, 2.584244551737274e-09], "7": [2.0351677980258766e-09, 2.955213035032816e-09], "8": 
 [2.379947213073259e-09, 3.5807287457548474e-09], "9": [2.9614387929489157e-09, 4.428794241338281e-09], "10": 
 [3.5808241708120252e-09, 5.247041226910477e-09], "11": [4.136034551200579e-09, 6.026277543624166e-09], "12": 
 [4.890752691340854e-09, 7.457985468634934e-09], "13": [6.2730994769101586e-09, 9.589732077862633e-09], "14":
 [7.893847904275436e-09, 1.1687225512360978e-08], "15": [9.345340174460752e-09, 1.3762043152233714e-08], "16":
 [1.0978493586549243e-08, 1.6262760201045273e-08], "17": [1.3526934344462764e-08, 2.157926406305932e-08], "18":
 [1.8045823361694697e-08, 2.711029182088656e-08], "19": [2.1550025677740064e-08, 3.1861088201059746e-08], "20": 
 [2.453350988841541e-08, 3.611626721404304e-08], "21": [2.7196984694271958e-08, 3.995278834492865e-08], "22": 
 [2.9641573529840143e-08, 4.3318065851782826e-08], "23": [3.187419245775203e-08, 4.5478238497900765e-08], "24":
 [3.3521538858575206e-08, 4.6046488264314935e-08]}

particles_pos_v14_lrnd={"0": [7.071203017729596e-10, 1.1005013232430319e-09], "1": [8.797742415579472e-10, 1.3341534055269684e-09], "2":
 [1.0249360413125393e-09, 1.5500231709924995e-09], "3": [1.1663524905165027e-09, 1.76746761608056e-09], "4":
 [1.3247911530621861e-09, 2.0313528986656205e-09], "5": [1.4993872332769961e-09, 2.3322201548214762e-09], "6":
 [1.7211922995795243e-09, 2.7237415257553574e-09], "7": [1.9709779821521803e-09, 3.127660894213944e-09], "8": 
 [2.2719521313715783e-09, 3.683906828365737e-09], "9": [2.6910604870370505e-09, 4.432425549741674e-09], "10": 
 [3.1273717566860052e-09, 5.176267101656894e-09], "11": [3.5214138488359515e-09, 5.88784346872649e-09], "12": 
 [4.344519392306214e-09, 8.075601998567308e-09], "13": [6.055647568060633e-09, 1.15707771780943e-08], "14": 
 [7.891354399664643e-09, 1.4822317894774355e-08], "15": [9.488534907698969e-09, 1.7987926934224552e-08], "16": 
 [1.1972289469619363e-08, 2.4352144139698158e-08], "17": [1.6205665582181165e-08, 3.529177218334523e-08], "18": 
 [2.1726991804436577e-08, 4.890078120743212e-08], "19": [2.747153199034406e-08, 5.80468336800726e-08], "20": 
 [3.181625187151363e-08, 6.282803902888346e-08]}

ions_neg_v141_lrnd_elm25_chv={"0": [7.162722829189713e-10, 1.069396256829605e-09], "1": [8.777038868584165e-10, 1.2946446240931807e-09], "2":
 [1.02759999667648e-09, 1.5009614340173589e-09], "3": [1.1745671577354954e-09, 1.7069296262514892e-09], "4":
 [1.330366866209285e-09, 1.936923662090671e-09], "5": [1.523569941517978e-09, 2.230584153641438e-09], "6":
 [1.763456724402997e-09, 2.5797080071606806e-09], "7": [2.0282822638164162e-09, 2.9544966841050904e-09], "8":
 [2.3814411392836232e-09, 3.585069910708154e-09], "9": [2.975412233075357e-09, 4.4586459608966074e-09], "10":
 [3.5950536764791088e-09, 5.264604651129057e-09], "11": [4.1451727671809455e-09, 6.0533707850852946e-09], "12":
 [4.908070335058664e-09, 7.554663305133556e-09], "13": [6.370251631603248e-09, 9.708647141925531e-09], "14":
 [7.99768693004782e-09, 1.1807392092933995e-08], "15": [9.421664599359865e-09, 1.3820815597006762e-08], "16":
 [1.1077372179456622e-08, 1.6484128380953205e-08], "17": [1.3672218205379754e-08, 2.166074698335169e-08], "18":
 [1.8004147704745442e-08, 2.705025422013022e-08], "19": [2.1526595943708422e-08, 3.186462197910958e-08], "20":
 [2.4479448411058165e-08, 3.602967602639687e-08], "21": [2.7167075804775312e-08, 3.995596608326538e-08], "22":
 [2.9709429581149162e-08, 4.347660778516122e-08], "23": [3.2017621669510734e-08, 4.541693553734151e-08], "24":
 [3.3681779306509456e-08, 4.5708942542429593e-08]}

ions_neg_v141_hrnd_elm25_chv={"0": [7.166088691273807e-10, 1.0680529042034586e-09], "1": [8.766077238618173e-10, 1.2948822433267513e-09], "2":
 [1.0299408040688298e-09, 1.499348034087609e-09], "3": [1.1734785792847775e-09, 1.7044494278972901e-09], "4":
 [1.3325133850369467e-09, 1.9382490836139435e-09], "5": [1.5246381007528337e-09, 2.225245739179207e-09], "6":
 [1.7639085256336676e-09, 2.5820656631195802e-09], "7": [2.031125366865139e-09, 2.9500866762907164e-09], "8":
 [2.3797421811002037e-09, 3.588570659248709e-09], "9": [2.974497863791539e-09, 4.454598568481367e-09], "10":
 [3.600685696291343e-09, 5.270248740213011e-09], "11": [4.147441219601736e-09, 6.040242224249017e-09], "12":
 [4.91677135667135e-09, 7.553001655458096e-09], "13": [6.368408114819548e-09, 9.720375884730652e-09], "14":
 [7.983853657878074e-09, 1.1805425434050108e-08], "15": [9.415625209974426e-09, 1.385125858049306e-08], "16":
 [1.1073898998765917e-08, 1.646956079843046e-08], "17": [1.3699649532490082e-08, 2.163235038140565e-08], "18":
 [1.802716724985027e-08, 2.7073230406563897e-08], "19": [2.150965694231408e-08, 3.18014983817787e-08], "20":
 [2.4492014342932346e-08, 3.605523207826252e-08], "21": [2.7196391232036758e-08, 3.9971877339499414e-08], "22":
 [2.9708672750890287e-08, 4.340340668907393e-08], "23": [3.1934336215501455e-08, 4.5516606356023787e-08], "24":
 [3.357606684507602e-08, 4.604620001912706e-08]}

particles_neg_v14_hrnd_elm25_chv={"0": [7.064042348613744e-10, 1.098448642629657e-09], "1": [8.804166704055082e-10, 1.3382314628239056e-09], "2":
 [1.029929687471369e-09, 1.5587130355258417e-09], "3": [1.1715316139497021e-09, 1.7855361801274118e-09], "4":
 [1.3375770105761571e-09, 2.0569498417960264e-09], "5": [1.5213418055613532e-09, 2.3650454557473994e-09], "6":
 [1.7508747954828767e-09, 2.7596483400473384e-09], "7": [2.0016660559184965e-09, 3.1683106443504784e-09], "8":
 [2.3537662142103147e-09, 3.892364103276864e-09], "9": [2.908132768324034e-09, 4.87724949455506e-09], "10":
 [3.4647588825814346e-09, 5.819021297995083e-09], "11": [3.9538012917605546e-09, 6.718817658651353e-09], "12":
 [4.699417872175317e-09, 8.46741114160911e-09], "13": [5.980472218281243e-09, 1.1075320194303405e-08], "14":
 [7.32788314110871e-09, 1.356178222960962e-08], "15": [8.546658915188655e-09, 1.5975204027458133e-08], "16":
 [1.003513710889909e-08, 1.9121314649296197e-08], "17": [1.2584497621235887e-08, 2.5586370462254533e-08], "18":
 [1.6193366252358348e-08, 3.375225311508183e-08], "19": [1.9296366093983868e-08, 4.135787670654157e-08], "20":
 [2.2205359456932222e-08, 4.775609758877974e-08], "21": [2.4949980719630763e-08, 5.296183331803598e-08], "22":
 [2.74591102013903e-08, 5.6781631708117125e-08], "23": [2.9663083522354196e-08, 5.960548252960392e-08], "24":
 [3.1563036449021557e-08, 6.184741780170188e-08]}

particles_pos_v14_hrnd_elm25_chv={"0": [7.078936573300082e-10, 1.1028124202267342e-09], "1": [8.820689074711288e-10, 1.338578480376556e-09], "2":
 [1.0282838004812897e-09, 1.5555733231186909e-09], "3": [1.1690201376888472e-09, 1.78149953583638e-09], "4":
 [1.3348577366613476e-09, 2.052450596702526e-09], "5": [1.5181851404773485e-09, 2.3598931859583523e-09], "6":
 [1.7509905816449534e-09, 2.7622818843531693e-09], "7": [2.005225058134512e-09, 3.173616142257642e-09], "8":
 [2.3525975248633675e-09, 3.883241010415373e-09], "9": [2.894779176315381e-09, 4.847916178622021e-09], "10":
 [3.4477052157384657e-09, 5.792051680298515e-09], "11": [3.944104923721565e-09, 6.701786044525521e-09], "12":
 [4.666951603466218e-09, 8.360587273279385e-09], "13": [5.8997089753584225e-09, 1.0919012148501147e-08], "14":
 [7.252706673752877e-09, 1.342019953209076e-08], "15": [8.486965748259635e-09, 1.586474588608095e-08], "16":
 [9.936765651010928e-09, 1.8874421710969515e-08], "17": [1.246739709406105e-08, 2.5483091135053706e-08], "18":
 [1.6208024471095697e-08, 3.3797206216952925e-08], "19": [1.9338854709118882e-08, 4.146124436417098e-08], "20":
 [2.2243072357824654e-08, 4.783315990811425e-08], "21": [2.4952099445624026e-08, 5.296484079233324e-08], "22":
 [2.7389130112705036e-08, 5.667341170880334e-08], "23": [2.9602350590226694e-08, 5.953977095698898e-08], "24":
 [3.149305285177872e-08, 6.176077534814207e-08]}

ions_neg_v14_lrnd={"0": [7.16444775804687e-10, 1.0700473216535486e-09], "1": [8.766005865635541e-10, 1.2912139078236106e-09], "2":
 [1.0233784015731513e-09, 1.494607599390042e-09], "3": [1.167004143869059e-09, 1.6953050539978397e-09], "4":
 [1.3171140158277396e-09, 1.9129953633709412e-09], "5": [1.5010400712091295e-09, 2.196448726880819e-09], "6":
 [1.7374102549917467e-09, 2.5397495919423145e-09], "7": [1.9987455846433743e-09, 2.909827835540643e-09], "8":
 [2.308391221999045e-09, 3.399628318414748e-09], "9": [2.743654460662328e-09, 4.064377206025429e-09], "10":
 [3.235105579106799e-09, 4.708980203333236e-09], "11": [3.657786198230896e-09, 5.320149444473465e-09], "12":
 [4.40743242351629e-09, 7.235758783698336e-09], "13": [6.341141170615947e-09, 1.0173608443214825e-08], "14":
 [8.61139257420043e-09, 1.2899482426374689e-08], "15": [1.0474248763637253e-08, 1.556121847426194e-08], "16":
 [1.2937549036927316e-08, 2.07105715562843e-08], "17": [1.778277482687919e-08, 2.8143993934982054e-08], "18":
 [2.3703396062208645e-08, 3.638374245002531e-08], "19": [2.9466514534877575e-08, 4.3834357529742775e-08], "20":
 [3.3648711881521194e-08, 4.601489588650497e-08]}

particles_neg_v14_lrnd_elm25_chv={"0": [7.069362615274498e-10, 1.0997213171324974e-09], "1": [8.808566990541904e-10, 1.3378548681525188e-09], "2":
 [1.0289107593642929e-09, 1.5567373372470626e-09], "3": [1.1735882698862885e-09, 1.7819165288950059e-09], "4":
 [1.33907446265167e-09, 2.0566964315737684e-09], "5": [1.5196769831618589e-09, 2.3693061857773854e-09], "6":
 [1.746857574391727e-09, 2.766483303572841e-09], "7": [1.9994652400593872e-09, 3.175124915389092e-09], "8":
 [2.356723607259348e-09, 3.898924433703989e-09], "9": [2.9132704671892387e-09, 4.885896237782877e-09], "10":
 [3.4706415007361796e-09, 5.829183030947836e-09], "11": [3.962593522847839e-09, 6.722985797091222e-09], "12":
 [4.7129989830246345e-09, 8.451324582668515e-09], "13": [5.988747492408655e-09, 1.1071809141751731e-08], "14":
 [7.327476301230385e-09, 1.3565350436407152e-08], "15": [8.540412648877852e-09, 1.5964663591630477e-08], "16":
 [1.0026869681645872e-08, 1.9124820128830518e-08], "17": [1.2598231836692178e-08, 2.5591353853168848e-08], "18":
 [1.6168464020950793e-08, 3.3780310510756356e-08], "19": [1.9217207258451423e-08, 4.1823547718041306e-08], "20":
 [2.2165154711387387e-08, 4.85721513696297e-08], "21": [2.4975870919608738e-08, 5.3775910543120685e-08], "22":
 [2.7516844715923574e-08, 5.7423919806126056e-08], "23": [2.9725447211466474e-08, 6.012495231939605e-08], "24":
 [3.160714147308623e-08, 6.222035590925216e-08]}

# Define standard conditions
temp_ref = 273.15 # K, 0C
pres_ref = 101325.0 # Pa, 1atm

def make_config_template(fn):
    """
    Make a configuration file template

    Parameters
    ----------

    fn : str
        full path to configuration file

        For example `/home/user/config.yml`

    Notes
    -----

    The default values are used to calculate the corrections
    in case the data is not available in the diagnostic data 
    either due missing or broken sensor.

    """
    with open(fn,"w") as f:
        f.write("measurement_location: # Name of the measurement site\n")
        f.write("data_folder: # Full paths to raw data folders\n")
        f.write("- # Data folder 1\n")
        f.write("- # Data folder 2, and so on...\n")
        f.write("processed_folder: # Full path to folder where procesed data is saved\n")
        f.write("database_file: # Full path to database file (will be created on first run) \n")
        f.write("start_date: # Format: yyyy-mm-dd\n")
        f.write("end_date: # Format: yyyy-mm-dd or '' for current day\n")
        f.write("inlet_length: # length of inlet in meters\n")
        f.write("do_inlet_loss_correction: # true or false\n")
        f.write("convert_to_standard_conditions: # true or false\n")
        f.write("do_wagner_ion_mode_correction: # true or false\n")
        f.write("remove_corona_ions: # true or false\n")
        f.write("remove_noisy_electrometers: # true or false\n")
        f.write("inverter_name: # hires_25, lores_25, lores_21 or '' (needed for noise removal, '' if noise not removed)\n")
        f.write("allow_reprocess: # true or false")
        f.write("choose_better_particle_polarity: # true or false")
        f.write("use_default_values: # true or false")
        f.write("default_temperature: # temperature in K used in corrections as fallback")
        f.write("default_pressure: # pressure in Pa used in corrections as fallback")
        f.write("default_flowrate: # flow rate in lpm used in corrections as fallback")
        f.write("include_flags: # include flags to the data file and make a separete flag file for each day, true or false")

def tubeloss(dpp,pflow,plength,temp,press):
    DPP,TEMP = np.meshgrid(dpp,temp)
    DPP,PRESS = np.meshgrid(dpp,press)
    DPP,PFLOW = np.meshgrid(dpp,pflow)
    rmuu = np.pi*af.particle_diffusivity(DPP,TEMP,PRESS)*plength/PFLOW
    pene = np.nan*np.ones(rmuu.shape)
    cond1=rmuu<0.02
    cond2=rmuu>=0.02
    pene[cond1] = 1. - 2.56*rmuu[cond1]**(2./3.) + 1.2*rmuu[cond1]+0.177*rmuu[cond1]**(4./3.)
    pene[cond2] = 0.819*np.exp(-3.657*rmuu[cond2]) + 0.097*np.exp(-22.3*rmuu[cond2]) + 0.032*np.exp(-57.0*rmuu[cond2])
    return pene

def read_file(fn,ftype):
    """
    Read NAIS raw data file into a pandas.DataFrame

    Parameters
    ----------

    fn : str
        Raw data filename with path

    ftype : str
        `"spectra"` (inverted size/mobility distribution) or
        `"records"` (diagnostic data and electrometer currents)

    Returns
    -------

    pandas.DataFrame
        Contents of the file

    str
        Explantions of flags, returned only if `ftype="records"`

    """

    with open(fn,'r') as f:

        header_found = False
        data_matrix = []
        flag_explanations = []
        lines = f.read().splitlines()
        
        for line in lines:

             # Skip empty and comments
             if (len(line)==0):
                 continue

             # Collect a list of flags and skip comments
             if line[:6]=="# flag":
                 # parse the line
                 diagnostic_comment_yaml = yaml.safe_load(line[7:].rstrip('\r\n'))
                 flag_name = list(diagnostic_comment_yaml.keys())[0]
                 flag_message = diagnostic_comment_yaml[flag_name]["message"]
                 flag_explanations.append([flag_name,flag_message])
             elif (line[0]=='#'):
                 continue
             else:
                 pass

             # Test if it is a header
             if (header_found==False):
                 if "opmode" in line:
                     delimiter = re.search('(.)opmode',line).group(1)
                     header = line.split(delimiter)
                     number_of_columns = len(header)
                     header_found = True
                     continue
                 else:
                     continue 
             else:
                 data_line = line.split(delimiter)
                 
                 if ((len(data_line)==number_of_columns) & ("opmode" not in data_line)):
                     data_matrix.append(data_line)
                 continue

    if len(data_matrix)==0:
        return None

    else:
        # Convert anything that can be converted to float and the rest is coerced to NaNs
        df = pd.DataFrame(columns = header, data = data_matrix)
        df_flags = pd.DataFrame(columns=["Flag","Message"], data = flag_explanations)

        # records: start_time, end_time, opmode, data..., flags
        # spectra: start_time, end_time, opmode, data...
        if ftype=="records":
            df.iloc[:,3:-1] = df.iloc[:,3:-1].apply(pd.to_numeric, errors='coerce').astype(float)
        if ftype=="spectra":
            df.iloc[:,3:] = df.iloc[:,3:].apply(pd.to_numeric, errors='coerce').astype(float)

        # Establish begin_time (first column) as index
        df = df.set_index(df.columns[0])
        df.index = pd.to_datetime(df.index)

        # if there is no tz information set the timezone to UTC
        df.index = [t.tz_localize('UTC') if (t.tzinfo is None) else t for t in df.index]
        
        if ftype=="records":
            return df, df_flags
        if ftype=="spectra":
            return df

def average_mob(y,h):
    data = pd.DataFrame([])
    
    for i in range(0,len(mob_ion_geomeans)):
        if i==0:
            y_block = y.iloc[:,h>mob_ion_geomeans[i]]
        else:
            y_block = y.iloc[:,((h>mob_ion_geomeans[i]) & (h<=mob_ion_geomeans[i-1]))]

        data[i] = y_block.median(axis=1)

    y_block = y.iloc[:,h<=mob_ion_geomeans[i]]
    data[i+1] = y_block.mean(axis=1)

    return data

def average_dp(y,h):
    data = pd.DataFrame([])

    for i in range(0,len(dp_par_geomeans)):
        if i==0:
            y_block = y.iloc[:,h<dp_par_geomeans[i]]
        else:
            y_block = y.iloc[:,((h<dp_par_geomeans[i]) & (h>=dp_par_geomeans[i-1]))]

        data[i] = y_block.median(axis=1)

    y_block = y.iloc[:,h>=dp_par_geomeans[i]]
    data[i+1] = y_block.mean(axis=1)

    return data

def find_diagnostic_names(diag_params):

    sampleflow_name=None
    temperature_name=None
    pressure_name=None

    for temp_name in possible_temperature_names:
         if temp_name in diag_params:
             temperature_name = temp_name
             break

    for pres_name in possible_pressure_names:
        if pres_name in diag_params:
            pressure_name = pres_name
            break

    # try single flow sensor
    for flow_name in possible_sampleflow_names1:
        if flow_name in diag_params:
            sampleflow_name = flow_name
            break

    if sampleflow_name is None:
        # try two flow sensors
        sf_name = []
        for flow_name in possible_sampleflow_names2:
            if flow_name in diag_params:
                sf_name.append(flow_name)
        if len(sf_name)==2:
            sampleflow_name=sf_name

    return temperature_name, pressure_name, sampleflow_name

def process_data(df,mode):
    """
    Convert spectra files to .sum format
    """

    if (df is None):
        return None, None

    elif not df.index.to_series().is_monotonic_increasing:
        return None, None

    else:
        df_columns = df.columns
        df_inverter_reso = int((len(df_columns)-2)/4)

        neg_df = df.iloc[:,2:2+df_inverter_reso]
        pos_df = df.iloc[:,2+2*df_inverter_reso:2+3*df_inverter_reso]

        if mode=="ions":
            mob_ion_inv = np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",y)[0])
                                    for y in df_columns[2:2+df_inverter_reso]])

            neg_df = average_mob(neg_df,mob_ion_inv)
            pos_df = average_mob(pos_df,mob_ion_inv)

            # Convert to number size distributions
            neg_df = neg_df * dlogmob_ion / dlogdp_ion
            pos_df = pos_df * dlogmob_ion / dlogdp_ion

        if mode=="particles":
            dp_par_inv = 2.0*np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",y)[0])
                                       for y in df_columns[2:2+df_inverter_reso]])
        
            neg_df = average_dp(neg_df,dp_par_inv)
            pos_df = average_dp(pos_df,dp_par_inv)

        # Construct the headers
        if mode=="ions":
            df_header = dp_ion*1e-9
        if mode=="particles":
            df_header = dp_par*1e-9

        negdf = pd.DataFrame(columns=df_header, index=df.index, data=neg_df.values)
        posdf = pd.DataFrame(columns=df_header, index=df.index, data=pos_df.values)

        negdf.index.name = "Time"
        posdf.index.name= "Time"

        if negdf.isna().all().all():
            negdf = None
        if posdf.isna().all().all():
            posdf = None

        return negdf, posdf

def get_environmental_data(
    df,
    rec,
    mode,
    use_default_values,
    default_pressure,
    default_temperature,
    default_flowrate):

    if ((rec is None) or (df is None)):
        return None,None,None

    else:        
        # Extract the records that match the mode
        if mode=="ions":
            df_rec = rec[rec.opmode=='ions']
        if mode=="particles":
            df_rec = rec[rec.opmode=='particles']

        if not df_rec.index.to_series().is_monotonic_increasing:
            return None,None,None
        
        df_rec = df_rec.reindex(df.index,method="nearest")

        # Check that the relevant diagnostic data is found
        t_name,p_name,sf_name = find_diagnostic_names(list(df_rec))

        if t_name is not None:
            t_df = 273.15 + pd.DataFrame(df_rec[t_name].astype(float))
            # Values may be missing: e.g. sensor is broken
            if (t_df.isna().all().all() and use_default_values):
                t_df = pd.DataFrame(index = df.index)
                t_df[0] = default_temperature
        elif use_default_values:
            t_df = pd.DataFrame(index = df.index)
            t_df[0] = default_temperature
        else:
            t_df = None

        if p_name is not None:
            p_df = 100.0 * pd.DataFrame(df_rec[p_name].astype(float))
            if (p_df.isna().all().all() and use_default_values):
                p_df = pd.DataFrame(index = df.index)
                p_df[0] = default_pressure
        elif use_default_values:
            p_df = pd.DataFrame(index = df.index)
            p_df[0] = default_pressure
        else:
            p_df = None

        if sf_name is not None:
            if len(sf_name)==2:
                flow_df = pd.DataFrame(df_rec[sf_name].sum(axis=1,min_count=2).astype(float))
            else:
                flow_df = pd.DataFrame(df_rec[sf_name].astype(float))
            # Test if the sampleflow is in cm3/s (old models) or
            # l/min and if necessary convert to l/min
            if (np.nanmedian(flow_df)>300):
                flow_df = (flow_df/1000.0) * 60.0
            else:
                pass
            if (flow_df.isna().all().all() and use_default_values):
                flow_df = pd.DataFrame(index = df.index)
                flow_df[0] = default_flowrate
        elif use_default_values:
            flow_df = pd.DataFrame(index = df.index)
            flow_df[0] = default_flowrate
        else:
            flow_df = None

        # Sanity check the values
        if t_df is not None:
            t_df = t_df.where(((t_df>=223.)|(t_df<=353.)),np.nan)
        
        if p_df is not None:
            p_df = p_df.where(((p_df>=37000.)|(p_df<=121000.)),np.nan)
        
        if flow_df is not None:
            flow_df = flow_df.where(((flow_df>=48.)|(flow_df<=65.)),np.nan)
     
        return t_df, p_df, flow_df

def bring_to_sealevel(
    df,
    t_df,
    p_df):
    """
    Notes
    -----

    NAIS keeps constant volumetric flowrate, However air expands
    and compresses depending on the pressure and temperature, chaging
    the number of particles per unit volume. In order to compare
    concentrations we need to transform the concentrations to standard
    conditions.
    """

    if ((df is None) or (t_df is None) or (p_df is None)):
        return None
    else:
        stp_corr_df = (pres_ref*t_df.values)/(temp_ref*p_df.values)
        df = stp_corr_df * df
        
        return df

def correct_inlet_losses(
    df,
    mode,
    pipe_length,
    t_df,
    p_df,
    flow_df):

    if ((df is None) or (t_df is None) or (p_df is None) or (flow_df is None)):
        return None

    # Diffusion loss correction
    if mode=="ions":
        throughput = tubeloss(dp_ion*1e-9,flow_df.values*1.667e-5,pipe_length,t_df.values,p_df.values)
    if mode=="particles":
        throughput = tubeloss(dp_par*1e-9,flow_df.values*1.667e-5,pipe_length,t_df.values,p_df.values)
    
    df = df / throughput

    return df
 
def wagner_ion_mode_correction(df):
    if df is None:
        return None
    else:
        roberts_corr = 0.713*dp_ion**0.120
        df = df / roberts_corr
     
        return df

def add_flags(
    df,
    rec,
    mode):

    if ((rec is None) or (df is None)):
        return None

    else:        
        # Extract the records that match the mode
        if mode=="ions":
            df_rec = rec[rec.opmode=='ions']
        if mode=="particles":
            df_rec = rec[rec.opmode=='particles']

        if not df_rec.index.to_series().is_monotonic_increasing:
            return None
        
        df_rec = df_rec.reindex(df.index,method="nearest")

        # Read the flags column from records and add it to 
        # the final data as the first column
        df.insert(0,"Flags",df_rec["flags"])

        return df


def clean_elem_noise(
    df,
    rec,
    mode,
    polarity,
    inverter_name):

    if ((df is None) or (rec is None)):
        return None

    if inverter_name =="hires_25":
        if mode=="ions":
            if polarity=="neg":
                elm2dp = ions_neg_v141_hrnd_elm25_chv
            if polarity=="pos":
                elm2dp = ions_pos_v141_hrnd_elm25_chv
        if mode=="particles":
            if polarity=="neg":
                elm2dp = particles_neg_v14_hrnd_elm25_chv
            if polarity=="pos":
                elm2dp = particles_pos_v14_hrnd_elm25_chv
    elif inverter_name == "lores_25":
        if mode=="ions":
            if polarity=="neg":
                elm2dp = ions_neg_v141_lrnd_elm25_chv
            if polarity=="pos":
                elm2dp = ions_pos_v141_lrnd_elm25_chv
        if mode=="particles":
            if polarity=="neg":
                elm2dp = particles_neg_v14_lrnd_elm25_chv
            if polarity=="pos": 
                elm2dp = particles_pos_v14_lrnd_elm25_chv
    elif inverter_name == "lores_21":
        if mode=="ions":
            if polarity=="neg":
                elm2dp = ions_neg_v14_lrnd
            if polarity=="pos":
                elm2dp = ions_pos_v14_lrnd
        if mode=="particles":
            if polarity=="neg":
                elm2dp = particles_neg_v14_lrnd
            if polarity=="pos": 
                elm2dp = particles_pos_v14_lrnd
    else:
        return df

    # Extract the records that match the mode
    if mode=="ions":
        df_rec = rec[rec.opmode=='ions']
    if mode=="particles":
        df_rec = rec[rec.opmode=='particles']

    df_rec = df_rec.reindex(df.index,method="nearest")

    elm2dp = {int(k):v for k,v in elm2dp.items()}
    number_of_elms = len(elm2dp)

    # Rolling time windows
    reso_in_seconds = (df.index[1]-df.index[0]).seconds
    small_window = int((10.*60.)/(reso_in_seconds))          # 10 minutes
    medium_window = int((4.*60.*60.)/(reso_in_seconds))      # 6 hours
    large_window = int((12.*60.*60.)/(reso_in_seconds))      # 12 hours

    # NOISE LEVEL FROM THE RECORDS
    if polarity == "neg":
        df_std = df_rec.iloc[:,2+2*number_of_elms:2+3*number_of_elms]
    if polarity == "pos":
        df_std = df_rec.iloc[:,2+3*number_of_elms:2+4*number_of_elms]
    else:
        return None

    # Set index to electrometer number
    elm_header = np.arange(0,number_of_elms).astype(int)
    df_std.columns = elm_header

    # Calculate noise level at each diameter
    df_std2 = df.copy()

    for d in df.columns.values:
        elms = []
        for elm in df_std.columns.values:
            if ((d >= elm2dp[elm][0]) & (d <= elm2dp[elm][1])):
                elms.append(elm)
        df_std2[d] = df_std[elms].mean(axis=1).values
    
    # Apply medium window to get rid of small fluctuations in electrometer noise
    df_std2 = df_std2.rolling(medium_window, min_periods=int((medium_window+1.)/2.), center=True).median()      
    
    # Get the median noise
    median_std2 = np.nanmedian(df_std2)

    # Then find where the noise is more than N times median 
    N = 500
    df_std3 = df_std2.where((df_std2>N*median_std2), np.nan)
    # NOISE LEVEL FROM THE INVERTED DATA

    # Calculate standard deviation in 10 min segments
    df2 = df.rolling(small_window, min_periods=int((small_window+1.)/2.), center=True).std()

    # In a bigger window (12 hours) calculate the 75th quantile of the standard deviations
    # (semi)continuous noise causes higher values compared to normal and rare sudden changes in conc
    df2 = df2.rolling(large_window, min_periods=int((large_window+1.)/2.), center=True).quantile(0.75)

    # find where the noise is more than M times the median
    M = 7
    threshold = M*np.nanmedian(df2)
    
    df3 = df2.where(df2 > threshold, np.nan)

    # REMOVE DATA FROM WHERE THE ELECTROMETER NOISE AND THE INVERTED DATA NOISE AGREE
    df = df[df3.isna() & df_std3.isna()]
    
    return df


def clean_corona_ions(
    df,
    rec,
    mode):

    if ((df is None) or (rec is None)):
        return None

    # Only consider likely limit range
    lower = 1.5e-9
    upper = 5.0e-9
    c = (lower <= df.columns.values) & (upper >= df.columns.values)
    df2 = df.loc[:, c]
 
    # Find maximum difference between size bin medians
    corona_lim = df2.columns.values[df2.median().diff().abs().argmax()]
 
    # Set values below corona ion limit to NaNs
    df.iloc[:,df.columns.values<=corona_lim]=np.nan

    return df

def choose_particle_polarity(negpar,pospar):

    if ((negpar is None) & (pospar is None)):
        return None
    elif ((negpar is None) & (pospar is not None)):
        return "pos"
    elif ((pospar is None) & (negpar is not None)):
        return "neg"
    else:
        pass

    # Calculate number concentration between 2-3 nm
    # and determine the better polarity based on that
    neg_conc = af.calc_conc(negpar,2e-9,3e-9)
    pos_conc = af.calc_conc(pospar,2e-9,3e-9)

    neg_med = float(neg_conc.median())
    pos_med = float(pos_conc.median())

    if (np.isnan(neg_med) & np.isnan(pos_med)):
        return None
    elif np.isnan(neg_med):
        return "pos"
    elif np.isnan(pos_med):
        return "neg"
    elif (pos_med<neg_med):
        return "pos"
    else:
        return "neg"

def nais_processor(config_file):
    """ Processes NAIS data

    Parameters
    ----------

    config_file : str
        full path to configuration file

    """

    with open(config_file,'r') as stream:
        config = yaml.safe_load(stream)
        load_path = config['data_folder']
        save_path = config['processed_folder']
        start_date = config['start_date']
        database = config['database_file']
        location = config['measurement_location']
        end_date = config['end_date']
        allow_reprocess = config["allow_reprocess"]
        pipelength = config['inlet_length']
        do_inlet_loss_correction = config['do_inlet_loss_correction']
        convert_to_standard_conditions = config['convert_to_standard_conditions']
        do_wagner_ion_mode_correction = config["do_wagner_ion_mode_correction"]
        choose_better_particle_polarity=config["choose_better_particle_polarity"]
        remove_noisy_electrometers = config["remove_noisy_electrometers"]
        remove_corona_ions = config["remove_corona_ions"]
        inverter_name = config["inverter_name"]
        use_default_values = config["use_default_values"]
        default_temperature = config["default_temperature"]
        default_pressure = config["default_pressure"]
        default_flowrate = config["default_flowrate"]
        include_flags = config["include_flags"]

    db = TinyDB(database)
    check = Query()

    assert isinstance(start_date,date)
    assert (end_date=='' or isinstance(end_date,date))
    assert os.path.exists(save_path)
    assert all([os.path.exists(x) for x in load_path])
    assert isinstance(allow_reprocess,bool)
    assert isinstance(remove_corona_ions,bool)
    assert isinstance(remove_noisy_electrometers,bool)
    assert isinstance(convert_to_standard_conditions,bool)
    assert isinstance(do_wagner_ion_mode_correction,bool)
    assert isinstance(do_inlet_loss_correction,bool)
    assert isinstance(choose_better_particle_polarity,bool)
    assert ((inverter_name=="hires_25") | (inverter_name=="lores_25") | (inverter_name=="lores_21") | (inverter_name==''))
    assert (isinstance(pipelength,(float, int)) & (not isinstance(pipelength,bool)))
    assert isinstance(use_default_values,bool)
    assert ((isinstance(default_temperature,(float, int)) & (not isinstance(default_temperature,bool))) | (default_temperature==''))
    assert ((isinstance(default_pressure,(float, int)) & (not isinstance(default_pressure,bool))) |  (default_pressure==''))
    assert ((isinstance(default_flowrate,(float, int)) & (not isinstance(default_flowrate,bool))) | (default_flowrate==''))
    assert isinstance(include_flags,bool)

    end_date = date.today() if end_date=='' else end_date

    db = TinyDB(database)
    check = Query()

    start_dt = pd.to_datetime(start_date)
    end_dt = pd.to_datetime(end_date)

    start_date_str = start_dt.strftime("%Y%m%d")
    end_date_str = end_dt.strftime("%Y%m%d")

    # list existing dates based on if diagnostic file was found
    list_of_existing_dates = [x["timestamp"] for x in db.search(check.diagnostics.exists())]

    if len(list_of_existing_dates)==0:
        print("building database...")
        list_of_datetimes = pd.date_range(start=start_date_str, end=end_date_str)
    else:
        last_existing_date = sorted(list_of_existing_dates)[-1]
        list_of_datetimes = pd.date_range(start=last_existing_date, end=end_date_str)
    
    # Add unprocessed datafiles to the database
    for x in list_of_datetimes:
        if (x.strftime("%Y%m%d") in list_of_existing_dates):
            continue
        else:
            files_found=False
            for z in load_path:
                for y in filename_formats:

                    ion_fn = os.path.join(z,x.strftime(y[0]))
                    particle_fn = os.path.join(z,x.strftime(y[1]))
                    diagnostic_fn = os.path.join(z,x.strftime(y[2]))

                    if ( (os.path.exists(ion_fn) | # ions
                         os.path.exists(particle_fn)) & # particles
                         os.path.exists(diagnostic_fn) # diagnostics
                       ):

                        dtstr = x.strftime("%Y%m%d")

                        db.insert(
                            {"timestamp":dtstr,
                            "diagnostics":diagnostic_fn}
                            )

                        if os.path.exists(ion_fn):
                            db.update(
                                {"ions":ion_fn},
                                check.timestamp==dtstr)

                        if os.path.exists(particle_fn):
                            db.update(
                                {"particles":particle_fn},
                                check.timestamp==dtstr)

                        files_found=True
                        break

                if files_found:
                    break

    # From the database find the last day with processed data
    processed_days = db.search(
        check.processed_neg_ion_file.exists() |
        check.processed_pos_ion_file.exists() |
        check.processed_neg_particle_file.exists() |
        check.processed_pos_particle_file.exists())

    if len(processed_days)!=0:
        last_day=np.max([datetime.strptime(x["timestamp"],"%Y%m%d") for x in processed_days]).strftime("%Y%m%d")
    else:
        last_day=end_date_str

    if allow_reprocess:
        iterator1 = iter(db.search(
         (check.diagnostics.exists() &
          (check.ions.exists() |
          check.particles.exists()) &
          (check.timestamp>=start_date_str) &
          (check.timestamp<=end_date_str))))
    else:
        iterator1 = iter(db.search(
            (check.diagnostics.exists() &
             (check.ions.exists() |
             check.particles.exists()) &
             (check.timestamp>=last_day) &
             (check.timestamp>=start_date_str) &
             (check.timestamp<=end_date_str))
            ))

    for x in iterator1:

        print("processing %s (%s)" % (x["timestamp"],location))

        ions_exist=bool(db.search(
            check.ions.exists() &
            (check.timestamp==x["timestamp"])))

        particles_exist=bool(db.search(
            check.particles.exists() &
            (check.timestamp==x["timestamp"])))

        records,flags = read_file(x["diagnostics"],"records")

        if include_flags:
            my_save_path_flags = os.path.join(save_path,"NAIS"+x["timestamp"]+".flags")
            flags.to_csv(my_save_path_flags,index=False)

        # ions
        if ions_exist:

            ions = read_file(x["ions"],"spectra")

            negion_datamatrix,posion_datamatrix = process_data(ions,"ions")

            if (convert_to_standard_conditions or do_inlet_loss_correction):
                temperature_ion_df,pressure_ion_df,flowrate_ion_df = get_environmental_data(
                    negion_datamatrix,
                    records,
                    "ions",
                    use_default_values,
                    default_pressure,
                    default_temperature,
                    default_flowrate)

            if convert_to_standard_conditions:
                 negion_datamatrix = bring_to_sealevel(negion_datamatrix,temperature_ion_df,pressure_ion_df)
                 posion_datamatrix = bring_to_sealevel(posion_datamatrix,temperature_ion_df,pressure_ion_df)

            if do_inlet_loss_correction:
                 negion_datamatrix = correct_inlet_losses(
                      negion_datamatrix,
                      "ions",
                      pipelength,
                      temperature_ion_df,
                      pressure_ion_df,
                      flowrate_ion_df)
                 posion_datamatrix = correct_inlet_losses(
                      posion_datamatrix,
                      "ions",
                      pipelength,
                      temperature_ion_df,
                      pressure_ion_df,
                      flowrate_ion_df)

            if do_wagner_ion_mode_correction:
                negion_datamatrix = wagner_ion_mode_correction(negion_datamatrix)
                posion_datamatrix = wagner_ion_mode_correction(posion_datamatrix)
 
            if remove_noisy_electrometers: 
                negion_datamatrix = clean_elem_noise(
                       negion_datamatrix,
                       records,
                       "ions",
                       "neg",
                       inverter_name)
 
                posion_datamatrix = clean_elem_noise(
                       posion_datamatrix,
                       records,
                       "ions",
                       "pos",
                       inverter_name)
            
            if (negion_datamatrix is not None):

                if include_flags:
                    negion_datamatrix = add_flags(negion_datamatrix,records,"ions")
                    
                my_save_path_neg=os.path.join(save_path,"NAISn"+x["timestamp"]+"nds.sum")
                negion_datamatrix.to_csv(my_save_path_neg)

                db.update({"processed_neg_ion_file": my_save_path_neg},
                    check.timestamp==x["timestamp"])

            if (posion_datamatrix is not None):

                if include_flags:
                    posion_datamatrix = add_flags(posion_datamatrix,records,"ions")

                my_save_path_pos = os.path.join(save_path,"NAISp"+x["timestamp"]+"nds.sum")
                posion_datamatrix.to_csv(my_save_path_pos)

                db.update({"processed_pos_ion_file": my_save_path_pos},
                    check.timestamp==x["timestamp"])

        # particles
        if particles_exist:

            particles = read_file(x["particles"], "spectra")

            negpar_datamatrix,pospar_datamatrix = process_data(particles,"particles")

            if (convert_to_standard_conditions or do_inlet_loss_correction):
                temperature_particle_df,pressure_particle_df,flowrate_particle_df = get_environmental_data(
                    negpar_datamatrix,
                    records,
                    "particles",
                    use_default_values,
                    default_pressure,
                    default_temperature,
                    default_flowrate)

            if convert_to_standard_conditions:
                 negpar_datamatrix = bring_to_sealevel(negpar_datamatrix,temperature_particle_df,pressure_particle_df)
                 pospar_datamatrix = bring_to_sealevel(pospar_datamatrix,temperature_particle_df,pressure_particle_df)

            if do_inlet_loss_correction:
                 negpar_datamatrix = correct_inlet_losses(
                      negpar_datamatrix,
                      "particles",
                      pipelength,
                      temperature_particle_df,
                      pressure_particle_df,
                      flowrate_particle_df)
                 pospar_datamatrix = correct_inlet_losses(
                      pospar_datamatrix,
                      "particles",
                      pipelength,
                      temperature_particle_df,
                      pressure_particle_df,
                      flowrate_particle_df)

            if remove_noisy_electrometers:
                negpar_datamatrix = clean_elem_noise(
                       negpar_datamatrix,
                       records,
                       "particles",
                       "neg",
                       inverter_name)
 
                pospar_datamatrix = clean_elem_noise(
                       pospar_datamatrix,
                       records,
                       "particles",
                       "pos",
                       inverter_name)
 
            if choose_better_particle_polarity:
                better_polarity = choose_particle_polarity(
                        negpar_datamatrix,
                        pospar_datamatrix)
            else:
                better_polarity = None
 
            if remove_corona_ions:
                negpar_datamatrix = clean_corona_ions(
                    negpar_datamatrix,
                    records,
                    "particles")
 
                pospar_datamatrix = clean_corona_ions(
                    pospar_datamatrix,
                    records,
                    "particles")

            if ((negpar_datamatrix is not None) & 
               ((choose_better_particle_polarity==True) & (better_polarity=="neg") |
               (choose_better_particle_polarity==False))):
                
                if include_flags:
                    negpar_datamatrix = add_flags(negpar_datamatrix,records,"particles")

                my_save_path_neg=os.path.join(save_path,"NAISn"+x["timestamp"]+"np.sum")
                negpar_datamatrix.to_csv(my_save_path_neg)

                db.update({"processed_neg_particle_file": my_save_path_neg},
                    check.timestamp==x["timestamp"])

            if ((pospar_datamatrix is not None) &
               ((choose_better_particle_polarity==True) & (better_polarity=="pos") |
               (choose_better_particle_polarity==False))): 

                if include_flags:
                    pospar_datamatrix = add_flags(pospar_datamatrix,records,"particles")

                my_save_path_pos=os.path.join(save_path,"NAISp"+x["timestamp"]+"np.sum")
                pospar_datamatrix.to_csv(my_save_path_pos)

                db.update({"processed_pos_particle_file": my_save_path_pos},
                    check.timestamp==x["timestamp"])

    print("Done!")


def combine_databases(database_list, combined_database):
    """Combine JSON databases

    If the measurement setup changes one may have to use multiple configuration files
    which results in multiple databases. With this function you can combine the databases
    into a single database after processing.

    Parameters
    ----------

    database_list : str
        List of full paths to databases that should be combined

        First database should have the earliest data, second database
        the second earliest and so on

    combined_database : str
        full path to combined database
    
    """

    DB = {}
    i = 0

    for database in database_list:

        fid=open(database)

        database_json=json.load(fid)

        for key in database_json["_default"]:
            DB[i] = database_json["_default"][key]
            i=i+1

    with open(combined_database, "w") as f:
        json.dump({"_default":DB},f)
