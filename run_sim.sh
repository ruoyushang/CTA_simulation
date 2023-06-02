#!/bin/bash

cd /nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/
export CTADATA=/a/data/serret/shang/1dc
. ${NevisAppBase}/adm/nevis-init.sh
echo "loading module..."
#module load common/miniconda3
module load ctools/2.0.0
#module load gammapy

#echo "conda activate env..."
#conda_setup 
#conda activate cta_sim_env

echo "running the job..."
#python3 hello.py 

exposure='10.0'

#layout='mult3_zd20_m4interp'

#layout='mult3_zd20_f4' # F4-14MSTs
layout='mult3_zd20_m2' # M2-14MSTs11SCTs
#layout='mult3_zd20_c0' # C0-25MSTs
python3 run_icrc2023.py $layout $exposure 'crab_bkg'
python3 run_icrc2023.py $layout $exposure 'crab_src'
python3 run_icrc2023.py $layout $exposure 'two_bkg'
python3 run_icrc2023.py $layout $exposure 'two_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p25kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p25kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p5kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p5kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_1p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_1p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_2p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_2p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_3p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_3p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_4p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_4p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'ic443_bkg'
python3 run_icrc2023.py $layout $exposure 'ic443_src'

layout='mult3_zd20_f4' # F4-14MSTs
#layout='mult3_zd20_m2' # M2-14MSTs11SCTs
#layout='mult3_zd20_c0' # C0-25MSTs
python3 run_icrc2023.py $layout $exposure 'crab_bkg'
python3 run_icrc2023.py $layout $exposure 'crab_src'
python3 run_icrc2023.py $layout $exposure 'two_bkg'
python3 run_icrc2023.py $layout $exposure 'two_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p25kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p25kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p5kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p5kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_1p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_1p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_2p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_2p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_3p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_3p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_4p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_4p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'ic443_bkg'
python3 run_icrc2023.py $layout $exposure 'ic443_src'

#layout='mult3_zd20_f4' # F4-14MSTs
#layout='mult3_zd20_m2' # M2-14MSTs11SCTs
layout='mult3_zd20_c0' # C0-25MSTs
python3 run_icrc2023.py $layout $exposure 'crab_bkg'
python3 run_icrc2023.py $layout $exposure 'crab_src'
python3 run_icrc2023.py $layout $exposure 'two_bkg'
python3 run_icrc2023.py $layout $exposure 'two_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p25kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p25kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p5kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_0p5kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_1p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_1p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_2p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_2p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_3p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_3p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_4p0kpc_bkg'
python3 run_icrc2023.py $layout $exposure 'rx_j1713_4p0kpc_src'
python3 run_icrc2023.py $layout $exposure 'ic443_bkg'
python3 run_icrc2023.py $layout $exposure 'ic443_src'
