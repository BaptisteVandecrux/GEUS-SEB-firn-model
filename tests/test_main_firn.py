import __init__
import numpy as np
import pandas as pd
from main_firn import run_main_firn_old, run_GEUS_model_opt, run_main_firn_parallel, run_GEUS_model_old, run_main_firn_new
from joblib import Parallel, delayed

def compare_old_opt_firn():
    print("--- Test firn old running ---")
    site_list = ["KAN_U"]

    for site in site_list:
        print(site)
        filename = "Firn viscosity/Input files/" + site + ".csv"
        (snowc_old, snic_old, slwc_old, tsoil_old, zrfrz_old, 
        rhofirn_old, zsupimp_old, dgrain_old, grndc_old, 
        grndd_old, compaction_old, zrogl_old, Tsurf_out_old,
        pgrndcapc_old, pgrndhflx_old, dH_comp_old, snowbkt_old) = run_GEUS_model_old(site, filename)

    output_old = [snowc_old, snic_old, slwc_old, tsoil_old, zrfrz_old, 
        rhofirn_old, zsupimp_old, dgrain_old, grndc_old, 
        grndd_old, compaction_old, zrogl_old, Tsurf_out_old,
        pgrndcapc_old, pgrndhflx_old, dH_comp_old, snowbkt_old]

    # Run new main_firn, parallel
    res = Parallel(n_jobs=8, verbose=10)(delayed(run_GEUS_model_opt)(site) for site in site_list)
    output_opt = res[0]

    diff_exist = False
    # Assert that the values from original and optimized code are the same
    for i in range(len(output_opt)):
        if (output_old[i] != output_opt[i]).all():
            print("Difference in objects for i: ", i)
            diff_exist = True

    if diff_exist == False:
        print("Test done. No differences exist.")
   

    return()


compare_old_opt_firn()