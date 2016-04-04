import os
import tempfile
import subprocess
import io
import pandas as pd
from ...correction import correct_with_percentage
import numpy as np
import atexit


rserve_started = False


def start_rserve():
    global rserve_started
    if not rserve_started:
        rserve_started = True
        os.system("R CMD Rserve --no-save\\\n")
        atexit.register(stop_rserve)


def stop_rserve():
    try:
        os.system('killall Rserve')
    except:
        pass


def Percentage_Cibersort():
    start_rserve()
    if not os.path.exists('cache/cibersort'):
        os.makedirs('cache/cibersort')

    def calc(df_tpm, genes_for_estimation):
        cibersort_input = df_tpm[['observed']].reset_index()
        cibersort_input.columns = ['index', 'observed']
        tf = tempfile.NamedTemporaryFile(
            suffix='.input.txt', dir='cache/cibersort')
        cibersort_input.sort('index').to_csv(tf, sep="\t", index=False)
        tf.flush()

        tf_sig = tempfile.NamedTemporaryFile(
            suffix='.sig.txt', dir='cache/cibersort')
        sig_df = df_tpm.ix[genes_for_estimation, [
            'reference_sample', 'reference_contamination']]
        sig_df = sig_df.reset_index()
        sig_df.columns = ['index', 'sample', 'contamination']
        sig_df.sort('index').to_csv(tf_sig, sep="\t", index=False)
        tf_sig.flush()

        p = subprocess.Popen("java -jar ../../../code/CIBERSORT/package/CIBERSORT.jar -M %s -B %s" % (os.path.abspath(tf.name),
                                                                                                      os.path.abspath(tf_sig.name)), shell=True, cwd='cache/cibersort', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        result = stdout[stdout.find("==================CIBERSORT====================") + len(
            "==================CIBERSORT===================="):].strip()
        result = io.BytesIO(result)
        print stdout, stderr
        df_result = pd.read_csv(result, sep="\t")
        calc_percentage = df_result.ix[0, 'contamination']
        p_value = df_result.ix[0, 'P-value']
        calc_percentage_std = np.nan
        percentage, corrected = correct_with_percentage(
            df_tpm,
            'observed',
            ['reference_sample'],
            ['reference_contamination'],
            calc_percentage)
        return {
            'percentage': calc_percentage,
            'percentage_std': calc_percentage_std,
            'corrected': corrected,
            'p-value': p_value
        }
    return 'cibersort', calc
