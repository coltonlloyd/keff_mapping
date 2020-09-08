import datetime
import os
from multiprocessing import Pool
from os.path import dirname, abspath
import numpy as np

import glob
import cobra
import cobrame
import pickle
import pandas as pd
from cobrame.io.json import load_json_me_model, save_json_me_model

from kappmax_prediction_scripts import new_update_keffs, run_simulations, compare_with_proteomics

media_dict = {'Glucose': 'EX_glc__D_e', 'Acetate': 'EX_ac_e',
              'Pyruvate': 'EX_pyr_e', 'Glycerol': 'EX_glyc_e',
              'Fumarate': 'EX_fum_e', 'Succinate': 'EX_succ_e',
              'LB': '', 'Glucosamine': 'EX_gam_e',
              'Mannose': 'EX_man_e', 'Xylose': 'EX_xyl__D_e',
              'Fructose': 'EX_fru_e', 'Glycerol + AA': '',
              'Galactose': 'EX_gal_e'}

keff_vectors = ["experim_kappmax_NULL",
                "experim_kcat_in_vitro_NULL",
                "kappmax.ensemble_MTF-FBA",
                "kcat.in.vitro.ensemble_MTF-FBA",
                "OLD"]
root = dirname(abspath(__file__))
name_suffix = 'keff_analysis'
if not name_suffix:
    name_suffix = datetime.datetime.now().strftime("%Y-%m-%d_at_%H:%M:%S")

save_prefix = 'kappmax_kcat_in_vitro_lm_enet_rf_dl_iML_MTF-FBA_bm_norm_FALSE_use_Dans_kcats_TRUEuse_met_dG_KmTRUE_iv_sub'
out_location = '%s/all_output_%s' % (root, name_suffix)
batch_sims_save_loc = '%s/batch_simulations/%s/' % (out_location, save_prefix)
validations_save_loc = '%s/validations/%s/' % (out_location, save_prefix)

#resource_dir = dirname(abspath(me_validation_resources.__file__))
proteomics_data_dir = '/'.join([root, 'data'])

os.makedirs(batch_sims_save_loc, exist_ok=True)
os.makedirs(validations_save_loc, exist_ok=True)


def set_keffs_media_and_solve(value):
    keff_column, media = value
    model = load_json_me_model('iJL1678b.json')
    try:
        model.reactions.EX_glc_e.id = 'EX_glc__D_e'
        model.repair()
    except:
        pass
    df = pd.read_csv('/'.join([root, 'data',
                               save_prefix + '.csv']),
                     index_col=0)

    transporters = pd.read_csv('/'.join([root, 'data',
                                         'transport_reactions.csv']),
                               index_col=0)['bigg_id'].values

    if keff_column != 'OLD':
        new_update_keffs.update_all_keffs(model, df[keff_column],
                                          transporters=transporters)
    simulation_savefile_name = \
        '%s/%s_%s_sim_flux.json' % (batch_sims_save_loc, keff_column, media)

    model.metabolites.biomass.compartment = 'mc'

    model.global_info['k_deg'] = 0
    for r in model.reactions:
        if isinstance(r, cobrame.TranslationReaction):
            r.update()
    model.unmodeled_protein_fraction = 0

    run_simulations.maximize_growth_rate(model, media,
                                         simulation_savefile_name,
                                         solver='gurobi',
                                         precision=1e-12)


def run_pool(function, values, processes=2):
    pool = Pool(processes=processes)
    pool.map(function, values)


if __name__ == '__main__':
    run_simulations_flag = True
    submit_list = []
    for keff in keff_vectors:
        for media in media_dict:
            submit_list.append([keff, media])

    if run_simulations_flag:
        run_pool(set_keffs_media_and_solve, submit_list, processes=1)

    model = load_json_me_model('iJL1678b.json')

    for keff in keff_vectors:
        for media in media_dict:

            proteomics_data_path = \
                '/'.join([proteomics_data_dir, 'Aebersold_copy_numbers.xlsx'])
            simulation_savefile_name = \
                '%s/%s_%s_sim_flux.json' % (batch_sims_save_loc, keff, media)

            dataframe = compare_with_proteomics.return_dataframe(model, media,
                                                                 proteomics_data_path,
                                                                 simulation_savefile_name)

            dataframe.to_csv('%s/spreadsheet_%s_%s.csv' % (validations_save_loc, media, keff))

            compare_with_proteomics.transform_df_to_mass_or_mol_fraction(
                dataframe, model, ['Measured',
                                   'Simulated'])


    # ########## Process output for david ###########
    ijo = cobra.io.load_json_model('%s/iJO1366_bigg.json' % proteomics_data_dir)

    tRNA_genes = set()
    for gene in ijo.genes:
        for r in gene.reactions:
            if r.subsystem == 'tRNA Charging':
                print(gene)
                tRNA_genes.add(gene.id)

    print(tRNA_genes)

    sim_loc = '%s/*.csv' % validations_save_loc
    i = 0
    df = pd.DataFrame()
    for f in glob.glob(sim_loc):
        csv = pd.read_csv(f, index_col=0)
        filename = f.split('/spreadsheet_')[-1].replace('.csv', '')
        split = filename.split('_')
        substrate = split[0]
        source = '_'.join(split[1:])
        for gene, v in csv.iterrows():
            if not v['Metabolic']:
                continue
            if gene in tRNA_genes:
                continue
            df.loc[i, 'gene_bnum'] = gene
            df.loc[i, 'abundance (mmol)'] = v['Simulated']
            df.loc[i, 'abundance(mg)'] = \
                v['Simulated'] * model.metabolites.get_by_id('protein_' + gene).formula_weight
            if v['Membrane_complex_associated']:
                df.loc[i, 'annotation'] = 'Membrane'
            elif v['Transport_gene']:
                df.loc[i, 'annotation'] = 'Transport'
            else:
                df.loc[i, 'annotation'] = ''

            df.loc[i, 'sum metabolic reactions'] = v['Metabolic flux']
            df.loc[i, 'growth_condition'] = substrate
            df.loc[i, 'source'] = source
            i += 1

    df.to_csv('simulation_validation_output.csv')
