from os.path import exists
from special_media import set_special_media

if __name__ == '__main__':
    keff_df = 'kappmax_kcat_in_vitro_lm_rf_dl_iJO_MOM_flux_iv_sub_w_enet'
    include_old = True

    if not include_old:
        cols = kapp_df.columns
    else:
        cols = ['OLD']
        cols.extend(list(kapp_df.columns))

    print(cols)
    for col in cols:
        if column and col not in column:
            continue
        print(col)
        if col != 'OLD':
            update_keffs(me, kapp_df[col], high_pdh=high_pdh)

        me_nlp = ME_NLP1(me, growth_key='mu')
        for media, uptake_rxn in media_list.items():

            print(media, col)
            if check_exist and exists(
                            '%s/%s_%s_sim.pickle' % (saveloc, col, media)):
                continue
            if sim_media and media not in sim_media:
                continue

            print(media)
            me.reactions.EX_glc__D_e.lower_bound = 0
            if media in ['LB', 'Glucose_AA']:
                set_special_media(me, media)
            else:
                me.reactions.get_by_id(uptake_rxn).lower_bound = -1000
            me.reactions.EX_o2_e.lower_bound = -1000

            x, status, hs = me_nlp.solvelp(.0001)

            if status != 'optimal':
                print(status)
                print(col, 'is infeasible')
                with open('%s/%s_%s_sim.pickle' % (saveloc, col, media),
                          'wb') as f:
                    pickle.dump(me.solution, f)
                continue
            print(me.solution.x_dict['biomass_dilution'])

            muopt, hs, xopt, cache = me_nlp.bisectmu(precision=1e-14,
                                                     basis=hs)
            sol = me.solution

            with open('%s/%s_%s_sim.pickle' % (saveloc, col, media),
                      'wb') as f:
                pickle.dump(sol, f)
            with open('%s/%s_%s_sim_flux.json' % (saveloc, col, media),
                      'w') as f:
                json.dump(sol.x_dict, f)
            with open('%s/%s_%s_sim_met_flux.json' % (saveloc, col, media),
                      'w') as f:
                json.dump(me.get_metabolic_flux(solution=sol), f)
            me.reactions.get_by_id(uptake_rxn).lower_bound = 0.
    saveloc = '%s/%s' % (parent, 'ME_sims_')
    run_with_simulations(model, saveloc, check_exist=True,
                         objective_rxn='PDH_FWD_PYRUVATEDEH-CPLX_mod_mg2_'
                                       'mod_fad_mod_thmpp_mod_lipo',
                         old=True, high_pdh=True)
