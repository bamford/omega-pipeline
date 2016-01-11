import numpy as np



def in_field(field):
    formato = ('%s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s \
                %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s\
		  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s')
    # This catalogue is created by cutting_catalogue_liners.in_field
    filename = 'F'+str(field)+'_catalogue_master_hst_liners.txt'
    f = np.loadtxt(filename,dtype='float')

    field_number,ID_5,Prob_posiew_5,prob_agn5,prob_sf5,prob_retired5,prob_weak_agn5,\
    cont_maxprob_5,cont_median_5,rel_error_cont_median_5,rel_error_cont_maxlik_5,\
    z_maxprob_5,z_median_5,rel_error_z_median_5,\
    fha_maxprob_5,fha_median_5,rel_error_fha_median_5,rel_error_fha_maxlik_5,lha_median_5,\
    fnii_maxprob_5,fnii_median_5,rel_error_fnii_median_5,rel_error_fnii_maxlik_5,\
    ewha_median_5,ewha_maxlik_5,error_ewha_median_5,error_ewha_maxlik_5,lineratio_median_5,lineratio_maxlik_5,error_lineratio_median_5,error_lineratio_maxlik_5,\
    chisq_5,twolines_5,oneline_5,\
    ID_total,Prob_posiew_total,prob_agn_total,prob_sf_total,prob_retiredtotal,prob_weak_agntotal,\
    cont_maxprob_total,cont_median_total,rel_error_cont_median_total,rel_error_cont_maxlik_total,\
    z_maxprob_total,z_median_total,rel_error_z_median_total,\
    fha_maxprob_total,fha_median_total,rel_error_fha_median_total,rel_error_fha_maxlik_total,lha_median_total,\
    fnii_maxprob_total,fnii_median_total,rel_error_fnii_median_total,rel_error_fnii_maxlik_total,\
    ewha_median_total,ewha_maxlik_total,error_ewha_median_total,error_ewha_maxlik_total,lineratio_median_total,lineratio_maxlik_total,error_lineratio_median_total,error_lineratio_maxlik_total,\
    chisq_total,twolines_total,oneline_total,\
    ST_A_IMAGE,ra,dec,Rmag,combo_flag,stages_flag,MC_z,logmass,logmass_cl,flux24,tir,tuv,tir_cl,tuv_cl,\
    sfr_det,sfr_lo,sfr_hi,sfr_det_cl,sfr_lo_cl,sfr_hi_cl,sed_type = f[:].T


    ### One line aper5
    cut_one_line_5 = np.where((oneline_5 == 1))[0]

    ### One line total aper
    cut_one_line_total = np.where((oneline_total  == 1))[0]

    ###Two lines aper5
    cut_two_lines_5 = np.where((twolines_5  == 1))[0]

    ###Two lines total aper
    cut_two_lines_total = np.where((twolines_total  == 1))[0]

    ###Two lines total aper more restricted
    cut_two_lines_total_01 = np.where((twolines_total  == 1))[0]

    ### Two lines and AGN
    cut_two_lines_agn_5 = np.where((twolines_5  == 1) & (prob_agn5  > 99.7))[0]

    ### Two lines and SF
    cut_two_lines_sf_5 = np.where((twolines_5  == 1) & (prob_sf5  > 99.7))[0]


    ### Catalogue with all the galaxies within the z cut in the total aperture measurement
    np.savetxt('F'+str(field)+'_catalogue_total_master_hst_liners_nocut.txt',f,fmt=formato)


    ### Catalogue with all the galaxies within the two lines cut in the total aperture measurement
    np.savetxt('F'+str(field)+'_catalogue_cut_two_lines_total_master_hst_liners_nocut.txt',f[cut_two_lines_total],fmt=formato)

    ### Catalogue with all the galaxies within the two lines cut in the total aperture measurement MORE RESTRICTED
    np.savetxt('F'+str(field)+'_catalogue_cut_two_lines_total_01_master_hst_liners_nocut.txt',f[cut_two_lines_total_01],fmt=formato)


    ### Catalogue with all the galaxies within the two lines cut in the 5 pixel aperture measurement
    np.savetxt('F'+str(field)+'_catalogue_cut_two_lines_aper5_master_hst_liners_nocut.txt',f[cut_two_lines_5],fmt=formato)

    ### Catalogue with all the galaxies within the one line cut in the total aperture measurement
    np.savetxt('F'+str(field)+'_catalogue_cut_one_line_total_master_hst_liners_nocut.txt',f[cut_one_line_total],fmt=formato)



    f1 = open('number_of_galaxies_hst_liners_nocut.txt','a+')
    f1.write('%2i %3i %3i %3i %3i %3i \n'%(field,len(field_number), len(field_number),len(cut_one_line_total),len(cut_two_lines_total), len(cut_two_lines_5)))
    f1.close()
