
# coding: utf-8

# In[ ]:

class CheckReplication(object):
    ext_ref_df = 'mortality_dataframe_cluster0'
    ext_ref_timegrid = 'relative_timegrid_cluster0'
    
    ext_df = '_mortality_dataframe.pickle'
    ext_meta='_analysis_metadata.pickle'
    ext_timeseries='_timeseries_dataframe.pickle'
    cutoff31 = 1.949
    cutoff32 = 1.858
    cutoff35 = 1.731
    cutoff21 = 1.628
    cutoff22 = 1.517
    cutoff25 = 1.358
    min_files1 = 2
    min_files2 = 2
    min_cells1 = 300
    min_cells2 = 150
    
    pre='CR-'
    endDendro ='dentrogram-all-files.eps'
    endDendroCluster = 'dentrogram-cluster.eps'
    endTsImage='colomap_all_files_cluster'
    extTs ='.png'
    
    
    def __init__(self, folderPath='',refFolderPath='',cutoff=-1,showProgress=True):
        self.folderPath = folderPath
        self.refFolderPath = refFolderPath
        self.cutoff_c = cutoff
        self.showProgress = showProgress
        
    
    def run(self):
        if len(self.refFolderPath):
            self.df_ref,self.timegrid_ref = CheckReplication.getRefDf_timegrid(path_folder=self.refFolderPath)
            self.getStatsFromRefDf_timegrid(self.df_ref,self.timegrid_ref)
            self.vsRef = True
        self.getDf_timegrid_timeseries(self.folderPath,showProgress = True)
        
  
    
    def getStatsFromRefDf_timegrid(self,df_ref, timegrid_ref):
        import sys,os
        sys.path.insert(1, os.path.join(sys.path[0], '..'))
        sys.path.insert(1, os.path.join(sys.path[0], '.'))
        from cohortMortalitySummary import KSm_2samples, KSm_test, generate_ts_image
        from cohortMortalitySummary import KaplanMeier, NelsonAalen, BSHazardR
        from cohortMortalitySummary import GGMfit,  KSm_gof, KSm_test
        from cohortMortalitySummary import loglik_test, GGM_test_2samples
        
        from cohortMortalitySummary import ro, globalenv, pandas2ri, Formula, rflexsurv, pth
        self.KM_ref = KaplanMeier(df_ref, timegrid_ref)
        self.NA_ref = NelsonAalen(df_ref, timegrid_ref)
        self.BSFit_ref = BSHazardR(df_ref)
        self.GGMfit_ref = GGMfit(df_ref)
        self.residues_ref = KSm_gof(df_ref, timegrid_ref, self.GGMfit_ref['ML_survivorship'])
        
    
    

    @staticmethod
    def getSignificance(files,dfs_by_pos, tss_by_pos):
        import sys,os
        sys.path.insert(1, os.path.join(sys.path[0], '..'))
        sys.path.insert(1, os.path.join(sys.path[0], '.'))
        from cohortMortalitySummary import KSm_2samples, KSm_test
        from scipy.special import binom
        import numpy as np
        distances = np.zeros((len(files),len(files)))
        significances = np.empty((len(files),len(files)))
        pdist = np.zeros(int(binom(len(files), int(2))))
        significances[:] = np.NAN
        significances_pd = np.zeros(int(binom(len(files), int(2))))
        k = 0
        for i in range(len(files)):
            for j in range(i+1,len(files)):
                print files[i], '------',files[j]
                residues = KSm_2samples(dfs_by_pos[files[i]], tss_by_pos[files[i]], 
                                        dfs_by_pos[files[j]], tss_by_pos[files[j]])
                test = KSm_test(residues[0],residues[1])
                distances[i][j] = test[0]
                distances[j][i] = test[0]
                if type(test[1]) is str:
                    significances[i][j] = 0
                elif test[1] > 0.05:
                    significances[i][j] = 0
                elif test[1] > 0.02:
                    significances[i][j] = 1
                else:
                    significances[i][j] = 2
                significances[j][i] = significances[i][j]
                pdist[k] = distances[i][j]
                significances_pd = significances[i][j]
                k+=1
        return significances, pdist, significances_pd
    
    
    def getDf_timegrid_timeseries(self,path_folder,showProgress=False):
        import pickle
        df_paths=CheckReplication.listFilePathsByExt(path_folder,CheckReplication.ext_df)
        meta_paths = CheckReplication.listFilePathsByExt(path_folder,CheckReplication.ext_meta)
        timeseries_paths= CheckReplication.listFilePathsByExt(path_folder,CheckReplication.ext_timeseries)
        if showProgress:
            print 'df_paths: ',df_paths
            print 'meta_paths: ',meta_paths
            print 'timeseries_paths: ',timeseries_paths            
        out_path_folder_s, out_name_s = CheckReplication.getSavePathsAndCheck(df_paths, meta_paths)
        if showProgress:
            print '-----------'
            print 'out_path_folder ',out_path_folder_s,'\n'
            print 'out_name_s ', out_name_s
        
            
            
        df_s, files_s, dfs_by_pos_s=CheckReplication.loadDf_timeseries_fromFiles(df_paths,df=True)
        profile_s,noname0,noname1 = CheckReplication.loadDf_timeseries_fromFiles(timeseries_paths,df=False)
        metadata_s = CheckReplication.loadMeta_fromFiles(meta_paths)    
        tss_by_pos_s = CheckReplication.getTss_by_pos(df_paths, meta_paths)
        significances_s, pdist_s, significances_pd_s = CheckReplication.getSignificanceFromList(files_s, dfs_by_pos_s, tss_by_pos_s)
        
        self.files_s = files_s
        self.pdist_s = pdist_s
        self.df_paths = df_paths
        self.profile_s = profile_s
        self.df_s = df_s
        self.out_path_folder_s = out_path_folder_s
        self.out_name_s = out_name_s
        CheckReplication.runDentrogramFromList(files_s,pdist_s,out_path_folder_s, out_name_s,self.cutoff_c)
        self.runClusterFromList(files_s, pdist_s,df_s, tss_by_pos_s,profile_s,out_path_folder_s, out_name_s)
    
    @staticmethod
    def plotStatsRefvsClusters2(GGMfit_ref,GGMfits,output_clusters,name, savePath):
        from matplotlib import pyplot as plt
        import seaborn as sns
        import numpy as np
        f, (ax1, ax3, ax2) = plt.subplots(1, 3, sharex=False, sharey=False)
        color_dict = [2,0,1]
        fontsize1 = 16
        fontsize2 = 14
        clr = sns.color_palette()[0]
        exp = len(output_clusters)-1
        try:
            ax1.bar(0, GGMfit_ref['model_paras'].loc['rate','est'], 1, color=clr, alpha=0.8,
                        yerr=GGMfits[exp]['model_paras'].loc['rate','se']*1.96)
            ax2.bar(0, GGMfit_ref['model_paras'].loc['s','est'], 1, color=clr, alpha=0.8,
                        yerr=GGMfits[exp]['model_paras'].loc['s','se']*1.96)
            ax3.bar(0, GGMfit_ref['model_paras'].loc['beta','est'], 1, color=clr, alpha=0.8,
                        yerr=GGMfits[exp]['model_paras'].loc['beta','se']*1.96)
        except:
            print "plotStatsRefvsClusters2 : ERROR"
        for (exp,ttl,row) in zip(range(len(output_clusters)),
                                    [name+' cluster'+str(i) for i in range(len(output_clusters))],
                                    range(1,len(output_clusters)+1)):
                if row < len(sns.color_palette()):
                    clr = sns.color_palette()[row]
                else:
                    flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71","#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
                    clr = sns.color_palette(flatui)[row-len(sns.color_palette())]
                    
                ax1.bar(row, GGMfits[exp]['model_paras'].loc['rate','est'], 1, color=clr, alpha=0.8,
                        yerr=GGMfits[exp]['model_paras'].loc['rate','se']*1.96)
                ax2.bar(row, GGMfits[exp]['model_paras'].loc['s','est'], 1, color=clr, alpha=0.8,
                        yerr=GGMfits[exp]['model_paras'].loc['s','se']*1.96)
                ax3.bar(row, GGMfits[exp]['model_paras'].loc['beta','est'], 1, color=clr, alpha=0.8,
                        yerr=GGMfits[exp]['model_paras'].loc['beta','se']*1.96)
        ax1.set_xticks(np.arange(0.25,len(output_clusters),1))
        #ax1.set_xticklabels([r'$\Delta$rpoS','Wildtype',r'$\Delta$rssB'],fontsize=fontsize2)
        ax2.set_xticks(np.arange(0.25,len(output_clusters),1))
        #ax2.set_xticklabels([r'$\Delta$rpoS','Wildtype',r'$\Delta$rssB'],fontsize=fontsize2)
        ax3.set_xticks(np.arange(0.25,len(output_clusters),1))
        #ax3.set_xticklabels([r'$\Delta$rpoS','Wildtype',r'$\Delta$rssB'],fontsize=fontsize2)
        ax2.set_yscale('log')
        ax3.set_yscale('log')
        ax3.set_ylim([10,5000])
        ax1.set_title('b: Ageing rate',fontsize=fontsize1)
        ax1.set_ylabel('$b$',fontsize=fontsize1)
        ax2.set_title(r'$s\equiv h_{max}/b$: hazard plateaux',fontsize=fontsize1)
        ax2.set_ylabel('$s$',fontsize=fontsize1)
        ax3.set_title(r'$\beta \equiv \frac{h_{max}}{h_0}$: ageing range',fontsize=fontsize1)
        ax3.set_ylabel(r'$\beta$',fontsize=fontsize1)
        f.set_size_inches(15,4.5)
        sns.despine(fig=plt.gcf())
        plt.tight_layout()
        f.savefig(savePath)
        
    @staticmethod
    def writeGGMtoFile(GGMfits,output_clusters,dfs_by_c,name,savePath):
        from cohortMortalitySummary import GGM_test_2samples
        fo = open(savePath,'w',0)
        for i in range(len(output_clusters)):
            print >> fo, name+' cluster'+str(i)+' #'+str(output_clusters[i][1])
            print >> fo, GGMfits[i]['model_paras']
            print >> fo, '\n'
        
        error = False
        if len(output_clusters) > 1:
            GGM_self_tests = dict()
            for (exp,ttl) in zip(range(1,len(output_clusters)),
                                    [name+' cluster'+str(i) for i in range(1,len(output_clusters))]):
                (c, n, nf, fns) = output_clusters[exp]
                try:
                    GGM_self_tests[ttl] = GGM_test_2samples(dfs_by_c[c],dfs_by_c[output_clusters[0][0]])
                except:
                    print "Error in writeGGMtoFile"
                    error = True
                    break
                    
            
            if not error:
                for ttl, test_res in GGM_self_tests.items():
                    print >> fo, ttl+':'
                    print >> fo, '-------'
                    print >> fo, test_res['bestHypothesisByAIC']
                    print >> fo, '-------'
                    print >> fo, test_res['HypothesesTested'].loc[:,['parameter','p-value']]
                    print >> fo, '\n'
        
        print >> fo, output_clusters
        
    @staticmethod
    def plotStatsRefvsClusters(timegrid_ref,KM_ref,NA_ref,BSFit_ref,GGMfit_ref,residues_ref,output_clusters,GGMfits,KMs,NAs,BSFits,tss_by_c,GGMgofs,name,savePath):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_context("paper")
        sns.set_style("white")
        f, (ax1, ax4, ax2, ax3) = plt.subplots(4, 1, sharex=True, sharey=False)
        fontsize1 = 16
        fontsize2 = 14

        clr = sns.color_palette()[0]; ttl = 'wildtype';
        ax1.fill_between(KM_ref.index,KM_ref.loc[:,'lower_ci'],KM_ref.loc[:,'upper_ci'],
                                    alpha=0.5,label=ttl,color=clr)
        ax2.fill_between(NA_ref.index,NA_ref.loc[:,'lower_ci'],NA_ref.loc[:,'upper_ci'],
                                    alpha=0.5,label=ttl,color=clr)
        ax3.fill_between(BSFit_ref['time'],BSFit_ref['lower.ci'],BSFit_ref['upper.ci'],
                                    alpha=0.5,label=ttl,color=clr)
        ax1.plot(timegrid_ref[1:],GGMfit_ref['ML_survivorship'](timegrid_ref[1:]),linewidth=2, linestyle='--',
                                label=ttl,color=clr)
        ax2.plot(timegrid_ref[1:],GGMfit_ref['ML_cumulative_hazard'](timegrid_ref[1:]),linewidth=2, linestyle='--',
                                label=ttl,color=clr)
        ax3.plot(timegrid_ref[1:],GGMfit_ref['ML_hazard'](timegrid_ref[1:]),linewidth=1.5, linestyle='--',
                                label=ttl,color=clr)
        ax4.plot(residues_ref[0].index, residues_ref[0].values, linewidth=2,
                                color=clr,alpha=0.6)


        for (exp,ttl,row) in zip(range(len(output_clusters)),
                                    [name+' cluster'+str(i) for i in range(len(output_clusters))],
                                    range(len(output_clusters))):
            #clr = sns.color_palette()[row+1]
            if row+1 < len(sns.color_palette()):
                clr = sns.color_palette()[row+1]
            else:
                flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71","#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
                clr = sns.color_palette(flatui)[row+1-len(sns.color_palette())]
                
            #ax1.plot(KMs[exp].index,KMs[exp].loc[:,'survivorship'],linewidth=0.75,label='',color=clr)
            ax1.fill_between(KMs[exp].index,KMs[exp].loc[:,'lower_ci'],KMs[exp].loc[:,'upper_ci'],
                                    alpha=0.5,label=ttl,color=clr)

            #ax2.plot(NAs[exp].index,NAs[exp].loc[:,'cumulative_hazard'],linewidth=0.75,label='',color=clr)
            ax2.fill_between(NAs[exp].index,NAs[exp].loc[:,'lower_ci'],NAs[exp].loc[:,'upper_ci'],
                                    alpha=0.5,label=ttl,color=clr)

            #ax3.plot(BSFits[exp]['time'],BSFits[exp]['hazard'],linewidth=0.75,color=clr)
            ax3.fill_between(BSFits[exp]['time'],BSFits[exp]['lower.ci'],BSFits[exp]['upper.ci'],
                                    alpha=0.5,label=ttl,color=clr)

            ax1.plot(tss_by_c[exp][1:],GGMfits[exp]['ML_survivorship'](tss_by_c[exp][1:]),linewidth=2, linestyle='--',
                                label=ttl,color=clr)

            ax2.plot(tss_by_c[exp][1:],GGMfits[exp]['ML_cumulative_hazard'](tss_by_c[exp][1:]),linewidth=2, linestyle='--',
                                label=ttl,color=clr)

            ax3.plot(tss_by_c[exp][1:],GGMfits[exp]['ML_hazard'](tss_by_c[exp][1:]),linewidth=1.5, linestyle='--',
                                label=ttl,color=clr)

            ax4.plot(GGMgofs[exp][0].index,GGMgofs[exp][0].values,linewidth=2,
                                color=clr,alpha=0.6)


        ax4.plot([0,125],[CheckReplication.cutoff22,CheckReplication.cutoff22],linewidth=2, linestyle='--',
                                color='k',label='alpha=0.01 thresholds',alpha=0.6)
        ax4.plot([0,125],[-CheckReplication.cutoff22,-CheckReplication.cutoff22],linewidth=2, linestyle='--',
                                color='k',alpha=0.6)
        ax4.plot([0,125],[CheckReplication.cutoff25,CheckReplication.cutoff25],linewidth=1, linestyle='--',
                                color='k',label='alpha=0.025 thresholds',alpha=0.6)
        ax4.plot([0,125],[-CheckReplication.cutoff25,-CheckReplication.cutoff25],linewidth=1, linestyle='--',
                                color='k',alpha=0.6)
        ax1.legend(loc=1,fontsize=fontsize2)
        ax2.legend(loc=4,fontsize=fontsize2)
        ax3.legend(loc=4,fontsize=fontsize2)
        ax1.set_ylim(0,1)
        ax2.set_ylim(0.003,10)
        ax3.set_xlim([0,125])
        ax3.set_xlabel('Time (h)',fontsize=fontsize1)
        ax2.set_yscale('log')
        ax3.set_yscale('log')
        ax1.set_ylabel('$S(t)$',fontsize=fontsize1)
        ax1.set_title('Survivorship',fontsize=fontsize1)
        ax2.set_ylabel('$H(t)$',fontsize=fontsize1)
        ax2.set_title('Cumulative hazard rate',fontsize=fontsize1)
        ax3.set_ylabel(r'$\hat{h}(t)$',fontsize=fontsize1)
        ax3.set_title('Smoothed hazard rate',fontsize=fontsize1)
        f.set_size_inches(f.get_size_inches()[0],f.get_size_inches()[1]*4)
        sns.despine(fig=f)
        plt.tight_layout()
        f.savefig(savePath)

    @staticmethod
    def getSavePathsAndCheck_abandoned(path_s_0, path_s_1):
        path_state = True
        folder_s = []
        name_s = []
        out_path_folder_s = []
        out_name_s = []
        import os
        for idx, val in enumerate(path_s_0):
            p_0 = path_s_0[idx][0]
            p_1 = path_s_1[idx][1]
            
            f_0 = os.path.split(p_0)[0]
            f_1 = os.path.split(p_1)[0]
            
            n_0 = os.path.split(p_0)[1]
            n_1 = os.path.split(p_1)[1]
            common = CheckReplication.commonStringFromstart(n_0,n_1)
            if f_0 != f_1 or len(common)==0:
                raise Exception("Invalid filepath!", p_0) 
                break
            out_path_folder_s.append(f_0)
            out_name_s.append(common)
        return out_path_folder_s, out_name_s
    
    @staticmethod
    def getSavePathsAndCheck(path_s_0, path_s_1):
        path_state = True
        folder_s = []
        name_s = []
        out_path_folder_s = []
        out_name_s = []
        import os
        for idx, val in enumerate(path_s_0):
            p_0 = path_s_0[idx][0]
            p_1 = path_s_1[idx][0]
            
            f_0 = os.path.split(p_0)[0]
            f_1 = os.path.split(p_1)[0]
            
            n_0 = os.path.split(p_0)[1]
            n_1 = os.path.split(p_1)[1]
            
            common = CheckReplication.commonStringFromstart(n_0,n_1)
            
            if f_0 != f_1 or len(common)==0:
                raise Exception("Invalid filepath!", p_0) 
                break
            out_path_folder_s.append(f_0)
            out_name_s.append(common)
        return out_path_folder_s, out_name_s
    
    @staticmethod          
    def commonStringFromstart(string1, string2):
        """
        Alert : non common utility,
        """
        
        len1, len2 = len(string1), len(string2)
        r1 = string1
        r2 = string2
        match=''
        for j in range(len(r2)):
            if r1[j]==r2[j]:
                match=match+r1[j]
            else:
                break
        return match
                          
    
    
    @staticmethod
    def runDentrogramFromList(files_s,pdist_s,path_folder_s, name_s,cutoff_c):
        import os
        for idx, val in enumerate(files_s):
            files = val
            folder_path = path_folder_s[idx]
            print "dendrogram will be saved to : ", folder_path
            pdist = pdist_s[idx]
            name = CheckReplication.pre + name_s[idx]+CheckReplication.endDendro
            save_path = os.path.join(folder_path, name)
            CheckReplication.runDentrogram(files, pdist,save_path,cutoff_c)
    
    
  
    @staticmethod
    def runDentrogram(files, pdist,savePath,cutoff_c):
        from cohortMortalitySummary import _criticalPartialBrownianBridge
        import matplotlib.pyplot as plt
        from scipy.cluster.hierarchy import complete, dendrogram
        import re
        fig = plt.figure()
        if cutoff_c < 0:
            cutoff = _criticalPartialBrownianBridge(0.05/float(len(pdist)))
        else:
            cutoff = cutoff_c
        dn = dendrogram(complete(pdist),labels=[re.findall('Position\s(\d+)',fn)[0] for fn in files], color_threshold=cutoff)
        fig.savefig(savePath)
        print "save eps to : ", savePath
        
    def runClusterFromList(self,files_s, pdist_s,df_s,tss_by_pos_s,profile_s,path_folder_s, name_s):
        for idx, val in enumerate(files_s):
            files = val
            pdist = pdist_s[idx]
            df = df_s[idx]
            tss_by_pos = tss_by_pos_s[idx]
            profile = profile_s[idx]
            out_folder = path_folder_s[idx]
            name = name_s[idx]
            self.runCluster(files,pdist,df, tss_by_pos,profile,out_folder,name)
            
    def runCluster(self,files, pdist,df, tss_by_pos,profile,out_folder,name,showProgress=True):
        from cohortMortalitySummary import _criticalPartialBrownianBridge
        from scipy.cluster.hierarchy import average, complete, dendrogram, fcluster
        from cohortMortalitySummary import KSm_2samples, KSm_test,generate_ts_image
        from scipy.cluster.hierarchy import dendrogram
        from cohortMortalitySummary import KaplanMeier, NelsonAalen, BSHazardR
        from cohortMortalitySummary import GGMfit,  KSm_gof, KSm_test
        import matplotlib.pyplot as plt
        import os
        import numpy as np
        import pandas as pd
        from scipy.special import binom
        
        if self.cutoff_c < 0:
            cutoff = _criticalPartialBrownianBridge(0.05/float(len(pdist)))
        else:
            cutoff = self.cutoff_c
        
        fclusters = fcluster(complete(pdist),cutoff,criterion='distance')
        nc = np.max(fclusters)
        fclusters_files = [[files[i] for i in range(len(files)) if fclusters[i]==c] for c in range(1,nc+1)]
        
        clustersFixed = False
        while clustersFixed is False:
            dfs_by_c = [pd.concat([df.loc[fn] for fn in cfiles]) for cfiles in fclusters_files]
            tss_by_c = [list(np.unique(np.sort(np.concatenate([tss_by_pos[fn] for fn in cfiles])))) for cfiles in fclusters_files]

            #tss_by_c = [list(np.sort(np.concatenate([tss_by_pos[fn] for fn in cfiles]))) for cfiles in fclusters_files]
            cdistances = np.zeros((nc,nc))
            csignificances = np.empty((nc,nc))
            cpdist = np.zeros(int(binom(nc, int(2))))
            csignificances[:] = np.NAN
            csignificances_pd = np.zeros(int(binom(nc, int(2))))
            k = 0
            for i in range(nc):
                for j in range(i+1,nc):
                    residues = KSm_2samples(dfs_by_c[i], tss_by_c[i], dfs_by_c[j], tss_by_c[j])
                    test = KSm_test(residues[0],residues[1])
                    cdistances[i][j] = test[0]
                    cdistances[j][i] = test[0]
                    if type(test[1]) is str:
                        csignificances[i][j] = 0
                    elif test[1] > 0.05:
                        csignificances[i][j] = 0
                    elif test[1] > 0.02:
                        csignificances[i][j] = 1
                    else:
                        csignificances[i][j] = 2
                    csignificances[j][i] = csignificances[i][j]
                    cpdist[k] = cdistances[i][j]
                    csignificances_pd = csignificances[i][j]
                    k+=1

            superclusters = fcluster(complete(cpdist),cutoff,criterion='distance')
            if nc == np.max(superclusters):
                clustersFixed = True
            else:
                print nc, superclusters
                nc = np.max(superclusters)
                fclusters_files = [[files[i] for i in range(len(files)) if fclusters[i]-1 < len(superclusters) and superclusters[fclusters[i]-1]==c] for c in range(1,nc+1)]
        if showProgress:
            print fclusters_files
            dn = dendrogram(average(cpdist), color_threshold=cutoff)
        dfs_by_c = [pd.concat([df.loc[fn].set_index(['fn','index'],drop=False) for fn in cfiles]) for cfiles in fclusters_files]
        #tss_by_c = [list(np.sort(np.concatenate([tss_by_pos[fn] for fn in cfiles]))) for cfiles in fclusters_files]
        tss_by_c = [list(np.unique(np.sort(np.concatenate([tss_by_pos[fn] for fn in cfiles])))) for cfiles in fclusters_files]

        for fn in np.unique(profile.index.get_level_values(0)):
            profile.loc[fn,'fn']=fn
            
        profile.loc[:,'index'] = profile.index.get_level_values(1)
        profiles_by_c = [pd.concat([profile.loc[fn].set_index(['fn','index'],drop=True) for fn in cfiles]) for cfiles in fclusters_files]
        
        output_clusters = []
        ln = ()
        for i in range(len(fclusters_files)):
            if len(fclusters_files[i]) < CheckReplication.min_files1:
                ln = ln + tuple(fclusters_files[i])
            elif len(dfs_by_c[i]) < CheckReplication.min_cells2:
                ln = ln + tuple(fclusters_files[i])
            else:
                output_clusters.append((i,len(dfs_by_c[i]), len(fclusters_files[i]), fclusters_files[i]))
        output_clusters = sorted(output_clusters,key=lambda tpl:tpl[1],reverse=True)  
            
        for i in range(len(output_clusters)):
            import pickle
            (c, n, nf, fns) = output_clusters[i]
            image = generate_ts_image(profiles_by_c[c], dfs_by_c[c])
            fig = plt.figure()
            plt.imshow(image, aspect=100/float(len(dfs_by_c[c])), interpolation='None')
            plt.title(CheckReplication.pre+name+' cluster'+str(i))
            plt.xlabel('Frame #')
            plt.ylabel('Cell #')
            fig.set_size_inches(10,10)
            plt.tight_layout()
            fileName = CheckReplication.pre+name+CheckReplication.endTsImage+str(i)+CheckReplication.extTs
            savePath = os.path.join(out_folder,fileName)
            print savePath
            fig.savefig(savePath)
            
            # save pickle original clusters into pickles
            fileName = CheckReplication.pre+name + 'mortality_orig_dataframe_cluster'+str(i)+'.pickle'
            savePath = os.path.join(out_folder,fileName)
            print savePath
            pickle.dump(dfs_by_c[c],open(savePath,'w'))
            
            fileName = CheckReplication.pre+name + 'relative_orig_timegrid_cluster'+str(i)+'.pickle'
            savePath = os.path.join(out_folder,fileName)
            print savePath
            pickle.dump(tss_by_c[c],open(savePath,'w'))
            
            #fig.savefig(os.path.join(strain_pth, exp_prefix+' colomap_all_files_cluster'+str(i)+'.png'))
            
        KMs = []
        NAs = []
        BSFits = []
        for i in range(len(output_clusters)):
            KMs.append(KaplanMeier(dfs_by_c[output_clusters[i][0]], tss_by_c[output_clusters[i][0]]) )
            NAs.append(NelsonAalen(dfs_by_c[output_clusters[i][0]], tss_by_c[output_clusters[i][0]]) )
            BSFits.append( BSHazardR(dfs_by_c[output_clusters[i][0]]) )

        GGMfits = []
        GGMgofs = []
        for i in range(len(output_clusters)):
            GGMfits.append(GGMfit(dfs_by_c[output_clusters[i][0]]))
            residues = KSm_gof(dfs_by_c[output_clusters[i][0]], tss_by_c[output_clusters[i][0]],GGMfits[i]['ML_survivorship'])
            test = KSm_test(residues[0],residues[1],alpha=0.02)
            GGMgofs.append((residues[0],residues[1],test[1],test[2])) 
            
        if self.vsRef:
            fileName = CheckReplication.pre+name+'vs-Ref-survivor.png'
            savePath = os.path.join(out_folder,fileName)
            CheckReplication.plotStatsRefvsClusters(self.timegrid_ref,self.KM_ref,self.NA_ref,self.BSFit_ref,self.GGMfit_ref,self.residues_ref,output_clusters,GGMfits,KMs,NAs,BSFits,tss_by_c,GGMgofs,name,savePath)
            
            fileName = CheckReplication.pre+name+'vs-Ref-bin.png' 
            savePath = os.path.join(out_folder,fileName)
            CheckReplication.plotStatsRefvsClusters2(self.GGMfit_ref,GGMfits,output_clusters,name, savePath)
            
            fileName = CheckReplication.pre+name+'vs-Ref-GGM.txt' 
            savePath = os.path.join(out_folder,fileName)
            CheckReplication.writeGGMtoFile(GGMfits,output_clusters, dfs_by_c,name,savePath)
        

        dfs_by_c_new,tss_by_c_new,profiles_by_c_new,output_clusters_new = CheckReplication.getNewDfsProfilesby_c(files, pdist,df, profile, tss_by_pos,output_clusters,name,out_folder,self.cutoff_c)
        KMs_new,NAs_new,BSFits_new,GGMfits_new,GGMgofs_new = CheckReplication.redefineStatsByCluster(dfs_by_c_new,tss_by_c_new,output_clusters_new)
        
       
        if self.vsRef:
            fileName = CheckReplication.pre+name+"vs-Ref-survivor-new.png"
            savePath = os.path.join(out_folder,fileName)
            CheckReplication.plotStatsRefvsClusters(self.timegrid_ref,self.KM_ref,self.NA_ref,self.BSFit_ref,self.GGMfit_ref,self.residues_ref,output_clusters_new,GGMfits_new,KMs_new,NAs_new,BSFits_new,tss_by_c_new,GGMgofs_new,name,savePath)
  
            fileName = CheckReplication.pre+name+'vs-Ref-bin-new' # without extension, to be added in function
            savePath = os.path.join(out_folder,fileName)
            CheckReplication.plotStatsRefvsClusters2(self.GGMfit_ref,GGMfits_new,output_clusters_new,name, savePath)

            fileName = CheckReplication.pre+name+'vs-Ref-GGM-new.txt'
            savePath = os.path.join(out_folder,fileName)
            lifespan_new = CheckReplication.getLifespan(KMs_new)
            CheckReplication.writeLifespan(GGMfits_new, output_clusters_new, lifespan_new,name,savePath)
        
    @staticmethod
    def getLifespan(KMs):
        import pandas as pd
        lifespan = dict();
        for i in range(len(KMs)):
            upper=KMs[i].loc[:,'upper_ci']
            lower = KMs[i].loc[:,'lower_ci']
            var =  KMs[i].loc[:,'variance']
            survivor=KMs[i].loc[:,'survivorship']
            s = survivor.values >= 0.5    
            slower = lower.values >= 0.5
            supper = upper.values >= 0.5
            ind = survivor.iloc[s].index[-1];  
            indlower = lower.iloc[slower].index[-1]
            indupper = upper.iloc[supper].index[-1]
            data = [{'est': ind, 'L95%': indlower, 'U95%':indupper, 'var':var.loc[ind]}]
            lifespan[i] = pd.DataFrame(data,index=['lifespan'],columns=['est','L95%','U95%','var'])
        print lifespan
        return lifespan
        
    @staticmethod
    def writeLifespan(GGMfits,output_clusters,lifespan,name,savePath):
        fo = open(savePath,'w',0)
        for i in range(len(output_clusters)):
            print >> fo, name+' cluster'+str(i)+' #'+str(output_clusters[i][1])
            print >> fo, GGMfits[i]['model_paras']
            print >> fo, lifespan[i]
            print >> fo, '\n'
    
    @staticmethod
    def redefineStatsByCluster(dfs_by_c, tss_by_c,output_clusters):
        import sys,os
        sys.path.insert(1, os.path.join(sys.path[0], '..'))
        sys.path.insert(1, os.path.join(sys.path[0], '.'))
        from cohortMortalitySummary import KaplanMeier, NelsonAalen, BSHazardR
        from cohortMortalitySummary import GGMfit,  KSm_gof, KSm_test
        KMs = []
        NAs = []
        BSFits = []
        for i in range(len(output_clusters)):
            KMs.append(KaplanMeier(dfs_by_c[output_clusters[i][0]], tss_by_c[output_clusters[i][0]]) )
            NAs.append(NelsonAalen(dfs_by_c[output_clusters[i][0]], tss_by_c[output_clusters[i][0]]) )
            BSFits.append( BSHazardR(dfs_by_c[output_clusters[i][0]]) )

        GGMfits = []
        GGMgofs = []
        for i in range(len(output_clusters)):
            GGMfits.append(GGMfit(dfs_by_c[output_clusters[i][0]]))
            residues = KSm_gof(dfs_by_c[output_clusters[i][0]], tss_by_c[output_clusters[i][0]],GGMfits[i]['ML_survivorship'])
            test = KSm_test(residues[0],residues[1],alpha=0.02)
            GGMgofs.append((residues[0],residues[1],test[1],test[2] )) 
        return KMs,NAs,BSFits,GGMfits,GGMgofs
        
        
        
    
    
    @staticmethod
    def getNewDfsProfilesby_c(files,pdist,df,profile, tss_by_pos,output_clusters,name,out_folder,cutoff_c):
        from cohortMortalitySummary import _criticalPartialBrownianBridge
        import matplotlib.pyplot as plt
        import re
        import pandas as pd
        import numpy as np
        from scipy.cluster.hierarchy import average, complete, dendrogram, fcluster
        from cohortMortalitySummary import generate_ts_image
        import os
        fig = plt.figure()
        
        if cutoff_c < 0:
            cutoff = _criticalPartialBrownianBridge(0.05/float(len(pdist)))
        else:
            cutoff = cutoff_c
        
        dn = dendrogram(average(pdist),labels=[re.findall('Position\s(\d+)',fn)[0] for fn in files], color_threshold=cutoff )
        fclusters = fcluster(average(pdist),cutoff,criterion='distance')
        nc = np.max(fclusters)
        fclusters_files = [[files[i] for i in range(len(files)) if fclusters[i]==c] for c in range(1,nc+1)]
        
        fileName = CheckReplication.pre+name + 'dendrogram-all-files-new.eps'
        savePath = os.path.join(out_folder,fileName)
        fig.savefig(savePath)
        
        dfs_by_c = [pd.concat([df.loc[fn].set_index(['fn','index'],drop=False) for fn in cfiles]) for cfiles in fclusters_files]
        tss_by_c = [list(np.unique(np.sort(np.concatenate([tss_by_pos[fn] for fn in cfiles])))) for cfiles in fclusters_files]
        for fn in np.unique(profile.index.get_level_values(0)):
            profile.loc[fn,'fn']=fn
        profile.loc[:,'index'] = profile.index.get_level_values(1)
        profiles_by_c = [pd.concat([profile.loc[fn].set_index(['fn','index'],drop=True) for fn in cfiles]) for cfiles in fclusters_files]
        
        
        fileName = CheckReplication.pre+name + "discarded-files.txt"
        savePath = os.path.join(out_folder,fileName)
        fo_discard = open(savePath,'w',0)
        output_clusters = []
        ln = ()
        for i in range(len(fclusters_files)):
            if len(fclusters_files[i]) < CheckReplication.min_files1:
                print >> fo_discard, "Discard orphan positions "+str(fclusters_files[i])
                ln = ln + tuple(fclusters_files[i])
            elif len(dfs_by_c[i]) < CheckReplication.min_cells2:
                print >> fo_discard,"Discard clusters with small sample size "+str(fclusters_files[i])
                ln = ln + tuple(fclusters_files[i])
            else:
                output_clusters.append((i,len(dfs_by_c[i]), len(fclusters_files[i]), fclusters_files[i]))
        
        if len(ln) > 0:
            fileName = CheckReplication.pre+name + "outlier_tabs.txt"
            savePath = os.path.join(out_folder,fileName)
            fl = open(savePath,'w')
            fl.write('\n'.join(ln))
            fl.close()

        output_clusters = sorted(output_clusters,key=lambda tpl:tpl[1],reverse=True) 
        
        for i in range(len(output_clusters)):
            (c, n, nf, fns) = output_clusters[i]
            import pickle
            fileName = CheckReplication.pre+name + 'mortality_dataframe_cluster'+str(i)+'.pickle'
            savePath = os.path.join(out_folder,fileName)
            print savePath
            print len(dfs_by_c),'--c--',c
            pickle.dump(dfs_by_c[c],open(savePath,'w'))
            
            fileName = CheckReplication.pre+name + 'relative_timegrid_cluster'+str(i)+'.pickle'
            savePath = os.path.join(out_folder,fileName)
            print savePath
            pickle.dump(tss_by_c[c],open(savePath,'w'))
            
            fileName = CheckReplication.pre+name + 'timeseries_dataframe_cluster'+str(i)+'.pickle'
            savePath = os.path.join(out_folder,fileName)
            print savePath
            pickle.dump(profiles_by_c[c],open(savePath,'w'))
            
            fileName = CheckReplication.pre+name + 'cluster'+str(i)+'_tabs.txt'
            savePath = os.path.join(out_folder,fileName)
            print savePath
            f_s = open(savePath,'w')
            f_s.write('\n'.join(fns))
            f_s.close()
            
            #pickle.dump(dfs_by_c[c],open(os.path.join(strain_pth, exp_prefix+' mortality_dataframe_cluster'+str(i)+'.pickle'),'w'))
            #pickle.dump(tss_by_c[c],open(os.path.join(strain_pth, exp_prefix+' relative_timegrid_cluster'+str(i)+'.pickle'),'w'))
            #pickle.dump(profiles_by_c[c],
            #            open(os.path.join(strain_pth, exp_prefix+' timeseries_dataframe_cluster'+str(i)+'.pickle'),'w'))
            #fl = open(os.path.join(strain_pth, exp_prefix+' cluster'+str(i)+'_tabs.txt'),'w')
            #fl.write('\n'.join(fns))
            #fl.close()

            image = generate_ts_image(profiles_by_c[c], dfs_by_c[c])
            fig = plt.figure()
            plt.imshow(image, aspect=100/float(len(dfs_by_c[c])), interpolation='None')
            plt.title(name+' cluster'+str(i))
            plt.xlabel('Frame #')
            plt.ylabel('Cell #')
            fig.set_size_inches(10,10)
            plt.tight_layout()
            fileName = CheckReplication.pre+name+'colomap_all_files_cluster'+str(i)+'_new.png'
            savePath = os.path.join(out_folder, fileName)
            fig.savefig(savePath)
            #fig.savefig(os.path.join(strain_pth, exp_prefix+' colomap_all_files_cluster'+str(i)+'.png'))

        return dfs_by_c, tss_by_c,profiles_by_c,output_clusters
    
    @staticmethod
    def getSignificanceFromList(files_s, dfs_by_pos_s,tss_by_pos_s):
        out_l_sig = []
        out_l_pdist = []
        out_l_sig_pd = []
        for idx, val in enumerate(files_s):
            files = files_s[idx]
            dfs_by_pos = dfs_by_pos_s[idx]
            tss_by_pos = tss_by_pos_s[idx]
            significances, pdist, significances_pd = CheckReplication.getSignificance(files,dfs_by_pos,tss_by_pos)
            out_l_sig.append(significances)
            out_l_pdist.append(pdist)
            out_l_sig_pd.append(significances_pd)
        return out_l_sig, out_l_pdist, out_l_sig_pd
    
    @staticmethod
    def loadMeta_fromFiles(f_list):
        import pickle
        import datetime
        out = []
        
        for l in f_list:
            data_list=[]
            key_set = None
            for p in l:
                data = pickle.load(open(p,'rb'))
                data_list.append(data)
                if key_set is None:
                    key_set = set(data)
                else:
                    key_set = key_set & set(data)
            
            c = dict()
            for a in key_set:
                print "check dataset type: ",a, type(data_list[0][a])
                if type(data_list[0][a]) is list:
                    c[a] = data_list[0][a]
                    for d in data_list:
                        if c[a] == d[a]:
                            continue
                        c[a] = c[a] + d[a]
                if type(data_list[0][a]) is dict:
                    c[a] = data_list[0][a].items()
                    for d in data_list:
                        if c[a] == d[a].items():
                            continue
                        c[a] = c[a] + d[a].items()                       
                if type(data_list[0][a]) is datetime:
                    c[a] = min(data_list[0][a],data_list[1][a])       
                    for d in data_list:
                        c[a] = min(c[a],d[a])
            out.append(c)
        return out
            
            
    @staticmethod
    def getTss_by_pos(df_list,meta_list):
        import pickle
        import pandas as pd
        import numpy as np
        out_l = []
        for idx_l_df, val_l_df in enumerate(df_list):
            tss_by_pos = None
            for idx_df_p,val_df_p in enumerate(val_l_df):
                print 'path df ',val_df_p
                print 'path meta ', meta_list[idx_l_df][idx_df_p]
                df = pickle.load(open(val_df_p,'rb'))
                files = list(np.unique(df.index.get_level_values(0)))
                meta = pickle.load(open(meta_list[idx_l_df][idx_df_p]))
                tss = {fn: [(abst - meta['time_grid'][0]).total_seconds()/3600 for abst in meta['time_indeces'][fn]] for fn in files}
                if tss_by_pos == None:
                    tss_by_pos = tss.items()
                else:
                    tss_by_pos = tss_by_pos+tss.items()
            tss_by_pos = dict(tss_by_pos)
            out_l.append(tss_by_pos)
        return out_l
    

        
    @staticmethod
    def fillNone(profile):
        # Non interpolation method, only fill the NaN by the previous value
        for row in profile.index:
            profile.loc[row].fillna(method = "ffill",inplace=True)
            profile.loc[row].fillna(method = "bfill",inplace=True)
    
    @staticmethod
    def loadDf_timeseries_fromFiles(f_list,df):
        import pandas as pd
        import pickle
        out = []
        out_files = []
        out_dfs_by_pos = []
        for l in f_list:
            data_list=[]
            for p in l:
                data = pickle.load(open(p,'rb'))
                data_list.append(data)
            c = pd.concat(data_list)            
            if df:
                import numpy as np
                files = list(np.unique(c.index.get_level_values(0)))
                dfs_by_pos = {fn:c.loc[fn] for fn in files}
                out_files.append(files)
                out_dfs_by_pos.append(dfs_by_pos)
            else:# means time_series or profile
                CheckReplication.fillNone(c)
            out.append(c)

        return out, out_files, out_dfs_by_pos
            
                
        
    @staticmethod
    def getRefDf_timegrid(path_df='', path_timegrid='',path_folder = ''):
        import pickle
        if len(path_df) == 0 and len(path_folder):
            path_df_list =CheckReplication.listFilePathsByKey(path_folder,CheckReplication.ext_ref_df)
            path_df = path_df_list[0][0]
        if len(path_timegrid)==0 and len(path_folder):
            path_timegrid_list = CheckReplication.listFilePathsByKey(path_folder, CheckReplication.ext_ref_timegrid)
            path_timegrid = path_timegrid_list[0][0]
        if len(path_folder):
            print "path ref_df: ",path_df
            print "path ref_timegrid: ", path_timegrid
            
        df_ref = pickle.load(open(path_df,'rb'))
        timegrid_ref= pickle.load(open(path_timegrid,'rb'))
        return df_ref,timegrid_ref
        
    
    @staticmethod
    def listFilePathsByExt(folderPath, ext=''):
        """
        -#-#-#
        Return files ended by ext, in sub folder list 
        :param str folderPath: folder path
        :return list files' absolute paths by (sub)folder
        """
        import path
        paths=[];
        import os
        if not os.path.isdir(folderPath):
            print "folder not existed, ",folderPath
            return paths
        
        #remove the last non alphanumeric char if user input / at the end
        if folderPath[-1]=='/':
            folderPath = folderPath[:-1]
        
        try:
            all_subfolders = next(os.walk(folderPath))[1]
        except:
            all_subfolders = []

        p=[]
        for fname in os.listdir(folderPath):
            if fname.endswith(ext):
                p.append(os.path.join(folderPath,fname))
        p.sort()
        paths.append(p)
        
        if len(all_subfolders):
            for sub in all_subfolders:
                p=[]
                subPath = os.path.join(folderPath,sub)
                for fname in os.listdir(subPath):
                    if fname.endswith(ext):
                        p.append(os.path.join(subPath,fname))
                p.sort()
                paths.append(p)
        paths_ = [x for x in paths if x]
                    
        return paths_    
    
    @staticmethod
    def listFilePathsByKey(folderPath, key=''):
        """
        -#-#-#
        Return files ended by ext, in sub folder list 
        :param str folderPath: folder path
        :return list files' absolute paths by (sub)folder
        """
        import path
        paths=[];
        import os
        if not os.path.isdir(folderPath):
            print "folder not existed"
            return paths
        
        #remove the last non alphanumeric char if user input / at the end
        if not folderPath[-1].isalnum():
            folderPath = folderPath[:-1]
            
        all_subfolders = next(os.walk(folderPath))[1]

        p=[]
        for fname in os.listdir(folderPath):
            if key in fname:
                p.append(os.path.join(folderPath,fname))
        paths.append(p)
        
        if len(all_subfolders):
            for sub in all_subfolders:
                p=[]
                subPath = os.path.join(folderPath,sub)
                for fname in os.listdir(subPath):
                    if key in fname:
                        p.append(os.path.join(subPath,fname))
                paths.append(p)
        paths_ = [x for x in paths if x]
                    
        return paths_
    

