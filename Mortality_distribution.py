
# coding: utf-8

# In[ ]:

class MortalityDistribution(object):
    
    ext_df_timegrid = 'df_relative_timegrid.pickle'
    ext_timeseries='timeseries_dataframe.pickle'
    
    ext_figPath = '_Mortality_summary & binning.png'
    ext_figPathM = '_Mortality_summary_binning poster.png'
    ext_figPathFitting = '_KM_NA_BSH_GGMfit.png'
    ext_figPathFittingLogScale='_KM_NA_BSH_GGMfit_LogScale.png'
    ext_ggmPath ='_GG_GGM_parameters.txt'
    pre='MD-'
    
    def __init__(self,folder_path='',ref_path = '',showProgress=True):
        """
        :param folder_path: contains df and timeseries pickles
        :param ref_path: reference file path of df_timegrid.pickle
        """
        if len(folder_path) > 0:
            self.folder_path = folder_path
        else:
            self.folder_path = MortalityDistribution.getCurrentPath()
        self.showProgress=showProgress
        self.ref_path = ref_path
        
    def run(self):
        import os
        files_df_timegrid = MortalityDistribution.listFilePathsByFolder(self.folder_path, ext = MortalityDistribution.ext_df_timegrid)
        files_timeseries= MortalityDistribution.listFilePathsByFolder(self.folder_path, ext = MortalityDistribution.ext_timeseries)
        if len(files_df_timegrid)==0 or len(files_timeseries)==0:
            print "No relative_timegrid.pickle or timeseries_dataframe.pickle found"
            return           
        self.comparisonList_df_timegrid = MortalityDistribution.getCombinationFromList(files_df_timegrid)
        if self.showProgress:
            print files_df_timegrid
            print self.comparisonList_df_timegrid
            print files_timeseries
        self.comparisonList_timeseries = MortalityDistribution.getCombinationFromList(files_timeseries)           
        #self.prefixName = MortalityDistribution.getPrefixSaveName(self.comparisonList_df_timegrid,self.showProgress)
        self.pairDf_tsPath = MortalityDistribution.getPairDF_TS_Path(self.comparisonList_df_timegrid, self.comparisonList_timeseries,showProgress=False)
        self.runCombinationPath(self.pairDf_tsPath,self.ref_path,files_df_timegrid, files_timeseries,self.showProgress)
        


    def runCombinationPath(self,pair_path_list,path_ref_df,files_df_timegrid, files_timeseries,showProgress=True):
        df_dict_ref=dict()
        timegrid_dict_ref=dict()
        lengths_df_list_ref=[]
        keys_list_ref=[]
        ref_loaded=False
        ref_Name=''
        if len(path_ref_df):
            ref_Name,xlimit_ref= MortalityDistribution.loadRefPickle(path_ref_df,df_dict_ref,timegrid_dict_ref,lengths_df_list_ref,keys_list_ref)
            if showProgress: 
                print 'load reference pickle, ', ref_Name, xlimit_ref

        import copy,os
        if len(pair_path_list)==0 and len(files_df_timegrid) and len(files_timeseries):
            for ind, val in enumerate(files_df_timegrid):
                df_1_path = val[0]
                ts_1_path = files_timeseries[ind][0]
                output_folder_path = os.path.split(df_1_path)[0]
                print "files will be written into : ",output_folder_path
                df_1,timegrid_1,ts1,KM_1, NA_1, BSFit_1, Image_1= MortalityDistribution.getStats(df_1_path,ts_1_path)
                xlimit = NA_1.index[-1]
                xlimit = int(xlimit+0.1*xlimit)
                pltdata_1 = MortalityDistribution.preparePlotData(df_1, timegrid_1,NA_1)
                biName = MortalityDistribution.getCombinedName(df_1_path,ts_1_path)
                print "biName is : ", biName
                
                if len(ref_Name):
                    df_dict = copy.deepcopy(df_dict_ref)
                    keys_list = copy.deepcopy(keys_list_ref)
                    timegrid_dict = copy.deepcopy(timegrid_dict_ref)
                    lengths_df_list = copy.deepcopy(lengths_df_list_ref)
                    
                    xlimit = max(xlimit_ref,xlimit) 
                    base = os.path.split(df_1_path)[1]
                    base.replace(MortalityDistribution.ext_df_timegrid,'')

                    df_dict[base]=df_1
                    timegrid_dict[base]=timegrid_1                
                    keys_list.append(base)
                    lengths_df_list.append(str(len(df_1)))


                    GGMfit_dict = MortalityDistribution.getGGMfit_dict(df_dict)
                    KM_dict,NA_dict,BSFit_dict=MortalityDistribution.getKMNABSfits_dict(df_dict,timegrid_dict)

                    figSavePath=os.path.join(output_folder_path,MortalityDistribution.pre+biName+'-vs-'+ref_Name+MortalityDistribution.ext_figPathFitting)
                    MortalityDistribution.runFitting(timegrid_dict,GGMfit_dict,KM_dict, NA_dict, BSFit_dict,figSavePath,xlimit,keys_list,lengths_df_list)

                    savePath = os.path.join(output_folder_path,MortalityDistribution.pre+biName+'-vs-'+ref_Name+MortalityDistribution.ext_ggmPath)
                    MortalityDistribution.saveGGMFit(GGMfit_dict,KM_dict,savePath)

                    figSavePath=os.path.join(output_folder_path,MortalityDistribution.pre+biName+'-vs-'+ref_Name+MortalityDistribution.ext_figPathFittingLogScale)
                    MortalityDistribution.plotLogScaleGGMFit(timegrid_dict,GGMfit_dict,KM_dict, NA_dict, BSFit_dict,figSavePath,xlimit,keys_list)

                    self.GGMfit_dict = GGMfit_dict
                    self.KM_dict = KM_dict
                    self.timegrid_dict = timegrid_dict
                    self.NA_dict = NA_dict
                    self.BSFit_dict = BSFit_dict
                    self.lengths_df_list = lengths_df_list
                    self.xlimit = xlimit
                    self.keys_list = keys_list
                print "working directory is : ",os.getcwd()
        
        
        for l in pair_path_list:
            if showProgress:
                print "--------PAIR-------"
            df_1_path = l[0][0]
            if showProgress:
                print 'df1 path : '+ df_1_path
            df_2_path = l[0][1]
            if showProgress:
                print 'df2 path : '+ df_2_path
            ts_1_path = l[1][0]
            if showProgress:
                print 'ts1 path : '+ ts_1_path
            ts_2_path = l[1][1]
            if showProgress:
                print 'ts2 path : '+ ts_2_path
            
            output_folder_path = os.path.split(df_1_path)[0]
            print "files will be written into : ",output_folder_path
            
            df_1,timegrid_1,ts1,KM_1, NA_1, BSFit_1, Image_1= MortalityDistribution.getStats(df_1_path,ts_1_path)
            df_2,timegrid_2,ts2,KM_2, NA_2, BSFit_2, Image_2= MortalityDistribution.getStats(df_2_path,ts_2_path)
            
            xlimit = max(NA_1.index[-1],NA_2.index[-1])
            xlimit = int(xlimit+0.1*xlimit)
            
            pltdata_1 = MortalityDistribution.preparePlotData(df_1, timegrid_1,NA_1)
            pltdata_2 = MortalityDistribution.preparePlotData(df_2, timegrid_2,NA_2)
            
            biName = MortalityDistribution.getCombinedName(df_1_path,df_2_path)
            
            figSavePath=os.path.join(output_folder_path,MortalityDistribution.pre+biName+MortalityDistribution.ext_figPath)
            MortalityDistribution.plotData(df_1,pltdata_1.copy(), pltdata_2.copy(),KM_1,KM_2,NA_1,NA_2,xlimit,figSavePath)
            
            figSavePath=os.path.join(output_folder_path,MortalityDistribution.pre+biName+MortalityDistribution.ext_figPathM)
            MortalityDistribution.plotReplicate(df_1,NA_1,NA_2,Image_1,pltdata_1.copy(),pltdata_2.copy(),xlimit,figSavePath)
            #update xlimit by NA Ref
            if len(ref_Name):
                df_dict = copy.deepcopy(df_dict_ref)
                keys_list = copy.deepcopy(keys_list_ref)
                timegrid_dict = copy.deepcopy(timegrid_dict_ref)
                lengths_df_list = copy.deepcopy(lengths_df_list_ref)
                
                xlimit = max(xlimit_ref,xlimit) 
                base = os.path.split(df_1_path)[1]
                base.replace(MortalityDistribution.ext_df_timegrid,'')

                
                df_dict[base]=df_1
                timegrid_dict[base]=timegrid_1                
                keys_list.append(base)
                lengths_df_list.append(str(len(df_1)))
                
                base = os.path.split(df_2_path)[1]
                base.replace(MortalityDistribution.ext_df_timegrid,'')
                df_dict[base]=df_2
                timegrid_dict[base]=timegrid_2
                keys_list.append(base)
                lengths_df_list.append(str(len(df_2)))
                
                GGMfit_dict = MortalityDistribution.getGGMfit_dict(df_dict)
                KM_dict,NA_dict,BSFit_dict=MortalityDistribution.getKMNABSfits_dict(df_dict,timegrid_dict)
                
                figSavePath=os.path.join(output_folder_path,MortalityDistribution.pre+biName+'-vs-'+ref_Name+MortalityDistribution.ext_figPathFitting)
                try:
                    MortalityDistribution.runFitting(timegrid_dict,GGMfit_dict,KM_dict, NA_dict, BSFit_dict,figSavePath,xlimit,keys_list,lengths_df_list)
                    savePath = os.path.join(output_folder_path,MortalityDistribution.pre+biName+'-vs-'+ref_Name+MortalityDistribution.ext_ggmPath)
                    MortalityDistribution.saveGGMFit(GGMfit_dict,KM_dict,savePath)
                    figSavePath=os.path.join(output_folder_path,MortalityDistribution.pre+biName+'-vs-'+ref_Name+MortalityDistribution.ext_figPathFittingLogScale)
                    MortalityDistribution.plotLogScaleGGMFit(timegrid_dict,GGMfit_dict,KM_dict, NA_dict, BSFit_dict,figSavePath,xlimit,keys_list)

                    self.GGMfit_dict = GGMfit_dict
                    self.KM_dict = KM_dict
                    self.timegrid_dict = timegrid_dict
                    self.NA_dict = NA_dict
                    self.BSFit_dict = BSFit_dict
                    self.lengths_df_list = lengths_df_list
                    self.xlimit = xlimit
                    self.keys_list = keys_list
                except:
                    print('error MortalityDistribution.runFitting')
                

            print "working directory is : ",os.getcwd()
    
    @staticmethod 
    def loadRefPickle(fPath, df_dict, timegrid_dict,lengths_df_list,keys_list):
        import sys,os
        sys.path.insert(1, os.path.join(sys.path[0], '..'))
        sys.path.insert(1, os.path.join(sys.path[0], '.'))
        from cohortMortalitySummary import KaplanMeier, NelsonAalen, Breslow, BSHazardR, generate_ts_image
        
        #try:
        df_timegrid=MortalityDistribution.loadPickle(fPath)
        df_ref = df_timegrid.df
        timegrid_ref = df_timegrid.time_grid_relative

        fName_base = os.path.split(fPath)[1]
        
        #print fName_base
        #fName_base.replace(MortalityDistribution.ext_df_timegrid,'')

        if not fName_base[-1].isalnum():
            fName_base = fName_base[:-1]    

        keys_list.append(fName_base);
        df_dict[fName_base] = df_ref
        timegrid_dict[fName_base] = timegrid_ref
        lengths_df_list.append(str(len(df_ref)));
        NARef = NelsonAalen(df_ref, timegrid_ref)
        if '.' in fName_base:
            separator_index = fName_base.index('.')
            fName_base = fName_base[:separator_index]
        print fName_base
        return fName_base,NARef.index[-1]
        #except:
            #print 'ERROR : Reference df pickle'
            #return ''
        
    @staticmethod
    def getGGMfit_dict(df_dict):
        import sys,os
        workingpath = os.getcwd()
        sys.path.insert(1, os.path.join(sys.path[0], '..'))
        sys.path.insert(1, os.path.join(sys.path[0], '.'))
        from cohortMortalitySummary import GGMfit
        print "getGGMfit_dict--(1)",os.getcwd()
        GGMfit_dicts = dict()
        
        try:
            for exp in df_dict.keys(): 
                GGMfit_dicts[exp] = GGMfit(df_dict[exp])
                print exp,GGMfit_dicts[exp].keys()
            os.chdir(workingpath)
            print "getGGMfit_dict--(2)",os.getcwd()
            return GGMfit_dicts
        except:
            return -1
    
    @staticmethod
    def getKMNABSfits_dict(df_dict, timegrid_dict):
        import sys,os
        print "getKMNABSfits_dict--",os.getcwd()
        sys.path.insert(1, os.path.join(sys.path[0], '..'))
        sys.path.insert(1, os.path.join(sys.path[0], '.'))
        from cohortMortalitySummary import KaplanMeier, NelsonAalen, BSHazardR
        KMs = dict()
        NAs = dict()
        BSFits = dict()
        for exp in df_dict.keys():
            print exp
            KMs[exp] = KaplanMeier(df_dict[exp], timegrid_dict[exp])
            NAs[exp] = NelsonAalen(df_dict[exp], timegrid_dict[exp])
            BSFits[exp] = BSHazardR(df_dict[exp])
        print "getKMNABSfits_dict--",os.getcwd()
        return KMs, NAs, BSFits
    
    @staticmethod
    def runFitting(tss,GGMfit_dicts,KMs, NAs, BSFits,figSavePath,xlimit,keys_list,lengths_df_list):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_context("paper")
        sns.set_style("ticks")
        f, axes = plt.subplots(3, 1, sharex=False, sharey=False)
        fontsize1 = 16
        fontsize2 = 14
        #print zip(['wt'],[baseRef],range(1))
        base ='strain'
        for (strain,ttl,row) in zip([base],[base],range(1)):
            axes[1].set_yscale('log')
            axes[2].set_yscale('log')
            axes[0].set_title(r'Survivorship $S(t)$',fontsize=fontsize1)
            axes[1].set_title(r'Cumulative hazard rate $H(t)$',fontsize=fontsize1)
            axes[2].set_title(r'Smoothed hazard rate $\hat{h}(t)$',fontsize=fontsize1)
            
            #axes[0].set_ylabel(ttl, fontsize=fontsize1)
            #axes[0].set_ylabel(ttl, fontsize=fontsize1)
            #axes[0].set_ylabel(ttl, fontsize=fontsize1)

            for ind,exp in enumerate(keys_list):
                #exp = strain + str(repeat)
                print ind, exp

                clr = sns.color_palette()[ind];
                nc = lengths_df_list[ind];
                #axes[row][0].plot(KMs[exp].index,KMs[exp].loc[:,'survivorship'],linewidth=0.75,
                #                  label=ttl+str(repeat),color=clr)
                axes[0].fill_between(KMs[exp].index,KMs[exp].loc[:,'lower_ci'],KMs[exp].loc[:,'upper_ci'],
                                          alpha=0.5,label=exp,color=clr)
                axes[0].plot(tss[exp][1:],GGMfit_dicts[exp]['ML_survivorship'](tss[exp][1:]),linewidth=2, linestyle='--',
                                label=exp+'-'+nc,color=clr)

                #axes[row][1].plot(NAs[exp].index,NAs[exp].loc[:,'cumulative_hazard'],linewidth=0.75,
                #                  label=ttl+str(repeat),color=clr)
                axes[1].fill_between(NAs[exp].index,NAs[exp].loc[:,'lower_ci'],NAs[exp].loc[:,'upper_ci'],
                                          alpha=0.5,label=exp,color=clr)
                axes[1].plot(tss[exp][1:],GGMfit_dicts[exp]['ML_cumulative_hazard'](tss[exp][1:]),linewidth=2, linestyle='--',
                                label=exp,color=clr)

                #axes[row][2].plot(BSFits[exp]['time'],BSFits[exp]['hazard'],linewidth=0.75,
                #                  label=ttl+str(repeat),color=clr)
                axes[2].fill_between(BSFits[exp]['time'],BSFits[exp]['lower.ci'],BSFits[exp]['upper.ci'],
                                          alpha=0.5,label=exp,color=clr)
                axes[2].plot(tss[exp][1:],GGMfit_dicts[exp]['ML_hazard'](tss[exp][1:]),linewidth=1.5, linestyle='--',
                                label=exp,color=clr)

            axes[0].legend(loc=1,fontsize=fontsize2)
            axes[1].legend(loc=4,fontsize=fontsize2)
            axes[2].legend(loc=4,fontsize=fontsize2)
            axes[0].set_xlim([0,xlimit])
            axes[1].set_xlim([0,xlimit])
            axes[2].set_xlim([0,xlimit])
            axes[0].set_xlabel('Time (h)',fontsize=fontsize1)
            axes[1].set_xlabel('Time (h)',fontsize=fontsize1)
            axes[2].set_xlabel('Time (h)',fontsize=fontsize1)
            axes[0].set_ylim(0,1)
            axes[1].set_ylim(0.003,10)
            axes[2].set_ylim(0.0002,0.5)

        f.set_size_inches(16,30)
        sns.despine(fig=f)
        plt.tight_layout()
        f.savefig(figSavePath)
        
    @staticmethod
    def getlifespan_dict(KMs):
        import pandas as pd
        lifespan = dict();
        for exp in sorted(KMs.keys()):
            print exp
            upper=KMs[exp].loc[:,'upper_ci']
            lower = KMs[exp].loc[:,'lower_ci']
            var =  KMs[exp].loc[:,'variance']
            survivor=KMs[exp].loc[:,'survivorship']
            s = survivor.values >= 0.5    
            slower = lower.values >= 0.5
            supper = upper.values >= 0.5
            ind_survivor = survivor.iloc[s].index[-1];  
            indlower = lower.iloc[slower].index[-1]
            indupper = upper.iloc[supper].index[-1]
            data = [{'est': ind_survivor, 'L95%': indlower, 'U95%':indupper, 'var':var.loc[ind_survivor]}]
            lifespan[exp] = pd.DataFrame(data,index=['lifespan'],columns=['est','L95%','U95%','var'])
            #print lifespan[exp]
        return lifespan
            
    @staticmethod
    def saveGGMFit(GGMfit_dict,KMs,savePath,showProgress= True):
        lifespan_dict = MortalityDistribution.getlifespan_dict(KMs)        
        if showProgress:
            print savePath
            for exp in sorted(GGMfit_dict.keys()):
                print exp
                print '-------------------------------\nGGM(AIC={:1f})'.format(GGMfit_dict[exp]['AIC'])
                print GGMfit_dict[exp]['model_paras']
                print lifespan_dict[exp]
                print '\n'
        
        fo = open(savePath,'w',0)
        print savePath
        for exp in sorted(GGMfit_dict.keys()):
            print >> fo, exp
            print >> fo, '-------------------------------\nGGM(AIC={:1f})'.format(GGMfit_dict[exp]['AIC'])
            print >> fo, GGMfit_dict[exp]['model_paras']
            print >> fo, lifespan_dict[exp]
            print >> fo, '\n'
            
    @staticmethod
    def plotLogScaleGGMFit(tss,GGMfit_dicts,KMs, NAs, BSFits,figSavePath,xlimit,keys_list):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_context("paper")
        sns.set_style("ticks")
        f, axes = plt.subplots(3, 1, sharex=False, sharey=False)
        fontsize1 = 16
        fontsize2 = 14
        #print zip(['wt'],[baseRef],range(1))
        base = 'strain'
        for (strain,ttl,row) in zip([base],[base],range(1)):
            axes[1].set_yscale('log')
            axes[2].set_yscale('log')
            axes[0].set_title(r'Survivorship $S(t)$',fontsize=fontsize1)
            axes[1].set_title(r'Cumulative hazard rate $H(t)$',fontsize=fontsize1)
            axes[2].set_title(r'Smoothed hazard rate $\hat{h}(t)$',fontsize=fontsize1)
            #axes[0].set_ylabel(ttl, fontsize=fontsize1)
            #axes[0].set_ylabel(ttl, fontsize=fontsize1)
            #axes[0].set_ylabel(ttl, fontsize=fontsize1)

            for ind,exp in enumerate(keys_list):
                #exp = strain + str(repeat)
                print ind, exp

                clr = sns.color_palette()[ind];
                #axes[row][0].plot(KMs[exp].index,KMs[exp].loc[:,'survivorship'],linewidth=0.75,
                #                  label=ttl+str(repeat),color=clr)
                axes[0].fill_between(KMs[exp].index,KMs[exp].loc[:,'lower_ci'],KMs[exp].loc[:,'upper_ci'],
                                          alpha=0.5,label=exp,color=clr)
                axes[0].plot(tss[exp][1:],GGMfit_dicts[exp]['ML_survivorship'](tss[exp][1:]),linewidth=2, linestyle='--',
                                label=exp,color=clr)

                #axes[row][1].plot(NAs[exp].index,NAs[exp].loc[:,'cumulative_hazard'],linewidth=0.75,
                #                  label=ttl+str(repeat),color=clr)
                axes[1].fill_between(NAs[exp].index,NAs[exp].loc[:,'lower_ci'],NAs[exp].loc[:,'upper_ci'],
                                          alpha=0.5,label=exp,color=clr)
                axes[1].plot(tss[exp][1:],GGMfit_dicts[exp]['ML_cumulative_hazard'](tss[exp][1:]),linewidth=2, linestyle='--',
                                label=exp,color=clr)

                #axes[row][2].plot(BSFits[exp]['time'],BSFits[exp]['hazard'],linewidth=0.75,
                #                  label=ttl+str(repeat),color=clr)
                axes[2].fill_between(BSFits[exp]['time'],BSFits[exp]['lower.ci'],BSFits[exp]['upper.ci'],
                                          alpha=0.5,label=exp,color=clr)
                axes[2].plot(tss[exp][1:],GGMfit_dicts[exp]['ML_hazard'](tss[exp][1:]),linewidth=1.5, linestyle='--',
                                label=exp,color=clr)

            axes[0].legend(loc=1,fontsize=fontsize2)
            axes[1].legend(loc=4,fontsize=fontsize2)
            axes[2].legend(loc=4,fontsize=fontsize2)
            axes[0].set_xlim([0,xlimit])
            axes[1].set_xscale('log')
            axes[2].set_xscale('log')
            axes[1].set_xlim([10,xlimit])
            axes[2].set_xlim([10,xlimit])
            axes[0].set_xlabel('Time (h)',fontsize=fontsize1)
            axes[1].set_xlabel('Time (h)',fontsize=fontsize1)
            axes[2].set_xlabel('Time (h)',fontsize=fontsize1)
            axes[0].set_ylim(0,1)
            axes[1].set_ylim(0.003,10)
            axes[2].set_ylim(0.0002,0.5)
        f.set_size_inches(16,30)
        sns.despine(fig=f)
        plt.tight_layout()
        f.savefig(figSavePath)
            
          
    
    @staticmethod
    def plotReplicate(df1,NA1,NA2,image1,pltdata1,pltdata2,xlimit,figSavePath):
        import seaborn as sns
        sns.set_context("paper")
        sns.set_style("ticks")
        import matplotlib.gridspec as gs
        import matplotlib.patches as patches
        import matplotlib.pyplot as plt
        import numpy as np


        dt = 2 #Hours
        figw = 12 #Inches
        figh = 24 #Inches
        gray1 = 13 #Hours
        gray2 = 93 #Hours

        f = plt.gcf()
        gs = gs.GridSpec(18, 2)
        ax1 = plt.subplot(gs[:8, :])
        image_asp = 0.47
        ax2 = plt.subplot(gs[8:13, :])
        ax3 = plt.subplot(gs[13:18, :])
        #f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
        ax1.imshow(image1, aspect=image_asp*figh/float(figw)*float(image1.shape[1])/float(image1.shape[0]), interpolation='None')
        ax1.set_ylabel('Cell #',fontsize=18)
        ax1.legend(loc=1,fontsize=12)
        ax1.set_ylim([image1.shape[0],0])
        ax2.plot(NA1.index,NA1.loc[:,'cumulative_hazard'],color='b',linewidth=0.75,label='Nelson-Aalen, exp1')
        ax2.fill_between(NA1.index,NA1.loc[:,'lower_ci'],NA1.loc[:,'upper_ci'],color=sns.color_palette()[0],alpha=0.5,label='95% CI, exp1')
        ax2.plot(NA2.index+dt,NA2.loc[:,'cumulative_hazard'],color='r',linewidth=0.75,label='Nelson-Aalen, exp2')
        ax2.fill_between(NA2.index+dt,NA2.loc[:,'lower_ci'],NA2.loc[:,'upper_ci'],color=sns.color_palette()[2],alpha=0.5,label='95% CI, exp2')
        ax2.legend(loc=4,fontsize=14)
        ax2.set_ylabel('Cumulative hazard rate H(t)',fontsize=18)
        ax2.set_yscale('log')
        ax2.set_ylim(-np.log((len(df1)-10)/float(len(df1))),-np.log(0.5/float(len(df1))))

        ax2s = ax2.twinx()
        ax2s.set_ylabel('Survivorship S(t)',fontsize=18,labelpad=15)
        ax2s.tick_params('y', colors='k')
        chf_lim = (-np.log((len(df1)-10)/float(len(df1))),-np.log(0.5/float(len(df1))))
        survival_ticks = list( np.arange(0,1,0.1) )
        survival_labels = ['%.2f'%(1/np.exp(s)) 
                           for s in np.exp(np.arange(np.log(chf_lim[0]),
                                                     np.log(chf_lim[1]),
                                                     0.1*np.log(chf_lim[1]/chf_lim[0]) ) 
                                          )]

        ax2s.set_yticks(survival_ticks)
        ax2s.set_yticklabels(survival_labels)
        ax3.errorbar(0.5*(pltdata1['tL']+pltdata1['tR']),pltdata1['dNA'] ,
                     yerr=np.array(pltdata1[['hr_eb_l','hr_eb_u']].transpose()),
                     xerr=0.5*pltdata1['dt'],
                    fmt='none',ecolor='b', elinewidth=1.2, capsize=5, capthick=0.8, alpha=0.6,
                    label='Binning estimates, exp1')
        ax3.errorbar(0.5*(pltdata2['tL']+pltdata2['tR'])+dt,pltdata2['dNA'] ,
                     yerr=np.array(pltdata2[['hr_eb_l','hr_eb_u']].transpose()),
                     xerr=0.5*pltdata2['dt'],
                    fmt='none',ecolor='r', elinewidth=1.2, capsize=5, capthick=0.8, alpha=0.6,
                    label='Binning estimates, exp2')
        ax3.set_ylim(3/float(len(df1)),0.3)
        ax3.legend(loc=4,fontsize=14)
        ax3.set_ylabel('Hazard rate h(t)',fontsize=18)
        ax3.set_yscale('log')
        ax3.set_xlabel('Time (h)',fontsize=18)
        ax1.set_xlim([0,xlimit])
        ax2.set_xlim([0,xlimit])
        ax3.set_xlim([0,xlimit])

        f.set_size_inches(figw,figh)
        sns.set(font_scale = 4)
        sns.despine(ax=ax1,bottom=True,left=True)
        sns.despine(ax=ax2,top=True,right=False)
        sns.despine(ax=ax2s,top=True,right=False)
        ax1.tick_params(axis='both', labelsize=12)
        ax2.tick_params(axis='both', labelsize=12)
        ax2s.tick_params(axis='both', labelsize=12)
        ax3.tick_params(axis='both', labelsize=12)
        sns.despine(ax=ax3)
        plt.tight_layout()
        print figSavePath
        #print 'f',f
        plt.savefig(figSavePath,dpi=300)
    
    @staticmethod
    def plotData(df,pltdata1, pltdata2,KM1,KM2,NA1,NA2,xlimit,figSavePath):
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set_context("paper")
        sns.set_style("white")
        dt = 2
        f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
        ax1.plot(KM1.index,KM1.loc[:,'survivorship'],color='b',linewidth=0.75,label='Kaplan-Meier, exp1')
        ax1.fill_between(KM1.index,KM1.loc[:,'lower_ci'],KM1.loc[:,'upper_ci'],color=sns.color_palette()[0],alpha=0.5,label='95% CI, exp1')
        ax1.plot(KM2.index+dt,KM2.loc[:,'survivorship'],color='r',linewidth=0.75,label='Kaplan-Meier, exp2')
        ax1.fill_between(KM2.index+dt,KM2.loc[:,'lower_ci'],KM2.loc[:,'upper_ci'],color=sns.color_palette()[2],alpha=0.5,label='95% CI, exp2')
        ax1.set_ylabel('Survivorship S(t)',fontsize=12)
        ax1.legend(loc=1,fontsize=10)
        ax2.plot(NA1.index,NA1.loc[:,'cumulative_hazard'],color='b',linewidth=0.75,label='Nelson-Aalen, exp1')
        ax2.fill_between(NA1.index,NA1.loc[:,'lower_ci'],NA1.loc[:,'upper_ci'],color=sns.color_palette()[0],alpha=0.5,label='95% CI, exp1')
        ax2.plot(NA2.index+dt,NA2.loc[:,'cumulative_hazard'],color='r',linewidth=0.75,label='Nelson-Aalen, exp2')
        ax2.fill_between(NA2.index+dt,NA2.loc[:,'lower_ci'],NA2.loc[:,'upper_ci'],color=sns.color_palette()[2],alpha=0.5,label='95% CI, exp2')
        ax2.legend(loc=4,fontsize=10)
        ax2.set_ylabel('Cumulative hazard rate H(t)',fontsize=12)
        ax2.set_yscale('log')
        ax2.set_ylim(-np.log((len(df)-1)/float(len(df))),-np.log(1/float(len(df))))
        ax3.errorbar(0.5*(pltdata1['tL']+pltdata1['tR']),pltdata1['dNA'] ,
                     yerr=np.array(pltdata1[['hr_eb_l','hr_eb_u']].transpose()),
                     xerr=0.5*pltdata1['dt'],
                    fmt='none',ecolor='b', elinewidth=1.2, capsize=5, capthick=0.8, alpha=0.6,
                    label='Binning estimates, exp1')
        ax3.errorbar(0.5*(pltdata2['tL']+pltdata2['tR'])+dt,pltdata2['dNA'] ,
                     yerr=np.array(pltdata2[['hr_eb_l','hr_eb_u']].transpose()),
                     xerr=0.5*pltdata2['dt'],
                    fmt='none',ecolor='r', elinewidth=1.2, capsize=5, capthick=0.8, alpha=0.6,
                    label='Binning estimates, exp2')
        ax3.legend(loc=4,fontsize=10)
        ax3.set_ylabel('Hazard rate h(t)',fontsize=12)
        ax3.set_yscale('log')
        ax3.set_xlabel('Time (h)',fontsize=12)
        ax3.set_xlim([0,xlimit])
        f.set_size_inches(5,12)
        sns.despine(fig=f)
        plt.tight_layout()
        f.savefig(figSavePath,dpi=300)
        #f.savefig('Diagnosis plots/Mortality summary wt1 & binning.png',dpi=300)
        
    @staticmethod
    def preparePlotData(data,relative_time_grid,NA):
        from sklearn.cluster import KMeans
        import numpy as np
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import pandas as pd
        from dateutil.parser import parse
        import os
        kmc = KMeans(n_clusters = int(np.floor(np.max(relative_time_grid)))).fit( np.reshape(relative_time_grid,(len(relative_time_grid),1)) )
        gap_mids = [ ((relative_time_grid[i]+relative_time_grid[i+1])/2,
                      relative_time_grid[i],relative_time_grid[i+1],      
                      i ) 
                    for i in range(len(relative_time_grid)-1) if kmc.labels_[i] != kmc.labels_[i+1] ]

        mortalities = pd.Series(data.loc[data['status']==1,'mortality'].value_counts(), index=relative_time_grid).fillna(0)
        censored = pd.Series(data.loc[data['status']==2,'mortality'].value_counts(), index=relative_time_grid).fillna(0)
        counts = mortalities + censored
        n = len(counts)
        # At Risk individuals at t_i are those who are at risk just before t_i, thus including the mortality & censoring events at t_i
        atRisk = pd.Series(data=np.dot(np.triu(np.ones((n,n)),k=0), counts), index=relative_time_grid)

        NA_limit = 0
        points = []
        for i in range(1,len(gap_mids)):
            (t,tm,tp,ind) = gap_mids[i]
            if NA.loc[tm,'lower_ci']>NA_limit:
                NA_limit = NA.loc[tm,'upper_ci']
                points.append(
                ( t,ind,tm,tp,
                     NA.loc[tm,'cumulative_hazard'],
                     NA.loc[tm,'lower_ci'],
                     NA.loc[tm,'upper_ci'],
                    )
                )

        from statsmodels.stats.proportion import proportion_confint
        txy = []
        for i in range(1,len(points)):
            current_point = points[i]
            previous_point = points[i-1]
            d = np.sum(mortalities.loc[previous_point[3]:current_point[2]])
            Y = atRisk.loc[previous_point[3]]
            txy.append(
                    (previous_point[2],current_point[3],current_point[0]-previous_point[0],
                     0.5*(previous_point[4]+current_point[4]),
                     0.5*(previous_point[5]+current_point[5]),
                     0.5*(previous_point[6]+current_point[6]),
                     (current_point[4]-previous_point[4])/(current_point[0]-previous_point[0]),
                     Y,d,
                     proportion_confint(d,Y,method='jeffrey',alpha=0.05)[0],
                     proportion_confint(d,Y,method='jeffrey',alpha=0.05)[1]
                    )
                    ) 

        pltdata = pd.DataFrame(data=txy,columns = ['tL','tR','dt','NA','lciNA','uciNA','dNA','Y','d','lcid','ucid'])
        pltdata['hr_eb_l'] = (pltdata['d']/pltdata['Y'] - pltdata['lcid'])/pltdata['dt']
        pltdata['hr_eb_u'] = (pltdata['ucid'] - pltdata['d']/pltdata['Y'])/pltdata['dt']
        pltdata['na_eb_l'] = [pltdata.loc[i,'NA']-pltdata.loc[i,'lciNA'] for i in pltdata.index]
        pltdata['na_eb_r'] = [pltdata.loc[i,'uciNA']-pltdata.loc[i,'NA'] for i in pltdata.index]
        return pltdata

    @staticmethod
    def getPrefixSaveName(comparisonList,showProgress=True):
        import os
        prefix=[]
        for c in comparisonList:
            if len(c) != 2:#need a pair in the tuple
                continue
            head_0, tail_0 = os.path.split(c[0])
            head_1, tail_1 = os.path.split(c[1])
            common = MortalityDistribution.commonStringFromlast(tail_0,tail_1)
            filenew_name = tail_0.replace(common, '')+'-'+tail_1.replace(common,'')
            prefix.append(filenew_name)
            if showProgress:
                print 'save file prefixed name is : '+ filenew_name
        return prefix
            
    @staticmethod
    def getCombinedName(path1,path2,showProgress=True):
        """
        remove common part of path file name, return combined one
        """
        import os
        head_0, tail_0 = os.path.split(path1)
        head_1, tail_1 = os.path.split(path2)
        
        common = MortalityDistribution.commonStringFromlast(tail_0,tail_1)
        filenew_name = tail_0.replace(common, '')+'-vs-'+tail_1.replace(common,'')
        
        if showProgress:
            print 'save file prefixed name is : '+ filenew_name
        if not filenew_name[-1].isalnum():
            return filenew_name[:-1]
        else:
            return filenew_name
  
        
    @staticmethod          
    def commonStringFromlast(string1, string2):
        """
        Alert : non common utility,
        """
        
        len1, len2 = len(string1), len(string2)
        r1 = string1[::-1]
        r2 = string2[::-1]
        match=''
        for j in range(len(r2)):
            if r1[j]==r2[j]:
                match=match+r1[j]
            else:
                break
        return match[::-1]
                    
            
        
    @staticmethod
    def getCombinationFromList(filePathsByFolder):
        """
        call by user
        :param list filePathsByFolder: [['a.pickle','b.pickle','c.pickle'],['d.pickle','e.pickle']]
        return [('a.pickle','b.pickle'),('a.pickle','c.pickle'),('d.pickle','e.pickle')]
        """
        pairFiles=[]
        for fps in filePathsByFolder:
            l = MortalityDistribution.getCombination(fps)
            for m in l:
                pairFiles.append(m)               
        return pairFiles
        
    
    @staticmethod
    def getCurrentPath():
        import os
        return os.getcwd()
    
    @staticmethod    
    def loadPickle(fPath):
        """
        load Df_timegrid (TabPostAnazlyer.py and Df_timegrid.py should be found in the current or parent path) or timeseries pickle
        """
        import sys,os        
        sys.path.insert(1, os.path.join(sys.path[0], '..'))
        sys.path.insert(1, os.path.join(sys.path[0], '.'))
        from TabPostAnalyzer import Df_timegrid
        
        import pickle
        import pandas as pd
        f = open(fPath, 'rb')
        #p = pd.read_pickle(f)
        p = pickle.load(f)
        return p
    
    @staticmethod
    def getPairDF_TS_Path(comparisonList_df_timegrid,comparisonList_timeseries,showProgress = True):
        """
        return list[tuple(dict(df), dict(ts))] each tuple content is a dict, each dict1 is df, dict2 is ts
        """
        import os
        output = []
        i = 0
        for c in comparisonList_df_timegrid:
            df_list = [None]*2
            ts_list = [None]*2
            
            if showProgress:
                print '----------PAIR('+str(i)+')----------'
            if len(c)!=2:
                continue
            df1_path = c[0]
            df2_path = c[1]
            base = os.path.split(df1_path)[0]
            for c1 in comparisonList_timeseries:
                if len(c)!=2:
                    continue     
                for c3 in c1:
                    if df1_path.replace(MortalityDistribution.ext_df_timegrid,'') == c3.replace(MortalityDistribution.ext_timeseries,''):
                        ts1_path = c3
                    if df2_path.replace(MortalityDistribution.ext_df_timegrid,'') == c3.replace(MortalityDistribution.ext_timeseries,''):
                        ts2_path = c3   
            
            df_list[0]=df1_path
            ts_list[0]=ts1_path
            
            df_list[1]=df2_path
            ts_list[1]=ts2_path
            i = i+1
            output.append([df_list,ts_list])
            if showProgress:
                print os.path.split(df1_path)[1], os.path.split(ts1_path)[1]
                print "-VS-"
                print os.path.split(df2_path)[1], os.path.split(ts2_path)[1]
        return output
    
    @staticmethod
    def getStats(path_df_timegrid, path_timeseries):
        df_timegrid = MortalityDistribution.loadPickle(path_df_timegrid)
        ts = MortalityDistribution.loadPickle(path_timeseries)
        KM,NA,BSFit,Image=MortalityDistribution.getStatistic(df_timegrid,ts)
        return df_timegrid.df, df_timegrid.time_grid_relative,ts, KM, NA, BSFit, Image
    
    @staticmethod
    def getStatistic(df_timegrid, ts):        
        import sys,os
        sys.path.insert(1, os.path.join(sys.path[0], '..'))
        sys.path.insert(1, os.path.join(sys.path[0], '.'))
        from cohortMortalitySummary import KaplanMeier, NelsonAalen, Breslow, BSHazardR, generate_ts_image
        df = df_timegrid.df
        timegrid = df_timegrid.time_grid_relative
        
        KM = KaplanMeier(df, timegrid)
        NA = NelsonAalen(df, timegrid)
        BSFit = BSHazardR(df)
        Image = generate_ts_image(ts, df)
        return KM,NA,BSFit,Image

        
    @staticmethod
    def listFilePathsByFolder(folderPath, ext=''):
        """
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
        if folderPath[-1]=='/':
            folderPath = folderPath[:-1]
            
        all_subfolders = next(os.walk(folderPath))[1]

        p=[]
        for fname in os.listdir(folderPath):
            if fname.endswith(ext):
                p.append(os.path.join(folderPath,fname))
        paths.append(p)
        
        if len(all_subfolders):
            for sub in all_subfolders:
                p=[]
                subPath = os.path.join(folderPath,sub)
                for fname in os.listdir(subPath):
                    if fname.endswith(ext):
                        p.append(os.path.join(subPath,fname))
                paths.append(p)
        paths_ = [x for x in paths if x]
                    
        return paths_

    @staticmethod
    def getCombination(l, c=2):
        """
        get combination of list l
        """
        import itertools
        com = list(itertools.combinations(l,c))
        return com
    
    #@staticmethod
    #def 
        

