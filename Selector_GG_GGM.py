
# coding: utf-8

# In[ ]:

class GG_GGM_Selector(object):
    
    ext_cluster0_df = 'mortality_dataframe_cluster0.pickle'
    ext_cluster0_timegrid='relative_timegrid_cluster0.pickle'
    pre = 'SEL_'
    
    def __init__(self, folder_path='',out_folder='', showProgress = True):
        import os
        if len(folder_path) > 0:
            self.folder_path = folder_path
        else:
            self.folder_path = GG_GGM_Selector.getCurrentPath()
        self.showProgress=showProgress
        if len(out_folder):
            self.out_folder=out_folder
        else:
            self.out_folder = self.folder_path
        if len(self.out_folder) and not os.path.exists(self.out_folder):
            os.makedirs(self.out_folder)
            print 'create output folder : ',self.out_folder
        self.paths_df_timegrid = GG_GGM_Selector.listFilePathsByExt(self.folder_path, GG_GGM_Selector.ext_cluster0_df,GG_GGM_Selector.ext_cluster0_timegrid)
        if len(self.paths_df_timegrid) == 0:
            print "ERROR, no cluster0 files found!"
        if showProgress:
            print self.paths_df_timegrid
     
    
    def run(self):
        import pickle
        import os
        res = None
        for path_pair in self.paths_df_timegrid:
            if len(path_pair) !=3:
                print "lack of files"
                continue
            strainName = path_pair[2]
            print strainName
            if GG_GGM_Selector.ext_cluster0_df in path_pair[0]:
                df = pickle.load(open(path_pair[0],'rU'))
                timegrid = pickle.load(open(path_pair[1],'rU'))
            else:
                df = pickle.load(open(path_pair[1],'rU'))
                timegrid = pickle.load(open(path_pair[0],'rU'))
            res = GG_GGM_Selector.GG_GGM_To_Table(strainName, df,timegrid,res,writeIndividual=True,out_folder_path=self.out_folder)
        res.to_csv(os.path.join(self.out_folder,GG_GGM_Selector.pre+'GG_GGM_all.csv'))
            
    @staticmethod
    def run2():
        # save all to one table
        import csv
        import pickle

        strainParentPath = '../Microfluidics Database/'
        res = None;
        i = 1;
        with open('microfluidics database information_07102017_UPDATED.csv','rU') as csvfile:
            reader = csv.reader(csvfile, dialect=csv.excel_tab )
            for row in reader:
                tup_strain = findStrainLine(row[0]);
                if tup_strain[0]:
                    bestClusterPath = getClusterFilePath(strainParentPath, tup_strain[1], tup_strain[2])
                    print bestClusterPath
                    df = pickle.load(open(bestClusterPath[0],'rU'))
                    df_ts = pickle.load(open(bestClusterPath[1],'rU'))
                    res = GG_GGM_To_Table(tup_strain[1], df, df_ts,res)
                    print 'save GG, GGM to file : ', tup_strain[1]
                    i=i+1;
                    #if i > 6:
                     #   break
        res.to_csv('./GG_GGM/'+'RES.csv')
        

                
    
    @staticmethod
    def getCurrentPath():
        import os
        return os.getcwd()
    
    @staticmethod
    def listFilePathsByExt(folderPath, ext='',ext2=''):
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
            
        all_subfolders = next(os.walk(folderPath))[1]

        p=[]
        for fname in os.listdir(folderPath):
            if len(ext)>0 and fname.endswith(ext):
                p.append(os.path.join(folderPath,fname))
            if len(ext2)>0 and fname.endswith(ext2):
                p.append(os.path.join(folderPath,fname))
                
        if len(p)==2:
            match = GG_GGM_Selector.commonStringFromstart(os.path.split(p[0])[1],os.path.split(p[1])[1])
            p.append(match)
            
        paths.append(p)
        
        if len(all_subfolders):
            for sub in all_subfolders:
                p=[]
                subPath = os.path.join(folderPath,sub)
                for fname in os.listdir(subPath):
                    if len(ext)>0 and fname.endswith(ext):
                        p.append(os.path.join(subPath,fname))
                    if len(ext2)>0 and fname.endswith(ext2):
                        p.append(os.path.join(subPath,fname))
                if len(p)==2:
                    match = GG_GGM_Selector.commonStringFromstart(os.path.split(p[0])[1],os.path.split(p[1])[1])
                    p.append(match)
                    
                paths.append(p)
        paths_ = [x for x in paths if x]
                    
        return paths_
    
    
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
    def AIC_ggm_gg(df):
        from cohortMortalitySummary import GGMfit, GGfit
        ggmFit = GGMfit(df);
        ggFit = GGfit(df);  
        return (ggmFit['AIC'], ggFit['AIC'])
    
    @staticmethod
    def GG_GGM_To_file(df,df_ts,strain_name,folderPath):
        import os
        savePath = os.path.join(folderPath,strain_name+'.txt')
        fo = open(savePath,'w',0)
        from cohortMortalitySummary import GGfit, GGMfit
        from cohortMortalitySummary import KaplanMeier, NelsonAalen, Breslow, BSHazardR, generate_ts_image
        import pandas as pd
        KM = KaplanMeier(df, df_ts)
        NA = NelsonAalen(df, df_ts)
        upper=KM.loc[:,'upper_ci']
        lower = KM.loc[:,'lower_ci']
        var =  KM.loc[:,'variance']
        survivor=KM.loc[:,'survivorship']
        s = survivor.values >= 0.5    
        slower = lower.values >= 0.5
        supper = upper.values >= 0.5
        ind = survivor.iloc[s].index[-1];  
        indlower = lower.iloc[slower].index[-1]
        indupper = upper.iloc[supper].index[-1]
        data = [{'est': ind, 'L95%': indlower, 'U95%':indupper, 'var':var.loc[ind]}]
        lifespan = pd.DataFrame(data,index=['lifespan'],columns=['est','L95%','U95%','var'])
        print lifespan
        survivorship = KM.iloc[-1,:]['survivorship'];
        GG = GGfit(df)
        print >> fo, 'GGFit_params'
        if True:
            print >> fo, tup_strain[1]
            print >> fo, '-------------------------------\nGG(AIC={:1f})'.format(GG['AIC'])
            print >> fo, GG['model_paras']
            print >> fo, '\n'


        GGM = GGMfit(df)
        print >> fo, 'GGMFit_params'
        if True:
            print >> fo, tup_strain[1]
            print >> fo, '-------------------------------\nGGM(AIC={:1f})'.format(GGM['AIC'])
            print >> fo, GGM['model_paras'],'\n'
            print >> fo, lifespan
            print >> fo, 'survivorship'
            print >> fo, survivorship
            print >> fo, '\n'
    
    # 1) Number of tab files in the cluster
    @staticmethod
    def countTabFiles(df):
        list_fn=df.index.get_level_values(0).tolist()
        uniq_list = list(set(list_fn))
        print len(uniq_list)
        return len(uniq_list);

    # 2) Number of dates in the cluster
    @staticmethod
    def getDatetime(s):
        from datetime import datetime
        for k in s.split(' '):
            #print k
            if len(k) >= 8:
                for m in k.split('_'):
                    if m.isdigit():
                        t = datetime(int(m[4:]), int(m[2:4]), int(m[0:2]), 0, 0)
                        return t;
            return None;
        
    @staticmethod
    def tabFileDate(df):
        list_fn=df.index.get_level_values(0).tolist()
        uniq_list = list(set(list_fn))
        date_list = [];
        for fn in uniq_list:
            s=fn
            date_list.append(str(GG_GGM_Selector.getDatetime(s)))

        uniq_date = list(set(date_list))
        print len(uniq_date), uniq_date
        return (len(uniq_date),uniq_date)
    

    #   3) Final survivorship: influence both K-S based testing, the reliability of GGM parameter estimation, and reliability of experimental estimation of lifespan. We distrust those strains with high final survivorship
    #   	finalKM = KaplanMeier(df, ts).iloc[-1,:].loc[" survivorship"]
    @staticmethod
    def finalKM(df, df_ts):
        from cohortMortalitySummary import KaplanMeier
        finalKM = KaplanMeier(df, df_ts).iloc[-1,:].loc["survivorship"]
        return finalKM

    # 4) GGM Goodness-of-fit statistic (K-S), p-value, and fit proportions:
    @staticmethod
    def ggmStatistic(df, df_ts):
        from cohortMortalitySummary import GGMfit, KSm_gof, KSm_test
        import numpy as np
        ggmfit = GGMfit(df)
        residues = KSm_gof(df,df_ts, ggmfit['ML_survivorship'])
        test = KSm_test(residues[0],residues[1],alpha=0.02)
        KSstatistic = test[0]
        pvalue = test[1]
        fit_prop = np.sum(np.abs(residues[0].iloc[::2]) < test[2])/float(len(df_ts))
        print pvalue, fit_prop
        return (pvalue, fit_prop)

    #5) Whether Makeham term is appropriate: Compare GG and GGM using AIC and Likelihood-ratio test:
    @staticmethod
    def ggm_vs_gg(df):
        from cohortMortalitySummary import GGMfit, GGfit
        from scipy.stats import chi2
        ggmfit = GGMfit(df)
        try:
            ggfit = GGfit(df)
        except:
            ggfit = GG_GGM_Selector.makeGG_empty()
            
        makeham_AIC = ggmfit['AIC'] - ggfit['AIC']
        D = -2*(ggfit['logLik'] - ggmfit['logLik'])
        likelihoodratio_pvalue = chi2.sf(D,1)
        return makeham_AIC, likelihoodratio_pvalue;

    #6) functions of GGM parameters: initial mortality (h0 = lambda+b*s/beta), the mortality plateaux (hmax = b*s+lambda), CV(b, beta, s, labmda)
    @staticmethod
    def ggmStats(df):
        from cohortMortalitySummary import GGMfit_flxsrv_fit, GGfit_flxsrv_fit
        import pandas as pd
        import numpy as np
        flxsrv_ggm_fit = GGMfit_flxsrv_fit(df)
        #flxsrv_gg_fit = GGfit_flxsrv_fit(df)
        import rpy2
        import os
        os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'
        import rpy2.robjects as ro
        import rpy2.robjects.packages as rpacks
        from rpy2.robjects import pandas2ri, Formula, globalenv
        #import pandas.rpy.common as com
        pandas2ri.activate()
        rflexsurv = rpacks.importr('flexsurv')
        flexsurvreg = rflexsurv.flexsurvreg
        #summary = rflexsurv.summary_flexsurvreg
        normboot = rflexsurv.normboot_flexsurvreg

        flxsrv_gg_fit = GGfit_flxsrv_fit(df)
        flxsrv_ggm_fit = GGMfit_flxsrv_fit(df)
        sample_ggm = pd.DataFrame(data=pandas2ri.ri2py( normboot(flxsrv_ggm_fit, B=1000000, transform=True)),columns=['rate','beta','s','lambda'])
        sample_gg = pd.DataFrame(data=pandas2ri.ri2py( normboot(flxsrv_gg_fit, B=1000000, transform=True)),columns=['beta','s','rate'])

        if flxsrv_gg_fit.rx2('AIC')[0] < flxsrv_ggm_fit.rx2('AIC')[0]:
            sample_hinit = sample_gg['rate']+sample_gg['s']-sample_gg['beta']
            sample_hmax = sample_gg['rate']+sample_gg['s']
        else:
            sample_hinit = np.log(np.exp(sample_ggm['rate']+sample_ggm['s']-sample_ggm['beta'])+np.exp(sample_ggm['lambda']))
            sample_hmax = np.log(np.exp(sample_ggm['rate']+sample_ggm['s'])+np.exp(sample_ggm['lambda']))

        from scipy.stats import norm
        stats = dict() 
        for (name,sample) in {'h_init':sample_hinit, 'h_max':sample_hmax}.items():
            mu = np.mean(sample)
            sigma = np.std(sample)
            stats[name]={'est':np.exp(mu), 
                        'L95%':np.exp(mu+norm.ppf(0.025)*sigma),
                        'U95%':np.exp(mu+norm.ppf(0.975)*sigma),
                        'se':np.exp(mu)*sigma}

        return pd.DataFrame(data=stats).T.loc[:,['est','L95%','U95%','se']]

    # 7) GammaGompertz parameters and their functions: b, beta, s, CV(b, beta, s)
    @staticmethod
    def ggParameterStats(df):
        from cohortMortalitySummary import GGMfit_flxsrv_fit, GGfit_flxsrv_fit
        import pandas as pd
        import numpy as np
        from mpmath import mp, hyp2f1, hyp3f2, power,sqrt
        mp.dps = 16; mp.pretty = True

        flxsrv_ggm_fit = GGMfit_flxsrv_fit(df)
        #flxsrv_gg_fit = GGfit_flxsrv_fit(df)
        import rpy2
        import os
        os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'
        import rpy2.robjects as ro
        import rpy2.robjects.packages as rpacks
        from rpy2.robjects import pandas2ri, Formula, globalenv
        #import pandas.rpy.common as com
        pandas2ri.activate()
        rflexsurv = rpacks.importr('flexsurv')
        flexsurvreg = rflexsurv.flexsurvreg
        #summary = rflexsurv.summary_flexsurvreg
        normboot = rflexsurv.normboot_flexsurvreg

        flxsrv_gg_fit = GGfit_flxsrv_fit(df)
        #flxsrv_ggm_fit = GGMfit_flxsrv_fit(df)
        #sample_ggm = pd.DataFrame(data=pandas2ri.ri2py( normboot(flxsrv_ggm_fit, B=10000, transform=True)),columns=['rate','beta','s','lambda'])
        sample_gg = pd.DataFrame(data=pandas2ri.ri2py( normboot(flxsrv_gg_fit, B=2000, transform=True)),columns=['beta','s','rate'])

        para_stats_gg = pd.DataFrame(columns=['mean','variance','cv'],index=sample_gg.index)
        for i in sample_gg.index:
            (beta, s, b) = np.exp(sample_gg.loc[i,:])
            mean = hyp2f1(s,1,s+1,(beta-1)/beta)/b/s
            var = 2*power(beta,s)*hyp3f2(s,s,s,s+1,s+1,(1-beta))/(b*b)/(s*s) - mean*mean
            para_stats_gg.loc[i,'mean'] = mean
            para_stats_gg.loc[i,'variance'] = var
            para_stats_gg.loc[i,'cv'] = sqrt(var)/mean

        from scipy.stats import norm
        stats = dict() 
        for (name,sample) in {'mean':para_stats_gg['mean'], 'variance':para_stats_gg['variance'], 'cv':para_stats_gg['cv']}.items():
            mu = np.mean(sample)
            sigma = np.std(sample)
            stats[name]={'est':mu, 
                        'L95%':mu+norm.ppf(0.025)*sigma,
                        'U95%':mu+norm.ppf(0.975)*sigma,
                        'se':sigma}
        stats = pd.DataFrame(data=stats).T.loc[:,['est','L95%','U95%','se']]
        samples = pd.concat([np.exp(sample_gg.loc[:,['rate','beta','s']]), para_stats_gg], axis=1)
        return (stats, samples)

    # 8) Experimental lifespan stats: mean, median, CV
    @staticmethod
    def stats_lifespan(df):
        import numpy as np
        m=np.mean(df['mortality'].values)
        std=np.std(df['mortality'].values)
        cv = std/m
        return [m, std, cv];

    @staticmethod
    def GG_GGM_To_file_More_Results(strain_name,df,df_ts):
        savePath = './GG_GGM/'+strain_name+'.txt'
        import os
        if os.path.isfile(savePath):
            fo = open(savePath,'a')
        else:
            fo = open(savePath, 'w',0)
        from cohortMortalitySummary import GGfit, GGMfit
        from cohortMortalitySummary import KaplanMeier, NelsonAalen, Breslow, BSHazardR, generate_ts_image
        import pandas as pd
        KM = KaplanMeier(df, df_ts)
        NA = NelsonAalen(df, df_ts)
        upper=KM.loc[:,'upper_ci']
        lower = KM.loc[:,'lower_ci']
        var =  KM.loc[:,'variance']
        survivor=KM.loc[:,'survivorship']
        s = survivor.values >= 0.5    
        slower = lower.values >= 0.5
        supper = upper.values >= 0.5
        ind = survivor.iloc[s].index[-1];  
        indlower = lower.iloc[slower].index[-1]
        indupper = upper.iloc[supper].index[-1]
        data = [{'est': ind, 'L95%': indlower, 'U95%':indupper, 'var':var.loc[ind]}]
        lifespan = pd.DataFrame(data,index=['lifespan'],columns=['est','L95%','U95%','var'])
        print lifespan
        survivorship = KM.iloc[-1,:]['survivorship'];

        GG = GGfit(df)
        print >> fo, 'GGFit_params'
        if True:
            print >> fo, tup_strain[1]
            print >> fo, '-------------------------------\nGG(AIC={:1f})'.format(GG['AIC'])
            print >> fo, GG['model_paras']
            print >> fo, '\n'


        GGM = GGMfit(df)
        print >> fo, 'GGMFit_params'
        if True:
            print >> fo, tup_strain[1]
            print >> fo, '-------------------------------\nGGM(AIC={:1f})'.format(GGM['AIC'])
            print >> fo, GGM['model_paras'],'\n'
            print >> fo, lifespan
            print >> fo, 'survivorship'
            print >> fo, survivorship
            print >> fo, '\n'

        print >> fo, '1) Number of tab files in the cluster'
        print >> fo, countTabFiles(df)
        print >> fo, '\n'


        print >> fo, '2) Number of dates in the cluster'
        dt = tabFileDate(df)
        print >> fo, dt[0],dt[1]
        print >> fo, '\n'

        print >> fo, '3) Final survivorship:'
        fKM = finalKM(df, df_ts)
        print >> fo, fKM
        print >> fo, '\n'

        print >> fo, '4) GGM Goodness-of-fit statistic (K-S), p-value, and fit proportions:'
        print >> fo, 'pvalue','\t','fit_proportion'
        print >> fo, ggmStatistic(df, df_ts)
        print >> fo, '\n'

        print >> fo, '5) AIC(GGM)-AIC(GG), likelihood-ratio test pvalue'
        dAIC_ggm_minus_gg , pvalue_lklh_r = ggm_vs_gg(df)
        print >> fo, dAIC_ggm_minus_gg , pvalue_lklh_r
        print >> fo, '\n'

        print >> fo, '6) initial mortality (h_init, h_max)'
        print >> fo, ggmStats(df)
        print >> fo, '\n'

        print >> fo, '7) GammaGompertz parameters and their functions: b, beta, s, CV(b, beta, s)'
        if dAIC_ggm_minus_gg > 0 and fKM<0.25:
            try:
                (stats, samples) = ggParameterStats(df)
                lastRow = samples.loc[len(samples)-1]
                print >> fo, stats
                print >> fo,'\n'
                print >> fo, lastRow
            except:
                print "Error in ggParameterStats()"
            #samples.to_csv('./GG_GGM/'+strain_name+'_GammaGompertz.csv')
        elif not dAIC_ggm_minus_gg > 0:
            print >> fo, "N/A, choose GGM due to AIC"
        else:
            print >> fo, "finalKM >=0.25, theoretical stats not meaningful"
        print >> fo, '\n'


        print >> fo, '8) lifespan stats'
        print >> fo, 'mean', '\t', 'std','\t','cv'
        print >> fo, stats_lifespan(df)
        print >> fo, '\n'

    @staticmethod
    def dfToOneLineDF(d,name,key_ext=''):
        import pandas as pd
        out = {};
        l = [];
        for ind in d.index:
            for key in d.keys():
                k_new = key_ext+ind+'_'+key
                out[k_new] = d.loc[ind][key]
                l.append(k_new)
                #print ind, key, k.loc[ind][key]
        df_new = pd.DataFrame(out,[name],columns=l)
        return df_new

    @staticmethod
    def concatTwoDf(df1, df2):
        import pandas as pd
        dfs= [df1, df2]
        df  = pd.concat(dfs,axis=1)
        return df

    @staticmethod
    def addColumnToDF(df, value, key):
        import pandas as pd
        df[key] = pd.Series(value, index=df.index)

    @staticmethod
    def makeGG_empty():
        import pandas as pd
        import numpy as np
        model_para_except = pd.DataFrame(np.full((3,4),-1), index=['beta','s','rate'], columns=['est', 'L95%', 'U95%', 'se'])
        gg_except = {'AIC': -1, 
                    'logLik': -1,
                    'events': -1,
                    'tRisk': -1, 
                    'model_paras': model_para_except,
                    'ML_survivorship': -1,
                    'ML_cumulative_hazard': -1,
                    'ML_hazard': -1}
        return gg_except
    
    @staticmethod
    def GG_GGM_To_Table(strain_name,df,df_ts, res, writeIndividual=True,out_folder_path=''):
        import sys,os
        sys.path.insert(1, os.path.join(sys.path[0], '../../..'))
        from cohortMortalitySummary import GGfit, GGMfit
        from cohortMortalitySummary import KaplanMeier, NelsonAalen, Breslow, BSHazardR, generate_ts_image
        import pandas as pd
        KM = KaplanMeier(df, df_ts)
        NA = NelsonAalen(df, df_ts)
        upper=KM.loc[:,'upper_ci']
        lower = KM.loc[:,'lower_ci']
        var =  KM.loc[:,'variance']
        survivor=KM.loc[:,'survivorship']
        s = survivor.values >= 0.5    
        slower = lower.values >= 0.5
        supper = upper.values >= 0.5
        ind = survivor.iloc[s].index[-1];  
        indlower = lower.iloc[slower].index[-1]
        indupper = upper.iloc[supper].index[-1]
        data = [{'est': ind, 'L95%': indlower, 'U95%':indupper, 'var':var.loc[ind]}]
        lifespan = pd.DataFrame(data,index=['lifespan'],columns=['est','L95%','U95%','var'])
        #print lifespan
        survivorship = KM.iloc[-1,:]['survivorship'];
        
        GG_run = True
        try:
            GG = GGfit(df)
            if True:
                df_GG= GG['model_paras']
                gg_oneline = GG_GGM_Selector.dfToOneLineDF(df_GG,strain_name,'gg_')
                res_oneline=gg_oneline;
        except:
            print "error in ggfit"
            GG = GG_GGM_Selector.makeGG_empty()
            res_oneline = GG_GGM_Selector.dfToOneLineDF(GG['model_paras'],strain_name,'gg_')
            GG_run = False
        
        if writeIndividual:
            import os
            savePath = os.path.join(out_folder_path,strain_name+'.txt')
            print '-----------',savePath
            fo = open(savePath, 'w',0)
         
        if writeIndividual:
            print >> fo, 'GGFit_params'
            print >> fo, strain_name
            print >> fo, '-------------------------------\nGG(AIC={:1f})'.format(GG['AIC'])
            print >> fo, GG['model_paras']
            print >> fo, '\n'


        GGM = GGMfit(df)
        if True:
            print GGM['AIC']
            ggm_oneline = GG_GGM_Selector.dfToOneLineDF(GGM['model_paras'],strain_name,'ggm_')
            res_oneline = GG_GGM_Selector.concatTwoDf(res_oneline,ggm_oneline);
            lifespan_oneline = GG_GGM_Selector.dfToOneLineDF(lifespan, strain_name)
            res_oneline = GG_GGM_Selector.concatTwoDf(res_oneline, lifespan_oneline);
            GG_GGM_Selector.addColumnToDF(res_oneline,survivorship,'suvrivorship')
            GG_GGM_Selector.addColumnToDF(res_oneline,GGM['AIC'],'AIC')
        

        if writeIndividual:
            print >> fo, 'GGMFit_params'
            print >> fo, strain_name
            print >> fo, '-------------------------------\nGGM(AIC={:1f})'.format(GGM['AIC'])
            print >> fo, GGM['model_paras'],'\n'
            print >> fo, lifespan
            print >> fo, 'survivorship'
            print >> fo, survivorship
            print >> fo, '\n'

        if True:
            count_tabs =GG_GGM_Selector.countTabFiles(df)
            GG_GGM_Selector.addColumnToDF(res_oneline,count_tabs,'count_tab_files');
            
        if writeIndividual:
            print >> fo, '1) Number of tab files in the cluster'
            print >> fo, count_tabs
            print >> fo, '\n'

        if True:
            dt = GG_GGM_Selector.tabFileDate(df);
            GG_GGM_Selector.addColumnToDF(res_oneline,dt[0],'count_dates')
        
        if writeIndividual:
            print >> fo, '2) Number of dates in the cluster'
            print >> fo, dt[0],dt[1]
            print >> fo, '\n'

        if True:
            fKM = GG_GGM_Selector.finalKM(df, df_ts);
            GG_GGM_Selector.addColumnToDF(res_oneline,fKM, 'final_survivorship')
        
        if writeIndividual:
            print >> fo, '3) Final survivorship:'
            print >> fo, fKM
            print >> fo, '\n'

        if True:
            ggmS = GG_GGM_Selector.ggmStatistic(df, df_ts);
            GG_GGM_Selector.addColumnToDF(res_oneline,ggmS[0], 'p-value_ggm_K-S');
            GG_GGM_Selector.addColumnToDF(res_oneline,ggmS[1], 'fit_proportion_ggm_K-S');
        
        if writeIndividual:
            print >> fo, '4) GGM Goodness-of-fit statistic (K-S), p-value, and fit proportions:'
            print >> fo, 'pvalue','\t','fit_proportion'
            print >> fo, ggmS
            print >> fo, '\n'

        if True:
            dAIC_ggm_minus_gg , pvalue_lklh_r = GG_GGM_Selector.ggm_vs_gg(df)
            GG_GGM_Selector.addColumnToDF(res_oneline,dAIC_ggm_minus_gg, 'diff_AIC_GGM-GG');
            GG_GGM_Selector.addColumnToDF(res_oneline,pvalue_lklh_r, 'likelihood_ration_pvalue');
        
        if writeIndividual:
            print >> fo, '5) AIC(GGM)-AIC(GG), likelihood-ratio test pvalue'
            print >> fo, dAIC_ggm_minus_gg , pvalue_lklh_r
            print >> fo, '\n'

        if True:
            try:
                ggm_stats=GG_GGM_Selector.ggmStats(df)
                initial_mortality = GG_GGM_Selector.dfToOneLineDF(ggm_stats,strain_name,'ggm_initial_mortality_')
                res_oneline = GG_GGM_Selector.concatTwoDf(res_oneline,initial_mortality)
        
                if writeIndividual:
                    print >> fo, '6) initial mortality (h_init, h_max)'
                    print >> fo, ggm_stats
                    print >> fo, '\n'
            except:
                print("\x1b[31m\"ERROR in step ggmStats(df)\"\x1b[0m")
                if writeIndividual:
                    print >> fo, '6) initial mortality (h_init, h_max)'
                    print >> fo, 'ERROR'
                    print >> fo, '\n'

        if True:
            error_ggParameterStats = False
            if GG_run and dAIC_ggm_minus_gg > 0 and fKM < 0.25:
                try:
                    (stats, samples) = GG_GGM_Selector.ggParameterStats(df)
                    stats_oneline_dAIC = GG_GGM_Selector.dfToOneLineDF(stats, strain_name,'gg_params_')
                    res_oneline = GG_GGM_Selector.concatTwoDf(res_oneline,stats_oneline_dAIC);
                    
                    lastRow = samples.loc[len(samples)-1]
                    lastRow=pd.DataFrame(lastRow.values,lastRow.index)
                    lastRow=lastRow.transpose()
                    lastRow = pd.DataFrame(lastRow.values, index=['gg_params_'],columns=lastRow.keys())
                    lastRow_oneline = GG_GGM_Selector.dfToOneLineDF(lastRow, strain_name)
                    res_oneline = GG_GGM_Selector.concatTwoDf(res_oneline, lastRow_oneline);
                except:
                    error_ggParameterStats = True
                    print "Error in ggParameterStats"

            elif not dAIC_ggm_minus_gg > 0:
                print "N/A, choose GGM due to AIC"
            else:
                print "finalKM >=0.25, theoretical stats not meaningful"
        
        if writeIndividual:
            print >> fo, '7) GammaGompertz parameters and their functions: b, beta, s, CV(b, beta, s)'
            if GG_run and dAIC_ggm_minus_gg > 0 and fKM<0.25:
                if not error_ggParameterStats:
                    print >> fo, stats
                    print >> fo,'\n'
                    print >> fo, lastRow
            elif not dAIC_ggm_minus_gg > 0:
                print >> fo, "N/A, choose GGM due to AIC"
            else:
                print >> fo, "finalKM >=0.25, theoretical stats not meaningful"
            print >> fo, '\n'
        

        if True:
            m,std,cv =GG_GGM_Selector.stats_lifespan(df);
            GG_GGM_Selector.addColumnToDF(res_oneline, m,'mean_lifespan');
            GG_GGM_Selector.addColumnToDF(res_oneline, std,'std_lifespan');
            GG_GGM_Selector.addColumnToDF(res_oneline, cv,'cv_lifespan');
            
        if writeIndividual:   
            print >> fo, '8) lifespan stats'
            print >> fo, 'mean', '\t', 'std','\t','cv'
            print >> fo, m,'\t',std,'\t',cv
            print >> fo, '\n'

        if res is None:
            res = res_oneline;
        else:
            res = res.append(res_oneline,ignore_index=False);

        return res


    @staticmethod
    def run0():
        import csv
        import pickle

        strainParentPath = '../Microfluidics Database/'
        res={}
        with open('microfluidics database information_07102017_UPDATED.csv','rU') as csvfile:
            reader = csv.reader(csvfile, dialect=csv.excel_tab )
            for row in reader:
                tup_strain = findStrainLine(row[0]);
                if tup_strain[0]:
                    bestClusterPath = getClusterFilePath(strainParentPath, tup_strain[1], tup_strain[2])
                    print bestClusterPath
                    df = pickle.load(open(bestClusterPath[0],'rU'))
                    df_ts = pickle.load(open(bestClusterPath[1],'rU'))
                    GG_GGM_To_file_More_Results(tup_strain[1], df, df_ts)
                    print 'save GG, GGM to file : ', tup_strain[1]
                    break

