
# coding: utf-8

# In[162]:

class Dendrogram_Replicates(object):
    ext_df_cluster0 = 'mortality_dataframe_cluster0.pickle'
    ext_timegride_cluster0 = 'relative_timegrid_cluster0.pickle'
    ext_df_timegrid ='df_relative_timegrid.pickle'
    ext_outlier = 'outlier_tabs.txt'
    pre = 'DR-'
    
    def __init__(self, folder_path='', file_replicates_path='',select_cluster0=True, showProgress=True):
        """
        :param select_cluster0: if False, all the df will be selected, if True, outlier is automatcally excluded
        :param outlier: if select_cluster0 is False, the outlier files will be excluded(if there is any outlier.txt)
        
        """
        if len(folder_path):
            self.folder_path = folder_path
        else:
            self.folder_path = Dendrogram_Replicates.getCurrentPath()
        
        if len(file_replicates_path) == 0:
            print "The replicate file list is not provided!"
            return
        self.file_replicates_path = file_replicates_path
        self.select_cluster0 = select_cluster0
        self.showProgress=showProgress
        
        if select_cluster0:
            path_list = Dendrogram_Replicates.listFilePathsByExt(self.folder_path, Dendrogram_Replicates.ext_df_cluster0,Dendrogram_Replicates.ext_timegride_cluster0)
        else:
            path_list = Dendrogram_Replicates.listFilePathsByExt(self.folder_path, Dendrogram_Replicates.ext_df_timegrid,Dendrogram_Replicates.ext_outlier)
        self.path_list = path_list
        self.select_cluster0
    
    def run(self):
        self.groupReplicate(self.path_list, self.file_replicates_path,self.select_cluster0)
        
        
    def groupReplicate(self,path_list,file_path,select_cluster0):
        import os
        dataset_base = {}
        data2add = dict()
        breakpoint = 0
        for l in path_list:
            #if breakpoint >=1:      
            #    break    
            breakpoint = breakpoint + 1
            print "\n finding strain ", l[-1]," in replicate list file"
            strain_name = l[-1]
            df = None
            timegrid = None
            l_group = Dendrogram_Replicates.groupByFilterFile(file_path, strain_name)
            outlier_l = []
            if select_cluster0:
                for m in l:
                    if Dendrogram_Replicates.ext_df_cluster0 in m:
                        df = Dendrogram_Replicates.loadPickle(m)
                    if Dendrogram_Replicates.ext_timegride_cluster0 in m:
                        timegrid = Dendrogram_Replicates.loadPickle(m)
                        timegrid_path = m
            else:
                df_s = []
                timegrid_s = []
                timegrid_path_s=[]
                for m in l:
                    if Dendrogram_Replicates.ext_df_timegrid in m:
                        df_timegrid = Dendrogram_Replicates.loadPickle(m)
                        df_s.append(df_timegrid.df)
                        timegrid_s.append(df_timegrid.time_grid_relative)
                        timegrid_path_s.append(m)
                    if Dendrogram_Replicates.ext_outlier in m:
                        outlier_l = Dendrogram_Replicates.getListFromFile(m)
                self.df_s = df_s
                self.timegrid_s = timegrid_s
                self.outlier_l = outlier_l
                
 
            if select_cluster0 and (df is None or timegrid is None):
                print " ERROR, df or timegrid not found"
                continue
            
            
            for idx, val in enumerate(l_group):
                if not select_cluster0:
                    df = None
                    timegrid = None
                    ind_df = 0
                    for df_loop in df_s:
                        fns = list(df_loop.index.levels[0])
                        if val[0] in fns:
                            df = df_loop
                            timegrid = timegrid_s[ind_df]
                            timegrid_path = timegrid_path_s[ind_df]
                            break
                        ind_df = ind_df + 1
                    for l_out in outlier_l:
                        df = df[df.fn != l_out]

                
                if idx==0:
                    if select_cluster0:
                        df_new = Dendrogram_Replicates.sliceDataFrame(val,df)
                    else: # No need to slice because the df_timegrid is already separated by division file
                        df_new = df
                    if df_new is None:
                        continue
                    dataset_base[strain_name+'-'+str(idx)] = dict();
                    dataset_base[strain_name+'-'+str(idx)]['fns'] = strain_name+'-'+str(idx)
                    dataset_base[strain_name+'-'+str(idx)]['dataframe'] = df_new
                    dataset_base[strain_name+'-'+str(idx)]['timegrid'] = list(timegrid);
                    dataset_base[strain_name+'-'+str(idx)]['strain'] = strain_name
                    dataset_base[strain_name+'-'+str(idx)]['exp'] = timegrid_path
                    dataset_base[strain_name+'-'+str(idx)]['date'] = '20140602'
                else:
                    if select_cluster0:
                        df_new = Dendrogram_Replicates.sliceDataFrame(val,df)
                    else:
                        df_new = df
                    if df_new is None:
                        continue                    
                    data2add[strain_name+'-'+str(idx)] = dict()
                    data2add[strain_name+'-'+str(idx)]['fns'] = strain_name+'-'+str(idx)
                    data2add[strain_name+'-'+str(idx)]['dataframe'] = df_new
                    data2add[strain_name+'-'+str(idx)]['timegrid'] = list(timegrid);
                    data2add[strain_name+'-'+str(idx)]['strain'] = strain_name
                    data2add[strain_name+'-'+str(idx)]['exp'] = timegrid_path
                    data2add[strain_name+'-'+str(idx)]['date'] = '20140602'
                
            

        #df_base = Dendrogram_Replicates.combineDataSet(dataset_base, data2add)
        saveName = Dendrogram_Replicates.pre + 'alldatasets.pickle'
        savePath = os.path.join(self.folder_path, saveName)
        self.alldatasets, self.names = Dendrogram_Replicates.updateDataset(dataset_base, data2add,savePath)
        #Dendrogram_Replicates.runPdist()
        
        saveName = Dendrogram_Replicates.pre + 'pdist.npy'
        savePath = os.path.join(self.folder_path, saveName)
        self.pdist = Dendrogram_Replicates.getPdist(self.alldatasets, self.names, savePath)
        
        saveName = Dendrogram_Replicates.pre + 'dendrogram_all.eps'
        savePath = os.path.join(self.folder_path,saveName)
        Dendrogram_Replicates.plotDendrogram(self.pdist,self.alldatasets,savePath)

        
    @staticmethod
    def getListFromFile(path):
        import csv
        l=[]
        with open(path, 'rU') as readf:
            csvreader = csv.reader(readf)
            for row in csvreader:
                if len(row) > 0:
                    l.append(row[0])  

        return l

    @staticmethod
    def updateDataset(dataset_base, data2add,savePath=''):
        alldatasets = dataset_base.copy()
        alldatasets.update(data2add)
        names = alldatasets.keys()
        import pickle
        if len(savePath):
            print 'save(alldagtasets) to : ', savePath
            with open(savePath, 'wb') as handle:
                pickle.dump(alldatasets, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return alldatasets, names
    
    @staticmethod
    def defRfunc():
        import rpy2,os
        os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'
        import rpy2.robjects as ro
        import rpy2.robjects.packages as rpacks
        from rpy2.robjects import pandas2ri, Formula, globalenv
        pandas2ri.activate()
        rsurvival = rpacks.importr('survival')
        rflexsurv = rpacks.importr('flexsurv')
        custom_env =  ro.r['new.env']()
        ro.r.source('Rmodels/GammaGompertzMakeham.R',local=custom_env, echo=False,verbose=False);
        ro.r.attach(custom_env)
        gammagompertzmakeham = custom_env['gammagompertzmakeham']
        flexsurvreg = rflexsurv.flexsurvreg
        
        from scipy.stats import chi2
        def loglik_test(fsf_null, fsf_hy):
            D = -2*(fsf_null.rx2('loglik')[0] - fsf_hy.rx2('loglik')[0])
            df = fsf_hy.rx2('npars')[0]- fsf_null.rx2('npars')[0]
            return  D, df, chi2.sf(D,df)
        
    @staticmethod
    def plotDendrogram(pdist,alldatasets, savePath=''):
        from cohortMortalitySummary import _criticalPartialBrownianBridge
        from scipy.cluster.hierarchy import dendrogram,complete
        cutoff31 = 1.949
        cutoff32 = 1.858
        cutoff35 = 1.731
        cutoff21 = 1.628
        cutoff22 = 1.517
        cutoff25 = 1.358
        cutoff = _criticalPartialBrownianBridge(0.01/float(len(pdist)))
        import matplotlib.pyplot as plt
        
        fig = plt.figure(figsize=(200,400))
        dn = dendrogram(complete(pdist),labels=alldatasets.keys(),orientation='right', color_threshold=cutoff22,leaf_font_size=58)
        fig.set_size_inches(105, 185)
        if len(savePath):
            print "save(dendrogram) to: ", savePath
            fig.savefig(savePath)
        
    @staticmethod
    def plotDendrogramFromFile(dataset_file_path,pdist_file_path, savePath=''):
        import pickle
        import numpy as np
        with open(dataset_file_path, 'rb') as handle:
            alldata = pickle.load(handle)
        pdistance=np.load(pdist_file_path)
        Dendrogram_Replicates.plotDendrogramFromFile(pdist,alldata,savePath)

        
    @staticmethod
    def getPdist(alldatasets,names,savePath=''):
        from cohortMortalitySummary import _criticalPartialBrownianBridge, KSm_2samples, KSm_test
        import numpy as np
        from scipy.special import binom
        distances = np.zeros((len(names),len(names)))
        significances = np.empty((len(names),len(names)))
        pdist = np.zeros(int(binom(len(names), int(2))))
        significances[:] = np.NAN
        significances_pd = np.zeros(int(binom(len(names), int(2))))
        k = 0
        for i in range(len(names)):
            for j in range(i+1,len(names)):
                print names[i], names[j]
                residues = KSm_2samples(alldatasets[names[i]]['dataframe'],alldatasets[names[i]]['timegrid'], 
                                        alldatasets[names[j]]['dataframe'],alldatasets[names[j]]['timegrid'])
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
        if len(savePath):
            print 'save(pdist) to: ',savePath
            np.save(savePath, pdist)
        return pdist
        
        
    @staticmethod
    def combineDataSet(dataset_base, data2add):
        import pandas as pd
        def creatSeries(index, v):
            l = [];
            for i in range(0,len(index)):
                l.append(v);
            return pd.Series(l, index)

        for key, dct in dataset_base.items():
            print key
            dct['dataframe']['exp'] = creatSeries(dct['dataframe']['x'].index,dct['exp'])
            dct['dataframe']['strain'] = creatSeries(dct['dataframe']['x'].index,dct['strain'])
            dct['dataframe']['date'] =creatSeries(dct['dataframe']['x'].index,dct['date'])


        for key, dct in data2add.items():
            print key
            dct['dataframe']['exp'] = creatSeries(dct['dataframe']['x'].index,dct['exp'])
            dct['dataframe']['strain'] = creatSeries(dct['dataframe']['x'].index,dct['strain'])
            dct['dataframe']['date'] =creatSeries(dct['dataframe']['x'].index,dct['date'])
        df_base = pd.concat( [ dct['dataframe'] for key, dct in dataset_base.items()] )
        return df_base
    
    
    
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
        f = open(fPath, 'rb')
        p = pickle.load(f)
        return p
    
    @staticmethod
    def sliceDataFrame(gList,df):
        """
        works only for the ordered index
        :params gList: list to slice df
        :params outlier_l: list to be excluded
        """
        import pandas as pd
        index_df = df.index;
        levels_0=df.index.levels[0]
        frames = []
        
        levels_0=df.index.levels[0]
        l_levels_0 = list(levels_0)
        
        for gLine in gList:
            if gLine in l_levels_0:
                mask, idx = index_df.get_loc_level(gLine,level=index_df.names[0])
                selectedDf = df.iloc[mask.start:mask.stop]
                frames.append(selectedDf)
        if len(frames) > 1:
            return pd.concat(frames)
        elif len(frames) == 1:
            return selectedDf   

    @staticmethod
    def groupByFilterFile(filterFile,name):
        '''
        get all folder names from filterFile, if files are separated by return label,
        multiple folder names will be generated,
        :param filterFile: contains file paths list
        :param name: get the name included list
        :return: list paths in replicates
        '''
        import glob,os
        listFiles = [];
        with open(filterFile, 'rU') as rf:
            for line in rf:
                listFiles.append(line.rstrip());
        listFiles.append('')#this line is important, if the line of file is not seprated by the next line
        previous = ''
        count = 0
        folderGroup = [];
        out = []
        for j in range(0,len(listFiles)):
            line_fn = listFiles[j]
            if len(line_fn):
                if name in line_fn:
                    if previous == '':
                        count=count+1
                        if len(folderGroup):
                            out.append(folderGroup)
                            folderGroup=[]
                    folderGroup.append(os.path.split(line_fn)[1])
                    #print count, len(folderGroup), line_fn
            previous = line_fn.rstrip()
            if previous == '':
                if len(folderGroup):
                    out.append(folderGroup)
                    folderGroup=[]
        return out
    
    @staticmethod
    def getCurrentPath():
        import os
        return os.getcwd()
    
    @staticmethod
    def listFilePathsByExt(folderPath, ext='',ext2='',ext3=''):
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
        if not folderPath[-1].isalnum():
            folderPath = folderPath[:-1]
            
        all_subfolders = next(os.walk(folderPath))[1]

        p=[]
        for fname in os.listdir(folderPath):
            if len(ext)>0 and fname.endswith(ext):
                p.append(os.path.join(folderPath,fname))
            if len(ext2)>0 and fname.endswith(ext2):
                p.append(os.path.join(folderPath,fname))
            if len(ext3)>0 and fname.endswith(ext3):
                p.append(os.path.join(folderPath,fname))
        if len(p):
            folder_name = os.path.split(folderPath)[1]
            p.append(folder_name)
        
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
                    if len(ext3)>0 and fname.endswith(ext3):
                        p.append(os.path.join(folderPath,fname))
                if len(p):
                    folder_name = os.path.split(subPath)[1]
                    p.append(folder_name)
                paths.append(p)                
        paths_ = [x for x in paths if x]
                    
        return paths_

