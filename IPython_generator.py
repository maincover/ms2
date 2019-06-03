class IPythonNotebookGenerator(object):
    text_PIP = """\
# PIProfile Analyzer Process (all .tab files in the folder pth)
TabPostAnalyzer.py is should be in the grand parent folder
This is an auto-generated notebook."""
    text_TAB = """\
# TAB Analyzer Process (only file filter listed .tab files (existed also in pth) )
TabPostAnalyzer.py is should be in the grand parent folder
This is an auto-generated notebook."""
    text_PRE = """\
# Preprocessing to detect duplication
This is an auto-generated notebook."""
    text_MD="""\
# Mortality Distribution, generate mortality curve, compare with reference(WT for example)
This is an auto-generated notebook."""
    text_CR="""\
# Check Replication and discard orphan files by clustering
    """
    text_SELECTOR="""\
# GG_GGM_SELECTOR Write the result to table
    """
    
    text_DEN="""\
# Check replicates by Dendrogram, result saved to the input folder,
  Consider Cluster0 or all dataframe
    """
    
    code_part_PIP = [None]*10
    code_part_TAB = [None]*10
    code_part_PIP[0]=['%matplotlib inline','\n',
                'import sys,os','\n',
                'python_class_path=\'','\'' ,'\n',
                'sys.path.insert(1, os.path.join(sys.path[0], python_class_path))','\n',
                'from TabPostAnalyzer import PIProfileAnalyzer','\n',
                'pth=\'','\'' ,'\n',
                'pA = PIProfileAnalyzer(pth,"",debug=True,redefineInCells=False)']
    
    code_part_PIP[1]=['step0PlotOriginalFiles = True','\n',
                     'step1TvReg = False','\n',
                     'step2Output = False','\n',
                     'pA.process(step0PlotOriginalFiles, step1TvReg, step2Output)']
    
    code_part_PIP[2]=['step0PlotOriginalFiles = False','\n',
                     'step1TvReg = True','\n',
                     'step2Output = False','\n',
                     'pA.process(step0PlotOriginalFiles, step1TvReg, step2Output)']
    
    code_part_PIP[3] = ['pA.denoisePlotAll(asp=0.005,onlyPlot=True)','\n',
                       'pA.defineAlcoholTime(more=20)']
    
    code_part_PIP[4] = ['step0PlotOriginalFiles = False','\n',
                       'step1TvReg = False','\n',
                       'step2Output = True','\n',
                       'pA.process(step0PlotOriginalFiles, step1TvReg, step2Output)']
    
    code_part_PIP[5] = ['pA.plotAllFiles(0.02)']
    
    code_part_PIP[6] = ['import pickle','\n',
                       'f0 = open(os.path.join(pA.path,pA.saveName, \'PIP_\'+pA.saveName+\' _df_relative_timegrid.pickle\'), \'rb\')','\n',
                       'df_ts = pickle.load(f0)','\n',
                       'df_ts.checkLength()']

    insert_position_1 = 7; #deprecated
    
    
    
    code_part_TAB[0] = ['%matplotlib inline','\n',
               'import sys,os','\n',
               'python_class_path=\'','\'' ,'\n',
               'sys.path.insert(1, os.path.join(sys.path[0], python_class_path))','\n',
               'from TabPostAnalyzer import TabPostAnalyzer','\n',
               'pth=\'','\'' , '\n',
                'fFn=\'','\'','\n',
                'partial=','\n',
                'pA = TabPostAnalyzer(pth,filterFile=fFn, group=partial)']
    
    code_part_TAB[1] = ['step0PlotOriginalFiles = True','\n',
                        'step1TvReg = False','\n',
                        'step2Output = False','\n',
                        'pA.process(step0PlotOriginalFiles, step1TvReg, step2Output)']
    
    code_part_TAB[2] = ['step0PlotOriginalFiles = False','\n',
                      'step1TvReg = True','\n',
                      'step2Output = False','\n',
                      'pA.process(step0PlotOriginalFiles, step1TvReg, step2Output)']
    
    code_part_TAB[3] = ['pA.denoisePlotAll(asp=0.01,onlyPlot=True)','\n',
                       'aIndex=','', '\n', 
                       'pA.defineAlcoholTime(more=20,aIndex=aIndex)']
    
    code_part_TAB[4] = ['step0PlotOriginalFiles = False','\n',
                        'step1TvReg = False','\n',
                        'step2Output = True','\n',
                        'pA.process(step0PlotOriginalFiles, step1TvReg, step2Output)']
    
    code_part_TAB[5] = ['pA.plotAllFiles(0.03)']
    
    code_part_TAB[6] = ['import pickle','\n',
                       'f0 = open(os.path.join(pA.path,pA.subfolderName, \'TA_\'+pA.saveName+\'_df_relative_timegrid.pickle\'), \'rb\')','\n',
                       'df_ts = pickle.load(f0)','\n',
                       'df_ts.checkLength()']
    
    
    code_PRE = [['import sys,os','\n',
               'python_class_path=\'','\'' ,'\n',
               'sys.path.insert(1, os.path.join(sys.path[0], python_class_path))','\n',
                'from Preprocessing import PreprocessingTabFile','\n',
                'current_path = \'','\'','\n',
                'prev_path = \'','\'','\n',
                'filter_file_path = \'','\'','\n',
                'selfFilter = ','\n',
                'pre = PreprocessingTabFile(current_path,prev_path,filterFile=filter_file_path,autoFilter = selfFilter,report_folder_name=\'report\',showProgress = True)','\n',
                'pre.run()']]
    
    code_MD=[['%matplotlib inline','\n',
              'import sys,os','\n',
              'python_class_path=\'','\'' ,'\n',
              'sys.path.insert(1, os.path.join(sys.path[0], python_class_path))','\n',
              'from Mortality_distribution import MortalityDistribution','\n',
              'folder_path=\'','\'','\n',
              'ref_df_path = \'','\'','\n',
              'M_Distribution = MortalityDistribution(folder_path,ref_df_path)','\n',
              'M_Distribution.run()']]
    
    
    code_SELECTOR=[['import sys,os','\n',
               'python_class_path=\'','\'' ,'\n',
               'sys.path.insert(1, os.path.join(sys.path[0], python_class_path))','\n',
                'from Selector_GG_GGM import GG_GGM_Selector','\n',
                'data_folder = \'','\'','\n',
                'out_folder = \'','\'','\n'
                'gg_ggm_selector = GG_GGM_Selector(data_folder,out_folder)','\n',
                'gg_ggm_selector.run()']]
    
    code_CR=[['%matplotlib inline','\n',
            'import sys,os','\n',
            'python_class_path=\'','\'' ,'\n',
            'sys.path.insert(1, os.path.join(sys.path[0], python_class_path))','\n',
            'from Check_Replication import CheckReplication','\n',
            'folderPath=\'','\'','\n',
            'path_ref_folder =\'','\'','\n',
            'cutoff = ','','\n',
            'CR = CheckReplication(folderPath, refFolderPath=path_ref_folder, cutoff=cutoff)','\n',
            'CR.run()']]

    code_DEN=[['%matplotlib inline','\n',
            'import sys,os','\n',
            'python_class_path=\'','\'' ,'\n',
            'sys.path.insert(1, os.path.join(sys.path[0], python_class_path))','\n',
            'from Dendrogram_replicates_check import Dendrogram_Replicates','\n',
            'data_folder = \'','\'','\n',
            'file_replicates_path = \'','\'','\n'
            'select_cluster0=','\n'
            'dr = Dendrogram_Replicates(folder_path=data_folder, file_replicates_path=file_replicates_path,select_cluster0=select_cluster0)','\n',
            'dr.run()']]
    
    
    def __init__(self,data_folder_path='',filterPath='',filterName = '', outputPath = '',showProgress = True):
        """
        Initialize the Notebook Generator,
        :param str data_folder_path: folder path containing .tab files on subfolders, if no input,
                              notebook will be created based on filterPath and strainName,
                              better to precise this path for checking the lack files of filterpath.
        :param str filterPath: predefined list files, same strain seperated by return label,
                               will be created autmotically,if not defined, the notebook will be 
                               created based on data_path (only absolute path files)
        :param str filterName: only speciefied strain name related notebook will be generated,
                               will work with data_path or filterPath
        
        """
        self.data_folder_path = data_folder_path;
        self.filterFilePath = filterPath
        self.filterName = filterName
        if len(outputPath): 
            self.outputPath = outputPath
        else:
            self.outputPath = '.'
        
        self.showProgress = showProgress
        self.callPIP = False;
        #parse the input
        if len(self.data_folder_path) == 0 and len(filterPath)>0:
            self.subFolderPaths = IPythonNotebookGenerator.getAllSubFolderPathsFromFilterFile(filterPath)
        elif len(self.data_folder_path) > 0 and len(filterPath) == 0:
            self.callPIP = True
            self.subFolderPaths = IPythonNotebookGenerator.getAllSubFolderPaths(data_folder_path)
        else:
            print "DatafolderPath and filterPath cann't be both defined, check your input args"
            return
        if len(filterName)>0:
            self.subFolderPaths = IPythonNotebookGenerator.filterList(self.subFolderPaths,filterName)
            
        if len(self.subFolderPaths) == 0:
            return

        if showProgress:
            print 'Initialization of Ipython auto generator'
            #print self.subFolderPaths
     
    @staticmethod
    def removeFile(folderPath, startWith=''):
        '''
        remove all files startwith input str in folder and subfolders
        '''
        if len(startWith)==0:
            return
        import os
        for parent, dirnames, filenames in os.walk(folderPath):
            for fn in filenames:
                if fn.startswith(startWith):
                    print 'remove file : ', os.path.join(parent, fn)
                    os.remove(os.path.join(parent, fn))
                
    def run_TAB_PIP(self,python_class_path,alcohol_time="-1", parent_path=''):
        """
        @param str parent_path: optional, used when filter file (not self.callPIP) has relative path, 
        """
        if len(self.subFolderPaths) == 0:
            return
        import nbformat as nbf
        import os
        import copy
       
        if self.showProgress:
            print "---Ipython code generation process---"
        
        for l in self.subFolderPaths:
            fileNames = l[1]
            dataFolderPath = l[0]
            partial = 1
            for fn in fileNames:
                if self.callPIP:
                    code_Group = copy.deepcopy(IPythonNotebookGenerator.code_part_PIP)
                    text = IPythonNotebookGenerator.text_PIP
                else:
                    code_Group = copy.deepcopy(IPythonNotebookGenerator.code_part_TAB)
                    IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'fFn',self.filterFilePath)
                    IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'partial',partial)
                    partial = partial + 1
                    if not os.path.isabs(l[0]): 
                        dataFolderPath = os.path.join(parent_path,l[0])
                    text = IPythonNotebookGenerator.text_TAB
                
                IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'pth',dataFolderPath)
                IPythonNotebookGenerator.removeFile(dataFolderPath, "TA_")
                IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'python_class_path',python_class_path)
                IPythonNotebookGenerator.insertCode_by_name(code_Group[3],'aIndex',alcohol_time)
                #print code_Group
                nb = IPythonNotebookGenerator.codeGroupToNb(code_Group, text)
                #print nb
                if self.showProgress:
                    print os.path.join(self.outputPath,fn+'.ipynb')
                f = IPythonNotebookGenerator.createFile(fn+'.ipynb',self.outputPath)
                nbf.write(nb, f) 
                f.close()
                #print nb   
    
    @staticmethod
    def run_PRE(python_class_path,folder_path='',current_path='',prev_path='',filterFile='',selfFilter = 'True'):
        import copy
        import nbformat as nbf
        import os
        text = IPythonNotebookGenerator.text_PRE
        code_Group = copy.deepcopy(IPythonNotebookGenerator.code_PRE)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'python_class_path',python_class_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'current_path',current_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'prev_path',prev_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'filter_file_path',filterFile)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'selfFilter',selfFilter)
        nb = IPythonNotebookGenerator.codeGroupToNb(code_Group, text)
        f = IPythonNotebookGenerator.createFile('Preprocessing.ipynb',folder_path)
        nbf.write(nb, f) 
        f.close()
        
    @staticmethod
    def run_MD(python_class_path,folder_path='',data_folder_path='',ref_df_path=''):
        import copy
        import nbformat as nbf
        import os
        IPythonNotebookGenerator.removeFile(data_folder_path, "MD-")
        text = IPythonNotebookGenerator.text_MD
        code_Group = copy.deepcopy(IPythonNotebookGenerator.code_MD)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'python_class_path',python_class_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'folder_path',data_folder_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'ref_df_path',ref_df_path)
        nb = IPythonNotebookGenerator.codeGroupToNb(code_Group, text)
        f = IPythonNotebookGenerator.createFile('Mortality_distribution.ipynb',folder_path)
        nbf.write(nb, f) 
        f.close()
        
    @staticmethod
    def run_SELECTOR(python_class_path,folder_path='',data_folder_path=''):
        import copy
        import nbformat as nbf
        import os
        text = IPythonNotebookGenerator.text_SELECTOR
        code_Group = copy.deepcopy(IPythonNotebookGenerator.code_SELECTOR)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'python_class_path',python_class_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'data_folder',data_folder_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'out_folder',os.path.join(data_folder_path,'GG_GGM_TABLE'))
        nb = IPythonNotebookGenerator.codeGroupToNb(code_Group, text)
        f = IPythonNotebookGenerator.createFile('GG_GGM_Selector.ipynb',folder_path)
        nbf.write(nb, f) 
        f.close()
        
    @staticmethod
    def run_DEN(python_class_path,folder_path='',data_folder_path='',filePath='',selectCluster0='True'):
        import copy
        import nbformat as nbf
        import os
        text = IPythonNotebookGenerator.text_DEN
        code_Group = copy.deepcopy(IPythonNotebookGenerator.code_DEN)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'python_class_path',python_class_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'data_folder',data_folder_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'file_replicates_path',filePath)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'select_cluster0',selectCluster0)
        nb = IPythonNotebookGenerator.codeGroupToNb(code_Group, text)
        f = IPythonNotebookGenerator.createFile('Dendrogram_check_replicates.ipynb',folder_path)
        nbf.write(nb, f) 
        f.close()
        
    @staticmethod
    def run_CR(python_class_path,folder_path='',folderPath='', path_ref_folder='', cutoff=''):
        import copy
        import nbformat as nbf
        import os
        from Check_Replication import CheckReplication
        IPythonNotebookGenerator.removeFile(folderPath, CheckReplication.pre)
        text = IPythonNotebookGenerator.text_CR
        code_Group = copy.deepcopy(IPythonNotebookGenerator.code_CR)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'python_class_path',python_class_path)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'folderPath',folderPath)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'path_ref_folder',path_ref_folder)
        IPythonNotebookGenerator.insertCode_by_name(code_Group[0],'cutoff = ',cutoff)
        nb = IPythonNotebookGenerator.codeGroupToNb(code_Group, text)
        f = IPythonNotebookGenerator.createFile('Check_replication.ipynb',folder_path)
        nbf.write(nb, f) 
        f.close()
            
    @staticmethod
    def execute(folder_path, filterNames=[]):
        """
        execute ipython automatically
        :param str folder_path: path contains .ipynb files
        :param str or list[str] filterNames: only .ipynb contains names will be executed, or ALL if no input
        """
        print "---Execute ipython code---"
        if isinstance(filterNames, str):
            fList = []
            fList.append(filterNames)
            filterNames = fList
            
        paths = IPythonNotebookGenerator.listFilesFromFolder(folder_path)
        IPythonNotebookGenerator.callPython(paths,filterNames)

    @staticmethod
    def createFile(fname, directory=''):
        import os
        if len(directory) and not os.path.exists(directory):
            os.makedirs(directory)
            print directory
        return open(os.path.join(directory,fname), 'w')
        
    @staticmethod
    def getAllSubFolderPaths(folderPath, ext='.tab'):
        """
        Return all child folders (only subfolder) included current folder which contain ext files 
        :param str folderPath: folder path
        :return tuple subfolders' absolute path and name
        """
        import path
        paths=[];
        import os
        if not os.path.isdir(folderPath):
            print "folder not existed"
            return paths
        all_subfolders = next(os.walk(folderPath))[1]
        if any(fname.endswith(ext) for fname in os.listdir(folderPath)):
            #remove the last non alphanumeric char if user input / at the end
            if not folderPath[-1].isalnum():
                folderPath = folderPath[:-1]
            paths.append((folderPath,[os.path.split(folderPath)[1]]))
        if len(all_subfolders):
            for sub in all_subfolders:
                subPath = os.path.join(folderPath,sub)
                if any(fname.endswith(ext) for fname in os.listdir(subPath)):
                    paths.append((subPath,[os.path.split(subPath)[1]]))  
        return paths

    
    @staticmethod
    def getAllSubFolderPathsFromFilterFile(filterFile,data_folderPath=''):
        '''
        get all folder names from filterFile, if files are separated by return label,
        multiple folder names will be generated,
        :param filterFile: contains file paths list
        :param data_folderPath: optional, if defined, check the file existency inside the filterFile
        :return: tuple subfolder path (can be relative) with folderGroup Names
        '''
        import glob,os
        existing_fns=[]
        if len(data_folderPath.rstrip()) > 0:
            existing_fns = glob.glob(os.path.join(data_folderPath, saveName,'*.tab'))
        # get all the files in the path corresponding the filterFile.txt

        listFiles = [];
        with open(filterFile, 'rU') as rf:
            for line in rf:
                listFiles.append(line.rstrip());
        listFiles.append('')
        saveName = '';
        subFoldersPath = [];

        for i in range(len(listFiles)):
            line = listFiles[i]  
            if len(line) > 0:
                if saveName == os.path.split(line)[0]:
                    continue;
                else:
                    saveName = os.path.split(line)[0]
                previous = ''
                count = 0
                folderGroup = [];
                for j in range(i,len(listFiles)):
                    line_fn = listFiles[j]
                    if len(line_fn) > 0:
                        if os.path.split(line_fn)[0]==saveName:
                            if len(existing_fns) > 0 and (line_fn not in existing_fns):
                                continue;
                            if previous == '':
                                count=count+1
                                folderGroup.append(os.path.split(saveName)[1]+'-'+str(count))
                    previous = line_fn.rstrip();
                if len(folderGroup) > 0:
                    subFoldersPath.append((saveName, folderGroup))
        return subFoldersPath
    
    @staticmethod
    def filterList(folderPathList,filterName):
        """
        :param:folderPathList: tuple ex:('appY', ['appY-1', 'appY-2'])
        optional, only return list containing filterName
        """
        lFiltered = []
        for l in folderPathList:
            if filterName in l[0]:
                lFiltered.append(l);
        return lFiltered;
    
                
    @staticmethod
    def generateNotebook():
        import nbformat as nbf
        nb = nbf.v4.new_notebook();
        return nb

    
    @staticmethod
    def appendCode(nb,code='',text=''):
        import nbformat as nbf
        if nb == None:
            nb = IPythonNotebookGenerator.generateNotebook()
            
        if len(text) > 0:
            nb['cells'].append(nbf.v4.new_markdown_cell(text))
        if len(code) > 0:
            nb['cells'].append(nbf.v4.new_code_cell(code))    
            
        return nb
    
    @staticmethod
    def appendCodeList(codeList,nb=None,text_head=''):
        """
        codeGroup should be firstly changed to str list by codeGroupToList()
        """
        if len(text_head): #run once
            nb = IPythonNotebookGenerator.appendCode(nb,text=text_head)
        for l in codeList:
            nb = IPythonNotebookGenerator.appendCode(nb,code=l)
        return nb
    
    @staticmethod
    def codeGroupToList(codeGroup):
        """
        code_part_PIP[0-end] should be firstly changed to str then added to the output list
        """
        codeList = [];
        for g in codeGroup:
            if g is None:
                continue
            codeList.append(IPythonNotebookGenerator.listToString(g))
        return codeList
    
    @staticmethod
    def codeGroupToNb(codeGroup,text=''):
        """
            high level api, call by user
        """ 
        codeList = IPythonNotebookGenerator.codeGroupToList(codeGroup)
        nb = IPythonNotebookGenerator.appendCodeList(codeList,text_head = text)
        return nb
        
                
            
    @staticmethod
    def writeNbToFile(nb, filePath):
        import nbformat as nbf
        with open(filePath, 'w') as f:
            nbf.write(nb, f)       
            
            
    @staticmethod
    def insertCode_by_pos(code_list,position,content):
        import warnings
        warnings.warn("deprecated", DeprecationWarning)
        code_list.insert(position,content);
    
    @staticmethod
    def insertCode_by_name(code_list,name,content):
        """
        only the first found name will be found, content will be inserted after the name
        """
        foundName = False
        for i in range(len(code_list)):
            l = code_list[i]
            if name in str(l):
                code_list.insert(i+1,str(content))
                foundName = True
                break;
                
        try:
            if not foundName:
                raise ValueError('Name not found in code list...')
        except ValueError as e:
            print e;
        
    
    @staticmethod
    def listToString(code_list):
        code = '';
        for l in code_list:
            code = code + str(l);
        return code
    
    
    @staticmethod
    def runIpython(path):
        import subprocess
        p = subprocess.check_call(['jupyter' ,'nbconvert', '--execute', '--inplace','--ExecutePreprocessor.timeout=-1',path])
        if p:
            print "Raise Call Error : ", path

                
    @staticmethod
    def filterPath(paths,filterList=[]):
        path_output = []
        if len(filterList) == 0:
            return paths
        for p in paths:
            found = False
            for f in filterList:
                if f in p:
                    found = True
                    path_output.append(p)
                    break
                    
        return path_output
    
    @staticmethod
    def callPython(paths,filterList=[]):
        path_filtered = IPythonNotebookGenerator.filterPath(paths, filterList)
        if len(path_filtered)==0:
            print "No ipython files found, please change the folder path or filter name"
        for p in path_filtered:
            print "--Execute-- ",p
            IPythonNotebookGenerator.runIpython(p)
            
    
    @staticmethod
    def listFilesFromFolder(folderPath,ext='.ipynb'):
        import os
        import glob
        return glob.glob(os.path.join(folderPath, '*'+ext))
                    
    
