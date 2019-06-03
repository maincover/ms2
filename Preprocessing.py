# coding: utf-8

# In[ ]:

class PreprocessingTabFile(object):
    """
    Duplication check independently or compared with existed database (Microfluidic database)
    To generate HTML documentation for this module issue the
    command:
        pydoc -w foo
    """

    def __init__(self, currDBPath='',prevDBpath='',filterFile = '',report_folder_name='',autoFilter = True,showProgress = False):     
        """
        :param currDBPath: folder path of the test database tab files in sub folder,
        :param prevDBpath: (OPT) folder path to existed database, to check with the current tabs
        :param filterFile: (OPT) pre existed .txt files' list with separation
        :param reportFolderPath: (OPT) The folder path where the report and result(pickle) will be saved to,
        :param autoFilter: (OPT) To select the ones with filter name to check duplications
        :param showProgress:(OPT) get more information of the process
        """
        if currDBPath=='':
            self.db_path = PreprocessingTabFile.getCurrentPath();
        else:
            self.db_path = currDBPath;
            
        if len(report_folder_name) == 0:
            self.report_fpath = os.path.join(self.db_path, "report")
        else:
            import os
            self.report_fpath = os.path.join(self.db_path, report_folder_name)
            if not os.path.exists(self.report_fpath):
                os.makedirs(self.report_fpath)
        
        self.filterFile = filterFile                
        self.prev_db_path = prevDBpath
        self.all_tabs = [];
        self.head_dups = [];
        self.md5_uniq = {};
        self.md5_dups = [];
        self.file_md5s = {};
        self.autoFilter = autoFilter;
        self.showProgress = showProgress;
        self.all_tabs_filter = []
    
    def run(self):
        #get all tab file in the db_path
        self.all_tabs = PreprocessingTabFile.updateTabFileByPath(self.all_tabs,self.db_path)
        if self.autoFilter:
            import copy
            self.all_tabs_filter = copy.deepcopy(self.all_tabs)
        #get unique file and dup file by md5
        self.md5_uniq,self.md5_dups,self.file_md5s = PreprocessingTabFile.checkMd5DupList(self.all_tabs)
        #get head duplication file
        PreprocessingTabFile.checkHeaderDubsToFile(self.all_tabs,self.report_fpath)
        if len(self.prev_db_path):
            #update all_tabs and md5 uniq, md5 dup files
            self.all_tabs = PreprocessingTabFile.updateTabFileByPath(self.all_tabs,self.prev_db_path)
            self.md5_uniq,self.md5_dups,self.file_md5s = PreprocessingTabFile.checkMd5DupList(self.all_tabs)       
        filteredFiles = PreprocessingTabFile.filterFiles(self.file_md5s.keys(), self.all_tabs_filter)
        if self.showProgress:
            print 'Total files length : ', len(self.file_md5s.keys()),' Filtered files length (after MD5 check) : ', len(filteredFiles)
            
        PreprocessingTabFile.checkContentDubsToFile(self.file_md5s.keys(),self.report_fpath,self.all_tabs_filter,self.showProgress)
        PreprocessingTabFile.filterDupFiles(self.report_fpath)
        self.report_file_path = PreprocessingTabFile.generateNoDupFiles(self.report_fpath);
        if len(self.filterFile)>0:
            PreprocessingTabFile.updateFinalFilesList(self.report_file_path,self.filterFile)
    
    @staticmethod
    def getCurrentPath():
        import os
        return os.getcwd()
    
    @staticmethod
    def updateTabFileByPath(all_tabs=[],path=''):
        """The tab located under the sub-folder of the path
        """
        import os,csv
        def checkDataFolder(fd):
            for fn in os.listdir(fd):
                if fn[-4:] == '.tab':
                    with open(os.path.join(fd, fn),'Ur') as f:
                        reader = csv.reader(f, dialect='excel',delimiter='\t')
                        header = reader.next()
                        for line in reader:
                            if len(line) != len(header):
                                sys.stdout.write('\t'.join([os.path.join(dir, fn),str(reader.line_num),str(len(header)-len(line))])+'\n')
                                break
                    all_tabs.append(os.path.abspath(os.path.join(fd, fn)))
                elif os.path.isdir(os.path.join(fd,fn)):
                    checkDataFolder(os.path.join(fd,fn))
             

        if path=='':
            return all_tabs
        else:#should include itself as well
            checkDataFolder(path)
        
        for fd in os.listdir(path):
            fd_full_path = os.path.join(path,fd)
            if os.path.isdir(fd_full_path):
                checkDataFolder(fd_full_path)
        
                
        return all_tabs

    @staticmethod
    def sort_nicely(l):
        def tryint(s):
            try:
                return int(s)
            except:
                return s
        def alphanum_key(s):
            import re
            """ Turn a string into a list of string and number chunks.
                "z23a" -> ["z", 23, "a"]
            """
            return [tryint(c) for c in re.split('([0-9]+)', s)]
        
        """ Sort the given list in the way that humans expect.
        """
        l.sort(key=alphanum_key) 


   
    
    @staticmethod
    def checkMd5DupList(all_tabs=[],uniq=None, dups=None):
        """return all uniq md5 files to uniq, all duplications in dups, and all files in file_md5s"""
        import hashlib
        file_md5s = {}
        if uniq is None:
            uniq={}
        if dups is None:
            dups=[]
        for fn in all_tabs:
            hash_md5 = hashlib.md5()
            with open(fn ,'Ur') as f:
                for line in f.readlines()[1:]:
                    hash_md5.update(line)
            file_md5s[fn] = hash_md5.hexdigest()
        for fn in file_md5s.keys():
            if file_md5s[fn] not in uniq.keys():
                uniq[file_md5s[fn]] = fn
            else:
                dups.append(tuple(sorted((fn, uniq[file_md5s[fn]]), reverse=True)) )

        sorted(dups)
        return uniq,dups,file_md5s
    
    @staticmethod
    def checkHeaderDubsToFile(all_tabs=[],report_folderPath='.'):
        """check header duplication for all_tabs files
        write the result to header_duplications.txt"""
        import hashlib
        import os
        header_md5s = {}
        reportf = open(os.path.join(report_folderPath,'header_duplications.txt'),'w',0)
        for fn in all_tabs:
            if fn[-4:]=='.tab':
                hash_md5 = hashlib.md5()
                with open(fn ,'Ur') as f:
                    hash_md5.update(f.readline())
                header_md5s[fn] = hash_md5.hexdigest()


        fn_by_h = {}
        dups_h_count = {}
        for fn in header_md5s.keys():
            if header_md5s[fn] not in fn_by_h.keys():
                fn_by_h[header_md5s[fn]] = [fn]
            else:
                fn_by_h[header_md5s[fn]].append(fn)
                dups_h_count[header_md5s[fn]] = len(fn_by_h[header_md5s[fn]])

        uniq_h = fn_by_h;

        for i in range(len(dups_h_count)):
            dups_md5 = sorted(dups_h_count.keys())[i]
            fns = uniq_h[dups_md5]
            reportf.write('\t'.join(fns)+'\n\n')
            #print '\t'.join(fns)+'\n'
            
            
    @staticmethod
    def filterFiles(files, files_filter):
        if len(files_filter)==0:
            return files
        filteredFiles = [None]*len(files);
        for f in files:
            if f in files_filter:
                filteredFiles.append(f);
        filteredFiles = filter(None, filteredFiles)
        return filteredFiles;
    
        
    @staticmethod
    def checkContentDubsToFile(files,report_folderPath='.',files_filter=[],showProgress = False,option=False):
        import os,csv
        """check content of files including those who have already header duplication
        write result to line_duplication_ratios.txt and no_duplications_files.txt"""
        
        def num(s):
            try:
                return str(int(s))
            except ValueError:
                try: 
                    return str(int(float(s)))
                except ValueError:
                    raise ValueError(s)

        import difflib
        import pandas as pd
        all_ratios = pd.DataFrame(index=files,columns=files)
        SeqMat = difflib.SequenceMatcher(lambda x: 0 )
        reportf = open(os.path.join(report_folderPath,'line_duplication_ratios.txt'),'w',buffering=1)
        checkedf = open(os.path.join(report_folderPath,'no_duplications_files.txt'),'w',buffering=1)
        checked = []
        data_rs = {}
        dups = {}
        
        m =1
        
        files = sorted(files)
        for fn in files:
            if showProgress and m%50 == 0:
                print "check line duplication(open file) : ",m,'/',len(files)
            with open(fn ,'Ur') as f:
                r = csv.reader(f,dialect='excel',delimiter='\t')
                header = r.next()
                data_rs[fn] = ['\t'.join([num(nb) for nb in row[4:]]) for row in r]
                f.close()
                del r
            dups[fn]=0
            m = m+1;
        

        for i in range(len(files)):
            if showProgress and i%50 == 0:
                print "check line duplication: ", i,'/',len(files)
            fn1 = files[i]  
            if len(files_filter)> 0 and fn1 not in files_filter:
                continue;
            data_rs1 = data_rs[fn1]
            for j in range(len(files)):                #range(i+1,len(files)):                
                fn2 = files[j]                
                if fn2==fn1:
                    print "exclude fn(itself) : ", fn1
                    continue
                data_rs2 = data_rs[fn2]
                SeqMat.set_seq1(data_rs1)
                SeqMat.set_seq2(data_rs2)
                all_ratios.loc[fn1,fn2] = SeqMat.real_quick_ratio()
                if all_ratios.loc[fn1,fn2] > 0:
                    all_ratios.loc[fn1,fn2] = SeqMat.quick_ratio()
                if all_ratios.loc[fn1,fn2] > 0:
                    all_ratios.loc[fn1,fn2] = SeqMat.ratio()
                if all_ratios.loc[fn1,fn2] > 0:
                    if showProgress:
                        print '\t'.join([fn1, fn2, str(all_ratios.loc[fn1,fn2])])+'\n';
                    reportf.write('\t'.join([fn1, fn2, str(all_ratios.loc[fn1,fn2])])+'\n')
                    print fn1,fn2,all_ratios.loc[fn1,fn2]
                    dups[fn1] += 1
                    dups[fn2] += 1
            if dups[fn1] == 0:
                checked.append(fn1)
                checkedf.write(fn1+'\n')
        all_ratios.to_pickle(os.path.join(report_folderPath,'line_duplication_dataframe.pickle'))
        
        reportf = open(os.path.join(report_folderPath,'line_duplication_count.txt'),'w',buffering=1)
        for fn in sorted(dups.keys()):
            if dups[fn]>0:
                reportf.write(fn+'\t'+str(dups[fn])+'\n')
                
        if option:
            import itertools
            for (fn1, fn2) in itertools.product(files,files):
                if fn1 != fn2:
                    if all_ratios.loc[fn1,fn2] != 0:
                        with open(fn1 ,'Ur') as f1:
                            lines1 = f1.readlines()[1:]
                        with open(fn2, 'Ur') as f2:
                            lines2 = f2.readlines()[1:]
                        linesSM = difflib.SequenceMatcher(lambda x: 0,lines1, lines2 )
                        all_ratios.loc[fn1,fn2] = linesSM.ratio()

            ldups_pct = []
            files = sorted(files)
            for i in range(len(files)):
                if len(files_filter)>0 and files[i] not in files_filter:
                    continue;
                for j in range(i+1,len(files)):
                    if float(all_ratios.loc[files[i],files[j]]) > 0:
                        ldups_pct.append((files[i], files[j], all_ratios.loc[files[i],files[j]])) 

            reportf = open(os.path.join(report_folderPath,'data_duplication_ratios.txt'),'w')
            for (fn1,fn2,pct) in sorted(ldups_pct, key=lambda tpl: tpl[0]):
                reportf.write('\t'.join([fn1, fn2, str(pct)])+'\n')
            all_ratios.to_pickle(os.path.join(report_pth,'data_duplication_dataframe.pickle'))

    @staticmethod
    def filterDupFiles(report_folderPath):
        import os, csv, sys, hashlib
        def getTreeFiles(baseFile, ldups, listFiles):
            listFiles.append(baseFile)
            for row in ldups:
                while row[0]==baseFile and (row[1] not in listFiles):
                    getTreeFiles(row[1], ldups, listFiles);
            return listFiles;

        def check2_lastFromSameFolder(listName):
            if len(listName) < 2:
                return True;
            firstName = listName[0].split("/")[:-1]
            for index, row in enumerate(listName):
                if firstName == row.split("/")[:-1]:
                    continue;
                else:
                    return False;
            return True;

        ldups = []
        import random;

        with open(os.path.join(report_folderPath,'line_duplication_ratios.txt'), 'rU') as readf:
            csvreader = csv.reader(readf, delimiter='\t')
            for row in csvreader:
                ldups.append((row[0],row[1],row[2]))

        uniq = {}
        name = "";
        fName = "";
        excludeFile =[]
        includeFile = []
        for index, row in enumerate(ldups):
            if row[0] in excludeFile:
                continue;
            if row[0] in includeFile:
                continue;
            if(index > len(ldups)):
                continue;
            #print "processing index : ", index
            fileName = row[0]
            listFiles = [];
            listFiles = getTreeFiles(fileName, ldups, listFiles)
            if check2_lastFromSameFolder(listFiles):
                randomFile = random.choice(listFiles)
                includeFile.append(randomFile);
                for r1 in listFiles:
                    if r1==randomFile:
                        continue;
                    excludeFile.append(r1);
            else:
                for r1 in listFiles:
                    excludeFile.append(r1);
        
        reportf = open(os.path.join(report_folderPath,'random_selected_files.txt'),'w')
        for fn1 in sorted(includeFile, key=lambda tpl: tpl[0]):
            reportf.write(fn1+'\n')


    @staticmethod
    def generateNoDupFiles(report_folderPath):
        import os,csv
        lFilter=[]
        if os.path.isfile(os.path.join(report_folderPath,'random_selected_files.txt')): 
            with open(os.path.join(report_folderPath,'random_selected_files.txt'), 'rU') as readf:
                csvreader = csv.reader(readf, delimiter='\t')
                for row in csvreader:
                    lFilter.append(row[0])
        else:
            print "random_selected_files.txt not existed"

        lNoDup=[]
        with open(os.path.join(report_folderPath,'no_duplications_files.txt'), 'rU') as readf:
            csvreader = csv.reader(readf, delimiter='\t')
            for row in csvreader:
                lNoDup.append(row[0])

        lHeader = []
        with open(os.path.join(report_folderPath,'header_duplications.txt'), 'rU') as readf:
            csvreader = csv.reader(readf, delimiter='\t')
            for row in csvreader:
                lHeader.append(row)

        newHeader = []
        for row in lHeader:
            if len(row) > 0:
                for t in row:
                    r = t.split('\t');
                    for r0 in r:
                        newHeader.append(r0)

        lFinal=[]
        for row in lNoDup:
            if row not in newHeader:
                lFinal.append(row);

        for row in lFilter:
            if row not in newHeader:
                if row not in lFinal:
                    lFinal.append(row);

        PreprocessingTabFile.sort_nicely(lFinal)
        reportf = open(os.path.join(report_folderPath,'finalFiles.txt'),'w')
        for fn1 in sorted(lFinal, key=lambda tpl: tpl[0]):
            reportf.write(fn1+'\n')
        return os.path.join(report_folderPath,'finalFiles.txt')
    
    @staticmethod
    def updateFinalFilesList(r0,r1):
        """1, independant call:
        r0: generated report file path
        r1: predefined filter file path
        
        if prev report updated, pre-clean the finalFiles of current report folder first
        if curr report updated, no need to clean,
        if both updated, pay attention to pre-clean carefully the prev report folder 
    
        2,if curr report updated
        Call the two static method,
        PreprocessingTabFile.filterDupFiles(report_fpath)
        PreprocessingTabFile.generateNoDupFiles(report_fpath);
        """
        import os, csv
        l1 = [];
        f1 = open(r0,'rU')
        try:
            with f1 as readf:
                csvreader = csv.reader(readf, delimiter='\t')
                for row in csvreader:
                    if len(row)>0:
                        l1.append(row[0])
                    else:
                        l1.append('\n')
        finally:
            f1.close()
        l2 = [];
        f2 = open(r1,'rU')
        try:
            with f2 as readf:
                csvreader = csv.reader(readf, delimiter='\t')
                for row in csvreader:
                    if len(row)>0:
                        l2.append(row[0])
                    else: 
                        l2.append('\n')
        finally:
            f2.close()
        if not PreprocessingTabFile.sublist(l2,l1):
            print "length of finalFiles before update : ",len(l1)
            l1.extend(l2)
            print "length finalFiles after update: ", len(l1)
            f3 = open(r0,'w')
            try:
                for fn1 in l1:
                    if fn1 == '\n':
                        f3.write(fn1)
                    else:
                        f3.write(fn1+'\n');
            finally:
                f3.close()
        else:
            print "files already found in finalFiles"
    
    @staticmethod
    def sublist(l2,l1):
        """
        check if l2 belonged to l1
        """
        arrayList = lambda x:any(l2 == x[offset:offset+len(l2)] for offset in range(len(x)))
        return arrayList(l1)