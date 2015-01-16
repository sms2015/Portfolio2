import pandas as pd
import tarfile
import os
import numpy as np
import itertools
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *
import shutil

def intro():
    #print introductory information to the console
    print ('About this Program\n')

    print ('This program uses data from the NOAA Global Historical Climatology '
           'Network (GCHN) to analyze historical snowfall and snowpack at North '
           'American weather stations and looks for correlations between these '
           'measures and historical El Nino/La Nina episodes between the years '
           '1950 and 2010.  More information about the data used in this '
           'program can be found at: http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt '
           'and more information about this program can be found in the readme.txt '
           'file and the project report file "IS 602 Final Project.docx", both of '
           'which are contained in the IS602_FinalProject.tar file\n')

def get_tar_file(filename,action):
    #create a TarFile object
    t = tarfile.open(filename, action)

    t.extractall()  #extracts all files from the tar files
    
    t.close

def get_station_file():
    #open the ghcnd-stations file which contains information about weather
    #stations, turn this file into a data frame
    stationfile ='ghcnd-stations.txt'

    #data file has fixed width columns: http://pandas.pydata.org/pandas-docs/stable/io.html
    widths = [11,9,10,7,3,31,4,4,6]

    df = pd.read_fwf(stationfile,widths=widths,header=None, index_col=0)
    df.columns =['Latitude','Longitude','Elevation','State','StationName',
                'GSNFlag','HCNFlag','WMOId']

    df.index.name ='StationId'

    return df

def NoAmStation():  #Narrow data set to just North American stations
    df=get_station_file()
    NoAmStations = df[pd.notnull(df['State'])]
    return NoAmStations

    #might want to show elevation information in output
    #might want to consider a googlemaps API for visualization

def get_inventory_file():
    #import inventory file return a data frame of the relevant data
    inventoryFile = 'ghcnd-inventory.txt'

    widthsInv = [11,9,10,5,5,5]

    dfI = pd.read_fwf(inventoryFile,widths=widthsInv,header=None, index_col=0)
    dfI.columns = ['Latitude','Longitude','Element','Firstyear','Lastyear']

    return dfI

def NoAM5010Station():
    #Identify the North American stations in the inventory file that have
    #data for 1950 and 2010. Data is not available in the inventory file for
    #interim years, just first and last year by element type.  A file must
    #be extracted and opened in a later step to check for interim years.
    dfI = get_inventory_file()
    df = get_station_file()
    NoAmStations = NoAmStation()
    
    SNOWfiles = dfI[dfI['Element']=='SNOW']
    SNWDfiles = dfI[dfI['Element']=='SNWD']
    
    #use only data from 1950 through 2010
    Snow5010 = SNOWfiles[(SNOWfiles['Firstyear']<=1950) & (SNOWfiles['Lastyear'] >= 2010)]
    Snwd5010 = SNWDfiles[(SNWDfiles['Firstyear']<=1950) & (SNWDfiles['Lastyear'] >= 2010)]

    Snow5010list = Snow5010.index.tolist()
    Snwd5010list = Snwd5010.index.tolist()

    #find the intersection of stations that have SNOW and SNWD
    list1 = Snow5010list 
    list2 = Snwd5010list
    SNlist = [val for val in list1 if val in list2]

    #find North American stations that have SNOW and SNWD codes
    dfNAsn=NoAmStations[NoAmStations.index.isin(SNlist)]
    #has 4465 rows

    NAsnlist = dfNAsn.index.tolist()

    #return list of North American stations that have SNOW and SNWD codes
    #where the first available record is from the year 1950 or earlier, and
    #where the last record is from 2010 or later.
    
    return NAsnlist

def get_specific_sttn_file(StationFile,filepath): #'ghcnd_all'
    #set up widths
    widths2= [11,4,2,4]+[5,1,1,1]*31
    dfs = pd.read_fwf(os.path.join(filepath, StationFile),
                      widths=widths2,header=None, index_col=0)
    dfs.index.name='StationId'
    return dfs
    
def check_year_month(StationFile, filepath, full_record):
    dfs = get_specific_sttn_file(StationFile,filepath)
    #get year and month data
    dfs1=dfs.iloc[:,[0,1,2]]
    dfs1.columns=['Year','Month','Element']
    dfs2= dfs1[(dfs1['Element']=='SNOW') &
               (dfs1['Month'].isin([1,2,3,4,10,11,12])) &
               dfs1['Year'].isin(range(1950,2011))]
    dfs3= dfs1[(dfs1['Element']=='SNWD') &
               (dfs1['Month'].isin([1,2,3,4,10,11,12])) &
               dfs1['Year'].isin(range(1950,2011))]

    if len(dfs2.index)==full_record and len(dfs3.index) == full_record:
        return True

def full_record_copy(listRows, NAsnlist,full_record,filepath):
    fullRecordlist = []
    for i in range(0,listRows):  #change to listrows
        StationFile = NAsnlist[i]+'.dly'
        if check_year_month(StationFile,filepath,full_record) == True:
            fullRecordlist=fullRecordlist + [NAsnlist[i]]
    return fullRecordlist

def print_sttn_details(StationFile,filename):
    df=get_station_file()
    #print the station number, name, and state
    wsState="".join(list(df['State'][df.index==str(StationFile[:-4])]))
    wsName="".join(list(df['StationName'][df.index==str(StationFile[:-4])]))
    ws=[wsName,wsState]
    return ws

def specific_station_transform(StationFile,filepath):
    #select specific columns from the data file needed for the analysis
    dfs=get_specific_sttn_file(StationFile,filepath)
    col_range = [3]
    for i in range(1,31):
        col_range=col_range+[col_range[i-1]+4]
    col_range = range(0,3)+col_range
    dfns = dfs.iloc[:,col_range]

    #add column names
    dfns.columns=['Year','Month','Element']+range(1,32)
    #range(1,32) is for day of the month 1 through 31

    #select rows with October thru April
    dfOA=dfns[dfns['Month'].isin([1,2,3,4,10,11,12])]

    #slect years from 1950 to 2010
    dfOA50=dfOA[(dfOA['Year']>=1950)&(dfOA['Year']<=2010)]

    return dfOA50

def element(StationFile,filepath,elem_type):
    #select rows where element = elem_type
    #elem_type is either SNOW (snowfall in mm) or SNWD (snowpack in mm))
    dfOA50 = specific_station_transform(StationFile,filepath)
    dfOA50e = dfOA50[dfOA50['Element']==elem_type]
    return dfOA50e

def NA5010_file_copy(NAsnlist,listRows):
    #copy North American files with data for 1950 and 2010 into new
    #directory
    mypath = 'new_ghcnd_all'
    if not os.path.isdir(mypath):
        os.makedirs(mypath)
    for i in range (0,listRows):
        filename = 'ghcnd_all/'+NAsnlist[i]+'.dly'
        shutil.copy2(filename, mypath)

def NA5010_inclusive_file_copy(newNAsnlist,newlistRows):
    mypath = 'new2_ghcnd_all'
    if not os.path.isdir(mypath):
        os.makedirs(mypath)
    for i in range (0,newlistRows):
        filename = 'new_ghcnd_all/'+newNAsnlist[i]+'.dly'
        shutil.copy2(filename, mypath)
    
def fill_missing(StationFile,filepath,elem_type):
    #in dataset -9999 represents missing data, replace -9999 with na
    #see note for this function at the bottom of this function.
    dfOA50e = element(StationFile,filepath, elem_type)

    #October through April by element type station data frame
    dfOA50elem = dfOA50e.replace(-9999, np.nan)
    
    #elemchk= dfOA50elem.count(axis=1,numeric_only=True)
    #need to check if any of these have less than 20 days of observations.

    return dfOA50elem

    #Note: This program ignores missing day values. Accounting for these missing
    #values is beyond the scope of this project. A systematic interpolation
    #process would be needed using Python methods to fill NA data for a true
    #scientific analysis. Or at least a systematic process for when these
    #should or should not be ignored.

def monthly_avg(StationFile,filepath,elem_type,col_name):
    #step 2, Monthly Average: : compute monthly total snowfall for each each
    #month for each year, and average monthly snowpack
    dfOA50elem=fill_missing(StationFile,filepath, elem_type)
    col_months = range(1,32)
    if elem_type=='SNOW':
        dfOA50elem[col_name]=dfOA50elem[col_months].sum(axis=1)
    else:
        if elem_type=='SNWD':
            dfOA50elem[col_name]=dfOA50elem[col_months].mean(axis=1)
    return dfOA50elem

def thirty_year_avg(StationFile,filepath, elem_type,col_name):
    #step 3, 30 Year Monthly Average: compute 30 year (1981-2010) average
    #monthly snowfall and snowpack
    dfOA50elem=monthly_avg(StationFile,filepath, elem_type,col_name)
    year30=range(1981,2011)
    
    #select the rows that correspond to this 30 year range.
    dfOAelem30=dfOA50elem[dfOA50elem['Year'].isin(year30)]

    #find the average for each month October through April
    elem30Avg = dfOAelem30.groupby('Month').mean()[col_name]
    return elem30Avg


def delta(StationFile,filepath, elem_type, col_name, col_name_30,col_name_delta):
    #step 4, Monthly Snowfall vs 30 Year Average: compute the delta between
    #monthly averages (step 2) and 30 year averages (step 3)
    #for years 1950 through 2010

    #create dataframes for snowfall and snowdepth with just year month and
    #the calculated fiels snowsum and snwdavg
    dfOA50elem = monthly_avg(StationFile,filepath,elem_type,col_name)
    elem30Avg = thirty_year_avg(StationFile,filepath,elem_type,col_name)
    
    dfMelem = dfOA50elem.loc[:,['Year','Month',col_name]]

    #create dataframe from SnowAvg30 and SnwdAvg30 series
    elem30Avgdf=pd.DataFrame(elem30Avg)
    elem30Avgdf.columns=[col_name_30]

    #merge 30 yr avg and monthly average dataframes
    elemMergedf=pd.merge(dfMelem,elem30Avgdf,how='outer', left_on='Month',
                    right_index=True)

    #compute the delta vs the 30 year montly average for each month
    elemMergedf[col_name_delta]=elemMergedf[col_name]-elemMergedf[col_name_30]
    return elemMergedf

def get_nino_file():
    #import the Nino Index
    ninoFile = 'Nino_index.csv'
    dfNino = pd.read_csv(ninoFile)

    dfNinoOA = dfNino.loc[:,['Year','1','2','3','4','10','11','12']]
    #Month number headers were read in as text, change to number
    dfNinoOA.columns=['Year',1,2,3,4,10,11,12]
    #print dfNinoOA

    #reshape to similar format as snowfall and snowdepth data
    rangeN=range(1,13)
    dfNinoM=pd.melt(dfNino, id_vars=['Year'], var_name='Month',value_name='NinoIndex')
    #change type from string to integer
    dfNinoM['Month']=dfNinoM['Month'].astype(int)
    return dfNinoM

def elem_nino_merge(StationFile,filepath,elem_type,col_name,col_name_30,
                    col_name_delta):
    #Step 5, Correlation Analysis
    elemMergedf = delta(StationFile,filepath,elem_type, col_name, col_name_30,
                        col_name_delta)
    dfNinoM = get_nino_file()

    #Merge Snow/Snwd and NinoIndex
    elemNinoMergedf=pd.merge(elemMergedf,dfNinoM,how='left',on=['Year','Month'])
    #SnowNinoMergedf.to_csv('snow_nino.csv')
    return elemNinoMergedf

def pearson_corr(StationFile,filepath, elem_type,col_name,col_name_30,
                 col_name_delta):
    #Calculate correlation
    elemNinoMergedf = elem_nino_merge(StationFile,filepath, elem_type,col_name,
                                      col_name_30,col_name_delta)
    elemCorr=stats.pearsonr(elemNinoMergedf[col_name_delta],
                            elemNinoMergedf['NinoIndex'])
    return elemCorr

def scatter(StationFile,filepath,elem_type,col_name,col_name_30,col_name_delta,
            xname,yname,titlename):

    #scatter plot + regression line
    #curve fit approach
    #this will only be shown if there is a significant correlation

    elemNinoMergedf = elem_nino_merge(StationFile,filepath,elem_type,col_name,
                                      col_name_30,col_name_delta)
    elem_arr = elemNinoMergedf[col_name_delta].values
    elemNino_arr = elemNinoMergedf['NinoIndex'].values
    plt.scatter(elemNino_arr, elem_arr)

    #check if p-value is less than or equal to 0.10, if it's not do
    #not draw regression line, since it will create an error
    elemCorr=stats.pearsonr(elemNinoMergedf[col_name_delta],
                            elemNinoMergedf['NinoIndex'])

    if elemCorr[1]<=0.10:
        def func(x, a, b):
            return a * x + b

        popt, pcov = curve_fit(func, elemNino_arr, elem_arr)

        x = array([min(elemNino_arr), max(elemNino_arr)])
        y = popt[1] + popt[0] * x
        plt.plot(x, y, 'r-')

        equation = "Y = X * " + str(round(popt[0],2)) + " + " + str(round(popt[1],2))
    else:
        equation = 'p>0.10, no regression line'

    plt.text(0,max(elem_arr)*0.9,equation)  #place equation at about 85% of the max point
    plt.ylabel(xname)
    plt.xlabel(yname)
    plt.title(titlename)

    plt.tight_layout()

    #source: http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html

def processinput():
    #Ask if the user wants to run through the file processing and extraction
    #process, or just use the pre-processed files
    while True:
        try:
            print ('\nDo you want to extract and process files from the '
                   'original ghcnd_all.tar.gz file, or do you prefer to '
                   'use the pre-processed file?\n')
            print ('Please note that processing and extracting the files '
                   'will require about 30GB of space and could take an '
                   'hour or more to run.\n')
            print 'Select 0 to extract and process files'
            print 'Select 1 to use the pre-processed files'
            input = int(raw_input('Enter 0 or 1 >>> '))
        except ValueError: 
            print 'That\'s not a number!, Try again'
        else:
            if  0 <= input <= 1:
                break
            else:
                print 'Out of range. Try again'
    return input

def siginput():
    #prompt for user to input the level of significance
    while True:
        try:
            print '\nAssign the level of significance'
            input = float(raw_input('You may enter 0.05 or 0.10 >>> '))
        except ValueError: 
            print 'That\'s not a number!, Try again'
        else:
            if  input == 0.05 or input == 0.10:
                break
            else:
                print 'Out of range. Try again'
    return input
       
def sttninput(indexlist):
    #prompt for user to input the station index
    while True:
       try:
           input = int(raw_input('\nSelect a station index from'
                                     ' the list above >>> '))
       except ValueError:
           print 'That\'s not a valid station index! Try again'
       else:
           if  input in indexlist:
               break
           else:
               print 'That\'s not a valid station index!. Try again'
    return input    
    #change this to use the index instead, easier for the user to enter an integ

def nextinput():
    #prompt user for next action
    while True:
        try:
            print '\nSelect 0 to Exit'
            print 'Select 1 to enter another weather station'
            print 'Select 2 to change the significance level'
            input = int(raw_input('Select 0,1,or 2 >>> '))
        except ValueError: 
            print 'That\'s not an integer! Try again'
        else:
            if  0 <= input <= 2:
                break
            else:
                print 'Out of range. Try again'
    return input

def corr_check(newlistRows,newNAsnlist,filepath,sig_level):
    #this module checks the subset of files for correlations and returns
    #a list of weather stations with their correlations, if significant as
    #defined by user input significance level.
    
    corr_list = []
    snow_corr=[]
    snwd_corr=[]
    snow_pval=0
    snwd_pval=0

    for i in range(0,newlistRows):
        #compute the pearson correlation between the Snowfall delta and the
        #Nino Index, and between the Snowpack delta and the Nino Index
        StationFile = newNAsnlist[i]+'.dly'
        ws = print_sttn_details(StationFile,filepath) #get weather station state and name
        snow_corr = pearson_corr(StationFile,filepath,'SNOW','Snowsum','30YearAvg',
                        'SnowDelta')
        snwd_corr = pearson_corr(StationFile,filepath,'SNWD','Snwdavg','30YearAvg',
                     'SnwdDelta')

        #If the p-value is greater than or equal to the significance level,
        #change correlation and p-value to NA
        if snow_corr[1]>=sig_level:
            snow_corr=[np.nan,np.nan]
        if snwd_corr[1]>=sig_level:
            snwd_corr=[np.nan,np.nan]
            
        #If either the snowfall or snowpack correlations are significant
        #return the station file and its correlation and p-value
        if snow_corr[1]<sig_level or snwd_corr[1]<sig_level:
            corr_list=corr_list+[[newNAsnlist[i],ws[0],ws[1],snow_corr[0],snow_corr[1],
                                   snwd_corr[0],snwd_corr[1]]]

    #Determine the number of weather stations with statistically significant
    #correlations
    corr_listRows=len(corr_list)

    return [corr_list,corr_listRows]

def sig_corr_file_copy(corr_list,corr_listRows):
    mypath = 'new3_ghcnd_all'
    if not os.path.isdir(mypath):
        os.makedirs(mypath)
    for i in range (0,corr_listRows):
        filename = 'new2_ghcnd_all/'+corr_list[i][0]+'.dly'
        shutil.copy2(filename, mypath)

def df_corr_list(corr_listRows,corr_list,filepath):
    #create a dataframe from the corr_list
    if corr_listRows>0:
        dfcorr=pd.DataFrame(corr_list)
        dfcorr.columns = ['StationCode','StationName','ST/Prov','SnowCorr','SnowPval',
                          'SnwdCorr','SnwdPval']
        return dfcorr

def print_corr_list(corr_listRows,dfcorr,filepath):
    #prints information to the console about the weather stations that have
    #a significant correlation
    if corr_listRows>0:
        print ('There are '+str(corr_listRows)+
               ' North American weather stations with snowfall (SNOW) and '
               'snowpack (SNWD) information for all years between 1950 and 2010 '
               'with statistically significant correlations* between the Snowfall'
               'delta and/or the Snowpack delta and the Nino Index. These values '
               'appear in the following table. '
               '\n\n*If only one of these correlations is significant for a '
               'particular weather station, then the station will appear in the table,'
               'but "NaN" will be shown for the non-statistically significant '
               'relationship.\n')
        print dfcorr 
        #dfcorr.to_csv('CorrFile')
    else:
        print ('There are no North American weather stations with statistically'
               ' significant correlations between either the Snowfall delta or the'
               ' Snowpack delta and the Nino Index at the selected level of '
               'significance.')

def sttn_plot(corr_listRows,dfcorr,filepath):
    #Create scatter plots and regression lines for the stations with
    #statistically signicant correlations.  

    if corr_listRows>0:
        #create station list of stations with significant correlations at
        #user selected level of significance
        #sttnlist = dfcorr['StationCode'].tolist()
        indexlist = dfcorr.index.tolist()
        wsindex = sttninput(indexlist) #user selects station index from list
        wstation = dfcorr['StationCode'][dfcorr.index == wsindex].tolist()[0]
        #use index to get station name 
        StationFile=str(wstation)+'.dly'  #file name

        wsinfo = print_sttn_details(StationFile,filepath) #get weather station state and name

        #print weather station code, name and state
        print ('You have selected ' + str(wstation)+' '+wsinfo[0]+' in '+wsinfo[1]+
               ', a scatter plot with regression line for this station will appear '
               'in a new window')
        print ('\nNote: You must close the plot window in order to continue')

        plt.close('all') 
            
        plt.figure(1)
        subplot(2,1,1)
        scatter(StationFile,filepath,'SNOW','Snowsum','30YearAvg','SnowDelta','Snowfall Delta',
                'Nino Index','Snowfall Delta vs Nino Index')
        subplot(2,1,2)
        scatter(StationFile,filepath,'SNWD','Snwdavg','30YearAvg','SnwdDelta','Snowpack Delta',
                'Nino Index','Snowpack Delta vs Nino Index')
        plt.show()

    else:
        print ('\nSince there are no North American weather stations with '
               'statistically significant correlations there are no valuable '
               'scatter plots or regression lines to show. \n\nYou can change '
               'the level of significance by entering 2 below and following '
               'the prompt to enter a different significance level')

def create_tar_file():
    #once the data has been fully processed create a tar file with just the
    #the files with significant correlations

    print 'creating archive...'
    out = tarfile.open('ghcnd_processed.tar', mode='w')
    try:
        print ('Creating tar file archive for files with p-value less than '
               'or equal to 0.10')
        out.add('new3_ghcnd_all',arcname='ghcnd_processed')
    finally:
        print 'closing'
        out.close()

def data_processing():
    #Extract all 93,000+ weather stations from the tar.gz file
    print '\nExtracting all weather station files. This may take several minutes'
    print 'Please wait...'
    get_tar_file('ghcnd_all.tar.gz',"r:gz")
    #the above takes about 10 minutes to run depending on processing speed

    #Part 1
    #Identify weather stations with sufficient data records and copy files
    #to new directories

    #Using the station and inventory files determine which North American
    #stations have SNOW and SNWD data for 1950 and 2010.

    print ('\nFinding North American stations with snowfall and snowpack data '
           'for 1950 and 2010 in the data inventory file')
    print 'Please wait...'

    #find North American stations
    NAsnlist = NoAM5010Station()  #list of relevant files
    listRows = len(NAsnlist)

    print ('\nThere are '+str(listRows)+
           ' North American weather stations with snowfall (SNOW) and '
           'snowpack (SNWD) information for the years 1950 and 2010.\n')
    print 'Copying files, this may take several minutes, please wait...'

    #Copy these files to a new folder
    NA5010_file_copy(NAsnlist,listRows)
    #the above takes about 8 minutes to run, depending on processing speed

    print '\nCompiling the list of weather stations with a complete record'
    print 'This could take as long as 45 minutes. Please wait...'
    
    #Check if stations have a complete record
    #full record is 61 years (1950-2010) * 7 months(Oct - Apr)
    full_record = 427
    newNAsnlist = full_record_copy(listRows, NAsnlist,full_record,
                                   'new_ghcnd_all')
    #the above takes about 45 minutes to run

    newlistRows = len(newNAsnlist)
    
    print ('\nThere are '+str(newlistRows)+
           ' North American weather stations with snowfall (SNOW) and '
           'snowpack (SNWD) information for ALL years between 1950 '
           'and 2010.\n')
    print 'Copying files, please wait...'
    
    #Copy files with full record to a new folder
    NA5010_inclusive_file_copy(newNAsnlist,newlistRows)
    #only 83 stations have a complete record

    print '\nCompiling the list of weather stations with a p-value <= 0.10'
    print 'Please wait...'

    #Part 2:
    #Perform correlation analysis on the reduced list of stations
    #Create a list of stations with either a snowdepth or snowpack p-value
    #less than or equal to 0.10
    [corr_list,corr_listRows] = corr_check(newlistRows,newNAsnlist,
                                           'new2_ghcnd_all',0.10)
    #the above takes about 10 minutes to run

    print ('\nThere are '+str(corr_listRows)+
           ' North American weather stations with p-value less than'
           'or equal to 0.10\n')
    print 'Copying files, please wait...'

    sig_corr_file_copy(corr_list,corr_listRows)

    #Create tar file 'ghcnd_processed.tar'
    create_tar_file()

    print '\nData extraction and processing complete\n'

def processed_list():
    newNAsnlist = []
    files_in_dir = os.listdir('ghcnd_processed')
    for file_in_dir in files_in_dir:
        filename = file_in_dir
        filename = filename.rstrip('.dly')
        newNAsnlist = newNAsnlist + [filename]
    return newNAsnlist

def run_corr_and_plot(newlistRows,newNAsnlist,filepath):
    #For each station calculate the Pearson Correlation by element type:
    #snowfall and snowpack.  The correlation is statistically significant only
    #if the p-value is less than 0.05. However, this program allows the user to
    #input any significance level anywhere between 0 and 1.

    nextAction = 2

    #user can select 0,1, or 2 at the end of each loop:
    # 0 [exit]
    # 1 [choose new station to plot]
    # 2 [choose new significance level]
    
    while nextAction != 0:
        #Return the list of stations with statistically significant correlations, with
        #its respective correlation and p-value.
        if nextAction == 2:
            sig_level=siginput()  #user enters significance level
            print ('\nGenerating list, this may take several minutes, please wait...\n')
            [corr_list,corr_listRows]=corr_check(newlistRows,newNAsnlist,filepath,sig_level)
            #the above takes about 5 minutes, a little longer for sig_level = 0.10
            
            #create a dataframe from the corr_list
            dfcorr = df_corr_list(corr_listRows,corr_list,filepath)   #create data frame
            print_corr_list(corr_listRows,dfcorr,filepath)   #print comments and data frame

            #Create scatter plots and regression lines for the stations with
            #statistically signicant correlations.
            sttn_plot(corr_listRows,dfcorr,filepath)
            
        else:
            if nextAction == 1:
                sttn_plot(corr_listRows,dfcorr,filepath)

        nextAction = nextinput()
    
def main():

    #show introductory text
    intro()
    
    #user will be prompted to decide if they want to go through the data
    #processing process, or use the pre-processed files

    #prompt user to enter 0 to process and extract data or 1 to use
    #pre-processed file:
    process = processinput()

    if process == 0:
        data_processing()

    #Extract the files from the pre-processed file
    #This will create a folder named 'ghcnd_processed', and extract the
    #pre-processed files to this folder

    get_tar_file('ghcnd_processed.tar','r')

    #For each station calculate the Pearson Correlation by element type:
    #snowfall and snowpack.  The correlation is statistically significant only
    #if the p-value is less than 0.05 or possibly up to 0.10. This program
    #allows the user to choose a significance level (p-value) of 0.05 or 0.10

    newNAsnlist = processed_list()
    newlistRows = len(newNAsnlist)
    run_corr_and_plot(newlistRows,newNAsnlist,'ghcnd_processed')


if __name__ == "__main__":
    main()

