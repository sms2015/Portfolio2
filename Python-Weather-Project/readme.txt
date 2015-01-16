ABOUT THIS PROGRAM
This program uses data from the NOAA Global Historical Climatology Network (GCHN) to analyze historical 
snowfall and snowpack at North American weather stations and looks for correlations between these measures 
and historical El Nino/La Nina episodes between the years 1950 and 2010.  The Nino Index is a measurement 
of the variation from a 30 year average. This program calculates similar variations from the snowfall and 
snowdepth 30 year averages and looks for correlations between these values and the Nino Index. These variations 
are referred to as "delta" in this program. More information about the data used in this analysis can be 
found here: http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt

More detailed information about the work done in this program can be found in the project report: 
'IS 602 Final Project.pdf' also found in the project tar file.

This program will provide prompts to the user to select the statistical significance level for the correlation 
analysis, and to then select a weather station to view a scatter plot of the data for that weather station.

This program was created as final project for a Data Science Python programming class. As such, the emphasis 
was mostly on the code and less on the scientific process, so there are some potential scientific issues with
this analysis. For example, in order to simplify the coding for this project this program ignores missing 
daily values when calculating monthly averages.

HOW TO USE:
download the IS602_FinalProject.tar, this file contains all files needed EXCEPT the ghcd_all.tar.gz file
which was too large to include in this file.

The user will be presented with two options
Select 0 to run all of the extraction and processing functions to create the processed data files
or
Select 1 to use the pre-processed file, which is the tar file: ghcnd_processed.tar


If 0 is selected, at least 30GB of free space will be required for extracting and copying files, and the
extraction and processing can take up to an hour to complete.

You must have all of the following files in your directory to run option 0:
ghcnd_all.tar, ghcnd-inventory.txt, ghcnd-stations.txt,Nino_index.csv

The ghcnd_all.tar.gz file must be downloaded from the NOAA ftp site:
ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/


If 1 is selected the ghcnd_processed.tar file will be used. This archive contains 37 files, which contain
complete data records (all years, all months between October and April) for snowfall and snowdepth for 
the years 1950 through 2010, which also have significant correlations to the Nino Index when considering 
average monthly variations from the 30 year monthly average.

You must have all of the following files in your directory to run option 1
ghcnd_processed.tar, ghcnd-inventory.txt, ghcnd-stations.txt, Nino_index.csv

METHODOLOGY AND CODING
Importing and Extracting the data: 
Data (“.dly” files) are extracted from the TAR file using the tarfile module.
In subsequent steps using the metadata files "ghcnd-stations" and "ghcnd-inventory" are used to pare down 
the list of all station files to just those that are North American stations where there there is data for
1950 and 2010.

To further reduce the files to just those that are relevant a station by station analysis is completed, 
to first look for stations with complete data records (see above) and finally correlations with the 
Nino Index as described briefly above.

STATION ANALYSIS
Note: Each “.dly” file has an identical structure

Monthly Average: 
Computes a total snowfall and a snowpack average for for the “snow months” defined in 
this project as October – April for all years (1950 to the present).

30 Year Monthly Average: 
Computes an average monthly snowfall (an average of the monthly snow totals) and 
snowpack for the snow months (October – April) for the 30 year period 1981 – 2010. This is the most recently 
completed 30 year period and is commonly used to calculate climatological “normals” in meteorology. It 
should be noted that the ONI data is based on multiple centered thirty year periods; however, that is not 
the standard practice in meteorology and it would significantly complicate the analysis for this project to 
follow the ONI method.  Therefore only the 1981-2010 thirty year period will be used to 
compute snowfall and snowpack averages.

Monthly Snowfall vs 30 Year Average: 
Calculates the snowfall/snowpack deltas between the values in step 1 (monthly averages) and the values in 
step 2 (30 year averages). This step will create an array of monthly delta values (October-April) for 
the weather station, for each year in the ONI record (1950 – present).

Correlation Analysis: 
Using scipy.stats.pearsonr to compute the correlation between each array (snowfall and snowpack)
in step 3 and the ONI data set.(Each data point in the ONI table is a three month rolling average so a 
January-February-March 2010 value would be compared to a February 2010 value from the snowfall/snowpack arrays).  

Scatter Plot and Regression analysis: :
A scatter plot with a linear regression line for snowfall and snowpack can be displayed for stations that show
statistically significant correlations. For the purpose of this project, Pearson correlation p-value = 0.10 or less.  
