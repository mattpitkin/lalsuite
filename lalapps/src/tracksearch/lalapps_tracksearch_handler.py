#
# Copyright (C) 2004, 2005 Cristina V. Torres
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
__author__ = 'Cristina Torres <cristina@phys.utb.edu>'
__date__ = '$Date$'
__version__ = ''

import sys
from optparse import OptionParser
import ConfigParser
import getopt
import math
import os
import string
import time
import copy
try:
    from lalapps.tracksearchutils import *
except:
    from tracksearchutils import *

disableGraphics=False
try:
    from pylab import *
except Exception, errorInfo:
    disableGraphics=True
    sys.stderr.write("Error trying to import pylab!\n")
    sys.stderr.write("Exception Instance :%s\n"%(str(type(errorInfo))))
    sys.stderr.write("Exception Args     :%s\n"%(str(errorInfo.args)))
    sys.stderr.write("Pylab functionality unavailable!\n")
    sys.stderr.write("Function calls requiring graphics may FAIL!\n")
    sys.stderr.flush()

#Begin MAIN part of this python code to actually manipulate the curve
#lists


#################### Internal Methods #################

def buildCandidateGlob(fileList,verboseSwitch=False):
    #Build glob var from filelist.  Memory friendly
    #Returns copy of results
    canList=fileList
    if canList.__len__() == 1:
        if verboseSwitch:
            print "Globbing is not possible only one file found!"
            print "Returning the result as that entry!"
        tmpObject=candidateList()
        tmpObject.__loadfileQuick__(canList.pop(0))
        newCandidateObject=copy.deepcopy(tmpObject)
    else:
        newCandidateObject=candidateList()
        newCandidateObject.__loadfileQuick__(canList.pop(0))
    for entry in canList:
        if verboseSwitch:
            print " "
            print "Loading candidate list file: ",entry
        tmpCandidate=candidateList()
        tmpCandidate.__loadfileQuick__(entry)
        newCandidateObject=newCandidateObject.globList(tmpCandidate,True)
        del tmpCandidate
    return newCandidateObject
#End method

def gnuplotScriptFile(filename):
    """
    Invokes gnuPlot in an unusual hackish was to create a PS scatter
    plot to view the pixels found as part of a curve(s). Add a .plt
    extension to the filename specified
    """
    txtTemplate='plot "%s" with lines\n set size 1.0,0.6\n set terminal postscript portrait enhanced mono dashed lw 1 "Helvetica" 14\n set output "%s.ps"\n replot\n set terminal x11\n set size 1,1\n set title "%s"\n'
    output_fp=open(filename+'.plt','w')
    output_fp.write(txtTemplate%(filename,filename,os.path.basename(filename)))
    output_fp.close()
#End method gnuplotScriptFile()

def wrapperGraphTriggers(candidateObject):
    """
    Optimized internal function
    """
    global graphGPSoffset
    global timescale
    global graph2file
    try:
        if candidateObject.curves==[]:
            if candidateObject.verboseMode:
                sys.stdout.write("Only trait summary for triggers in memory.\n")
                sys.stdout.write("Loading entire candidate list.\n")
                sys.stdout.flush()
            candidateObject.__loadfileQuick__(candidateObject.filename[0])
        candidateObject.graphTriggers(graph2file,graphGPSoffset,timescale)
    except Exception, errorInfo:
        sys.stderr.write("Problem trying to create plot, %s !\n"%(os.path.basename(graph2file)))
        sys.stderr.write("Exception Instance :%s\n"%(str(type(errorInfo))))
        sys.stderr.write("Exception Args     :%s\n"%(str(errorInfo.args)))
        sys.stderr.write("Skipping this plot. Sorry.\n")
        sys.stderr.flush()
#End wrapperGraphTriggers(candidateObject)

def wrapperShowHistogram(candidateObject):
    """
    Optimized internal function
    """
    global triggerTrait
    global verboseSwitch
    minBinValue=0
    for singleTriggerTrait in triggerTrait:
        newHistCol=histogramcolumns
        if singleTriggerTrait.__contains__(":"):
            property=singleTriggerTrait.strip().split(":")
            myTrait=str(property[0])
            newHistCol=int(property[1])
            if property.__len__() > 2:
                if str(property[2]).__contains__('!'):
                    property[2]=property[2].strip('!')
                    minBinValue=-1*float(property[2])
                maxBinValue=float(property[2])
                #Since we have a max bin and bin count
                #make newHistCol not a number but vector of
                #bins to histogram with.
                resolution=(maxBinValue-minBinValue)/newHistCol
                binList=[]
                for index in range(0,newHistCol,1):
                    binList.append((index*resolution)+minBinValue)
                newHistCol=binList
            else:
                maxBinValue=float(-1)
        else:
            myTrait=singleTriggerTrait

        if verboseSwitch:
            property=candidateObject.__getCurveField__('',myTrait)
            sys.stdout.write("Creating histogram for %s \n"%
                             (property[1]))
            sys.stdout.flush()
        try:
            candidateObject.setHistogramType(histogramType)
            candidateObject.showHistogram(myTrait,
                                          graph2file,
                                          newHistCol,
                                          topTol)
        except Exception, errorInfo:
            sys.stderr.write("Problem trying to create histogram plot ::: %s :::\n"%(os.path.basename(graph2file)))
            sys.stderr.write("Exception Instance :%s\n"%(str(type(errorInfo))))
            sys.stderr.write("Exception Args     :%s\n"%(str(errorInfo.args)))
            sys.stderr.write("Skipping this plot. Sorry.\n")
            sys.stderr.flush()
#End wrapperShowHistogram(candidateObject)

def wrapperShowScatterPlot(candidateObject):
    """
    Creates a scatter plot of paired variables and prints
    them to the screen of a file if requested.
    """
    global triggerTrait
    global verboseSwitch
    for singleTriggerTrait in triggerTrait:
        try:
            (colX,colY)=singleTriggerTrait.strip().split(":")
        except:
            sys.stderr.write("Error plotting this scatter plot!\n")
            sys.stderr.write("Check trigger property string.\n")
            sys.stderr.write("%s\n"%(triggerTrait))
            sys.stderr.flush()
            if (triggerTrait=="p"):
                colX="p"
                colY="v"
            else:
                return
        myTraitX=str(colX)
        myTraitY=str(colY)
        if verboseSwitch:
            propertyX=candidateObject.__getCurveField__('',myTraitX)[1]
            propertyY=candidateObject.__getCurveField__('',myTraitY)[1]
            sys.stdout.write("Creating scatter plot for %s,%s\n"%
                             (propertyX,propertyY))
            sys.stdout.flush()
        try:
            formatString=''
            candidateObject.showScatterPlot(myTraitX,
                                            myTraitY,
					    graph2file,
                                            formatString)
        except Exception, errorInfo:
            sys.stderr.write("Problem trying to create scatter plot ::: %s :::\n"%(os.path.basename(graph2file)))
            sys.stderr.write("Exception Instance :%s\n"%(str(type(errorInfo))))
            sys.stderr.write("Exception Args     :%s\n"%(str(errorInfo.args)))
            sys.stderr.write("Skipping this plot. Sorry.\n")
            sys.stderr.flush()
            
#End wrapperShowScatterPlot(candidateObject)

###################End internal methods #################
usage=''
if disableGraphics:
    usage = "X windows related methods available.\n tshandler [args]"
else:
    usage = "tshandler [args]"
parser = OptionParser(usage)
parser.add_option("-f","--file",dest="filename",
                  default="",
                  help="This specifies the list of curve files to process contained in FILE, either rethresholding, globbing.  If the input is a directory then it is assumed that we will be working on all the seen files.  /PATH/*.ext is allowed.  Only one of the other following options should be specified at run-time.  Combining the following options will yield unpredictable results!",
                  metavar="FILE"
                  )
parser.add_option("-g","--glob",dest="glob",
                  default=False,
                  action="store_true",
                  help="This option tells the code to glob all the lists together into one larger curve candidate listing.  The output be a single or sequence of files with GLOB:Start:00000000,0:Stop:00000000,0:TF:000:000:.candidates. This option can not be used with clobber. This is a RAM and CPU intensive task! A quicker option assuming the candidate files are generated with the same paramters is --glue.  To use this then invoke --glob --glue",
                  )
parser.add_option("-k","--show-legend",dest="legend",
                  default=False,
                  action="store_true",
                  help="This shows the properties that can be accessed from the triggers in the trigger file.  These are typically useful for plotting properties relationships."
                  )
parser.add_option("-l","--glue",dest="glue",
                  default=False,
                  action="store_true",
                  help="This only works if --glob has been invoked also.  It forces all found candidate files to be concatenated together.  This does not account for differently constructed candidate libraries.  In that case a measure of track frequency and time resolution must be made before trying to join multiple libraries."
                  )
parser.add_option("-o","--outfile",dest="outfile",
                  default="",
                  help="This overrides the automagic output file nameing characteristics of this utility.  Please invoke the option to write to the output to a definite file location of your choice."
                  )
parser.add_option("-c","--clobber",dest="clobberFilename",
                   default="",
                   help="This is a file or directory/mask of all files that should be globbed together then checked against the globbed list from --file option.  Entries candidates in FILE with reside in list CLOBBER will be deleted from list FILE and written to disk along with an globbed list which is not clobbered.The filename will be like CLOBBER:Start:00000000,0:Stop:00000000,0:TF:000:000:AND:Start:00000000,0:Stop:00000000,0:TF:000:000:.candidates.  This option can not be used with glob, we assume specified --file are globs of candidates.   We glob the --clobber option if neccessary to perform the clobber."
                   )
parser.add_option("-K","--percentile-cut",dest="percentileCut",
                  default="",
                  help="A string formated by property letter ex: P or \
!P (negate) then a colon then the decimal percentage \
P:0.05 if making several different cuts separete by commas, EX: P:0.01,L:0.50")
parser.add_option("-T","--expression-threshold",dest="expThreshold", 
                  default="",
                  help="This is a flexible thresholding interface.  We are allowed to manipulate variables to build expressions. Variables allowed are:\n P ,power \n L ,curve length \n D , time duration in seconds \n B , frequency bandwith of kurve \n  T, start time float\n S, stop time float\n F, start Freq \n, G, stop Freq \n V, bright pixel frequency\n H, bright pixel time(float)\n J, power of bright pixel \n M, mean pixel power \n C, variance pixel power \n A, trigger angle \n See the help for the candidateList class for better explaination of the valid syntax allowed. You can impose multiple threshold expresssions for one file operation.  This saves time lost to IO of running candidateHandler.py multiple times.  Use the comma as a delimiter B>10,L<100 would apply to thresholds and save them in two appropriately marked files. ENCLOSE THE EXPRESSION IN DOUBLE QUOTES",
                  )
parser.add_option("-s","--write-summary",dest="dumpSummaryDisk",
                  default=False,
                  action="store_true",
                  help="This will write the summary information for the candidate file(s) opened. The filename scheme replaces candidate with summary. The data store is Length Pixels,Power(actually energy),Duration Sec,Bandwidth Hz.  This option can be invoked along with the thresholding option to write out summary files of the results immediately after thresholding the data."
                  )
parser.add_option("-d","--display-summary",dest="dumpSummaryScreen",
                  default=False,
                  action="store_true",
                  help="This will display the summary information for the candidate file(s) opened."
                  )
parser.add_option("-p","--print",dest="imageFilemask",
                  default='',
                  help="This option states that files specified with the --file option should be printed to a PNG file.  This method will name the file according to the directive given on the command line.  If more than one image will be created then depending on the --file option then the AUTO argument will be invoked by force to create multiple graphic files.  If the AUTO option is given then the image will have the same name as the candidate library file with the extension changed to PNG."
                  )
parser.add_option("-L","--line-plot",dest="linePlot",
                  default=False,
                  action="store_true",
                  help="This flag must be set to invoke creating the line plot of triggers from the candidate object."
                  )
parser.add_option("-i","--image",dest="imageCount",
                  default='0',
                  help="This switch allows you to tell the candidateHandler to graph the data, if the input library is clobbered together from many candidate lists then specify the maximum number of figures you are willing to create.  The default is currently 1")
parser.add_option("-r","--graphoffset",dest="graphGPSoffset",
                  default=float(0),
                  help="Specify GPS offset if graphing."
                  )
parser.add_option("-t","--timescale",dest="timescale",
                  default='second',
                  help="Specify an alternative plot time scale.  The valid options are second, hour, and day. Anything else will result in a scale of seconds in the plot."
                  )
parser.add_option("-x","--stats",dest="stats",
                  default='',
                  help="Determine the specified stats for all candidates in input file for variable in question.  See the -T options variable labels for explaination.  This result is written out to the screen as mean,min,max,std.  To measure stats of more than one trait use a comma seperated list P,L as input args. FUNCTIONALITY NOT IMPLEMENTED YET!"
                  )
parser.add_option("-b","--glitchdb",dest="glitchdb",
                  default=False,
                  action="store_true",
                  help="Invokes building a glitch database file for trigger auto-classification. Output is writting to file from INPUT.file the name of the output file will be called FILE.glitchDB"
                  )
parser.add_option("-v","--verbose",dest="verbose",
                  default=False,
                  action="store_true",
                  help="Sets up the verbose features. Prints diagnostic messages and progress meters where appropriate.")
parser.add_option("-V","--mute-spinner",dest="muteSpinner",
                  default=False,
                  action="store_true",
                  help="Setting this option leaves all diagnostic messages, if they are set.  It turns off the progress meter.")
parser.add_option("-n","--bin-count",
                  default=50,
                  dest="histCol",
                  help="Sets up,for example, the integrated power histogram column count so we can see significant outliers in the trigger distributions.")
parser.add_option("-m","--histogram",
                  default=False,
                  dest="histogram",
                  action="store_true",
                  help="This option will show a histogram. For the candidate trigger library.")
parser.add_option("-e","--scatter",
                  default=False,
                  dest="scatter",
                  action="store_true",
                  help="This option will show a scatter plot of paired trigger properties. These properties should be specified as the argument to --trigger-property as D:P,A:F. So refer to available property tags in this help menu. The default scatter plot if --trigger-property is not specified will be P:V.  So for the default the X-axis is Integrated Power and the Y-axis is Bright Frequency."
                  )
parser.add_option("-q","--trigger-property",
                  default="p",
                  dest="triggerTrait",
                  help="Please give single character field as argument.  See the explaination for thresholding for allowable fields. You can specify multiple histogram to create via using comma delimited entries IE p,d,v.  To override the default value of histCol for any particular plot use a : for example p,d:10,v means plot the histogram for p with histCol columns but for d use 10 columns instead. If you want to set a maximum value to histogram up to use, try instead d:10:1274, which means take property d create 10 bins equally spaced up to value 1274"
                  )
parser.add_option("-P","--top-percentile",
                  default=0.05,
                  dest="topPercentile",
                  help="This options set the histogram top percentile option for the histogram. In floating point style 1% --> 0.01")
parser.add_option("-z","--histogram-type",
                  default='default',
                  dest="histogramType",
                  help="This sets the type of histrogram to create if a histogram is invoked for creation either to screen or disk.  The options are default, logy logxy." )
parser.add_option("-D","--dqlist",
                  default=-1,
                  dest="makeDQlist",
                  help="This creates a file of DQ flags from candidateList object using a padsize specified as the argument.")

(options,args)=parser.parse_args()
filename=str(options.filename)
glitchdb=bool(options.glitchdb)
glob=options.glob
clobberFilename=str(options.clobberFilename)
expThreshold=str(options.expThreshold)
percentileCut=str(options.percentileCut)
graph2file=options.imageFilemask
linePlot=bool(options.linePlot)
createHistogram=bool(options.histogram)
createScatterPlot=bool(options.scatter)
triggerTrait=str(options.triggerTrait).split(',')
histogramcolumns=int(options.histCol)
graph2screen=int(options.imageCount)
graphGPSoffset=float(options.graphGPSoffset)
outfile=str(options.outfile)
dumpSummaryScreen=bool(options.dumpSummaryScreen)
dumpSummaryDisk=bool(options.dumpSummaryDisk)
verboseSwitch=bool(options.verbose)
muteSpinner=bool(options.muteSpinner)
measureStats=str(options.stats)
glue=bool(options.glue)
timescale=str(options.timescale).lower()
topTol=float(options.topPercentile)
histogramType=str(options.histogramType).lower()
dqList=int(options.makeDQlist)
showLegend=bool(options.legend)

if verboseSwitch:
    print "Setting verbose mode to candidateHandler call."

if (filename == "" and not showLegend):
    print "Filename argument either not specified or invalid!"
    os.abort()

#Load the file/dir/single file specified by --file option
canList=[]
canList=generateFileList(filename)

if showLegend:
    x=candidateList()
    myLegend=x.showQualityLegend()
    for line in myLegend:
        sys.stdout.write(line)
#SECTION TO DO THE GLOBBING OF MANY CANDIDATE FILES
elif (glob and (canList.__len__() >= 1)):
    if outfile != "":
        outName=outfile
        #If file preexists erase it first!
        if os.path.isfile(outName):
            if verboseSwitch:
                print "Prexisting file found:",outName
                print "Removing!"
            os.unlink(outName)
    #Depending on the presence of the --glue switch we actually read
    #the data or just concatenate the files together.
    #Thu-Oct-11-2007:200710111752 
    if glue:
        #Just join the files
        nullLibrary=candidateList()
        for libraryEntry in canList:
            if verboseSwitch:
                sys.stdout.write("Processing file : %s"%libraryEntry)
                sys.stdout.flush()
            nullLibrary.globListFile(outName,libraryEntry)
    else:
        #Glob them through RAM allowing the proper adjustments to be made
        newGlobFile=buildCandidateGlob(canList,verboseSwitch)
        #Set output file name
        if outfile !="":
            outName=outfile
        else:
            outName="Glob:"+newGlobFile.__filemaskGlob__()
        #Write this globbed file to disk.
        newGlobFile.writefile(outName)

    #SECTION TO DO THE CLOBBERING OF A CANDIDATE FILE WITH ANOTHER
elif (clobberFilename != '') and (canList.__len__() == 1):
    #Create clobberer
    #Load file(s) to clobber with
    clobberList=generateFileList(clobberFilename)
    #If there is more than one file we need to glob them first
    newCandidateClobberObject=buildCandidateGlob(clobberList,verboseSwitch)
    #Create candidate list to clobber.
    clobberVictim=candidateList(verboseSwitch)
    if muteSpinner:
        clobberVictim.spinnerVerboseMode=False
    clobberVictim.__loadfileQuick__(canList[0])
    #Clobber the Victim
    newClobberedList=clobberVictim.clusterClobberWith(newCandidateClobberObject)
    #Write the results to the disk!
    if outfile != "":
        newClobberedList.writefile(outfile)
    else:
        newClobberedList.writefile(clobberVictim.__filemaskClob__(newCandidateClobberObject))
    #
    # Add logic branch to do a report graphic...  Try original trigger
    # curve plots, show power histogram, then show a cut trigger plot
    # 
    #Create a lineplot on screen of the curves in the candidateList
elif (linePlot and not createHistogram and (canList.__len__() >=1) and ((graph2screen>0) or (graph2file!=''))):
    if ((graph2screen < canList.__len__()) and (graph2file=='')):
        print "Number of images to create exceeds your request!"
        print "Estimate images in this data structure are ",canList.__len__()
        #Needs continue anyway option here!
        os.abort()
    if ((graph2file != 'AUTO') and (canList.__len__() > 1) and (graph2screen==0)):
        graph2file='AUTO'
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        if muteSpinner:
            candidateObject.spinnerVerboseMode=False
        candidateObject.__loadfileQuick__(entry)
        wrapperGraphTriggers(candidateObject)
        del candidateObject
    #SECTION CREATE HISTOGRAMS
elif (createHistogram and (canList.__len__() >=1)):
    if ((graph2file != 'AUTO') and (canList.__len__() > 1) and (graph2screen==0)):
        graph2file='AUTO'
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        if muteSpinner:
            candidateObject.spinnerVerboseMode=False
        candidateObject.__propertyunpickler__(entry)
        wrapperShowHistogram(candidateObject)
        if linePlot:
            wrapperGraphTriggers(candidateObject)
        del candidateObject
    #SECTION CREATE SCATTER PLOTS
elif (createScatterPlot and (canList.__len__() >=1)):
    if ((graph2file != 'AUTO') and (canList.__len__() > 1) and (graph2screen==0)):
        graph2file='AUTO'
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        if muteSpinner:
            candidateObject.spinnerVerboseMode=False
        candidateObject.__propertyunpickler__(entry)
        wrapperShowScatterPlot(candidateObject)
        if linePlot:
            wrapperGraphTriggers(candidateObject)
        del candidateObject
    
    #SECTION APPLY ABITRARY THRESHOLDS
    #Adjust logic to do expThres or percentileCut together!!!
    #Use to local defs to clean up this branch in script!
    # A def for -T and a def for -K
elif (((expThreshold != "") or (percentileCut != ""))and (canList.__len__() >=1)):
    #Treshold will always go before percentile Cuts
    #Carry out thresholding on all entries in canList
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        if muteSpinner:
            candidateObject.spinnerVerboseMode=False        
        candidateObject.__loadfileQuick__(entry)
        #Build Threshold list
        if (expThreshold != ""):
            if expThreshold.__contains__('"'):
                if verboseSwitch:
                    sys.stdout.write("Found extra parenthesis in threshold expressions!")
                    sys.stdout.flush()
            expThreshold=expThreshold.strip('"')
            expThresholdLIST=str(expThreshold).split(',')
        else:
            expThresholdLIST=[None]
        #Build Percentile Cut list
        if (percentileCut != ""):
            perCut=percentileCut.strip('"')
            perCutLIST=str(perCut).split(',')
        else:
            perCutLIST=[None]
        #Process both lists
        for singleThreshold in expThresholdLIST:
            for singleCut in perCutLIST:
                singleThresholdName=str(singleThreshold).replace('(','--').replace(')','--')
                singleCutName=str(singleCut).replace('(','--').replace(')','--')
                #Fix up name schema and apply cuts as appropriate
                if ((singleCut != None) and \
                    (singleThreshold != None)):
                    singleCompleteName=singleThresholdName+":"+singleCutName
                    candidateResults=candidateObject.applyArbitraryThresholds(singleThreshold)
                    argProp,argVal=singleCut.split(':')
                    candidateResults2=candidateResults.applyPercentageThreshold(argProp,argVal)
                    candidateResults=candidateResults2
                elif singleCut != None:
                    singleCompleteName=singleCutName
                    argProp,argVal=singleCut.split(':')
                    candidateResults=candidateObject.applyPercentageThreshold(argProp,argVal)
                elif singleThreshold != None:
                    singleCompleteName=singleThresholdName
                    candidateResults=candidateObject.applyArbitraryThresholds(singleThreshold)
                else:
                    singleCompleteName="FUNCERR"
                    candidateResults=candidateObject
                if (outfile != "") and (canList.__len__() == 1):
                    candidateResults.writefile(outfile)
                    if verboseSwitch:
                        print "Wrote file :",outfile
                    if dumpSummaryDisk:
                        candidateResults.writeSummary(outfile)
                    if glitchdb:
                        candidateResults.writeGlitchDatabase('',outfile)
                    if (dqList>=0):
                        candidateObject.writeDQtable(dqList,outfile)
                else:
                    pathName=os.path.dirname(entry)
                    saveFiles=pathName+'/Threshold:'+str(singleCompleteName)+':'+os.path.basename(entry)
                    candidateResults.writefile(saveFiles)
                    if verboseSwitch:
                        print "Wrote file :",saveFiles
                    if dumpSummaryDisk:
                            candidateResults.writeSummary(saveFiles)
                    if glitchdb:
                            candidateResults.writeGlitchDatabase('',saveFiles)
                    if (dqList>=0):
                        candidateObject.writeDQtable(dqList,outfile)
                    candidateResults.__setfilename__(saveFiles)
                    if createHistogram:
                        wrapperShowHistogram(candidateResults)
                    if linePlot:
                        wrapperGraphTriggers(candidateResults)
                del candidateResults
        del entry
    del candidateObject
# #Do percentile cuts
# elif ((percentileCut != "") and (canList.__len__() >=1)):
#     #Carry out thresholding on all entries in canList
#     for entry in canList:
#         candidateObject=candidateList(verboseSwitch)
#         if muteSpinner:
#             candidateObject.spinnerVerboseMode=False        
#         candidateObject.__loadfileQuick__(entry)
#         perCut=percentileCut.strip('"')
#         perCutLIST=str(perCut).split(',')
#         for singleThreshold in perCutLIST:
#             argProp,argVal=singleThreshold.split(':')
#             candidateResults=candidateObject.applyPercentageThreshold(argProp,argVal)
#             singleThresholdName=str(singleThreshold).replace('(','--').replace(')','--')
#             pathName=''
#             if (outfile != "") and (canList.__len__() == 1):
#                 #Prepend Threshold label to outfile
#                 if perCutLIST.__len__()>1:
#                     outfile=singleThreshold+"_PercentileCut_"+outfile
#                 candidateResults.writefile(outfile)
#                 if verboseSwitch:
#                     print "Wrote file :",outfile
#                 if dumpSummaryDisk:
#                     candidateResults.writeSummary(outfile)
#                 if glitchdb:
#                     candidateResults.writeGlitchDatabase('',outfile)
#                 if (dqList>=0):
#                     candidateObject.writeDQtable(dqList,outfile)
#             else:
#                 pathName=os.path.dirname(entry)
#                 saveFiles=pathName+'/Threshold:'+str(singleThresholdName)+':'+os.path.basename(entry)
#                 candidateResults.writefile(saveFiles)
#                 if verboseSwitch:
#                     print "Wrote file :",saveFiles
#             	if dumpSummaryDisk:
#                 	candidateResults.writeSummary(saveFiles)
#             	if glitchdb:
#                 	candidateResults.writeGlitchDatabase('',saveFiles)
#                 if (dqList>=0):
#                     candidateObject.writeDQtable(dqList,outfile)
#                 candidateResults.__setfilename__(saveFiles)
#                 if createHistogram:
#                     wrapperShowHistogram(candidateResults)
#                 if linePlot:
#                     wrapperGraphTriggers(candidateResults)
#             del candidateResults
#         del entry
#     del candidateObject
    
    #SECTION TO DUMP SUMMARY TO DISK WHEN NO THRESHOLDING REQD
    #Want to rely on traitSummary to speed things up if available.
elif ((canList.__len__() >=1) and (dumpSummaryDisk or glitchdb) and (expThreshold =="") ):
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        if muteSpinner:
            candidateObject.spinnerVerboseMode=False        
        if verboseSwitch:
            print "Processing:",entry
        candidateObject.summaryRead(entry)
        if ((outfile != "") and dumpSummaryDisk):
            candidateObject.writeSummary(outfile)
        elif dumpSummaryDisk:
            candidateObject.writeSummary(entry)
        if ((outfile != "") and glitchdb):
            candidateObject.writeGlitchDatabase('',outfile)
        elif glitchdb:
            candidateObject.writeGlitchDatabase('',entry)
        if ((outfile != "") and (dqList>=0)):
            candidateObject.writeDQtable(dqList,outfile)
        elif (dqList>=0):
            candidateObject.writeDQtable(dqList,entry)
        del candidateObject
    #SECTION TO DUMP SUMMARY TO SCREEN        
elif ((canList.__len__() >=1) and dumpSummaryScreen):
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        if muteSpinner:
            candidateObject.spinnerVerboseMode=False
        candidateObject.summaryRead(entry)
        candidateObject.printSummary()
        del entry

elif ((canList.__len__() >=1) and dqList >= 0):
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        candidateObject.summaryRead(entry)
        if ((outfile != "") and (canList.__len__() == 1)):
            candidateObject.writeDQtable(dqList,outfile)
        else:
            candidateObject.writeDQtable(dqList)

    #IF there are no files to work on found.    
elif (canList.__len__() < 0):
      print "It appears there are no files to process!"

    #THIS SECTION SHOULD NEVER HAPPEN
else:
    print "Error with combination of arguments given!"
    print "Legitimate command line arguments don't get this error!"
    print "Options in parser data structure."
    print options
    print "Corresponding argument values."
    print args
    print "Candidates files found :",canList.__len__()
    os.abort()
 
