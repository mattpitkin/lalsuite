from __future__ import division
"""
coh_PTF_post_processing.in - post processing driver script for coh_PTF externally triggered pipeline

This script generated the condor DAG necessary to run the post processing for the coh_PTF GRB pipeline.
"""

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__date__ = '$Date$'
__version__ = '$Revision$'

# =============================================================================
# import standard modules
# =============================================================================

import os,sys,ConfigParser,glob,shutil
from glue import pipeline,lal
from optparse import OptionParser
import tempfile

# =============================================================================
# Initialise subDAG
# =============================================================================

def init_subDAG(tag, outdir, logdir, cp):

  # generate job
  dag  = pipeline.CondorDAG('%s/%s_%s.log' % (logdir,tag,\
                                           tempfile.mktemp(dir='',prefix='')))
  dag.set_dag_file(os.path.join(outdir,tag))
  # set max jobs args
  if cp.has_option('condor-max-jobs', tag):
    dag.add_maxjobs_category(tag, cp.getint('condor-max-jobs', tag))

  return dag

# =============================================================================
# Write arguments for CondorJobs
# =============================================================================

def init_job(exe, universe, tag, subdir, logdir, cp, memory=None):

  logtag = '$(cluster)-$(process)'

  subdir = os.path.abspath(subdir)
  logdir = os.path.abspath(logdir)
  nfslogdir = '%s/logs' % subdir
  if not os.path.isdir(nfslogdir):
    os.mkdir(nfslogdir)

  job = pipeline.CondorDAGJob(universe,exe)
  job.set_sub_file(os.path.join(subdir,'%s.sub' % (tag)))
  job.set_stderr_file(os.path.join(nfslogdir,'%s-%s.err' % (tag,logtag)))
  job.set_stdout_file(os.path.join(nfslogdir,'%s-%s.out' % (tag,logtag)))
  job.add_condor_cmd('getenv','True')

  if cp.has_section(tag):
    job.add_ini_opts(cp, tag)
  if memory:
    job.add_condor_cmd('requirements', 'memory > %s' % memory)

  return job

# =============================================================================
# Write arguments for trig_combiner job
# =============================================================================

def trig_combiner_setup(job, category, ifotag, usertag, grb, onoffcache,\
                        grbdir, numtrials, outdir, timeslidecache=None,\
                        slidetag = None):

  # setup node
  node = pipeline.CondorDAGNode(job)
  node.set_category(category)

  node.add_var_opt('ifo-tag', ifotag)
  node.add_var_opt('user-tag', usertag)
  if slidetag:
    node.add_var_opt('slide-tag',slidetag)
  node.add_var_opt('grb-name', grb)
  node.add_var_opt('segment-dir', os.path.abspath(grbdir))
  if onoffcache:
    node.add_var_opt('cache', os.path.abspath(onoffcache))
  node.add_var_opt('num-trials', numtrials)
  node.add_var_opt('output-dir', os.path.abspath(outdir))
  if timeslidecache:
    node.add_var_opt('slide-cache',timeslidecache)

  return node

# =============================================================================
# Write arguments for trig_cluster job
# =============================================================================

def trig_cluster_setup(job, category, trigfile, outdir):

  node = pipeline.CondorDAGNode(job)
  node.set_category(category)
  node.add_var_opt('trig-file', trigfile)
  node.add_var_opt('output-dir', os.path.abspath(outdir))

  return node

# =============================================================================
# Write arguments for injfind job
# =============================================================================

def injfind_setup(job, category, injdir, injrun, ifotag, grb,\
                  datastart, dataduration):

  # construct cache file
  injcachefile = '%s/HL-INJECTION_GRB%s_%s-%s-%s.cache'\
                 % (injdir, grb, injrun, datastart, dataduration)
  # if cache does not exist, make one
  if not os.path.isfile(injcachefile):
    injfiles = glob.glob('%s/HL-INJECTION_GRB%s_%s_*-%s-%s.xml'\
                         % (injdir, grb, injrun, datastart, dataduration))
    injcache = lal.Cache.from_urls(injfiles)
    injcache.tofile(open(injcachefile, 'w'))

  # construct trigger cache file
  trigcachefile = '%s/%s-INSPIRAL_HIPE_GRB%s_%s-%s-%s.cache'\
                  % (injdir, ifotag, grb, injrun, datastart, dataduration)

  # initialise node
  node = pipeline.CondorDAGNode(job)
  node.set_category(category)
  node.add_var_opt('cache', trigcachefile)
  node.add_var_opt('inj-cache', injcachefile)

  return node

# =============================================================================
# Write arguments for injcombiner job
# =============================================================================

def injcombiner_setup(job, category, outdir, fmcache, injpattern, inclination):

  node = pipeline.CondorDAGNode(job)
  node.set_category(category)
  node.add_var_opt('inj-cache', fmcache)
  node.add_var_opt('output-dir', outdir)
  node.add_var_opt('inj-string', injpattern)
  node.add_var_opt('max-inclination', inclination)

  return node

# =============================================================================
# write sbv_plotter args
# =============================================================================

def sbv_setup(job, category, trigfile, grb, outdir, grbdir, vetodir=None,\
              injfile=None):

  node = pipeline.CondorDAGNode(job)
  node.set_category(category)
  node.add_var_opt('trig-file', trigfile)
  node.add_var_opt('grb-name', grb)
  node.add_var_opt('segment-dir', grbdir)
  node.add_var_opt('output-path', os.path.abspath(outdir))
  if vetodir:
    node.add_var_opt('veto-directory',vetodir)
  if injfile:
    node.add_var_arg('--inj-file %s ' % injfile)

  return node

# =============================================================================
# write efficiency args
# =============================================================================

def onoff_efficiency_setup(job, category, outdir, segdir, offsource, onsource,\
                           vetodir=None):

  node = pipeline.CondorDAGNode(job)
  node.set_category(category)
  node.add_var_opt('output-path', outdir)
  node.add_var_opt('offsource-file', offsource)
  node.add_var_opt('onsource-file', onsource)
  node.add_var_opt('segment-dir', segdir)
  if vetodir:
    node.add_var_opt('veto-directory',vetodir)

  return node

def injection_efficiency_setup(job, category, outdir, segdir, offsource,\
                               onsource, injrun, cp, found, missed,\
                               vetodir=None): 

  node = pipeline.CondorDAGNode(job)
  node.set_category(category)
  node.add_var_opt('output-path', outdir)
  node.add_var_opt('upper-inj-dist', cp.getfloat(injrun, 'max-distance')/1000.)
  node.add_var_opt('lower-inj-dist', cp.getfloat(injrun, 'min-distance')/1000.)
  node.add_var_opt('offsource-file', offsource)
  node.add_var_opt('onsource-file', onsource)
  node.add_var_opt('found-file', found)
  node.add_var_opt('missed-file', missed)
  node.add_var_opt('segment-dir', segdir)
  if vetodir:
    node.add_var_opt('veto-directory',vetodir)

  return node

# =============================================================================
# Setup for horizon distance plot
# =============================================================================

def horizon_distance_setup(job, category, ifotag, grb, onoffcache,\
                        grbdir, outdir):

  # setup node
  node = pipeline.CondorDAGNode(job)
  node.set_category(category)

  node.add_var_opt('ifo-tag', ifotag)
  node.add_var_opt('grb-xml', "%s/grb%s.xml" %(os.path.abspath(grbdir),grb))
  node.add_var_opt('cache', os.path.abspath(onoffcache))
  node.add_var_opt('output-dir', os.path.abspath(outdir))

  return node

# =============================================================================
# Finalise DAG object
# =============================================================================

def finalise_DAG(dag, parents=[]):

  dag.write_sub_files()
  dag.write_dag()
  dag.write_script()
  dagfile    = os.path.split(dag.get_dag_file())
  DAGManJob  = pipeline.CondorDAGManJob(dagfile[1],dagfile[0])
  DAGManNode = pipeline.CondorDAGManNode(DAGManJob)
  for node in parents:
    DAGManNode.add_parent(node)

  return DAGManNode

# =============================================================================
# parse command line
# =============================================================================

def parse_command_line():

  usage = """usage: %prog [options]

lalapps_coh_PTF_post_processing will set up a DAG to run post processing for the coh_PTF triggered pipeline end to end. It will set up six sub DAGs for each of the following jobs:

1) trig_combiner: combines triggers from onoff coh_PTF_inspiral jobs and
   separates into 'ONSOURCE' 'OFF_TRIAL_X' etc

2) trig_cluster: clusters each of the files from trig_combiner

3) injfinder: tests found/missed status of all injections

4) injcombiner: combines GRB injection sets based on inclination angles

5) sbv_plotter: calculates signal-based-vetoes and final detection statistic,
   and plots a whole bunch of stuff

6) efficiency: calculates inverse false-alarm rate (IFAR) for all injections
   and detection candidate events, and calculates search 90% efficiency
   distances. This is split into two submit files, separating onoff jobs and
   injection jobs.

The default is to run the full post processing pipeline, various steps can be skipp by specifying the --skip options. The requires arguments are:

--run-dir
--grb-name
--config-file
"""

  parser = OptionParser(usage=usage, version= "%prog CVS\n$Id$\n$Name$\n")

  parser.add_option( "--verbose", action="store_true", default=False,\
                    help="verbose output")

  parser.add_option("-r", "--run-dir", type="string", action="store",\
                    default=None,\
                    help="run directory, parent of the GRBXXXXXX directory")

  parser.add_option("-o", "--output-dir", type="string", action="store",\
                    default=os.getcwd(),\
                    help="output directory for post processing, "+\
                         "default: %default")

  parser.add_option("-V", "--veto-directory", type="string", action="store",\
                    default=None, help="Directory of veto files.")

  parser.add_option("-i", "--ifo-tag", type="string", action="store",\
                    default="H1L1V1", help="IFO tag, default: %default")

  parser.add_option("-n", "--grb-name", type="string", action="store",\
                    default=None, help="GRB identifier, e.g. 090802")

  parser.add_option("-f", "--config-file", type="string", action="store",\
                    default=None, help="post processing inifile for analsis")

  parser.add_option("-a", "--inj-config-file", type="string", action="store",\
                    default=None, help="injection inifile for analysis")

  parser.add_option("-p", "--log-path", action="store", type="string",\
                    default=None, help="directory to write condor log file")

  parser.add_option("-T", "--skip-trig-combiner", action="store_true",\
                    default=False,\
                    help="skip trig combiner, default: %default")

  parser.add_option("-C", "--skip-clustering", action="store_true",\
                    default=False, help="skip clustering, default: %default")

  parser.add_option("-I", "--skip-injfind", action="store_true",\
                    default=False,\
                    help="skip injection finding, default: %default")

  parser.add_option("-J", "--skip-injcombiner", action="store_true",\
                    default=False,\
                    help="skip injection combining by inclination angle, "+\
                         "default: %default")

  parser.add_option("-S", "--skip-sbv-plotter", action="store_true",\
                    default=False, help="skip SBV plotter, default: %default")

  parser.add_option("-E", "--skip-efficiency", action="store_true",\
                    default=False, help="skip efficiency, default: %default") 

  parser.add_option("-Z", "--skip-horizon-plot", action="store_true",\
                    default=False, help="skip horizon plot, default: %default")

  parser.add_option("-t", "--time-slides", action="store_true",\
                    default=False, help="Include time slides.")

  (opts, args) = parser.parse_args()

  if not opts.run_dir:
    parser.error('Must give --run-dir')

  if not opts.grb_name:
    parser.error('Must give --grb-name') 

  if not opts.config_file:
    parser.error('Must give --config-file')

  if not opts.inj_config_file:
    print >>sys.stdout,\
         "Injection config file not given. Running with no injections."
    opts.skip_injfind=True
    opts.skip_injcombiner=True
  #  parser.error('Must give --inj-config-file')

  return opts, args

# =============================================================================
# Main function
# =============================================================================

def main(rundir, outdir, ifotag, grb, inifile, injfile, verbose=False,\
         logdir=None, vetoDir=None,run_combiner=True, run_clustering=True,\
         run_injfind=True, run_sbvplotter=True, run_efficiency=True,\
         run_injcombiner=True,run_horizon_dist_plot=True,timeSlides=False):

  # load ini files
  if verbose:
    print >>sys.stdout
    print >>sys.stdout, 'Initialising post processing driver, '+\
                        'loading configuration files...'

  # get directory
  grbdir = os.path.abspath('%s/GRB%s' % (rundir, grb))
  if not os.path.isdir(grbdir):
    raise ValueError, 'Cannot find directory GRB%s in %s' % (grb, rundir)

  # generate post processing directory
  if not os.path.isdir(outdir):
    os.makedirs(outdir)
  # generat subdirectories
  os.chdir(outdir)
  plotdir = 'output'
  exedir = 'executables'
  if not logdir:
    logdir = '%s/%s' % (outdir, 'logs')
  for d in [plotdir, exedir, logdir]:
    if not os.path.isdir(d):
      os.mkdir(d)

  # load ini file
  cp = ConfigParser.ConfigParser()
  cp.optionxform = str
  cp.read(inifile)

  # Add veto directory if needed
  if vetoDir:
    if run_sbvplotter:
      cp.set('sbv_plotter','veto-directory',vetoDir)
    if run_efficiency:
      cp.set('efficiency','veto-directory',vetoDir)
      cp.set('injection-efficiency','veto-directory',vetoDir)

  # load inj file
  if injfile:
    injcp = ConfigParser.ConfigParser()
    injcp.optionxform = str
    injcp.read(injfile)

    # find injection runs
    injruns = injcp.sections()
  else:
    injruns = []

  usertag = cp.get('input', 'user-tag')

  # =========
  # get times
  # =========

  # get times from datafind cache
  datafindstr = '%s/datafind/%s-INSPIRAL_HIPE_GRB%s_DATAFIND-*-*.cache'\
                % (grbdir, ifotag, grb)
  datafindglob = glob.glob(datafindstr)
  if len(datafindglob)!=1:
    raise ValueError, 'Cannot find single datafind cache matching %s' % datafindstr
  datafindcache = datafindglob[0]

  datastart, dataduration = map(int, os.path.splitext(datafindcache)[0]\
                                          .split('-')[-2:])
  
  pad = cp.getint('data', 'pad-data') 
  start = datastart+pad
  duration = dataduration-2*pad

  # ================
  # find onoff and timeslides cache
  # ================

  onoffcache = '%s/onoff/' % grbdir
  onoffcache += '%s-INSPIRAL_HIPE_GRB%s_ZERO_LAG_CATEGORY_1-%s-%s.cache'\
                % (ifotag, grb, datastart, dataduration)
  # Get the appropriate tag
  zlCache = lal.Cache.fromfile(open(onoffcache, 'r'))
  zlCache = zlCache.sieve(description=usertag)
  fileDescription = zlCache[0].description.split('_')
  grbIndex = fileDescription.index('GRB%s' %(grb)) + 1
  zlString = '_'.join(fileDescription[grbIndex:])

  if timeSlides:
    timeSlidesCache = '%s/timeslides/' % grbdir
    timeSlidesCache += \
                '%s-INSPIRAL_HIPE_GRB%s_TIME_SLIDES_CATEGORY_1-%s-%s.cache'\
                % (ifotag, grb, datastart, dataduration)
    # Count and store how many long slides they are, and their names

    # Open cache and remove the TMPLTBANK file
    slideCache = lal.Cache.fromfile(open(timeSlidesCache, 'r'))
    slideCache = slideCache.sieve(description=usertag)

    # Identify and store unique desciriptions
    slideStrings = [zlString]
    for entry in slideCache:
      fileDescription = entry.description.split('_')
      # Remove user tag and split bank number
      grbIndex = fileDescription.index('GRB%s' %(grb)) + 1
      redDescription = '_'.join(fileDescription[grbIndex:])
      if redDescription not in slideStrings:
        slideStrings.append(redDescription)
    numLongSlides = len(slideStrings)
 
  else:
    timeSlidesCache = None

  # ==============
  # set parameters
  # ==============

  universe = cp.get('condor', 'universe')

  for (job, executable) in cp.items('condor'):
    if job=='universe':  continue
    # replace tilde in executable
    executable = os.path.expanduser(executable)
    # replace environment variables in executable
    executable = os.path.expandvars(executable)
    # copy executable to exedir
    executable2 = os.path.join(outdir, exedir, os.path.basename(executable))
    if not os.path.isfile(executable2) or\
       not os.path.samefile(executable, executable2):
      shutil.copy(executable, executable2)
    cp.set('condor', job, executable2)
 
  dirsList = ['trig_combiner','efficiency','injection-efficiency']
  if timeSlides:
    dirsList.append('trig_combiner_part2')

  for dir in dirsList:
    cp.set(dir,'pad-data',cp.get('data','pad-data'))
    cp.set(dir,'segment-length',cp.get('data','segment-length'))

  # ==========
  # write dags
  # ==========

  if verbose:
    print >>sys.stdout
    print >>sys.stdout, "Generating dag..."

  # initialise uberdag
  dagtag  = os.path.splitext(os.path.basename(inifile))[0]
  uberdag = pipeline.CondorDAG("%s/%s_uberdag.log" % (logdir,dagtag))
  uberdag.set_dag_file('%s_uberdag' % (dagtag))

  DAGManNode = {}

  # ==================
  # generate time tags
  # ==================

  numtrials = cp.getint('input', 'num-trials')
  timetags = [ 'ALL_TIMES', 'ONSOURCE', 'OFFSOURCE' ] +\
             [ 'OFFTRIAL_%d' % t for t in xrange(1,numtrials+1) ]

  if timeSlides:
    timeSlideTags = [ 'ALL_TIMES', 'OFFSOURCE' ]

  minmemory = 0
  if cp.has_option('pipeline', 'minimum-memory'):
    minmemory = cp.getfloat('pipeline', 'minimum-memory')

  # =============
  # trig_combiner
  # =============

  if run_combiner:

    tag    = 'trig_combiner'
    exe    = cp.get('condor', tag)
  
    # generate job
    dag = init_subDAG(tag, outdir, logdir, cp)
    job = init_job(exe, universe, tag, outdir, logdir, cp, minmemory)

    # setup onoff node
    if timeSlides:
      node = trig_combiner_setup(job, tag, ifotag, usertag, grb, None,\
                                 grbdir, numtrials, outdir, slidetag=zlString,\
                                 timeslidecache=onoffcache)
    else:
      node = trig_combiner_setup(job, tag, ifotag, usertag, grb, onoffcache,\
                                 grbdir, numtrials, outdir)

    dag.add_node(node)

    if timeSlides:
      # setup long slide nodes
      for slideTag in slideStrings:
        if slideTag == zlString:
          continue
        node = trig_combiner_setup(job, tag, ifotag, usertag, grb, None,\
                                   grbdir, numtrials, outdir,slidetag=slideTag,\
                                   timeslidecache=timeSlidesCache,)
        dag.add_node(node)
  
    # finalise DAG
    DAGManNode[tag] = finalise_DAG(dag)
    uberdag.add_node(DAGManNode[tag])

  # ============
  # trig_cluster
  # ============

  trigfile = {}
  clstfile = {}
  for timetype in timetags:
    trigfile[timetype] = '%s-%s_GRB%s_%s-%s-%s.xml.gz'\
                         % (ifotag, usertag, grb, timetype, start, duration)
    clstfile[timetype] = trigfile[timetype].replace(timetype,\
                                                     '%s_CLUSTERED' % timetype)

  if timeSlides:
    trigslidefile = {}
    clstslidefile = {}
    fileList = {}
    for timetype in timeSlideTags:
      fileList[timetype] = []
      for slideString in slideStrings:
        trigslidefile[timetype+slideString] = \
            '%s-%s_TIMESLIDES_GRB%s_%s_%s-%s-%s.xml.gz'\
            % (ifotag,usertag,grb,slideString,timetype,start,duration)
        clstslidefile[timetype+slideString] = \
            trigslidefile[timetype+slideString].replace(\
            timetype,'%s_CLUSTERED' % timetype)
        fileList[timetype].append(clstslidefile[timetype+slideString])

  if run_clustering:

    tag = 'trig_cluster'
    exe    = cp.get('condor', tag)
  
    # generate job
    dag = init_subDAG(tag, outdir, logdir, cp)
    job = init_job(exe, universe, tag, outdir, logdir, cp, minmemory)

    # loop time tags
    for timetype in timetags:
      node = trig_cluster_setup(job, tag, trigfile[timetype], outdir)
      dag.add_node(node)

    if timeSlides:
      for timetype in timeSlideTags:
        for slideString in slideStrings:
          node = trig_cluster_setup(job, tag,\
              trigslidefile[timetype+slideString], outdir)
          dag.add_node(node)

    # finalise DAG
    parents = []
    if run_combiner: parents.append(DAGManNode['trig_combiner'])
    DAGManNode[tag] = finalise_DAG(dag, parents)
    uberdag.add_node(DAGManNode[tag])

  # =========
  # trigcombiner (part 2)
  # =========

  # Only needed if doing time slides and output files are split up at this point
  # This only combines clustered triggers, as there are too many unclustered

  if timeSlides:
    tag    = 'trig_combiner_part2'
    exe    = cp.get('condor', 'trig_combiner')
  
    # generate job
    dag = init_subDAG(tag, outdir, logdir, cp)
    job = init_job(exe, universe, tag, outdir, logdir, cp, minmemory)

    # Create cache file
    currTag = 'COMBINED_CLUSTERED_COHPTF_%s' %(usertag)
    cache_fname = '%s-%s_GRB%s_%s-%s-%s.cache'\
        % (ifotag,currTag,grb,'ALL_TIMES',start,duration)
    output_data_cache = lal.Cache.from_urls(fileList['ALL_TIMES'])
    output_data_cache.tofile(open(cache_fname, "w"))

    for timetype in timeSlideTags:
      clstslidefile[timetype] = \
            '%s-%s_TIMESLIDES_GRB%s_%s-%s-%s.xml.gz'\
            % (ifotag,usertag,grb,timetype,\
               start,duration)

    node = trig_combiner_setup(job, tag, ifotag, usertag, grb, None,\
                               grbdir, numtrials, outdir,\
                               timeslidecache=cache_fname)
  
    dag.add_node(node)

    # finalise DAG
    parents = []
    if run_combiner: parents.append(DAGManNode['trig_cluster'])
    DAGManNode[tag] = finalise_DAG(dag,parents)
    uberdag.add_node(DAGManNode[tag])


  # =========
  # injfinder
  # =========

  if run_injfind:

    # find buffer segments
    buffseg = '%s/%s' % (grbdir, 'bufferSeg.txt')
    if not os.path.isfile(buffseg):
      raise ValueError, 'Cannot find buffer segment file as %s' % buffseg
  
    tag = 'injfinder'
    exe = cp.get('condor', tag)
  
    # generate job
    dag = init_subDAG(tag, outdir, logdir, cp)
    job = init_job(exe, universe, tag, outdir, logdir, cp, minmemory)

    # construct arguments
    job.add_opt('output-dir', outdir)
    job.add_opt('exclude-segments', buffseg)

    for injrun in injruns:
      node = injfind_setup(job, tag, '%s/%s' %(grbdir,injrun), injrun,\
                           ifotag, grb, datastart, dataduration)
      dag.add_node(node)

    # finalise DAG
    parents = []
    DAGManNode[tag] = finalise_DAG(dag, parents)
    uberdag.add_node(DAGManNode[tag])

  # ===========
  # injcombiner
  # ===========

  tag = 'injcombiner'
  # get necessary configuration variables
  injpatterns = cp.get('%s-meta' % tag, 'injection-patterns').split(',')
  inclinations = map(int, cp.get('%s-meta' % tag,\
                                 'injection-inclinations').split(','))

  if run_injcombiner:

    exe = cp.get('condor', tag)

    # generate job
    dag = init_subDAG(tag, outdir, logdir, cp)
    job = init_job(exe, universe, tag, outdir, logdir, cp, minmemory)

    # generate cache of FOUND/MISSED files
    fmcachefile = '%s/%s-GRB%s_FOUNDMISSED-%s-%s.cache'\
                  % (outdir, ifotag, grb, start, duration)
    fmfiles = ['%s/%s-INJECTION_GRB%s_%s_FOUND-%s-%s.xml'\
               % (outdir, ifotag, grb, injrun, start, duration)\
               for injrun in injruns] + \
              ['%s/%s-INJECTION_GRB%s_%s_MISSED-%s-%s.xml'\
               % (outdir, ifotag, grb, injrun, start, duration)\
               for injrun in injruns]
    fmcache = lal.Cache.from_urls(fmfiles)
    fmcache.tofile(open(fmcachefile, 'w'))

    for injpattern in injpatterns:
      if injpattern =='NONE':
        continue
      for inc in inclinations:
        node = injcombiner_setup(job, tag, outdir, fmcachefile, injpattern, inc)
        dag.add_node(node)

    # finalise DAG
    parents = []
    if run_injfind: parents.append(DAGManNode['injfinder'])
    DAGManNode[tag] = finalise_DAG(dag, parents)
    uberdag.add_node(DAGManNode[tag])

  # ===========
  # sbv_plotter
  # ===========

  filteredInjRuns = []
  filteredDistRuns = {}
  # Append runs that have not been split
  if injruns:
    for run in injruns:
      append = True
      for injpattern in injpatterns:
        if injpattern in run:
          append = False
          break
      if append:
        filteredInjRuns.append(run)
        filteredDistRuns[run]=run

    # and append the runs that have been split
    for injpattern in injpatterns:
      if injpattern == 'NONE':
        continue
      distRun = None
      for run in injruns:
        if injpattern in run:
          distRun = run
          break
      if not distRun:
        raise BrokenError, "Cannot find any injections matching %s in ini file"\
                           % (injpattern)
      for inc in inclinations:
        injrun = 'injectionsAstro%s_FILTERED_%d' % (injpattern,inc)
        filteredInjRuns.append(injrun)
        filteredDistRuns[injrun] = distRun


  if run_sbvplotter:

    tag = 'sbv_plotter'
    exe = cp.get('condor', tag)

    # generate job
    dag = init_subDAG(tag, outdir, logdir, cp)
    job = init_job(exe, universe, tag, outdir, logdir, cp, minmemory)

    for timetype in [ 'ALL_TIMES', 'OFFSOURCE' ]:
      # setup SBV plotting job
      sbvoutpath = '%s/%s/%s/plots' % (outdir, plotdir, timetype)
      node = sbv_setup(job, tag, trigfile[timetype], grb, sbvoutpath, grbdir,\
                       vetodir=rundir)
      dag.add_node(node)
      # setup SBV clustered plotting job
      sbvoutpath = '%s/%s/%s/plots_clustered' % (outdir, plotdir, timetype)
      node = sbv_setup(job, tag, clstfile[timetype], grb, sbvoutpath, grbdir,\
                       vetodir=rundir)
      dag.add_node(node)

    if timeSlides:
      for timetype in timeSlideTags:
        # setup SBV clustered plotting job
        sbvoutpath = '%s/%s/%s_slides/plots_clustered'\
                     % (outdir, plotdir, timetype)
        node = sbv_setup(job, tag, clstslidefile[timetype], grb, sbvoutpath,\
                         grbdir, vetodir=rundir)
        dag.add_node(node)

    # Make the injection plot with whatever background is available
    if timeSlides:
      injclstfile = clstslidefile['OFFSOURCE']
    else:
      injclstfile = clstfile['OFFSOURCE']
    
    for injrun in filteredInjRuns: 
      # setup SBV clusterd injection plots
      injfile = '%s/%s-INJECTION_GRB%s_%s_FOUND-%d-%d.xml'\
                % (outdir, ifotag, grb, injrun, start, duration)
      sbvoutpath = '%s/%s/%s/plots_clustered' % (outdir, plotdir, injrun)
      node = sbv_setup(job, tag, injclstfile, grb, sbvoutpath, grbdir,\
                       vetodir=rundir,injfile=injfile)
      dag.add_node(node)

    # finalise DAG
    parents = []
    if run_clustering: parents.append(DAGManNode['trig_cluster'])
    if run_injcombiner:    parents.append(DAGManNode['injcombiner'])
    if timeSlides: parents.append(DAGManNode['trig_combiner_part2'])
    DAGManNode[tag] = finalise_DAG(dag, parents)
    uberdag.add_node(DAGManNode[tag])
 
  # ==========
  # efficiency
  # ==========

  if run_efficiency:

    tag = 'efficiency'
    exe = cp.get('condor', tag)

    # generate job
    dag = init_subDAG(tag, outdir, logdir, cp)
    job = init_job(exe, universe, tag, outdir, logdir, cp, minmemory)

    # construct arguments
    job.add_opt('segment-dir', grbdir)

    # Select the appropriate offsource file
    if timeSlides:
      offfile = clstslidefile['OFFSOURCE']
    else:
      offfile = clstfile['OFFSOURCE']


    for timetype in timetags:
      if timetype in [ 'OFFSOURCE', 'ALL_TIMES' ]:  continue
      # setup onoff efficiency jobs
      effoutdir = '%s/%s/%s/efficiency' % (outdir, plotdir, timetype)
      node = onoff_efficiency_setup(job, tag, effoutdir, grbdir,\
                                    offfile, clstfile[timetype],\
                                    vetodir=rundir)
      dag.add_node(node)

    # define new job for inejction efficiency
    tag = 'injection-efficiency'
    job = init_job(exe, universe, tag, outdir, logdir, cp, minmemory)

    for timetype in timetags:
      if timetype in [ 'OFFSOURCE', 'ALL_TIMES' ]:  continue

      for injrun in filteredInjRuns:
        # setup injection efficiency jobs
        found = '%s/%s-INJECTION_GRB%s_%s_FOUND-%d-%d.xml'\
                % (outdir, ifotag, grb, injrun, start, duration)
        missed = found.replace('FOUND', 'MISSED')
        effoutdir = '%s/%s/%s/efficiency_%s' % \
                    (outdir, plotdir, injrun, timetype)
        node = injection_efficiency_setup(job, tag, effoutdir, grbdir,\
                           offfile,clstfile[timetype],\
                           filteredDistRuns[injrun], injcp,found, missed,\
                           vetodir=rundir)
        dag.add_node(node)

    # finalise DAG
    parents = []
    if run_clustering: parents.append(DAGManNode['trig_cluster'])
    if run_injcombiner:    parents.append(DAGManNode['injcombiner'])
    if timeSlides: parents.append(DAGManNode['trig_combiner_part2'])
    DAGManNode[tag] = finalise_DAG(dag, parents)
    uberdag.add_node(DAGManNode[tag])

  # =================
  # Horizon dist plot
  # =================

  if run_horizon_dist_plot:

    tag    = 'horizon_dist'
    exe    = cp.get('condor', tag)

    # generate job
    dag = init_subDAG(tag, outdir, logdir, cp)         
    job = init_job(exe, universe, tag, outdir, logdir, cp, minmemory)

    # setup single node
    node = horizon_distance_setup(job, tag, ifotag, grb, onoffcache,\
                               rundir, '%s/%s' %(outdir,plotdir))
    dag.add_node(node)

    # finalise DAG
    parents = []
    if run_combiner: parents.append(DAGManNode['trig_combiner'])
    DAGManNode[tag] = finalise_DAG(dag,parents=parents)
    uberdag.add_node(DAGManNode[tag])
    
  # =============
  # write uberdag
  # =============

  uberdag.write_sub_files()
  uberdag.write_dag()

  # print message
  print >>sys.stdout
  print >>sys.stdout, '------------------------------------'
  print >>sys.stdout, 'Ready. To submit, run:'
  print >>sys.stdout
  subcmd = 'condor_submit_dag '
  if cp.has_option('pipeline', 'maxjobs'):
    subcmd += '-maxjobs %s ' % cp.getint('pipeline', 'maxjobs')
  subcmd += os.path.abspath(uberdag.get_dag_file())
  print >>sys.stdout, subcmd
  print >>sys.stdout

  print >>sys.stdout, 'Once submitted, to monitor status, run:'
  print >>sys.stdout
  print >>sys.stdout, 'lalapps_ihope_status --dag-file %s'\
                      % (os.path.abspath(uberdag.get_dag_file()))
  print >>sys.stdout, '------------------------------------'
  print >>sys.stdout


if __name__=='__main__':

  opts, args = parse_command_line()
  verbose = opts.verbose
  rundir  = os.path.abspath(opts.run_dir)
  outdir  = os.path.abspath(opts.output_dir)
  logdir  = opts.log_path
  ifotag  = opts.ifo_tag
  grb     = opts.grb_name
  inifile = os.path.abspath(opts.config_file)
  injfile = opts.inj_config_file
  if injfile:
    injfile = os.path.abspath(injfile)
  vetoDir = None
  if opts.veto_directory:
    vetoDir = opts.veto_directory

  main(rundir, outdir, ifotag, grb, inifile, injfile, verbose=verbose,\
       logdir=logdir,vetoDir=vetoDir,\
       run_combiner=not opts.skip_trig_combiner,\
       run_clustering=not opts.skip_clustering,\
       run_injfind=not opts.skip_injfind,\
       run_sbvplotter=not opts.skip_sbv_plotter,\
       run_efficiency=not opts.skip_efficiency,\
       run_injcombiner=not opts.skip_injcombiner,\
       run_horizon_dist_plot = not opts.skip_horizon_plot,\
       timeSlides = opts.time_slides)
