/* cmdline.h */

/* File autogenerated by gengetopt version 2.11  */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
#define CMDLINE_PARSER_PACKAGE PACKAGE
#endif

#ifndef CMDLINE_PARSER_VERSION
#define CMDLINE_PARSER_VERSION VERSION
#endif

struct gengetopt_args_info
{
  char * source_name_arg;	/* Name of source (default='NONAME_SOURCE').  */
  double right_ascension_arg;	/* Right Ascension of source, in rad.  */
  double declination_arg;	/* Declination of source, in rad.  */
  double orientation_arg;	/* Orientation of source, in rad.  */
  char * detector_arg;	/* Detector name; legal names: lho, llo, virgo, geo, tama, cit, test.  */
  int start_time_sec_arg;	/* GPS seconds field of start time of observation.  */
  int start_time_nanosec_arg;	/* GPS nanoseconds field of start time of observation (default='0').  */
  int nsample_arg;	/* number of samples.  */
  double sampling_interval_arg;	/* sampling time interval, in seconds.  */
  int n_ra_arg;	/* Number of grid points in RA (default='256').  */
  int n_dec_arg;	/* Number of grid points in Dec (default='64').  */
  int start_ra_arg;	/* Starting RA index (>= 0) (default='0').  */
  int count_ra_arg;	/* Number of sub-grid points in RA.  */
  int start_dec_arg;	/* Starting Dec index (>= 0) (default='0').  */
  int count_dec_arg;	/* Number of sub-grid points in Dec.  */
  char * format_arg;	/* output format (default='mam').  */
  char * output_dir_arg;	/* Output directory (default='.').  */
  int verbosity_arg;	/* verbosity level for debugging (default='0').  */
  int debug_arg;	/* debug level (default='0').  */
  char * earth_ephemeris_arg;	/* Earth ephemeris file.  */
  char * sun_ephemeris_arg;	/* Sun ephemeris file.  */

  int help_given ;	/* Whether help was given.  */
  int version_given ;	/* Whether version was given.  */
  int single_source_given ;	/* Whether single-source was given.  */
  int whole_sky_given ;	/* Whether whole-sky was given.  */
  int snapshot_given ;	/* Whether snapshot was given.  */
  int timeseries_given ;	/* Whether timeseries was given.  */
  int average_given ;	/* Whether average was given.  */
  int source_name_given ;	/* Whether source-name was given.  */
  int right_ascension_given ;	/* Whether right-ascension was given.  */
  int declination_given ;	/* Whether declination was given.  */
  int orientation_given ;	/* Whether orientation was given.  */
  int detector_given ;	/* Whether detector was given.  */
  int start_time_sec_given ;	/* Whether start-time-sec was given.  */
  int start_time_nanosec_given ;	/* Whether start-time-nanosec was given.  */
  int nsample_given ;	/* Whether nsample was given.  */
  int sampling_interval_given ;	/* Whether sampling-interval was given.  */
  int n_ra_given ;	/* Whether n-ra was given.  */
  int n_dec_given ;	/* Whether n-dec was given.  */
  int start_ra_given ;	/* Whether start-ra was given.  */
  int count_ra_given ;	/* Whether count-ra was given.  */
  int start_dec_given ;	/* Whether start-dec was given.  */
  int count_dec_given ;	/* Whether count-dec was given.  */
  int format_given ;	/* Whether format was given.  */
  int output_dir_given ;	/* Whether output-dir was given.  */
  int verbosity_given ;	/* Whether verbosity was given.  */
  int debug_given ;	/* Whether debug was given.  */
  int earth_ephemeris_given ;	/* Whether earth-ephemeris was given.  */
  int sun_ephemeris_given ;	/* Whether sun-ephemeris was given.  */

} ;

int cmdline_parser (int argc, char * const *argv, struct gengetopt_args_info *args_info);

void cmdline_parser_print_help(void);
void cmdline_parser_print_version(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
