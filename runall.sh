#!/bin/bash
# Ruth Gómez Graciani
# 01 04 2020

###############################################################################
# Description:                                                                
# Run all the analyses in this project                                        
###############################################################################

# SAVE HISTORY
# Save date and command 
# =========================================================================== #
# echo "$(date +"%Y%m%d") : runall.sh $*" >> project/history.txt

# USAGE 
# Help message. Don't forget to mention new options in README.md!!  
# =========================================================================== #

usage()
{
  echo "Usage: 

  $(basename $0) <COMMAND> [OPTIONS]
  "
  echo "Commands:

  all [-s]                Run all commands.
  download [-s]           Download publicly available data.
  preprocess [-s]         Prepare raw data to use it in the analysis commands.
  "
  echo "Options:

    -h                    Show help.
    -s <screens>          Set a number of screens to use. 
  "

}

# PARSE OPTIONS
# Parse options and get help 
# =========================================================================== #

# Parse help message
while getopts "h" OPT; do
   case "" in
      h)  usage; exit 1 ;;
   esac
done
shift $((OPTIND -1))

# Save command, if any
COMMAND=$1; shift

# Set default optional variables
SCREENS=1;

# Set empty mandatory variables

# Parse command optons
case "$COMMAND" in
  #Run all commands
  all ) 
    while getopts "s:" OPTIONS ; do
      case "$OPTIONS" in
        s)  SCREENS=${OPTARG} ;;
      esac
    done
    shift $((OPTIND -1)) 
    ;;
  #Download public data
  download ) 
    while getopts "s:" OPTIONS ; do
      case "$OPTIONS" in
        s)  SCREENS=${OPTARG} ;;
      esac
    done
    shift $((OPTIND -1))
    ;;
  #Prepare raw data
  preprocess ) 
    while getopts "s:" OPTIONS ; do
      case "$OPTIONS" in
        s)  SCREENS=${OPTARG} ;;
      esac
    done
    shift $((OPTIND -1))
    ;;
esac

# # Check that empty mandatory variables are full
# if [ "$COMMAND" == "all"]; then
#   if [-z "${VAR}"] ||[] ; then ; ;fi
# elif []; then
# fi

# DOWNLOAD RAW DATA
# All publicly available data is now downloaded to data/raw.
# =========================================================================== #

# if [ "$COMMAND" == "all" ] || [ "$COMMAND" == "download" ]; then
# # Download file x
# # Origin, content and format in data/raw/README.md
# fi

# PREPROCESS RAW DATA
# Files common for any analysis are created only once and stored in data/use.
# =========================================================================== #

# if [ "$COMMAND" == "all" ] || [ "$COMMAND" == "preprocess" ]; then
# # Process file x to y
# # Origin, content and format in data/use/README.md
# fi


