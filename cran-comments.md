## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

There was one NOTE that is found only on Windows Server 2022, R-devel 64 bit:

	* checking for detritus in the temp directory ... NOTE
	Found the following files/directories:
  	'lastMiKTeXException'

As noted in R-hub issue #503, this could be due to a bug/crash in MiKTeX and can likely be ignored.

As a general note, this is my first submission to CRAN.

## Downstream dependencies
There are no downstream dependencies yet