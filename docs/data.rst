ULTRACAM / ULTRASPEC data
=========================

It is helpful in dealing with ULTRACAM and ULTRASPEC to have some 
understanding of the nature and purposes of the various data files
you may encounter. There are two main types which are now described.

raw data
--------

ULTRACAM and ULTRASPEC data comes in two main forms. During observations the
"raw" data is generated, consisting of two files per run, one an xml file,
e.g. "run012.xml", and a corresponding binary data file, i.e.
"run012.dat". The xml files contain information about the windows set, the
readout mode used, along with whatever other information is defined at the
start of a run, while the data file contains all the timing and data frames.
The xml files are text files and can be printed to screen and even edited if
you take care. It is worth looking inside some to see what they contain. The
data files by contrast are efficient binary files are not easily understood.
During observing the data files are added to frame by frame, first some timing
bytes, then the data, 2 bytes per pixel. The ordering is not simple, and as
binary files, no standard software can know their structure. These files are
handled by the C++ pipeline (hereafter just the "pipeline"), and by
trm.ultracam.

ucm files
---------

The pipeline often needs to write and read single frames, e.g. to represent
flat-fields or bias frames. It does so in files with extension ".ucm". These
are also a specialised binary format. trm.ultracam can also read and write
these files. 

exporting data
--------------

Both the raw and ucm files can be exported to FITS for examination with other 
software.
