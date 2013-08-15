An introduction to ULTRACAM and ULTRASPEC
=========================================

ULTRACAM and ULTRASPEC are high-speed astronomical cameras designed to acquire
low noise, time-series images of the sky in order to study time-variable 
astronomical objects. They use 1024x1024 CCDs. ULTRACAM uses three normal CCDs
(e2v 4720s) to observed in three wavebands simultaneously. ULTRASPEC has
single CCD but it has both normal and avalanche-gain readouts, the latter
giving better performance under some conditions.

ULTRACAM
--------

ULTRACAM has been running since May 2002 and has had over 350 nights on
telescopes. An archive of around 9TB has been built up. Its CCDs have readout
noise levels around 3 electrons. It has been on the 4.2m William Herschel
Telescope, 3.5m New Technology Telescope and the 8m Very Large Telescope.

The exposure tims of astronomical CCDs have to be slow to keep a low noise
level on readout. This limits the rate a which pixels can be read. To keep the
noise down on ULTRACAM, we read out at a rate of about 200,000 pixels/second,
thus a full read takes 6 seconds. One of the best ways to speed things is
therefore to read out sub-sections of the CCD which we call "windows". Most 
ULTRACAM data is of this windowed form. The windows of all three CCDs are
identical.

ULTRACAM uses dichroics to split the light into different wavebands. Typically
the CCDs work in the ultraviolet (370nm), green (470nm) and red(670nm). These
are often referred to as blue, green and red.

ULTRASPEC
---------

ULTRASPEC, despite its name (which is a hold-over from an original use on a
spectrogaph), is also an imaging camera. It has one CCD and will be used
on the 2.4m Thai National Telescope. Unlike ULTRACAM, it does not have 
dichroics which means that it can make use of "white light" filters to accept
more spectrum.
