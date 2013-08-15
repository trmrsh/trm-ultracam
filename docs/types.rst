ULTRACAM run classification
===========================

Before ULTRACAM data can be used, they must be \"reduced". In the simplest
case each exposure (up to 12 MB of data) may end up as just 32 bytes of data:
8 bytes for the time, plus three sets of 8 bytes to store the data for the
target star. The reduction is currently accomplished with a fairly hands-on
reduction \"pipeline" but use of these data would be greatly facilitated by
an automatic reduction of all data just to see what is present.  Unfortunately
this is essentially impossible at the moment because we do not have a reliable
classifcation of the various types of runs. I have come to the conclusion that
the only way forward is a classification by eye of all the runs. This
represents a considerable amount of effort, so the more people who carry it
out the better. The purpose of this document is to help relatively
inexperienced classifiers get started.

Key principles
--------------

The two most important principles of the classification are:

* *It is* **much** *better to be sure than to be fast!* i .e. if
  you are unsure of a run's type, mark it *unsure* for later checking, do not
  simply choose a category in order to move on. Take your time in
  deciding. This cannot be emphasized enough. *Unreliable classification is
  worse than none at all.*

* *If you notice something unusual, comment upon it.* You can type  
  a comment for each run. If you don't, to some extent the effort 
  of looking at the images is going to waste.

Data types
--------------------

ULTRACAM's basic purpose is to measure the amount of light coming from
astronomical objects. It does this by accumulating \"counts" in \"picture
elements" (hereafter pixels). Each pixel corresponds to a place on the sky
determined by the pointing of the telescope during the observation. ULTRACAM's
field of view is small, around 6 arcminutes or 1/5\ :sup:`th` of the
angular diameter of the Moon. A count is directly related to the number of
electrons captured in a pixels following the detection of photons (particles
of light) which generate electron/hole pairs in the detector. In ULTRACAM, the
relation is close to 1-to-1 counts vs electrons.

In addition to the \"science" data, there are a variety of other types,
including some that simply have to be classified as \"junk". The main purpose
of the classification is to decide whether a frame is science, junk etc. This
is facilitated by a script called *rchecker.py*. The complete list of
data types recognised by *rchecker.py* and their meanings is as follows:

* **acquisition**: when first moving to a target it usually takes at least a
  few frames, and sometimes more, to define its precise position and
  orientation. One often ends up with sequence of images in which the targets
  change position. These are useful for later working out where the target is,
  but not for scientific data because the flux will vary in unpredictable
  ways. Such frames should be classified as *acquisition*. Note that one
  should prefer if possible to classify data frames as *science* (see below),
  so *acquisition* should be reserved for runs where there is clearly
  something wrong compared to a classic science run. Usually this will be a
  change is position, possibly multiple times.


* **bias**: even a zero second exposure on a CCD will return \"counts"
  but these are an artefact of the electronics and known as the \"bias
  level". One of the first steps in data reduction is to remove this bias, in
  order to measure the true number of counts in each pixel. This requires
  measuring the bias level by taking short exposure frames with as little
  light getting to the detectors as possible. This means all lights off in the
  telescope dome, and trying to block the light path. A well-taken bias should
  have no obvious images on it. There are many such frames in the ULTRACAM
  archive. Only classify a run as *bias* if there are no evident problems
  with it. This is one type of frame where a certain level of automated
  checking is possible, but it is still important to make the classification.


* **dark**: if one exposes a CCD with zero light, counts still pile up on
  top of the bias to to thermal excitation. These are called \"dark
  counts". They have the nasty property of varying considerably from 
  pixel-to-pixel. To calibrate these, we need exposures with zero light. These
  are classified as *dark*. They are in fact very difficult to get right
  because it is hard to ensure a complete absence of light. Look for comments
  in the logs: \"dome closed", \"lights off", \"mirror in beam" etc, which
  show that the observers have been careful. A short exposure *dark* will
  look very similar to a *bias*. If in doubt, mark *unsure*, and make a
  comment about what you were unsure about.


* **flat**: the pixels of CCDs vary slightly from each other in
  sensitivity. To calibrate this one takes images of a uniformly lit
  scene. The usual \"object" of choice is the twilight blue sky. These
  produce images known as \"flat fields", with type *flat*. Ideally they 
  have a large number of counts to reduce the noise. The maximum number of
  counts = 2\ :sup:`16` - 1 = 65535, set by the 16-bits used to store 
  each pixel, however the ULTRACAM chips are not actually reliable up to 
  this level. Instead a phenomenon we call \"peppering" affects the 
  green and blue CCDs at levels of around 28,000 or more. 


* **flux**: to calibrate the overall system throughput, observations of
  stars called \"flux standards" are taken. These are constant stars whose
  brightnesses have been carefully measured by others. Flux standard runs are
  often marked as such, but otherwise they typically happen during the period
  immediately after and before twilight in the evening and morning, and
  involve a bright star moved into the right-hand half of the CCDs. 


* **junk**: data of **no** use at all. Examples are runs with all CCDs
  saturated, ones affected by readout problems (and not even of any use
  as examples), all zeroes, or just no stars visible but some light
  getting on to the CCDs. Please be absolutely sure before deciding to 
  call a run *junk* because it carries the implication that we would 
  lose nothing by deleting the data.


* **publicity**: \"pretty pictures" for publicity. Galaxies, nebulae and
  the like, which are not of interest from a time-variability perspective.


* **science**: the runs which actually produce what we are after. These
  should generally be as dull as possible -- long sequence of images, all in
  the same place. Sometimes the telescope will drift off position, and you
  might then see the objects move at the end of a run. These should still be
  classified as **science** if you think there is any worthwhile data there.
  If you see any translation or rotation during a run, **comment on it**!
  You may also see the objects disappear, perhaps because of clouds. If you do
  **comment on it**! Such comments will be invaluable for diagnosing
  problems with the automated pipeline.


* **technical**: a miscellaneous category to cover test data that does
  not fit into any other bin. e.g. noise tests where successively brighter
  frames are taken to calibrate the noise characteristics of the CCDs.


* **unsure**: this one is an easy one. If you can't decide on the data
  type, mark it as *unsure*. Write a comment to explain which categories
  you were trying to decide between.

Filters
-------

As well as the data type, you will need to decide upon the filters in use.
This is tricky. *rchecker.py* pulls filter names from the
`night logs 
<http://deneb.astro.warwick.ac.uk/phsaap/ultracam/logs/index.html>`_,
which would be great if they were 100% reliable. Unfortunately
they are not. Practically the only way to judge this is to look through the
`logs <http://deneb.astro.warwick.ac.uk/phsaap/ultracam/logs/index.html>`_,
of a night to see the conistency of the filter names. If there appear to
be changes of filter, look at the runs before and after, and see how long the
filter changes took. Any apparent filter change in less than about 3 minutes
is suspiciously quick.

Targets
-------

Usually the target name sugested (which is derived from the night logs) is 
correct, but keep an eye out for dramatic changes in the star pattern without
the target apparently having changed -- there might have been an error. Such
errors are important to track down because the target name is the key way in
which one tracks down runs in the archive.

Comments
--------

The final item of information required is a commnt. This is your chance to
point out anything unusual in the data. Try to avoid saying things that can be
deduced from the logs such as \"four windows" or \"fast readout". Similarly,
don't repeat what is said in the night logs. Useful things to comment on are

* any rotation of the field

* any large movement of the field

* any disappearance of the targetsm owing to cloud perhaps.

* very variable blurring due to the atmosphere (known as \"seeing")

* moving objects (meteors, satellites, asteroids)

* excessive noise patterns

* readout problems


