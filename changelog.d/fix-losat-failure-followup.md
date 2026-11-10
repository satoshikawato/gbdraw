- Fixed: Web app: a Generate with LOSAT execution set to Threaded, whose threaded engine could not start or keep
  its threads running in the browser (a thread that did not start in time, a worker that could not be made or
  failed, shared memory that could not be allocated), now says so and names Serial and Auto, instead of an
  unrecognized failure; Sequence search and Find in other files say to click Search again. Auto still continues with
  the serial engine (OV-3316).
- Fixed: Web app: when LOSAT stops a comparison search and says why (an option it refuses), the Generation Error
  panel says that Generate with the same settings stops the same way and no longer offers Retry Generate. In a
  Linear diagram, Show in Comparison opens Comparison Settings at the LOSAT mode in use instead of the Execution
  setting; in a Circular diagram, which had no button for it, Show in Comparison now opens Pairwise Comparisons
  (OV-3317).
- Fixed: Web app: a comparison search in which LOSAT stops unexpectedly (a LOSAT crash, for example on a genetic code
  it does not know) now says so on both the serial and the threaded engine, with Show in Comparison and without
  Retry Generate, instead of an unrecognized failure (OV-3327).
