- Fixed: Web app: a Generate with LOSAT execution set to Threaded, whose threaded engine could not start or keep
  its threads running in the browser (a thread that did not start in time, a worker that could not be made), now says
  so and names Serial and Auto, instead of an unrecognized failure. Auto still continues with the serial engine
  (OV-3316).
- Fixed: Web app: when LOSAT stops a comparison search and says why (an option it refuses), the Generation Error
  panel no longer offers Retry Generate, which would stop the same way; Show in Comparison now opens the LOSAT
  options in Comparison Settings instead of the Execution setting (OV-3317).
