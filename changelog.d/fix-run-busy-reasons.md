- Fixed: Web app: a Review alignment opened before a Generate stays open after it. While the Generate runs, its
  Apply waits with the reason shown; once the Generate draws another Result, Apply is disabled with
  "The diagram changed after this review opened. Start alignment again." instead of failing with "The operation
  failed without recognized diagnostic information".
- Fixed: Web app: Rotate record using this feature, refused while a Generate runs, shows the reason ("Generating
  diagram. Retry after generation finishes.") instead of "Record rotation failed.".
- Fixed: Web app: buttons disabled while a Generate runs (Apply now, Align…, Review alignment options…, the review's
  Apply, Align records, Reset record positions, Reset positions, Reset Settings) state why next to them or to
  screen readers.
- Fixed: Web app: Align on a Similarity Group with several candidates is refused while a Generate runs, as it is for
  one candidate, instead of opening a review.
- Fixed: Web app: Align… or Review alignment options… on a feature that is not exactly one member of its Similarity
  Group, on a Similarity Group the shown diagram no longer has, or whose final geometry could not be checked, says
  what to do instead of "The operation failed without recognized diagnostic information"; a failed one-click Align
  no longer leaves the alignment in progress (OV-3849).
