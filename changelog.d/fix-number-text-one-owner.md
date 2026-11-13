- Fixed: Web app: in a browser set to a language other than English, some numbers followed the browser's format
  (4.641.652 in German, ৪৬,৪১,৬৫২ in Bengali) instead of 4,641,652: the sequence search's query length and hit
  Start–End, a feature's Length in its popup, the span lengths of a comparison match's sequences, the "Coordinates
  exceed the sequence length" message, the count of selected features, and the feature and span lengths of an
  exported interactive SVG. They now use the English digits and grouping in every language, as every other number
  does. The sequence search's "No hits" note also groups the query length (1,300 bp, not 1300 bp) (OV-2777).
