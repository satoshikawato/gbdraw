SURVEY: feature-label placement on circular maps (9 tools, plus section 10 on plasmid-map conventions). Researched 2026-10-12. "fetched" = I read the document or source in this session. "snippet" = text from a WebSearch result only, because the page returned 403. "recollection" = my memory, not checked.

## 1. Geneious (Prime), circular sequence view
Placement model: I found no Geneious page that describes the circular-view label algorithm (could not verify). The manual's viewer controls give an "Annotations > Labels" option with Inside / Outside / Inside or Outside / None. It is not stated whether this option applies to the circular view. (recollection) The circular view draws labels horizontally outside the ring, with short leader lines.
When labels do not fit: the manual has a "Hide excessive labels" checkbox: "removing annotation labels which are too frequent". It lists the option without saying which view it governs or how labels are chosen.
User controls: Labels placement drop-down, "Hide excessive labels", label font size under Sizes, and the annotation-types filter. The manual says the circular overview is unaffected by these controls except the annotation types.
Sources:
- https://manual.geneious.com/en/latest/Sequences.html (fetched; sequence viewer controls and the circular-sequences section)
- https://help.geneious.com/hc/en-us/articles/360045071831 (403, not read)
- Several UGENE circular-view pages came up in search. They describe a different product and are not Geneious evidence.

## 2. SnapGene, plasmid map
Placement model: a circular map has external labels (outside the ring) plus internal labels (inside). The SnapGene 2.8.0 release notes say Map view can show "either external plus internal feature labels, or only internal feature labels, or no feature labels". The 5.1.0 notes say "long map labels" can span two lines in circular Map view. (recollection) External labels sit on short connector lines, horizontal, and are auto-arranged.
When labels do not fit: from a search-result snippet of the support article "How do I ensure that feature labels are always displayed below or inside a DNA map?": "If feature labels are crowded, SnapGene may display some of them above a linear map or outside a circular map". A preference forces them inside, which "may cause some of the feature labels to be hidden". Labels that do not fit are hidden. To see more, reduce the font size or enlarge the window.
User controls: the Map-view label button (all / internal only / none), the DNA preference "No feature names outside or above a map except point feature names", font size, and window size. 5.1.0 also fixed the centre "map label" so it never overlaps features.
Sources:
- https://support.snapgene.com/hc/en-us/articles/10383722725524-Display-Feature-Labels-Below-or-Inside-a-Map (snippet; fetch returned 403)
- https://www.snapgene.com/updates/snapgene-2-8-0-release-notes (fetched)
- https://www.snapgene.com/updates/snapgene-5-1-0-release-notes (fetched)

## 3. Proksee / CGView.js (read from source)
Placement model: labels are horizontal text boxes outside the ring, joined to the feature by a straight leader line. The line is perpendicular (radial) in 'default' placement and angled in 'angled' placement. Each label gets a clock position (1 to 12) for where the line meets its box. Defaults: labelLineLength 20 px, priorityMax 50, font monospace 12.
When labels do not fit: in `Annotation.draw()`, labels are sorted by priority (favorite first, then longer feature first). The first priorityMax (50) get full collision-free placement. The rest get a plain box at the default line length and are dropped if they overlap any already accepted rect. Labels never shrink. The set is recomputed at each zoom level, so zooming in reveals more labels.
User controls: `annotation.labelPlacement` ('default' or 'angled'; separate Fast and Full variants used during interaction), `labelLineLength`, `priorityMax`, `onlyDrawFavorites`, per-feature `favorite`, per-label `lineAttachment`, `font`, `color`.
Algorithm, 'default' (`LabelPlacementDefault.placeLabels`):
1. For each label, in priority order, find the point at radius `outerOffset + margins + lineLength` at the label's bp angle.
2. Build the label rect from its measured text width and the font height, anchored at that point by its clock attachment.
3. If the rect overlaps any rect already placed (`Rect.overlap`), add `label.height` to lineLength and retry (do/while). Labels stack outward in radial tiers on the same ray. Nothing is dropped here.
4. Store the rect and the line attachment point.
Algorithm, 'angled' (`LabelPlacementAngled.js`, about 860 lines):
1. Sort labels by bp. Compute each label's default rect.
2. Group consecutive overlapping labels into "islands" (not if more than half the sequence apart).
3. Place each island by fanning label lines away from each other. Boundary labels go to a maximum angle (80 degrees of line tilt, converted to a bp shift; half that for 6 and 12 o'clock labels). Labels fill inward from both ends of the island.
4. Labels that cannot fit without exceeding the max angle are "popped". They are spread evenly over the bp range between the neighbouring non-popped labels and their line length grows until nothing clashes.
5. Adjacent islands that still clash are merged and re-placed, repeating up to the number of islands. A final pass lengthens lines for any remaining overlap.
Sources:
- https://raw.githubusercontent.com/stothard-group/cgview-js/main/src/LabelPlacementDefault.js (fetched)
- https://raw.githubusercontent.com/stothard-group/cgview-js/main/src/LabelPlacementAngled.js (fetched)
- https://raw.githubusercontent.com/stothard-group/cgview-js/main/src/Annotation.js (fetched)
- https://raw.githubusercontent.com/stothard-group/cgview-js/main/src/Label.js (fetched)

## 3b. Java CGView (older)
Placement model: horizontal labels outside the backbone circle with connector lines. Reverse-strand labels can go inside the circle (`useInnerLabels`).
When labels do not fit: `labelsToKeep` (default 1000) randomly discards labels before placement if there are more. `moveInnerLabelsToOuter` (default true) moves inner labels outside when the inside has no room. `labelPlacementQuality` (good / better / best, default better) trades placement quality for speed. The page does not describe the algorithm.
User controls: `labelLineLength` (initial connector length), `labelLineThickness`, `labelPlacementQuality`, `labelsToKeep`, `useInnerLabels` (true / false / auto; auto means inside only on a zoomed map), `moveInnerLabelsToOuter`, font.
Sources:
- https://sites.ualberta.ca/~stothard/cgview/xml_cgview.html (fetched)

## 4. OGDRAW (OrganellarGenomeDRAW)
Placement model: not documented in the 2019 NAR paper, the 2013 NAR paper, or the FAQ/help pages I fetched (could not verify). (recollection) Gene boxes sit in the ring, inside or outside by strand. Names are horizontal text next to the boxes. I am not sure whether lines connect them.
When labels do not fit: the FAQ says over-long gene names can "overlap with the GC graph or the image border". The fix is the "tidy up" option (shortens or reformats names) or editing the GenBank names. There is no collision avoidance in the docs. The FAQ also says names can be duplicated or overlap when the annotation has duplicate features. Users can export PS and edit by hand.
User controls: "tidy up", "Label intron-containing genes with *", output format and resolution, custom feature-class config XML, a FAQ item on label tags (tag priority: gene, label, note, locus_tag, product, operon). No label-placement options.
Sources:
- https://chlorobox.mpimp-golm.mpg.de/OGDraw-FAQ.html (fetched)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC6602502/ (fetched; mentions only tidy-up and the intron asterisk)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC3692101/ (fetched; nothing on label placement)

## 5. DNAPlotter / Artemis
Placement model: labels are opt-in per feature. The 2009 paper says "Labels and arrows can be added to a feature". The Artemis circular code `Block.java` (the circular-plot package I could read; I did not find DNAPlotter's own repository) draws the label as unrotated horizontal text. It is placed inside the ring at the feature's mid-angle, offset from the stroke.
When labels do not fit: the `Block.java` label code has no overlap test, so labels collide. Features can be dragged and edited by hand.
User controls: per-feature "Show Label" checkbox and label text in the feature dialog, the track manager (filter, size, position), and a list of labels or gene names can be pasted into a track filter.
Sources:
- https://pmc.ncbi.nlm.nih.gov/articles/PMC2612626/ (fetched; DNAPlotter paper)
- https://raw.githubusercontent.com/sanger-pathogens/Artemis/master/src/main/java/uk/ac/sanger/artemis/circular/Block.java (fetched; `if(drawLabel)` block)

## 6. Circos
Placement model: a text track is an annulus between r0 and r1. Each label is placed at the angular position of its data point, left-aligned at the track radius and rotated so the baseline is radial. `label_rotate=no` turns that off. Labels point outward only. Optional link lines (`show_links`, `link_dims` with 5 segments) connect the label to its true position, because snuggling can shift labels.
When labels do not fit: labels are tiled so they do not overlap. A label that would overlap is drawn at the same angle but radially shifted out. With `label_snuggle=yes` the label may also be shifted in the angular direction to use fewer layers. The docs say snuggling is "heuristic, not based on any kind of global optimization". A label that cannot fit inside r0..r1 is not drawn. `-debug_group textplace` lists the `not_placed` ones. The fix is to increase r1. With hundreds of labels, optimization can take minutes.
User controls: `label_snuggle`, `max_snuggle_distance` (e.g. 2r, in pixels or relative to label size), `snuggle_sampling` (pixels between tested positions), `snuggle_tolerance`, `snuggle_link_overlap_test` and `snuggle_link_overlap_tolerance`, `snuggle_refine` (swaps label pairs whose links cross, only with links shown), `padding` and `rpadding`, `r0` and `r1`, `label_size`, `label_font`.
Sources:
- https://circos.ca/documentation/tutorials/2d_tracks/text_1/ (fetched with `curl -k`; the site certificate failed WebFetch)
- https://circos.ca/documentation/tutorials/2d_tracks/text_2/ (fetched with `curl -k`)

## 7. pyCirclize
Placement model: `track.annotate(x, label, ...)` draws a radial line from `min_r` (default the track outer edge) to `max_r` (default `min_r + 5`) and a horizontal text label at its end, with `label_size` 8 and `shorten=20` characters by default. `track.text()` is separate. It places text at a position with `adjust_rotation=True` and `orientation` "horizontal" or "vertical", using a rotation chosen from the angle.
When labels do not fit: `adjust_annotation()` runs once at plot time in `circos.py`, for annotations only. The docstring calls the algorithm "experimental". Steps:
1. Sort annotations in four quadrants in reading order.
2. For each one, test it against the bboxes (expanded by 1.2 x 1.2) of all earlier ones. If it overlaps and sits in the side-facing angle range (30 to 150 or 210 to 330 degrees), try angular shifts within a window around its angle, nearest first, up to `max_rad_shift` (3 degrees) and the neighbours' angles, in 0.1 degree steps.
3. If it still overlaps, step the label outward by `dr` = 0.1 radius units per iteration.
4. Gives up at `max_iter` (1000). With more than `limit` (200) annotations, no adjustment is done and a warning is issued, so labels overlap.
User controls: `config.ann_adjust` fields `enable`, `limit`, `max_iter`, `drad`, `dr`, `expand`, `max_rad_shift`. Per-call `min_r`, `max_r`, `label_size`, `shorten`, `line_kws`, `text_kws`. Note that `annotate()` has no `annotation_adjust` argument; adjustment is global config.
Sources:
- https://raw.githubusercontent.com/moshi4/pyCirclize/main/src/pycirclize/annotation.py (fetched)
- https://raw.githubusercontent.com/moshi4/pyCirclize/main/src/pycirclize/config.py (fetched)
- https://raw.githubusercontent.com/moshi4/pyCirclize/main/src/pycirclize/track.py (fetched)
- https://raw.githubusercontent.com/moshi4/pyCirclize/main/src/pycirclize/circos.py (fetched; call site `adjust_annotation(ax)`)

## 8. Benchling, plasmid map
Placement model: could not verify. The help-centre pages returned 403 to every fetch. (recollection) Labels are drawn outside the ring, horizontal, with short connector lines. I am not sure of the exact behaviour.
When labels do not fit: could not verify. A search snippet says the 2023 release notes added a way to show overlapping parts on the plasmid map as named annotations. I could not read the entry, and it is about overlapping features rather than label collision.
User controls: could not verify.
Sources:
- https://help.benchling.com/hc/en-us/articles/9684225694093-Sequence-Annotation (snippet; 403; lists Sequence map, Linear map, Plasmid map)
- https://help.benchling.com/hc/en-us/articles/16691595907469-Release-Notes-Volume-5-2023 (snippet; 403)

## 9. Biopython GenomeDiagram
Placement model: `_CircularDrawer.get_feature_sigil` places each label as radial text at the feature's start, middle or end, chosen by `label_position`. It sits just outside the feature (`top`) for strand +1 and just inside (`btm`) for strand -1. `label_strand` can override which side. Text is flipped by 180 degrees on the left half so it reads upright, with text anchor "end" there where needed.
When labels do not fit: there is no collision avoidance and no dropping. Labels overlap each other and the neighbouring features. Users limit labels by setting `label=True` only on chosen features and by reducing `label_size`. The default `label_size` is 6 and the default font is Helvetica.
User controls: `label`, `label_position`, `label_strand`, `label_size`, `label_color`, `label_font`, `name_qualifiers`. `label_angle` is documented as linear only. The circular drawer makes text radial regardless.
Sources:
- https://raw.githubusercontent.com/biopython/biopython/master/Bio/Graphics/GenomeDiagram/_CircularDrawer.py (fetched; roughly lines 353 to 431)
- https://raw.githubusercontent.com/biopython/biopython/master/Bio/Graphics/GenomeDiagram/_Feature.py (fetched; the attribute list)

## Patterns across tools
(a) Radial or tangent text on one ring plus drop or overlap: GenomeDiagram (radial text, always overlaps), Circos (radial text, label dropped if it cannot fit the annulus, but stacks first), DNAPlotter (horizontal text inside, overlaps).
(b) Staggered tiers on the same ray: CGView.js 'default' (grow line length by one label height until the rect is free), Circos stacking (same angle, shifted out). Both extend the leader line.
(c) Labels pushed along the ring with leader lines: Circos snuggle (angular search within `max_snuggle_distance`, then refine crossings), CGView.js 'angled' (islands, max line angle, popped labels), pyCirclize `adjust_annotation` (angle shift up to 3 degrees, then radial push).
(d) Side-ordered callout columns (left and right, ordered leaders): none of the nine, from what I could verify. The closest is pyCirclize, which sorts per quadrant in reading order and uses horizontal text. Section 10 finds this layout in plasmid-map tools (CLC "Stacked", SimVector "Bracket View").
(e) Priority and zoom hiding: CGView.js (priority sort, only the first 50 get full placement, the rest dropped if they overlap, recomputed per zoom; Java CGView `labelsToKeep` discards randomly), Geneious "Hide excessive labels" (rule not documented), SnapGene (hides labels that do not fit; asks the user to shrink the font or enlarge the window).
(f) Inside versus outside as an option: SnapGene (external plus internal, or internal only), Geneious (Inside / Outside / Inside or Outside), Java CGView (`useInnerLabels`, with fallback to outside).
(g) User fixes it by hand or by shortening names: OGDRAW (tidy up, PS export), DNAPlotter (per-feature label, drag features), pyCirclize (global adjust limits).
Take-away for 100 to 200 labels: only CGView.js and Circos have documented, tunable algorithms that cope with that density. Both combine a leader-line stack with priority or bounded displacement. pyCirclize disables adjustment above 200 annotations by default.
Not verified: Geneious circular label behaviour, Benchling, OGDRAW placement, and how DNAPlotter's own code (not Artemis's `Block.java`) places labels.

## 10. Plasmid-map convention: cluster callout lists
Scope: the look the owner described is (a) neighbouring labels, typically restriction sites, grouped into a left-aligned vertical list beside the cluster; (b) leaders fanning to the list rows, or one shared leader or bracket for a boxed list such as an MCS; (c) gene names along or beside their arcs. I found documentation for (a) and (b) in three tools. I found none for (c). No tool I could read documents the grouping rule (distance threshold, row order, when a single leader is used).

### 10.1 SimVector (PREMIER Biosoft): "Bracket View" for the MCS
Placement model: the vendor page says "Vector catalogs generally display MCS's in a unique format, with enzymes one above the other in a vertical list" and calls SimVector's version the "Bracket View". SimVector finds the MCS from GenBank header annotations or the user specifies it manually. The page does not say how the bracket or leader is drawn.
When labels do not fit: the rule is a filter. It "list[s] only the enzymes that cut within the MCS and nowhere else", so the list is limited to unique cutters inside the MCS. Behaviour for a very long MCS list is not documented (could not verify).
User controls: manual MCS definition; enzyme sets from more than 20 commercial providers; EPS and SVG export for hand editing. Other options could not be verified.
Sources:
- https://premierbiosoft.com/plasmid_maps/VectorNTI.html (fetched; the MCS / Bracket View paragraph; the same text is on https://premierbiosoft.com/plasmid_maps/featuressv/cloning.html, fetched)

### 10.2 CLC Genomics Workbench (QIAGEN): "Stacked" versus "Radial" enzyme labels
Placement model: for restriction sites in the circular view the manual offers two label layouts. "Radial ... will place the restriction site labels as close to the cut site as possible." "Stacked ... will stack the labels so that all enzymes are shown. For circular views, it will align all the labels on each side of the circle. This can be useful for clearly seeing the order of the cut sites when they are located closely together." The manual adds that in a circular view, Stacked and Radial "also affect the layout of annotations".
When labels do not fit: Stacked shows every enzyme by moving labels into aligned columns on the left and right of the circle. The manual warns that the linear-view "Flag" option "will make it hard to see when several cut sites are located close to each other". Behaviour for a very long stack is not documented.
User controls: label mode No labels / Flag (linear) / Radial and Stacked (circular); "Show name flags"; per-enzyme list with colours; enzyme filters.
Caution: my search tool attributed this text to Geneious. It is CLC text. I did not find any Geneious page that documents the same layout.
Sources:
- https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/700/Dynamic_restriction_sites.html (fetched)

### 10.3 SnapGene and Addgene maps (SnapGene Server)
Placement model: the documentation I could read says little about layout. Addgene's announcement says SnapGene "positions restriction enzyme annotations throughout the map at locations better aligned to their precise positions within the plasmid". The 5.1.0 release notes (snippet) say enzyme sites and annotations on a linear map may overlap to reduce map height. A search summary of Addgene's SnapGene-generated PDF maps reports enzyme labels such as "PaeR7I - XhoI(3767)", where a dash joins isoschizomers cutting at one position. I did not read those PDFs myself (snippet). My memory of the vertical-list look for clustered sites in SnapGene and Addgene maps is recollection only.
When labels do not fit: the support pages that describe label hiding are for feature labels (section 2). For enzyme labels, the documented control is choosing fewer enzymes (enzyme sets such as "Unique 6+ Cutters", individual check boxes), not a layout rule. Could not verify any automatic grouping or MCS-bracket rule.
User controls: Show Enzymes button, Choose Enzymes (sets or manual), Show locations (coordinates) checkbox, font size.
Sources:
- https://blog.addgene.org/improved-plasmid-maps-powered-by-snapgene (fetched)
- https://support.snapgene.com/hc/en-us/articles/24818378518036-Show-or-Hide-All-Enzymes-Sites (snippet; 403)
- https://support.snapgene.com/hc/en-us/articles/10383760625172-Show-or-Hide-Restriction-Site-or-Primer-Coordinates-on-a-Map (snippet)
- https://media.addgene.org/data/plasmids/62/62665/62665-map_xjMk28tbER4p.pdf (snippet; an example SnapGene-generated map, not read)

### 10.4 Other tools asked about
- PlasMapper (2004 NAR): fetched. The map is drawn by the Java CGView. Feature names are placed "using an iterative collision detection and shifting process that results in a visually pleasing label arrangement and no label overlap". Restriction enzyme labels and feature labels can be switched on or off. The paper says SVG output lets users reposition labels by hand. No grouping or list rule. Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC441548/ (fetched).
- MacVector: fetched (tutorial PDF, text extracted). Preferences include "Enable Label Hiding" (speeds layout; "laying out a lot of labels, particularly on circular maps, can be fairly slow"), "Simplify layout for large sequences", and "Maximum Levels" (feature rows). No list or bracket layout described. Source: https://www.macvector.com/Resources/Creating%20Vector%20Maps.pdf (fetched).
- Older GCG-era PlasmidMap: fetched. Tick labels are printed next to a tick around the outside of the circle. The restrictions section says "The displays become very crowded if many labels are plotted." Source: https://rothlab.ucdavis.edu/genhelp/plasmidmap.html (fetched).
- Geneious: the Cloning manual page (https://manual.geneious.com/en/latest/Cloning.html, fetched) says restriction annotations show the recognition site, cut site and methylation sensitivity, and adds Enzymes and Fragments tabs. It says nothing about label layout (could not verify any circular grouping).
- GENETYX, pDRAW32, Vector NTI, SeqBuilder/Lasergene, ApE: could not verify. Searches returned no page describing their circular label layout. Vector NTI: the PREMIER Biosoft page only mentions display options in general. ApE: the Frontiers paper page was fetched but the search for label text found nothing useful.
- Addgene map viewer: it is SnapGene Server output (blog above); it has an interactive Map tab with hover and a sortable enzyme table that highlights a site on the map. No layout rule documented.

### 10.5 Points relevant to the design
1. The side-column layout exists in shipped tools (CLC "Stacked": labels aligned on each side of the circle to show cut-site order; SimVector "Bracket View": vertical enzyme list for the MCS). Both pair it with a filter or a mode switch, not with an automatic distance rule.
2. SimVector's rule is selection-based: the grouped list is only the unique cutters inside a user- or file-defined MCS region. The group is defined by an annotation, not by label proximity.
3. For long lists the documented answers are filtering (enzyme sets, unique cutters, "Hide excessive labels", "Enable Label Hiding"), not scrolling or shrinking.
4. Gene names along arcs (c): no documentation found (could not verify). GenomeDiagram's radial text at feature start, middle or end (section 9) is the only documented tangent-or-radial gene-name placement I read.
