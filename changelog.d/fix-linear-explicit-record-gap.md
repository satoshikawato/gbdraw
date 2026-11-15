- Changed: The Linear record gap now defaults to Auto. `LinearLayout.record_gap` and
  `LinearMultiRecordOptions.record_gap_px` default to `None` instead of `24.0`, and `--linear_record_gap` has no
  default value; code that reads these defaults as numbers (such as `LinearLayout().record_gap + 10`) must allow
  `None`. Auto draws 24 px, narrowed in a row whose gaps would take more than 25% of its width. A gap you set is
  drawn as set; when a row cannot hold it, the diagram is refused with the reason. The Web app's **Record gap (px)**
  field is blank (**Auto**) by default, and a Session saved with a gap keeps it.
