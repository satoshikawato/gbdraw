- Fixed: After the Undo of a feature color edit, a new Feature visibility rule
  hides its features in the preview at once, as Generate does; before, the
  preview kept drawing them until the next Generate.
- Fixed: Undo of an emptied number setting, such as Legend font size, shows
  the earlier value again and draws it; before, the field stayed empty and the
  figure kept the Auto size. Redo of the emptying leaves the field empty, and
  loading a Session or configuration file while a number field is empty
  applies its stored value.
- Fixed: After a Generate that brings an older Session's saved protein
  comparison results up to date, a new Feature visibility rule hides its
  features in the preview at once, as Generate does; before, the preview kept
  drawing them until the next Generate.
