- Fixed: Web app: after you delete a Legend row, press **Generate Diagram**, and then remove or recaption its color
  rule or give its feature another Legend name, **Restore** no longer lists a row the figure does not draw, and the
  next **Generate Diagram** no longer fails with "The generated result could not be accepted", on Circular and
  Linear diagrams (OV-2412, OV-3851). The 0.14.5 known issue said that a row without a stroke does not fail; a row
  with a Legend color of its own failed too, and a row without either was listed until the next Generate. A
  **Restore** pressed while a Generate or a diagram update runs now applies when it finishes, so it lists the rows
  that update draws.
