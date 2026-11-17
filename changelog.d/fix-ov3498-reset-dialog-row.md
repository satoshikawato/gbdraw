- Fixed: the **Reset fill color** dialog now asks by the feature's Legend row,
  the row its **Reset all** choice resets: it names that row, counts the row's
  other features, and opens whenever the row draws another feature. Before, it
  named and counted the features with the same label, so it could offer to
  reset two features while **Reset all** reset four, ask for a feature alone
  in its row, or not ask for a feature whose row draws others. A feature alone
  in its Legend row resets at once, as **This feature only**.
