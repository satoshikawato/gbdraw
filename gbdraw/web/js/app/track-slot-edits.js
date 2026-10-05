// R10, Q3, R13: a track stack editor routes each of its edits of a feature-slot
// input through `changeTrackLayout`, the feature placement owner's transition
// that the composition root injects as a port. A DOM event the template passes
// last names the control for Cancel and focus (its listener's element, not an
// icon inside a button); the edit itself runs without the event.
export const featureSlotEdits = (changeTrackLayout, actions) => Object.fromEntries(
  Object.entries(actions).map(([name, action]) => [name, (...args) => {
    const event = typeof Event === 'function' && args.at(-1) instanceof Event ? args.pop() : null;
    return changeTrackLayout(() => action(...args), event?.currentTarget);
  }])
);
