export const bypassWithAliasClone = (input) => {
  const borrowedOwner = input.features;
  const copy = cloneJsonValue;
  return copy(borrowedOwner, []);
};
