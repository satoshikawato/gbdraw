export const bypassMutationJournalProperty = (element) => {
  element.style.fill = '#abcdef';
  return true;
};
