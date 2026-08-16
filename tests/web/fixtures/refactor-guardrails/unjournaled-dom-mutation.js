export const bypassMutationJournal = (element) => {
  element['setAttribute']('fill', '#abcdef');
  return true;
};
