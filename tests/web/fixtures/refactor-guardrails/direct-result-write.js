export const bypassResultOwner = (state, index, content) => {
  const resultItems = state.results.value;
  resultItems[index] = { ...resultItems[index], content };
};
