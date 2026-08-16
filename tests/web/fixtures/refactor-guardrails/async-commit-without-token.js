export const bypassAsyncTargetToken = async (previewRuntime, gate) => {
  await gate;
  return previewRuntime.commitDomEdit({
    mutate: ({ mutation, svg }) => mutation.setAttribute(svg, 'data-stale', 'true')
  });
};
