export const bypassAsyncCommitToken = async (previewRuntime, gate) => {
  const targetToken = previewRuntime.captureDomEditToken();
  await gate;
  if (!previewRuntime.isDomEditTokenCurrent(targetToken)) return false;
  return previewRuntime.commitDomEdit({
    mutate: ({ mutation, svg }) => mutation.setAttribute(svg, 'data-stale', 'true')
  });
};
