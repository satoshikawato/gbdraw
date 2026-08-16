export const bypassAsyncRevalidation = async (previewRuntime, gate) => {
  const targetToken = previewRuntime.captureDomEditToken();
  await gate;
  return previewRuntime.commitDomEdit({
    targetToken,
    mutate: ({ mutation, svg }) => mutation.setAttribute(svg, 'data-stale', 'true')
  });
};
