export const bypassWithPlaceholderMutation = (previewRuntime) => previewRuntime.commitDomEdit({
  mutate: function reportMutationWithoutOwningIt() {
    return true;
  }
});
