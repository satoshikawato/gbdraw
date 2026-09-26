// Presentation only: generation equivalence and operation facts belong to config.
export const describeGenerationApplication = ({ status, unknown = [], error, operations }) => {
  const messages = {
    clean: ['Applied', 'The generation settings match the current Result.'],
    pending: ['Pending', 'Generation settings have changes for the next successful Generate. The current Result keeps its applied settings.'],
    unknown: ['Unknown', 'Some settings cannot be compared with the current Result. A successful Generate establishes the applied settings.'],
    invalid: ['Invalid settings', 'Check the generation inputs before generating. A Result is replaced only after a successful Generate.'],
    ungenerated: ['Not generated', 'There is no Result yet. Generate creates a Result from the current settings.']
  };
  const [label, message] = messages[status];
  const liveMessage = operations.liveApplying
    ? 'Live edit applying: geometry is being updated. Direct edits already applied are kept.'
    : operations.liveError
      ? 'Live edit failed: direct edits already applied are kept; geometry may still need updating. Retry the live edit or use Generate.'
      : '';
  const unknownMessage = status === 'pending' && unknown.length
    ? 'Some settings also have unknown application status.' : '';
  const generationMessage = operations.generating
    ? operations.cancelRequested ? 'Canceling Generate. The current Result is kept.'
      : 'Generating. The current Result stays until a successful replacement.'
    : '';
  return {
    label, message, unknownMessage, liveMessage, generationMessage,
    error: status === 'invalid' ? error : '',
    // No paths, counts, or changing validation details in the live region.
    // Equal messages produce no repeated DOM announcement on each keystroke.
    announcement: [label, message, unknownMessage, generationMessage, liveMessage].filter(Boolean).join(' ')
  };
};
