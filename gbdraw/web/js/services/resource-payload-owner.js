const payloadOwners = new WeakMap();

const objectKey = (value) => (
  value && (typeof value === 'object' || typeof value === 'function') ? value : null
);

export const setResourcePayloadOwner = (descriptor, owner) => {
  const descriptorKey = objectKey(descriptor);
  const ownerKey = objectKey(owner);
  if (descriptorKey && ownerKey) payloadOwners.set(descriptorKey, ownerKey);
  return descriptor;
};

export const getResourcePayloadOwner = (descriptor) => {
  const key = objectKey(descriptor);
  return (key && payloadOwners.get(key)) || key;
};
