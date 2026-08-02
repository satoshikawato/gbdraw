import assert from 'node:assert/strict';

globalThis.location = { href: 'https://example.test/gbdraw/web/' };

const capabilityModule = await import(
  '../../gbdraw/web/js/services/runtime-capabilities.js'
);
const {
  DIAGRAM_ENGINE_COMPATIBILITY_MESSAGE,
  DiagramRuntimeCompatibilityError,
  EXPECTED_WEB_RUNTIME_CAPABILITIES,
  validateWebRuntimeCapabilities
} = capabilityModule;

const cloneExpected = () =>
  JSON.parse(JSON.stringify(EXPECTED_WEB_RUNTIME_CAPABILITIES));

assert.equal(EXPECTED_WEB_RUNTIME_CAPABILITIES.rendering.optionSchema, 3);
assert.deepEqual(
  EXPECTED_WEB_RUNTIME_CAPABILITIES.rendering.featureRenderings,
  ['arrow', 'rectangle', 'underlay']
);

{
  const source = cloneExpected();
  const validated = validateWebRuntimeCapabilities(source);
  assert.deepEqual(validated, EXPECTED_WEB_RUNTIME_CAPABILITIES);
  assert.notEqual(validated, source);
  assert.equal(Object.isFrozen(validated), true);
  assert.equal(Object.isFrozen(validated.request), true);
  assert.equal(Object.isFrozen(validated.request.supportedSchemas), true);
}

for (const mutate of [
  (manifest) => {
    manifest.schema += 1;
  },
  (manifest) => {
    manifest.request.supportedSchemas.reverse();
  },
  (manifest) => {
    delete manifest.rendering.optionSchema;
  },
  (manifest) => {
    manifest.rendering.optionSchema = 1;
  },
  (manifest) => {
    manifest.rendering.featureRenderings = ['arrow', 'arrowhead', 'rectangle', 'underlay'];
  },
  (manifest) => {
    manifest.unexpected = true;
  }
]) {
  const manifest = cloneExpected();
  mutate(manifest);
  assert.throws(
    () => validateWebRuntimeCapabilities(manifest),
    (error) => {
      assert.equal(error instanceof DiagramRuntimeCompatibilityError, true);
      assert.equal(error.message, DIAGRAM_ENGINE_COMPATIBILITY_MESSAGE);
      assert.match(error.diagnostic, /^capabilities/);
      assert.doesNotMatch(error.message, /schema|protocol|wheel|deploy/i);
      return true;
    }
  );
}

class FakeWorker {
  static instances = [];
  static capabilities = cloneExpected();

  constructor() {
    this.listeners = new Map();
    this.messages = [];
    this.terminated = false;
    FakeWorker.instances.push(this);
  }

  addEventListener(type, listener) {
    const listeners = this.listeners.get(type) || new Set();
    listeners.add(listener);
    this.listeners.set(type, listeners);
  }

  removeEventListener(type, listener) {
    this.listeners.get(type)?.delete(listener);
  }

  postMessage(message) {
    this.messages.push(message);
    if (message.type !== 'init') return;
    queueMicrotask(() => {
      this.emit('message', {
        data: {
          id: message.id,
          type: 'init',
          ok: true,
          capabilities: structuredClone(FakeWorker.capabilities)
        }
      });
    });
  }

  emit(type, event) {
    this.listeners.get(type)?.forEach((listener) => listener(event));
  }

  terminate() {
    this.terminated = true;
  }
}

globalThis.Worker = FakeWorker;
const generationModule = await import(
  '../../gbdraw/web/js/services/diagram-generation.js'
);
const {
  disposeDiagramGenerationWorker,
  getDiagramGenerationWorkerCapabilities,
  preinitializeDiagramGenerationWorker
} = generationModule;

{
  const [first, second] = await Promise.all([
    preinitializeDiagramGenerationWorker(),
    preinitializeDiagramGenerationWorker()
  ]);
  assert.equal(FakeWorker.instances.length, 1);
  assert.equal(
    FakeWorker.instances[0].messages.filter((message) => message.type === 'init').length,
    1
  );
  assert.equal(first, second);
  assert.equal(first, getDiagramGenerationWorkerCapabilities());
  assert.equal(Object.isFrozen(first), true);
  disposeDiagramGenerationWorker();
}

{
  const incompatible = cloneExpected();
  incompatible.renderProtocol += 1;
  FakeWorker.capabilities = incompatible;
  await assert.rejects(
    preinitializeDiagramGenerationWorker(),
    (error) => {
      assert.equal(error instanceof DiagramRuntimeCompatibilityError, true);
      assert.equal(error.message, DIAGRAM_ENGINE_COMPATIBILITY_MESSAGE);
      assert.match(error.diagnostic, /renderProtocol/);
      return true;
    }
  );
  const worker = FakeWorker.instances.at(-1);
  assert.equal(worker.terminated, true);
  assert.equal(getDiagramGenerationWorkerCapabilities(), null);
}

{
  const oldWheel = cloneExpected();
  oldWheel.rendering.optionSchema = 2;
  oldWheel.rendering.featureRenderings = ['arrow', 'arrowhead', 'rectangle', 'underlay'];
  FakeWorker.capabilities = oldWheel;
  await assert.rejects(
    preinitializeDiagramGenerationWorker(),
    (error) => {
      assert.equal(error instanceof DiagramRuntimeCompatibilityError, true);
      assert.equal(error.message, DIAGRAM_ENGINE_COMPATIBILITY_MESSAGE);
      assert.match(error.diagnostic, /rendering\.(?:featureRenderings|optionSchema)/);
      return true;
    }
  );
  const worker = FakeWorker.instances.at(-1);
  assert.equal(worker.terminated, true);
  assert.equal(getDiagramGenerationWorkerCapabilities(), null);
}

console.log('runtime capability tests passed');
