const REMOTE_ASSET_MANIFEST_PATH = '/gallery/remote-assets.json';
const REMOTE_GALLERY_PREFIXES = [
  '/gallery/examples/',
  '/gallery/sessions/',
  '/gallery/sources/',
  '/gallery/media/',
];
const GALLERY_ROOT_PATH = '/gallery/';
const PALETTE_EXPLORER_PATH = '/gallery/palettes';
const PALETTE_EXPLORER_ASSET_PATH = '/gallery/palettes/';
const GALLERY_STATIC_SEGMENTS = new Set([
  'examples',
  'files',
  'media',
  'palettes',
  'sessions',
  'sources',
  'thumbnails',
  'tutorials',
]);
const REMOTE_CACHE_CONTROL = 'public, max-age=86400';
const SVG_CONTENT_SECURITY_POLICY = [
  "default-src 'self' data: blob:",
  "script-src 'unsafe-inline'",
  "style-src 'unsafe-inline'",
  "img-src 'self' data: blob:",
  "frame-ancestors 'self'",
].join('; ');

let remoteAssetsPromise;

const isRemoteGalleryCandidate = (pathname) =>
  REMOTE_GALLERY_PREFIXES.some((prefix) => pathname.startsWith(prefix));

const isGalleryViewRoute = (pathname) => {
  const match = pathname.match(/^\/gallery\/([^/.]+)\/?$/);
  return Boolean(match && !GALLERY_STATIC_SEGMENTS.has(match[1].toLowerCase()));
};

const redirectToPath = (request, pathname) => {
  const redirectUrl = new URL(request.url);
  redirectUrl.pathname = pathname;
  redirectUrl.hash = '';
  return Response.redirect(redirectUrl.toString(), 308);
};

const fetchAssetAtPath = (request, env, pathname) => {
  const assetUrl = new URL(request.url);
  assetUrl.pathname = pathname;
  assetUrl.search = '';
  assetUrl.hash = '';
  return env.ASSETS.fetch(new Request(assetUrl, request));
};

const fetchGalleryShell = (request, env) => fetchAssetAtPath(request, env, GALLERY_ROOT_PATH);

const fetchPaletteExplorer = (request, env) =>
  fetchAssetAtPath(request, env, PALETTE_EXPLORER_ASSET_PATH);

const inferContentType = (pathname, upstreamHeaders) => {
  if (pathname.endsWith('.svg')) return 'image/svg+xml; charset=utf-8';
  if (pathname.endsWith('.json')) return 'application/json; charset=utf-8';
  if (pathname.endsWith('.webp')) return 'image/webp';
  if (pathname.endsWith('.mp4')) return 'video/mp4';
  if (pathname.endsWith('.webm')) return 'video/webm';
  if (pathname.endsWith('.ogg')) return 'video/ogg';
  return upstreamHeaders.get('content-type') || 'application/octet-stream';
};

const getRemoteAssets = async (request, env) => {
  if (!remoteAssetsPromise) {
    const manifestUrl = new URL(REMOTE_ASSET_MANIFEST_PATH, request.url);
    remoteAssetsPromise = env.ASSETS.fetch(new Request(manifestUrl, { method: 'GET' }))
      .then(async (response) => {
        if (!response.ok) return {};
        const payload = await response.json();
        return payload && typeof payload === 'object' && !Array.isArray(payload) ? payload : {};
      })
      .catch(() => ({}));
  }
  return remoteAssetsPromise;
};

const buildUpstreamHeaders = (request) => {
  const headers = new Headers();
  for (const name of ['accept', 'range', 'if-none-match', 'if-modified-since']) {
    const value = request.headers.get(name);
    if (value) headers.set(name, value);
  }
  return headers;
};

const fetchRemoteAsset = async (request, pathname, remoteUrl) => {
  const upstream = await fetch(
    new Request(remoteUrl, {
      method: request.method,
      headers: buildUpstreamHeaders(request),
    }),
    {
      cf: {
        cacheEverything: true,
        cacheTtl: 86400,
      },
    }
  );
  const headers = new Headers(upstream.headers);
  headers.set('Cache-Control', REMOTE_CACHE_CONTROL);
  headers.set('Content-Type', inferContentType(pathname, upstream.headers));
  headers.set('Cross-Origin-Embedder-Policy', 'require-corp');
  headers.set('Cross-Origin-Opener-Policy', 'same-origin');
  headers.set('Cross-Origin-Resource-Policy', 'same-origin');
  headers.delete('Content-Security-Policy-Report-Only');
  headers.delete('X-Frame-Options');

  if (pathname.endsWith('.svg')) {
    headers.set('Content-Security-Policy', SVG_CONTENT_SECURITY_POLICY);
  } else {
    headers.delete('Content-Security-Policy');
  }

  return new Response(upstream.body, {
    status: upstream.status,
    statusText: upstream.statusText,
    headers,
  });
};

export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    const isReadRequest = request.method === 'GET' || request.method === 'HEAD';
    if (isReadRequest && url.pathname === PALETTE_EXPLORER_PATH) {
      return redirectToPath(request, `${PALETTE_EXPLORER_PATH}/`);
    }
    if (isReadRequest && url.pathname === PALETTE_EXPLORER_ASSET_PATH) {
      return fetchPaletteExplorer(request, env);
    }
    if (isReadRequest && isGalleryViewRoute(url.pathname)) {
      if (url.pathname.endsWith('/')) {
        return redirectToPath(request, url.pathname.slice(0, -1));
      }
      return fetchGalleryShell(request, env);
    }
    if (isReadRequest && isRemoteGalleryCandidate(url.pathname)) {
      const remoteAssets = await getRemoteAssets(request, env);
      const remoteUrl = remoteAssets[url.pathname.replace(/^\/+/, '')];
      if (remoteUrl) {
        return fetchRemoteAsset(request, url.pathname, remoteUrl);
      }
    }
    return env.ASSETS.fetch(request);
  },
};
