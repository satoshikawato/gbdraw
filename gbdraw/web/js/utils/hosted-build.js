// @ts-check
// What the build that produced this page says about where it is hosted. The
// gbdraw.app bundle (`tools/prepare_cloudflare_pages.py`, which owns the fact
// and the meta name) loads Google Analytics and writes
// `<meta name="gbdraw-hosted-analytics" content="google-analytics-4">`; every
// other build (the package's offline GUI, a local server) has no such element.
/**
 * @typedef {object} HostedBuild
 * @property {boolean} googleAnalytics the page loads Google Analytics (the hosted site)
 */

/**
 * Reads the build flags once, at app setup.
 * @param {Pick<ParentNode, 'querySelector'>} root the page's document
 * @returns {Readonly<HostedBuild>}
 */
export const readHostedBuild = (root) => Object.freeze({
  googleAnalytics: root.querySelector('meta[name="gbdraw-hosted-analytics"]') !== null
});
