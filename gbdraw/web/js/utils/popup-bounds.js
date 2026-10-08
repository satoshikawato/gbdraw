// @ts-check

// The lowest y a popup or floating panel may reach: the top of the app footer
// (`[data-app-footer]`) when it is on screen, else the viewport bottom. The
// feature and match popups, their drags and resizes, their first placement,
// and the alignment palette read this one limit.
/** @returns {number} */
export const popupViewportBottom = () => {
  const viewportBottom = Math.max(1, window.innerHeight || 1);
  const footerTop = document.querySelector('[data-app-footer]')?.getBoundingClientRect().top;
  return typeof footerTop === 'number' && footerTop > 0 && footerTop < viewportBottom
    ? footerTop
    : viewportBottom;
};
