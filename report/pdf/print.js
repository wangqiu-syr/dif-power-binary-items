// print.js: print the report's HTML to PDF with headless Chromium through Playwright.
//   node report/pdf/print.js <in.html> <out.pdf> [title]
// US Letter, the margins from style.css (@page), page numbers in the footer.
const path = require('path');
const { chromium } = require('playwright');

(async () => {
  const [html, pdf, title = ''] = process.argv.slice(2);
  if (!html || !pdf) { console.error('usage: node print.js in.html out.pdf [title]'); process.exit(1); }
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto('file://' + path.resolve(html), { waitUntil: 'load' });
  await page.emulateMedia({ media: 'print' });
  const esc = s => s.replace(/&/g, '&amp;').replace(/</g, '&lt;');
  await page.pdf({
    path: pdf,
    format: 'Letter',
    printBackground: true,
    preferCSSPageSize: true,
    displayHeaderFooter: true,
    headerTemplate: '<div></div>',
    footerTemplate:
      '<div style="width:100%;font:8px Liberation Sans,Arial,sans-serif;color:#555;padding:0 0.9in;display:flex;justify-content:space-between">' +
      `<span>${esc(title)}</span><span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>`,
  });
  await browser.close();
})().catch(e => { console.error(e); process.exit(1); });
