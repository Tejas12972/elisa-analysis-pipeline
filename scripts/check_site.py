"""Browser smoke test for the built site (``_site/``), including the in-browser app.

It serves the site locally, then uses a headless browser to check that:

* the landing page and every linked local page and image resolve;
* the stlite app boots, runs the pipeline on each built-in dataset, and shows the
  expected plate-QC verdict;
* no Python exception is rendered in the app.

Optionally it saves a screenshot for the README / landing page.

Usage::

    pip install playwright && playwright install chromium
    python scripts/check_site.py [--channel chrome] [--screenshot docs/img/app_screenshot.png]
"""

from __future__ import annotations

import argparse
import functools
import http.server
import sys
import threading
from pathlib import Path
from urllib.parse import urljoin, urlparse

from playwright.sync_api import Page, sync_playwright

REPO = Path(__file__).resolve().parents[1]
BOOT_TIMEOUT_MS = 300_000
EXPECTED = {  # dataset radio label prefix -> "Plates passing QC" metric
    "Demo study": "4/4",
    "Real: mouse IL-6": "3/4",
    "Real: SoftMax": "5/5",
}


class _QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        pass


def serve(root: Path) -> tuple[http.server.ThreadingHTTPServer, str]:
    handler = functools.partial(_QuietHandler, directory=str(root))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}/"


def check_links(page: Page, base: str) -> list[str]:
    """Every same-origin <a href> and <img src> on the landing page must return 200."""
    page.goto(base)
    urls = page.eval_on_selector_all(
        "a[href], img[src]", "els => els.map(e => e.getAttribute('href') || e.getAttribute('src'))"
    )
    problems = []
    for u in sorted(set(urls)):
        full = urljoin(base, u)
        if urlparse(full).netloc != urlparse(base).netloc:
            continue
        status = page.request.get(full).status
        if status != 200:
            problems.append(f"{u} -> HTTP {status}")
    return problems


def metric_value(page: Page) -> str:
    return page.locator("[data-testid='stMetricValue']").first.inner_text().strip()


def check_app(page: Page, base: str, screenshot: Path | None) -> list[str]:
    problems: list[str] = []
    page.goto(urljoin(base, "app/"))
    page.locator("[data-testid='stMetricValue']").first.wait_for(timeout=BOOT_TIMEOUT_MS)
    runtime = page.locator("[data-testid='stSidebar'] [data-testid='stCaptionContainer']").last
    print(f"  app booted · {runtime.inner_text()}")
    sidebar = page.locator("[data-testid='stSidebar']")
    for prefix, expected in EXPECTED.items():
        sidebar.get_by_text(prefix, exact=False).first.click()
        page.wait_for_function(
            "exp => { const m = document.querySelector(\"[data-testid='stMetricValue']\");"
            " return m && m.innerText.trim() === exp; }",
            arg=expected,
            timeout=120_000,
        )
        page.wait_for_timeout(1500)
        if page.locator("[data-testid='stException']").count():
            problems.append(f"{prefix}: Python exception rendered in the app")
        got = metric_value(page)
        print(f"  app · {prefix:<18} plates passing QC = {got}")
        if got != expected:
            problems.append(f"{prefix}: expected {expected}, got {got}")
        if screenshot is not None and prefix == "Demo study":
            screenshot.parent.mkdir(parents=True, exist_ok=True)
            page.screenshot(path=str(screenshot))
    page.get_by_role("tab", name="Download").click()
    try:
        with page.expect_download(timeout=30_000) as download:
            page.get_by_role("button", name="Download results.csv").click()
        print(f"  app · download works ({download.value.suggested_filename})")
    except Exception as exc:  # any failure mode is a check failure
        problems.append(f"download of results.csv failed: {type(exc).__name__}")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=REPO / "_site")
    ap.add_argument("--channel", default=None, help="e.g. 'chrome' to use an installed Chrome")
    ap.add_argument("--screenshot", type=Path, default=None)
    args = ap.parse_args()

    server, base = serve(args.site)
    errors: list[str] = []
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(channel=args.channel)
            page = browser.new_page(accept_downloads=True, viewport={"width": 1400, "height": 900})
            page.on("pageerror", lambda e: errors.append(f"page error: {e}"))
            problems = check_links(page, base)
            print(f"  landing page: {len(problems)} broken link(s)")
            problems += check_app(page, base, args.screenshot)
            browser.close()
    finally:
        server.shutdown()
    problems += errors
    for p in problems:
        print(f"FAIL {p}")
    print("site check:", "FAILED" if problems else "OK")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
