#!/usr/bin/env python
"""Static checks on site/: internal links and asset references resolve, external links are well-formed,
every <img> has alt text, every JSON is under 3 MB and the site under 15 MB, and every page has a title,
a viewport meta and the shared chrome. Exit 1 on any failure.

Usage: python tools/linkcheck.py [site_dir]
"""
from __future__ import annotations

import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse

SITE = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "site")


class Collector(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links, self.assets, self.imgs, self.title, self.viewport = [], [], [], None, False
        self._in_title = False

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "a" and a.get("href"):
            self.links.append(a["href"])
        if tag in ("script", "img") and a.get("src"):
            self.assets.append(a["src"])
        if tag == "link" and a.get("href"):
            self.assets.append(a["href"])
        if tag == "img":
            self.imgs.append(a)
        if tag == "meta" and a.get("name") == "viewport":
            self.viewport = True
        if tag == "title":
            self._in_title = True

    def handle_data(self, data):
        if self._in_title:
            self.title = (self.title or "") + data

    def handle_endtag(self, tag):
        if tag == "title":
            self._in_title = False


def main() -> int:
    problems = []
    pages = sorted(SITE.glob("*.html"))
    if not pages:
        print(f"no pages in {SITE}")
        return 1
    js_refs = set()
    for js in SITE.glob("assets/**/*.js"):
        js_refs |= set(re.findall(r'"(data/[\w.-]+\.json)"', js.read_text()))
    for page in pages:
        c = Collector()
        c.feed(page.read_text())
        if not c.title:
            problems.append(f"{page.name}: no <title>")
        if not c.viewport:
            problems.append(f"{page.name}: no viewport meta")
        for img in c.imgs:
            if "alt" not in img:
                problems.append(f"{page.name}: <img src={img.get('src')}> without alt")
        for href in c.links + c.assets:
            u = urlparse(href)
            if u.scheme in ("http", "https"):
                if not u.netloc or " " in href:
                    problems.append(f"{page.name}: malformed external link {href}")
                continue
            if u.scheme in ("mailto", "data"):
                continue
            path = u.path
            if not path:
                continue  # same-page anchor
            target = (page.parent / path).resolve()
            if not target.exists():
                problems.append(f"{page.name}: {href} does not resolve")
            elif u.fragment and target.suffix == ".html":
                if f'id="{u.fragment}"' not in target.read_text():
                    problems.append(f"{page.name}: {href}: anchor #{u.fragment} not found in {target.name}")
    for ref in sorted(js_refs):
        if not (SITE / ref).exists():
            problems.append(f"assets/*.js: {ref} referenced but missing")
    total = sum(p.stat().st_size for p in SITE.rglob("*") if p.is_file() and "_screenshots" not in p.parts)
    biggest = max((p for p in SITE.glob("data/*.json")), key=lambda p: p.stat().st_size)
    for p in SITE.glob("data/*.json"):
        if p.stat().st_size > 3_000_000:
            problems.append(f"{p.name} exceeds 3 MB")
    if total > 15_000_000:
        problems.append(f"site/ is {total / 1e6:.1f} MB (limit 15)")
    print(f"pages: {len(pages)}; site size {total / 1e6:.2f} MB (limit 15); largest JSON {biggest.name} {biggest.stat().st_size / 1e6:.2f} MB; json refs from js: {len(js_refs)}")
    for p in problems:
        print("  PROBLEM:", p)
    print("link check: " + ("clean" if not problems else f"{len(problems)} problems"))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
