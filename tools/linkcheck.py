#!/usr/bin/env python
"""Static checks on site/: internal links and asset references resolve, external links are well-formed,
every <img> has alt text, every JSON is under 3 MB and the site under 15 MB, every page has a title and
a viewport meta, every PAGES entry in assets/site.js has its file, favicon.ico exists, and every page carries the
head tags (og:image resolving to a local file, og:title, theme-color, apple-touch-icon; missing ones are WARNs
until every page has them). Exit 1 on any failure.

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
        self.metas, self.link_tags = [], []
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
        if tag == "meta":
            self.metas.append(a)
            if a.get("name") == "viewport":
                self.viewport = True
        if tag == "link":
            self.link_tags.append(a)
        if tag == "title":
            self._in_title = True

    def handle_data(self, data):
        if self._in_title:
            self.title = (self.title or "") + data

    def handle_endtag(self, tag):
        if tag == "title":
            self._in_title = False


HEAD_TAGS = ("og:image", "og:title", "theme-color", "apple-touch-icon")


def head_tags(c: Collector) -> dict:
    """The social / browser head tags of a page: the four values, empty when absent."""
    props = {m.get("property"): m.get("content") or "" for m in c.metas if m.get("property")}
    names = {m.get("name"): m.get("content") or "" for m in c.metas if m.get("name")}
    rels = {r for link in c.link_tags for r in (link.get("rel") or "").split()}
    return {"og:image": props.get("og:image", ""), "og:title": props.get("og:title", ""),
            "theme-color": names.get("theme-color", ""), "apple-touch-icon": "apple-touch-icon" in rels}


def og_local_path(content: str) -> str | None:
    """The site-relative path of an og:image URL, relative ("assets/...") or absolute (".../assets/...")."""
    path = urlparse(content).path
    i = path.find("assets/")
    return path[i:] if i >= 0 else (path.lstrip("/") or None)


def main() -> int:
    problems, warnings = [], []
    pages = sorted(SITE.glob("*.html"))
    if not pages:
        print(f"no pages in {SITE}")
        return 1
    js_refs = set()
    for js in SITE.glob("assets/**/*.js"):
        js_refs |= set(re.findall(r'"(data/[\w.-]+\.json)"', js.read_text()))
    head_complete = 0
    site_js = SITE / "assets" / "site.js"
    nav = re.findall(r'\["([\w.-]+\.html)",\s*"([^"]+)"\]', site_js.read_text()) if site_js.exists() else []
    if not nav:
        problems.append("assets/site.js: no PAGES entries found")
    for href, label in nav:
        if not (SITE / href).exists():
            problems.append(f"assets/site.js: PAGES entry {href} ({label}) has no file")
    if not (SITE / "favicon.ico").exists():
        problems.append("favicon.ico missing (run `make brand`)")
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
        tags = head_tags(c)
        missing = [k for k in HEAD_TAGS if not tags[k]]
        if missing:
            warnings.append(f"{page.name}: head tags missing: {', '.join(missing)}")
        else:
            head_complete += 1
        if tags["og:image"]:
            local = og_local_path(tags["og:image"])
            if not local or not (SITE / local).exists():
                problems.append(f"{page.name}: og:image {tags['og:image']} does not resolve to a file under site/")
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
    print(f"nav: {len(nav)} PAGES entries in assets/site.js; head tags: {head_complete} of {len(pages)} pages complete")
    for w in warnings:
        print("  WARN:", w)
    for p in problems:
        print("  PROBLEM:", p)
    print("link check: " + ("clean" if not problems else f"{len(problems)} problems"))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
