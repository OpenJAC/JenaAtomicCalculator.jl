#!/usr/bin/env python3
"""
html-to-flatodf.py  --  turn a simple HTML report into flat ODF (.fodt), for LibreOffice to make a PDF of.

    python3 tools/html-to-flatodf.py  in.html  out.fodt  [line-height]  [image-width-cm]
    soffice --headless --convert-to pdf out.fodt

WHY THIS EXISTS, so that nobody re-learns it:  LibreOffice can convert HTML to PDF directly, but its HTML
filter SILENTLY IGNORES line spacing and several other CSS properties, and it routes the file through the
Writer/Web document type, whose page setup is not the one the stylesheet asks for.  A 10-page report came
out as 19 pages that way.  Flat ODF is plain XML, converts as a real Writer document, and honours the
point size, the margins, the line spacing and the page breaks exactly.

WHAT IT HANDLES:  h1/h2/h3, p, ul/li, pre, tables with th/td (class "r" right-aligns a cell), figure with
img and figcaption, div class "box" (a bordered note), b/i/tt/sup/sub, br, and class "brk" to force a page
break.  Images are referenced by file name and must sit beside the .fodt.  Anything else is ignored rather
than mangled.

IT LIVES IN tools/ AND NOT IN AN apps/ FOLDER.  It was copied into two application folders before this was
noticed;  a generic format converter carries no application's information and belongs where both can reach
it.  See the directory contract in CLAUDE.md:  tools/ is the assistant's own, and is tracked.
"""
import re, sys, html
from html.parser import HTMLParser

SRC, DST = sys.argv[1], sys.argv[2]
LH = sys.argv[3] if len(sys.argv) > 3 else "100%"
IMGW = float(sys.argv[4]) if len(sys.argv) > 4 else 11.4

HEAD = '''<?xml version="1.0" encoding="UTF-8"?>
<office:document xmlns:office="urn:oasis:names:tc:opendocument:xmlns:office:1.0"
 xmlns:style="urn:oasis:names:tc:opendocument:xmlns:style:1.0"
 xmlns:text="urn:oasis:names:tc:opendocument:xmlns:text:1.0"
 xmlns:table="urn:oasis:names:tc:opendocument:xmlns:table:1.0"
 xmlns:fo="urn:oasis:names:tc:opendocument:xmlns:xsl-fo-compatible:1.0"
 xmlns:draw="urn:oasis:names:tc:opendocument:xmlns:drawing:1.0"
 xmlns:xlink="http://www.w3.org/1999/xlink"
 xmlns:svg="urn:oasis:names:tc:opendocument:xmlns:svg-compatible:1.0"
 office:version="1.3" office:mimetype="application/vnd.oasis.opendocument.text">
<office:font-face-decls>
 <style:font-face style:name="LS" svg:font-family="&apos;Liberation Serif&apos;" style:font-family-generic="roman"/>
 <style:font-face style:name="LM" svg:font-family="&apos;Liberation Mono&apos;" style:font-family-generic="modern"/>
</office:font-face-decls>
<office:automatic-styles>
 <style:page-layout style:name="pm1"><style:page-layout-properties fo:page-width="21cm"
   fo:page-height="29.7cm" fo:margin-top="1.9cm" fo:margin-bottom="1.6cm"
   fo:margin-left="1.9cm" fo:margin-right="1.9cm" style:print-orientation="portrait"/>
   <style:footer-style><style:header-footer-properties fo:min-height="0.5cm" fo:margin-top="0.4cm"/></style:footer-style>
 </style:page-layout>
 %STYLES%
</office:automatic-styles>
<office:master-styles><style:master-page style:name="Standard" style:page-layout-name="pm1">
 <style:footer><text:p text:style-name="FOOT">b26.01-average-summary-claude &#8212; page <text:page-number text:select-page="current">1</text:page-number></text:p></style:footer>
</style:master-page></office:master-styles>
<office:body><office:text>
%BODY%
</office:text></office:body></office:document>
'''

def para(name, size, align="justify", mb="0.10cm", mt="0cm", bold=False, italic=False,
         font="LS", brk=None, border=None, pad=None, keep=False, lh=None):
    pp = ['fo:text-align="%s"' % align, 'fo:margin-top="%s"' % mt, 'fo:margin-bottom="%s"' % mb,
          'fo:line-height="%s"' % (lh or LH)]
    if brk:   pp.append('fo:break-before="%s"' % brk)
    if border: pp.append('fo:border="%s"' % border)
    if pad:    pp.append('fo:padding="%s"' % pad)
    if keep:   pp.append('fo:keep-with-next="always"')
    tp = ['style:font-name="%s"' % font, 'fo:font-size="%s"' % size]
    if bold:   tp.append('fo:font-weight="bold"')
    if italic: tp.append('fo:font-style="italic"')
    return ('<style:style style:name="%s" style:family="paragraph">'
            '<style:paragraph-properties %s/><style:text-properties %s/></style:style>'
            % (name, " ".join(pp), " ".join(tp)))

STYLES = [
  para("P",    "11pt"),
  para("PB",   "11pt", brk="page"),
  para("H1",   "18pt", align="left", mb="0.05cm", bold=True, keep=True),
  para("SUB",  "12pt", align="left", mb="0.05cm", italic=True, keep=True),
  para("BY",   "9.5pt", align="left", mb="0.45cm"),
  para("H2",   "12.5pt", align="left", mt="0.42cm", mb="0.10cm", bold=True, keep=True),
  para("H2B",  "12.5pt", align="left", mt="0cm", mb="0.10cm", bold=True, keep=True, brk="page"),
  para("H3",   "11pt", align="left", mt="0.28cm", mb="0.04cm", bold=True, italic=True, keep=True),
  para("LI",   "11pt", mb="0.06cm"),
  para("PRE",  "8.5pt", align="left", mb="0.02cm", mt="0.02cm", font="LM", lh="100%"),
  para("CAP",  "9.5pt", mt="0.08cm", mb="0.30cm", lh="100%"),
  para("FIG",  "11pt", align="center", mt="0.22cm", mb="0.02cm"),
  para("BOX",  "10.5pt", mt="0.22cm", mb="0.28cm", border="0.5pt solid #777777", pad="0.15cm"),
  para("SMALL","9.5pt", mb="0.08cm", lh="100%"),
  para("TH",   "9.5pt", align="left", mb="0cm", mt="0cm", bold=True, lh="100%"),
  para("TD",   "9.5pt", align="left", mb="0cm", mt="0cm", lh="100%"),
  para("TDR",  "9.5pt", align="right", mb="0cm", mt="0cm", lh="100%"),
  para("FOOT", "8.5pt", align="right", mb="0cm"),
  '<style:style style:name="fr1" style:family="graphic"><style:graphic-properties '
  'style:wrap="none" style:vertical-pos="top" style:vertical-rel="text" '
  'style:horizontal-pos="center" style:horizontal-rel="paragraph"/></style:style>',
  '<style:style style:name="T-b" style:family="text"><style:text-properties fo:font-weight="bold"/></style:style>',
  '<style:style style:name="T-i" style:family="text"><style:text-properties fo:font-style="italic"/></style:style>',
  '<style:style style:name="T-m" style:family="text"><style:text-properties style:font-name="LM" fo:font-size="10pt"/></style:style>',
  '<style:style style:name="T-sup" style:family="text"><style:text-properties style:text-position="super 58%"/></style:style>',
  '<style:style style:name="T-sub" style:family="text"><style:text-properties style:text-position="sub 58%"/></style:style>',
  '<style:style style:name="tbl" style:family="table"><style:table-properties style:width="17.2cm" '
  'table:align="left" fo:margin-top="0.12cm" fo:margin-bottom="0.22cm"/></style:style>',
  '<style:style style:name="cellH" style:family="table-cell"><style:table-cell-properties '
  'fo:border="0.5pt solid #999999" fo:padding="0.035cm 0.09cm" fo:background-color="#eeeeee"/></style:style>',
  '<style:style style:name="cell" style:family="table-cell"><style:table-cell-properties '
  'fo:border="0.5pt solid #999999" fo:padding="0.035cm 0.09cm"/></style:style>',
]

class P(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.out = [];  self.buf = [];  self.mode = None;  self.pend_break = False
        self.in_tr = False;  self.ncols = 0;  self.tbl_rows = []; self.row = []
        self.inbox = False
        self.cellbuf = None;  self.cell_is_h = False; self.cell_align = "TD"
    # ---- helpers
    def esc(self, t):
        return (t.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;"))
    def emit(self, style, text):
        text = text.strip()
        if not text:  return
        if self.pend_break and style in ("P","H2"):
            style = "PB" if style == "P" else "H2B"
            self.pend_break = False
        self.out.append('<text:p text:style-name="%s">%s</text:p>' % (style, text))
    def flush(self):
        if self.mode and self.buf:
            self.emit(self.mode, "".join(self.buf))
        self.buf = [];  self.mode = None
    # ---- tags
    def handle_starttag(self, tag, attrs):
        a = dict(attrs);  cls = a.get("class","")
        if tag in ("h1","h2","h3","p","li","pre","figcaption","div"):
            self.flush()
            if tag == "div" and "box" in cls:       self.mode = "BOX"; self.inbox = True
            elif tag == "div":                       self.mode = None
            elif tag == "p" and self.inbox:          self.mode = "BOX"
            elif tag == "p" and "sub" in cls:        self.mode = "SUB"
            elif tag == "p" and "by" in cls:         self.mode = "BY"
            elif tag == "p" and "small" in cls:      self.mode = "SMALL"
            elif tag == "figcaption":                self.mode = "CAP"
            elif tag == "li":                        self.mode = "LI"
            elif tag == "pre":                       self.mode = "PRE"
            else:                                    self.mode = tag.upper()
            if "brk" in cls:  self.pend_break = True
        elif tag == "b":    (self.cellbuf if self.cellbuf is not None else self.buf).append('<text:span text:style-name="T-b">')
        elif tag == "i":    (self.cellbuf if self.cellbuf is not None else self.buf).append('<text:span text:style-name="T-i">')
        elif tag == "tt":   (self.cellbuf if self.cellbuf is not None else self.buf).append('<text:span text:style-name="T-m">')
        elif tag == "sup":  (self.cellbuf if self.cellbuf is not None else self.buf).append('<text:span text:style-name="T-sup">')
        elif tag == "sub":  (self.cellbuf if self.cellbuf is not None else self.buf).append('<text:span text:style-name="T-sub">')
        elif tag == "br":   (self.cellbuf if self.cellbuf is not None else self.buf).append('<text:line-break/>')
        elif tag == "img":
            self.flush()
            h = IMGW * 420.0/620.0
            self.out.append('<text:p text:style-name="FIG"><draw:frame draw:style-name="fr1" '
                'text:anchor-type="as-char" svg:width="%.2fcm" svg:height="%.2fcm">'
                '<draw:image xlink:href="%s" xlink:type="simple" xlink:show="embed" '
                'xlink:actuate="onLoad"/></draw:frame></text:p>' % (IMGW, h, a.get("src","")))
        elif tag == "table":
            self.flush();  self.tbl_rows = []
        elif tag == "tr":
            self.row = []
        elif tag in ("td","th"):
            self.cellbuf = [];  self.cell_is_h = (tag == "th")
            self.cell_align = "TDR" if "r" in a.get("class","") else "TD"
    def handle_endtag(self, tag):
        if tag in ("h1","h2","h3","p","li","pre","figcaption","div"):
            self.flush()
            if tag == "div":  self.inbox = False
        elif tag in ("b","i","tt","sup","sub"):
            (self.cellbuf if self.cellbuf is not None else self.buf).append('</text:span>')
        elif tag in ("td","th"):
            self.row.append((self.cell_is_h, self.cell_align, "".join(self.cellbuf).strip()))
            self.cellbuf = None
        elif tag == "tr":
            if self.row:  self.tbl_rows.append(self.row)
            self.row = []
        elif tag == "table":
            if not self.tbl_rows:  return
            n = max(len(r) for r in self.tbl_rows)
            s = ['<table:table table:name="t%d" table:style-name="tbl">' % len(self.out),
                 '<table:table-column table:number-columns-repeated="%d"/>' % n]
            for r in self.tbl_rows:
                s.append('<table:table-row>')
                for (ish, al, txt) in r:
                    s.append('<table:table-cell table:style-name="%s" office:value-type="string">'
                             '<text:p text:style-name="%s">%s</text:p></table:table-cell>'
                             % ("cellH" if ish else "cell", "TH" if ish else al, txt))
                for _ in range(n - len(r)):
                    s.append('<table:table-cell table:style-name="cell"/>')
                s.append('</table:table-row>')
            s.append('</table:table>')
            self.out.append("".join(s))
            self.tbl_rows = []
    def handle_data(self, d):
        d = self.esc(d)
        if self.cellbuf is not None:
            self.cellbuf.append(d)
        elif self.mode == "PRE":
            self.buf.append(d.replace("\n", "</text:p><text:p text:style-name=\"PRE\">")
                             .replace("  ", "<text:s text:c=\"2\"/>"))
        elif self.mode:
            self.buf.append(re.sub(r"\s+", " ", d))

src = open(SRC, encoding="utf-8").read()
body = src[src.index("<body>")+6 : src.index("</body>")]
p = P();  p.feed(body);  p.flush()
out = HEAD.replace("%STYLES%", "\n ".join(STYLES)).replace("%BODY%", "\n".join(p.out))
open(DST, "w", encoding="utf-8").write(out)
print("fodt written:", len(out), "bytes;", len(p.out), "blocks")
