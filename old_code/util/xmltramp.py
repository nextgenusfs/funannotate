#!/usr/bin/env python

"""xmltramp: Make XML documents easily accessible."""

__version__ = "2.18"
__author__ = "Aaron Swartz"
__credits__ = "Many thanks to pjz, bitsko, and DanC."
__copyright__ = "(C) 2003-2006 Aaron Swartz. GNU GPL 2."

__all__ = ['Element', 'Namespace', 'seed', 'load', 'parse', 'indent']

if not hasattr(__builtins__, 'True'):
    True, False = 1, 0


def isstr(f):
    return isinstance(f, (str, unicode))


def islst(f):
    return isinstance(f, (tuple, list))


def isint(f):
    return isinstance(f, (int, long))


def iselement(x, n):
    return isinstance(x, Element) and x._name == n


empty = {
    'http://www.w3.org/1999/xhtml': [
        'img', 'br', 'hr', 'meta', 'link', 'base', 'param', 'input', 'col', 'area'
    ],
    'http://www.w3.org/2005/Atom': ['link', 'category'],
}

indent = '\t'


def quote(x, elt=True):
    if elt and '<' in x and len(x) > 24 and x.find(']]>') == -1:
        return "<![CDATA[" + x + "]]>"
    else:
        x = x.replace('&', '&amp;').replace('<', '&lt;').replace(']]>', ']]&gt;')
    if not elt:
        x = x.replace('"', '&quot;')
    return x


def transpose(d):
    return dict((v, k) for k, v in d.items())


def delns(name):
    return name[1] if islst(name) and name[0] is None else name


class Element:

    def __init__(self, name, attrs=None, children=None, prefixes=None, value=None):
        self._name = delns(name)
        self._attrs = dict((delns(k), v) for k, v in attrs.items()) if attrs else {}
        self._dir = children or []

        prefixes = prefixes or {}
        self._prefixes = transpose(prefixes)
        self._dNS = prefixes.get(None, None) if prefixes else None

        if value is not None:
            if islst(value):
                self(*value)
            elif isinstance(value, dict):
                self(**value)
            elif value:
                self._dir.append(value)

    def __repr__(self, recursive=0, multiline=0, inner=0, inprefixes=None):
        def qname(name, inprefixes):
            if islst(name):
                if inprefixes[name[0]] is not None:
                    return inprefixes[name[0]] + ':' + name[1]
                else:
                    return name[1]
            else:
                return name

        def arep(a, inprefixes, addns=1):
            out = ''

            for p in self._prefixes:
                if p not in inprefixes:
                    if addns:
                        out += ' xmlns'
                        if self._prefixes[p]:
                            out += ':' + self._prefixes[p]
                        out += '="' + quote(p, False) + '"'
                    inprefixes[p] = self._prefixes[p]

            for k in a:
                out += ' ' + qname(k, inprefixes) + '="' + quote(a[k], False) + '"'

            return out

        inprefixes = inprefixes or {u'http://www.w3.org/XML/1998/namespace': 'xml'}

        if not inner:
            # need to call arep() first to set inprefixes
            attributes = arep(self._attrs, inprefixes, recursive)
            qn = qname(self._name, inprefixes)
            out = '<' + qn + attributes
        else:
            out = ''

        if not self._dir and islst(self._name) and self._name[1] in empty.get(self._name[0], []):
            if not inner:
                out += ' />'
            return out

        if not inner:
            out += '>'

        if recursive:
            indent_content = multiline and any(isinstance(x, Element) for x in self._dir)
            for x in self._dir:
                if indent_content:
                    out += '\n' + indent * recursive
                if isstr(x):
                    out += quote(x)
                elif isinstance(x, Element):
                    out += x.__repr__(recursive + 1, multiline, 0, inprefixes.copy())
                else:
                    raise TypeError("I wasn't expecting '%r'" % x)
            if indent_content:
                out += '\n' + indent * (recursive - 1)
        else:
            if self._dir:
                out += '...'

        if not inner:
            out += '</' + qn + '>'

        return out

    def __unicode__(self):
        return ''.join(unicode(x) for x in self._dir)

    def __str__(self):
        return self.__unicode__().encode('ascii', 'xmlcharrefreplace')

    def _mkns(self, n):
        return (self._dNS, n) if self._dNS and not islst(n) else n

    def _first(self, n, default=None):
        n = self._mkns(n)
        for x in self._dir:
            if iselement(x, n):
                return x
        return default

    def __getattr__(self, n):
        if n[0] == '_':
            raise AttributeError("Use foo['" + n + "'] to access the child element")
        child = self._first(n)
        if child is None:
            raise AttributeError('No child element named %r' % n)
        return child

    def __hasattr__(self, n):
        return bool(self._first(n))

    def __setattr__(self, n, v):
        if n[0] == '_':
            self.__dict__[n] = v
        else:
            self[n] = v

    def _getchild(self, n, default=None):
        return self._first(n, default)

    def __getitem__(self, n):
        if isint(n):  # d[1] == d._dir[1]
            return self._dir[n]
        elif isinstance(n, slice):
            # numerical slices
            if isint(n.start):
                return self._dir[n.start:n.stop]
            # d['foo':] == all <foo>s
            n = self._mkns(n.start)
            return [x for x in self._dir if iselement(x, n)]
        else:  # d['foo'] == first <foo>
            child = self._getchild(n)
            if child is None:
                raise KeyError(n)
            return child

    def __setitem__(self, n, v):
        if isint(n):  # d[1]
            self._dir[n] = v
        elif n is None:  # d[None] appends v
            self._dir.append(v)
        elif isinstance(n, slice):  # d['foo':] adds a new foo
            n = self._mkns(n.start)
            nv = Element(n, prefixes=transpose(self._prefixes), value=v)
            self._dir.append(nv)
        else:  # d["foo"] replaces first <foo> and dels rest
            n = self._mkns(n)
            nv = Element(n, prefixes=transpose(self._prefixes), value=v)
            idx = [i for (i, x) in enumerate(self._dir) if iselement(x, n)]
            if idx:
                self._dir[idx[0]] = nv
                for x in idx[1:]:
                    del self._dir[x]
            else:
                self._dir.append(nv)

    def __delitem__(self, n):
        if isint(n):
            del self._dir[n]
        elif isinstance(n, slice):  # delete all <foo>s
            n = self._mkns(n.start)
            for i in range(len(self)):
                if self[i]._name == n:
                    del self[i]
        elif isinstance(n, Element):  # delete exactly element n
            try:
                self._dir.remove(n)
            except ValueError:
                pass
        else:  # delete first foo
            n = self._mkns(n)
            for i in range(len(self)):
                if self._dir[i]._name == n:
                    del self._dir[i]
                    break

    def __call__(self, *_pos, **_set):
        if _set:  # e(a=1, b=2) sets attrs
            for k in _set:
                self._attrs[k] = _set[k]
            if not _pos:  # without positional args allows chaining
                return self
        if len(_pos) > 1:  # e('a', 1, 'b', 2) sets attrs pairwise
            for i in range(0, len(_pos), 2):
                self._attrs[_pos[i]] = _pos[i + 1]
            return self
        if len(_pos) == 1:  # e('attr') returns the attr value
            return self._attrs[_pos[0]]
        if len(_pos) == 0:  # e() returns the complete attr dict
            return self._attrs

    def __len__(self):
        return len(self._dir)

    def __contains__(self, n):
        # true if self has a child element n
        return self._getchild(n) is not None

    def _new(self, n, v=None):
        # appends element n with content v and returns it
        self[n:] = v
        return self._dir[-1]


class Namespace:

    def __init__(self, uri):
        self.__uri = uri

    def __getattr__(self, n):
        return (self.__uri, n)

    def __getitem__(self, n):
        return (self.__uri, n)

    def _prefix(self, p):
        return {p: self.__uri}


from collections import defaultdict
from xml.sax.handler import EntityResolver, DTDHandler, ContentHandler, ErrorHandler


class Seeder(EntityResolver, DTDHandler, ContentHandler, ErrorHandler):
    def __init__(self):
        self.stack = []
        self.ch = ''
        self.prefixes = defaultdict(list)
        ContentHandler.__init__(self)

    def startPrefixMapping(self, prefix, uri):
        self.prefixes[prefix].append(uri)

    def endPrefixMapping(self, prefix):
        self.prefixes[prefix].pop()

    def startElementNS(self, name, qname, attrs):
        ch = self.ch
        self.ch = ''
        if ch and not ch.isspace():
            self.stack[-1]._dir.append(ch)

        attrs = dict(attrs)
        newprefixes = dict((k, v[-1]) for k, v in self.prefixes.items())

        self.stack.append(Element(name, attrs, prefixes=newprefixes))

    def characters(self, ch):
        self.ch += ch

    def endElementNS(self, name, qname):
        ch = self.ch
        self.ch = ''
        if ch and not ch.isspace():
            self.stack[-1]._dir.append(ch)

        element = self.stack.pop()
        if self.stack:
            self.stack[-1]._dir.append(element)
        else:
            self.result = element

from xml.sax import make_parser
from xml.sax.handler import feature_namespaces


def seed(fileobj):
    seeder = Seeder()
    parser = make_parser()
    parser.setFeature(feature_namespaces, 1)
    parser.setContentHandler(seeder)
    parser.parse(fileobj)
    return seeder.result


def parse(text):
    from StringIO import StringIO
    return seed(StringIO(text))


def load(url):
    import urllib
    return seed(urllib.urlopen(url))


def unittest():
    parse('<doc>a<baz>f<b>o</b>ob<b>a</b>r</baz>a</doc>').__repr__(1, 1) == \
      '<doc>\n\ta<baz>\n\t\tf<b>o</b>ob<b>a</b>r\n\t</baz>a\n</doc>'

    assert str(parse("<doc />")) == ""
    assert str(parse("<doc>I <b>love</b> you.</doc>")) == "I love you."
    assert parse("<doc>\nmom\nwow\n</doc>")[0].strip() == "mom\nwow"
    assert str(parse('<bing>  <bang> <bong>center</bong> </bang>  </bing>')) == "center"
    assert unicode(parse('<doc>\xcf\x80</doc>')).encode('utf-8') == '\xcf\x80'

    d = Element('foo', attrs={'foo': 'bar'}, children=['hit with a', Element('bar'), Element('bar')])

    try:
        d._doesnotexist
        raise "ExpectedError", "but found success. Damn."
    except AttributeError:
        pass
    assert d.bar._name == 'bar'
    try:
        d.doesnotexist
        raise "ExpectedError", "but found success. Damn."
    except AttributeError:
        pass

    assert hasattr(d, 'bar')
    assert 'bar' in d
    assert 'foo' not in d
    assert d._getchild('baz') is None
    assert d._getchild('baz', 1) == 1

    assert d('foo') == 'bar'
    d(silly='yes')
    assert d('silly') == 'yes'
    assert d() == d._attrs

    assert d[0] == 'hit with a'
    d[0] = 'ice cream'
    assert d[0] == 'ice cream'
    del d[0]
    assert d[0]._name == "bar"
    assert len(d[:]) == len(d._dir)
    assert len(d[1:]) == len(d._dir) - 1
    assert len(d['bar':]) == 2
    d['bar':] = 'baz'
    assert len(d['bar':]) == 3
    assert d['bar']._name == 'bar'
    d[-1](sillier='no')
    assert d[-1].__repr__(1) == '<bar sillier="no">baz</bar>'
    del d[d[-1]]
    assert len(d['bar':]) == 2
    assert d[-1].__repr__(1) == '<bar></bar>'

    d = Element('foo')

    doc = Namespace("http://example.org/bar")
    bbc = Namespace("http://example.org/bbc")
    dc = Namespace("http://purl.org/dc/elements/1.1/")
    d = parse("""<doc version="2.7182818284590451"
      xmlns="http://example.org/bar"
      xmlns:dc="http://purl.org/dc/elements/1.1/"
      xmlns:bbc="http://example.org/bbc">
        <author>John Polk and John Palfrey</author>
        <dc:creator>John Polk</dc:creator>
        <dc:creator>John Palfrey</dc:creator>
        <bbc:show bbc:station="4">Buffy</bbc:show>
    </doc>""")

    assert repr(d) == '<doc version="2.7182818284590451">...</doc>'
    assert d.__repr__(1) == '<doc xmlns:bbc="http://example.org/bbc" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns="http://example.org/bar" version="2.7182818284590451"><author>John Polk and John Palfrey</author><dc:creator>John Polk</dc:creator><dc:creator>John Palfrey</dc:creator><bbc:show bbc:station="4">Buffy</bbc:show></doc>'
    assert d.__repr__(1, inner=1) == '<author xmlns:bbc="http://example.org/bbc" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns="http://example.org/bar">John Polk and John Palfrey</author><dc:creator xmlns:bbc="http://example.org/bbc" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns="http://example.org/bar">John Polk</dc:creator><dc:creator xmlns:bbc="http://example.org/bbc" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns="http://example.org/bar">John Palfrey</dc:creator><bbc:show xmlns:bbc="http://example.org/bbc" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns="http://example.org/bar" bbc:station="4">Buffy</bbc:show>'
    assert d.__repr__(1, 1) == '''<doc xmlns:bbc="http://example.org/bbc" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns="http://example.org/bar" version="2.7182818284590451">
	<author>John Polk and John Palfrey</author>
	<dc:creator>John Polk</dc:creator>
	<dc:creator>John Palfrey</dc:creator>
	<bbc:show bbc:station="4">Buffy</bbc:show>
</doc>'''

    assert repr(parse("<doc xml:lang='en' />")) == '<doc xml:lang="en"></doc>'

    assert str(d.author) == str(d['author']) == "John Polk and John Palfrey"
    assert d.author._name == doc.author
    assert str(d[dc.creator]) == "John Polk"
    assert d[dc.creator]._name == dc.creator
    assert str(d[dc.creator:][1]) == "John Palfrey"
    d[dc.creator] = "Me!!!"
    assert str(d[dc.creator]) == "Me!!!"
    assert len(d[dc.creator:]) == 1
    d[dc.creator:] = "You!!!"
    assert len(d[dc.creator:]) == 2

    assert d[bbc.show](bbc.station) == "4"
    d[bbc.show](bbc.station, "5")
    assert d[bbc.show](bbc.station) == "5"

    e = Element('e')
    e.c = '<img src="foo">'
    assert e.__repr__(1) == '<e><c>&lt;img src="foo"></c></e>'
    e.c = '2 > 4'
    assert e.__repr__(1) == '<e><c>2 > 4</c></e>'
    assert e.__repr__(1, inner=1) == '<c>2 > 4</c>'
    e.c = 'CDATA sections are <em>closed</em> with ]]>.'
    assert e.__repr__(1) == '<e><c>CDATA sections are &lt;em>closed&lt;/em> with ]]&gt;.</c></e>'
    e.c = parse('<div xmlns="http://www.w3.org/1999/xhtml">i<br /><span></span>love<br />you</div>')
    assert e.__repr__(1) == '<e><c><div xmlns="http://www.w3.org/1999/xhtml">i<br /><span></span>love<br />you</div></c></e>'

    e = Element('e')
    e('c', 'that "sucks"')
    assert e.__repr__(1) == '<e c="that &quot;sucks&quot;"></e>'

    assert quote("]]>") == "]]&gt;"
    assert quote('< dkdkdsd dkd sksdksdfsd fsdfdsf]]> kfdfkg >') == '&lt; dkdkdsd dkd sksdksdfsd fsdfdsf]]&gt; kfdfkg >'

    assert parse('<x a="&lt;"></x>').__repr__(1) == '<x a="&lt;"></x>'
    assert parse('<a xmlns="http://a"><b xmlns="http://b"/></a>').__repr__(1) == '<a xmlns="http://a"><b xmlns="http://b"></b></a>'

if __name__ == '__main__':
    unittest()