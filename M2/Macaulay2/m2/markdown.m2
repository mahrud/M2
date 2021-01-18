--  Copyright 2020 by Mahrud Sayrafi
-----------------------------------------------------------------------------
-- markdown output
-- See https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet
-- TODO: task lists
-----------------------------------------------------------------------------

-* Note: users can modify the behavior of markdown on individual types in
-- order to adjust its output to different markdown processors.
-- For example, to specify the layout and excerpt:
markdown HEAD := x -> concatenate("---\n",
    "layout: package", newline,
    "excerpt: [headline]", newline,
    apply(x, markdown), "---\n")
-- or to skip wrapping BODY in Liquid raw tags:
markdown BODY := x -> concatenate(shorten apply(x, markdown))
*-

needs "format.m2"

-- YAML Ain't Markup Language, but it is a MarkUpType! See yaml.org
-- TODO: add parser and generator
YAML = new IntermediateMarkUpType of HypertextContainer

-----------------------------------------------------------------------------
-- Local utilities
-----------------------------------------------------------------------------

-- Trim extra newlines
shorten := s -> replace("\n+$", "\n", concatenate(s, newline))

-- Parse attributes
parseAttr := (T, x) -> first override(T.Options, toSequence x)

-- FIXME
markdownLiteral := identity

-----------------------------------------------------------------------------
-- Setup default rendering
-----------------------------------------------------------------------------

-- Rendering by concatenation of rendered inputs
setupRenderer(markdown, concatenate, Hypertext)

-----------------------------------------------------------------------------
-- (markdown, MarkUpType) methods
-----------------------------------------------------------------------------

-- Markdown-native hypertext objects must begin with a YAML object rather than HEAD
markdown YAML  := x -> concatenate("---\n", between_newline apply(toList x, markdown), "\n---\n")

-- See hypertext.m2 for the full list
markdown HTML  := x -> concatenate(apply(x, markdown))
markdown HEAD  := x -> concatenate("---\n",   apply(x, markdown), "\n---\n")
markdown TITLE := x -> concatenate("title: ", apply(x, markdown))

markdown BODY  := x -> concatenate(
    "{% raw %}", shorten apply(x, markdown), "{% endraw %}")
markdown SPAN  := x -> concatenate(apply(x, markdown), newline) -- TODO: should this have a newline?
markdown PARA  := x -> concatenate(apply(x, markdown), 2:newline)
markdown DIV   := x -> (
    c := (parseAttr(DIV, x))#"class";
    concatenate(shorten apply(x, markdown), if c =!= null then ("{:.", c, "}"), 2:newline))

markdown LABEL :=
markdown SMALL :=
markdown SUB   :=
markdown SUP   := html

markdown LINK  :=
markdown MENU  :=
markdown META  :=
markdown STYLE :=
markdown Option  :=
markdown Nothing := x -> ""
markdown COMMENT := x -> concatenate("<!--", apply(x, markdown), "-->")

-- TODO: implement other literal character changes based on this reference:
--       https://daringfireball.net/projects/markdown/syntax#code
markdown LITERAL :=
markdown String  := identity
markdown TEX     := x -> concatenate apply(x, markdown)

-- Headers
markdown HEADER1 := x -> concatenate(newline, 1:"#", " ", apply(x, markdown), newline)
markdown HEADER2 := x -> concatenate(newline, 2:"#", " ", apply(x, markdown), newline)
markdown HEADER3 := x -> concatenate(newline, 3:"#", " ", apply(x, markdown), newline)
markdown HEADER4 := x -> concatenate(newline, 4:"#", " ", apply(x, markdown), newline)
markdown HEADER5 := x -> concatenate(newline, 5:"#", " ", apply(x, markdown), newline)
markdown HEADER6 := x -> concatenate(newline, 6:"#", " ", apply(x, markdown), newline)

-- Emphasis
-- TODO: markdown STRIKE -* ~~Scratch this.~~ *-
markdown EM     :=
markdown ITALIC := x -> concatenate("_",  apply(x, markdown), "_")  -- *asterisks* or _underscores_
markdown BOLD   :=
markdown STRONG := x -> concatenate("**", apply(x, markdown), "**") -- **asterisks** or __underscores__.

-- Lists
lvl := -1              -- nesting level
ctr := new MutableList -- counters for ordered list
pushLvl :=  n     -> (lvl = lvl + n; n)
popLvl  := (n, s) -> (lvl = lvl - n; s)
markdown OL := x -> popLvl(pushLvl 1, (ctr#lvl =  0; shorten concatenate(newline, apply(x, markdown))))
markdown UL := x -> popLvl(pushLvl 1, (ctr#lvl = -1; shorten concatenate(newline, apply(x, markdown))))
markdown LI := x -> (
    pre := concatenate(lvl:"  ", if ctr#lvl < 0 then  "- " else (ctr#lvl = ctr#lvl + 1; toString ctr#lvl | ". "));
    shorten concatenate(pre, apply(x, markdown)))

-- Description lists
markdown DL := x -> shorten concatenate(apply(x, markdown));
markdown DT := x -> shorten concatenate(apply(x, markdown));
markdown DD := x -> shorten concatenate(lvl:"   ", ": ", apply(x, markdown))

-- Links
-- (markdown, TOH) defined in format.m2
markdown TO     := x -> markdown TO2{tag := x#0, format tag | if x#?1 then x#1 else ""}
markdown TO2    := x -> (
    (tag, name) := (getPrimaryTag fixup x#0, x#1);
    title := if (excerpt := headline tag) =!= null then (1, format markdownLiteral excerpt) else "";
    if isMissingDoc tag or isUndocumented tag then concatenate(markdown TT name, " (missing<!-- tag: ", toString tag.Key, " -->)") else
    concatenate("[", name, "](/", markdownBasename tag, title, ")"))
markdown ANCHOR := html
markdown HREF   := x -> concatenate("[", markdown last x, "](", toURL first x, ")")

-- Images
markdown IMG := x -> (
    (o, cn) := override(IMG.Options, toSequence x);
    if o#"alt" === null then error ("IMG item is missing alt attribute");
    concatenate("![", format o#"alt", "](", toURL o#"src", ")"))

-- Code and Syntax Highlighting
markdown PRE  := x -> markdown concatenate x -* TODO: syntax highlighting *-
markdown TT   := x -> concatenate("`", apply(x, markdown), "`")
markdown CODE := x -> concatenate(
    newline, "```M2", newline, shorten apply(x, markdown), "```", newline)

-- Blockquotes
-- TODO: new lines should also start with >
markdown BLOCKQUOTE := x -> concatenate("> ", apply(x, markdown), 2:newline)

-- Horizontal Rule
markdown HR := x -> "---"

-- Line Breaks
markdown BR := x -> "<br/>\n"

-- Tables
-- | Tables        | Are           | Cool  |
-- | ------------- |:-------------:| -----:|
-- | col 3 is      | right-aligned | $1600 |
-- | col 2 is      | centered      |   $12 |
-- | zebra stripes | are neat      |    $1 |
markdown TABLE := x -> (
    c := (parseAttr(TABLE, x))#"class";
    concatenate(shorten markdown new CODE from toSequence x, if c =!= null then ("{:.", c, "}\n")))
markdown TR    :=
markdown TD    := x -> concatenate(apply(x, markdown), newline)
