-- packages
restart
debug Core
needsPackage "Graphs"

L = separate(" ", version#"packages") -- allPackages()
P = new MutableHashTable from {}
for pkg in L do P#pkg = readPackage pkg

G = digraph new HashTable from {}

for pkg in L do (
    print pkg;
    G = addVertices(G, {pkg});
    deps = P#pkg.PackageExports | P#pkg.PackageImports;
    if #deps == 0 then continue;
    G = addVertices(G, deps);
    G = addEdges'(G, (v -> {pkg, v}) \ deps);
    )

displayGraph G
q#vertices G

-- core sources
restart
debug Core
needsPackage "Graphs"

loads = get(Core#"source directory" | "loadsequence")
L = select("^\\w+\\.m2", "$&", loads)
--L = select(values loadedFiles, src -> match(regexQuote Core#"source directory", src));

P = new MutableHashTable from {}
for src in L do P#src = select("needs \"(.*\\.m2)\"", "$1", get(Core#"source directory"|src))

G = digraph(pairs P, EntryMode => "neighbors")

displayGraph G
#vertices G

S = new HashTable from for src in L list src => #P#src


-- interpreter
restart
needsPackage "Graphs"

p = Core#"source directory" | "../d/"
L = select("(\\w+\\.d+)", "$1", get concatenate("!ls ", p, "*.d*"))

P = new MutableHashTable from {}
for src in L do P#(first select("[^.]+", src)) = select("^use +(.*);", "$1", get(p | src))

G = digraph(pairs P, EntryMode => "neighbors")

displayGraph G
#vertices G

S = new HashTable from for src in L list src => #P#src
