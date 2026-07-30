doc ///
Node
  Key
    HolonomicSystems
  Headline
    Examples of Holonomic D-modules
  Description
    Text
      This package contains various classes of examples of holonomic systems represented as D-modules.
      It also includes an implementation of the canonical series solution algorithm for a regular
      holonomic system by [@HREF("https://mathscinet.ams.org/mathscinet/pdf/1734566.pdf","SST")@].
    Tree
      :Some examples of D-modules
        @TOH "gkz"@
	@TOH "eulerOperators"@
	@TOH "toricIdealPartials"@
	@TOH "AppellF1"@
      :@TOH "Canonical Series Tutorial"@
        @TOH "distraction"@
	@TOH "cssExpts"@
	@TOH "cssExptsMult"@
	@TOH "isTorusFixed"@
	@TOH "solveFrobeniusIdeal"@
      :Differential Operators
        @TOH "diffOps"@
  References
    [@HREF("https://link.springer.com/book/10.1007/978-3-662-04112-3","SST")@]
    M. Saito, B. Sturmfels, and N. Takayama. {\em Gröbner Deformations of Hypergeometric Differential Equations}.
    Volume 6 of {\em Algorithms and Computation in Mathematics}. Springer, 2000.
  Subnodes
   "gkz"
   "eulerOperators"
   "toricIdealPartials"
   "AppellF1"
   "Canonical Series Tutorial"
   "diffOps"
///

end--

restart
uninstallPackage "HolonomicSystems"
installPackage "HolonomicSystems"
