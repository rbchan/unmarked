(TeX-add-style-hook "distsamp"
 (lambda ()
    (LaTeX-add-bibliographies
     "unmarked")
    (LaTeX-add-labels
     "eq:1"
     "eq:2"
     "fig:umfhist"
     "fig:predplot"
     "fig:detplot")
    (TeX-add-symbols
     '("package" 1)
     '("code" 1))
    (TeX-run-style-hooks
     "fullpage"
     "natbib"
     "Sweave"
     "fontenc"
     "OT1"
     "latex2e"
     "art10"
     "article"
     "a4paper")))

