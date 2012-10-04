(TeX-add-style-hook "spp-dist"
 (lambda ()
    (LaTeX-add-bibliographies
     "unmarked")
    (LaTeX-add-labels
     "fig:swiss"
     "fig:ef"
     "fig:psi1"
     "fig:predict"
     "fig:issj")
    (TeX-run-style-hooks
     "geometry"
     "2cm}"
     "color"
     "verbatim"
     "fullpage"
     "natbib"
     "authoryear"
     "round"
     "Sweave"
     "fontenc"
     "OT1"
     "latex2e"
     "art10"
     "article"
     "a4paper")))

