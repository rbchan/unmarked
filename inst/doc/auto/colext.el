(TeX-add-style-hook "colext"
 (lambda ()
    (LaTeX-add-labels
     "fig:sim"
     "fig:yearlysim"
     "fig:gof"
     "fig:cov")
    (TeX-add-symbols
     "rf")
    (TeX-run-style-hooks
     "graphicx"
     "setspace"
     "indentfirst"
     "lineno"
     "float"
     "framed"
     "url"
     "amssymb"
     "amsmath"
     "fullpage"
     "Sweave"
     "fontenc"
     "OT1"
     "latex2e"
     "art12"
     "article"
     "12pt")))

