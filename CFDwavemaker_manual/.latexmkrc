# latexmk configuration for the CFDwavemaker manual.
# minted (code listings) shells out to pygmentize, which requires
# -shell-escape. This file makes plain `latexmk` (and editors that call it,
# e.g. VS Code LaTeX Workshop) build correctly out of the box.
$pdf_mode = 1;
set_tex_cmds('-shell-escape -interaction=nonstopmode %O %S');
