import shutil
from pathlib import Path

source = Path("logs_v2/plots/combined_figure.svg")
latex_dir = Path("latex_exports")
latex_dir.mkdir(exist_ok=True)

# Copy + output LaTeX include code
shutil.copy(source, latex_dir / "figure.svg")

with open(latex_dir / "figure_tex_snippet.txt", "w") as f:
    f.write(r"""
% LaTeX snippet
\usepackage{graphicx}
\begin{figure}[htbp]
    \centering
    \includegraphics[width=0.9\textwidth]{figure.svg}
    \caption{Multi-panel visualization of peptide–membrane simulation.}
\end{figure}
    """.strip())

print("[✔] LaTeX export complete — see 'latex_exports/'")

