
import subprocess



latex_str = r"""
\documentclass{article}
\usepackage{multirow}
\begin{document}
\begin{center}
\begin{tabular}{ |c|c|c| } 
\hline
g / n & 0 & 1 \\
\hline
\multirow{2}{4em}{0} 
& k=12: 1 & k=14: 1 \\ 
& k=13: 2 &  \\ 
\hline
\multirow{2}{4em}{1} 
& k=12: 1 & 0 \\ 
& k=13: 2 &  \\ 
\hline
\multirow{2}{4em}{1} 
& k=12: 1 & ? \\ 
& k=13: 2 &  \\ 
\hline
\end{tabular}
\end{center}
\end{document}
"""

# Write LaTeX string to file
with open("table.tex", "w") as f:
    f.write(latex_str)


subprocess.run([
    "pdflatex",
    "-output-directory", "/mnt/c/Users/ps200/Desktop/Semesterarbeit/Table",
    "table.tex"
])
