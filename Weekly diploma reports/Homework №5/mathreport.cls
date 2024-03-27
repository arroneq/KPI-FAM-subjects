\ProvidesClass{mathreport}[2023/11/09 university assignment report]
\NeedsTeXFormat{LaTeX2e}

% ------------------------------------------------------------------------------------------------------------
% Define The Document Class Options
% ------------------------------------------------------------------------------------------------------------

\LoadClass[a4paper,14pt]{extarticle}

% ------------------------------------------------------------------------------------------------------------
% Title Options
% ------------------------------------------------------------------------------------------------------------

% Title Name -------------------------------------------------------------------------------------------------

\newcommand{\ReportPreamble}[1]{
  \def\ReportPreamble{#1}
}

\newcommand{\ReportName}[1]{
  \def\ReportName{#1}
}

\newcommand{\ReportSubject}[1]{
  \def\ReportSubject{#1}
}

% Author's Info ----------------------------------------------------------------------------------------------

\newcommand{\AuthorName}[1]{
  \def\AuthorName{#1}
}

\newcommand{\AuthorInfo}[1]{
  \def\AuthorInfo{#1}
}

% Supervisor's Info ------------------------------------------------------------------------------------------

\newcommand{\SupervisorName}[1]{
  \def\SupervisorName{#1}
}

\newcommand{\SupervisorInfo}[1]{
  \def\SupervisorInfo{#1}
}

% The Longest Label ------------------------------------------------------------------------------------------

\newcommand{\TheLongestLabel}[1]{
  \def\TheLongestLabel{#1}
}

% ------------------------------------------------------------------------------------------------------------
% Load Necessary Packages
% ------------------------------------------------------------------------------------------------------------

\RequirePackage{fontspec}
\setmainfont{CMU serif}
% \setmainfont{Times New Roman}
\RequirePackage[english,ukrainian]{babel}

\RequirePackage{cmap} % make LaTeX PDF output copy-and-pasteable
\RequirePackage{indentfirst} % set an additional space before a paragraph at the begining of new section
\RequirePackage{microtype}
\parskip=2mm % space between paragraphs

\RequirePackage{mathtools,amsmath,amssymb,amsfonts}

\RequirePackage{amsthm}

\theoremstyle{plain}
\newtheorem{theorem}{\indent Теорема}
\newtheorem{claim}{\indent Твердження}
\newtheorem{lemma}{\indent Лема}
\newtheorem{corollary}{\indent Наслідок}

\theoremstyle{definition}
\newtheorem{definition}{\indent Означення}
\newtheorem{algorithm}{\indent Алгоритм}
\newtheorem{problem}{\indent Задача}
\newtheorem{example}{\indent Приклад}

\theoremstyle{remark}
\newtheorem*{remark}{\indent\textbf{Коментар}}

\renewenvironment{proof}{{\indent\bfseries Доведення.~}}{\qed}

\RequirePackage{setspace}
\RequirePackage{textcomp}
\RequirePackage{enumitem}
\RequirePackage{tabularray}

\RequirePackage{leftidx} % this package enables left subscripts and superscripts in math mode.

\RequirePackage{import} % for adding a file by path https://tex.stackexchange.com/questions/246/when-should-i-use-input-vs-include

\RequirePackage{geometry} 
\geometry{
    left=1.25cm,
    right=1.25cm,
    top=2cm,
    bottom=2cm
}

\RequirePackage[table,xcdraw,dvipsnames]{xcolor}
\RequirePackage{color}
% 1) tutorial about xcolor:  https://www.overleaf.com/learn/latex/Using_colours_in_LaTeX
% 2) huge tutorial about xcolor: https://latex-tutorial.com/color-latex/ 
% 3) RGB calculator: https://www.w3schools.com/colors/colors_rgb.asp

\RequirePackage[hidelinks]{hyperref} % [hidelinks] for clickable non-colored links
% \definecolor{linkcolor}{HTML}{0000FF}
% \definecolor{urlcolor}{HTML}{0000FF} 
% \definecolor{citecolor}{HTML}{0000FF} 
% \hypersetup{
%     pdfstartview=FitH, 
%     unicode=true, 
%     linkcolor=linkcolor, 
%     urlcolor=urlcolor, 
%     citecolor=citecolor, 
%     colorlinks=true
% }

\RequirePackage{graphicx}
\RequirePackage{wrapfig}
\RequirePackage{float}
\graphicspath{{Images/}} % path to images

\RequirePackage[bibstyle=gost-numeric,sorting=none]{biblatex} % bibliography settings
\addbibresource{citations.bib}
\RequirePackage[autostyle=false]{csquotes}

\RequirePackage{tikz} % package for creating graphics
\usetikzlibrary{arrows.meta}

\RequirePackage{pgfplots} % package for plotting data
\RequirePackage{pgfplotstable} % package for plotting data from tables
\pgfplotsset{compat=1.3}

\RequirePackage{lipsum}

\RequirePackage{listings}
\lstset{
    frame             = single,
    language          = Python,
    aboveskip         = 3mm,
    belowskip         = 3mm,
    columns           = flexible,
    basicstyle        = \linespread{1}\small\ttfamily,
    numbers           = left,
    numberstyle       = \tiny\color{gray},
    commentstyle      = \color{OliveGreen},
    % stringstyle      = \color{Mahogany},
    morestring        = [b]''',
    showstringspaces  = false,
    morekeywords      = {import, as, for, while, in, return, range},
    keywordstyle      = \bfseries\color{blue!50!black},
    emph              = {[1]import, as, for, while, in, return},
    emphstyle         = {[1]\bfseries\color{magenta!60!black}},
    breaklines        = true,
    breakatwhitespace = true,
    tabsize           = 4,
}