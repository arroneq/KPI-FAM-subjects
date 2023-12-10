\ProvidesClass{englishreport}[university english assignment report]
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

% The Longest Label ------------------------------------------------------------------------------------------

\newcommand{\TheLongestLabel}[1]{
  \def\TheLongestLabel{#1}
}

% ------------------------------------------------------------------------------------------------------------
% Load Necessary Packages
% ------------------------------------------------------------------------------------------------------------

\RequirePackage{fontspec}
% \setmainfont{CMU serif}
\setmainfont{Times New Roman}
\RequirePackage[ukrainian,english]{babel}

\RequirePackage{cmap} % make LaTeX PDF output copy-and-pasteable
\RequirePackage{indentfirst} % set an additional space before a paragraph at the begining of new section
\RequirePackage{microtype}

\linespread{1.5} % space between text lines
\setlength{\parindent}{1.25cm}
\RequirePackage{titletoc}

\RequirePackage{mathtools,amsmath,amssymb,amsfonts}

\RequirePackage{setspace}
\RequirePackage{textcomp}
\RequirePackage{enumitem}
\RequirePackage{tabularray}

\RequirePackage{leftidx} % this package enables left subscripts and superscripts in math mode.

\RequirePackage{import} % for adding a file by path https://tex.stackexchange.com/questions/246/when-should-i-use-input-vs-include

\RequirePackage{geometry} 
\geometry{
    left=2.5cm,
    right=1cm,
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

% prohibit the use of transfers
\tolerance=1
\emergencystretch=\maxdimen
\hyphenpenalty=10000
\hbadness=10000

\RequirePackage[style=apa]{biblatex} % bibliography settings
% \DeclareLanguageMapping{english}{english-apa}
\RequirePackage[autostyle=false]{csquotes}
\addbibresource{citations.bib}