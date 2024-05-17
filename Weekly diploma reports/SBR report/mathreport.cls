\ProvidesClass{mathreport}[2024/04/20 university assignment report]
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

\RequirePackage[onehalfspacing]{setspace}
\setlength{\parindent}{1.25cm}

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
    left=2cm,
    right=2cm,
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

\RequirePackage[
    bibstyle=gost-numeric,
    sorting=none,
    maxbibnames=99
]{biblatex} % bibliography settings

\addbibresource{citations.bib}
\RequirePackage[autostyle=false]{csquotes}

\RequirePackage{tikz} % package for creating graphics
\usetikzlibrary{arrows.meta}

\RequirePackage{pgfplots} % package for plotting data
\RequirePackage{pgfplotstable} % package for plotting data from tables
\pgfplotsset{compat=1.3}

\definecolor{forestgreen}{RGB}{34,139,34}
\definecolor{orangered}{RGB}{239,134,64}
\definecolor{lightred}{rgb}{1,0.4,0.5}
\definecolor{orange}{rgb}{1,0.45,0.13}	
\definecolor{darkblue}{rgb}{0.0,0.0,0.6}
\definecolor{lightblue}{rgb}{0.1,0.57,0.7}
\definecolor{gray}{rgb}{0.4,0.4,0.4}
\definecolor{lightgray}{rgb}{0.95,0.95,0.95}
\definecolor{editorGray}{rgb}{0.95,0.95,0.95}
\definecolor{editorOcher}{rgb}{1,0.5,0}
\definecolor{chaptergrey}{rgb}{0.6,0.6,0.6}
\definecolor{editorGreen}{rgb}{0,0.5,0}
\definecolor{olive}{rgb}{0.17,0.59,0.20}
\definecolor{brown}{rgb}{0.69,0.31,0.31}
\definecolor{purple}{rgb}{0.38,0.18,0.81}

\colorlet{commentcolour}{green!50!black}
\colorlet{stringcolour}{red!60!black}
\colorlet{keywordcolour}{blue!90!black}
\colorlet{keywordcolour2}{magenta!80!black}
\colorlet{commandcolour}{blue!60!black}
\colorlet{numpycolour}{yellow!60!black}
\colorlet{exceptioncolour}{yellow!50!red}
\colorlet{literatecolour}{gray!80!white}
\colorlet{promptcolour}{green!50!black}
\colorlet{specmethodcolour}{violet}

\RequirePackage{listings}

\lstdefinestyle{mypython}{
    frame            = single,
    keepspaces       = true,
    language         = python,
    numbers          = left,
    numberstyle      = \tiny,
    % numbersep        = 5pt,
    showtabs         = true,
    tabsize          = 4,
    basicstyle       = \ttfamily\footnotesize,
    stringstyle      = \color{stringcolour},
    showstringspaces = false,
    upquote          = true,
    morecomment      = [s]{"""}{"""},
    commentstyle     = \color{commentcolour}\slshape,
    rulecolor        = \color{black!40},
    backgroundcolor  = \color{lightgray!60!white},
    breakindent      = .5\textwidth,
    breaklines       = true,
    % prebreak         = \raisebox{0ex}[0ex][0ex]{\ensuremath{\hookleftarrow}},
    % otherkeywords    = {\%, \}, \{, \&, \|},
    alsoletter       = {1234567890},
    keywordstyle     = \color{keywordcolour}\bfseries,
    emph             = {[1]and, class, def, global, in, lambda, not, or, pass, print, raise, try, assert, with},
    emphstyle        = {[1]\color{keywordcolour}\bfseries},
    emph             = {[2]True, False, None},
    emphstyle        = {[2]\color{keywordcolour}},
    emph             = {[3]object, type, isinstance, copy, deepcopy, zip, enumerate, reversed, list, set, dict, tuple, xrange, append, add, execfile, real, imag, reduce, str, repr},
    emphstyle        = {[3]\color{commandcolour}},
    emph             = {[4]Exception, NameError, IndexError, SyntaxError, TypeError, ValueError, OverflowError, ZeroDivisionError},
    emphstyle        = {[4]\color{exceptioncolour}\bfseries},
    emph             = {[5]range, ode, fsolve, sqrt, exp, sin, cos,arctan, arctan2, arccos, pi,  array, norm, solve, dot, range, len, max, sum, flatten, shape, reshape, find, any, all, abs, plot, linspace, legend, quad, polyval, polyfit, hstack, concatenate, vstack, column_stack, empty, zeros, ones, rand, vander, grid, pcolor, eig, eigs, eigvals, svd, qr, tan, det, logspace, roll, min, mean, cumsum, cumprod, diff, vectorize, lstsq, cla, eye, xlabel, ylabel, squeeze},
    emphstyle        = {[5]\color{numpycolour}},
    emph             = {[6]__init__,__add__,__mul__,__div__,__sub__,__call__,__getitem__,__setitem__,__eq__,__ne__,__nonzero__,__rmul__,__radd__,__repr__,__str__,__get__,__truediv__,__pow__,__name__,__future__,__all__},
    emphstyle        = {[6]\color{specmethodcolour}},
    emph             = {[7]break, continue, yield, del, elif, else, except, exec, finally, for, from, if, import, return},
    emphstyle        = {[7]\color{keywordcolour2}},
}

\lstdefinestyle{myR}{
    frame            = single,
    keepspaces       = true,
    language         = R,
    numbers          = left,
    numberstyle      = \tiny,
    % numbersep        = 5pt,
    showtabs         = true,
    tabsize          = 4,
    basicstyle       = \ttfamily\footnotesize,
    stringstyle      = \color{stringcolour},
    showstringspaces = false,
    upquote          = true,
    morecomment      = [s]{"""}{"""},
    commentstyle     = \color{commentcolour}\slshape,
    rulecolor        = \color{black!40},
    backgroundcolor  = \color{lightgray!60!white},
    breakindent      = .5\textwidth,
    breaklines       = true,
    % prebreak         = \raisebox{0ex}[0ex][0ex]{\ensuremath{\hookleftarrow}},
    % otherkeywords    = {\%, \}, \{, \&, \|},
    alsoletter       = {1234567890},
    keywordstyle     = {},
    emph             = {[1]and, or, not, in, function, print, raise, try, with},
    emphstyle        = {[1]\color{keywordcolour}\bfseries},
    emph             = {[2]True, False, None},
    emphstyle        = {[2]\color{keywordcolour}},
    emph             = {[3]list, set, seed, c, array, read, csv},
    emphstyle        = {[3]\color{commandcolour}},
    emph             = {[4]Exception, NameError, IndexError, SyntaxError, TypeError, ValueError, OverflowError, ZeroDivisionError},
    emphstyle        = {[4]\color{exceptioncolour}\bfseries},
    emph             = {[5]summary, head, tail, ncol, nrow, cbind, rbind, sqrt, exp, sin, cos, arctan, arccos, pi, dnorm, dpois, dexp, dgamma, rnorm, rpois, rexp, rgamma, range, len, max, min, sum, abs, mean},
    emphstyle        = {[5]\color{numpycolour}},
    emph             = {[9]return, if, else, while, for, repeat, next, break, return, switch, try, catch},
    emphstyle        = {[9]\color{keywordcolour2}},
}